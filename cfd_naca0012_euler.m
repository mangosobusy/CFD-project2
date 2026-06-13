%% cfd_naca0012_euler.m
% 基于第一题生成的 NACA0012 O 型贴体网格，求解二维可压缩 Euler 方程。
% 方法：半离散迎风有限体积 + Rusanov 通量 + MUSCL-TVD minmod 限制器 + 三阶 SSP Runge-Kutta。
% 功能：计算翼型跨声速绕流，对比不同时间步收敛效果，输出流场与气动力结果。
clear; clc; close all;

%% 1. 图形与基本参数
set(0,'DefaultFigureColor','w');
try
    set(0,'DefaultAxesFontName','SimSun');
    set(0,'DefaultTextFontName','SimSun');
catch
end
set(0,'DefaultAxesFontSize',10.5);

gamma = 1.4;
Minf  = 0.80;              % 可改为 0.40 或 0.80
alphaDeg = 1.25;           % 攻角，单位：deg
rhoInf = 1.0;
pInf   = 1.0;
CFL    = 0.5;             % 稳定 CFL 数。若发散，可降至 0.20~0.35
maxIter = 1500;            % 计算可先取，正式结果建议继续增大
tol     = 1.0e-6;
useLimiter = true;         % true：MUSCL-TVD；false：一阶迎风 Rusanov
mainTimeMode = 'local';    % 'local' 或 'global'

doTimeStepComparison = true;   % 是否比较局部/全局时间步的残差下降
compareIter = 350;
plotEvery = 100;

%% 2. 读取或生成网格
if ~exist('naca0012_grid.mat','file')
    fprintf('未发现 naca0012_grid.mat，先执行第一题网格生成程序。\n');
    run('cfd_naca0012_grid.m');
end
load('naca0012_grid.mat','X','Y','params');
geom = build_geometry(X,Y);

fprintf('读取网格：nI=%d, nJ(cell)=%d, 远场半径约 %.1f c\n', ...
    geom.nI,geom.nJ,params.Rfar/params.c);

%% 3. 远场状态与初值
alpha = alphaDeg*pi/180;
aInf = sqrt(gamma*pInf/rhoInf);
Vinf = Minf*aInf;
uInf = Vinf*cos(alpha);
vInf = Vinf*sin(alpha);
EInf = pInf/(gamma-1) + 0.5*rhoInf*(uInf^2+vInf^2);
UInf = [rhoInf, rhoInf*uInf, rhoInf*vInf, EInf];
free = struct('rho',rhoInf,'p',pInf,'u',uInf,'v',vInf,'E',EInf, ...
              'U',UInf,'M',Minf,'alphaDeg',alphaDeg,'gamma',gamma);

U0 = zeros(geom.nI,geom.nJ,4);
for k = 1:4
    U0(:,:,k) = UInf(k);
end

%% 4. 比较局部时间步长与全局时间步长的收敛速度
if doTimeStepComparison
    fprintf('\n开始局部/全局时间步残差比较，比较步数=%d。\n',compareIter);
    [~,histLocal]  = advance_euler(U0,geom,free,gamma,CFL,compareIter,tol,useLimiter,'local',plotEvery,false);
    [~,histGlobal] = advance_euler(U0,geom,free,gamma,CFL,compareIter,tol,useLimiter,'global',plotEvery,false);

    figC = figure('Name','局部与全局时间步收敛比较');
    semilogy(histGlobal.iter,histGlobal.res,'LineWidth',1.2); hold on;
    semilogy(histLocal.iter,histLocal.res,'LineWidth',1.2);
    grid on; xlabel('迭代步数'); ylabel('密度相对变化残差');
    legend('全局时间步','局部时间步','Location','southwest');
    title('局部时间步与全局时间步收敛比较');
    safe_export(figC,'naca0012_timestep_comparison.png');
end

%% 5. 主计算
fprintf('\n开始主计算：M=%.2f, alpha=%.2f deg, CFL=%.2f, timeMode=%s\n', ...
    Minf,alphaDeg,CFL,mainTimeMode);
[U,history] = advance_euler(U0,geom,free,gamma,CFL,maxIter,tol,useLimiter,mainTimeMode,plotEvery,true);

%% 6. 后处理：流场云图、压力系数、气动力系数
[rho,u,v,p,Ma] = primitive_from_U_field(U,gamma);
qInf = 0.5*rhoInf*Vinf^2;
CpCell = (p-pInf)/qInf;

fig1 = figure('Name','马赫数云图');
pcolor(geom.xc,geom.yc,Ma); shading interp; axis equal tight;
colorbar; xlabel('x/c'); ylabel('y/c');
title(sprintf('NACA0012 马赫数云图，M_\\infty=%.2f，\\alpha=%.2f°',Minf,alphaDeg));
safe_export(fig1,'naca0012_mach_contour.png');

fig2 = figure('Name','压力系数云图');
pcolor(geom.xc,geom.yc,CpCell); shading interp; axis equal tight;
colorbar; xlabel('x/c'); ylabel('y/c');
title('NACA0012 压力系数云图');
safe_export(fig2,'naca0012_cp_contour.png');

% 壁面压力系数。采用近壁单元压力近似壁面压力。
xWall = X(:,1); yWall = Y(:,1);
pWall = p(:,1);
CpWall = (pWall-pInf)/qInf;

fig3 = figure('Name','翼型表面压力系数');
idxUpper = yWall >= -1.0e-10;
idxLower = yWall <  -1.0e-10;
plot(xWall(idxUpper),CpWall(idxUpper),'o-','LineWidth',1.0,'MarkerSize',3); hold on;
plot(xWall(idxLower),CpWall(idxLower),'s-','LineWidth',1.0,'MarkerSize',3);
set(gca,'YDir','reverse'); grid on; xlabel('x/c'); ylabel('C_p');
legend('上翼面','下翼面','Location','best');
title('NACA0012 表面压力系数分布');
safe_export(fig3,'naca0012_wall_cp.png');

% 压力积分估算气动力系数。nWall 为流体单元指向壁面的外法向。
Fx = 0.0; Fy = 0.0;
for i = 1:geom.nI
    nvec = squeeze(geom.nWall(i,:))';
    Fx = Fx + pWall(i)*nvec(1);
    Fy = Fy + pWall(i)*nvec(2);
end
CD = ( Fx*cos(alpha) + Fy*sin(alpha) )/qInf;
CL = (-Fx*sin(alpha) + Fy*cos(alpha) )/qInf;

fprintf('\n计算完成：\n');
fprintf('  最终残差 = %.6e\n',history.res(end));
fprintf('  迭代步数 = %d\n',history.iter(end));
fprintf('  压力积分估算：CL = %.6f, CD = %.6f\n',CL,CD);

save('naca0012_euler_result.mat','U','rho','u','v','p','Ma','CpCell','CpWall', ...
     'history','free','CL','CD','geom');
fprintf('已保存流场结果：naca0012_euler_result.mat\n');

fig4 = figure('Name','主计算残差曲线');
semilogy(history.iter,history.res,'LineWidth',1.3); grid on;
xlabel('迭代步数'); ylabel('密度相对变化残差');
title('Euler 方程时间推进收敛曲线');
safe_export(fig4,'naca0012_euler_residual.png');

%% ======================== 局部函数 ========================
function geom = build_geometry(X,Y)
    nI = size(X,1);
    nJ = size(X,2)-1;
    xc = zeros(nI,nJ); yc = zeros(nI,nJ); area = zeros(nI,nJ);

    for j = 1:nJ
        for i = 1:nI
            ip = i+1; if ip>nI, ip=1; end
            xv = [X(i,j), X(ip,j), X(ip,j+1), X(i,j+1)];
            yv = [Y(i,j), Y(ip,j), Y(ip,j+1), Y(i,j+1)];
            xc(i,j) = mean(xv); yc(i,j) = mean(yv);
            area(i,j) = abs(0.5*sum(xv.*yv([2:end,1])-yv.*xv([2:end,1])));
        end
    end

    nXi = zeros(nI,nJ,2);       % 周向内部面，方向：cell(i,j) -> cell(i+1,j)
    lenXi = zeros(nI,nJ);
    for j = 1:nJ
        for i = 1:nI
            ip = i+1; if ip>nI, ip=1; end
            p1 = [X(ip,j),Y(ip,j)];
            p2 = [X(ip,j+1),Y(ip,j+1)];
            e = p2-p1;
            n = [e(2),-e(1)];
            toR = [xc(ip,j)-xc(i,j), yc(ip,j)-yc(i,j)];
            if dot(n,toR) < 0, n = -n; end
            nXi(i,j,:) = n;
            lenXi(i,j) = norm(n);
        end
    end

    nEta = zeros(nI,nJ-1,2);    % 法向内部面，方向：cell(i,j) -> cell(i,j+1)
    lenEta = zeros(nI,nJ-1);
    for j = 1:nJ-1
        for i = 1:nI
            ip = i+1; if ip>nI, ip=1; end
            p1 = [X(i,j+1),Y(i,j+1)];
            p2 = [X(ip,j+1),Y(ip,j+1)];
            e = p2-p1;
            n = [e(2),-e(1)];
            toR = [xc(i,j+1)-xc(i,j), yc(i,j+1)-yc(i,j)];
            if dot(n,toR) < 0, n = -n; end
            nEta(i,j,:) = n;
            lenEta(i,j) = norm(n);
        end
    end

    nWall = zeros(nI,2); lenWall = zeros(nI,1);
    for i = 1:nI
        ip = i+1; if ip>nI, ip=1; end
        p1 = [X(i,1),Y(i,1)];
        p2 = [X(ip,1),Y(ip,1)];
        mid = 0.5*(p1+p2);
        e = p2-p1;
        n = [e(2),-e(1)];
        toFace = mid - [xc(i,1),yc(i,1)];
        if dot(n,toFace) < 0, n = -n; end
        nWall(i,:) = n;
        lenWall(i) = norm(n);
    end

    nFar = zeros(nI,2); lenFar = zeros(nI,1);
    for i = 1:nI
        ip = i+1; if ip>nI, ip=1; end
        p1 = [X(i,end),Y(i,end)];
        p2 = [X(ip,end),Y(ip,end)];
        mid = 0.5*(p1+p2);
        e = p2-p1;
        n = [e(2),-e(1)];
        toFace = mid - [xc(i,end),yc(i,end)];
        if dot(n,toFace) < 0, n = -n; end
        nFar(i,:) = n;
        lenFar(i) = norm(n);
    end

    geom = struct('X',X,'Y',Y,'nI',nI,'nJ',nJ,'xc',xc,'yc',yc,'area',area, ...
        'nXi',nXi,'lenXi',lenXi,'nEta',nEta,'lenEta',lenEta, ...
        'nWall',nWall,'lenWall',lenWall,'nFar',nFar,'lenFar',lenFar);
end

function [U,history] = advance_euler(U0,geom,free,gamma,CFL,maxIter,tol,useLimiter,timeMode,plotEvery,printFlag)
    U = U0;
    history.iter = zeros(maxIter,1);
    history.res  = zeros(maxIter,1);

    for iter = 1:maxIter
        Uold = U;

        [L0,sigma0] = compute_rhs(U,geom,free,gamma,useLimiter);
        dt = compute_dt(CFL,geom.area,sigma0,timeMode);
        U1 = add_update(U,L0,dt);
        U1 = repair_nonphysical(U1,free,gamma);

        [L1,~] = compute_rhs(U1,geom,free,gamma,useLimiter);
        Ua = add_update(U1,L1,dt);
        U2 = 0.75*U + 0.25*Ua;
        U2 = repair_nonphysical(U2,free,gamma);

        [L2,~] = compute_rhs(U2,geom,free,gamma,useLimiter);
        Ub = add_update(U2,L2,dt);
        U  = (1.0/3.0)*U + (2.0/3.0)*Ub;
        U = repair_nonphysical(U,free,gamma);

        drho = U(:,:,1)-Uold(:,:,1);
        res = sqrt(mean(drho(:).^2))/(sqrt(mean(Uold(:,:,1).^2,'all'))+eps);
        history.iter(iter) = iter;
        history.res(iter)  = res;

        if printFlag && (mod(iter,plotEvery)==0 || iter==1)
            fprintf('  iter=%6d, residual=%.3e\n',iter,res);
        end
        if res < tol
            history.iter = history.iter(1:iter);
            history.res  = history.res(1:iter);
            if printFlag
                fprintf('达到收敛判据：iter=%d, residual=%.3e\n',iter,res);
            end
            return;
        end
    end
    history.iter = history.iter(1:maxIter);
    history.res  = history.res(1:maxIter);
end

function [L,sigma] = compute_rhs(U,geom,free,gamma,useLimiter)
    nI = geom.nI; nJ = geom.nJ;
    L = zeros(size(U));
    sigma = zeros(nI,nJ);

    % 1) 周向内部面，周期边界
    for j = 1:nJ
        for i = 1:nI
            ip = i+1; if ip>nI, ip=1; end
            im = i-1; if im<1, im=nI; end
            ipp = ip+1; if ipp>nI, ipp=1; end

            UL0 = squeeze(U(i,j,:))';
            UR0 = squeeze(U(ip,j,:))';
            UL = UL0; UR = UR0;

            if useLimiter
                slopeL = minmod_vec(UL0-squeeze(U(im,j,:))',UR0-UL0);
                slopeR = minmod_vec(UR0-UL0,squeeze(U(ipp,j,:))'-UR0);
                ULc = UL0 + 0.5*slopeL;
                URc = UR0 - 0.5*slopeR;
                if is_physical(ULc,gamma) && is_physical(URc,gamma)
                    UL = ULc; UR = URc;
                end
            end

            nvec = squeeze(geom.nXi(i,j,:))';
            [flux,smax] = rusanov_flux(UL,UR,nvec,gamma);
            for k = 1:4
                L(i,j,k)  = L(i,j,k)  - flux(k)/geom.area(i,j);
                L(ip,j,k) = L(ip,j,k) + flux(k)/geom.area(ip,j);
            end
            sigma(i,j)  = sigma(i,j)  + smax;
            sigma(ip,j) = sigma(ip,j) + smax;
        end
    end

    % 2) 法向内部面
    for j = 1:nJ-1
        for i = 1:nI
            jm = max(j-1,1);
            jp = j+1;
            jpp = min(j+2,nJ);

            UL0 = squeeze(U(i,j,:))';
            UR0 = squeeze(U(i,j+1,:))';
            UL = UL0; UR = UR0;

            if useLimiter
                slopeL = minmod_vec(UL0-squeeze(U(i,jm,:))',UR0-UL0);
                slopeR = minmod_vec(UR0-UL0,squeeze(U(i,jpp,:))'-UR0);
                ULc = UL0 + 0.5*slopeL;
                URc = UR0 - 0.5*slopeR;
                if is_physical(ULc,gamma) && is_physical(URc,gamma)
                    UL = ULc; UR = URc;
                end
            end

            nvec = squeeze(geom.nEta(i,j,:))';
            [flux,smax] = rusanov_flux(UL,UR,nvec,gamma);
            for k = 1:4
                L(i,j,k)   = L(i,j,k)   - flux(k)/geom.area(i,j);
                L(i,j+1,k) = L(i,j+1,k) + flux(k)/geom.area(i,j+1);
            end
            sigma(i,j)   = sigma(i,j)   + smax;
            sigma(i,j+1) = sigma(i,j+1) + smax;
        end
    end

    % 3) 壁面滑移边界：无穿透，质量与能量通量为零，仅保留压力通量
    for i = 1:nI
        Ui = squeeze(U(i,1,:))';
        [~,~,~,p] = primitive_from_U(Ui,gamma);
        nvec = geom.nWall(i,:);
        flux = [0, p*nvec(1), p*nvec(2), 0];
        for k = 1:4
            L(i,1,k) = L(i,1,k) - flux(k)/geom.area(i,1);
        end
        [~,~,~,~,c] = primitive_from_U(Ui,gamma);
        sigma(i,1) = sigma(i,1) + c*norm(nvec);
    end

    % 4) 远场边界：采用自由来流外状态的 Rusanov 通量
    for i = 1:nI
        UL = squeeze(U(i,nJ,:))';
        UR = free.U;
        nvec = geom.nFar(i,:);
        [flux,smax] = rusanov_flux(UL,UR,nvec,gamma);
        for k = 1:4
            L(i,nJ,k) = L(i,nJ,k) - flux(k)/geom.area(i,nJ);
        end
        sigma(i,nJ) = sigma(i,nJ) + smax;
    end
end

function dt = compute_dt(CFL,area,sigma,timeMode)
    localDt = CFL*area./max(sigma,1.0e-14);
    if strcmpi(timeMode,'global')
        dt = min(localDt(:));
    else
        dt = localDt;
    end
end

function Unew = add_update(U,L,dt)
    Unew = U;
    for k = 1:4
        Unew(:,:,k) = U(:,:,k) + dt.*L(:,:,k);
    end
end

function [flux,smax] = rusanov_flux(UL,UR,nvec,gamma)
    [FL,unL,cL] = normal_flux(UL,nvec,gamma);
    [FR,unR,cR] = normal_flux(UR,nvec,gamma);
    len = norm(nvec) + eps;
    smax = max(abs(unL)/len + cL, abs(unR)/len + cR)*len;
    flux = 0.5*(FL+FR) - 0.5*smax*(UR-UL);
end

function [F,un,c] = normal_flux(U,nvec,gamma)
    [rho,u,v,p,c] = primitive_from_U(U,gamma);
    E = U(4);
    un = u*nvec(1) + v*nvec(2);
    F = [rho*un, rho*u*un+p*nvec(1), rho*v*un+p*nvec(2), (E+p)*un];
end

function [rho,u,v,p,c] = primitive_from_U(U,gamma)
    rho = max(U(1),1.0e-12);
    u = U(2)/rho;
    v = U(3)/rho;
    kinetic = 0.5*rho*(u^2+v^2);
    p = (gamma-1)*(U(4)-kinetic);
    p = max(p,1.0e-12);
    c = sqrt(gamma*p/rho);
end

function [rho,u,v,p,Ma] = primitive_from_U_field(U,gamma)
    rho = U(:,:,1);
    u   = U(:,:,2)./rho;
    v   = U(:,:,3)./rho;
    p   = (gamma-1)*(U(:,:,4)-0.5*rho.*(u.^2+v.^2));
    p   = max(p,1.0e-12);
    Ma  = sqrt(u.^2+v.^2)./sqrt(gamma*p./rho);
end

function s = minmod_vec(a,b)
    s = zeros(size(a));
    idx = (a.*b) > 0;
    s(idx) = sign(a(idx)).*min(abs(a(idx)),abs(b(idx)));
end

function tf = is_physical(U,gamma)
    if any(~isfinite(U)) || U(1) <= 1.0e-10
        tf = false; return;
    end
    rho = U(1); u = U(2)/rho; v = U(3)/rho;
    p = (gamma-1)*(U(4)-0.5*rho*(u^2+v^2));
    tf = p > 1.0e-10;
end

function U = repair_nonphysical(U,free,gamma)
    [nI,nJ,~] = size(U);
    for j = 1:nJ
        for i = 1:nI
            q = squeeze(U(i,j,:))';
            if ~is_physical(q,gamma)
                for k = 1:4
                    U(i,j,k) = free.U(k);
                end
            end
        end
    end
end

function safe_export(figHandle,fileName)
    try
        exportgraphics(figHandle,fileName,'Resolution',300);
    catch
        saveas(figHandle,fileName);
    end
end
