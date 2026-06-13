%% cfd_naca0012_grid.m
% NACA0012 尖尾缘翼型 O 型贴体网格生成
% 方法：先代数插值形成初始 O 网格，再求解椭圆型网格生成方程进行光顺。
% 功能：控制壁面法向间距与网格正交性，输出网格文件与质量统计结果。

clear; clc; close all;

%% 1. 参数设置
set(0,'DefaultFigureColor','w');
try
    set(0,'DefaultAxesFontName','SimSun');
    set(0,'DefaultTextFontName','SimSun');
catch
end
set(0,'DefaultAxesFontSize',10.5);

c      = 1.0;       % 弦长
thick  = 0.12;      % NACA0012 厚度比
nI     = 160;       % 周向网格点数，机翼表面点数大于 80
nJ     = 65;        % 法向网格点数
Rfar   = 40.0*c;    % 远场圆半径，满足大于 25 倍弦长；第二题要求可取 40 倍
xc     = 0.5*c;     % 远场圆心 x 坐标
yc     = 0.0;       % 远场圆心 y 坐标

wallSpacing    = 8.0e-4*c;   % 第一层法向间距，用于控制壁面附近加密程度
stretchStrength = 3.5;       % 法向拉伸系数，越大表示靠近壁面越密
maxIter        = 2500;        % 椭圆光顺最大迭代次数
tol            = 1.0e-9;      % 光顺收敛判据
omega          = 1.45;        % SOR 松弛因子，1.0 为 Gauss-Seidel

% 控制源项参数。sourceStrength=0 时退化为 Laplace 型椭圆光顺；
% 取 0.02~0.08 可进一步加强靠壁面网格吸引，但过大可能导致网格扭曲。
sourceStrength = 0.00;
sourceDecay    = 12.0;

%% 2. NACA0012 尖尾缘内边界与圆形远场边界
[xWall,yWall] = naca0012_sharp_te(nI,c,thick);

% 外边界采用相对翼型中心的圆形远场
ang = atan2(yWall-yc,xWall-xc);
xFar = xc + Rfar*cos(ang);
yFar = yc + Rfar*sin(ang);

%% 3. 构造初始代数网格，并固定第一层法向网格以增强正交性
X = zeros(nI,nJ);
Y = zeros(nI,nJ);
X(:,1)   = xWall(:);
Y(:,1)   = yWall(:);
X(:,end) = xFar(:);
Y(:,end) = yFar(:);

% 壁面外法向方向：由切向量旋转得到，并通过远场方向修正符号
[nxWall,nyWall] = wall_outward_normal(xWall,yWall,xFar,yFar);
X(:,2) = X(:,1) + wallSpacing*nxWall(:);
Y(:,2) = Y(:,1) + wallSpacing*nyWall(:);

% 从第一层法向网格到远场作指数拉伸插值
for j = 3:nJ
    eta = (j-2)/(nJ-2);
    s   = (exp(stretchStrength*eta)-1)/(exp(stretchStrength)-1);
    X(:,j) = (1-s)*X(:,2) + s*X(:,end);
    Y(:,j) = (1-s)*Y(:,2) + s*Y(:,end);
end

%% 4. 求解椭圆型网格生成方程进行光顺
% 采用非线性椭圆方程：
% alpha*r_xixi - 2*beta*r_xieta + gamma*r_etaeta + J^2(P*r_xi+Q*r_eta)=0
% 其中 r=(x,y)。默认 P=0，Q 可在靠壁面区域设置指数型源项以控制法向间距。
fprintf('开始椭圆型贴体网格光顺：nI=%d, nJ=%d, Rfar=%.1f c\n',nI,nJ,Rfar/c);

for iter = 1:maxIter
    maxMove = 0.0;
    Xold = X;
    Yold = Y;

    for j = 3:nJ-1
        qControl = -sourceStrength*exp(-((j-2)/sourceDecay)^2);
        pControl = 0.0;

        for i = 1:nI
            ip = i + 1; if ip > nI, ip = 1; end
            im = i - 1; if im < 1,  im = nI; end

            x_xi  = 0.5*(Xold(ip,j)-Xold(im,j));
            y_xi  = 0.5*(Yold(ip,j)-Yold(im,j));
            x_eta = 0.5*(Xold(i,j+1)-Xold(i,j-1));
            y_eta = 0.5*(Yold(i,j+1)-Yold(i,j-1));

            alpha = x_eta^2 + y_eta^2;
            beta  = x_xi*x_eta + y_xi*y_eta;
            gamma = x_xi^2 + y_xi^2;
            jac   = x_xi*y_eta - x_eta*y_xi;

            crossX = Xold(ip,j+1)-Xold(ip,j-1)-Xold(im,j+1)+Xold(im,j-1);
            crossY = Yold(ip,j+1)-Yold(ip,j-1)-Yold(im,j+1)+Yold(im,j-1);

            srcX = jac^2*(pControl*x_xi + qControl*x_eta);
            srcY = jac^2*(pControl*y_xi + qControl*y_eta);

            denom = 2*(alpha+gamma) + eps;
            xNew = (alpha*(X(ip,j)+X(im,j)) + gamma*(X(i,j+1)+X(i,j-1)) ...
                   -0.5*beta*crossX + srcX)/denom;
            yNew = (alpha*(Y(ip,j)+Y(im,j)) + gamma*(Y(i,j+1)+Y(i,j-1)) ...
                   -0.5*beta*crossY + srcY)/denom;

            xSor = (1-omega)*X(i,j) + omega*xNew;
            ySor = (1-omega)*Y(i,j) + omega*yNew;

            maxMove = max(maxMove,hypot(xSor-X(i,j),ySor-Y(i,j)));
            X(i,j) = xSor;
            Y(i,j) = ySor;
        end
    end

    % 固定三条边界：壁面、第一层法向控制线、远场
    X(:,1)   = xWall(:);       Y(:,1)   = yWall(:);
    X(:,2)   = xWall(:) + wallSpacing*nxWall(:);
    Y(:,2)   = yWall(:) + wallSpacing*nyWall(:);
    X(:,end) = xFar(:);        Y(:,end) = yFar(:);

    if mod(iter,100)==0 || iter==1
        fprintf('  iter=%5d, maxMove=%.3e\n',iter,maxMove);
    end
    if maxMove < tol
        fprintf('椭圆光顺收敛：iter=%d, maxMove=%.3e\n',iter,maxMove);
        break;
    end
end

%% 5. 网格质量统计
[cellArea,minArea,maxArea] = structured_cell_area(X,Y);
normalSpacing = hypot(X(:,2)-X(:,1),Y(:,2)-Y(:,1));
orthDev = wall_orthogonality_deviation(X,Y);

fprintf('\n网格质量统计：\n');
fprintf('  最小单元面积 = %.6e\n',minArea);
fprintf('  最大单元面积 = %.6e\n',maxArea);
fprintf('  第一层法向间距：min=%.6e, mean=%.6e, max=%.6e\n', ...
        min(normalSpacing),mean(normalSpacing),max(normalSpacing));
fprintf('  壁面正交性偏差：mean=%.3f deg, max=%.3f deg\n', ...
        mean(orthDev),max(orthDev));

if minArea <= 0
    warning('检测到非正面积单元，建议降低 sourceStrength 或 omega，或增大 nJ。');
end

%% 6. 保存结果与绘图
params = struct('c',c,'thick',thick,'nI',nI,'nJ',nJ,'Rfar',Rfar, ...
    'wallSpacing',wallSpacing,'stretchStrength',stretchStrength, ...
    'sourceStrength',sourceStrength,'sourceDecay',sourceDecay);
save('naca0012_grid.mat','X','Y','params','cellArea','normalSpacing','orthDev');
fprintf('\n已保存网格文件：naca0012_grid.mat\n');

fig1 = figure('Name','NACA0012 O 型贴体网格');
plot_structured_grid(X,Y,4,4);
axis equal; xlabel('x/c'); ylabel('y/c');
title('NACA0012 O 型贴体网格');
safe_export(fig1,'naca0012_grid_full.png');

fig2 = figure('Name','NACA0012 壁面附近网格');
plot_structured_grid(X,Y,2,1);
axis equal; xlim([-0.15,1.15]); ylim([-0.25,0.25]);
xlabel('x/c'); ylabel('y/c');
title('NACA0012 壁面附近网格');
safe_export(fig2,'naca0012_grid_zoom.png');

fig3 = figure('Name','壁面第一层间距与正交性');
yyaxis left; plot(1:nI,normalSpacing,'LineWidth',1.2); ylabel('第一层法向间距');
yyaxis right; plot(1:nI,orthDev,'LineWidth',1.2); ylabel('正交性偏差/(°)');
xlabel('壁面离散点序号'); grid on;
title('壁面附近网格控制效果');
safe_export(fig3,'naca0012_grid_quality.png');

%% ======================== 局部函数 ========================
function [x,y] = naca0012_sharp_te(nI,c,t)
% NACA 四位数对称翼型尖尾缘厚度公式，尾缘系数为 -0.1036。
    if mod(nI,2) ~= 0
        error('nI 应取偶数，例如 160、200、240。');
    end
    nHalf = nI/2 + 1;
    beta  = linspace(0,pi,nHalf);
    xx    = 0.5*c*(1+cos(beta));       % TE -> LE，前后缘加密
    xxn   = xx/c;
    yt    = 5*t*c*(0.2969*sqrt(xxn) - 0.1260*xxn - 0.3516*xxn.^2 ...
             + 0.2843*xxn.^3 - 0.1036*xxn.^4);

    xUpper = xx;
    yUpper = yt;
    xLower = xx(end-1:-1:2);
    yLower = -yt(end-1:-1:2);

    x = [xUpper, xLower];
    y = [yUpper, yLower];
end

function [nx,ny] = wall_outward_normal(x,y,xFar,yFar)
    n = numel(x);
    nx = zeros(size(x)); ny = zeros(size(y));
    for i = 1:n
        ip = i+1; if ip>n, ip=1; end
        im = i-1; if im<1, im=n; end
        tx = x(ip)-x(im);
        ty = y(ip)-y(im);
        nnx = ty;
        nny = -tx;
        len = hypot(nnx,nny) + eps;
        nnx = nnx/len; nny = nny/len;
        toFarX = xFar(i)-x(i);
        toFarY = yFar(i)-y(i);
        if nnx*toFarX + nny*toFarY < 0
            nnx = -nnx; nny = -nny;
        end
        nx(i) = nnx; ny(i) = nny;
    end
end

function [area,minArea,maxArea] = structured_cell_area(X,Y)
    nI = size(X,1); nJ = size(X,2)-1;
    area = zeros(nI,nJ);
    for j = 1:nJ
        for i = 1:nI
            ip = i+1; if ip>nI, ip=1; end
            xv = [X(i,j), X(ip,j), X(ip,j+1), X(i,j+1)];
            yv = [Y(i,j), Y(ip,j), Y(ip,j+1), Y(i,j+1)];
            area(i,j) = abs(0.5*sum(xv.*yv([2:end,1])-yv.*xv([2:end,1])));
        end
    end
    minArea = min(area(:));
    maxArea = max(area(:));
end

function orthDev = wall_orthogonality_deviation(X,Y)
    nI = size(X,1);
    orthDev = zeros(nI,1);
    for i = 1:nI
        ip = i+1; if ip>nI, ip=1; end
        im = i-1; if im<1, im=nI; end
        tx = X(ip,1)-X(im,1);
        ty = Y(ip,1)-Y(im,1);
        rx = X(i,2)-X(i,1);
        ry = Y(i,2)-Y(i,1);
        cs = (tx*rx+ty*ry)/(hypot(tx,ty)*hypot(rx,ry)+eps);
        cs = max(-1,min(1,cs));
        angle = acosd(abs(cs));
        orthDev(i) = abs(90-angle);
    end
end

function plot_structured_grid(X,Y,skipI,skipJ)
    nI = size(X,1); nJ = size(X,2);
    hold on;
    for i = 1:skipI:nI
        plot(X(i,:),Y(i,:),'k-','LineWidth',0.35);
    end
    for j = 1:skipJ:nJ
        plot([X(:,j);X(1,j)],[Y(:,j);Y(1,j)],'k-','LineWidth',0.35);
    end
    plot([X(:,1);X(1,1)],[Y(:,1);Y(1,1)],'k-','LineWidth',1.2);
    box on; grid off;
end

function safe_export(figHandle,fileName)
    try
        exportgraphics(figHandle,fileName,'Resolution',300);
    catch
        saveas(figHandle,fileName);
    end
end
