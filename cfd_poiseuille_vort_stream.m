%% cfd_poiseuille_vort_stream.m
% 二维平面泊肃叶流动求解，涡量-流函数法
% 方法：涡量输运方程+泊松方程，对流项二阶迎风格式，扩散项中心差分，SOR迭代
% 功能：计算低雷诺数层流流动，对比数值解与解析解，输出流场云图与误差结果
clear; clc; close all;

%% 1. 基础参数与网格设置
Re      = 200;
Uin     = 1.0;
H       = 1.0;
L       = 20.0;        % 长槽道保证下游充分发展
Nx      = 141;
Ny      = 51;
dx = L/(Nx-1);
dy = H/(Ny-1);
x  = linspace(0, L, Nx);
y  = linspace(0, H, Ny);

% CFL对流+扩散双重稳定时间步
Umax_est = 1.5*Uin;
dt_cfl   = min(dx, dy)/Umax_est;
dt_diff  = 0.5*Re*(dx^2*dy^2)/(dx^2+dy^2);
dt       = min(0.4*dt_cfl, 0.4*dt_diff);
fprintf('稳定时间步 dt = %.4e\n', dt);

% 全局迭代参数
max_iter = 100000;
tol_res  = 1e-6;

% SOR泊松迭代参数
omega_sor   = 1.4;
max_sor     = 300;
sor_tol     = 1e-6;

%% 2. 场变量与边界初始化
psi   = zeros(Nx, Ny);
omega = zeros(Nx, Ny);
u     = zeros(Nx, Ny);
v     = zeros(Nx, Ny);

% 流函数固定边界条件
psi(:,1)    = 0;
psi(:,Ny)   = Uin*H;
psi(1,:)    = Uin*y;
psi(Nx,:)   = psi(Nx-1,:);

% 进口涡量=0，出口一阶外推
omega(1,:)  = 0;
omega(Nx,:) = omega(Nx-1,:);

% Thompson公式计算壁面涡量（不含出口）
for i = 1:Nx-1
    omega(i,1)  = -2*(psi(i,2)-psi(i,1))/dy^2;
    omega(i,Ny) = -2*(psi(i,Ny-1)-psi(i,Ny))/dy^2;
end

% 初始速度场赋值
for i = 2:Nx-1
    for j = 2:Ny-1
        u(i,j) = (psi(i,j+1)-psi(i,j-1))/(2*dy);
        v(i,j) = -(psi(i+1,j)-psi(i-1,j))/(2*dx);
    end
end
u(1,:) = Uin;
u(Nx,:) = u(Nx-1,:);  % 出口初值速度外推

%% 3. 主迭代循环（时间推进+SOR求解）
iter_cnt = 0;
res_omg  = 1e9;

while iter_cnt < max_iter && res_omg > tol_res
    iter_cnt = iter_cnt + 1;

    % 涡量输运方程：对流项二阶迎风格式离散
    RHS = zeros(Nx, Ny);
    for i = 2:Nx-1
        for j = 2:Ny-1
            ui = u(i,j); vi = v(i,j);
            % x方向二阶迎风，近边界降级一阶迎风兜底
            if i>=3 && i<=Nx-2
                if ui >= 0
                    dwdx = (3*omega(i,j)-4*omega(i-1,j)+omega(i-2,j))/(2*dx);
                else
                    dwdx = (-3*omega(i,j)+4*omega(i+1,j)-omega(i+2,j))/(2*dx);
                end
            else
                if ui >= 0
                    dwdx = (omega(i,j)-omega(i-1,j))/dx;
                else
                    dwdx = (omega(i+1,j)-omega(i,j))/dx;
                end
            end
            % y方向二阶迎风
            if j>=3 && j<=Ny-2
                if vi >= 0
                    dwdy = (3*omega(i,j)-4*omega(i,j-1)+omega(i,j-2))/(2*dy);
                else
                    dwdy = (-3*omega(i,j)+4*omega(i,j+1)-omega(i,j-2))/(2*dy);
                end
            else
                if vi >= 0
                    dwdy = (omega(i,j)-omega(i,j-1))/dy;
                else
                    dwdy = (omega(i,j+1)-omega(i,j))/dy;
                end
            end

            conv = ui*dwdx + vi*dwdy;
            diff_x = (omega(i+1,j)-2*omega(i,j)+omega(i-1,j))/dx^2;
            diff_y = (omega(i,j+1)-2*omega(i,j)+omega(i,j-1))/dy^2;
            diff = (diff_x+diff_y)/Re;
            RHS(i,j) = -conv + diff;
        end
    end

    % 显式更新涡量场
    omega_new = omega + dt*RHS;
    omega_new(1,:)  = 0;
    omega_new(Nx,:) = omega_new(Nx-1,:); % 出口∂ω/∂x=0一阶外推

    % ========== SOR迭代求解泊松方程 ∇²ψ = -ω ==========
    psi_new = psi;
    rhs_psi = -omega_new;
    sor_err = 1e9;
    for sor = 1:max_sor
        sor_err = 0;
        for i = 2:Nx-1
            for j = 2:Ny-1
                psi_old = psi_new(i,j);
                psi_new(i,j) = (1-omega_sor)*psi_old + omega_sor*( ...
                    (psi_new(i+1,j)+psi_new(i-1,j))/dx^2 + ...
                    (psi_new(i,j+1)+psi_new(i,j-1))/dy^2 - rhs_psi(i,j) ...
                    ) / (2/dx^2+2/dy^2);
            end
        end
        % 统一施加流函数边界
        psi_new(1,:)  = Uin*y;
        psi_new(:,1)  = 0;
        psi_new(:,Ny) = Uin*H;
        psi_new(Nx,:) = psi_new(Nx-1,:);
        if sor_err < sor_tol
            break;
        end
    end

    % 壁面涡量更新：跳过出口，保护出口边界
    for i = 1:Nx-1
        omega_new(i,1)  = -2*(psi_new(i,2)-psi_new(i,1))/dy^2;
        omega_new(i,Ny) = -2*(psi_new(i,Ny-1)-psi_new(i,Ny))/dy^2;
    end

    % 更新全场速度场
    for i = 2:Nx-1
        for j = 2:Ny-1
            u(i,j) = (psi_new(i,j+1)-psi_new(i,j-1))/(2*dy);
            v(i,j) = -(psi_new(i+1,j)-psi_new(i-1,j))/(2*dx);
        end
    end
    u(1,:)=Uin;
    v(:,:)=0;
    u(Nx,:) = u(Nx-1,:);  % 出口速度零梯度外推
    v(Nx,:) = v(Nx-1,:);

    % 迭代变量更新与残差计算
    d_omg = abs(omega_new-omega);
    omega = omega_new;
    psi   = psi_new;
    res_omg = max(d_omg(:));

    % 定时输出迭代信息
    if mod(iter_cnt,500)==0
        fprintf('迭代%6d步, 涡量最大变化量=%.3e, SOR内迭代=%d\n',iter_cnt,res_omg,sor);
    end
end
fprintf('\n二维计算收敛！总迭代步数=%d, 最终残差=%.3e\n',iter_cnt,res_omg);

%% 4. 速度剖面对比与误差分析
u_exact = 6*Uin*(y/H).*(1-y/H);
x_idx = [1, round(Nx/4), round(Nx/2), round(3*Nx/4), Nx];
titles = {'进口x=0','x=L/4','x=L/2','x=3L/4','出口x=L'};

figure('Position',[80,80,1250,520]);
for k=1:length(x_idx)
    ii = x_idx(k);
    subplot(2,3,k);
    h1 = plot(y, u(ii,:), 'bo-','LineWidth',1.5,'MarkerSize',3); hold on;
    h2 = plot(y, u_exact, 'r--','LineWidth',2);
    legend([h1,h2],{'数值解','解析解'},'Location','best');
    xlabel('y/H'); ylabel('u/U_{in}');
    title(titles{k}); grid on;
end
sgtitle(sprintf('Re=%d 二维泊肃叶流动沿程速度剖面演化（二阶迎风格式）',Re));

% 自动保存速度剖面图
saveas(gcf,'poiseuille_velocity_profile.png');

% 出口截面误差定量计算
u_out = u(Nx,:);
rel_err = abs(u_out-u_exact)/Uin;
fprintf('出口最大相对误差 = %.3f%%\n',max(rel_err)*100);
fprintf('出口中心线最大速度数值解=%.4f，理论值=1.5\n',max(u_out));

%% 5. 流函数与涡量云图绘制
figure('Position',[100,220,1500,420]);

% 流函数云图
subplot(1,2,1);
contourf(x, y, psi', 30);
colorbar;
caxis([0, Uin*H]);
title('流函数 \psi 云图');
xlabel('x'); ylabel('y/H');
axis tight;

% 涡量云图
subplot(1,2,2);
contourf(x, y, omega', 30);
colorbar;
caxis([-6, 6]);
title('涡量 \omega 云图');
xlabel('x'); ylabel('y/H');
axis tight;

colormap(jet);

% 自动保存流函数+涡量云图
saveas(gcf,'poiseuille_psi_omega_contour.png');
