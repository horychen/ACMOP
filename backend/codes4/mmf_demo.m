%% 计算按照槽口分布的电流密度
%%%%%%%%%%% 输入参数 %%%%%%%%
Qs = 12;                 % 定子槽数
p = 4;                   % 极对数（这里没用到）
m = 3;                   % 相数，后面计算绕组系数会用到
alpha_u = 2*pi/Qs;       % 槽间距角 [rad]
N = 30000;               % 采样步数
alpha = linspace(0,2*pi,N);  % 机械角度（共N个采样点）
alpha_ref = 0:alpha_u:2*pi;  % 每个槽的参考角，注意最后一个值 2pi 与 0 重合

% 三相电流幅值
[i_u, i_v, i_w] = deal(1,-0.5,-0.5);

% U相上层导体的电流密度（标幺）及每槽导体数，这里只考虑上层
I_U = [-1  1  0  -1  1  0  -1  1  0  -1  1  0  ];  
ZQ_U = ones(1, Qs);   % 每槽导体数
A_U = I_U .* ZQ_U;    % U相上层导体电流密度



%%%%%%%%%%% 参数：齿宽与槽开口 %%%%%%%%%%
b_t =0.1745;
b_0 =    0.3491;
tooth_width = b_t;    % 齿宽 10°
slot_open   = b_0;    % 槽开口 20°
% 注：对于Qs=12时，槽间距alpha_u为30°，满足 tooth_width+slot_open = 30°。

%%%%%%%%%%% 计算 U 相 MMF（同时考虑上层与下层） %%%%%%%%%%

% 初始化两个矩阵：
MMF_U_upper = zeros(N, Qs);   % 上层导体的 MMF 贡献
MMF_U_lower = zeros(N, Qs);   % 下层导体的 MMF 贡献（返向电流）

for i = 1:Qs
    %% 上层导体贡献
    % 对于槽 i，上层导体所在槽的导通区间定义为：
    % [ theta_start_upper, theta_start_upper + slot_open )
    % 其中 theta_start_upper = mod(alpha_ref(i) + tooth_width, 2*pi)
    theta_start_upper = mod(alpha_ref(i) + tooth_width, 2*pi);
    theta_end_upper   = mod(theta_start_upper + slot_open, 2*pi);

    for j = 1:N
        theta = mod(alpha(j), 2*pi);
        if theta_start_upper < theta_end_upper
            % 不跨越2pi的情况
            if (theta >= theta_start_upper) && (theta < theta_end_upper)
                MMF_U_upper(j,i) = A_U(i);
            else
                MMF_U_upper(j,i) = 0;
            end
        else
            % 如果跨越2pi（例如从350°到10°）
            if (theta >= theta_start_upper) || (theta < theta_end_upper)
                MMF_U_upper(j,i) = A_U(i);
            else
                MMF_U_upper(j,i) = 0;
            end
        end
    end

    %% 下层导体贡献
    % 对于下层导体，假定其位置在与上层对应的下一槽内，即槽 index: i_lower = mod(i, Qs) + 1
    % 下层导体电流方向与上层相反，所以贡献为 -A_U(i)。
    lower_slot = mod(i, Qs) + 1;  
    % 利用该槽的参考角计算下层导体的导通区间
    theta_start_lower = mod(alpha_ref(lower_slot) + tooth_width, 2*pi);
    theta_end_lower   = mod(theta_start_lower + slot_open, 2*pi);

    for j = 1:N
        theta = mod(alpha(j), 2*pi);
        if theta_start_lower < theta_end_lower
            if (theta >= theta_start_lower) && (theta < theta_end_lower)
                MMF_U_lower(j,i) = -A_U(i);  % 负号代表返向电流
            else
                MMF_U_lower(j,i) = 0;
            end
        else
            if (theta >= theta_start_lower) || (theta < theta_end_lower)
                MMF_U_lower(j,i) = -A_U(i);
            else
                MMF_U_lower(j,i) = 0;
            end
        end
    end
end

% 总的 U 相 MMF 为上层和下层的贡献之和（各槽独立叠加）
MMF_U_pha = sum(MMF_U_upper + MMF_U_lower, 2);
% 注意：此时 MMF_U_pha 已经包含了成对线圈上侧和下侧的作用。

%%%%%%%%%%% 将U相平移得到V、W相MMF %%%%%%%%%%%
Delta_N = round(N/3);   % 移相点数
MMF_V_pha = circshift(MMF_U_pha, Delta_N);
MMF_W_pha = circshift(MMF_V_pha, Delta_N);

% 计算三相总MMF（各相乘以对应电流）
MMF_tot = MMF_U_pha*i_u + MMF_V_pha*i_v + MMF_W_pha*i_w;

% 归一化到 [-1,1]，利用 mapminmax（注意转置）
% [MMF_tot_norm, PS] = mapminmax(MMF_tot');
MMF_tot_norm = MMF_tot;

%%%%%%%%%%% 频谱分析 %%%%%%%%%%%
N_fft = length(MMF_tot_norm);
y_fft = fft(MMF_tot_norm);
P2 = abs(y_fft/N_fft);
P1 = P2(1:N_fft/2+1);
P1(2:end-1) = 2*P1(2:end-1);
n_harm = length(P1);
Order = 0:(n_harm-1);

%%%%%%%%%%% 绘图 %%%%%%%%%%%
figure('Position',[400, 130, 700, 480]);
subplot(2,1,1);
plot(alpha, MMF_tot_norm, 'r', 'LineWidth', 1);
grid on;
set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridAlpha', 1);
set(gcf, 'color', 'w');
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('Mechanical position [rad]');
ylabel('Current linkage [A]');
xlim([0, 2*pi]);
set(gca, 'XTick', 0:pi/6:2*pi);
ylim([-2, 3]);
set(gca, 'YTick', -2:0.5:3);

subplot(2,1,2);
stem(Order, P1, 'r', 'LineWidth', 1);
grid on;
set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridAlpha', 1);
set(gcf, 'color', 'w');
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('Order');
ylabel('Current linkage [A]');
xlim([0, 30]);
set(gca, 'XTick', 0:1:30);
ylim([0, 2]);
set(gca, 'YTick', 0:0.2:2);
MMF_tot_trans = MMF_tot.'./5; % 没用

