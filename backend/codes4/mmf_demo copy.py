"""
计算按照槽口分布的电流密度和MMF（磁动势）
从MATLAB脚本 mmf_demo.m 转换而来
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from typing import Tuple, Optional, Dict, List



# 配置matplotlib中文字体支持
def setup_chinese_font():
        """配置matplotlib以支持中文字体显示"""
        chinese_fonts = [
            'Microsoft YaHei',      # 微软雅黑
            'SimHei',               # 黑体
            'SimSun',               # 宋体
            'KaiTi',                # 楷体
            'FangSong',             # 仿宋
        ]
        
        font_found = False
        for font_name in chinese_fonts:
            try:
                plt.rcParams['font.sans-serif'] = [font_name]
                plt.rcParams['axes.unicode_minus'] = False
                font_found = True
                print(f"Chinese font configured: {font_name}")
                break
            except Exception:
                continue
        
        if not font_found:
            try:
                plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
                plt.rcParams['axes.unicode_minus'] = False
                print("Warning: Chinese font not found, using default font")
            except:
                pass

setup_chinese_font()


def compute_u_phase_mmf_core(
    I_U_upper: np.ndarray,
    I_U_lower: Optional[np.ndarray] = None,
    Qs: int = 12,
    N: int = 30000,
    tooth_width: float = 0.1745,
    slot_open: float = 0.3491
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    核心MMF计算函数，严格按照MATLAB逻辑，不绘图
    
    Args:
        I_U_upper: U相上层导体的激励模式数组（长度为Qs）
        I_U_lower: U相下层导体的激励模式数组（长度为Qs），如果为None则自动计算为上层反向
        Qs: 定子槽数
        N: 采样步数
        tooth_width: 齿宽 [rad]
        slot_open: 槽开口 [rad]
    
    Returns:
        (alpha, alpha_ref, MMF_U_upper, MMF_U_lower, MMF_U_pha)
        - alpha: 机械角度数组 [rad]
        - alpha_ref: 槽参考角数组
        - MMF_U_upper: 上层导体MMF贡献矩阵 (N, Qs)
        - MMF_U_lower: 下层导体MMF贡献矩阵 (N, Qs)
        - MMF_U_pha: U相MMF数组
    """
    # ========== 输入参数 ==========
    alpha_u = 2 * np.pi / Qs  # 槽间距角 [rad]
    alpha = np.linspace(0, 2*np.pi, N)  # 机械角度（共N个采样点）
    # 每个槽的参考角，注意最后一个值 2pi 与 0 重合
    # MATLAB: alpha_ref = 0:alpha_u:2*pi
    alpha_ref = np.linspace(0, 2*np.pi, Qs + 1)
    
    # 处理激励模式
    if I_U_lower is None:
        # 如果未提供下层激励，则自动计算为上层反向（MATLAB逻辑：-A_U(i)）
        I_U_lower = -I_U_upper
    
    # 确保是numpy数组
    I_U_upper = np.array(I_U_upper)
    I_U_lower = np.array(I_U_lower)
    
    # 每槽导体数（MATLAB: ZQ_U = ones(1, Qs)）
    ZQ_U = np.ones(Qs)
    A_U_upper = I_U_upper * ZQ_U  # U相上层导体电流密度
    A_U_lower = I_U_lower * ZQ_U  # U相下层导体电流密度
    
    # ========== 计算 U 相 MMF（同时考虑上层与下层） ==========
    # 初始化两个矩阵：
    MMF_U_upper = np.zeros((N, Qs))   # 上层导体的 MMF 贡献
    MMF_U_lower = np.zeros((N, Qs))   # 下层导体的 MMF 贡献（返向电流）
    
    for i in range(Qs):
        # 上层导体贡献
        # MATLAB: theta_start_upper = mod(alpha_ref(i) + tooth_width, 2*pi)
        theta_start_upper = (alpha_ref[i] + tooth_width) % (2 * np.pi)
        theta_end_upper = (theta_start_upper + slot_open) % (2 * np.pi)
        
        for j in range(N):
            theta = alpha[j] % (2 * np.pi)
            if theta_start_upper < theta_end_upper:
                # 不跨越2pi的情况
                if (theta >= theta_start_upper) and (theta < theta_end_upper):
                    MMF_U_upper[j, i] = A_U_upper[i]
                else:
                    MMF_U_upper[j, i] = 0
            else:
                # 如果跨越2pi（例如从350°到10°）
                if (theta >= theta_start_upper) or (theta < theta_end_upper):
                    MMF_U_upper[j, i] = A_U_upper[i]
                else:
                    MMF_U_upper[j, i] = 0
        
        # 下层导体贡献
        # MATLAB: lower_slot = mod(i, Qs) + 1
        # Python索引从0开始，所以是 (i + 1) % Qs
        lower_slot = (i + 1) % Qs
        theta_start_lower = (alpha_ref[lower_slot] + tooth_width) % (2 * np.pi)
        theta_end_lower = (theta_start_lower + slot_open) % (2 * np.pi)
        
        for j in range(N):
            theta = alpha[j] % (2 * np.pi)
            if theta_start_lower < theta_end_lower:
                if (theta >= theta_start_lower) and (theta < theta_end_lower):
                    # MATLAB: MMF_U_lower(j,i) = -A_U(i)
                    MMF_U_lower[j, i] = -A_U_upper[i]  # 注意：使用-A_U_upper，不是A_U_lower
                else:
                    MMF_U_lower[j, i] = 0
            else:
                if (theta >= theta_start_lower) or (theta < theta_end_lower):
                    MMF_U_lower[j, i] = -A_U_upper[i]
                else:
                    MMF_U_lower[j, i] = 0
    
    # 总的 U 相 MMF 为上层和下层的贡献之和（各槽独立叠加）
    # MATLAB: MMF_U_pha = sum(MMF_U_upper + MMF_U_lower, 2)
    MMF_U_pha = np.sum(MMF_U_upper + MMF_U_lower, axis=1)
    
    return alpha, alpha_ref, MMF_U_upper, MMF_U_lower, MMF_U_pha


def compute_and_plot_u_phase_mmf(
    I_U_upper: np.ndarray,
    I_U_lower: Optional[np.ndarray] = None,
    Qs: int = 12,
    N: int = 30000,
    tooth_width: float = 0.1745,
    slot_open: float = 0.3491,
    plot_title: str = "U相MMF",
    show_plot: bool = True,
    compute_three_phase: bool = False,
    i_u: float = 1.0,
    i_v: float = -0.5,
    i_w: float = -0.5,
    use_hardcoded_I_U: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    计算并绘制U相的MMF（磁动势）
    
    Args:
        I_U_upper: U相上层导体的激励模式数组（长度为Qs），例如 [-1, 0, 0, -1, ...]
        I_U_lower: U相下层导体的激励模式数组（长度为Qs），如果为None则自动计算为上层反向
        Qs: 定子槽数，默认12
        N: 采样步数，默认30000
        tooth_width: 齿宽 [rad]，默认0.1745 (约10°)
        slot_open: 槽开口 [rad]，默认0.3491 (约20°)
        plot_title: 图表标题
        show_plot: 是否显示图形
        compute_three_phase: 是否计算三相合成MMF（需要设置三相电流）
        i_u, i_v, i_w: 三相电流幅值（仅在compute_three_phase=True时使用）
        use_hardcoded_I_U: 是否使用硬编码的I_U模式（MATLAB模式）
    
    Returns:
        (alpha, MMF_U_pha, MMF_tot, Order, P1)
        - alpha: 机械角度数组 [rad]
        - MMF_U_pha: U相MMF数组
        - MMF_tot: 总MMF数组（如果compute_three_phase=True则为三相合成，否则等于MMF_U_pha）
        - Order: 谐波阶次数组
        - P1: 频谱幅值数组
    """
    # 如果使用硬编码的I_U模式（MATLAB模式）
    if use_hardcoded_I_U:
        I_U_upper = np.array([-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0])
        I_U_lower = None  # 将自动计算为-I_U_upper
    
    # 使用核心计算函数
    alpha, alpha_ref, MMF_U_upper, MMF_U_lower, MMF_U_pha = compute_u_phase_mmf_core(
        I_U_upper=I_U_upper,
        I_U_lower=I_U_lower,
        Qs=Qs,
        N=N,
        tooth_width=tooth_width,
        slot_open=slot_open
    )
    
    # ========== 计算总MMF ==========
    if compute_three_phase:
        # 将U相平移得到V、W相MMF
        Delta_N = round(N / 3)   # 移相点数
        MMF_V_pha = np.roll(MMF_U_pha, Delta_N)
        MMF_W_pha = np.roll(MMF_V_pha, Delta_N)
        # 计算三相总MMF（各相乘以对应电流）
        MMF_tot = MMF_U_pha * i_u + MMF_V_pha * i_v + MMF_W_pha * i_w
    else:
        # 只计算U相MMF
        MMF_tot = MMF_U_pha
    
    MMF_tot_norm = MMF_tot
    
    # ========== 频谱分析 ==========
    N_fft = len(MMF_tot_norm)
    y_fft = np.fft.fft(MMF_tot_norm)
    P2 = np.abs(y_fft / N_fft)
    P1 = P2[:N_fft//2 + 1]
    P1[1:-1] = 2 * P1[1:-1]
    n_harm = len(P1)
    Order = np.arange(n_harm)
    
    # ========== 绘图 ==========
    if show_plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        # 子图1：MMF空间分布
        ax1.plot(alpha, MMF_tot_norm, 'r', linewidth=1)
        ax1.grid(True, linestyle=':', color='k', alpha=1)
        ax1.set_xlabel('Mechanical position [rad]', fontsize=14)
        ax1.set_ylabel('Current linkage [A]', fontsize=14)
        ax1.set_xlim([0, 2*np.pi])
        ax1.set_xticks(np.arange(0, 2*np.pi + np.pi/6, np.pi/6))
        # 自动调整y轴范围
        y_max = np.max(MMF_tot_norm)
        y_min = np.min(MMF_tot_norm)
        y_range = y_max - y_min
        ax1.set_ylim([y_min - 0.1*y_range, y_max + 0.1*y_range])
        ax1.tick_params(labelsize=14)
        ax1.set_title(plot_title, fontsize=16, fontweight='bold')
        
        # 子图2：频谱分析
        markerline, stemlines, baseline = ax2.stem(Order, P1, linefmt='r-', markerfmt='ro', basefmt='k-')
        plt.setp(stemlines, 'linewidth', 1)
        plt.setp(markerline, 'markersize', 3)
        ax2.grid(True, linestyle=':', color='k', alpha=1)
        ax2.set_xlabel('Order', fontsize=14)
        ax2.set_ylabel('Current linkage [A]', fontsize=14)
        ax2.set_xlim([0, 30])
        ax2.set_xticks(np.arange(0, 31, 1))
        # 自动调整y轴范围
        P1_max = np.max(P1)
        ax2.set_ylim([0, P1_max * 1.1])
        ax2.tick_params(labelsize=14)
        
        plt.tight_layout()
        plt.show()
        
        print(f"{plot_title}计算完成！")
        print(f"MMF最大值: {np.max(MMF_tot_norm):.4f}")
        print(f"MMF最小值: {np.min(MMF_tot_norm):.4f}")
        print(f"MMF平均值: {np.mean(MMF_tot_norm):.4f}")
    
    return alpha, MMF_U_pha, MMF_tot, Order, P1


def compute_phase_excitation_patterns(
    coil_pitch_y: int,
    phase_winding_info: Dict[str, List[int]],
    Q: int = 12
) -> Tuple[List[int], List[int], List[int], List[int]]:
    """
    根据线圈节距和绕组信息计算相的激励模式
    
    该函数生成四个激励列表，用于后续绘图：
    - 扭矩激励上层
    - 扭矩激励下层
    - 悬挂激励上层
    - 悬挂激励下层
    
    Args:
        coil_pitch_y: 线圈节距
        phase_winding_info: 相的绕组信息字典，包含：
            - 'upper_slot_conductor_list': 上层导体槽号列表（槽号从1开始）
            - 'lower_slot_conductor_list': 下层导体槽号列表（槽号从1开始）
            - 'suspension_slot_conductor_list': 悬挂绕组槽号列表，负号表示反向（槽号从1开始）
        Q: 总槽数
    
    Returns:
        (torque_upper, torque_lower, suspension_upper, suspension_lower)
        四个长度为Q的列表，每个元素对应一个槽的激励值
    """
    # 初始化四个激励列表（长度为Q，初始值为0）
    torque_upper = [0] * Q
    torque_lower = [0] * Q
    suspension_upper = [0] * Q
    suspension_lower = [0] * Q
    
    # 处理扭矩激励 - 上层
    for slot in phase_winding_info['upper_slot_conductor_list']:
        if 1 <= slot <= Q:
            # 槽号从1开始，转换为Python索引（从0开始）
            torque_upper[slot - 1] = -1
    
    # 处理扭矩激励 - 下层
    for slot in phase_winding_info['lower_slot_conductor_list']:
        if 1 <= slot <= Q:
            torque_lower[slot - 1] = 1
    
    # 处理悬挂激励 - 上层
    for slot_value in phase_winding_info['suspension_slot_conductor_list']:
        slot = abs(slot_value)  # 获取槽号（绝对值）
        sign = 1 if slot_value > 0 else -1  # 获取符号
        
        if 1 <= slot <= Q:
            suspension_upper[slot - 1] = sign
    
    # 处理悬挂激励 - 下层
    # 下层导体的位置 = 上层导体位置 + coil_pitch_y（考虑周期性）
    for slot_value in phase_winding_info['suspension_slot_conductor_list']:
        upper_slot = abs(slot_value)  # 上层槽号
        sign = 1 if slot_value > 0 else -1  # 获取符号
        
        if 1 <= upper_slot <= Q:
            # 计算下层槽号：槽号从1开始，所以是 (upper_slot + coil_pitch_y - 1) % Q + 1
            lower_slot = ((upper_slot - 1 + coil_pitch_y) % Q) + 1
            
            # 下层导体的激励与上层相反（因为电流方向相反）
            suspension_lower[lower_slot - 1] = -sign
    
    return torque_upper, torque_lower, suspension_upper, suspension_lower


def compute_and_plot_mmf_diagnostic(
    I_U_upper: Optional[np.ndarray] = None,
    I_U_lower: Optional[np.ndarray] = None,
    phase_excitation_dict: Optional[Dict[str, Dict[str, np.ndarray]]] = None,
    excitation_type: str = 'torque',
    use_hardcoded_I_U: bool = True,
    use_circshift: bool = True,
    Qs: int = 12,
    N: int = 30000,
    tooth_width: float = 0.1745,
    slot_open: float = 0.3491,
    i_u: float = 1.0,
    i_v: float = -0.5,
    i_w: float = -0.5
) -> Dict:
    """
    完整的MMF诊断函数，展示所有中间步骤的计算结果
    
    Args:
        I_U_upper: U相上层导体的激励模式数组（如果use_hardcoded_I_U=False时使用）
        I_U_lower: U相下层导体的激励模式数组
        phase_excitation_dict: 各相激励模式字典（用于不同激励模式方法）
        excitation_type: 激励类型，'torque' 或 'suspension'
        use_hardcoded_I_U: 是否使用硬编码的I_U模式（MATLAB模式）
        use_circshift: 是否使用circshift方法计算B、C相（MATLAB方法）
        Qs: 定子槽数
        N: 采样步数
        tooth_width: 齿宽 [rad]
        slot_open: 槽开口 [rad]
        i_u, i_v, i_w: 三相电流幅值
    
    Returns:
        包含所有中间结果的字典
    """
    # 确定激励模式
    if use_hardcoded_I_U:
        I_U_upper = np.array([-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0])
        I_U_lower = None  # 将自动计算为-I_U_upper
    
    # ========== 步骤1：计算U相MMF ==========
    alpha, alpha_ref, MMF_U_upper, MMF_U_lower, MMF_U_pha = compute_u_phase_mmf_core(
        I_U_upper=I_U_upper,
        I_U_lower=I_U_lower,
        Qs=Qs,
        N=N,
        tooth_width=tooth_width,
        slot_open=slot_open
    )
    
    # 计算上层和下层MMF的叠加结果（用于可视化）
    MMF_U_upper_sum = np.sum(MMF_U_upper, axis=1)
    MMF_U_lower_sum = np.sum(MMF_U_lower, axis=1)
    
    # ========== 步骤2：计算V、W相MMF（两种方法） ==========
    # 方法1：circshift方法（MATLAB方法）
    if use_circshift:
        Delta_N = round(N / 3)  # MATLAB: Delta_N = round(N/3)
        MMF_V_pha_circshift = np.roll(MMF_U_pha, Delta_N)  # MATLAB: circshift
        MMF_W_pha_circshift = np.roll(MMF_V_pha_circshift, Delta_N)
    else:
        MMF_V_pha_circshift = None
        MMF_W_pha_circshift = None
    
    # 方法2：使用不同激励模式计算B、C相
    MMF_B_pha_excitation = None
    MMF_C_pha_excitation = None
    if phase_excitation_dict is not None:
        if 'B' in phase_excitation_dict:
            _, _, _, _, MMF_B_pha_excitation = compute_u_phase_mmf_core(
                I_U_upper=phase_excitation_dict['B']['upper'],
                I_U_lower=phase_excitation_dict['B']['lower'],
                Qs=Qs, N=N, tooth_width=tooth_width, slot_open=slot_open
            )
        if 'C' in phase_excitation_dict:
            _, _, _, _, MMF_C_pha_excitation = compute_u_phase_mmf_core(
                I_U_upper=phase_excitation_dict['C']['upper'],
                I_U_lower=phase_excitation_dict['C']['lower'],
                Qs=Qs, N=N, tooth_width=tooth_width, slot_open=slot_open
            )
    
    # ========== 步骤3：计算三相合成MMF ==========
    # 方法1：circshift方法
    if use_circshift:
        MMF_total_circshift = MMF_U_pha * i_u + MMF_V_pha_circshift * i_v + MMF_W_pha_circshift * i_w
    else:
        MMF_total_circshift = None
    
    # 方法2：不同激励模式方法
    if MMF_B_pha_excitation is not None and MMF_C_pha_excitation is not None:
        MMF_total_excitation = MMF_U_pha * i_u + MMF_B_pha_excitation * i_v + MMF_C_pha_excitation * i_w
    else:
        MMF_total_excitation = None
    
    # ========== 创建诊断图 ==========
    fig = plt.figure(figsize=(20, 24))
    gs = fig.add_gridspec(7, 4, hspace=0.4, wspace=0.3)
    
    # ========== 步骤1：输入参数验证 ==========
    # Subplot 1.1: 激励模式条形图
    ax1_1 = fig.add_subplot(gs[0, 0])
    I_U_upper_plot = I_U_upper if I_U_upper is not None else np.array([-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0])
    I_U_lower_plot = -I_U_upper_plot if I_U_lower is None else I_U_lower
    x_slots = np.arange(1, Qs + 1)
    ax1_1.bar(x_slots - 0.2, I_U_upper_plot, width=0.4, label='上层', alpha=0.7)
    ax1_1.bar(x_slots + 0.2, I_U_lower_plot, width=0.4, label='下层', alpha=0.7)
    ax1_1.set_xlabel('槽号', fontsize=10)
    ax1_1.set_ylabel('激励值', fontsize=10)
    ax1_1.set_title('U相激励模式', fontsize=11, fontweight='bold')
    ax1_1.legend(fontsize=9)
    ax1_1.grid(True, alpha=0.3)
    ax1_1.set_xticks(x_slots)
    
    # Subplot 1.2: 槽参考角分布
    ax1_2 = fig.add_subplot(gs[0, 1])
    ax1_2.plot(alpha_ref[:-1], np.ones(Qs), 'o', markersize=8, label='槽参考角')
    ax1_2.set_xlabel('角度 [rad]', fontsize=10)
    ax1_2.set_ylabel('', fontsize=10)
    ax1_2.set_title('槽参考角分布', fontsize=11, fontweight='bold')
    ax1_2.set_xlim([0, 2*np.pi])
    ax1_2.set_xticks(np.arange(0, 2*np.pi + np.pi/6, np.pi/6))
    ax1_2.grid(True, alpha=0.3)
    ax1_2.legend(fontsize=9)
    
    # ========== 步骤2：上层导体MMF贡献 ==========
    # Subplot 2.1: 上层导体MMF贡献矩阵（热图）
    ax2_1 = fig.add_subplot(gs[1, 0])
    # 为了可视化，对矩阵进行降采样
    downsample_factor = max(1, N // 1000)
    MMF_U_upper_vis = MMF_U_upper[::downsample_factor, :]
    alpha_vis = alpha[::downsample_factor]
    im = ax2_1.imshow(MMF_U_upper_vis.T, aspect='auto', origin='lower', 
                      extent=[0, 2*np.pi, 0, Qs], cmap='RdBu_r', interpolation='nearest')
    ax2_1.set_xlabel('机械角度 [rad]', fontsize=10)
    ax2_1.set_ylabel('槽号', fontsize=10)
    ax2_1.set_title('上层导体MMF贡献矩阵', fontsize=11, fontweight='bold')
    plt.colorbar(im, ax=ax2_1, label='MMF贡献')
    
    # Subplot 2.2: 上层导体MMF贡献叠加
    ax2_2 = fig.add_subplot(gs[1, 1])
    ax2_2.plot(alpha, MMF_U_upper_sum, 'r-', linewidth=1.5, label='上层MMF')
    ax2_2.set_xlabel('机械角度 [rad]', fontsize=10)
    ax2_2.set_ylabel('MMF [A]', fontsize=10)
    ax2_2.set_title('上层导体MMF贡献叠加', fontsize=11, fontweight='bold')
    ax2_2.set_xlim([0, 2*np.pi])
    ax2_2.grid(True, alpha=0.3)
    ax2_2.legend(fontsize=9)
    
    # ========== 步骤3：下层导体MMF贡献 ==========
    # Subplot 3.1: 下层导体MMF贡献矩阵（热图）
    ax3_1 = fig.add_subplot(gs[2, 0])
    MMF_U_lower_vis = MMF_U_lower[::downsample_factor, :]
    im = ax3_1.imshow(MMF_U_lower_vis.T, aspect='auto', origin='lower',
                      extent=[0, 2*np.pi, 0, Qs], cmap='RdBu_r', interpolation='nearest')
    ax3_1.set_xlabel('机械角度 [rad]', fontsize=10)
    ax3_1.set_ylabel('槽号', fontsize=10)
    ax3_1.set_title('下层导体MMF贡献矩阵', fontsize=11, fontweight='bold')
    plt.colorbar(im, ax=ax3_1, label='MMF贡献')
    
    # Subplot 3.2: 下层导体MMF贡献叠加
    ax3_2 = fig.add_subplot(gs[2, 1])
    ax3_2.plot(alpha, MMF_U_lower_sum, 'b-', linewidth=1.5, label='下层MMF')
    ax3_2.set_xlabel('机械角度 [rad]', fontsize=10)
    ax3_2.set_ylabel('MMF [A]', fontsize=10)
    ax3_2.set_title('下层导体MMF贡献叠加', fontsize=11, fontweight='bold')
    ax3_2.set_xlim([0, 2*np.pi])
    ax3_2.grid(True, alpha=0.3)
    ax3_2.legend(fontsize=9)
    
    # ========== 步骤4：U相MMF ==========
    # Subplot 4.1: U相MMF
    ax4_1 = fig.add_subplot(gs[3, 0])
    ax4_1.plot(alpha, MMF_U_pha, 'k-', linewidth=2, label='U相MMF')
    ax4_1.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax4_1.set_xlabel('机械角度 [rad]', fontsize=10)
    ax4_1.set_ylabel('MMF [A]', fontsize=10)
    ax4_1.set_title('U相MMF（上层+下层）', fontsize=11, fontweight='bold')
    ax4_1.set_xlim([0, 2*np.pi])
    ax4_1.grid(True, alpha=0.3)
    ax4_1.legend(fontsize=9)
    
    # Subplot 4.2: U相MMF频谱
    ax4_2 = fig.add_subplot(gs[3, 1])
    N_fft = len(MMF_U_pha)
    y_fft = np.fft.fft(MMF_U_pha)
    P2 = np.abs(y_fft / N_fft)
    P1 = P2[:N_fft//2 + 1]
    P1[1:-1] = 2 * P1[1:-1]
    Order = np.arange(len(P1))
    markerline, stemlines, baseline = ax4_2.stem(Order, P1, linefmt='k-', markerfmt='ko', basefmt='k-')
    plt.setp(stemlines, 'linewidth', 1)
    plt.setp(markerline, 'markersize', 3)
    ax4_2.set_xlabel('Order', fontsize=10)
    ax4_2.set_ylabel('MMF Amplitude [A]', fontsize=10)
    ax4_2.set_title('U相MMF频谱', fontsize=11, fontweight='bold')
    ax4_2.set_xlim([0, 30])
    ax4_2.grid(True, alpha=0.3)
    
    # ========== 步骤5：V、W相MMF对比 ==========
    # Subplot 5.1: circshift方法（MATLAB方法）
    ax5_1 = fig.add_subplot(gs[4, 0])
    if use_circshift:
        ax5_1.plot(alpha, MMF_U_pha, 'r-', linewidth=1.5, label='U相', alpha=0.7)
        ax5_1.plot(alpha, MMF_V_pha_circshift, 'b-', linewidth=1.5, label='V相(circshift)', alpha=0.7)
        ax5_1.plot(alpha, MMF_W_pha_circshift, 'g-', linewidth=1.5, label='W相(circshift)', alpha=0.7)
    ax5_1.set_xlabel('机械角度 [rad]', fontsize=10)
    ax5_1.set_ylabel('MMF [A]', fontsize=10)
    ax5_1.set_title('V、W相MMF (circshift方法)', fontsize=11, fontweight='bold')
    ax5_1.set_xlim([0, 2*np.pi])
    ax5_1.grid(True, alpha=0.3)
    ax5_1.legend(fontsize=9)
    
    # Subplot 5.2: 不同激励模式方法
    ax5_2 = fig.add_subplot(gs[4, 1])
    if MMF_B_pha_excitation is not None and MMF_C_pha_excitation is not None:
        ax5_2.plot(alpha, MMF_U_pha, 'r-', linewidth=1.5, label='A相', alpha=0.7)
        ax5_2.plot(alpha, MMF_B_pha_excitation, 'b-', linewidth=1.5, label='B相(激励)', alpha=0.7)
        ax5_2.plot(alpha, MMF_C_pha_excitation, 'g-', linewidth=1.5, label='C相(激励)', alpha=0.7)
    ax5_2.set_xlabel('机械角度 [rad]', fontsize=10)
    ax5_2.set_ylabel('MMF [A]', fontsize=10)
    ax5_2.set_title('B、C相MMF (不同激励模式)', fontsize=11, fontweight='bold')
    ax5_2.set_xlim([0, 2*np.pi])
    ax5_2.grid(True, alpha=0.3)
    ax5_2.legend(fontsize=9)
    
    # Subplot 5.3: 两种方法的差异
    ax5_3 = fig.add_subplot(gs[4, 2])
    if use_circshift and MMF_B_pha_excitation is not None:
        diff_V = MMF_V_pha_circshift - MMF_B_pha_excitation
        diff_W = MMF_W_pha_circshift - MMF_C_pha_excitation
        ax5_3.plot(alpha, diff_V, 'b-', linewidth=1.5, label='V/B相差异', alpha=0.7)
        ax5_3.plot(alpha, diff_W, 'g-', linewidth=1.5, label='W/C相差异', alpha=0.7)
        ax5_3.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    ax5_3.set_xlabel('机械角度 [rad]', fontsize=10)
    ax5_3.set_ylabel('MMF差异 [A]', fontsize=10)
    ax5_3.set_title('两种方法差异', fontsize=11, fontweight='bold')
    ax5_3.set_xlim([0, 2*np.pi])
    ax5_3.grid(True, alpha=0.3)
    ax5_3.legend(fontsize=9)
    
    # ========== 步骤6：三相合成MMF ==========
    # Subplot 6.1: 各相MMF对比（circshift方法）
    ax6_1 = fig.add_subplot(gs[5, 0])
    if use_circshift:
        ax6_1.plot(alpha, MMF_U_pha * i_u, 'r-', linewidth=1.5, label=f'A相 (i={i_u})', alpha=0.7)
        ax6_1.plot(alpha, MMF_V_pha_circshift * i_v, 'b-', linewidth=1.5, label=f'V相 (i={i_v})', alpha=0.7)
        ax6_1.plot(alpha, MMF_W_pha_circshift * i_w, 'g-', linewidth=1.5, label=f'W相 (i={i_w})', alpha=0.7)
    ax6_1.set_xlabel('机械角度 [rad]', fontsize=10)
    ax6_1.set_ylabel('MMF [A]', fontsize=10)
    ax6_1.set_title('各相MMF (circshift, 乘以电流)', fontsize=11, fontweight='bold')
    ax6_1.set_xlim([0, 2*np.pi])
    ax6_1.grid(True, alpha=0.3)
    ax6_1.legend(fontsize=9)
    
    # Subplot 6.2: 三相合成MMF（circshift方法）
    ax6_2 = fig.add_subplot(gs[5, 1])
    if use_circshift:
        ax6_2.plot(alpha, MMF_total_circshift, 'k-', linewidth=2, label='合成MMF(circshift)')
        ax6_2.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax6_2.set_xlabel('机械角度 [rad]', fontsize=10)
    ax6_2.set_ylabel('MMF [A]', fontsize=10)
    ax6_2.set_title('三相合成MMF (circshift方法)', fontsize=11, fontweight='bold')
    ax6_2.set_xlim([0, 2*np.pi])
    ax6_2.grid(True, alpha=0.3)
    ax6_2.legend(fontsize=9)
    
    # Subplot 6.3: 三相合成MMF（不同激励模式方法）
    ax6_3 = fig.add_subplot(gs[5, 2])
    if MMF_total_excitation is not None:
        ax6_3.plot(alpha, MMF_total_excitation, 'm-', linewidth=2, label='合成MMF(激励)')
        ax6_3.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax6_3.set_xlabel('机械角度 [rad]', fontsize=10)
    ax6_3.set_ylabel('MMF [A]', fontsize=10)
    ax6_3.set_title('三相合成MMF (不同激励模式)', fontsize=11, fontweight='bold')
    ax6_3.set_xlim([0, 2*np.pi])
    ax6_3.grid(True, alpha=0.3)
    ax6_3.legend(fontsize=9)
    
    # Subplot 6.4: 两种合成MMF的差异
    ax6_4 = fig.add_subplot(gs[5, 3])
    if MMF_total_circshift is not None and MMF_total_excitation is not None:
        diff_total = MMF_total_circshift - MMF_total_excitation
        ax6_4.plot(alpha, diff_total, 'r-', linewidth=1.5, label='合成MMF差异')
        ax6_4.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    ax6_4.set_xlabel('机械角度 [rad]', fontsize=10)
    ax6_4.set_ylabel('MMF差异 [A]', fontsize=10)
    ax6_4.set_title('两种合成MMF差异', fontsize=11, fontweight='bold')
    ax6_4.set_xlim([0, 2*np.pi])
    ax6_4.grid(True, alpha=0.3)
    ax6_4.legend(fontsize=9)
    
    # ========== 步骤7：频谱分析 ==========
    # Subplot 7.1: 合成MMF频谱（circshift方法）
    ax7_1 = fig.add_subplot(gs[6, 0])
    if use_circshift:
        N_fft = len(MMF_total_circshift)
        y_fft = np.fft.fft(MMF_total_circshift)
        P2 = np.abs(y_fft / N_fft)
        P1 = P2[:N_fft//2 + 1]
        P1[1:-1] = 2 * P1[1:-1]
        Order = np.arange(len(P1))
        markerline, stemlines, baseline = ax7_1.stem(Order, P1, linefmt='k-', markerfmt='ko', basefmt='k-')
        plt.setp(stemlines, 'linewidth', 1)
        plt.setp(markerline, 'markersize', 3)
        ax7_1.set_xlabel('Order', fontsize=10)
        ax7_1.set_ylabel('MMF Amplitude [A]', fontsize=10)
        ax7_1.set_title('合成MMF频谱 (circshift)', fontsize=11, fontweight='bold')
        ax7_1.set_xlim([0, 30])
        ax7_1.grid(True, alpha=0.3)
    
    # Subplot 7.2: 合成MMF频谱（不同激励模式方法）
    ax7_2 = fig.add_subplot(gs[6, 1])
    if MMF_total_excitation is not None:
        N_fft = len(MMF_total_excitation)
        y_fft = np.fft.fft(MMF_total_excitation)
        P2 = np.abs(y_fft / N_fft)
        P1 = P2[:N_fft//2 + 1]
        P1[1:-1] = 2 * P1[1:-1]
        Order = np.arange(len(P1))
        markerline, stemlines, baseline = ax7_2.stem(Order, P1, linefmt='m-', markerfmt='mo', basefmt='k-')
        plt.setp(stemlines, 'linewidth', 1)
        plt.setp(markerline, 'markersize', 3)
        ax7_2.set_xlabel('Order', fontsize=10)
        ax7_2.set_ylabel('MMF Amplitude [A]', fontsize=10)
        ax7_2.set_title('合成MMF频谱 (不同激励模式)', fontsize=11, fontweight='bold')
        ax7_2.set_xlim([0, 30])
        ax7_2.grid(True, alpha=0.3)
    
    # 添加数值验证信息
    info_text = f"""
    数值验证信息:
    U相MMF: max={np.max(MMF_U_pha):.4f}, min={np.min(MMF_U_pha):.4f}, mean={np.mean(MMF_U_pha):.6f}
    """
    if use_circshift:
        info_text += f"""
    合成MMF(circshift): max={np.max(MMF_total_circshift):.4f}, min={np.min(MMF_total_circshift):.4f}, mean={np.mean(MMF_total_circshift):.6f}
    电流验证: i_u+i_v+i_w = {i_u + i_v + i_w:.6f}
    """
    if MMF_total_excitation is not None:
        info_text += f"""
    合成MMF(激励): max={np.max(MMF_total_excitation):.4f}, min={np.min(MMF_total_excitation):.4f}, mean={np.mean(MMF_total_excitation):.6f}
    """
    
    # 在空白区域显示信息
    ax_info = fig.add_subplot(gs[6, 2:])
    ax_info.axis('off')
    ax_info.text(0.1, 0.5, info_text, fontsize=10, family='monospace', verticalalignment='center')
    
    plt.suptitle(f'{excitation_type.capitalize()} Excitation - 完整MMF诊断', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.show()
    
    # 返回所有中间结果
    return {
        'alpha': alpha,
        'alpha_ref': alpha_ref,
        'MMF_U_upper': MMF_U_upper,
        'MMF_U_lower': MMF_U_lower,
        'MMF_U_pha': MMF_U_pha,
        'MMF_V_pha_circshift': MMF_V_pha_circshift if use_circshift else None,
        'MMF_W_pha_circshift': MMF_W_pha_circshift if use_circshift else None,
        'MMF_B_pha_excitation': MMF_B_pha_excitation,
        'MMF_C_pha_excitation': MMF_C_pha_excitation,
        'MMF_total_circshift': MMF_total_circshift if use_circshift else None,
        'MMF_total_excitation': MMF_total_excitation
    }


def compute_and_plot_three_phase_mmf_comparison(
    phase_excitation_dict: Dict[str, Dict[str, np.ndarray]],
    excitation_type: str = 'suspension',  # 'torque' or 'suspension'
    Qs: int = 12,
    N: int = 30000,
    tooth_width: float = 0.1745,
    slot_open: float = 0.3491,
    show_plot: bool = True
) -> Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]:
    """
    计算并绘制三相MMF对比图（包括U、B、C相和合成MMF）
    
    Args:
        phase_excitation_dict: 字典，格式为：
            {
                'A': {'upper': array, 'lower': array},
                'B': {'upper': array, 'lower': array},
                'C': {'upper': array, 'lower': array}
            }
        excitation_type: 激励类型，'torque' 或 'suspension'
        Qs: 定子槽数
        N: 采样步数
        tooth_width: 齿宽 [rad]
        slot_open: 槽开口 [rad]
        show_plot: 是否显示图形
    
    Returns:
        (alpha, phase_mmf_dict, MMF_total)
        - alpha: 机械角度数组
        - phase_mmf_dict: 各相MMF字典 {'A': array, 'B': array, 'C': array}
        - MMF_total: 三相合成MMF
    """
    # 计算各相MMF
    phase_mmf_dict = {}
    alpha = None
    
    for phase_name in ['A', 'B', 'C']:
        if phase_name in phase_excitation_dict:
            I_upper = phase_excitation_dict[phase_name]['upper']
            I_lower = phase_excitation_dict[phase_name]['lower']
            
            # 计算该相MMF（不绘图）
            alpha, MMF_U_pha, _, _, _ = compute_and_plot_u_phase_mmf(
                I_U_upper=I_upper,
                I_U_lower=I_lower,
                Qs=Qs,
                N=N,
                tooth_width=tooth_width,
                slot_open=slot_open,
                plot_title="",
                show_plot=False,
                compute_three_phase=False
            )
            phase_mmf_dict[phase_name] = MMF_U_pha
    
    # 计算三相电流（满足 ia+ib+ic=0）
    # 使用对称三相电流：ia = I*cos(wt), ib = I*cos(wt-120°), ic = I*cos(wt-240°)
    # 在某个时刻，例如 wt=0: ia=1, ib=-0.5, ic=-0.5 (满足 ia+ib+ic=0)
    i_a = 1.0
    i_b = -0.5
    i_c = -0.5
    
    # 计算三相合成MMF
    # 注意：B相和C相的空间相位差已经通过不同的槽位配置体现在各自的MMF中了
    # 因此不需要额外的空间偏移
    MMF_A = phase_mmf_dict['A']
    MMF_B = phase_mmf_dict['B']
    MMF_C = phase_mmf_dict['C']
    
    # 三相合成MMF（各相MMF乘以对应电流）
    MMF_total = MMF_A * i_a + MMF_B * i_b + MMF_C * i_c
    
    # 绘图
    if show_plot:
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        
        # 子图1：各相MMF对比
        ax1 = axes[0, 0]
        ax1.plot(alpha, phase_mmf_dict['A'], 'r-', linewidth=1.5, label='A相')
        ax1.plot(alpha, phase_mmf_dict['B'], 'b-', linewidth=1.5, label='B相')
        ax1.plot(alpha, phase_mmf_dict['C'], 'g-', linewidth=1.5, label='C相')
        ax1.grid(True, linestyle=':', color='k', alpha=0.5)
        ax1.set_xlabel('Mechanical position [rad]', fontsize=12)
        ax1.set_ylabel('Current linkage [A]', fontsize=12)
        ax1.set_xlim([0, 2*np.pi])
        ax1.set_xticks(np.arange(0, 2*np.pi + np.pi/6, np.pi/6))
        ax1.legend(fontsize=11)
        ax1.set_title(f'{excitation_type.capitalize()} Excitation - 各相MMF', fontsize=14, fontweight='bold')
        ax1.tick_params(labelsize=11)
        
        # 子图2：三相合成MMF
        ax2 = axes[0, 1]
        ax2.plot(alpha, MMF_total, 'k-', linewidth=2, label='合成MMF')
        ax2.grid(True, linestyle=':', color='k', alpha=0.5)
        ax2.set_xlabel('Mechanical position [rad]', fontsize=12)
        ax2.set_ylabel('Current linkage [A]', fontsize=12)
        ax2.set_xlim([0, 2*np.pi])
        ax2.set_xticks(np.arange(0, 2*np.pi + np.pi/6, np.pi/6))
        ax2.legend(fontsize=11)
        ax2.set_title(f'{excitation_type.capitalize()} Excitation - 三相合成MMF', fontsize=14, fontweight='bold')
        ax2.tick_params(labelsize=11)
        # 添加零线
        ax2.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        
        # 子图3：合成MMF频谱
        ax3 = axes[1, 0]
        N_fft = len(MMF_total)
        y_fft = np.fft.fft(MMF_total)
        P2 = np.abs(y_fft / N_fft)
        P1 = P2[:N_fft//2 + 1]
        P1[1:-1] = 2 * P1[1:-1]
        Order = np.arange(len(P1))
        
        markerline, stemlines, baseline = ax3.stem(Order, P1, linefmt='k-', markerfmt='ko', basefmt='k-')
        plt.setp(stemlines, 'linewidth', 1)
        plt.setp(markerline, 'markersize', 3)
        ax3.grid(True, linestyle=':', color='k', alpha=0.5)
        ax3.set_xlabel('Order', fontsize=12)
        ax3.set_ylabel('Current linkage [A]', fontsize=12)
        ax3.set_xlim([0, 30])
        ax3.set_xticks(np.arange(0, 31, 5))
        P1_max = np.max(P1)
        ax3.set_ylim([0, P1_max * 1.1])
        ax3.set_title(f'{excitation_type.capitalize()} Excitation - 合成MMF频谱', fontsize=14, fontweight='bold')
        ax3.tick_params(labelsize=11)
        
        # 子图4：各相MMF频谱对比
        ax4 = axes[1, 1]
        colors = ['r', 'b', 'g']
        phase_names = ['A', 'B', 'C']
        for idx, (phase_name, color) in enumerate(zip(phase_names, colors)):
            if phase_name in phase_mmf_dict:
                mmf_phase = phase_mmf_dict[phase_name]
                
                N_fft = len(mmf_phase)
                y_fft = np.fft.fft(mmf_phase)
                P2 = np.abs(y_fft / N_fft)
                P1_phase = P2[:N_fft//2 + 1]
                P1_phase[1:-1] = 2 * P1_phase[1:-1]
                Order = np.arange(len(P1_phase))
                
                ax4.plot(Order, P1_phase, color=color, linewidth=1.5, marker='o', markersize=3, label=f'{phase_name}相')
        
        ax4.grid(True, linestyle=':', color='k', alpha=0.5)
        ax4.set_xlabel('Order', fontsize=12)
        ax4.set_ylabel('Current linkage [A]', fontsize=12)
        ax4.set_xlim([0, 30])
        ax4.set_xticks(np.arange(0, 31, 5))
        ax4.legend(fontsize=11)
        ax4.set_title(f'{excitation_type.capitalize()} Excitation - 各相MMF频谱', fontsize=14, fontweight='bold')
        ax4.tick_params(labelsize=11)
        
        plt.tight_layout()
        plt.show()
        
        print(f"\n{excitation_type.capitalize()} Excitation MMF计算完成！")
        print(f"合成MMF最大值: {np.max(MMF_total):.4f}")
        print(f"合成MMF最小值: {np.min(MMF_total):.4f}")
        print(f"合成MMF平均值: {np.mean(MMF_total):.6f} (应该接近0)")
        print(f"电流验证: ia+ib+ic = {i_a + i_b + i_c:.6f}")
    
    return alpha, phase_mmf_dict, MMF_total


def animate_three_phase_mmf_time_evolution(
    phase_excitation_dict: Dict[str, Dict[str, np.ndarray]],
    excitation_type: str = 'suspension',
    Qs: int = 12,
    N: int = 30000,
    tooth_width: float = 0.1745,
    slot_open: float = 0.3491,
    num_frames: int = 100,
    current_amplitude: float = 1.0,
    animation_interval: int = 50
):
    """
    创建动画，展示三相对称交流电下各相MMF和合成MMF随时间的变化
    
    Args:
        phase_excitation_dict: 各相激励模式字典
        excitation_type: 激励类型
        Qs: 定子槽数
        N: 采样步数
        tooth_width: 齿宽 [rad]
        slot_open: 槽开口 [rad]
        num_frames: 动画帧数
        current_amplitude: 电流幅值
        animation_interval: 动画帧间隔（毫秒）
    """
    # 计算各相的单位MMF（单位电流下的MMF）
    phase_unit_mmf_dict = {}
    alpha = None
    
    for phase_name in ['A', 'B', 'C']:
        if phase_name in phase_excitation_dict:
            I_upper = phase_excitation_dict[phase_name]['upper']
            I_lower = phase_excitation_dict[phase_name]['lower']
            
            # 计算该相单位MMF（不绘图）
            alpha, MMF_U_pha, _, _, _ = compute_and_plot_u_phase_mmf(
                I_U_upper=I_upper,
                I_U_lower=I_lower,
                Qs=Qs,
                N=N,
                tooth_width=tooth_width,
                slot_open=slot_open,
                plot_title="",
                show_plot=False,
                compute_three_phase=False
            )
            phase_unit_mmf_dict[phase_name] = MMF_U_pha
    
    # 时间数组（一个电周期）
    t_array = np.linspace(0, 2*np.pi, num_frames)
    
    # 创建图形和子图
    fig = plt.figure(figsize=(16, 10))
    
    # 创建子图布局
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # 子图1：三相电流
    ax1 = fig.add_subplot(gs[0, 0])
    line_ia, = ax1.plot([], [], 'r-', linewidth=2, label='ia')
    line_ib, = ax1.plot([], [], 'b-', linewidth=2, label='ib')
    line_ic, = ax1.plot([], [], 'g-', linewidth=2, label='ic')
    ax1.set_xlim([0, 2*np.pi])
    ax1.set_ylim([-1.2*current_amplitude, 1.2*current_amplitude])
    ax1.set_xlabel('Time (electrical angle [rad])', fontsize=11)
    ax1.set_ylabel('Current [A]', fontsize=11)
    ax1.set_title('三相电流', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # 子图2：各相MMF（随时间变化）
    ax2 = fig.add_subplot(gs[0, 1])
    line_mmf_a, = ax2.plot([], [], 'r-', linewidth=1.5, label='A相MMF', alpha=0.8)
    line_mmf_b, = ax2.plot([], [], 'b-', linewidth=1.5, label='B相MMF', alpha=0.8)
    line_mmf_c, = ax2.plot([], [], 'g-', linewidth=1.5, label='C相MMF', alpha=0.8)
    ax2.set_xlim([0, 2*np.pi])
    y_max = max([np.max(np.abs(phase_unit_mmf_dict[ph])) for ph in ['A', 'B', 'C']])
    ax2.set_ylim([-1.2*y_max*current_amplitude, 1.2*y_max*current_amplitude])
    ax2.set_xlabel('Mechanical position [rad]', fontsize=11)
    ax2.set_ylabel('MMF [A]', fontsize=11)
    ax2.set_title('各相MMF（随时间变化）', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # 子图3：合成MMF（随时间变化）
    ax3 = fig.add_subplot(gs[1, :])
    line_mmf_total, = ax3.plot([], [], 'k-', linewidth=2, label='合成MMF')
    ax3.set_xlim([0, 2*np.pi])
    ax3.set_ylim([-1.5*y_max*current_amplitude, 1.5*y_max*current_amplitude])
    ax3.set_xlabel('Mechanical position [rad]', fontsize=11)
    ax3.set_ylabel('MMF [A]', fontsize=11)
    ax3.set_title('三相合成MMF（随时间变化）', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=10)
    ax3.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # 子图4：合成MMF频谱（随时间变化）
    ax4 = fig.add_subplot(gs[2, :])
    ax4.set_xlim([0, 30])
    ax4.set_ylim([0, 1.5*y_max*current_amplitude])
    ax4.set_xlabel('Order', fontsize=11)
    ax4.set_ylabel('MMF Amplitude [A]', fontsize=11)
    ax4.set_title('合成MMF频谱（随时间变化）', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # 添加时间显示文本
    time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12, fontweight='bold')
    
    # 预计算所有帧的数据
    mmf_a_time = []
    mmf_b_time = []
    mmf_c_time = []
    mmf_total_time = []
    spectrum_time = []
    
    for t in t_array:
        # 计算三相电流
        ia = current_amplitude * np.cos(t)
        ib = current_amplitude * np.cos(t - 2*np.pi/3)
        ic = current_amplitude * np.cos(t - 4*np.pi/3)
        
        # 计算各相MMF（单位MMF × 电流）
        mmf_a = phase_unit_mmf_dict['A'] * ia
        mmf_b = phase_unit_mmf_dict['B'] * ib
        mmf_c = phase_unit_mmf_dict['C'] * ic
        
        # 合成MMF
        mmf_total = mmf_a + mmf_b + mmf_c
        
        # 计算频谱
        N_fft = len(mmf_total)
        y_fft = np.fft.fft(mmf_total)
        P2 = np.abs(y_fft / N_fft)
        P1 = P2[:N_fft//2 + 1]
        P1[1:-1] = 2 * P1[1:-1]
        Order = np.arange(len(P1))
        
        mmf_a_time.append(mmf_a)
        mmf_b_time.append(mmf_b)
        mmf_c_time.append(mmf_c)
        mmf_total_time.append(mmf_total)
        spectrum_time.append((Order, P1))
    
    # 动画更新函数
    def animate(frame):
        t = t_array[frame]
        
        # 计算三相电流
        ia = current_amplitude * np.cos(t)
        ib = current_amplitude * np.cos(t - 2*np.pi/3)
        ic = current_amplitude * np.cos(t - 4*np.pi/3)
        
        # 更新电流图
        line_ia.set_data(t_array[:frame+1], [current_amplitude * np.cos(tt) for tt in t_array[:frame+1]])
        line_ib.set_data(t_array[:frame+1], [current_amplitude * np.cos(tt - 2*np.pi/3) for tt in t_array[:frame+1]])
        line_ic.set_data(t_array[:frame+1], [current_amplitude * np.cos(tt - 4*np.pi/3) for tt in t_array[:frame+1]])
        
        # 更新各相MMF
        mmf_a = mmf_a_time[frame]
        mmf_b = mmf_b_time[frame]
        mmf_c = mmf_c_time[frame]
        line_mmf_a.set_data(alpha, mmf_a)
        line_mmf_b.set_data(alpha, mmf_b)
        line_mmf_c.set_data(alpha, mmf_c)
        
        # 更新合成MMF
        mmf_total = mmf_total_time[frame]
        line_mmf_total.set_data(alpha, mmf_total)
        
        # 更新频谱
        Order, P1 = spectrum_time[frame]
        # 清除旧的stem图
        ax4.clear()
        if len(Order) > 0 and len(P1) > 0:
            markerline, stemlines, baseline = ax4.stem(Order, P1, linefmt='k-', markerfmt='ko', basefmt='k-')
            plt.setp(stemlines, 'linewidth', 1)
            plt.setp(markerline, 'markersize', 3)
        ax4.set_xlim([0, 30])
        P1_max = np.max(P1) if len(P1) > 0 else 1
        ax4.set_ylim([0, P1_max * 1.1])
        ax4.set_xlabel('Order', fontsize=11)
        ax4.set_ylabel('MMF Amplitude [A]', fontsize=11)
        ax4.set_title('合成MMF频谱（随时间变化）', fontsize=12, fontweight='bold')
        ax4.grid(True, alpha=0.3)
        
        # 更新时间文本
        time_text.set_text(f'Time: t = {t:.3f} rad (ωt = {np.degrees(t):.1f}°) | '
                          f'ia = {ia:.3f} A, ib = {ib:.3f} A, ic = {ic:.3f} A | '
                          f'ia+ib+ic = {ia+ib+ic:.6f} A')
        
        return (line_ia, line_ib, line_ic, line_mmf_a, line_mmf_b, line_mmf_c, 
                line_mmf_total, time_text)
    
    # 创建动画
    anim = animation.FuncAnimation(
        fig, animate, frames=num_frames, interval=animation_interval, 
        blit=False, repeat=True
    )
    
    plt.suptitle(f'{excitation_type.capitalize()} Excitation - 三相MMF时间演化动画', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.show()
    
    return anim


if __name__ == '__main__':

    coil_pitch_y = 1
    Q = 12
    torque_winding_info = {
        'A': {
            'upper_slot_conductor_list': [1,4,7,10],
            'lower_slot_conductor_list': [2,5,8,11],
            'suspension_slot_conductor_list': [1,10,-4,-7],
        },
        'B': {
            'upper_slot_conductor_list': [2,5,8,11],
            'lower_slot_conductor_list': [3,6,9,12],
            'suspension_slot_conductor_list': [2,5,-11,-8],
        },
        'C': {
            'upper_slot_conductor_list': [3,6,9,12],
            'lower_slot_conductor_list': [4,7,10,13],
            'suspension_slot_conductor_list': [3,12,-6,-9],
        }
    }

    # 计算各相的激励模式（用于不同激励模式方法对比）
    torque_upper_A, torque_lower_A, suspension_upper_A, suspension_lower_A = compute_phase_excitation_patterns(
        coil_pitch_y, torque_winding_info['A'], Q)
    torque_upper_B, torque_lower_B, suspension_upper_B, suspension_lower_B = compute_phase_excitation_patterns(
        coil_pitch_y, torque_winding_info['B'], Q)
    torque_upper_C, torque_lower_C, suspension_upper_C, suspension_lower_C = compute_phase_excitation_patterns(
        coil_pitch_y, torque_winding_info['C'], Q)
    
    # 使用硬编码的值（与计算结果一致）
    A_phase_torque_excitation_upper_layer = [-1,  0,  0,  -1,  0,  0,  -1,  0,  0,  -1,  0,  0]
    A_phase_torque_excitation_lower_layer = [ 0,  1,  0,   0,  1,  0,   0,  1,  0,   0,  1,  0]
    A_phase_suspension_excitation_upper_layer = [ 1,  0,  0,  -1,  0,  0,  -1,  0,  0,   1,  0,  0]
    A_phase_suspension_excitation_lower_layer = [ 0, -1,  0,   0,  1,  0,   0,  1,  0,   0, -1,  0]
    
    # 准备激励模式字典（用于不同激励模式方法）
    torque_excitation_dict = {
        'A': {
            'upper': np.array(A_phase_torque_excitation_upper_layer),
            'lower': np.array(A_phase_torque_excitation_lower_layer)
        },
        'B': {
            'upper': np.array(torque_upper_B),
            'lower': np.array(torque_lower_B)
        },
        'C': {
            'upper': np.array(torque_upper_C),
            'lower': np.array(torque_lower_C)
        }
    }
    
    suspension_excitation_dict = {
        'A': {
            'upper': np.array(A_phase_suspension_excitation_upper_layer),
            'lower': np.array(A_phase_suspension_excitation_lower_layer)
        },
        'B': {
            'upper': np.array(suspension_upper_B),
            'lower': np.array(suspension_lower_B)
        },
        'C': {
            'upper': np.array(suspension_upper_C),
            'lower': np.array(suspension_lower_C)
        }
    }
    
    # ========== 转矩激励完整诊断 ==========
    print("\n" + "="*80)
    print("转矩激励 - 完整MMF诊断（与MATLAB对比）")
    print("="*80)
    print("使用硬编码的I_U模式（MATLAB模式）：I_U = [-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0]")
    print("使用circshift方法计算B、C相MMF（MATLAB方法）")
    print("电流值：i_u=1.0, i_v=-0.5, i_w=-0.5")
    
    torque_results = compute_and_plot_mmf_diagnostic(
        I_U_upper=None,  # 将使用硬编码模式
        I_U_lower=None,
        phase_excitation_dict=torque_excitation_dict,  # 用于对比
        excitation_type='torque',
        use_hardcoded_I_U=True,  # 使用MATLAB的硬编码模式
        use_circshift=True,  # 使用MATLAB的circshift方法
        Qs=Q,
        N=30000,
        i_u=1.0,
        i_v=-0.5,
        i_w=-0.5
    )
    
    # 数值验证：与MATLAB预期值对比
    print("\n" + "-"*80)
    print("数值验证（转矩激励，circshift方法）：")
    print("-"*80)
    print(f"U相MMF: max={np.max(torque_results['MMF_U_pha']):.4f}, "
          f"min={np.min(torque_results['MMF_U_pha']):.4f}, "
          f"mean={np.mean(torque_results['MMF_U_pha']):.6f}")
    if torque_results['MMF_total_circshift'] is not None:
        print(f"合成MMF: max={np.max(torque_results['MMF_total_circshift']):.4f}, "
              f"min={np.min(torque_results['MMF_total_circshift']):.4f}, "
              f"mean={np.mean(torque_results['MMF_total_circshift']):.6f}")
        print(f"MATLAB预期范围: ylim([-2, 3])")
        print(f"电流验证: i_u+i_v+i_w = {1.0 + (-0.5) + (-0.5):.6f}")
    
    # ========== 悬浮激励完整诊断 ==========
    print("\n" + "="*80)
    print("悬浮激励 - 完整MMF诊断")
    print("="*80)
    print("使用从winding_info计算的激励模式")
    print("使用circshift方法计算B、C相MMF")
    
    suspension_results = compute_and_plot_mmf_diagnostic(
        I_U_upper=np.array(A_phase_suspension_excitation_upper_layer),
        I_U_lower=np.array(A_phase_suspension_excitation_lower_layer),
        phase_excitation_dict=suspension_excitation_dict,
        excitation_type='suspension',
        use_hardcoded_I_U=False,  # 使用计算的激励模式
        use_circshift=True,  # 使用circshift方法
        Qs=Q,
        N=30000,
        i_u=1.0,
        i_v=-0.5,
        i_w=-0.5
    )
    
    # 数值验证
    print("\n" + "-"*80)
    print("数值验证（悬浮激励，circshift方法）：")
    print("-"*80)
    print(f"U相MMF: max={np.max(suspension_results['MMF_U_pha']):.4f}, "
          f"min={np.min(suspension_results['MMF_U_pha']):.4f}, "
          f"mean={np.mean(suspension_results['MMF_U_pha']):.6f}")
    if suspension_results['MMF_total_circshift'] is not None:
        print(f"合成MMF: max={np.max(suspension_results['MMF_total_circshift']):.4f}, "
              f"min={np.min(suspension_results['MMF_total_circshift']):.4f}, "
              f"mean={np.mean(suspension_results['MMF_total_circshift']):.6f}")
    
    print("\n" + "="*80)
    print("诊断完成！所有中间步骤已显示在subplots中。")
    print("="*80)



