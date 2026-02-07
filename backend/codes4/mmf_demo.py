"""
计算按照槽口分布的电流密度和MMF（磁动势）
从MATLAB脚本 mmf_demo.m 转换而来
简化版本：使用硬编码值，显示每相MMF和合成MMF
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple


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
    I_U_lower: np.ndarray = None,
    Qs: int = 12,
    N: int = 30000,
    tooth_width: float = 0.1745,
    slot_open: float = 0.3491
) -> Tuple[np.ndarray, np.ndarray]:
    """
    核心MMF计算函数，严格按照MATLAB逻辑
    
    Args:
        I_U_upper: U相上层导体的激励模式数组（长度为Qs）
        I_U_lower: U相下层导体的激励模式数组（长度为Qs），如果为None则自动计算为上层反向
        Qs: 定子槽数
        N: 采样步数
        tooth_width: 齿宽 [rad]
        slot_open: 槽开口 [rad]
    
    Returns:
        (alpha, MMF_U_pha)
        - alpha: 机械角度数组 [rad]
        - MMF_U_pha: U相MMF数组
    """
    # ========== 输入参数 ==========
    alpha = np.linspace(0, 2*np.pi, N)  # 机械角度（共N个采样点）
    # 每个槽的参考角，注意最后一个值 2pi 与 0 重合
    alpha_ref = np.linspace(0, 2*np.pi, Qs + 1)
    
    # 确保是numpy数组
    I_U_upper = np.array(I_U_upper)
    
    # 处理下层激励
    if I_U_lower is None:
        # 如果未提供下层激励，则自动计算为上层反向（MATLAB逻辑：-A_U(i)）
        I_U_lower = -I_U_upper
    else:
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
        lower_slot = (i + 1) % Qs
        theta_start_lower = (alpha_ref[lower_slot] + tooth_width) % (2 * np.pi)
        theta_end_lower = (theta_start_lower + slot_open) % (2 * np.pi)
        
        for j in range(N):
            theta = alpha[j] % (2 * np.pi)
            if theta_start_lower < theta_end_lower:
                if (theta >= theta_start_lower) and (theta < theta_end_lower):
                    # 使用指定的下层激励值
                    MMF_U_lower[j, i] = A_U_lower[i]
                else:
                    MMF_U_lower[j, i] = 0
            else:
                if (theta >= theta_start_lower) or (theta < theta_end_lower):
                    MMF_U_lower[j, i] = A_U_lower[i]
                else:
                    MMF_U_lower[j, i] = 0
    
    # 总的 U 相 MMF 为上层和下层的贡献之和（各槽独立叠加）
    MMF_U_pha = np.sum(MMF_U_upper + MMF_U_lower, axis=1)
    
    return alpha, MMF_U_pha


def plot_u_phase_mmf_check():
    """
    绘制转矩和悬浮力的U相磁动势对比图，供检查
    """
    # 硬编码参数
    Qs = 12
    N = 30000
    tooth_width = 0.1745
    slot_open = 0.3491
    
    # 硬编码的激励模式
    I_U_torque_upper = np.array([-1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0])
    I_U_torque_lower = np.array([0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0])
    
    # 悬浮激励：硬编码值
    I_U_suspension_upper = np.array([1, 0, 0, -1, 0, 0, -1, 0, 0, 1, 0, 0])
    I_U_suspension_lower = np.array([0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0])
    
    # 计算转矩激励的U相MMF
    alpha, MMF_U_torque = compute_u_phase_mmf_core(
        I_U_upper=I_U_torque_upper,
        I_U_lower=I_U_torque_lower,
        Qs=Qs, N=N, tooth_width=tooth_width, slot_open=slot_open
    )
    
    # 计算悬浮激励的U相MMF
    alpha, MMF_U_suspension = compute_u_phase_mmf_core(
        I_U_upper=I_U_suspension_upper,
        I_U_lower=I_U_suspension_lower,
        Qs=Qs, N=N, tooth_width=tooth_width, slot_open=slot_open
    )
    
    # 创建图形
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # 转矩激励U相MMF
    ax = axes[0]
    ax.plot(alpha, MMF_U_torque, 'r-', linewidth=2, label='转矩激励U相MMF')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=12)
    ax.set_ylabel('MMF [A]', fontsize=12)
    ax.set_title('转矩激励 - U相磁动势', fontsize=14, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=11)
    
    # 添加数值信息到标题
    ax.text(0.02, 0.98, 
            f'max={np.max(MMF_U_torque):.4f}, min={np.min(MMF_U_torque):.4f}, mean={np.mean(MMF_U_torque):.6f}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 悬浮激励U相MMF
    ax = axes[1]
    ax.plot(alpha, MMF_U_suspension, 'b-', linewidth=2, label='悬浮激励U相MMF')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=12)
    ax.set_ylabel('MMF [A]', fontsize=12)
    ax.set_title('悬浮激励 - U相磁动势', fontsize=14, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=11)
    
    # 添加数值信息到标题
    ax.text(0.02, 0.98, 
            f'max={np.max(MMF_U_suspension):.4f}, min={np.min(MMF_U_suspension):.4f}, mean={np.mean(MMF_U_suspension):.6f}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.show()
    
    # 打印详细信息
    print("\n" + "="*80)
    print("U相磁动势检查结果")
    print("="*80)
    
    print("\n转矩激励 - U相MMF:")
    print("-"*80)
    print(f"上层激励: {I_U_torque_upper}")
    print(f"下层激励: {I_U_torque_lower}")
    print(f"最大值: {np.max(MMF_U_torque):.6f} A")
    print(f"最小值: {np.min(MMF_U_torque):.6f} A")
    print(f"平均值: {np.mean(MMF_U_torque):.6f} A")
    print(f"峰峰值: {np.max(MMF_U_torque) - np.min(MMF_U_torque):.6f} A")
    
    print("\n悬浮激励 - U相MMF:")
    print("-"*80)
    print(f"上层激励: {I_U_suspension_upper}")
    print(f"下层激励: {I_U_suspension_lower}")
    print(f"最大值: {np.max(MMF_U_suspension):.6f} A")
    print(f"最小值: {np.min(MMF_U_suspension):.6f} A")
    print(f"平均值: {np.mean(MMF_U_suspension):.6f} A")
    print(f"峰峰值: {np.max(MMF_U_suspension) - np.min(MMF_U_suspension):.6f} A")


def plot_mmf_comparison():
    """
    绘制MMF对比图：左边列是转矩激励，右边列是悬浮激励
    每列显示：U相、V相、W相MMF和合成MMF
    """
    # 硬编码参数
    Qs = 12
    N = 30000
    tooth_width = 0.1745
    slot_open = 0.3491
    i_u = 1.0
    i_v = -0.5
    i_w = -0.5
    
    # 硬编码的激励模式（MATLAB模式）
    I_U_torque_upper = np.array([-1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0])
    I_U_torque_lower = np.array([0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0])
    
    # 悬浮激励：硬编码值
    I_U_suspension_upper = np.array([1, 0, 0, -1, 0, 0, -1, 0, 0, 1, 0, 0])
    I_U_suspension_lower = np.array([0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0])

    # 计算转矩激励的U相MMF
    alpha, MMF_U_torque = compute_u_phase_mmf_core(
        I_U_upper=I_U_torque_upper,
        I_U_lower=I_U_torque_lower,
        Qs=Qs, N=N, tooth_width=tooth_width, slot_open=slot_open
    )
    
    # 计算悬浮激励的U相MMF
    alpha, MMF_U_suspension = compute_u_phase_mmf_core(
        I_U_upper=I_U_suspension_upper,
        I_U_lower=I_U_suspension_lower,
        Qs=Qs, N=N, tooth_width=tooth_width, slot_open=slot_open
    )
    
    # 使用circshift计算V、W相MMF（MATLAB方法）
    Delta_N = round(N / 3)
    
    # 转矩激励的V、W相
    MMF_V_torque = np.roll(MMF_U_torque, Delta_N)
    MMF_W_torque = np.roll(MMF_V_torque, Delta_N)
    MMF_total_torque = MMF_U_torque * i_u + MMF_V_torque * i_v + MMF_W_torque * i_w
    
    # 悬浮激励的V、W相
    MMF_V_suspension = np.roll(MMF_U_suspension, Delta_N)
    MMF_W_suspension = np.roll(MMF_V_suspension, Delta_N)
    MMF_total_suspension = MMF_U_suspension * i_u + MMF_V_suspension * i_v + MMF_W_suspension * i_w
    
    # 创建图形：左边列是转矩激励，右边列是悬浮激励
    fig, axes = plt.subplots(4, 2, figsize=(16, 14))
    
    # 左边列：转矩激励
    # U相MMF
    ax = axes[0, 0]
    ax.plot(alpha, MMF_U_torque, 'r-', linewidth=1.5, label='U相')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('转矩激励 - U相MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # V相MMF
    ax = axes[1, 0]
    ax.plot(alpha, MMF_V_torque, 'b-', linewidth=1.5, label='V相')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('转矩激励 - V相MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # W相MMF
    ax = axes[2, 0]
    ax.plot(alpha, MMF_W_torque, 'g-', linewidth=1.5, label='W相')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('转矩激励 - W相MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # 合成MMF
    ax = axes[3, 0]
    ax.plot(alpha, MMF_total_torque, 'k-', linewidth=2, label='合成MMF')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('转矩激励 - 合成MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # 右边列：悬浮激励
    # U相MMF
    ax = axes[0, 1]
    ax.plot(alpha, MMF_U_suspension, 'r-', linewidth=1.5, label='U相')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('悬浮激励 - U相MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # V相MMF
    ax = axes[1, 1]
    ax.plot(alpha, MMF_V_suspension, 'b-', linewidth=1.5, label='V相')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('悬浮激励 - V相MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # W相MMF
    ax = axes[2, 1]
    ax.plot(alpha, MMF_W_suspension, 'g-', linewidth=1.5, label='W相')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('悬浮激励 - W相MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # 合成MMF
    ax = axes[3, 1]
    ax.plot(alpha, MMF_total_suspension, 'k-', linewidth=2, label='合成MMF')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('悬浮激励 - 合成MMF', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    plt.show()
    
    # 打印数值信息
    print("\n" + "="*80)
    print("转矩激励MMF数值信息：")
    print("-"*80)
    print(f"U相MMF: max={np.max(MMF_U_torque):.4f}, min={np.min(MMF_U_torque):.4f}, mean={np.mean(MMF_U_torque):.6f}")
    print(f"V相MMF: max={np.max(MMF_V_torque):.4f}, min={np.min(MMF_V_torque):.4f}, mean={np.mean(MMF_V_torque):.6f}")
    print(f"W相MMF: max={np.max(MMF_W_torque):.4f}, min={np.min(MMF_W_torque):.4f}, mean={np.mean(MMF_W_torque):.6f}")
    print(f"合成MMF: max={np.max(MMF_total_torque):.4f}, min={np.min(MMF_total_torque):.4f}, mean={np.mean(MMF_total_torque):.6f}")
    print(f"电流: i_u={i_u}, i_v={i_v}, i_w={i_w}, 和={i_u+i_v+i_w:.6f}")
    
    print("\n" + "="*80)
    print("悬浮激励MMF数值信息：")
    print("-"*80)
    print(f"U相MMF: max={np.max(MMF_U_suspension):.4f}, min={np.min(MMF_U_suspension):.4f}, mean={np.mean(MMF_U_suspension):.6f}")
    print(f"V相MMF: max={np.max(MMF_V_suspension):.4f}, min={np.min(MMF_V_suspension):.4f}, mean={np.mean(MMF_V_suspension):.6f}")
    print(f"W相MMF: max={np.max(MMF_W_suspension):.4f}, min={np.min(MMF_W_suspension):.4f}, mean={np.mean(MMF_W_suspension):.6f}")
    print(f"合成MMF: max={np.max(MMF_total_suspension):.4f}, min={np.min(MMF_total_suspension):.4f}, mean={np.mean(MMF_total_suspension):.6f}")
    print(f"电流: i_u={i_u}, i_v={i_v}, i_w={i_w}, 和={i_u+i_v+i_w:.6f}")


def compute_mmf_matlab_style():
    """
    严格按照MATLAB脚本 mmf_demo.m 的逻辑计算MMF
    使用MATLAB的硬编码值：I_U = [-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0]
    """
    # ========== 输入参数（与MATLAB完全一致） ==========
    Qs = 12
    N = 30000
    alpha_u = 2 * np.pi / Qs
    alpha = np.linspace(0, 2*np.pi, N)
    alpha_ref = np.linspace(0, 2*np.pi, Qs + 1)  # 对应MATLAB: 0:alpha_u:2*pi
    
    # 三相电流幅值（MATLAB: [i_u, i_v, i_w] = deal(1,-0.5,-0.5)）
    i_u = 1.0
    i_v = -0.5
    i_w = -0.5
    
    # U相上层导体的电流密度（MATLAB: I_U = [-1  1  0  -1  1  0  -1  1  0  -1  1  0]）
    I_U = np.array([-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0])
    ZQ_U = np.ones(Qs)  # 每槽导体数
    A_U = I_U * ZQ_U    # U相上层导体电流密度
    
    # 参数：齿宽与槽开口
    tooth_width = 0.1745
    slot_open = 0.3491
    
    # ========== 计算 U 相 MMF（严格按照MATLAB逻辑） ==========
    MMF_U_upper = np.zeros((N, Qs))
    MMF_U_lower = np.zeros((N, Qs))
    
    for i in range(Qs):
        # 上层导体贡献
        theta_start_upper = (alpha_ref[i] + tooth_width) % (2 * np.pi)
        theta_end_upper = (theta_start_upper + slot_open) % (2 * np.pi)
        
        for j in range(N):
            theta = alpha[j] % (2 * np.pi)
            if theta_start_upper < theta_end_upper:
                if (theta >= theta_start_upper) and (theta < theta_end_upper):
                    MMF_U_upper[j, i] = A_U[i]
                else:
                    MMF_U_upper[j, i] = 0
            else:
                if (theta >= theta_start_upper) or (theta < theta_end_upper):
                    MMF_U_upper[j, i] = A_U[i]
                else:
                    MMF_U_upper[j, i] = 0
        
        # 下层导体贡献（MATLAB: MMF_U_lower(j,i) = -A_U(i)）
        lower_slot = (i + 1) % Qs  # MATLAB: mod(i, Qs) + 1
        theta_start_lower = (alpha_ref[lower_slot] + tooth_width) % (2 * np.pi)
        theta_end_lower = (theta_start_lower + slot_open) % (2 * np.pi)
        
        for j in range(N):
            theta = alpha[j] % (2 * np.pi)
            if theta_start_lower < theta_end_lower:
                if (theta >= theta_start_lower) and (theta < theta_end_lower):
                    MMF_U_lower[j, i] = -A_U[i]  # MATLAB: 负号代表返向电流
                else:
                    MMF_U_lower[j, i] = 0
            else:
                if (theta >= theta_start_lower) or (theta < theta_end_lower):
                    MMF_U_lower[j, i] = -A_U[i]
                else:
                    MMF_U_lower[j, i] = 0
    
    # 总的 U 相 MMF（MATLAB: MMF_U_pha = sum(MMF_U_upper + MMF_U_lower, 2)）
    MMF_U_pha = np.sum(MMF_U_upper + MMF_U_lower, axis=1)
    
    # ========== 将U相平移得到V、W相MMF（MATLAB方法） ==========
    Delta_N = round(N / 3)  # MATLAB: Delta_N = round(N/3)
    MMF_V_pha = np.roll(MMF_U_pha, Delta_N)  # MATLAB: circshift
    MMF_W_pha = np.roll(MMF_V_pha, Delta_N)
    
    # 计算三相总MMF
    MMF_tot = MMF_U_pha * i_u + MMF_V_pha * i_v + MMF_W_pha * i_w
    
    return alpha, MMF_U_pha, MMF_V_pha, MMF_W_pha, MMF_tot


def compare_matlab_vs_current():
    """
    对比MATLAB风格的计算结果和当前Python脚本的计算结果
    """
    # MATLAB风格的计算
    alpha_matlab, MMF_U_matlab, MMF_V_matlab, MMF_W_matlab, MMF_tot_matlab = compute_mmf_matlab_style()
    
    # 当前Python脚本的计算（使用用户给定的硬编码值）
    Qs = 12
    N = 30000
    tooth_width = 0.1745
    slot_open = 0.3491
    
    I_U_torque_upper = np.array([-1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0])
    I_U_torque_lower = np.array([0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0])
    
    alpha_current, MMF_U_current = compute_u_phase_mmf_core(
        I_U_upper=I_U_torque_upper,
        I_U_lower=I_U_torque_lower,
        Qs=Qs, N=N, tooth_width=tooth_width, slot_open=slot_open
    )
    
    # 使用circshift计算V、W相
    Delta_N = round(N / 3)
    MMF_V_current = np.roll(MMF_U_current, Delta_N)
    MMF_W_current = np.roll(MMF_V_current, Delta_N)
    
    i_u = 1.0
    i_v = -0.5
    i_w = -0.5
    MMF_tot_current = MMF_U_current * i_u + MMF_V_current * i_v + MMF_W_current * i_w
    
    # 创建对比图
    fig, axes = plt.subplots(3, 2, figsize=(16, 12))
    
    # U相MMF对比
    ax = axes[0, 0]
    ax.plot(alpha_matlab, MMF_U_matlab, 'r-', linewidth=2, label='MATLAB风格 (I_U=[-1,1,0,...])', alpha=0.7)
    ax.plot(alpha_current, MMF_U_current, 'b--', linewidth=2, label='当前Python (I_U_upper=[-1,0,0,...])', alpha=0.7)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('U相MMF对比', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # V相MMF对比
    ax = axes[1, 0]
    ax.plot(alpha_matlab, MMF_V_matlab, 'r-', linewidth=2, label='MATLAB风格', alpha=0.7)
    ax.plot(alpha_current, MMF_V_current, 'b--', linewidth=2, label='当前Python', alpha=0.7)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('V相MMF对比', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # W相MMF对比
    ax = axes[2, 0]
    ax.plot(alpha_matlab, MMF_W_matlab, 'r-', linewidth=2, label='MATLAB风格', alpha=0.7)
    ax.plot(alpha_current, MMF_W_current, 'b--', linewidth=2, label='当前Python', alpha=0.7)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('W相MMF对比', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # 合成MMF对比
    ax = axes[0, 1]
    ax.plot(alpha_matlab, MMF_tot_matlab, 'r-', linewidth=2, label='MATLAB风格合成MMF', alpha=0.7)
    ax.plot(alpha_current, MMF_tot_current, 'b--', linewidth=2, label='当前Python合成MMF', alpha=0.7)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF [A]', fontsize=11)
    ax.set_title('合成MMF对比', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    # 差异分析
    diff_U = MMF_U_matlab - MMF_U_current
    diff_V = MMF_V_matlab - MMF_V_current
    diff_W = MMF_W_matlab - MMF_W_current
    diff_tot = MMF_tot_matlab - MMF_tot_current
    
    ax = axes[1, 1]
    ax.plot(alpha_matlab, diff_U, 'r-', linewidth=1.5, label='U相差异', alpha=0.7)
    ax.plot(alpha_matlab, diff_V, 'b-', linewidth=1.5, label='V相差异', alpha=0.7)
    ax.plot(alpha_matlab, diff_W, 'g-', linewidth=1.5, label='W相差异', alpha=0.7)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF差异 [A]', fontsize=11)
    ax.set_title('各相MMF差异', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    ax = axes[2, 1]
    ax.plot(alpha_matlab, diff_tot, 'k-', linewidth=2, label='合成MMF差异')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('机械角度 [rad]', fontsize=11)
    ax.set_ylabel('MMF差异 [A]', fontsize=11)
    ax.set_title('合成MMF差异', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 2*np.pi])
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    plt.show()
    
    # 打印详细对比信息
    print("\n" + "="*80)
    print("MATLAB风格 vs 当前Python脚本 - MMF对比分析")
    print("="*80)
    
    print("\n【MATLAB风格】（I_U = [-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0]）")
    print("-"*80)
    print(f"U相MMF: max={np.max(MMF_U_matlab):.6f}, min={np.min(MMF_U_matlab):.6f}, mean={np.mean(MMF_U_matlab):.6f}")
    print(f"V相MMF: max={np.max(MMF_V_matlab):.6f}, min={np.min(MMF_V_matlab):.6f}, mean={np.mean(MMF_V_matlab):.6f}")
    print(f"W相MMF: max={np.max(MMF_W_matlab):.6f}, min={np.min(MMF_W_matlab):.6f}, mean={np.mean(MMF_W_matlab):.6f}")
    print(f"合成MMF: max={np.max(MMF_tot_matlab):.6f}, min={np.min(MMF_tot_matlab):.6f}, mean={np.mean(MMF_tot_matlab):.6f}")
    
    print("\n【当前Python脚本】（I_U_upper = [-1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0]）")
    print("-"*80)
    print(f"U相MMF: max={np.max(MMF_U_current):.6f}, min={np.min(MMF_U_current):.6f}, mean={np.mean(MMF_U_current):.6f}")
    print(f"V相MMF: max={np.max(MMF_V_current):.6f}, min={np.min(MMF_V_current):.6f}, mean={np.mean(MMF_V_current):.6f}")
    print(f"W相MMF: max={np.max(MMF_W_current):.6f}, min={np.min(MMF_W_current):.6f}, mean={np.mean(MMF_W_current):.6f}")
    print(f"合成MMF: max={np.max(MMF_tot_current):.6f}, min={np.min(MMF_tot_current):.6f}, mean={np.mean(MMF_tot_current):.6f}")
    
    print("\n【差异分析】")
    print("-"*80)
    print(f"U相MMF最大差异: {np.max(np.abs(diff_U)):.6f} A")
    print(f"V相MMF最大差异: {np.max(np.abs(diff_V)):.6f} A")
    print(f"W相MMF最大差异: {np.max(np.abs(diff_W)):.6f} A")
    print(f"合成MMF最大差异: {np.max(np.abs(diff_tot)):.6f} A")
    print(f"\n原因：激励模式不同！")
    print(f"  MATLAB: I_U = [-1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0]")
    print(f"  当前Python: I_U_upper = [-1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0]")
    print(f"            I_U_lower = [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0]")


if __name__ == '__main__':
    # 对比MATLAB风格和当前Python脚本的结果
    compare_matlab_vs_current()
