

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
import sys
import os
import matplotlib

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

# ========== 绕组参数 ==========
m = 3  # 相数
coil_pitch_y = 1  # 线圈节距
Q = 12  # 总槽数
p = 4  # 极对数
ps = 5  # 每极每相槽数（可选参数）

# ========== 匝函数偏置 ==========
# 注意：匝函数偏置不是空间相位差，而是初始累积匝数值（通常为0）
# 空间相位差已经通过connection_star_raw_dict中的角度实现
turn_func_bias_phase_a = 0.0
turn_func_bias_phase_b = 0.0
turn_func_bias_phase_c = 0.0

# ========== 连接星形图 ==========
# 格式：{'相名': [(电角度(度), 槽号), ...]}
# 空间相位差通过角度体现：
# - A相：0°, 360°, 720°, 1080° (槽 1, 4, 7, 10)
# - B相：120°, 480°, 840°, 1200° (槽 2, 5, 8, 11) - 120°空间偏移
# - C相：240°, 600°, 960°, 1320° (槽 3, 6, 9, 12) - 240°空间偏移
connection_star_raw_dict = {
    'A': [(0.0, 1), (360.0, 4), (720.0, 7), (1080.0, 10)],
    'B': [(120.0, 2), (480.0, 5), (840.0, 8), (1200.0, 11)],  # 120°空间偏移
    'C': [(240.0, 3), (600.0, 6), (960.0, 9), (1320.0, 12)]   # 240°空间偏移
}

# 注意：coil_pitch_y已在文件开头定义
torque_winding_info = {
    'A': {
        'upper_slot_conductor_list': [1,4,7,10],
        'lower_slot_conductor_list': [2,5,8,11],
        'suspension_slot_conductor_list': [1,10,-4,-7],
    },
    'B': {
        'upper_slot_conductor_list': [2,5,8,11],
        'lower_slot_conductor_list': [3,6,9,12],
        'suspension_slot_conductor_list': [2,11,-5,-8],
    },
    'C': {
        'upper_slot_conductor_list': [3,6,9,12],
        'lower_slot_conductor_list': [4,7,10,13],
        'suspension_slot_conductor_list': [3,12,-6,-9],
    }
}

torque_upper_A, torque_lower_A, suspension_upper_A, suspension_lower_A = compute_phase_excitation_patterns(coil_pitch_y, torque_winding_info['A'], Q)
print(torque_upper_A, torque_lower_A, suspension_upper_A, suspension_lower_A)

A_phase_torque_excitation_upper_layer = [-1,  0,  0,  -1,  0,  0,  -1,  0,  0,  -1,  0,  0  ]
A_phase_torque_excitation_lower_layer = [ 0,  1,  0,   0,  1,  0,   0,  1,  0,   0,  1,  0  ]
A_phase_suspension_excitation_upper_layer = [ 1,  0,  0,  -1,  0,  0,  -1,  0,  0,   1,  0,  0  ]
A_phase_suspension_excitation_lower_layer = [ 0, -1,  0,   0,  1,  0,   0,  1,  0,   0, -1,  0  ]

# 配置matplotlib中文字体支持
def setup_chinese_font():
    """
    配置matplotlib以支持中文字体显示
    """
    # Windows系统常见中文字体列表（按优先级排序）
    chinese_fonts = [
        'Microsoft YaHei',      # 微软雅黑
        'SimHei',               # 黑体
        'SimSun',               # 宋体
        'KaiTi',                # 楷体
        'FangSong',             # 仿宋
    ]
    
    # 尝试设置中文字体
    font_found = False
    for font_name in chinese_fonts:
        try:
            # 设置字体
            plt.rcParams['font.sans-serif'] = [font_name]
            # 解决负号显示问题
            plt.rcParams['axes.unicode_minus'] = False
            font_found = True
            # 使用英文输出，避免终端编码问题
            print(f"Chinese font configured: {font_name}")
            break
        except Exception as e:
            continue
    
    if not font_found:
        # 如果找不到中文字体，尝试使用系统默认字体
        try:
            plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
            plt.rcParams['axes.unicode_minus'] = False
            print("Warning: Chinese font not found, using default font (Chinese may not display correctly)")
        except:
            pass

# 立即配置中文字体
setup_chinese_font()


# 添加路径以导入其他模块的函数
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from winding_mmf_moving_wave_three_phase import (
    normalize_slot_number,
    extract_phase_winding_data,
    create_detailed_turn_function,
    interpolate_turn_function
)


def compute_winding_function(
    connection_star_raw_dict: Dict[str, List[Tuple[float, int]]],
    phase_name: str,
    coil_pitch_y: int,
    Q: int,
    p: int,
    turn_func_bias: float = 0.0,
    turns_per_coil_side: int = 5,
    phase_dpnv_grouping_list: Optional[List] = None,
    num_points: int = 1000
) -> Tuple[np.ndarray, np.ndarray]:
    """
    计算绕组函数（Winding Function）
    
    绕组函数 = 匝函数 - 匝函数在整个圆周上的平均值
    
    Args:
        connection_star_raw_dict: 连接星形图字典
        phase_name: 相名，如 'Aa'
        coil_pitch_y: 线圈节距
        Q: 总槽数
        p: 极对数
        turn_func_bias: 匝函数初始偏置值
        turns_per_coil_side: 每个线圈边的匝数
        phase_dpnv_grouping_list: 可选的DPNV分组列表
        num_points: 采样点数
    
    Returns:
        (电角度数组, 绕组函数值数组)
    """
    # 提取该相的绕组数据
    upper_layer, lower_layer, conn_list, rev_upper, rev_lower = extract_phase_winding_data(
        connection_star_raw_dict, phase_name, coil_pitch_y, Q, phase_dpnv_grouping_list
    )
    
    # 创建详细的匝函数（在槽位置空间）
    detailed_tf = create_detailed_turn_function(
        upper_layer, lower_layer, conn_list, rev_upper, rev_lower,
        turn_func_bias, Q, turns_per_coil_side
    )
    
    # 电角度范围：0 到 p*360 度
    alpha_deg = np.linspace(0, p * 360, num_points)
    
    # 将电角度转换为槽位置
    # 相邻槽之间的电角度差 = 360° * p / Q
    # 槽位置 = 电角度 / (360 * p / Q) = 电角度 * Q / (360 * p)
    slot_positions = alpha_deg * Q / (360 * p)
    
    # 处理周期性：确保槽位置在 [0, Q] 范围内
    slot_positions = slot_positions % Q
    
    # 插值匝函数到电角度对应的槽位置
    turn_function_values = interpolate_turn_function(detailed_tf, slot_positions)
    
    # 计算匝函数在整个圆周上的平均值
    # 平均值 = 积分(匝函数) / 周长
    # 周长 = p * 360 度（电角度）
    # 使用数值积分
    circumference = p * 360  # 电角度
    # 使用trapezoid替代已弃用的trapz
    if hasattr(np, 'trapezoid'):
        turn_function_mean = np.trapezoid(turn_function_values, alpha_deg) / circumference
    else:
        turn_function_mean = np.trapz(turn_function_values, alpha_deg) / circumference
    
    # 计算绕组函数
    winding_function = turn_function_values - turn_function_mean
    
    return alpha_deg, winding_function


def plot_winding_function(
    alpha_deg: np.ndarray,
    winding_function: np.ndarray,
    phase_name: str = 'A',
    save_path: Optional[str] = None
):
    """
    绘制绕组函数
    
    Args:
        alpha_deg: 电角度数组（度）
        winding_function: 绕组函数值数组
        phase_name: 相名
        save_path: 保存路径（可选）
    """
    plt.figure(figsize=(12, 6))
    plt.plot(alpha_deg, winding_function, 'b-', linewidth=2, label=f'{phase_name}相绕组函数')
    plt.xlabel('电角度 (度)', fontsize=12)
    plt.ylabel('绕组函数值', fontsize=12)
    plt.title(f'{phase_name}相绕组函数 (Winding Function)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # 添加零线
    plt.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    
    # 标记极对数周期
    for i in range(1, p + 1):
        plt.axvline(x=i * 360, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图形已保存到: {save_path}")
    
    plt.show()


# ========== 主程序：计算并绘制A相绕组函数 ==========
if __name__ == '__main__':
    print("=" * 80)
    print("计算A相绕组函数")
    print("=" * 80)
    
    # ========== 示例：计算A相激励模式 ==========
    print("\n计算A相激励模式...")
    torque_upper, torque_lower, suspension_upper, suspension_lower = compute_phase_excitation_patterns(
        coil_pitch_y=coil_pitch_y,
        phase_winding_info=torque_winding_info['A'],
        Q=Q
    )
    
    print(f"扭矩激励上层: {torque_upper}")
    print(f"扭矩激励下层: {torque_lower}")
    print(f"悬挂激励上层: {suspension_upper}")
    print(f"悬挂激励下层: {suspension_lower}")
    
    # 验证结果
    print("\n验证结果:")
    print(f"扭矩激励上层是否匹配: {torque_upper == A_phase_torque_excitation_upper_layer}")
    print(f"扭矩激励下层是否匹配: {torque_lower == A_phase_torque_excitation_lower_layer}")
    print(f"悬挂激励上层是否匹配: {suspension_upper == A_phase_suspension_excitation_upper_layer}")
    print(f"悬挂激励下层是否匹配: {suspension_lower == A_phase_suspension_excitation_lower_layer}")
    print("=" * 80)
    
    # 计算A相绕组函数
    alpha_deg, winding_function_a = compute_winding_function(
        connection_star_raw_dict=connection_star_raw_dict,
        phase_name='Aa',
        coil_pitch_y=coil_pitch_y,
        Q=Q,
        p=p,
        turn_func_bias=turn_func_bias_phase_a,
        turns_per_coil_side=5,
        num_points=2000
    )
    
    # 打印一些统计信息
    print(f"\n绕组函数统计信息（A相）:")
    print(f"  电角度范围: {alpha_deg[0]:.1f}° 到 {alpha_deg[-1]:.1f}°")
    print(f"  绕组函数最大值: {np.max(winding_function_a):.4f}")
    print(f"  绕组函数最小值: {np.min(winding_function_a):.4f}")
    print(f"  绕组函数平均值: {np.mean(winding_function_a):.6f} (应该接近0)")
    print(f"  绕组函数标准差: {np.std(winding_function_a):.4f}")
    
    # 绘制A相绕组函数
    plot_winding_function(alpha_deg, winding_function_a, phase_name='A')
    
    print("\n绘图完成！")
