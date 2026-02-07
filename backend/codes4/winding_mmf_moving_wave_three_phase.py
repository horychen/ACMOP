"""
三相绕组合成MMF（磁动势）计算与可视化脚本

本脚本从头开始实现三相绕组的合成MMF计算，包括：
1. 匝函数（Turn Function）计算
2. 单相MMF计算
3. 三相合成MMF计算
4. 静态分布图可视化
5. 行波动画可视化
6. 基波分量提取

MMF计算原理：
MMF(θ) = ia * turn_function_a(θ) + ib * turn_function_b(θ) + ic * turn_function_c(θ)

其中：
- turn_function(θ) 是匝函数，表示累积匝数随空间位置的变化
- ia, ib, ic 是三相电流的瞬时值
- θ 是空间位置（以槽为单位）
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from typing import Dict, List, Tuple, Optional
import math

# ============================================================================
# 核心计算函数
# ============================================================================

def normalize_slot_number(slot_number: int, Q: int) -> int:
    """
    将槽号归一化到[1, Q]范围内（处理周期性）
    
    Args:
        slot_number: 槽号
        Q: 总槽数
    
    Returns:
        归一化后的槽号
    """
    if slot_number > Q:
        return slot_number - Q
    elif slot_number < 1:
        return slot_number + Q
    else:
        return slot_number


def extract_phase_winding_data(
    connection_star_raw_dict: Dict[str, List[Tuple[float, int]]],
    phase_name: str,
    coil_pitch_y: int,
    Q: int,
    phase_dpnv_grouping_list: Optional[List] = None
) -> Tuple[List[int], List[int], List[int], Optional[List[int]], Optional[List[int]]]:
    """
    从connection_star_raw_dict中提取指定相的绕组数据
    
    Args:
        connection_star_raw_dict: 连接星形图字典，格式为 {'A': [(角度, 槽号), ...], ...}
        phase_name: 相名，如 'Aa' 或 'Bb' 或 'Cc'
        coil_pitch_y: 线圈节距
        Q: 总槽数
        phase_dpnv_grouping_list: 可选的DPNV分组列表（用于反向激励）
    
    Returns:
        (上层槽列表, 下层槽列表, 连接方向列表, 反向激励上层列表, 反向激励下层列表)
    """
    # 处理反向激励
    if phase_dpnv_grouping_list is not None:
        reversed_excitation_upper_layer = [abs(int(el)) for el in phase_dpnv_grouping_list]
        reversed_excitation_lower_layer = [
            normalize_slot_number(abs(int(el)) + coil_pitch_y, Q) 
            for el in phase_dpnv_grouping_list
        ]
    else:
        reversed_excitation_upper_layer = None
        reversed_excitation_lower_layer = None
    
    # 提取该相的槽位信息
    phase_upper_layer_list = []
    phase_lower_layer_list = []
    connection_list = []
    
    # phase_name的第一个字符对应connection_star_raw_dict的键（如'A'、'B'、'C'）
    phase_key = phase_name[0]
    
    if phase_key in connection_star_raw_dict:
        phase_data = connection_star_raw_dict[phase_key]
        
        for angle, slot_number in phase_data:
            # 上层槽
            phase_upper_layer_list.append(slot_number)
            
            # 下层槽（根据线圈节距计算）
            lower_slot = normalize_slot_number(slot_number + coil_pitch_y, Q)
            phase_lower_layer_list.append(lower_slot)
            
            # 连接方向：如果phase_name的最后一个字符是小写（如'a'），表示反向连接
            if len(phase_name) > 1 and phase_name[-1].islower():
                connection_list.append(-1)
            else:
                connection_list.append(1)
    
    return (
        phase_upper_layer_list,
        phase_lower_layer_list,
        connection_list,
        reversed_excitation_upper_layer,
        reversed_excitation_lower_layer
    )


def compute_turn_function(
    upper_layer_list: List[int],
    lower_layer_list: List[int],
    connection_list: List[int],
    rev_upper: Optional[List[int]],
    rev_lower: Optional[List[int]],
    turn_func_bias: float,
    Q: int,
    turns_per_coil_side: int = 5
) -> Dict[int, float]:
    """
    计算匝函数（Turn Function）
    
    匝函数表示累积匝数随空间位置的变化。当经过一个槽时：
    - 如果该槽有上层导体，匝数会跳跃
    - 如果该槽有下层导体，匝数会再次跳跃
    
    Args:
        upper_layer_list: 上层槽号列表
        lower_layer_list: 下层槽号列表
        connection_list: 连接方向列表（1或-1）
        rev_upper: 反向激励的上层槽列表
        rev_lower: 反向激励的下层槽列表
        turn_func_bias: 匝函数初始偏置值
        Q: 总槽数
        turns_per_coil_side: 每个线圈边的匝数
    
    Returns:
        字典：{槽号: 累积匝数}
    """
    accumulated_turns = turn_func_bias
    turn_function = {}
    
    # 初始值（槽0，即槽1之前）
    turn_function[0] = accumulated_turns
    
    # 遍历所有槽
    for slot_number in range(1, Q + 1):
        # 槽开始时的匝数
        turn_function[slot_number] = accumulated_turns
        
        # 处理上层导体
        if slot_number in upper_layer_list:
            sign = connection_list[upper_layer_list.index(slot_number)]
            # 如果该槽在反向激励列表中，符号取反
            if rev_upper is not None and slot_number in rev_upper:
                sign *= -1
            accumulated_turns += sign * turns_per_coil_side
        
        # 处理下层导体
        if slot_number in lower_layer_list:
            # 下层导体的符号与上层相反
            sign = -1 * connection_list[lower_layer_list.index(slot_number)]
            # 如果该槽在反向激励列表中，符号取反
            if rev_lower is not None and slot_number in rev_lower:
                sign *= -1
            accumulated_turns += sign * turns_per_coil_side
    
    return turn_function


def create_detailed_turn_function(
    upper_layer_list: List[int],
    lower_layer_list: List[int],
    connection_list: List[int],
    rev_upper: Optional[List[int]],
    rev_lower: Optional[List[int]],
    turn_func_bias: float,
    Q: int,
    turns_per_coil_side: int = 5
) -> List[Tuple[float, float]]:
    """
    创建详细的匝函数，包含槽内的所有跳跃点
    
    这对于精确插值很重要，因为匝函数在槽内会有跳跃
    
    Args:
        参数同compute_turn_function
    
    Returns:
        列表：[(位置, 匝数值), ...]
    """
    accumulated_turns = turn_func_bias
    detailed_tf = []
    
    # 初始点
    detailed_tf.append((0.0, accumulated_turns))
    
    for slot_number in range(1, Q + 1):
        slot_start = float(slot_number - 1)
        slot_end = float(slot_number)
        
        # 槽开始前（如果不是第一个槽）
        if slot_number > 1:
            detailed_tf.append((slot_start - 1e-6, accumulated_turns))
        detailed_tf.append((slot_start, accumulated_turns))
        
        # 处理上层导体（在槽内约0.3位置）
        if slot_number in upper_layer_list:
            sign = connection_list[upper_layer_list.index(slot_number)]
            if rev_upper is not None and slot_number in rev_upper:
                sign *= -1
            accumulated_turns += sign * turns_per_coil_side
            detailed_tf.append((slot_start + 0.3, accumulated_turns))
        
        # 处理下层导体（在槽内约0.7位置）
        if slot_number in lower_layer_list:
            sign = -1 * connection_list[lower_layer_list.index(slot_number)]
            if rev_lower is not None and slot_number in rev_lower:
                sign *= -1
            accumulated_turns += sign * turns_per_coil_side
            detailed_tf.append((slot_start + 0.7, accumulated_turns))
        
        # 槽结束前
        detailed_tf.append((slot_end - 1e-6, accumulated_turns))
        detailed_tf.append((slot_end, accumulated_turns))
    
    return detailed_tf


def interpolate_turn_function(
    detailed_tf: List[Tuple[float, float]],
    positions: np.ndarray
) -> np.ndarray:
    """
    从详细匝函数插值到指定位置
    
    使用阶梯插值（前向填充），因为匝函数是阶梯函数
    
    Args:
        detailed_tf: 详细匝函数列表 [(位置, 值), ...]
        positions: 需要插值的位置数组
    
    Returns:
        插值后的匝函数值数组
    """
    pos_array = np.array([p[0] for p in detailed_tf])
    val_array = np.array([p[1] for p in detailed_tf])
    
    values = []
    for pos in positions:
        # 找到最后一个位置 <= pos 的点
        idx = np.searchsorted(pos_array, pos, side='right') - 1
        idx = max(0, min(idx, len(val_array) - 1))
        values.append(val_array[idx])
    
    return np.array(values)


def compute_mmf_spatial_distribution(
    connection_star_raw_dict: Dict[str, List[Tuple[float, int]]],
    coil_pitch_y: int,
    Q: int,
    ia: float,
    ib: float,
    ic: float,
    phase_Aa_dpnv_grouping_list: Optional[List] = None,
    phase_Bb_dpnv_grouping_list: Optional[List] = None,
    phase_Cc_dpnv_grouping_list: Optional[List] = None,
    turn_func_bias_phase_a: float = 0.0,
    turn_func_bias_phase_b: float = 0.0,
    turn_func_bias_phase_c: float = 0.0,
    bool_double_layer_winding: bool = True,
    num_points_per_slot: int = 10,
    turns_per_coil_side: int = 5
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    计算三相MMF的空间分布
    
    Args:
        connection_star_raw_dict: 连接星形图字典
        coil_pitch_y: 线圈节距
        Q: 总槽数
        ia, ib, ic: 三相电流瞬时值
        phase_Aa_dpnv_grouping_list: A相DPNV分组（可选）
        phase_Bb_dpnv_grouping_list: B相DPNV分组（可选）
        phase_Cc_dpnv_grouping_list: C相DPNV分组（可选）
        turn_func_bias_phase_a: A相匝函数偏置
        turn_func_bias_phase_b: B相匝函数偏置
        turn_func_bias_phase_c: C相匝函数偏置
        bool_double_layer_winding: 是否为双层绕组
        num_points_per_slot: 每个槽的采样点数
        turns_per_coil_side: 每个线圈边的匝数
    
    Returns:
        (位置数组, 合成MMF数组, A相MMF数组, B相MMF数组, C相MMF数组)
    """
    if not bool_double_layer_winding:
        raise NotImplementedError("单层绕组暂未实现")
    
    # 提取各相绕组数据
    upper_a, lower_a, conn_a, rev_upper_a, rev_lower_a = extract_phase_winding_data(
        connection_star_raw_dict, 'Aa', coil_pitch_y, Q, phase_Aa_dpnv_grouping_list
    )
    upper_b, lower_b, conn_b, rev_upper_b, rev_lower_b = extract_phase_winding_data(
        connection_star_raw_dict, 'Bb', coil_pitch_y, Q, phase_Bb_dpnv_grouping_list
    )
    upper_c, lower_c, conn_c, rev_upper_c, rev_lower_c = extract_phase_winding_data(
        connection_star_raw_dict, 'Cc', coil_pitch_y, Q, phase_Cc_dpnv_grouping_list
    )
    
    # 创建详细的匝函数
    detailed_tf_a = create_detailed_turn_function(
        upper_a, lower_a, conn_a, rev_upper_a, rev_lower_a,
        turn_func_bias_phase_a, Q, turns_per_coil_side
    )
    detailed_tf_b = create_detailed_turn_function(
        upper_b, lower_b, conn_b, rev_upper_b, rev_lower_b,
        turn_func_bias_phase_b, Q, turns_per_coil_side
    )
    detailed_tf_c = create_detailed_turn_function(
        upper_c, lower_c, conn_c, rev_upper_c, rev_lower_c,
        turn_func_bias_phase_c, Q, turns_per_coil_side
    )

    # ========== 调试：可视化匝函数 ==========
    # 设置 DEBUG_TURN_FUNCTION = True 来调试匝函数
    # 设置 DEBUG_AND_QUIT = True 来在调试后退出（不继续计算MMF）
    # 注意：如果 DEBUG_AND_QUIT=True，此函数会抛出异常，调用者需要处理
    DEBUG_TURN_FUNCTION = True
    DEBUG_AND_QUIT = False  # 设为True时，调试后退出，不继续计算MMF
    
    if DEBUG_TURN_FUNCTION:
        debug_turn_functions(
            detailed_tf_a, detailed_tf_b, detailed_tf_c,
            upper_a, lower_a, conn_a, rev_upper_a, rev_lower_a,
            upper_b, lower_b, conn_b, rev_upper_b, rev_lower_b,
            upper_c, lower_c, conn_c, rev_upper_c, rev_lower_c,
            Q
        )
        if DEBUG_AND_QUIT:
            print("\n调试完成。已退出（DEBUG_AND_QUIT=True）。")
            raise SystemExit("调试模式：已查看匝函数，退出程序。")
    
    # 创建空间位置网格
    num_total_points = Q * num_points_per_slot + 1
    slot_positions = np.linspace(0, Q, num_total_points)
    
    # 插值匝函数到网格
    turn_func_values_a = interpolate_turn_function(detailed_tf_a, slot_positions)
    turn_func_values_b = interpolate_turn_function(detailed_tf_b, slot_positions)
    turn_func_values_c = interpolate_turn_function(detailed_tf_c, slot_positions)
    
    # 计算各相MMF（电流 × 匝函数）
    phase_mmf_a = turn_func_values_a * ia
    phase_mmf_b = turn_func_values_b * ib
    phase_mmf_c = turn_func_values_c * ic
    
    # 计算合成MMF
    mmf_values = phase_mmf_a + phase_mmf_b + phase_mmf_c
    
    return slot_positions, mmf_values, phase_mmf_a, phase_mmf_b, phase_mmf_c


# ============================================================================
# 调试函数
# ============================================================================

def debug_turn_functions(
    detailed_tf_a: List[Tuple[float, float]],
    detailed_tf_b: List[Tuple[float, float]],
    detailed_tf_c: List[Tuple[float, float]],
    upper_a: List[int],
    lower_a: List[int],
    conn_a: List[int],
    rev_upper_a: Optional[List[int]],
    rev_lower_a: Optional[List[int]],
    upper_b: List[int],
    lower_b: List[int],
    conn_b: List[int],
    rev_upper_b: Optional[List[int]],
    rev_lower_b: Optional[List[int]],
    upper_c: List[int],
    lower_c: List[int],
    conn_c: List[int],
    rev_upper_c: Optional[List[int]],
    rev_lower_c: Optional[List[int]],
    Q: int
):
    """
    重新设计的调试函数：清晰展示匝函数与导体的关系
    
    每个图包含：
    1. 槽编号（1到Q，高度与y=0对齐）
    2. 导体标注：⊗（流入）或⊙（流出），黑色=上层，蓝色=下层
    3. 连接线：从1号槽开始，每遇到导体就发生匝函数跳变
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    
    # 设置字体为Times New Roman
    plt.rcParams['font.family'] = 'Times New Roman'
    
    # 创建图形：3个子图，每个显示一相的匝函数
    fig, axs = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
    fig.suptitle("Turn Functions with Conductor Positions", 
                 fontsize=24, fontweight='bold', fontfamily='Times New Roman')
    
    # 处理每一相
    turn_functions = [detailed_tf_a, detailed_tf_b, detailed_tf_c]
    phase_names = ['A', 'B', 'C']
    phase_colors = ['red', 'blue', 'green']
    
    phase_data = [
        (upper_a, lower_a, conn_a, rev_upper_a, rev_lower_a),
        (upper_b, lower_b, conn_b, rev_upper_b, rev_lower_b),
        (upper_c, lower_c, conn_c, rev_upper_c, rev_lower_c)
    ]
    
    # 定义颜色：上层黑色，下层蓝色
    upper_color = 'black'
    lower_color = 'blue'
    
    for idx, (tf, phase_name, phase_color, (upper_list, lower_list, conn_list, rev_upper, rev_lower)) in enumerate(
        zip(turn_functions, phase_names, phase_colors, phase_data)
    ):
        # 提取位置和值
        positions = np.array([p[0] for p in tf])
        values = np.array([p[1] for p in tf])
        
        # ========== 1. 绘制槽编号（与y=0对齐） ==========
        for slot_num in range(1, Q + 1):
            slot_center = slot_num - 0.5
            axs[idx].text(slot_center, 0, str(slot_num), 
                         fontsize=16, ha='center', va='center',
                         color='gray', weight='bold', fontfamily='Times New Roman',
                         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.7),
                         zorder=5)
        
        # ========== 2. 收集所有导体信息并按槽号排序 ==========
        conductors = []  # [(slot_num, layer, sign, color), ...]
        
        # 上层导体
        for slot_num in upper_list:
            sign = conn_list[upper_list.index(slot_num)]
            if rev_upper is not None and slot_num in rev_upper:
                sign *= -1
            conductors.append((slot_num, 'upper', sign, upper_color))
        
        # 下层导体
        for slot_num in lower_list:
            sign = -1 * conn_list[lower_list.index(slot_num)]
            if rev_lower is not None and slot_num in rev_lower:
                sign *= -1
            conductors.append((slot_num, 'lower', sign, lower_color))
        
        # 按槽号排序
        conductors.sort(key=lambda x: x[0])
        
        # ========== 3. 绘制匝函数阶梯图 ==========
        axs[idx].step(positions, values, where='post', linewidth=3, 
                     color=phase_color, label=f'Phase {phase_name} Turn Function', 
                     alpha=0.8, zorder=1)
        
        # ========== 4. 绘制导体和连接线 ==========
        # 从位置0开始，沿着匝函数曲线绘制连接线
        # 每遇到一个导体就绘制连接线到该导体
        
        # 创建导体位置字典，便于查找
        conductor_dict = {}  # {slot_num: [(layer, sign, color), ...]}
        for slot_num, layer, sign, color in conductors:
            if slot_num not in conductor_dict:
                conductor_dict[slot_num] = []
            conductor_dict[slot_num].append((layer, sign, color))
        
        # 从位置0开始，沿着匝函数曲线前进
        prev_x = 0.0
        prev_y = values[0] if len(values) > 0 else 0.0
        
        # 遍历所有槽（从1到Q）
        for slot_num in range(1, Q + 1):
            slot_start = slot_num - 1
            slot_center = slot_num - 0.5
            slot_end = slot_num
            
            # 找到该槽开始时的匝函数值
            pos_idx_start = np.searchsorted(positions, slot_start, side='right') - 1
            pos_idx_start = max(0, min(pos_idx_start, len(values) - 1))
            tf_value_at_start = values[pos_idx_start]
            
            # 如果该槽有导体，绘制连接线和导体
            if slot_num in conductor_dict:
                # 先绘制从上一个位置到槽开始的线
                axs[idx].plot([prev_x, slot_start], [prev_y, tf_value_at_start],
                            color='gray', linewidth=2, alpha=0.4, linestyle='--', zorder=2)
                
                # 处理该槽的所有导体（可能有上层和下层）
                for layer, sign, color in conductor_dict[slot_num]:
                    # 确定导体y位置
                    if layer == 'upper':
                        conductor_y = tf_value_at_start + 1.2  # 上层在匝函数上方
                    else:
                        conductor_y = tf_value_at_start - 1.2  # 下层在匝函数下方
                    
                    # 绘制从槽开始到导体的垂直线
                    axs[idx].plot([slot_center, slot_center], [tf_value_at_start, conductor_y],
                                color=color, linewidth=2.5, alpha=0.7, zorder=2)
                    
                    # 绘制导体符号
                    # ⊕ (U+2295) 表示流入，⊙ (U+2299) 表示流出
                    if sign > 0:  # 流入
                        symbol = '⊕'  # 带圈的加号（流入）
                    else:  # 流出
                        symbol = '⊙'  # 带点的圈（流出）
                    
                    # 使用scatter绘制圆，确保是正圆（不受坐标轴比例影响）
                    # scatter会自动处理aspect ratio，保持圆是正圆
                    axs[idx].scatter([slot_center], [conductor_y], 
                                   s=800,  # 圆的大小（以点为单位）
                                   facecolor='white', edgecolor=color, linewidth=3,
                                   zorder=3, marker='o')
                    
                    # 绘制符号
                    axs[idx].text(slot_center, conductor_y, symbol,
                                 fontsize=28, ha='center', va='center',
                                 color=color, weight='bold', fontfamily='Times New Roman',
                                 zorder=4)
                
                # 找到该槽结束时的匝函数值（导体导致跳变后）
                pos_idx_end = np.searchsorted(positions, slot_end, side='right') - 1
                pos_idx_end = max(0, min(pos_idx_end, len(values) - 1))
                tf_value_at_end = values[pos_idx_end]
                
                # 更新到槽结束位置
                prev_x = slot_end
                prev_y = tf_value_at_end
            else:
                # 该槽没有导体，直接绘制到槽结束
                pos_idx_end = np.searchsorted(positions, slot_end, side='right') - 1
                pos_idx_end = max(0, min(pos_idx_end, len(values) - 1))
                tf_value_at_end = values[pos_idx_end]
                
                # 绘制从上一个位置到槽结束的线
                axs[idx].plot([prev_x, slot_end], [prev_y, tf_value_at_end],
                            color='gray', linewidth=2, alpha=0.4, linestyle='--', zorder=2)
                
                prev_x = slot_end
                prev_y = tf_value_at_end
        
        # ========== 5. 绘制槽边界线 ==========
        for slot_num in range(Q + 1):
            axs[idx].axvline(x=slot_num, color='lightgray', linestyle=':', 
                           linewidth=1, alpha=0.5, zorder=0)
        
        # ========== 6. 设置图形属性 ==========
        axs[idx].set_ylabel(f'Accumulated Turns\nPhase {phase_name}', 
                           fontsize=18, fontweight='bold', fontfamily='Times New Roman')
        axs[idx].legend(loc='upper right', fontsize=15, prop={'family': 'Times New Roman'})
        axs[idx].grid(True, alpha=0.3, zorder=0)
        axs[idx].axhline(y=0, color='k', linestyle='-', linewidth=1.5, zorder=0)
        axs[idx].set_xlim(0, Q)
        
        # 设置刻度标签字体
        for label in axs[idx].get_xticklabels() + axs[idx].get_yticklabels():
            label.set_fontfamily('Times New Roman')
    
    # 设置x轴标签（只在最后一个子图）
    axs[2].set_xlabel('Slot Position', fontsize=19.5, fontweight='bold', fontfamily='Times New Roman')
    
    # 添加图例说明
    legend_text = '⊗ = Into (Black=Upper, Blue=Lower)  |  ⊙ = Out'
    fig.text(0.5, 0.02, legend_text, fontsize=14, ha='center', 
            fontfamily='Times New Roman',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.97])
    plt.show()
    
    # 打印详细信息到控制台
    print("\n" + "="*80)
    print("Turn Function Debug Info")
    print("="*80)
    
    for tf, phase_name, (upper_list, lower_list, conn_list, rev_upper, rev_lower) in zip(
        turn_functions, phase_names, phase_data
    ):
        print(f"\n[Phase {phase_name}]")
        print(f"  Upper layer slots: {sorted(upper_list)}")
        print(f"  Lower layer slots: {sorted(lower_list)}")
        print(f"  Turn function points: {len(tf)}")
        print(f"  Turn function range: [{min(v for _, v in tf):.2f}, {max(v for _, v in tf):.2f}]")
        
        # 打印每个导体的电流方向
        print(f"  Upper conductor current directions:")
        for slot_num in sorted(upper_list):
            sign = conn_list[upper_list.index(slot_num)]
            if rev_upper is not None and slot_num in rev_upper:
                sign *= -1
            direction = "Into (×)" if sign > 0 else "Out (·)"
            print(f"    Slot {slot_num}: {direction} (sign={sign:+d})")
        
        print(f"  Lower conductor current directions:")
        for slot_num in sorted(lower_list):
            sign = -1 * conn_list[lower_list.index(slot_num)]
            if rev_lower is not None and slot_num in rev_lower:
                sign *= -1
            direction = "Into (×)" if sign > 0 else "Out (·)"
            print(f"    Slot {slot_num}: {direction} (sign={sign:+d})")
    
    print("\n" + "="*80)


# ============================================================================
# 工具函数
# ============================================================================

def infer_pole_pairs_from_connection_star(
    connection_star_raw_dict: Dict[str, List[Tuple[float, int]]],
    Q: int
) -> int:
    """
    从connection_star_raw_dict推断极对数
    
    原理：相邻同相槽之间的电角度差 = 360° * p / Q
    
    Args:
        connection_star_raw_dict: 连接星形图字典
        Q: 总槽数
    
    Returns:
        极对数 p
    """
    if 'A' not in connection_star_raw_dict:
        raise ValueError("无法推断极对数：连接星形图中没有A相数据")
    
    phase_a_data = connection_star_raw_dict['A']
    if len(phase_a_data) < 2:
        raise ValueError("无法推断极对数：A相至少需要2个槽")
    
    # 获取A相的前两个槽
    angle1, slot1 = phase_a_data[0]
    angle2, slot2 = phase_a_data[1]
    
    # 计算电角度差和槽差
    angle_diff = abs(angle2 - angle1)
    slot_diff = abs(slot2 - slot1)
    
    if slot_diff == 0:
        raise ValueError("无法推断极对数：两个槽号相同")
    
    # 计算极对数：angle_diff / slot_diff = 360° * p / Q
    # 因此：p = (angle_diff / slot_diff) * Q / 360°
    p = (angle_diff / slot_diff) * Q / 360.0
    
    return int(round(p))


def extract_fundamental_component(
    slot_positions: np.ndarray,
    mmf_values: np.ndarray,
    p: int,
    Q: int
) -> np.ndarray:
    """
    使用FFT提取基波分量（p对极）
    
    Args:
        slot_positions: 空间位置数组
        mmf_values: MMF值数组
        p: 极对数
        Q: 总槽数
    
    Returns:
        基波分量数组
    """
    N = len(mmf_values)
    if N < 2:
        return mmf_values
    
    # 执行FFT
    fft_result = np.fft.fft(mmf_values)
    
    # 基波谐波索引为p
    fundamental_harmonic_idx = max(1, min(p, N // 2))
    
    # 创建滤波器，只保留基波分量
    fft_filtered = np.zeros_like(fft_result, dtype=complex)
    
    # 保留基波及其复共轭（用于实信号）
    idx_pos = fundamental_harmonic_idx
    idx_neg = (N - fundamental_harmonic_idx) % N
    
    if idx_pos < N:
        fft_filtered[idx_pos] = fft_result[idx_pos]
    if idx_neg < N and idx_neg != idx_pos:
        fft_filtered[idx_neg] = fft_result[idx_neg]
    
    # 逆FFT得到基波分量
    mmf_fundamental = np.real(np.fft.ifft(fft_filtered))
    
    return mmf_fundamental


# ============================================================================
# 可视化函数
# ============================================================================

def plot_mmf_static(
    connection_star_raw_dict: Dict[str, List[Tuple[float, int]]],
    coil_pitch_y: int,
    Q: int,
    ia: float,
    ib: float,
    ic: float,
    phase_Aa_dpnv_grouping_list: Optional[List] = None,
    phase_Bb_dpnv_grouping_list: Optional[List] = None,
    phase_Cc_dpnv_grouping_list: Optional[List] = None,
    turn_func_bias_phase_a: float = 0.0,
    turn_func_bias_phase_b: float = 0.0,
    turn_func_bias_phase_c: float = 0.0,
    bool_double_layer_winding: bool = True,
    num_points_per_slot: int = 10,
    save_path: Optional[str] = None
):
    """
    绘制静态MMF分布图
    
    Args:
        参数同compute_mmf_spatial_distribution
        save_path: 保存路径（可选）
    """
    # 计算MMF分布
    slot_positions, mmf_values, mmf_a, mmf_b, mmf_c = compute_mmf_spatial_distribution(
        connection_star_raw_dict, coil_pitch_y, Q,
        ia, ib, ic,
        phase_Aa_dpnv_grouping_list, phase_Bb_dpnv_grouping_list, phase_Cc_dpnv_grouping_list,
        turn_func_bias_phase_a, turn_func_bias_phase_b, turn_func_bias_phase_c,
        bool_double_layer_winding, num_points_per_slot
    )
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 绘制各相MMF
    ax.plot(slot_positions, mmf_a, 'r-', linewidth=1.5, alpha=0.7, label=f'Phase A (i={ia:.2f} A)')
    ax.plot(slot_positions, mmf_b, 'b-', linewidth=1.5, alpha=0.7, label=f'Phase B (i={ib:.2f} A)')
    ax.plot(slot_positions, mmf_c, 'g-', linewidth=1.5, alpha=0.7, label=f'Phase C (i={ic:.2f} A)')
    
    # 绘制合成MMF
    ax.plot(slot_positions, mmf_values, 'k-', linewidth=2.5, label='Composite MMF')
    
    # 设置图形属性
    ax.set_xlabel('Slot Position', fontsize=12)
    ax.set_ylabel('MMF (Ampere-turns)', fontsize=12)
    ax.set_title(f'Three-Phase MMF Distribution\nQ={Q}, y={coil_pitch_y}, ia={ia:.2f}, ib={ib:.2f}, ic={ic:.2f}', 
                 fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    
    # 标记槽位置
    for i in range(Q + 1):
        ax.axvline(x=i, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f'MMF分布图已保存至: {save_path}')
    else:
        plt.show()


def animate_mmf_traveling_wave(
    connection_star_raw_dict: Dict[str, List[Tuple[float, int]]],
    coil_pitch_y: int,
    Q: int,
    p: int,
    phase_Aa_dpnv_grouping_list: Optional[List] = None,
    phase_Bb_dpnv_grouping_list: Optional[List] = None,
    phase_Cc_dpnv_grouping_list: Optional[List] = None,
    turn_func_bias_phase_a: float = 0.0,
    turn_func_bias_phase_b: float = 0.0,
    turn_func_bias_phase_c: float = 0.0,
    bool_double_layer_winding: bool = True,
    current_amplitude: float = 1.0,
    frequency: float = 1.0,
    num_frames: int = 200,
    num_points_per_slot: int = 10
):
    """
    动画显示MMF行波
    
    Args:
        connection_star_raw_dict: 连接星形图字典
        coil_pitch_y: 线圈节距
        Q: 总槽数
        p: 极对数
        其他参数同compute_mmf_spatial_distribution
        current_amplitude: 电流幅值
        frequency: 电频率 (Hz)
        num_frames: 动画帧数
    """
    # 时间数组
    time_period = 1.0 / frequency if frequency > 0 else 1.0
    times = np.linspace(0, time_period, num_frames)
    
    # 创建图形（3个子图）
    fig = plt.figure(figsize=(15, 10))
    
    # 子图1：完整MMF（包含谐波）
    ax1 = plt.subplot(3, 1, 1)
    line_composite, = ax1.plot([], [], 'k-', linewidth=2.5, label='Composite MMF')
    line_phase_a, = ax1.plot([], [], 'r-', linewidth=1.5, alpha=0.7, label='Phase A')
    line_phase_b, = ax1.plot([], [], 'b-', linewidth=1.5, alpha=0.7, label='Phase B')
    line_phase_c, = ax1.plot([], [], 'g-', linewidth=1.5, alpha=0.7, label='Phase C')
    ax1.set_xlabel('Slot Position')
    ax1.set_ylabel('MMF (Ampere-turns)')
    ax1.set_title('Full MMF Distribution (with Harmonics)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    
    # 子图2：基波分量
    ax2 = plt.subplot(3, 1, 2)
    line_fundamental_composite, = ax2.plot([], [], 'k-', linewidth=2.5, label='Composite Fundamental')
    line_fundamental_a, = ax2.plot([], [], 'r-', linewidth=1.5, alpha=0.7, label='Phase A Fundamental')
    line_fundamental_b, = ax2.plot([], [], 'b-', linewidth=1.5, alpha=0.7, label='Phase B Fundamental')
    line_fundamental_c, = ax2.plot([], [], 'g-', linewidth=1.5, alpha=0.7, label='Phase C Fundamental')
    ax2.set_xlabel('Slot Position')
    ax2.set_ylabel('MMF (Ampere-turns)')
    ax2.set_title(f'Fundamental Component (p={p} pole pairs)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    
    # 子图3：三相电流
    ax3 = plt.subplot(3, 1, 3)
    line_ia, = ax3.plot([], [], 'r-', linewidth=2, label='i_a')
    line_ib, = ax3.plot([], [], 'b-', linewidth=2, label='i_b')
    line_ic, = ax3.plot([], [], 'g-', linewidth=2, label='i_c')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Current (A)')
    ax3.set_title('Three-Phase Currents')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    
    # 预计算MMF范围以设置合适的y轴范围
    omega = 2 * np.pi * frequency
    sample_times = np.linspace(0, time_period, 20)
    max_mmf_abs = 0
    max_fundamental_abs = 0
    
    for t in sample_times:
        ia_sample = current_amplitude * np.cos(omega * t)
        ib_sample = current_amplitude * np.cos(omega * t - 2*np.pi/3)
        ic_sample = current_amplitude * np.cos(omega * t - 4*np.pi/3)
        
        slot_pos_sample, mmf_sample, _, _, _ = compute_mmf_spatial_distribution(
            connection_star_raw_dict, coil_pitch_y, Q,
            ia_sample, ib_sample, ic_sample,
            phase_Aa_dpnv_grouping_list, phase_Bb_dpnv_grouping_list, phase_Cc_dpnv_grouping_list,
            turn_func_bias_phase_a, turn_func_bias_phase_b, turn_func_bias_phase_c,
            bool_double_layer_winding, num_points_per_slot
        )
        
        max_mmf_abs = max(max_mmf_abs, np.max(np.abs(mmf_sample)))
        
        mmf_fundamental_sample = extract_fundamental_component(slot_pos_sample, mmf_sample, p, Q)
        max_fundamental_abs = max(max_fundamental_abs, np.max(np.abs(mmf_fundamental_sample)))
    
    # 设置y轴范围
    mmf_range = max_mmf_abs * 1.3 if max_mmf_abs > 1e-6 else 10.0
    fundamental_range = max_fundamental_abs * 1.3 if max_fundamental_abs > 1e-6 else mmf_range
    
    ax1.set_xlim(0, Q)
    ax1.set_ylim(-mmf_range, mmf_range)
    ax2.set_xlim(0, Q)
    ax2.set_ylim(-fundamental_range, fundamental_range)
    ax3.set_xlim(0, time_period)
    ax3.set_ylim(-current_amplitude * 1.2, current_amplitude * 1.2)
    
    # 存储电流历史
    current_history = {'ia': [], 'ib': [], 'ic': [], 'time': []}
    
    # 预计算位置数组（只需要计算一次）
    slot_positions = np.linspace(0, Q, Q * num_points_per_slot + 1)
    
    def animate(frame):
        # 计算三相电流
        t = times[frame]
        omega = 2 * np.pi * frequency
        ia = current_amplitude * np.cos(omega * t)
        ib = current_amplitude * np.cos(omega * t - 2*np.pi/3)
        ic = current_amplitude * np.cos(omega * t - 4*np.pi/3)
        
        # 计算MMF分布
        _, mmf_values, mmf_a, mmf_b, mmf_c = compute_mmf_spatial_distribution(
            connection_star_raw_dict, coil_pitch_y, Q,
            ia, ib, ic,
            phase_Aa_dpnv_grouping_list, phase_Bb_dpnv_grouping_list, phase_Cc_dpnv_grouping_list,
            turn_func_bias_phase_a, turn_func_bias_phase_b, turn_func_bias_phase_c,
            bool_double_layer_winding, num_points_per_slot
        )
        
        # 提取基波分量
        mmf_fundamental_composite = extract_fundamental_component(slot_positions, mmf_values, p, Q)
        mmf_fundamental_a = extract_fundamental_component(slot_positions, mmf_a, p, Q)
        mmf_fundamental_b = extract_fundamental_component(slot_positions, mmf_b, p, Q)
        mmf_fundamental_c = extract_fundamental_component(slot_positions, mmf_c, p, Q)
        
        # 更新子图1
        line_composite.set_data(slot_positions, mmf_values)
        line_phase_a.set_data(slot_positions, mmf_a)
        line_phase_b.set_data(slot_positions, mmf_b)
        line_phase_c.set_data(slot_positions, mmf_c)
        
        # 更新子图2
        line_fundamental_composite.set_data(slot_positions, mmf_fundamental_composite)
        line_fundamental_a.set_data(slot_positions, mmf_fundamental_a)
        line_fundamental_b.set_data(slot_positions, mmf_fundamental_b)
        line_fundamental_c.set_data(slot_positions, mmf_fundamental_c)
        
        # 更新子图3
        current_history['ia'].append(ia)
        current_history['ib'].append(ib)
        current_history['ic'].append(ic)
        current_history['time'].append(t)
        
        # 只保留最后一个周期的数据
        if len(current_history['time']) > num_frames:
            current_history['ia'] = current_history['ia'][-num_frames:]
            current_history['ib'] = current_history['ib'][-num_frames:]
            current_history['ic'] = current_history['ic'][-num_frames:]
            current_history['time'] = current_history['time'][-num_frames:]
        
        line_ia.set_data(current_history['time'], current_history['ia'])
        line_ib.set_data(current_history['time'], current_history['ib'])
        line_ic.set_data(current_history['time'], current_history['ic'])
        
        return (line_composite, line_phase_a, line_phase_b, line_phase_c,
                line_fundamental_composite, line_fundamental_a, line_fundamental_b, line_fundamental_c,
                line_ia, line_ib, line_ic)
    
    # 创建动画
    anim = animation.FuncAnimation(
        fig, animate, frames=num_frames, interval=50, blit=True, repeat=True
    )
    
    plt.tight_layout()
    plt.show()
    
    return anim


# ============================================================================
# 主程序
# ============================================================================

if __name__ == '__main__':
    # ========== 绕组参数 ==========
    m = 3  # 相数
    coil_pitch_y = 1  # 线圈节距
    Q = 12  # 总槽数
    
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
    
    # ========== 推断极对数 ==========
    try:
        p = infer_pole_pairs_from_connection_star(connection_star_raw_dict, Q)
        print(f'推断的极对数 p = {p}')
    except Exception as e:
        print(f'警告：无法从connection_star_raw_dict推断极对数: {e}')
        print('使用默认值 p = 4。如果这不正确，请手动指定。')
        p = 4
    
    # ========== 选择运行模式 ==========
    # 模式1：绘制静态MMF分布图
    RUN_STATIC_PLOT = False
    if RUN_STATIC_PLOT:
        print('\n=== 绘制静态MMF分布图 ===')
        # 指定三相电流值
        ia_real = 1.0
        ib_real = -0.5
        ic_real = -0.5
        
        plot_mmf_static(
            connection_star_raw_dict,
            coil_pitch_y,
            Q,
            ia_real,
            ib_real,
            ic_real,
            turn_func_bias_phase_a=turn_func_bias_phase_a,
            turn_func_bias_phase_b=turn_func_bias_phase_b,
            turn_func_bias_phase_c=turn_func_bias_phase_c,
            bool_double_layer_winding=True,
            save_path=None  # 设置为路径字符串以保存图片
        )
    
    # 模式2：动画显示MMF行波（推荐用于验证）
    RUN_ANIMATION = True
    if RUN_ANIMATION:
        print('\n=== 启动MMF行波动画 ===')
        print('合成MMF应该形成沿气隙传播的行波。')
        print('关闭动画窗口以退出。')
        
        anim = animate_mmf_traveling_wave(
            connection_star_raw_dict,
            coil_pitch_y,
            Q,
            p,  # 极对数（用于基波分量提取）
            phase_Aa_dpnv_grouping_list=None,
            phase_Bb_dpnv_grouping_list=None,
            phase_Cc_dpnv_grouping_list=None,
            turn_func_bias_phase_a=turn_func_bias_phase_a,
            turn_func_bias_phase_b=turn_func_bias_phase_b,
            turn_func_bias_phase_c=turn_func_bias_phase_c,
            bool_double_layer_winding=True,
            current_amplitude=1.0,  # 电流幅值 (A)
            frequency=1.0,  # 电频率 (Hz)
            num_frames=200,  # 动画帧数
            num_points_per_slot=10  # 每个槽的采样点数
        )
