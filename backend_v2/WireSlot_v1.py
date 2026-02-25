import math
from typing import Dict, Tuple, Any

# 线规标准数据库 (裸铜直径, 含漆外径) 单位: mm
WIRE_DATA = {
    30: {'bare': 0.254, 'coated': 0.290},
    31: {'bare': 0.226, 'coated': 0.260},
    32: {'bare': 0.203, 'coated': 0.235},
}

def calc_stator_geometry(OD: float, ID: float, tooth_width: float, tooth_depth: float, yoke: float, liner: float, slots: int = 12) -> Dict[str, float]:
    """
    计算定子槽部的净几何空间（扣除槽纸）。
    """
    r_outer = OD / 2 - yoke
    
    # 平行齿下的槽口和槽底净宽度
    slot_opening_net = ((math.pi * ID) / slots) - tooth_width - (2 * liner)
    slot_bottom_net = ((math.pi * (r_outer * 2)) / slots) - tooth_width - (2 * liner)
    
    # 近似梯形净槽面积
    net_area = ((slot_opening_net + slot_bottom_net) / 2) * (tooth_depth - liner)
    gross_area = ((math.pi * ID / slots - tooth_width) + (math.pi * r_outer * 2 / slots - tooth_width)) / 2 * tooth_depth
    
    return {
        "slot_opening_net": max(0.0, slot_opening_net),
        "slot_bottom_net": max(0.0, slot_bottom_net),
        "net_area": max(0.0, net_area),
        "gross_area": max(0.0, gross_area)
    }

def calc_winding_capacity(net_area: float, gross_area: float, awg: int, orthocyclic_factor: float = 0.88, fill_heuristic: float = 0.7) -> Dict[str, Any]:
    """
    根据给定的线规和净槽面积，计算槽内能容纳的最大物理导线数及满槽率。
    """
    if awg not in WIRE_DATA:
        raise ValueError(f"不支持的 AWG 线规: {awg}")
        
    d_coated = WIRE_DATA[awg]['coated']
    d_bare = WIRE_DATA[awg]['bare']
    
    # 单根导线截面积
    a_bare = math.pi * (d_bare / 2) ** 2
    a_coated = math.pi * (d_coated / 2) ** 2
    
    # 基于整列系数的理论最大容量
    raw_capacity = (net_area * orthocyclic_factor) / a_coated
    # 结合双层排布与飞叉绕线工程极限的启发式折算
    z_slot = math.floor(raw_capacity * fill_heuristic)
    
    # 确保双层绕组下，导线数为偶数（两边线圈平分）
    if z_slot % 2 != 0:
        z_slot -= 1
        
    total_copper_area = z_slot * a_bare
    fill_factor = total_copper_area / gross_area
    
    return {
        "z_slot": z_slot,
        "fill_factor": fill_factor,
        "a_bare": a_bare
    }

def calc_thermal_load(ID: float, z_slot: int, a_bare: float, J: float, slots: int = 12) -> Dict[str, float]:
    """
    计算热负荷、单根电流与安匝数。
    """
    I = a_bare * J
    NI = z_slot * I
    
    # 线负荷 A (A/mm) = (总导线数 * 单根电流) / (pi * 内径)
    total_conductors = slots * z_slot
    A_loading = (total_conductors * I) / (math.pi * ID) 
    
    # 热负荷 AJ = A_loading * 10 (转化为 A/cm) * J (A/mm^2)
    AJ = (A_loading * 10) * J 
    
    return {
        "current_per_wire": I,
        "ampere_turns": NI,
        "thermal_load_AJ": AJ
    }

def calc_bemf_constants(La: float, yoke: float, z_slot: int, B_sat: float = 1.8, poles: int = 10, slots: int = 12) -> Dict[str, float]:
    """
    推算反电动势常数 (考虑轭部饱和瓶颈)。
    """
    # 每相串联线圈数 (3相电机)
    coils_per_phase = slots / 3
    # 每相总匝数 (单槽线圈边匝数为 z_slot / 2)
    turns_per_phase = coils_per_phase * (z_slot / 2)
    
    kw1 = 0.933 # 12S10P 基波绕组系数
    pole_pairs = poles / 2
    
    # 轭部截面积限制极限磁通 (La 和 yoke 单位需从 mm 转为 m 计算磁通)
    a_yoke_m2 = (yoke * 1e-3) * (La * 1e-3)
    phi_max = 2 * B_sat * a_yoke_m2 # 极最大有效磁通 Wb
    
    # 相电压峰值常数 V/(rad/s)
    Ke_ph_pk = turns_per_phase * kw1 * phi_max * pole_pairs
    
    # 线间电压有效值常数 V/(rad/s)
    Ke_LL_rms = math.sqrt(3) * (Ke_ph_pk / math.sqrt(2))
    
    # 线间有效值 V/krpm
    Ke_V_krpm = Ke_LL_rms * (1000 * 2 * math.pi / 60)
    
    # KV 值 (RPM/V) 基于线间峰值
    E_LL_pk = math.sqrt(3) * Ke_ph_pk
    KV = 1 / (E_LL_pk * (2 * math.pi / 60)) if E_LL_pk > 0 else 0
    
    return {
        "turns_per_phase": turns_per_phase,
        "phi_max_Wb": phi_max,
        "Ke_V_krpm": Ke_V_krpm,
        "KV_rpm_V": KV
    }

def calc_motor_losses(La: float, tooth_width: float, turns_per_phase: float, current: float, rpm: float, awg: int, poles: int = 10) -> Dict[str, float]:
    """
    预估 100°C 下的绕组铜损与简化的经验铁损。
    """
    a_bare = math.pi * (WIRE_DATA[awg]['bare'] / 2) ** 2
    rho_100C = 0.023 # 铜在100度时的电阻率 Ohm*mm^2/m
    
    # 估算平均匝长 (端部跨距简化计算)
    L_end = tooth_width * 2.5 + 0.5 
    L_turn_m = (2 * La + 2 * L_end) * 1e-3 
    
    # 计算相电阻
    L_ph_m = turns_per_phase * L_turn_m
    R_ph = rho_100C * (L_ph_m / a_bare)
    
    # 总铜损 (3相)
    P_cu = 3 * (current ** 2) * R_ph
    
    # 铁损估算 (基于之前的经验模型与轭部惩罚)
    freq = (rpm / 60) * (poles / 2)
    # 此处的铁损模型为一个高度简化的经验拟合，仅用于反映频率上升带来的非线性爆炸
    # 基础损耗系数对应 10000 RPM 时的 0.3W
    base_loss_coeff = 0.3 / (833 ** 1.5) 
    P_fe = base_loss_coeff * (freq ** 1.5)
    
    # 额外的高频饱和惩罚项 (针对极薄轭部在高速时的涡流恶化)
    if freq > 1000:
        P_fe *= 1.5
        
    return {
        "phase_resistance_100C": R_ph,
        "frequency_Hz": freq,
        "copper_loss_W": P_cu,
        "iron_loss_W": P_fe,
        "total_loss_W": P_cu + P_fe
    }

# ==========================================
# 后端调用示例
# ==========================================
if __name__ == "__main__":
    # 1. 初始输入参数
    params = {
        "OD": 13.0, "ID": 8.3, "tooth_width": 1.2, "tooth_depth": 2.1, 
        "yoke": 0.25, "liner": 0.1, "slots": 12, "poles": 10,
        "awg": 31, "J": 14.0, "La": 10.0, "B_sat": 1.8
    }
    
    # 2. 依次调用各模块执行计算管线
    geo = calc_stator_geometry(params["OD"], params["ID"], params["tooth_width"], params["tooth_depth"], params["yoke"], params["liner"], params["slots"])
    
    wind = calc_winding_capacity(geo["net_area"], geo["gross_area"], params["awg"])
    
    thermal = calc_thermal_load(params["ID"], wind["z_slot"], wind["a_bare"], params["J"], params["slots"])
    
    bemf = calc_bemf_constants(params["La"], params["yoke"], wind["z_slot"], params["B_sat"], params["poles"], params["slots"])
    
    # 计算 20000 RPM 时的损耗
    losses = calc_motor_losses(params["La"], params["tooth_width"], bemf["turns_per_phase"], thermal["current_per_wire"], rpm=20000, awg=params["awg"], poles=params["poles"])
    
    # 打印结果供核对
    print(f"--- 12S10P {params['OD']}mm 微电机后端验算结果 ---")
    print(f"满槽率: {wind['fill_factor']*100:.1f}%, 槽内总导线: {wind['z_slot']} 根")
    print(f"热负荷 (AJ): {thermal['thermal_load_AJ']:.0f} A^2/(cm*mm^2)")
    print(f"KV 值估算: {bemf['KV_rpm_V']:.0f} RPM/V")
    print(f"相电阻 (100°C): {losses['phase_resistance_100C']:.3f} Ohm")
    print(f"20k RPM 总损耗: {losses['total_loss_W']:.2f} W (铜损 {losses['copper_loss_W']:.2f}W, 铁损 {losses['iron_loss_W']:.2f}W)")
