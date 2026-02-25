import math

class MotorThermalAnalyzer:
    """
    电机绕组空间利用率与热负荷计算分析类
    专门针对微型电机（如 13mm 外径）的定子设计
    """

    def __init__(self, stator_od=13.0, rotor_od=8.0, air_gap=0.15, 
                 tooth_depth=2.1, tooth_width=1.2, slots=12, liner=0.1):
        """
        初始化电机几何参数
        :param stator_od: 定子外径 (mm)
        :param rotor_od: 转子外径 (mm)
        :param air_gap: 单边气隙 (mm)
        :param tooth_depth: 齿深 (mm)
        :param tooth_width: 平行齿齿宽 (mm)
        :param slots: 槽数
        :param liner: 槽绝缘纸厚度 (mm)
        """
        self.stator_od = stator_od
        self.rotor_od = rotor_od
        self.air_gap = air_gap
        self.stator_id = rotor_od + 2 * air_gap
        self.tooth_depth = tooth_depth
        self.tooth_width = tooth_width
        self.slots = slots
        self.liner = liner
        
        # 计算定子轭部厚度
        self.yoke_thickness = (stator_od - (self.stator_id + 2 * tooth_depth)) / 2

    def calculate_slot_area(self):
        """
        计算单个槽的几何面积 (mm^2)
        公式: (定子内部环形面积 - 总齿部面积) / 槽数
        """
        r_in = self.stator_id / 2
        r_out = r_in + self.tooth_depth
        total_annulus_area = math.pi * (r_out**2 - r_in**2)
        total_teeth_area = self.slots * self.tooth_width * self.tooth_depth
        return (total_annulus_area - total_teeth_area) / self.slots

    def get_wire_properties(self, awg_size):
        """
        获取标准 AWG 线规的物理特性
        :return: (裸铜直径, 含漆膜直径, 裸铜截面积)
        """
        awg_db = {
            30: {"d_bare": 0.254, "d_od": 0.290},
            31: {"d_bare": 0.226, "d_od": 0.260},
            32: {"d_bare": 0.203, "d_od": 0.235}
        }
        if awg_size not in awg_db:
            raise ValueError("Unsupported AWG size. Use 30, 31, or 32.")
            
        data = awg_db[awg_size]
        area_bare = math.pi * (data["d_bare"] / 2)**2
        return data["d_bare"], data["d_od"], area_bare

    def estimate_max_wires_in_slot(self, d_od):
        """
        基于几何干涉模拟槽内能放下的最大导线根数 (双层绕组)
        逻辑: 模拟两个相邻齿的线圈在槽中轴线汇合的过程
        """
        r_in = self.stator_id / 2
        slot_angle = (2 * math.pi) / self.slots
        tooth_angles = [-slot_angle / 2, slot_angle / 2]
        
        total_wires = 0
        
        # 模拟两个齿侧的布线
        for t_angle in tooth_angles:
            is_left_tooth = t_angle < 0
            side = 1 if is_left_tooth else -1
            
            # 遍历层 (向槽中心方向)
            for layer in range(6):
                # 距离齿中心线的横向距离
                h_offset = (self.tooth_width / 2 + self.liner + d_od / 2 + layer * d_od * 0.88)
                
                # 遍历行 (从槽口向槽底)
                for row in range(15):
                    v_offset = self.liner + d_od / 2 + row * d_od
                    if v_offset > self.tooth_depth - self.liner:
                        break
                    
                    r = r_in + v_offset
                    if h_offset >= r:
                        continue
                        
                    # 计算该点相对于齿中心的角度偏移
                    theta_offset = math.asin(h_offset / r) * side
                    final_angle = t_angle + theta_offset
                    
                    # 碰撞检测: 不得越过槽中心线 (0度)
                    if side == 1 and final_angle > -0.002: continue
                    if side == -1 and final_angle < 0.002: continue
                    
                    total_wires += 1
                    
        return total_wires

    def analyze_thermal_performance(self, awg_size, target_j=14.0):
        """
        执行完整的电机热负荷与激励分析
        :param awg_size: AWG 线规号
        :param target_j: 目标电流密度 (A/mm^2)
        :return: 包含所有关键指标的字典
        """
        d_bare, d_od, area_bare = self.get_wire_properties(awg_size)
        slot_area_geo = self.calculate_slot_area()
        
        # 1. 物理排布计算
        n_wires = self.estimate_max_wires_in_slot(d_od)
        
        # 2. 槽满率计算
        cu_area_total = n_wires * area_bare
        fill_factor_cu = (cu_area_total / slot_area_geo) * 100
        
        # 3. 电流指标
        current_per_wire = target_j * area_bare
        total_ampere_turns = n_wires * current_per_wire # 单槽总安匝
        
        # 4. 热负荷指标 (AJ)
        # 线负荷 A (A/cm) = (总导线数 * 单线电流) / (PI * Dsi_cm)
        stator_id_cm = self.stator_id / 10.0
        total_conductors = n_wires * self.slots
        line_load_a = (total_conductors * current_per_wire) / (math.pi * stator_id_cm)
        
        aj_value = line_load_a * target_j
        
        return {
            "awg": awg_size,
            "wires_per_slot": n_wires,
            "turns_per_tooth": n_wires / 2,
            "slot_fill_factor_cu": round(fill_factor_cu, 2),
            "current_per_wire_a": round(current_per_wire, 4),
            "ampere_turns_per_slot": round(total_ampere_turns, 2),
            "line_load_a_per_cm": round(line_load_a, 2),
            "thermal_load_aj": round(aj_value, 0),
            "yoke_thickness_mm": round(self.yoke_thickness, 3)
        }

# --- 调用示例 ---
if __name__ == "__main__":
    # 使用你的 13mm 电机参数
    analyzer = MotorThermalAnalyzer(
        stator_od=13.0, 
        rotor_od=8.0, 
        air_gap=0.15, 
        tooth_depth=2.1, 
        tooth_width=1.2, 
        slots=12
    )
    
    print(f"{'线规':<8} | {'单槽导线':<8} | {'安匝(NI)':<8} | {'槽满率%':<8} | {'AJ热负荷':<10}")
    print("-" * 60)
    
    for awg in [30, 31, 32]:
        res = analyzer.analyze_thermal_performance(awg)
        print(f"AWG {res['awg']:<4} | {res['wires_per_slot']:<12} | {res['ampere_turns_per_slot']:<12} | {res['slot_fill_factor_cu']:<10} | {res['thermal_load_aj']:<10}")

    # 针对 AWG 31 的详细诊断
    final_res = analyzer.analyze_thermal_performance(31)
    if final_res['thermal_load_aj'] > 1000:
        print(f"\n警告: 当前热负荷 AJ={final_res['thermal_load_aj']} 极高，需强化散热。")
    if final_res['yoke_thickness_mm'] < 0.3:
        print(f"警告: 轭部厚度 {final_res['yoke_thickness_mm']}mm 过薄，磁路可能严重饱和。")
