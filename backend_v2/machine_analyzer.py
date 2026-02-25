import math
from collections import OrderedDict

# 线规标准数据库 (裸铜直径, 含漆外径) 
# WIRE_DATA包含了不同AWG规格导线的：
# bare: 裸铜心直径 (mm)
# coated: 包含绝缘漆皮后的最大外径 (mm)
WIRE_DATA = {
    30: {'bare': 0.254, 'coated': 0.290},
    31: {'bare': 0.226, 'coated': 0.260},
    32: {'bare': 0.203, 'coated': 0.235},
}

class MotorPerformanceAnalyzer:
    """
    电机绕组空间利用率与热负荷计算分析类
    结合了几何干涉法 (v2) 与 完整的电气/热负荷/反电动势分析 (v1)。
    通过该类，我们可以从物理干涉上推断槽内导线最大数量，并在此基础上评估电机整体电磁性能。
    """
    def __init__(self, machine_dict: OrderedDict):
        self.machine_dict = machine_dict
        
        # gp (Geometry Parameters): 包含所有的几何尺寸参数，如外径、内径、齿宽、气隙等
        self.gp = machine_dict['geometry']
        
        # wp (Winding Parameters): 包含绕组配置信息，如极槽配合、线规、目标电流密度、转速等
        self.wp = machine_dict['winding']
        
        # 添加默认的槽绝缘衬套厚度（liner），如果未指定，默认给 0.1mm
        if 'liner' not in self.gp:
            self.gp['liner'] = 0.1
            
        # 如果未指定使用的线规AWG，默认采用AWG 31
        if 'awg' not in self.wp:
            self.wp['awg'] = 31

    def get_wire_properties(self, awg_size):
        """
        根据给定的AWG尺寸获取导线物理参数。
        返回: (裸铜直径 mm, 含漆外径 mm, 裸铜截面积 mm^2)
        """
        if awg_size not in WIRE_DATA:
            raise ValueError(f"不支持的 AWG 线规: {awg_size}")
        data = WIRE_DATA[awg_size]
        area_bare = math.pi * (data["bare"] / 2)**2
        return data["bare"], data["coated"], area_bare

    def calculate_slot_area(self):
        """
        计算单个槽的纯几何面积 (单位: mm^2)
        算法说明:
        1. 算出包含整个槽深在内的环形总面积 (Total Annulus Area)。
        2. 算出所有定子齿所占的矩形面积之和 (Total Teeth Area)。
        3. 用总环形面积减去所有定子齿的面积，得到所有槽的总面积。
        4. 除以槽数 slot_count 得到单槽几何面积。
        """
        # 计算定子内径 (转子外径加上气隙)
        stator_id = (self.gp['r_rotor_outer'] + self.gp['d_air_gap']) * 2
        r_in = stator_id / 2
        
        # r_out 为定子槽底的半径 (内径 + 齿深)
        r_out = r_in + self.gp['d_tooth']
        
        # 定子槽所在环形区域的总面积
        total_annulus_area = math.pi * (r_out**2 - r_in**2)
        
        # 所有齿的面积 (假设齿为矩形: 槽数 * 齿宽 * 齿深)
        total_teeth_area = self.wp['slot_count'] * self.gp['w_tooth'] * self.gp['d_tooth']
        
        # 计算单槽有效面积
        return (total_annulus_area - total_teeth_area) / self.wp['slot_count']

    def estimate_max_wires_in_slot(self, d_od):
        """
        基于几何干涉原理，模拟槽内能放下的最大导线根数。(主要针对双层绕组)
        参数:
        - d_od: 导线的含漆外径 (带绝缘皮的直径)
        
        算法逻辑:
        1. 从定子两个相邻的齿壁开始，一圈一圈（layer）往槽中心线排布导线。
        2. 检测每根导线的坐标 (圆周坐标: r, theta)。
        3. 如果导线与定子槽底/槽口/绝缘纸(liner)发生干涉，则该位置放置失败。
        4. 统计直到左右两侧排满或者到达槽中心汇合为止，槽内所能容纳的总导线数。
        """
        stator_id = (self.gp['r_rotor_outer'] + self.gp['d_air_gap']) * 2
        r_in = stator_id / 2
        
        # 确定每一个槽占据的弧度角，以及当前分析的左右槽壁所在的初始角度
        slot_angle = (2 * math.pi) / self.wp['slot_count']
        tooth_angles = [-slot_angle / 2, slot_angle / 2]
        
        total_wires = 0
        w_tooth = self.gp['w_tooth']
        d_tooth = self.gp['d_tooth']
        liner = self.gp['liner']
        
        # 遍历左侧齿壁和右侧齿壁，分别从两侧往中间填线
        for t_angle in tooth_angles:
            is_left_tooth = t_angle < 0
            side = 1 if is_left_tooth else -1
            
            # 假设一个槽径向最多可以排布 6 层 (Layer)
            for layer in range(6):
                # 计算当前层的水平偏移 (考虑半个齿宽、绝缘厚度、半径及层间互相嵌入的 0.88 排布紧密度系数)
                h_offset = (w_tooth / 2 + liner + d_od / 2 + layer * d_od * 0.88)
                
                # 在垂直方向上逐排（Row）放置导线，最大假设允许15排
                for row in range(15):
                    # 径向高度的垂直偏移 (沿半径方向)
                    v_offset = liner + d_od / 2 + row * d_od
                    
                    # 超过齿的深度，说明碰到定子背铁边缘，本层排不下了，直接跳到下一层
                    if v_offset > d_tooth - liner:
                        break
                    
                    r = r_in + v_offset
                    
                    # 几何干涉: 如果水平偏移导致的直线距离已经大于了弧线半径，这种位置是不物理存在的，跳过。
                    if h_offset >= r:
                        continue
                        
                    # 计算导线中心所在的绝对角度 (相对于齿壁的角度偏移)
                    theta_offset = math.asin(h_offset / r) * side
                    final_angle = t_angle + theta_offset
                    
                    # 碰撞检测：确保左侧排布的导线不会越过槽的几何中垂线 (越过了说明左右导线重叠干涉)
                    if side == 1 and final_angle > -0.002: continue
                    if side == -1 and final_angle < 0.002: continue
                    
                    # 成功放下一根导线
                    total_wires += 1
                    
        return total_wires

    def calc_bemf_constants(self, z_slot: int, B_sat: float = 1.8):
        """
        计算反电动势 (Back-EMF) 相关常数，包括每相反电动势与电机KV值。
        参数:
        - z_slot: 每槽平均导线数
        - B_sat: 假设的定子轭部饱和磁密设定值 (默认为 1.8 Tesla)
        
        重要公式与物理法则:
        1. 每相线圈匝数 = (槽数/3) * (每槽导线数/2)
        2. 最大磁通量 phi_max = 2 * 饱和磁通密度 * 轭部截面积
        3. 相反电动势峰值 Ke_ph_pk = 匝数 * 绕组系数 * 磁通 * 极对数
        """
        coils_per_phase = self.wp['slot_count'] / 3
        # 假设双层绕组，所以每个线圈边的匝数是 z_slot / 2
        turns_per_phase = coils_per_phase * (z_slot / 2)
        
        kw1 = 0.933 # 固定使用 12S10P 基波绕组系数
        pole_pairs = self.wp['pole_count'] / 2
        
        # 轭部厚度 = 定子外半径 - (内半径 + 气隙 + 齿深)
        yoke_thickness = self.gp['r_stator_outer'] - (self.gp['r_rotor_outer'] + self.gp['d_air_gap'] + self.gp['d_tooth'])
        
        # 定子轭部横截面积(m^2) = 轭部厚度 * 叠长
        a_yoke_m2 = (yoke_thickness * 1e-3) * (self.wp['l_stack'] * 1e-3)
        # 每个极通过的最大磁通 (韦伯)
        phi_max = 2 * B_sat * a_yoke_m2 
        
        # 反电动势常数 (Ke), Peak值，每相
        Ke_ph_pk = turns_per_phase * kw1 * phi_max * pole_pairs
        
        # 线反电动势 RMS值 = sqrt(3) * (相电压峰值 / sqrt(2))
        Ke_LL_rms = math.sqrt(3) * (Ke_ph_pk / math.sqrt(2))
        
        # 每1000 RPM对应的反电动势电压 (Volts / krpm)
        Ke_V_krpm = Ke_LL_rms * (1000 * 2 * math.pi / 60)
        
        # 线反电动势峰值 
        E_LL_pk = math.sqrt(3) * Ke_ph_pk
        # KV值 (RPM/V) 定义为 使得线反电动势峰值等于 1V 时对应的转速RPM
        KV = 1 / (E_LL_pk * (2 * math.pi / 60)) if E_LL_pk > 0 else 0
        
        return {
            "turns_per_phase": turns_per_phase,
            "phi_max_Wb": phi_max,
            "Ke_V_krpm": Ke_V_krpm,
            "KV_rpm_V": KV
        }

    def calc_motor_losses(self, turns_per_phase: float, current: float, rpm: float):
        """
        计算电机的定子铜损与粗略的铁损评估。
        参数:
        - turns_per_phase: 每相匝数
        - current: RMS 电流或设定的目标直流线排电流
        - rpm: 评估的工作转速
        """
        d_bare, _, area_bare = self.get_wire_properties(self.wp['awg'])
        rho_100C = 0.023 # 铜在 100°C 左右的近似电阻率 (Ohm·mm²/m)
        
        # 估算端部绕组伸出长度: 通常以2.5倍齿宽再加上0.5mm裕量作为端部悬置长度一半
        L_end = self.gp['w_tooth'] * 2.5 + 0.5 
        
        # 单匝线圈的总长度 (米): 涵盖两次进入定子铁芯的叠长 + 两端端部绕组伸出长度之和
        L_turn_m = (2 * self.wp['l_stack'] + 2 * L_end) * 1e-3 
        
        # 每相总线束长度
        L_ph_m = turns_per_phase * L_turn_m
        
        # 计算 100°C 时的相电阻
        R_ph = rho_100C * (L_ph_m / area_bare)
        
        # 三相总铜损 = 3 * I^2 * R
        P_cu = 3 * (current ** 2) * R_ph
        
        # 定子工作频率 (Electrical Frequency Hz)
        freq = (rpm / 60) * (self.wp['pole_count'] / 2)
        
        # 利用经验常数对最高转速/工作频率点做铁损 (Iron Loss) 的初步逼近预估 (Steinmetz 简化公式)
        base_loss_coeff = 0.3 / (833 ** 1.5) 
        P_fe = base_loss_coeff * (freq ** 1.5)
        
        # 当电频率过高时，额外增加损耗余量，模拟集肤效应以及更高的涡流损耗
        if freq > 1000:
            P_fe *= 1.5
            
        return {
            "phase_resistance_100C": R_ph,
            "frequency_Hz": freq,
            "copper_loss_W": P_cu,
            "iron_loss_W": P_fe,
            "total_loss_W": P_cu + P_fe
        }

    def analyze(self):
        """
        主执行函数，协调并集成所有的分析项。
        整合了几何约束、填充率、安匝数计算以及综合性能指标(如KV、铜损及其热负荷)。
        返回详细的评测指标字典，可直接输出至前端。
        """
        # 获取线规细节
        d_bare, d_od, area_bare = self.get_wire_properties(self.wp['awg'])
        
        # 计算纯几何槽面积
        slot_area_geo = self.calculate_slot_area()
        
        # 几何干涉推导槽内可装载最大导线数
        n_wires = self.estimate_max_wires_in_slot(d_od)
        
        # 满槽率判断: (导线裸铜面积 * 根数) / 槽几何面积
        cu_area_total = n_wires * area_bare
        fill_factor_cu = (cu_area_total / slot_area_geo) * 100
        
        # 计算当前电流设定下的热负荷: 安匝数 (Ampere-Turns)
        target_j = self.wp['rated_current_density']
        current_per_wire = target_j * area_bare
        total_ampere_turns = n_wires * current_per_wire
        
        # 定子内径换算为厘米，计算线负荷 Line Load (A/cm)
        stator_id = (self.gp['r_rotor_outer'] + self.gp['d_air_gap']) * 2
        stator_id_cm = stator_id / 10.0
        
        # 所有导线产生的总电枢安倍 / 定子内圆周长
        total_conductors = n_wires * self.wp['slot_count']
        line_load_a = (total_conductors * current_per_wire) / (math.pi * stator_id_cm)
        
        # AJ 值 (热负荷参数): 用于快速衡量发热能力，AJ = 线负荷(A/cm) * 电流密度(A/mm^2)
        aj_value = line_load_a * target_j
        
        # 结合以上的计算结果进行电气与反电动势求值
        bemf = self.calc_bemf_constants(n_wires)
        
        # 结合电阻率，转速，电流，计算各类损耗
        losses = self.calc_motor_losses(bemf["turns_per_phase"], current_per_wire, self.wp['rated_speed'])
        
        # 输出给UI展示的定子轭部厚度
        yoke_thickness = self.gp['r_stator_outer'] - (self.gp['r_rotor_outer'] + self.gp['d_air_gap'] + self.gp['d_tooth'])
        
        return {
            "awg": self.wp['awg'],                                     # 选定线规
            "wires_per_slot": n_wires,                                 # 单槽总导线根数
            "turns_per_tooth": n_wires / 2,                            # 每个齿(通常包含双层一边)上的匝数
            "slot_fill_factor_cu": round(fill_factor_cu, 2),           # 槽满率 (单纯铜占比)
            "current_per_wire_a": round(current_per_wire, 4),          # 单根导线相电流
            "ampere_turns_per_slot": round(total_ampere_turns, 2),     # 单槽安匝数
            "line_load_a_per_cm": round(line_load_a, 2),               # 线负荷指标
            "thermal_load_aj": round(aj_value, 0),                     # 综合发热载荷 AJ 值
            "yoke_thickness_mm": round(yoke_thickness, 3),             # 轭部厚度
            "KV_rpm_V": round(bemf["KV_rpm_V"], 0),                    # 预估KV值
            "phase_resistance_100C": round(losses["phase_resistance_100C"], 3), # 相电阻
            "total_loss_W": round(losses["total_loss_W"], 2),          # 总损耗
            "copper_loss_W": round(losses["copper_loss_W"], 2),        # 铜损耗
            "iron_loss_W": round(losses["iron_loss_W"], 2)             # 铁损耗预估
        }
