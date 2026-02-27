from collections import OrderedDict

class LegacyCodes:
    def PracticalInitialDesign(self, fea_config_dict, SI, GP, EX):

        # 定子外径 # this is related to the stator current density and should be determined by Js and power.
        stator_outer_diameter_Dso = GP['mm_r_so'].value*2 * 1e-3 # [m]
        # 定子内径
        # 定子内半径
        stator_inner_diameter_Dsi = stator_outer_diameter_Dso * GP['split_ratio'].value # [m]
        stator_inner_radius_r_is = stator_inner_diameter_Dsi * 0.5
        # 机械气隙长度
        # 护套长度（相当于气隙）
        # 总等效气隙长度（不导磁也不是永磁体的部分都算作气隙；如果把永磁体作为气隙怪怪的，我们不要这么做）
        mech_air_gap_length = GP['mm_d_mech_air_gap'].value *1e-3 # [m]
        sleeve_length = GP['mm_d_sleeve'].value * 1e-3 # [m]
        equivalent_air_gap_length = sleeve_length + mech_air_gap_length
        # 包含永磁体（如果有）在内的转子外径
        rotor_outer_radius_r_or = stator_inner_radius_r_is - equivalent_air_gap_length

        # 基于磁负荷去分配定子齿和轭的尺寸
        Bg = SI['guess_air_gap_flux_density_Bg'] # 0.9 T
        Bst = SI['guess_stator_tooth_flux_density_Bst'] # 1.5 T
        Bsy = SI['guess_stator_yoke_flux_density_Bsy'] # 1.2 T

        def get_alpha_rm_over_alpha_rp(p):
            if SI['p'] >= 2:
                ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.75
                alpha_rm_over_alpha_rp = 1.0
                # stator_yoke_flux_density_Bsy = 1.2
            else:
                # penalty for p=1 motor, i.e., large yoke height
                ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.5
                alpha_rm_over_alpha_rp = 0.75
                # stator_yoke_flux_density_Bsy = 1.5
            return alpha_rm_over_alpha_rp

        Q = SI['Qs']
        p = SI['p']

        # 单个永磁体极面下，气隙中的磁通全部进入两侧的定子轭部（的深度）
        # 定子除了轭部以外的部分全部为齿部
        stator_yoke_depth_d_sy  = Bg * math.pi * stator_inner_diameter_Dsi * get_alpha_rm_over_alpha_rp(p) / (2*Bsy * 2*p)

        # 定子齿部深度依赖于轭部深度
        stator_tooth_depth_d_st = (stator_outer_diameter_Dso - stator_inner_diameter_Dsi) *0.5 - stator_yoke_depth_d_sy
        stator_slot_depth_d_ss = stator_tooth_depth_d_st

        # 单个永磁体极面下，气隙中的磁通全部进入Qs个定子齿部（的宽度）
        stator_tooth_width_w_st = Bg * math.pi * stator_inner_diameter_Dsi / (Bst* Q)

        # 计算槽面积
        def get_stator_slot_area(Q, stator_outer_diameter_Dso, stator_yoke_depth_d_sy, stator_inner_diameter_Dsi, stator_tooth_width_w_st, stator_tooth_depth_d_st):
            return math.pi/(4*Q) * ((stator_outer_diameter_Dso - 2*stator_yoke_depth_d_sy)**2 - stator_inner_diameter_Dsi**2) - stator_tooth_width_w_st * stator_tooth_depth_d_st
        EX['stator_slot_area'] = get_stator_slot_area(Q, stator_outer_diameter_Dso, stator_yoke_depth_d_sy, stator_inner_diameter_Dsi, stator_tooth_width_w_st, stator_tooth_depth_d_st)

        # 计算端部绕组长度
        slot_pitch_pps = math.pi * (stator_inner_diameter_Dsi + stator_slot_depth_d_ss) / Q
        kov = EX['end_winding_length_factor_kov']
        EX['end_winding_length_Lew'] = math.pi*0.5 * (slot_pitch_pps + stator_tooth_width_w_st) + slot_pitch_pps*kov * (SI['coil_pitch_y'] - 1)

        # STATOR
        GP['mm_r_si'].value              = 1e3*stator_inner_radius_r_is # mm
        GP['mm_r_so'].value              = 1e3*stator_outer_diameter_Dso*0.5 # mm
        GP['mm_d_st'].value              = 1e3*stator_tooth_depth_d_st # mm
        GP['mm_d_sy'].value              = 1e3*stator_yoke_depth_d_sy # mm
        GP['mm_w_st'].value              = 1e3*stator_tooth_width_w_st # mm
        # ROTOR
        GP['mm_d_sleeve'].value          = 1e3*sleeve_length
        GP['mm_d_mech_air_gap'].value    = 1e3*mech_air_gap_length
        GP['mm_r_ro'].value              = 1e3*rotor_outer_radius_r_or
        GP['mm_d_ri'].value              = 1e3*(stator_inner_radius_r_is - equivalent_air_gap_length) - GP['mm_d_pm'].value - GP['mm_r_ri'].value

