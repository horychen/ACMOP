import inner_rotor_motor
from collections import OrderedDict

# template 有点类似analytical的电机（由几何尺寸组成）
# variant则有点像是具体的电机实现类（由各个局部类，比如转子、定子等组成）

class vernier_motor_VShapePM_template(inner_rotor_motor.template_machine_as_numbers):
    def __init__(self, fea_config_dict=None, spec_input_dict=None):
        # 基本信息
        self.machine_type = 'PMVM'
        self.name = '__PMVM'

        # 仿真输入
        self.fea_config_dict = fea_config_dict
        self.spec_input_dict = spec_input_dict
        self.SD = SD = spec_input_dict # short name

        # 初始化搜索空间
        childGP = OrderedDict({
            # Vernier specific
            "deg_alpha_vspm"    : acmop_parameter("free",     "v-shape_magnet_tilt_angle",     None, [None, None]),
            "mm_d_bg_air"       : acmop_parameter("free",     "magnet_bridge_depth_to_air",    None, [None, None]),
            "mm_d_bg_magnet"    : acmop_parameter("free",     "magnet_bridge_depth_to_magnet", None, [None, None]),
        })
        self.d.update( {"GP": childGP} )

        # all in one place
        self.d = {
            "which_filter": fea_config_dict['which_filter'],
            "GP": geometric_parameters,
            "OP": OrderedDict()
        }

        if 'FixedSleeveLength' in self.d['which_filter']:
            self.d['GP']['mm_d_sleeve'].type = "fixed"
        elif 'VariableSleeveLength' in self.d['which_filter']:
            self.d['GP']['mm_d_sleeve'].type = "free"
        else:
            raise Exception('Not defined', self.d['which_filter'])

        # Get Analytical Design
        self.ModifiedBianchi2006(fea_config_dict, SD)

        # 定义搜索空间，determine bounds
        original_template_neighbor_bounds = self.get_template_neighbor_bounds()
        self.bounds_denorm = []
        for key, val in original_template_neighbor_bounds.items():
            parameter = self.d['GP'][key]
            parameter.bounds = val
            if parameter.type == 'free':
                self.bounds_denorm.append(val)
        print(f'[template] BOUNDS_denorm {len(self.bounds_denorm)}', self.bounds_denorm)

        # debug
        # from pprint import pprint
        # for key, val in self.d['GP'].items():
        #     if val.type=='free':
        #         print(key, '\t', val)
        # print('--------------------')
        # for key, val in self.d['GP'].items():
        #     if val.type!='free':
        #         print(key, '\t', val)
        # quit()

        # Template's Other Properties (Shared by the swarm)
        OP = self.d['OP']
        if True:
            # WINDING Layout
            OP['wily'] = wily = winding_layout.winding_layout_v2(SD['DPNV_or_SEPA'], SD['Qs'], SD['p'], SD['ps'], SD['coil_pitch_y'])
            # STACK LENGTH
            OP['mm_stack_length'] = pyrhonen_procedure_as_function.get_mm_stack_length(SD) # mm TODO:
            OP['mm_mechanical_air_gap_length'] = SD['minimum_mechanical_air_gap_length_mm']
            # THERMAL Properties
            OP['Js']                = SD['Js'] # Arms/mm^2 im_OP['Js'] 
            OP['fill_factor']       = SD['space_factor_kCu'] # im_OP['fill_factor'] 
            # MOTOR Winding Excitation Properties
            OP['DriveW_zQ']         =            pyrhonen_procedure_as_function.get_zQ(SD, self.d['GP']['mm_r_si'].value*2*1e-3, self.d['GP']['mm_r_or'].value*2*1e-3) # TODO:
            OP['DriveW_CurrentAmp'] = np.sqrt(2)*pyrhonen_procedure_as_function.get_stator_phase_current_rms(SD) # TODO:
            OP['DriveW_Freq']       = SD['ExcitationFreqSimulated']
            OP['DriveW_Rs']         = 1.0 # TODO: Must be greater than zero to let JMAG work
            OP['DriveW_poles']      = SD['p']*2
        # self.d.update( {"OP": OP} )

    def ModifiedBianchi2006(self, fea_config_dict, spec_input_dict):
        # ease of writing
        template = self
        GP = self.d['GP']
        OP = self.d['OP']
        SD = spec_input_dict

        air_gap_flux_density_B = SD['guess_air_gap_flux_density_B']
        # air_gap_flux_density_B = 0.9
        stator_tooth_flux_density_B_ds = SD['guess_stator_tooth_flux_density_B_ds']
        # stator_tooth_flux_density_B_ds = 1.5

        stator_yoke_flux_density_Bys = SD['guess_stator_yoke_flux_density_Bys']

        if SD['p'] >= 2:
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.75
            alpha_rm_over_alpha_rp = 1.0
            # stator_yoke_flux_density_Bys = 1.2
        else:
            # penalty for p=1 motor, i.e., large yoke height
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.5
            alpha_rm_over_alpha_rp = 0.75
            # stator_yoke_flux_density_Bys = 1.5

        stator_outer_diameter_Dse = 0.128 # m
        stator_outer_radius_r_os  = 0.5*stator_outer_diameter_Dse

        speed_rpm = SD['ExcitationFreqSimulated'] * 60 / SD['p'] # rpm

        split_ratio = 0.5 # refer to 2020-MLMS-0953@Fig. 5
        rotor_outer_radius_r_or = pyrhonen_procedure_as_function.eric_specify_tip_speed_get_radius(SD['tip_speed'], speed_rpm)
        rotor_outer_diameter_Dr = rotor_outer_radius_r_or*2
        sleeve_length = 0.5
        stator_inner_radius_r_is  = rotor_outer_radius_r_or + (sleeve_length+SD['minimum_mechanical_air_gap_length_mm'])*1e-3 # m (sleeve 3 mm, air gap 0.75 mm)
        stator_inner_diameter_Dis = stator_inner_radius_r_is*2
        split_ratio = stator_inner_diameter_Dis / stator_outer_diameter_Dse

        stator_yoke_height_h_ys = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis * alpha_rm_over_alpha_rp / (2*stator_yoke_flux_density_Bys * 2*SD['p'])
        stator_tooth_height_h_ds = (stator_outer_diameter_Dse - stator_inner_diameter_Dis) / 2 - stator_yoke_height_h_ys
        stator_slot_height_h_ss = stator_tooth_height_h_ds
        stator_tooth_width_b_ds = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis / (stator_tooth_flux_density_B_ds* SD['Qs'])

        OP['stator_slot_area'] = stator_slot_area = np.pi/(4*SD['Qs']) * ((stator_outer_diameter_Dse - 2*stator_yoke_height_h_ys)**2 - stator_inner_diameter_Dis**2) - stator_tooth_width_b_ds * stator_tooth_height_h_ds

        slot_pitch_pps = np.pi * (stator_inner_diameter_Dis + stator_slot_height_h_ss) / SD['Qs']
        kov = 1.8 # \in [1.6, 2.0]
        OP['end_winding_length_Lew'] = end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_b_ds) + slot_pitch_pps*kov * (SD['coil_pitch_y'] - 1)

        Q = SD['Qs']
        p = SD['p']
        # STATOR
        GP['deg_alpha_st'].value         = 360/Q - 2 # deg
        GP['deg_alpha_so'].value         = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
        GP['mm_r_si'].value              = 1e3*stator_inner_radius_r_is # mm
        GP['mm_r_os'].value              = 1e3*stator_outer_diameter_Dse/2 # mm
        GP['mm_d_so'].value              = 1 # mm
        GP['mm_d_sp'].value              = 1.5*GP['mm_d_so'].value
        GP['mm_d_st'].value              = 1e3*(0.5*stator_outer_diameter_Dse - stator_yoke_height_h_ys) - GP['mm_r_si'].value - GP['mm_d_sp'].value  # mm
        GP['mm_d_sy'].value              = 1e3*stator_yoke_height_h_ys # mm
        GP['mm_w_st'].value              = 1e3*stator_tooth_width_b_ds # mm
        # ROTOR
        GP['mm_d_sleeve'].value          = 2 # 3 # mm
        GP['mm_d_fixed_air_gap'].value   = SD['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_d_pm'].value              = 4  # mm
        GP['mm_d_ri'].value              = 1e3*ROTOR_STATOR_YOKE_HEIGHT_RATIO*stator_yoke_height_h_ys # TODO：This ratio (0.75) is randomly specified
        GP['mm_r_or'].value              = 1e3*rotor_outer_radius_r_or
        GP['mm_r_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_d_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_fixed_air_gap'].value
        # Vernier specific
        GP["deg_alpha_vspm"].value       = 20.3
        GP["mm_d_bg_air"].value          = 1.5
        GP["mm_d_bg_magnet"].value       = 1.5


    def get_template_neighbor_bounds(self):
        Q = self.SD['Qs']
        p = self.SD['p']
        pr = self.SD['pr']

        GP = self.d['GP']

        original_template_neighbor_bounds = {
            # STATOR
            "deg_alpha_st": [ 0.35*360/Q, 0.9*360/Q],
            "mm_d_so":      [  0.5,   5],                                                       
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.2*GP['mm_d_st'].value],                
            "mm_r_os":      [1.0*GP['mm_r_os'].value, 1.2*GP['mm_r_os'].value], 
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],                
            # ROTOR
            "mm_d_sleeve":  [0.4,   2], # this is air gap length for vernier (note min. mech. air gap length is 0.1mm)
            "split_ratio":  [0.2, 0.6], # Binder-2020-MLMS-0953@Fig.7
            "mm_d_pm":      [2.5, 7],
            "mm_d_ri":      [0.8*GP['mm_d_ri'].value,  1.2*GP['mm_d_ri'].value],                              
            # Vernier Specific
            "deg_alpha_vspm"    : [10, 40],
            "mm_d_bg_air"       : [1, 5],
            "mm_d_bg_magnet"    : [1, 5],
        }
        return original_template_neighbor_bounds

    def get_rotor_volume(self, stack_length=None):
        if stack_length is None:
            return np.pi*(self.d['GP']['mm_r_or'].value*1e-3)**2 * (self.d['OP']['mm_stack_length']*1e-3)
        else:
            return np.pi*(self.d['GP']['mm_r_or'].value*1e-3)**2 * (stack_length*1e-3)

    def get_rotor_weight(self, gravity=9.8, stack_length=None):
        material_density_rho = pyrhonen_procedure_as_function.get_material_data()[0]
        if stack_length is None:
            return gravity * self.get_rotor_volume() * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg
        else:
            return gravity * self.get_rotor_volume(stack_length=stack_length) * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg

    def build_x_denorm(self):
        # design_parameters = build_design_parameters_list(GP, SD) # those member variables are defined in Pyrhonen's procedure

        GP = self.d['GP']

        self.x_denorm_dict = self.get_x_denorm_dict_from_geometric_parameters(GP)

        x_denorm = [val for key, val in self.x_denorm_dict.items()]
        if False:
            from pprint import pprint
            pprint(x_denorm_dict)
            pprint(x_denorm)
            pprint(GP)
            quit()
        return x_denorm

    def build_design_parameters_list(self):
        GP = self.d['GP']
        SD = self.SD
        # obsolete feature
        design_parameters = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0, #'mm_r_st',
            0.0, #'mm_r_sf',
            0.0, #'mm_r_sb',
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
            ]
        return design_parameters

    def get_x_denorm_dict_from_geometric_parameters(self, GP):
        x_denorm_dict = OrderedDict()
        for key, parameter in GP.items():
            if parameter.type == 'free':
                x_denorm_dict[key] = parameter.value
        return x_denorm_dict
    def get_x_denorm_dict_from_x_denorm_list(self, x_denorm):
        # 先拿个模板来，但是几何尺寸的变量值是旧的
        x_denorm_dict = self.get_x_denorm_dict_from_geometric_parameters(self.d['GP'])

        # 对模板进行遍历，挨个把新的几何尺寸的值从x_denorm中读取出来并更新x_denorm_dict
        for key, new_val in zip(x_denorm_dict.keys(), x_denorm):
            x_denorm_dict[key] = new_val
        return x_denorm_dict

    def update_geometric_parametes_using_x_denorm_dict(self, x_denorm_dict):
        for key, val in x_denorm_dict.items():
            self.d['GP'][key].value = val
        return self.d['GP']

