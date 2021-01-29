import CrossSectInnerNotchedRotor
import CrossSectStator
import Location2D
import pyrhonen_procedure_as_function, winding_layout
import numpy as np
import logging
import utility
from time import time as clock_time
from collections import OrderedDict
from recordtype import recordtype
EPS = 1e-2 # unit: mm

# template 有点类似analytical的电机（由几何尺寸组成）
# variant则有点像是具体的电机实现类（由各个局部类，比如转子、定子等组成）

acmop_parameter = recordtype('acmop_parameter', 'type, name, value, bounds')
class bearingless_spmsm_template(object):
    def __init__(self, fea_config_dict=None, spec_input_dict=None):
        # 基本信息
        self.machine_type = 'SPMSM'
        self.name = '__SPMSM'

        # 仿真输入
        self.fea_config_dict = fea_config_dict
        self.spec_input_dict = spec_input_dict
        self.SD = SD = spec_input_dict # short name

        # 初始化搜索空间
        self.bound_filter = None
        geometric_parameters = OrderedDict({
            # STATOR
            "deg_alpha_st" : acmop_parameter("free", "stator_tooth_span_angle"       , None, [None, None]),
            "mm_d_so"      : acmop_parameter("free", "stator_tooth_open_depth"       , None, [None, None]),
            "mm_d_st"      : acmop_parameter("free", "stator_tooth_depth"            , None, [None, None]),
            "mm_r_os"      : acmop_parameter("free", "outer_stator_radius"           , None, [None, None]),
            "mm_w_st"      : acmop_parameter("free", "stator_tooth_width"            , None, [None, None]),
            "mm_r_si"      : acmop_parameter("derived", "inner_stator_radius"        , None, [None, None]),
            "deg_alpha_so" : acmop_parameter("derived", "stator_tooth_open_angle"    , None, [None, None]),
            "mm_d_sp"      : acmop_parameter("derived", "stator_tooth_tip_depth"     , None, [None, None]),
            "mm_d_sy"      : acmop_parameter("derived", "stator_yoke_depth"          , None, [None, None]),
            # ROTOR
            "mm_d_fixed_air_gap": acmop_parameter("fixed",                                                "sleeve_length",                 None, [None, None]),
            "mm_d_sleeve"       : acmop_parameter("fixed",                                                "sleeve_length",                 None, [None, None]),
            "split_ratio"       : acmop_parameter("free",                                                 "split_ratio_r_is_slash_r_os",   None, [None, None]),
            "mm_d_pm"           : acmop_parameter("free",                                                 "magnet_depth",                  None, [None, None]),
            "deg_alpha_rm"      : acmop_parameter("free",                                                 "magnet_pole_span_angle",        None, [None, None]),
            "mm_d_ri"           : acmop_parameter("free",                                                 "rotor_iron (back iron) depth",  None, [None, None]),
            "mm_d_rp"           : acmop_parameter("free",                                                 "inter_polar_iron_thickness",    None, [None, None]),
            "mm_r_or"           : acmop_parameter("derived",                                              "outer_rotor_radius",            None, [None, None]),
            "mm_r_ri"           : acmop_parameter("derived",                                              "inner_rotor_radius",            None, [None, None]),
            "deg_alpha_rs"      : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "magnet_segment_span_angle",     None, [None, None]),
            "mm_d_rs"           : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "inter_segment_iron_thickness",  None, [None, None]),
        })
        # all in one place
        self.d = {
            "which_filter": fea_config_dict['which_filter'],
            "GP": geometric_parameters,
            "OP": None
        }

        if 'FixedSleeveLength' in self.d['which_filter']:
            self.d['GP']['mm_d_sleeve'].type = "fixed"
        elif 'VariableSleeveLength' in self.d['which_filter']:
            self.d['GP']['mm_d_sleeve'].type = "free"
        else:
            raise Exception('Not defined', self.d['which_filter'])

        # Get Analytical Design
        self.Bianchi2006(fea_config_dict, SD)

        # 定义搜索空间，determine bounds
        original_template_neighbor_bounds = self.get_template_neighbor_bounds()
        self.bounds_denorm = []
        for key, val in original_template_neighbor_bounds.items():
            self.d['GP'][key].bounds = val
            self.bounds_denorm.append(val)
        print('BOUNDS_denorm', self.bounds_denorm)

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
        OP = OrderedDict()
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
            # BEARING Winding Excitation Properties
            OP['BeariW_zQ']         = OP['DriveW_zQ']
            OP['BeariW_CurrentAmp'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] * (OP['DriveW_CurrentAmp'] / fea_config_dict['TORQUE_CURRENT_RATIO'])
            OP['BeariW_Freq']       = OP['DriveW_Freq']
            OP['BeariW_Rs']         = OP['DriveW_Rs'] * OP['BeariW_zQ'] / OP['DriveW_zQ']
            OP['BeariW_poles']      = SD['ps']*2
            OP['slot_current_utilizing_ratio'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] + fea_config_dict['TORQUE_CURRENT_RATIO'] # will be less than 1 for separate winding
        self.d.update( {"OP": OP} )

    def Bianchi2006(self, fea_config_dict, spec_input_dict):
        # ease of writing
        template = self
        GP = self.d['GP']
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

        stator_outer_diameter_Dse = 0.225 # m

        speed_rpm = SD['ExcitationFreqSimulated'] * 60 / SD['p'] # rpm

        rotor_outer_radius_r_or = pyrhonen_procedure_as_function.eric_specify_tip_speed_get_radius(SD['tip_speed'], speed_rpm)
        rotor_outer_diameter_Dr = rotor_outer_radius_r_or*2
        sleeve_length = 3
        stator_inner_radius_r_is  = rotor_outer_radius_r_or + (sleeve_length+SD['minimum_mechanical_air_gap_length_mm'])*1e-3 # m (sleeve 3 mm, air gap 0.75 mm)
        stator_inner_diameter_Dis = stator_inner_radius_r_is*2
        split_ratio = stator_inner_diameter_Dis / stator_outer_diameter_Dse

        stator_yoke_height_h_ys = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis * alpha_rm_over_alpha_rp / (2*stator_yoke_flux_density_Bys * 2*SD['p'])
        stator_tooth_height_h_ds = (stator_outer_diameter_Dse - stator_inner_diameter_Dis) / 2 - stator_yoke_height_h_ys
        stator_slot_height_h_ss = stator_tooth_height_h_ds
        stator_tooth_width_b_ds = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis / (stator_tooth_flux_density_B_ds* SD['Qs'])

        stator_slot_area = np.pi/(4*SD['Qs']) * ((stator_outer_diameter_Dse - 2*stator_yoke_height_h_ys)**2 - stator_inner_diameter_Dis**2) - stator_tooth_width_b_ds * stator_tooth_height_h_ds

        slot_pitch_pps = np.pi * (stator_inner_diameter_Dis + stator_slot_height_h_ss) / SD['Qs']
        kov = 1.8 # \in [1.6, 2.0]
        end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_b_ds) + slot_pitch_pps*kov * (SD['coil_pitch_y'] - 1)

        Q = SD['Qs']
        p = SD['p']
        GP['deg_alpha_st'].value         = 360/Q - 2 # deg
        GP['deg_alpha_so'].value         = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
        GP['mm_r_si'].value              = 1e3*stator_inner_radius_r_is # mm
        GP['mm_r_os'].value              = 1e3*stator_outer_diameter_Dse/2 # mm
        GP['mm_d_so'].value              = 1 # mm
        GP['mm_d_sp'].value              = 1.5*GP['mm_d_so'].value
        GP['mm_d_st'].value              = 1e3*(0.5*stator_outer_diameter_Dse - stator_yoke_height_h_ys) - GP['mm_r_si'].value - GP['mm_d_sp'].value  # mm
        GP['mm_d_sy'].value              = 1e3*stator_yoke_height_h_ys # mm
        GP['mm_w_st'].value              = 1e3*stator_tooth_width_b_ds # mm
        GP['mm_d_sleeve'].value          = 2 # 3 # mm
        GP['mm_d_fixed_air_gap'].value   = SD['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_r_or'].value              = 1e3*rotor_outer_radius_r_or
        GP['mm_d_pm'].value              = 4  # mm
        GP['deg_alpha_rm'].value         = 0.95*360/(2*p) # deg
        GP['deg_alpha_rs'].value         = 0.975*GP['deg_alpha_rm'].value / SD['no_segmented_magnets']
        GP['mm_d_ri'].value              = 1e3*ROTOR_STATOR_YOKE_HEIGHT_RATIO*stator_yoke_height_h_ys # TODO：This ratio (0.75) is randomly specified
        GP['mm_r_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_d_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_fixed_air_gap'].value
        GP['mm_d_rp'].value              = 3  # mm
        GP['mm_d_rs'].value              = 0.20*GP['mm_d_rp'].value # d_pm > d_rp and d_pm > d_rs

        # Those are some obsolete variables that are convenient to have.
        # template.Radius_OuterStatorYoke = spec_geometry_dict['Radius_OuterStatorYoke'] = 1e3*0.5*stator_outer_diameter_Dse # mm
        # template.Radius_OuterRotor      = spec_geometry_dict['Radius_OuterRotor'] = 1e3*rotor_outer_radius_r_or # mm

        # Those variables are for PMSM convenience
        # template.rotor_steel_outer_radius = spec_geometry_dict['rotor_steel_outer_radius'] = template.Radius_OuterRotor

        # required_torque = SD['mec_power']/(2*np.pi*speed_rpm)*60
        # rotor_volume_Vr = required_torque/(2*SD['TangentialStress'])
        # template.required_torque = required_torque

    def get_template_neighbor_bounds(self):
        Q = self.SD['Qs']
        p = self.SD['p']
        s = self.SD['no_segmented_magnets']

        GP = self.d['GP']

        original_template_neighbor_bounds = {
            "deg_alpha_st": [ 0.35*360/Q, 0.9*360/Q],
            "mm_d_so":      [  0.5,   5],                                                       
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.2*GP['mm_d_st'].value],                
            "mm_r_os":      [1.0*GP['mm_r_os'].value, 1.2*GP['mm_r_os'].value], 
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],                
            "mm_d_sleeve":  [3,   6],
            "mm_d_pm":      [2.5, 7],
            "deg_alpha_rm": [0.6*360/(2*p),          1.0*360/(2*p)],                                     
            "deg_alpha_rs": [0.8*360/(2*p)/s,        0.975*360/(2*p)/s],                               
            "mm_d_ri":      [0.8*GP['mm_d_ri'].value,  1.2*GP['mm_d_ri'].value],                              
            "split_ratio":  [0.2, 0.6], # Binder-2020-MLMS-0953@Fig.7
            "mm_d_rp":      [2.5,   6],                                                         
            "mm_d_rs":      [2.5,   6]
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
            GP['deg_alpha_st'].value,
            GP['deg_alpha_so'].value,
            GP['mm_r_si'].value,
            GP['mm_d_so'].value,
            GP['mm_d_sp'].value,
            GP['mm_d_st'].value,
            GP['mm_d_sy'].value,
            GP['mm_w_st'].value,
            0.0, #'mm_r_st',
            0.0, #'mm_r_sf',
            0.0, #'mm_r_sb',
            SD['Qs'],
            GP['mm_d_sleeve'].value,
            GP['mm_d_fixed_air_gap'].value,
            GP['mm_d_pm'].value,
            GP['deg_alpha_rm'].value,
            GP['deg_alpha_rs'].value,
            GP['mm_d_ri'].value,
            GP['mm_r_ri'].value,
            GP['mm_d_rp'].value,
            GP['mm_d_rs'].value,
            SD['p'],
            SD['no_segmented_magnets']
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

class bearingless_spmsm_design_variant(object):

    def __init__(self, spmsm_template=None, x_denorm=None, counter=None, counter_loop=None):

        self.template = spmsm_template

        #00 Settings
        # self.template.fea_config_dict = spmsm_template.fea_config_dict
        # self.template.spec_input_dict = spmsm_template.spec_input_dict
        # self.spec_geometry_dict = spmsm_template.spec_geometry_dict

        #01 Model ID
        # self.model_name_prefix
        self.counter = counter
        self.counter_loop = counter_loop
        if counter is not None:
            if counter_loop == 1:
                self.name = f"p{self.template.spec_input_dict['p']}ps{self.template.spec_input_dict['ps']}-Q{self.template.spec_input_dict['Qs']}y{self.template.spec_input_dict['coil_pitch_y']}-{counter:04d}"
            else:
                self.name = f"p{self.template.spec_input_dict['p']}ps{self.template.spec_input_dict['ps']}-Q{self.template.spec_input_dict['Qs']}y{self.template.spec_input_dict['coil_pitch_y']}-{counter:04d}-redo{counter_loop}"
        else:
            self.name = 'SPMSM_InitialDesign'
        self.ID = 'Q%dp%ds%d'%(self.template.SD['Qs'], self.template.SD['p'],self.template.SD['no_segmented_magnets'])

        #02 Geometry Data
        if x_denorm is None:
            GP = self.template.d['GP'] # do nothing, use template's GP
        else:
            x_denorm_dict = self.template.get_x_denorm_dict(x_denorm)
            GP = self.template.update_geometric_parametes_using_x_denorm_dict(x_denorm_dict)


        #test: yes, it is by reference!
        # print(self.template.d['GP'])
        # GP['aaaaaa'] = 1239129
        # print(self.template.d['GP'])
        # quit()

        # 不合理的变量选择（mm_d_rp）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        if GP['mm_d_rp'].value > GP['mm_d_pm'].value:
            GP['mm_d_rp'].value            = GP['mm_d_pm'].value
            free_variables[11] = free_variables[6]

            msg = '[Warning from bearingless_spmsm_design.py]: Inter-pole notch depth mm_d_rp cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rp is set to mm_d_pm.'
            print(msg)
            logger = logging.getLogger(__name__).warn(msg)

        # 不合理的变量选择（mm_d_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        if GP['mm_d_rs'].value > GP['mm_d_pm'].value:
            GP['mm_d_rs'].value            = GP['mm_d_pm'].value
            free_variables[12] = free_variables[6]

            msg = '[Warning from bearingless_spmsm_design.py]: Inter-segment notch depth mm_d_rs cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rs is set to mm_d_pm.'
            print(msg)
            logger = logging.getLogger(__name__).warn(msg)

        # 不合理的变量选择（deg_alpha_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        if not (GP['deg_alpha_rs'].value > GP['deg_alpha_rm'].value/spmsm_template.SD['no_segmented_magnets']):
            GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value/spmsm_template.SD['no_segmented_magnets']
            msg = '[Warning from bearingless_spmsm_design.py]: deg_alpha_rs cannot be larger than deg_alpha_rm/s deg_alpha_rs is set to deg_alpha_rm/s.'
            print(msg)
            logger = logging.getLogger(__name__).warn(msg)

        # 如果没有永磁体分段，那么alpha_rs应该等于alpha_rm。
        if spmsm_template.SD['no_segmented_magnets'] == 1:
            GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value # raise Exception('Invalid alpha_rs. Check that it is equal to alpha_rm for s=1')
            GP['mm_d_rs'].value = 0 # raise Exception('Invalid d_rs. Check that it is equal to 0 for s =1')
        # This is all we need



        # Parts
        self.rotorCore = CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                            name = 'NotchedRotor',
                            mm_d_pm      = GP['mm_d_pm'].value,
                            deg_alpha_rm = GP['deg_alpha_rm'].value, # angular span of the pole: class type DimAngular
                            deg_alpha_rs = GP['deg_alpha_rs'].value, # segment span: class type DimAngular
                            mm_d_ri      = GP['mm_d_ri'].value, # rotor iron thickness: class type DimLinear
                            mm_r_ri      = GP['mm_r_ri'].value, # inner radius of rotor: class type DimLinear
                            mm_d_rp      = GP['mm_d_rp'].value, # interpolar iron thickness: class type DimLinear
                            mm_d_rs      = GP['mm_d_rs'].value, # inter segment iron thickness: class type DimLinear
                            p = spmsm_template.SD['p'], # Set pole-pairs to 2
                            s = spmsm_template.SD['no_segmented_magnets'], # Set magnet segments/pole to 4
                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0))

        self.shaft = CrossSectInnerNotchedRotor.CrossSectShaft(name = 'Shaft',
                                                      notched_rotor = self.rotorCore
                                                    )

        self.rotorMagnet = CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
                                                      notched_rotor = self.rotorCore
                                                    )

        self.stator_core = CrossSectStator.CrossSectInnerRotorStator( name = 'StatorCore',
                                            deg_alpha_st = GP['deg_alpha_st'].value, #40,
                                            deg_alpha_so = GP['deg_alpha_so'].value, #20,
                                            mm_r_si = GP['mm_r_si'].value,
                                            mm_d_so = GP['mm_d_so'].value,
                                            mm_d_sp = GP['mm_d_sp'].value,
                                            mm_d_st = GP['mm_d_st'].value,
                                            mm_d_sy = GP['mm_d_sy'].value,
                                            mm_w_st = GP['mm_w_st'].value,
                                            mm_r_st = 0.0, # =0
                                            mm_r_sf = 0.0, # =0
                                            mm_r_sb = 0.0, # =0
                                            Q = spmsm_template.SD['Qs'],
                                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                            )

        self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                    stator_core = self.stator_core)

        self.sleeve = CrossSectInnerNotchedRotor.CrossSectSleeve(
                            name = 'Sleeve',
                            notched_magnet = self.rotorMagnet,
                            d_sleeve = GP['mm_d_sleeve'].value
                            )

        #03 Inherit properties
        self.template.d['OP']
            # self.DriveW_CurrentAmp = None ########### will be assisned when drawing the coils
            # self.DriveW_poles      = spmsm_template.DriveW_poles

            # self.Js                = spmsm_template.Js 
            # self.fill_factor       = spmsm_template.fill_factor 

            # self.template.d['OP']['mm_stack_length']      = spmsm_template.stack_length 

        #03 Mechanical Parameters
        self.update_mechanical_parameters()

        # #04 Material Condutivity Properties
        # if self.template.fea_config_dict is not None:
        #     self.End_Ring_Resistance = fea_config_dict['End_Ring_Resistance']
        #     self.Bar_Conductivity = fea_config_dict['Bar_Conductivity']
        # self.Copper_Loss = self.DriveW_CurrentAmp**2 / 2 * self.DriveW_Rs * 3
        # # self.Resistance_per_Turn = 0.01 # TODO

        #05 Winidng Excitation
        self.CurrentAmp_per_phase = None # will be used in copper loss calculation

        #06 Meshing & Solver Properties
        self.max_nonlinear_iteration = 50 # 30 for transient solve
        self.meshSize_Magnet = 2 # mm

    def update_mechanical_parameters(self, syn_freq=None):
        OP = self.template.d['OP']
        if syn_freq is None:
            OP['the_speed'] = OP['DriveW_Freq']*60. / (0.5*OP['DriveW_poles']) # rpm
            OP['Omega']     = OP['the_speed'] / 60. * 2*np.pi
            # self.omega = None # This variable name is devil! you can't tell its electrical or mechanical! #+ self.DriveW_Freq * (1-self.the_slip) * 2*pi
        else:
            raise Exception('Not implemented.')

    def draw_spmsm(self, toolJd, bool_pyx=False):

        # blue
        # color_rgb_A = np.array([113, 142, 164])/255
        # color_rgb_B = np.array([73, 109, 137])/255

        # yellow
        # color_rgb_A = np.array([255, 252, 170])/255
        # color_rgb_B = np.array([212, 208, 166])/255

        # gray
        color_rgb_A = np.array([236,236,236])/255
        color_rgb_B = np.array([226,226,226])/255

        # Rotor Core
        list_regions_1 = self.rotorCore.draw(toolJd)
        toolJd.bMirror = False
        toolJd.iRotateCopy = self.rotorCore.p*2
        region1 = toolJd.prepareSection(list_regions_1, color=color_rgb_A)

        # Shaft
        if not bool_pyx:
            list_regions = self.shaft.draw(toolJd)
            toolJd.bMirror = False
            toolJd.iRotateCopy = 1
            region0 = toolJd.prepareSection(list_regions)

        # Rotor Magnet
        list_regions = self.rotorMagnet.draw(toolJd)
        toolJd.bMirror = False
        toolJd.iRotateCopy = self.rotorMagnet.notched_rotor.p*2
        region2 = toolJd.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

        # This is only for post-processing and it is for handle a un-fixable filling bug with PyX.
        if bool_pyx:
            region1 = toolJd.prepareSection(list_regions_1, color=color_rgb_A)


        # Sleeve
        if not bool_pyx:
            list_regions = self.sleeve.draw(toolJd)
            toolJd.bMirror = False
            toolJd.iRotateCopy = self.rotorMagnet.notched_rotor.p*2
            regionS = toolJd.prepareSection(list_regions)

        # Stator Core
        list_regions = self.stator_core.draw(toolJd)
        toolJd.bMirror = True
        toolJd.iRotateCopy = self.stator_core.Q
        region3 = toolJd.prepareSection(list_regions, color=color_rgb_A)

        if not bool_pyx:
            # Stator Winding
            list_regions = self.coils.draw(toolJd)
            toolJd.bMirror = False
            toolJd.iRotateCopy = self.coils.stator_core.Q
            region4 = toolJd.prepareSection(list_regions)

            # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
            OP = self.template.d['OP']
            CurrentAmp_in_the_slot = self.coils.mm2_slot_area * OP['fill_factor'] * OP['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / OP['DriveW_zQ']
            CurrentAmp_per_phase = CurrentAmp_per_conductor * OP['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
            # Maybe there is a bug here... regarding the excitation for suspension winding...
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
            OP['CurrentAmp_per_phase'] = CurrentAmp_per_phase
            OP['DriveW_CurrentAmp'] = self.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            OP['BeariW_CurrentAmp'] = self.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            # self.spec_geometry_dict['DriveW_CurrentAmp'] = self.DriveW_CurrentAmp

            slot_current_utilizing_ratio = (OP['DriveW_CurrentAmp'] + OP['BeariW_CurrentAmp']) / OP['CurrentAmp_per_phase']
            print('---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio)

            # print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
            # print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
            # print('---self.DriveW_CurrentAmp =', self.DriveW_CurrentAmp)
            # print('---self.BeariW_CurrentAmp =', self.BeariW_CurrentAmp)
            # print('---TORQUE_CURRENT_RATIO:', self.template.fea_config_dict['TORQUE_CURRENT_RATIO'])
            # print('---SUSPENSION_CURRENT_RATIO:', self.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            # Import Model into Designer
            toolJd.save(self.name, self.show(toString=True))

        return True

    def pre_process(self, app, model):
        # pre-process : you can select part by coordinate!
        ''' Group '''
        def group(name, id_list):
            model.GetGroupList().CreateGroup(name)
            for the_id in id_list:
                model.GetGroupList().AddPartToGroup(name, the_id)
                # model.GetGroupList().AddPartToGroup(name, name) #<- this also works

        part_ID_list = model.GetPartIDs()
        # print(part_ID_list)
        # quit()

        # view = app.View()
        # view.ClearSelect()
        # sel = view.GetCurrentSelection()
        # sel.SelectPart(123)
        # sel.SetBlockUpdateView(False)
                                #  轴 转子    永磁体         套 定子      绕组
        SD = self.template.SD
        p = SD['p']
        s = SD['no_segmented_magnets']
        Q = SD['Qs']
        if len(part_ID_list) != int(1 + 1 + p*2*s + 1 + 1 + Q*2):
            msg = 'Number of Parts is unexpected. Should be %d but get %d.\n'%(int(1 + 1 + p*2*s + 1 + 1 + Q*2), len(part_ID_list)) + self.show(toString=True)
            logger = logging.getLogger(__name__)
            logger.error(msg)
            raise utility.ExceptionBadNumberOfParts(msg)

        self.id_backiron = id_backiron = part_ID_list[0]
        id_shaft = part_ID_list[1]
        partIDRange_Magnet = part_ID_list[2:int(2+p*s*2)]
        id_sleeve = part_ID_list[int(2+p*s*2)]
        id_statorCore = part_ID_list[int(2+p*s*2)+1]
        partIDRange_Coil = part_ID_list[int(2+p*s*2)+2 : int(2+p*s*2)+2 + int(Q*2)]

        # debug
        # print(id_backiron)
        # print(id_shaft)
        # print(partIDRange_Magnet)
        # print(id_sleeve)
        # print(id_statorCore)
        # print(partIDRange_Coil)

        model.SuppressPart(id_sleeve, 1)

        group("Magnet", partIDRange_Magnet)
        group("Coils", partIDRange_Coil)

        ''' Add Part to Set for later references '''
        def add_part_to_set(name, x, y, ID=None):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection()
            if ID is None:
                # print x,y
                sel.SelectPartByPosition(x,y,0) # z=0 for 2D
            else:
                sel.SelectPart(ID)
            model.GetSetList().GetSet(name).AddSelected(sel)

        # def edge_set(name,x,y):
        #     model.GetSetList().CreateEdgeSet(name)
        #     model.GetSetList().GetSet(name).SetMatcherType(u"Selection")
        #     model.GetSetList().GetSet(name).ClearParts()
        #     sel = model.GetSetList().GetSet(name).GetSelection()
        #     sel.SelectEdgeByPosition(x,y,0) # sel.SelectEdge(741)
        #     model.GetSetList().GetSet(name).AddSelected(sel)
        # edge_set(u"AirGapCoast", 0, self.template.d['GP']['mm_r_or'].value+0.5*self.Length_AirGap)

        # Shaft
        add_part_to_set('ShaftSet', 0.0, 0.0, ID=id_shaft) # 坐标没用，不知道为什么，而且都给了浮点数了

        # Create Set for 4 poles Winding
        Angle_StatorSlotSpan = 360/Q
        # R = self.mm_r_si + self.mm_d_sp + self.mm_d_st *0.5 # this is not generally working (JMAG selects stator core instead.)
        # THETA = 0.25*(Angle_StatorSlotSpan)/180.*np.pi
        R = np.sqrt(self.coils.PCoil[0]**2 + self.coils.PCoil[1]**2)
        THETA = np.arctan(self.coils.PCoil[1]/self.coils.PCoil[0])
        X = R*np.cos(THETA)
        Y = R*np.sin(THETA)
        countXL = 0
        wily = self.template.d['OP']['wily']
        for UVW, UpDown in zip(wily.layer_X_phases,wily.layer_X_signs):
            countXL += 1 
            add_part_to_set("CoilLX%s%s %d"%(UVW,UpDown,countXL), X, Y)

            # print(X, Y, THETA)
            THETA += Angle_StatorSlotSpan/180.*np.pi
            X = R*np.cos(THETA)
            Y = R*np.sin(THETA)

        # Create Set for 2 poles Winding
        # THETA = 0.75*(Angle_StatorSlotSpan)/180.*np.pi # 这里这个角度的选择，决定了悬浮绕组产生悬浮力的方向！！！！！
        THETA = np.arctan(-self.coils.PCoil[1]/self.coils.PCoil[0]) + (2*np.pi)/Q
        X = R*np.cos(THETA)
        Y = R*np.sin(THETA)
        countYL = 0
        for UVW, UpDown in zip(wily.layer_Y_phases,wily.layer_Y_signs):
            countYL += 1 
            add_part_to_set("CoilLY%s%s %d"%(UVW,UpDown,countYL), X, Y)

            THETA += Angle_StatorSlotSpan/180.*np.pi
            X = R*np.cos(THETA)
            Y = R*np.sin(THETA)

        # Create Set for Magnets
        GP = self.template.d['GP']
        R = GP['mm_r_si'].value - GP['mm_d_sleeve'].value - GP['mm_d_fixed_air_gap'].value - 0.5*GP['mm_d_pm'].value
        alpha_rs = GP['deg_alpha_rs'].value /180*np.pi
        deg_pole_span = 360 / (p*2)

        if s>1:
            deg_alpha_notch  = (GP['deg_alpha_rm'].value - s*GP['deg_alpha_rs'].value) / (s-1) # inter-segment notch占的角度
            alpha_notch = deg_alpha_notch /180*np.pi

        list_xy_magnets = []
        # list_xy_airWithinRotorSlot = []
        for ind in range(int(p*2)):
            natural_ind = ind + 1

            if s==1:
                      # v---This negative sign means we walk CCW to assign sets.
                THETA = - (180/p-GP['deg_alpha_rm'].value + 0.5*GP['deg_alpha_rm'].value + deg_pole_span*ind) /180.*np.pi
                X = R*np.cos(THETA)
                Y = R*np.sin(THETA)

                add_part_to_set("Magnet %d"%(natural_ind), X, Y)
                list_xy_magnets.append([X,Y])
            else:     # v---This negative sign means we walk CCW to assign sets.
                THETA = - ( 180/p-GP['deg_alpha_rm'].value + 0.5*GP['deg_alpha_rs'].value + deg_pole_span*ind ) /180*np.pi # initial position
                # THETA = ( 0.5*self.deg_alpha_rs + deg_pole_span*ind ) /180*np.pi # initial position
                for s in range(s):
                    X = R*np.cos(THETA)
                    Y = R*np.sin(THETA)
                    add_part_to_set("Magnet %d s%d"%(natural_ind, s), X, Y)
                    list_xy_magnets.append([X,Y])
                    THETA -= alpha_notch + alpha_rs
                        # ^---This negative sign means we walk CCW to assign sets.

        # Create Set for Motion Region
        def part_list_set(name, list_xy, list_part_id=None, prefix=None):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection() 
            for xy in list_xy:
                sel.SelectPartByPosition(xy[0],xy[1],0) # z=0 for 2D
            if list_part_id is not None:
                for ID in list_part_id:
                    sel.SelectPart(ID)
            model.GetSetList().GetSet(name).AddSelected(sel)
        part_list_set('Motion_Region', list_xy_magnets, list_part_id=[id_backiron, id_shaft])

        part_list_set('MagnetSet', list_xy_magnets)
        return True

    def add_magnetic_transient_study(self, app, model, dir_csv_output_folder, study_name):
        logger = logging.getLogger(__name__)
        spmsm_variant = self

        model.CreateStudy("Transient2D", study_name)
        app.SetCurrentStudy(study_name)
        study = model.GetStudy(study_name)

        # SS-ATA
        # study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
        # study.GetStudyProperties().SetValue("SpecifySlip", 0)
        # study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)
        # study.GetStudyProperties().SetValue(u"TimePeriodicType", 2) # This is for TP-EEC but is not effective

        # misc
        GP = self.template.d['GP']
        OP = self.template.d['OP']
        wily = OP['wily']
        study.GetStudyProperties().SetValue("ConversionType", 0)
        study.GetStudyProperties().SetValue("NonlinearMaxIteration", self.max_nonlinear_iteration)
        study.GetStudyProperties().SetValue("ModelThickness", OP['mm_stack_length']) # [mm] Stack Length

        # Material
        self.add_material(study)

        # Conditions - Motion
        study.CreateCondition("RotationMotion", "RotCon") # study.GetCondition(u"RotCon").SetXYZPoint(u"", 0, 0, 1) # megbox warning
        # print('the_speed:', self.the_speed)
        study.GetCondition("RotCon").SetValue("AngularVelocity", int(OP['the_speed']))
        study.GetCondition("RotCon").ClearParts()
        study.GetCondition("RotCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)
        # Implementation of id=0 control:
        #   After rotate the rotor by half the inter-pole notch span, The d-axis initial position is at pole pitch angle divided by 2.
        #   The U-phase current is sin(omega_syn*t) = 0 at t=0 and requires the d-axis to be at the winding phase axis (to obtain id=0 control)
        deg_pole_span = 180/self.template.SD['p']
        #                                                              inter-pole notch (0.5 for half)         rotate to x-axis    winding placing bias (half adjacent slot angle)      reverse north and south pole to make torque positive.
        print(                     '[PMSM JMAG] InitialRotationAngle :',(deg_pole_span-GP['deg_alpha_rm'].value)*0.5, - deg_pole_span*0.5, + wily.deg_winding_U_phase_phase_axis_angle,     + deg_pole_span)
        print(                     '[PMSM JMAG] InitialRotationAngle =',(deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + wily.deg_winding_U_phase_phase_axis_angle     + deg_pole_span, 'deg')
        study.GetCondition("RotCon").SetValue(u"InitialRotationAngle",  (deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + wily.deg_winding_U_phase_phase_axis_angle     + deg_pole_span) 


        study.CreateCondition("Torque", "TorCon") # study.GetCondition(u"TorCon").SetXYZPoint(u"", 0, 0, 0) # megbox warning
        study.GetCondition("TorCon").SetValue("TargetType", 1)
        study.GetCondition("TorCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("TorCon").ClearParts()

        study.CreateCondition("Force", "ForCon")
        study.GetCondition("ForCon").SetValue("TargetType", 1)
        study.GetCondition("ForCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("ForCon").ClearParts()


        # Conditions - FEM Coils & Conductors (i.e. stator/rotor winding)
        self.add_circuit(app, model, study, bool_3PhaseCurrentSource=wily.bool_3PhaseCurrentSource)


        # True: no mesh or field results are needed
        study.GetStudyProperties().SetValue("OnlyTableResults", self.template.fea_config_dict['designer.OnlyTableResults'])

        # Linear Solver
        if False:
            # sometime nonlinear iteration is reported to fail and recommend to increase the accerlation rate of ICCG solver
            study.GetStudyProperties().SetValue("IccgAccel", 1.2) 
            study.GetStudyProperties().SetValue("AutoAccel", 0)
        else:
            # this can be said to be super fast over ICCG solver.
            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

        if self.template.fea_config_dict['designer.MultipleCPUs'] == True:
            # This SMP(shared memory process) is effective only if there are tons of elements. e.g., over 100,000.
            # too many threads will in turn make them compete with each other and slow down the solve. 2 is good enough for eddy current solve. 6~8 is enough for transient solve.
            study.GetStudyProperties().SetValue("UseMultiCPU", True)
            study.GetStudyProperties().SetValue("MultiCPU", 2) 

        # 上一步的铁磁材料的状态作为下一步的初值，挺好，但是如果每一个转子的位置转过很大的话，反而会减慢非线性迭代。
        # 我们的情况是：0.33 sec 分成了32步，每步的时间大概在0.01秒，0.01秒乘以0.5*497 Hz = 2.485 revolution...
        # study.GetStudyProperties().SetValue(u"NonlinearSpeedup", 0) # JMAG17.1以后默认使用。现在后面密集的步长还多一点（32步），前面16步慢一点就慢一点呗！

        # two sections of different time step
        if True:
            number_cycles_in_1stTSS = self.template.fea_config_dict['designer.number_cycles_in_1stTSS']
            number_cycles_in_2ndTSS = self.template.fea_config_dict['designer.number_cycles_in_2ndTSS']
            number_cycles_in_3rdTSS = self.template.fea_config_dict['designer.number_cycles_in_3rdTSS']
            number_cycles_prolonged = self.template.fea_config_dict['designer.number_cycles_prolonged']
            number_of_steps_1stTSS = self.template.fea_config_dict['designer.number_of_steps_1stTSS'] 
            number_of_steps_2ndTSS = self.template.fea_config_dict['designer.number_of_steps_2ndTSS'] 
            number_of_steps_3rdTSS = number_cycles_in_3rdTSS*self.template.fea_config_dict['designer.StepPerCycle_3rdTSS']
            DM = app.GetDataManager()
            DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
            if number_cycles_prolonged == 0:
                if number_cycles_in_3rdTSS == 0:
                    refarray = [[0 for i in range(3)] for j in range(3)]
                    refarray[0][0] = 0
                    refarray[0][1] =    1
                    refarray[0][2] =        50
                    refarray[1][0] = number_cycles_in_1stTSS/OP['DriveW_poles']
                    refarray[1][1] =    number_of_steps_1stTSS
                    refarray[1][2] =        50
                    refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/OP['DriveW_poles']
                    refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                    refarray[2][2] =        50
                else:
                    refarray = [[0 for i in range(3)] for j in range(4)]
                    refarray[0][0] = 0
                    refarray[0][1] =    1
                    refarray[0][2] =        50
                    refarray[1][0] = number_cycles_in_1stTSS/OP['DriveW_poles']
                    refarray[1][1] =    number_of_steps_1stTSS
                    refarray[1][2] =        50
                    refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/OP['DriveW_poles']
                    refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                    refarray[2][2] =        50
                    refarray[3][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS+number_cycles_in_3rdTSS)/OP['DriveW_poles']
                    refarray[3][1] =    number_of_steps_3rdTSS
                    refarray[3][2] =        50
            else:
                refarray = [[0 for i in range(3)] for j in range(4)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = number_cycles_in_1stTSS/OP['DriveW_poles']
                refarray[1][1] =    number_of_steps_1stTSS
                refarray[1][2] =        50
                refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/OP['DriveW_poles']
                refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                refarray[2][2] =        50
                refarray[3][0] = refarray[2][0] + number_cycles_prolonged/OP['DriveW_poles'] 
                refarray[3][1] =    number_cycles_prolonged*self.template.fea_config_dict['designer.TranRef-StepPerCycle'] 
                refarray[3][2] =        50
                # refarray[4][0] = refarray[3][0] + 0.5/OP['DriveW_poles'] # 最后来一个超密的半周期400步
                # refarray[4][1] =    400
                # refarray[4][2] =        50
            number_of_total_steps = 1 + number_of_steps_1stTSS + number_of_steps_2ndTSS + number_of_steps_3rdTSS + number_cycles_prolonged*self.template.fea_config_dict['designer.TranRef-StepPerCycle'] # [Double Check] don't forget to modify here!
            DM.GetDataSet("SectionStepTable").SetTable(refarray)
            study.GetStep().SetValue("Step", number_of_total_steps)
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

        # add equations
        study.GetDesignTable().AddEquation("freq")
        study.GetDesignTable().AddEquation("speed")
        study.GetDesignTable().GetEquation("freq").SetType(0)
        study.GetDesignTable().GetEquation("freq").SetExpression("%g"%((OP['DriveW_Freq'])))
        study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency")
        study.GetDesignTable().GetEquation("speed").SetType(1)
        study.GetDesignTable().GetEquation("speed").SetExpression("freq * %d"%(60/(OP['DriveW_poles']/2)))
        study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed of four pole")

        # speed, freq, slip
        study.GetCondition("RotCon").SetValue("AngularVelocity", 'speed')
        if self.template.spec_input_dict['DPNV_or_SEPA']==False:
            app.ShowCircuitGrid(True)
            study.GetCircuit().GetComponent("CS4").SetValue("Frequency", "freq")
            study.GetCircuit().GetComponent("CS2").SetValue("Frequency", "freq")

        # max_nonlinear_iteration = 50
        # study.GetStudyProperties().SetValue(u"NonlinearMaxIteration", max_nonlinear_iteration)
        # study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
        # study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)

        # # add other excitation frequencies other than 500 Hz as cases
        # for case_no, DriveW_Freq in enumerate([50.0, slip_freq_breakdown_torque]):
        #     slip = slip_freq_breakdown_torque / DriveW_Freq
        #     study.GetDesignTable().AddCase()
        #     study.GetDesignTable().SetValue(case_no+1, 0, DriveW_Freq)
        #     study.GetDesignTable().SetValue(case_no+1, 1, slip)

        # 你把Tran2TSS计算周期减半！
        # 也要在计算铁耗的时候选择1/4或1/2的数据！（建议1/4）
        # 然后，手动添加end step 和 start step，这样靠谱！2019-01-09：注意设置铁耗条件（iron loss condition）的Reference Start Step和End Step。

        # Iron Loss Calculation Condition
        # Stator 
        if True:
            cond = study.CreateCondition("Ironloss", "IronLossConStator")
            cond.SetValue("RevolutionSpeed", "freq*60/%d"%(0.5*(OP['DriveW_poles'])))
            cond.ClearParts()
            sel = cond.GetSelection()
            sel.SelectPartByPosition(self.template.d['GP']['mm_r_si'].value + EPS, 0 ,0)
            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results and to have a FFT plot
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3) # 3:Custom
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTSS*0.5) # 1/4 period <=> number_of_steps_2ndTSS*0.5
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
        # Check CSV reults for iron loss (You cannot check this for Freq study) # CSV and save space
        study.GetStudyProperties().SetValue("CsvOutputPath", dir_csv_output_folder) # it's folder rather than file!
        study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss")
        study.GetStudyProperties().SetValue("DeleteResultFiles", self.template.fea_config_dict['delete_results_after_calculation'])
        # Terminal Voltage/Circuit Voltage: Check for outputing CSV results 
        study.GetCircuit().CreateTerminalLabel("TerminalGroupACU", 8, -13)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupACV", 8, -11)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupACW", 8, -9)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupBDU", 23, -13)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupBDV", 23, -11)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupBDW", 23, -9)
        # Export Stator Core's field results only for iron loss calculation (the csv file of iron loss will be clean with this setting)
            # study.GetMaterial(u"Rotor Core").SetValue(u"OutputResult", 0) # at least one part on the rotor should be output or else a warning "the jplot file does not contains displacement results when you try to calc. iron loss on the moving part." will pop up, even though I don't add iron loss condition on the rotor.
        # study.GetMeshControl().SetValue(u"AirRegionOutputResult", 0)
        # study.GetMaterial("Shaft").SetValue("OutputResult", 0)
        # study.GetMaterial("Cage").SetValue("OutputResult", 0)
        # study.GetMaterial("Coil").SetValue("OutputResult", 0)
        # Rotor
        if True:
            cond = study.CreateCondition("Ironloss", "IronLossConRotor")
            cond.SetValue("BasicFrequencyType", 2)
            cond.SetValue("BasicFrequency", "freq")
                # cond.SetValue(u"BasicFrequency", u"slip*freq") # this require the signal length to be at least 1/4 of slip period, that's too long!
            cond.ClearParts()
            sel = cond.GetSelection()
            # sel.SelectPartByPosition(self.mm_r_ri + EPS, 0 ,0) # Why this is not working??? Because it is integer.... you must use 0.0 instead of 0!!!
            sel.SelectPart(self.id_backiron)

            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3)
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTSS*0.5) # 1/4 period <=> number_of_steps_2ndTSS*0.5
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
        self.study_name = study_name
        return study

    def add_structural_static_study(self):
        pass

    def add_mesh(self, study, model):
        # this is for multi slide planes, which we will not be using
        refarray = [[0 for i in range(2)] for j in range(1)]
        refarray[0][0] = 3
        refarray[0][1] = 1
        study.GetMeshControl().GetTable("SlideTable2D").SetTable(refarray) 

        study.GetMeshControl().SetValue("MeshType", 1) # make sure this has been exe'd: study.GetCondition(u"RotCon").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)
        study.GetMeshControl().SetValue("RadialDivision", 4) # for air region near which motion occurs
        study.GetMeshControl().SetValue("CircumferentialDivision", 720) #1440) # for air region near which motion occurs 这个数足够大，sliding mesh才准确。
        study.GetMeshControl().SetValue("AirRegionScale", 1.05) # [Model Length]: Specify a value within the following area. (1.05 <= value < 1000)
        study.GetMeshControl().SetValue("MeshSize", 4) # mm
        study.GetMeshControl().SetValue("AutoAirMeshSize", 0)
        study.GetMeshControl().SetValue("AirMeshSize", 4) # mm
        study.GetMeshControl().SetValue("Adaptive", 0)

        # This is not neccessary for whole model FEA. In fact, for BPMSM simulation, it causes mesh error "The copy target region is not found".
        # study.GetMeshControl().CreateCondition("RotationPeriodicMeshAutomatic", "autoRotMesh") # with this you can choose to set CircumferentialDivision automatically

        study.GetMeshControl().CreateCondition("Part", "MagnetMeshCtrl")
        study.GetMeshControl().GetCondition("MagnetMeshCtrl").SetValue("Size", self.meshSize_Magnet)
        study.GetMeshControl().GetCondition("MagnetMeshCtrl").ClearParts()
        study.GetMeshControl().GetCondition("MagnetMeshCtrl").AddSet(model.GetSetList().GetSet("MagnetSet"), 0)

        study.GetMeshControl().CreateCondition("Part", "ShaftMeshCtrl")
        study.GetMeshControl().GetCondition("ShaftMeshCtrl").SetValue("Size", 10) # 10 mm
        study.GetMeshControl().GetCondition("ShaftMeshCtrl").ClearParts()
        study.GetMeshControl().GetCondition("ShaftMeshCtrl").AddSet(model.GetSetList().GetSet("ShaftSet"), 0)

        def mesh_all_cases(study):
            numCase = study.GetDesignTable().NumCases()
            for case in range(0, numCase):
                study.SetCurrentCase(case)
                if study.HasMesh() == False:
                    study.CreateMesh()
                # if case == 0:
                #     app.View().ShowAllAirRegions()
                #     app.View().ShowMeshGeometry()
                #     app.View().ShowMesh()

        mesh_all_cases(study)

    # TranFEAwi2TSS
    def add_material(self, study):
        if 'M19' in self.template.spec_input_dict['Steel']:
            study.SetMaterialByName("StatorCore", "M-19 Steel Gauge-29")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 95)
                # study.GetMaterial(u"Stator Core").SetValue(u"UserConductivityValue", 1900000)

            study.SetMaterialByName("NotchedRotor", "M-19 Steel Gauge-29")
            study.GetMaterial("NotchedRotor").SetValue("Laminated", 1)
            study.GetMaterial("NotchedRotor").SetValue("LaminationFactor", 95)

        elif 'M15' in self.template.spec_input_dict['Steel']:
            study.SetMaterialByName("StatorCore", "M-15 Steel")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 98)

            study.SetMaterialByName("NotchedRotor", "M-15 Steel")
            study.GetMaterial("NotchedRotor").SetValue("Laminated", 1)
            study.GetMaterial("NotchedRotor").SetValue("LaminationFactor", 98)

        elif self.template.spec_input_dict['Steel'] == 'Arnon5':
            study.SetMaterialByName("StatorCore", "Arnon5-final")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 96)

            study.SetMaterialByName("NotchedRotor", "Arnon5-final")
            study.GetMaterial("NotchedRotor").SetValue("Laminated", 1)
            study.GetMaterial("NotchedRotor").SetValue("LaminationFactor", 96)

        else:
            msg = 'Warning: default material is used: DCMagnetic Type/50A1000.'
            print(msg)
            logging.getLogger(__name__).warn(msg)
            study.SetMaterialByName("StatorCore", "DCMagnetic Type/50A1000")
            study.GetMaterial("StatorCore").SetValue("UserConductivityType", 1)
            study.SetMaterialByName("NotchedRotor", "DCMagnetic Type/50A1000")
            study.GetMaterial("NotchedRotor").SetValue("UserConductivityType", 1)

        study.SetMaterialByName("Coils", "Copper")
        study.GetMaterial("Coils").SetValue("UserConductivityType", 1)

        # study.SetMaterialByName("Cage", "Aluminium")
        # study.GetMaterial("Cage").SetValue("EddyCurrentCalculation", 1)
        # study.GetMaterial("Cage").SetValue("UserConductivityType", 1)
        # study.GetMaterial("Cage").SetValue("UserConductivityValue", self.Bar_Conductivity)

        # N40H Reversible
        study.SetMaterialByName(u"Magnet", u"Arnold/Reversible/N40H")
        study.GetMaterial(u"Magnet").SetValue(u"EddyCurrentCalculation", 1)
        available_temperature_list = [-40, 20, 60, 80, 100, 120, 150, 180, 200, 220] # according to JMAG
        magnet_temperature = min(available_temperature_list, key=lambda x:abs(x-self.template.spec_input_dict['Temperature']))        
        print('magnet_temperature is', magnet_temperature, 'deg C.')
        study.GetMaterial(u"Magnet").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)
        study.GetMaterial(u"Magnet").SetValue(u"Poles", self.template.d['OP']['DriveW_poles'])

        study.GetMaterial(u"Magnet").SetDirectionXYZ(1, 0, 0)
        study.GetMaterial(u"Magnet").SetAxisXYZ(0, 0, 1)
        study.GetMaterial(u"Magnet").SetOriginXYZ(0, 0, 0)
        study.GetMaterial(u"Magnet").SetPattern(u"RadialCircular")
        study.GetMaterial(u"Magnet").SetValue(u"StartAngle", 0.5* 360/(2*self.template.SD['p']) ) # 半个极距


        # add_carbon_fiber_material(app)

    def add_circuit(self, app, model, study, bool_3PhaseCurrentSource=True):
        # Circuit - Current Source
        app.ShowCircuitGrid(True)
        study.CreateCircuit()

        # 4 pole motor Qs=24 dpnv implemented by two layer winding (6 coils). In this case, drive winding has the same slot turns as bearing winding
        def circuit(Grouping,turns,Rs,ampD,ampB,freq,phase=0, CommutatingSequenceD=0, CommutatingSequenceB=0, x=10,y=10, bool_3PhaseCurrentSource=True):
            study.GetCircuit().CreateSubCircuit("Star Connection", "Star Connection %s"%(Grouping), x, y) # è¿™äº›æ•°å­—æŒ‡çš„æ˜¯gridçš„ä¸ªæ•°ï¼Œç¬¬å‡ è¡Œç¬¬å‡ åˆ—çš„æ ¼ç‚¹å¤„
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil1").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil1").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil2").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil2").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil3").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil3").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil1").SetName("CircuitCoil%sU"%(Grouping))
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil2").SetName("CircuitCoil%sV"%(Grouping))
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil3").SetName("CircuitCoil%sW"%(Grouping))
            # Star Connection_2 is GroupAC
            # Star Connection_4 is GroupBD

            if bool_3PhaseCurrentSource == True: # must use this for frequency analysis

                study.GetCircuit().CreateComponent("3PhaseCurrentSource", "CS%s"%(Grouping))
                study.GetCircuit().CreateInstance("CS%s"%(Grouping), x-4, y+1)
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("Amplitude", ampD+ampB)
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("Frequency", "freq") # this is not needed for freq analysis # "freq" is a variable
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("PhaseU", phase)
                # Commutating sequence is essencial for the direction of the field to be consistent with speed: UVW rather than UWV
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("CommutatingSequence", CommutatingSequenceD) 
            else: 
                I1 = "CS%s-1"%(Grouping)
                I2 = "CS%s-2"%(Grouping)
                I3 = "CS%s-3"%(Grouping)
                study.GetCircuit().CreateComponent("CurrentSource", I1)
                study.GetCircuit().CreateInstance(                   I1, x-4, y+3)
                study.GetCircuit().CreateComponent("CurrentSource", I2)
                study.GetCircuit().CreateInstance(                   I2, x-4, y+1)
                study.GetCircuit().CreateComponent("CurrentSource", I3)
                study.GetCircuit().CreateInstance(                   I3, x-4, y-1)

                phase_shift_drive = -120 if CommutatingSequenceD == 1 else 120
                phase_shift_beari = -120 if CommutatingSequenceB == 1 else 120

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 0*phase_shift_drive) # "freq" variable cannot be used here. So pay extra attension here when you create new case of a different freq.
                f2 = app.FunctionFactory().Sin(ampB, freq, 0*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I1).SetFunction(func)

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 1*phase_shift_drive)
                f2 = app.FunctionFactory().Sin(ampB, freq, 1*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I2).SetFunction(func)

                func = app.FunctionFactory().Composite()
                f1 = app.FunctionFactory().Sin(ampD, freq, 2*phase_shift_drive)
                f2 = app.FunctionFactory().Sin(ampB, freq, 2*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I3).SetFunction(func)

            study.GetCircuit().CreateComponent("Ground", "Ground")
            study.GetCircuit().CreateInstance("Ground", x+2, y+1)
        # 这里电流幅值中的0.5因子源自DPNV导致的等于2的平行支路数。没有考虑到这一点，是否会对initial design的有效性产生影响？
        # 仔细看DPNV的接线，对于转矩逆变器，绕组的并联支路数为2，而对于悬浮逆变器，绕组的并联支路数为1。

        OP = self.template.d['OP']
        wily = OP['wily']

        npb = wily.number_parallel_branch
        nwl = wily.number_winding_layer # number of windign layers 
        # if self.template.fea_config_dict['DPNV_separate_winding_implementation'] == True or self.template.spec_input_dict['DPNV_or_SEPA'] == False:
        if self.template.spec_input_dict['DPNV_or_SEPA'] == False:
            # either a separate winding or a DPNV winding implemented as a separate winding
            ampD =  0.5 * (OP['DriveW_CurrentAmp']/npb + self.BeariW_CurrentAmp) # 为了代码能被四极电机和二极电机通用，代入看看就知道啦。
            ampB = -0.5 * (OP['DriveW_CurrentAmp']/npb - self.BeariW_CurrentAmp) # 关于符号，注意下面的DriveW对应的circuit调用时的ampB前还有个负号！
            if bool_3PhaseCurrentSource != True:
                raise Exception('Logic Error Detected.')
        else:
            # case: DPNV as an actual two layer winding
            ampD = OP['DriveW_CurrentAmp']/npb
            ampB = OP['BeariW_CurrentAmp']
            if bool_3PhaseCurrentSource != False:
                raise Exception('Logic Error Detected.')

        circuit('GroupAC',  OP['DriveW_zQ']/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=OP['DriveW_Rs'],ampD= ampD,
                              ampB=-ampB, freq=self.template.d['OP']['DriveW_Freq'], phase=0,
                              CommutatingSequenceD=wily.CommutatingSequenceD,
                              CommutatingSequenceB=wily.CommutatingSequenceB)
        circuit('GroupBD',  OP['BeariW_zQ']/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=OP['BeariW_Rs'],ampD= ampD,
                              ampB=+ampB, freq=self.template.d['OP']['BeariW_Freq'], phase=0,
                              CommutatingSequenceD=wily.CommutatingSequenceD,
                              CommutatingSequenceB=wily.CommutatingSequenceB,x=25) # CS4 corresponds to uauc (conflict with following codes but it does not matter.)

        # Link FEM Coils to Coil Set     
        # if self.template.fea_config_dict['DPNV_separate_winding_implementation'] == True or self.template.spec_input_dict['DPNV_or_SEPA'] == False:
        if self.template.spec_input_dict['DPNV_or_SEPA'] == False:
            def link_FEMCoils_2_CoilSet(Grouping,l1,l2):
                # link between FEM Coil Condition and Circuit FEM Coil
                for UVW in ['U','V','W']:
                    which_phase = "%s%s-Phase"%(Grouping,UVW)
                    study.CreateCondition("FEMCoil", which_phase)
                    condition = study.GetCondition(which_phase)
                    condition.SetLink("CircuitCoil%s%s"%(Grouping,UVW))
                    condition.GetSubCondition("untitled").SetName("Coil Set 1")
                    condition.GetSubCondition("Coil Set 1").SetName("delete")
                count = 0
                dict_dir = {'+':1, '-':0, 'o':None}
                # select the part to assign the FEM Coil condition
                for UVW, UpDown in zip(l1,l2):
                    count += 1 
                    if dict_dir[UpDown] is None:
                        # print 'Skip', UVW, UpDown
                        continue
                    which_phase = "%s%s-Phase"%(Grouping,UVW)
                    condition = study.GetCondition(which_phase)
                    condition.CreateSubCondition("FEMCoilData", "Coil Set %d"%(count))
                    subcondition = condition.GetSubCondition("Coil Set %d"%(count))
                    subcondition.ClearParts()
                    subcondition.AddSet(model.GetSetList().GetSet("Coil%s%s%s %d"%(Grouping,UVW,UpDown,count)), 0) # 未修改，这里的和Set对不上的，Grouping要换成LX或LY
                    subcondition.SetValue("Direction2D", dict_dir[UpDown])
                # clean up
                for UVW in ['U','V','W']:
                    which_phase = "%s%s-Phase"%(Grouping,UVW)
                    condition = study.GetCondition(which_phase)
                    condition.RemoveSubCondition("delete")
            link_FEMCoils_2_CoilSet('GroupAC', # 这边还是古老的，还未修改啊。。。
                                    wily.dict_coil_connection['layer X phases'], # 40 for 4 poles, +1 for UVW, 
                                    wily.dict_coil_connection['layer X signs'])                 # +2 for up or down, 
            link_FEMCoils_2_CoilSet('GroupBD', 
                                    wily.dict_coil_connection['layer Y phases'], # 20 for 2 poles, +1 for UVW, .
                                    wily.dict_coil_connection['layer Y signs'])                 # +2 for up or down,  这里的2和4等价于leftlayer和rightlayer。
        else:
            # 两个改变，一个是激励大小的改变（本来是200A 和 5A，现在是205A和195A），
            # 另一个绕组分组的改变，现在的A相是上层加下层为一相，以前是用俩单层绕组等效的。

            # Link FEM Coils to Coil Set as double layer short pitched winding
            # Create FEM Coil Condition
            # here we map circuit component `Coil2A' to FEM Coil Condition 'phaseAuauc
            # here we map circuit component `Coil4A' to FEM Coil Condition 'phaseAubud
            for suffix, poles in zip(['GroupAC', 'GroupBD'], [OP['DriveW_poles'], OP['BeariW_poles']]): # 仍然需要考虑poles，是因为为Coil设置Set那里的代码还没有更新。这里的2(self.DriveW_poles)和4(self.BeariW_poles)等价于leftlayer和rightlayer。
                for UVW in ['U','V','W']:
                    study.CreateCondition("FEMCoil", 'phase'+UVW+suffix)
                    # link between FEM Coil Condition and Circuit FEM Coil
                    condition = study.GetCondition('phase'+UVW+suffix)
                    condition.SetLink("CircuitCoil%s%s"%(suffix,UVW))
                    condition.GetSubCondition("untitled").SetName("delete")
            countXL = 0 # countXL indicates which slot the current rightlayer is in.
            index = 0
            dict_dir = {'+':1, '-':0}
            coil_pitch = wily.coil_pitch_y #wily.dict_coil_connection[0]
            # select the part (via `Set') to assign the FEM Coil condition
            for UVW, UpDown in zip(wily.layer_X_phases, wily.layer_X_signs):

                countXL += 1 
                if wily.grouping_AC[index] == 1:
                    suffix = 'GroupAC'
                else:
                    suffix = 'GroupBD'
                condition = study.GetCondition('phase'+UVW+suffix)

                # right layer
                # print (countXL, "Coil Set %d"%(countXL), end=' ')
                condition.CreateSubCondition("FEMCoilData", "Coil Set Layer X %d"%(countXL))
                subcondition = condition.GetSubCondition("Coil Set Layer X %d"%(countXL))
                subcondition.ClearParts()
                subcondition.AddSet(model.GetSetList().GetSet("CoilLX%s%s %d"%(UVW,UpDown,countXL)), 0) # poles=4 means right layer, rather than actual poles
                subcondition.SetValue("Direction2D", dict_dir[UpDown])

                # left layer
                Q = self.template.SD['Qs']
                if coil_pitch > 0:
                    if countXL+coil_pitch <= Q:
                        count_leftlayer = countXL+coil_pitch
                        index_leftlayer = index+coil_pitch
                    else:
                        count_leftlayer = int(countXL+coil_pitch - Q)
                        index_leftlayer = int(index+coil_pitch - Q)
                else:
                    if countXL+coil_pitch > 0:
                        count_leftlayer = countXL+coil_pitch
                        index_leftlayer = index+coil_pitch
                    else:
                        count_leftlayer = int(countXL+coil_pitch + Q)
                        index_leftlayer = int(index+coil_pitch + Q)

                # Check if it is a distributed windg???
                if wily.distributed_or_concentrated == False:
                    print('Concentrated winding!')
                    UVW    = wily.layer_Y_phases[index_leftlayer]
                    UpDown = wily.layer_Y_signs[index_leftlayer]
                else:
                    # print('Distributed winding.')
                    if wily.layer_Y_phases[index_leftlayer] != UVW:
                        print('[Warn] Potential bug in your winding layout detected.')
                        raise Exception('Bug in winding layout detected.')
                    # 右层导体的电流方向是正，那么与其串联的一个coil_pitch之处的左层导体就是负！不需要再检查l_leftlayer2了~
                    if UpDown == '+': 
                        UpDown = '-'
                    else:
                        UpDown = '+'
                # print (count_leftlayer, "Coil Set %d"%(count_leftlayer))
                condition.CreateSubCondition("FEMCoilData", "Coil Set Layer Y %d"%(count_leftlayer))
                subcondition = condition.GetSubCondition("Coil Set Layer Y %d"%(count_leftlayer))
                subcondition.ClearParts()
                subcondition.AddSet(model.GetSetList().GetSet("CoilLY%s%s %d"%(UVW,UpDown,count_leftlayer)), 0) # poles=2 means left layer, rather than actual poles
                subcondition.SetValue("Direction2D", dict_dir[UpDown])
                # print 'coil_pitch=', coil_pitch
                # print layer_X_phases[index], UVW
                # print l_leftlayer1[index_leftlayer]
                # print layer_X_phases
                # print l_leftlayer1
                index += 1
            # clean up
            for suffix in ['GroupAC', 'GroupBD']:
                for UVW in ['U','V','W']:
                    condition = study.GetCondition('phase'+UVW+suffix)
                    condition.RemoveSubCondition("delete")
            # raise Exception('Test DPNV PE.')

    ''' JMAG Description
    '''
    def show(self, toString=False):
        def get_tuple_list(object):
            attrs = list(vars(object).items())
            key_list = [el[0] for el in attrs]
            val_list = [el[1] for el in attrs]
            the_dict = dict(list(zip(key_list, val_list)))
            sorted_key = sorted(key_list, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item)) # this is also useful for string beginning with digiterations '15 Steel'.
            tuple_list = [(key, the_dict[key]) for key in sorted_key]
            return tuple_list
        variant_tuple_list = get_tuple_list(self)
        template_tuple_list = get_tuple_list(self.template)
        if toString==False:
            print('- Bearingless PMSM Individual #%s\n\t' % (self.name), end=' ')
            print(', \n\t'.join("%s = %s" % item for item in tuple_list))
            return ''
        else:
            return '\n- Bearingless PMSM Individual #%s\n\t' % (self.name) \
                    + ', \n\t'.join("%s = %s" % item for item in variant_tuple_list) \
                    + '\n' + '--'*10 + '\n- Template:\n\t' \
                    + ', \n\t'.join("%s = %s" % item for item in template_tuple_list)

    def get_individual_name(self):
        if self.template.fea_config_dict['flag_optimization'] == True:
            return "ID%s" % (self.ID)
        else:
            return "%s_ID%s" % (self.name, self.ID)

    # Produce JMAG Project
    def build_jmag_project(self, project_meta_data=None, bool_re_evaluate=False):
        if project_meta_data is None:
            # this object is reloaded, so retrieve saved meta data
            project_meta_data = self.project_meta_data
        else:
            # save for reload this object from json
            self.project_meta_data = project_meta_data
        expected_project_file = project_meta_data["expected_project_file"]
        study_name            = project_meta_data["study_name"]
        dir_csv_output_folder = project_meta_data["dir_csv_output_folder"]
        output_dir            = project_meta_data["output_dir"]

        # use alias
        spmsm_variant = self

        # Leave the solving task to JMAG
        if bool_re_evaluate==False:

            # def draw_jmag_bpmsm():
            import JMAG
            toolJd = JMAG.JMAG(self.template.fea_config_dict, self.template.spec_input_dict)

            toolJd.open(expected_project_file)
            DRAW_SUCCESS = spmsm_variant.draw_spmsm(toolJd)
            if DRAW_SUCCESS != 1:
                raise Exception('Drawer failed.')

            app = toolJd.app

            # JMAG
            if app.NumModels()>=1:
                model = app.GetModel(spmsm_variant.name)
            else:
                logger.error('there is no model yet for %s'%(spmsm_variant.name))
                raise Exception('why is there no model yet? %s'%(spmsm_variant.name))

            spmsm_variant.pre_process(app, model)
            study = spmsm_variant.add_magnetic_transient_study(app, model, dir_csv_output_folder, study_name) # Change here and there 
            self.mesh_study(spmsm_variant, app, model, study, output_dir=output_dir)
            # raise KeyboardInterrupt
            self.run_study(spmsm_variant, app, study, self.template.fea_config_dict, clock_time())


            # export Voltage if field data exists.
            if self.template.fea_config_dict['delete_results_after_calculation'] == False:
                # Export Circuit Voltage
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(dir_csv_output_folder + study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
        else:

            THE_mm2_magnet_area = self.spmsm_variant.rotorMagnet.draw(None, bool_re_evaluate=True)
            THE_mm2_slot_area = self.spmsm_variant.coils.draw(None, bool_re_evaluate=True)

            CurrentAmp_in_the_slot = THE_mm2_slot_area * spmsm_variant.fill_factor * spmsm_variant.Js*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / spmsm_variant.template.d['OP']['DriveW_zQ']
            CurrentAmp_per_phase = CurrentAmp_per_conductor * spmsm_variant.template.d['OP']['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            spmsm_variant.CurrentAmp_per_phase = CurrentAmp_per_phase
            spmsm_variant.DriveW_CurrentAmp = spmsm_variant.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            spmsm_variant.BeariW_CurrentAmp = spmsm_variant.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            slot_area_utilizing_ratio = (spmsm_variant.DriveW_CurrentAmp + spmsm_variant.BeariW_CurrentAmp) / spmsm_variant.CurrentAmp_per_phase
            print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)

            print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
            print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
            print('---spmsm_variant.DriveW_CurrentAmp =', spmsm_variant.DriveW_CurrentAmp)
            print('---spmsm_variant.BeariW_CurrentAmp =', spmsm_variant.BeariW_CurrentAmp)
            print('---TORQUE_CURRENT_RATIO:', spmsm_variant.fea_config_dict['TORQUE_CURRENT_RATIO'])
            print('---SUSPENSION_CURRENT_RATIO:', spmsm_variant.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

    @classmethod
    def run_study(cls, im_variant, app, study, fea_config_dict, toc):
        logger = logging.getLogger(__name__)
        if fea_config_dict['designer.JMAG_Scheduler'] == False:
            print('Run jam.exe...')
            # if run_list[1] == True:
            try:
                study.RunAllCases()
            except Exception as error:
                raise error
            msg = 'Time spent on %s is %g s.'%(study.GetName() , clock_time() - toc)
            logger.debug(msg)
            print(msg)
        else:
            print('Submit to JMAG_Scheduler...')
            job = study.CreateJob()
            job.SetValue("Title", study.GetName())
            job.SetValue("Queued", True)
            job.Submit(False) # Fallse:CurrentCase, True:AllCases
            logger.debug('Submit %s to queue (Tran2TSS).'%(im_variant.individual_name))
            # wait and check
            # study.CheckForCaseResults()
        app.Save()
        # if the jcf file already exists, it pops a msg window
        # study.WriteAllSolidJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Solid', True) # True : Outputs cases that do not have results 
        # study.WriteAllMeshJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Mesh', True)

        # # run
        # if self.template.fea_config_dict['JMAG_Scheduler'] == False:
        #     study.RunAllCases()
        #     app.Save()
        # else:
        #     job = study.CreateJob()
        #     job.SetValue(u"Title", study.GetName())
        #     job.SetValue(u"Queued", True)
        #     job.Submit(True)
        #     logger.debug('Submit %s to queue (Freq).'%(im_variant.individual_name))
        #     # wait and check
        #     # study.CheckForCaseResults()

    @classmethod
    def mesh_study(cls, im_variant, app, model, study, output_dir):

        # this `if' judgment is effective only if JMAG-DeleteResultFiles is False 
        # if not study.AnyCaseHasResult(): 
        # mesh
        im_variant.add_mesh(study, model)

        # Export Image
        app.View().ShowAllAirRegions()
        # app.View().ShowMeshGeometry() # 2nd btn
        app.View().ShowMesh() # 3rn btn
        app.View().Zoom(3)
        app.View().Pan(-im_variant.template.d['GP']['mm_r_or'].value, 0)
        app.ExportImageWithSize(output_dir + model.GetName() + '.png', 2000, 2000)
        app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.

# circumferential segmented rotor 
if __name__ == '__main__':
    import JMAG
    import Location2D
    import CrossSectInnerNotchedRotor
    import CrossSectStator
    # from pylab import np

    if True:
        from utility import my_execfile
        my_execfile('./default_setting.py', g=globals(), l=locals())
        fea_config_dict
        toolJd = JMAG.JMAG(fea_config_dict)

        project_name          = 'proj%d'%(0)
        expected_project_file_path = './' + "%s.jproj"%(project_name)
        toolJd.open(expected_project_file_path)

    spmsm_template = bearingless_spmsm_template()
    spmsm_template.fea_config_dict = fea_config_dict
    Q = 6

    spmsm_template.deg_alpha_st = 360/Q*0.8   # deg_alpha_st # span angle of tooth: class type DimAngular
    spmsm_template.deg_alpha_so = 0                          # deg_alpha_so # angle of tooth edge: class type DimAngular
    spmsm_template.mm_r_si      = 50     # mm_r_si           # inner radius of stator teeth: class type DimLinear
    spmsm_template.mm_d_so      = 5      # mm_d_so           # tooth edge length: class type DimLinear
    spmsm_template.mm_d_sp      = 1.5*spmsm_template.mm_d_so # mm_d_sp      # tooth tip length: class type DimLinear
    spmsm_template.mm_d_st      = 15     # mm_d_st      # tooth base length: class type DimLinear
    spmsm_template.mm_d_sy      = 15     # mm_d_sy      # back iron thickness: class type DimLinear
    spmsm_template.mm_w_st      = 13     # mm_w_st      # tooth base width: class type DimLinear
    spmsm_template.mm_r_st      = 0         # mm_r_st      # fillet on outter tooth: class type DimLinear
    spmsm_template.mm_r_sf      = 0         # mm_r_sf      # fillet between tooth tip and base: class type DimLinear
    spmsm_template.mm_r_sb      = 0         # mm_r_sb      # fillet at tooth base: class type DimLinear
    spmsm_template.Q            = 6      # number of stator slots (integer)
    spmsm_template.sleeve_length        = 2 # mm
    spmsm_template.fixed_air_gap_length = 0.75 # mm
    spmsm_template.mm_d_pm      = 6      # mm_d_pm          # manget depth
    spmsm_template.deg_alpha_rm = 60     # deg_alpha_rm     # angular span of the pole: class type DimAngular
    spmsm_template.deg_alpha_rs = 10     # deg_alpha_rs     # segment span: class type DimAngular
    spmsm_template.mm_d_ri      = 8      # mm_d_ri          # rotor iron thickness: class type DimLinear
    spmsm_template.mm_r_ri      = 40     # mm_r_ri          # inner radius of rotor: class type DimLinear
    spmsm_template.mm_d_rp      = 5      # mm_d_rp          # interpolar iron thickness: class type DimLinear
    spmsm_template.mm_d_rs      = 3      # mm_d_rs          # inter segment iron thickness: class type DimLinear
    spmsm_template.p = 2     # p     # number of pole pairs
    spmsm_template.s = 3     # s     # number of segments  

    build_design_parameters_list(spmsm_template)

    spmsm_template.DriveWinding_Freq       = 1000
    spmsm_template.DriveWinding_Rs         = 0.1 # TODO
    spmsm_template.DriveWinding_zQ         = 1
    spmsm_template.DriveWinding_CurrentAmp = None # this means it depends on the slot area
    spmsm_template.DriveWinding_poles = 2*spmsm_template.p

    spmsm_template.Js = 4e6 # Arms/m^2
    spmsm_template.fill_factor = 0.45

    spmsm_template.stack_length = 100 # mm

    # logger = logging.getLogger(__name__) 
    # logger.info('spmsm_variant ID %s is initialized.', self.name)

    spmsm = bearingless_spmsm_design(   spmsm_template=spmsm_template,
                                        free_variables=None,
                                        counter=None, 
                                        counter_loop=None
                                        )
    # Rotor Core
    list_segments = spmsm.rotorCore.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = 0 #spmsm.rotorCore.p*2
    region1 = toolJd.prepareSection(list_segments)

    # Rotor Magnet    
    list_regions = spmsm.rotorMagnet.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = 0 #spmsm.rotorMagnet.notched_rotor.p*2
    region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')


# notched rotor 
if __name__ == '__main__':
    import JMAG
    import Location2D
    import CrossSectInnerNotchedRotor
    import CrossSectStator
    # from pylab import np

    if True:
        from utility import my_execfile
        my_execfile('./default_setting.py', g=globals(), l=locals())
        fea_config_dict
        toolJd = JMAG.JMAG(fea_config_dict, spec_input_dict=None)

        project_name          = 'proj%d'%(0)
        expected_project_file_path = './' + "%s.jproj"%(project_name)
        toolJd.open(expected_project_file_path)

    spmsm_template = bearingless_spmsm_template()
    spmsm_template.fea_config_dict = fea_config_dict
    Q = 6

    spmsm_template.deg_alpha_st = 360/Q*0.8   # deg_alpha_st # span angle of tooth: class type DimAngular
    spmsm_template.deg_alpha_so = 0                          # deg_alpha_so # angle of tooth edge: class type DimAngular
    spmsm_template.mm_r_si      = 50     # mm_r_si           # inner radius of stator teeth: class type DimLinear
    spmsm_template.mm_d_so      = 5      # mm_d_so           # tooth edge length: class type DimLinear
    spmsm_template.mm_d_sp      = 1.5*spmsm_template.mm_d_so # mm_d_sp      # tooth tip length: class type DimLinear
    spmsm_template.mm_d_st      = 15     # mm_d_st      # tooth base length: class type DimLinear
    spmsm_template.mm_d_sy      = 15     # mm_d_sy      # back iron thickness: class type DimLinear
    spmsm_template.mm_w_st      = 13     # mm_w_st      # tooth base width: class type DimLinear
    spmsm_template.mm_r_st      = 0         # mm_r_st      # fillet on outter tooth: class type DimLinear
    spmsm_template.mm_r_sf      = 0         # mm_r_sf      # fillet between tooth tip and base: class type DimLinear
    spmsm_template.mm_r_sb      = 0         # mm_r_sb      # fillet at tooth base: class type DimLinear
    spmsm_template.Q            = 6      # number of stator slots (integer)
    spmsm_template.sleeve_length        = 2 # mm
    spmsm_template.fixed_air_gap_length = 0.75 # mm
    spmsm_template.mm_d_pm      = 6      # mm_d_pm          # manget depth
    spmsm_template.deg_alpha_rm = 60     # deg_alpha_rm     # angular span of the pole: class type DimAngular
    spmsm_template.deg_alpha_rs = 60     # deg_alpha_rs     # segment span: class type DimAngular
    spmsm_template.mm_d_ri      = 8      # mm_d_ri          # inner radius of rotor: class type DimLinear
    spmsm_template.mm_r_ri      = 40     # mm_r_ri          # rotor iron thickness: class type DimLinear
    spmsm_template.mm_d_rp      = 5      # mm_d_rp          # interpolar iron thickness: class type DimLinear
    spmsm_template.mm_d_rs      = 0*3      # mm_d_rs          # inter segment iron thickness: class type DimLinear
    spmsm_template.p = 2     # p     # number of pole pairs
    spmsm_template.s = 1     # s     # number of segments  

    build_design_parameters_list(spmsm_template)

    spmsm_template.DriveWinding_Freq       = 1000
    spmsm_template.DriveWinding_Rs         = 0.1 # TODO
    spmsm_template.DriveWinding_zQ         = 1
    spmsm_template.DriveWinding_CurrentAmp = None # this means it depends on the slot area
    spmsm_template.DriveWinding_poles = 2*spmsm_template.p

    spmsm_template.Js = 4e6 # Arms/m^2
    spmsm_template.fill_factor = 0.45

    spmsm_template.stack_length = 100 # mm

    # logger = logging.getLogger(__name__) 
    # logger.info('spmsm_variant ID %s is initialized.', self.name)

    spmsm = bearingless_spmsm_design(   spmsm_template=spmsm_template,
                                        free_variables=None,
                                        counter=None, 
                                        counter_loop=None
                                        )
    # Rotor Core
    list_segments = spmsm.rotorCore.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = spmsm.rotorCore.p*2
    region1 = toolJd.prepareSection(list_segments)

    # Rotor Magnet    
    list_regions = spmsm.rotorMagnet.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = spmsm.rotorMagnet.notched_rotor.p*2
    region2 = toolJd.prepareSection(list_regions)

    # Rotor Magnet
    sleeve = CrossSectInnerNotchedRotor.CrossSectSleeve(
                    name = 'Sleeve',
                    notched_magnet = spmsm.rotorMagnet,
                    d_sleeve = spmsm_template.sleeve_length
                    )

    list_regions = sleeve.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = spmsm.rotorMagnet.notched_rotor.p*2
    regionS = toolJd.prepareSection(list_regions)


    # # Stator Core
    # list_regions = spmsm.stator_core.draw(toolJd)
    # toolJd.bMirror = True
    # toolJd.iRotateCopy = spmsm.stator_core.Q
    # region1 = toolJd.prepareSection(list_regions)

    # # Stator Winding
    # list_regions = spmsm.coils.draw(toolJd)
    # toolJd.bMirror = False
    # toolJd.iRotateCopy = spmsm.coils.stator_core.Q
    # region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')


def add_carbon_fiber_material(app):
    app.GetMaterialLibrary().CreateCustomMaterial(u"CarbonFiber", u"Custom Materials")
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"Density", 1.6)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"CoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"DemagnetizationCoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"MagnetizationSaturated", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulus", 110000)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulus", 5000)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"PoissonRatio", 0.1)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"Thermal Expansion", 8.4e-06)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulusX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulusY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulusZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulusXY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulusYZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulusZX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G11", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G12", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G13", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G14", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G15", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G16", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G22", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G23", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G24", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G25", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G26", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G33", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G34", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G35", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G36", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G44", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G45", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G46", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G55", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G56", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G66", 0)





    # -*- coding: utf-8 -*-
    app = designer.GetApplication()
    app.GetMaterialLibrary().CopyMaterial(u"Arnold/Reversible/N40H")
    app.GetMaterialLibrary().GetUserMaterial(u"N40H(reversible) copy").SetValue(u"Name", u"MyN40H(reversible)")
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"Density", 7.5)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"CoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"DemagnetizationCoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"MagnetizationSaturated", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"MagnetizationSaturatedMakerValue", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulus", 160000)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"PoissonRatio", 0.24)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulusX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulusY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulusZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"ShearModulusXY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"ShearModulusYZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"ShearModulusZX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G11", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G12", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G13", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G14", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G15", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G16", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G22", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G23", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G24", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G25", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G26", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G33", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G34", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G35", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G36", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G44", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G45", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G46", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G55", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G56", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G66", 0)



