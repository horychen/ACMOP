
# 修改此文件为连续极电机：

import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
from pprint import pprint

import CrossSectInnerNotchedRotor, CrossSectInnerConsequentPoleRotor
import CrossSectStator
import Location2D

# import pin

def derive_mm_r_ri(GP,SI):
    GP       ['mm_r_ri'].value = GP['mm_r_ro'].value - GP['mm_d_pm'].value - GP['mm_d_ri'].value
    if GP    ['mm_r_ri'].value<=0:
        print()
        print(GP       ['mm_r_ri'].value, GP['mm_r_ro'].value , GP['mm_d_pm'].value , GP['mm_d_ri'].value)
        print('背铁太厚了 或 split_ration太小了！建议增大split_ratio是下限')
        print()
    return GP['mm_r_ri'].value 

class bearingless_consequentPole_template(inner_rotor_motor.template_machine_as_numbers):
    ''' This is a surface mounted PM motor but it might have saliency on q-axis if alpha_rm is less than 180/p.
        就是说，允许永磁体陷入转子铁芯，只是永磁体外没有铁包裹防止永磁体飞出，而是需要额外增加碳纤维套。
    '''
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(bearingless_consequentPole_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'CPPM'
        self.name = '__CPPM'

       # 初始化搜索空间
        GP = self.d['GP'] # Geometry Parameter
        EX = self.d['EX'] # EXcitations (was OP: Other Property)
        SI = self.SI      # Specification Input dictionary (was SD)
        childGP = OrderedDict({
            # CPPM Peculiar
            "mm_d_pm"           : acmop_parameter("free",     "magnet_depth",                  None, [None, None], lambda GP,SI:None),
            "mm_d_ri"           : acmop_parameter("free",     "rotor_iron (back iron) depth",  None, [None, None], lambda GP,SI:None),
            "deg_alpha_rm"      : acmop_parameter("free",     "magnet_pole_span_angle",        None, [None, None], lambda GP,SI:None),
            "mm_d_rp"           : acmop_parameter("free",     "inter_polar_iron_thickness",    None, [None, None], lambda GP,SI:None),
        #    "deg_alpha_rs"      : acmop_parameter("free",     "magnet_segment_span_angle",     None, [None, None], lambda GP,SI:None),
        #    "mm_d_rs"           : acmop_parameter("free",     "inter_segment_iron_thickness",  None, [None, None], lambda GP,SI:None),
            "mm_r_ri"           : acmop_parameter("derived",  "rotor_inner_radius",            None, [None, None], lambda GP,SI:derive_mm_r_ri(GP,SI)),
        })
        GP.update(childGP)

        # Get Analytical Design
        self.Bianchi2006(fea_config_dict, SI, GP, EX)

        # 定义搜索空间，determine bounds
        original_template_neighbor_bounds = self.get_template_neighbor_bounds()
        self.bounds_denorm = self.define_search_space(GP, original_template_neighbor_bounds)

        # Template's Other Properties (Shared by the swarm)
        EX = self.get_other_properties_after_geometric_parameters_are_initialized(GP, SI)
        # BEARING Winding Excitation Properties
        if True:
            EX['BeariW_zQ']         = EX['DriveW_zQ']
            EX['BeariW_CurrentAmp'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] * (EX['DriveW_CurrentAmp'] / fea_config_dict['TORQUE_CURRENT_RATIO'])
            EX['BeariW_Freq']       = EX['DriveW_Freq']
            EX['BeariW_Rs']         = EX['DriveW_Rs'] * EX['BeariW_zQ'] / EX['DriveW_zQ']
            EX['BeariW_poles']      = SI['ps']*2
            EX['slot_current_utilizing_ratio'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] + fea_config_dict['TORQUE_CURRENT_RATIO'] # will be less than 1 for separate winding

    def Bianchi2006(self, fea_config_dict, SI, GP, EX):

        # inputs: air gap flux density
        air_gap_flux_density_Bg = SI['guess_air_gap_flux_density_Bg'] # 0.9 T
        stator_tooth_flux_density_Bst = SI['guess_stator_tooth_flux_density_Bst'] # 1.5 T
        stator_yoke_flux_density_Bsy = SI['guess_stator_yoke_flux_density_Bsy']

        if SI['p'] >= 2:
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.75
            alpha_rm_over_alpha_rp = 1.0
            # stator_yoke_flux_density_Bsy = 1.2
        else:
            # penalty for p=1 motor, i.e., large yoke height
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.5
            alpha_rm_over_alpha_rp = 0.75
            # stator_yoke_flux_density_Bsy = 1.5

        # ureg = pint.UnitRegistry()  # 0.225* ureg.meter
        stator_outer_diameter_Dse = 0.390 # this is related to the stator current density and should be determined by Js and power.
        sleeve_length = 1

        speed_rpm = SI['ExcitationFreqSimulated'] * 60 / SI['p'] # rpm

        rotor_outer_radius_r_or = pyrhonen_procedure_as_function.eric_specify_tip_speed_get_radius(SI['tip_speed'], speed_rpm)
        rotor_outer_diameter_Dr = rotor_outer_radius_r_or*2
        stator_inner_radius_r_is  = rotor_outer_radius_r_or + (sleeve_length+SI['minimum_mechanical_air_gap_length_mm'])*1e-3 # m (sleeve 3 mm, air gap 0.75 mm)
        stator_inner_diameter_Dis = stator_inner_radius_r_is*2
        split_ratio = stator_inner_diameter_Dis / stator_outer_diameter_Dse

        stator_yoke_height_h_ys = air_gap_flux_density_Bg * np.pi * stator_inner_diameter_Dis * alpha_rm_over_alpha_rp / (2*stator_yoke_flux_density_Bsy *2 *SI['p'])
        # print(stator_outer_diameter_Dse, stator_inner_diameter_Dis, stator_yoke_height_h_ys)
        stator_tooth_height_h_ds = (stator_outer_diameter_Dse - stator_inner_diameter_Dis) / 2 - stator_yoke_height_h_ys
        stator_slot_height_h_ss = stator_tooth_height_h_ds
        stator_tooth_width_b_ds = air_gap_flux_density_Bg * np.pi * stator_inner_diameter_Dis / (stator_tooth_flux_density_Bst* SI['Qs'])

        EX['stator_slot_area'] = stator_slot_area = np.pi/(4*SI['Qs']) * ((stator_outer_diameter_Dse - 2*stator_yoke_height_h_ys)**2 - stator_inner_diameter_Dis**2) - stator_tooth_width_b_ds * stator_tooth_height_h_ds

        slot_pitch_pps = np.pi * (stator_inner_diameter_Dis + stator_slot_height_h_ss) / SI['Qs']
        kov = 1.8 # \in [1.6, 2.0]
        EX['end_winding_length_Lew'] = end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_b_ds) + slot_pitch_pps*kov * (SI['coil_pitch_y'] - 1)

        Q = SI['Qs']
        p = SI['p']
        # STATOR
        GP['deg_alpha_st'].value         = 360/Q - 2 # deg
        GP['deg_alpha_sto'].value         = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
        GP['mm_r_si'].value              = 1e3*stator_inner_radius_r_is # mm
        GP['mm_r_so'].value              = 1e3*stator_outer_diameter_Dse/2 # mm
        GP['mm_d_sto'].value              = 3 # mm
        GP['mm_d_stt'].value              = 1.5*GP['mm_d_sto'].value
        GP['mm_d_st'].value              = 1e3*(0.5*stator_outer_diameter_Dse - stator_yoke_height_h_ys) - GP['mm_r_si'].value - GP['mm_d_stt'].value  # mm
        # print(f"{GP['mm_d_st'].value=}")
        # print (f"{1e3*stator_outer_diameter_Dse=}")
        # print(f"{1e3*stator_yoke_height_h_ys=}")
        # print(f"{speed_rpm=}")
        # print(f"{rotor_outer_radius_r_or=}")
        # print(f"{stator_tooth_height_h_ds=}")
        # print(f"{GP['mm_r_si'].value=}")
        # print(f"{GP['mm_r_si'].value=}")
        # print(f"{GP['mm_r_si'].value=}")
        # print(f"{air_gap_flux_density_Bg =}")
        # print(f"{stator_inner_diameter_Dis=}")
        # print(f"{alpha_rm_over_alpha_rp =}")
        # print(f"{2*stator_yoke_flux_density_Bsy =}")
        # print(f"{2*SI['p']=}")
        # print (f"{GP['mm_d_stt'].value =}")
        # quit()
        GP['mm_d_sy'].value              = 1e3*stator_yoke_height_h_ys # mm
        GP['mm_w_st'].value              = 1e3*stator_tooth_width_b_ds # mm
        # ROTOR
        GP['mm_d_sleeve'].value          = sleeve_length
        GP['mm_d_mech_air_gap'].value    = SI['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_d_pm'].value              = 10  # mm
        GP['mm_d_ri'].value              = 1e3*ROTOR_STATOR_YOKE_HEIGHT_RATIO*stator_yoke_height_h_ys # TODO：This ratio (0.75) is epirically specified
        GP['mm_r_ro'].value              = 1e3*rotor_outer_radius_r_or
        GP['mm_r_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_d_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value
        # SPMSM specific
        GP['deg_alpha_rm'].value         = 0.95*360/(2*p) # deg
        GP['mm_d_rp'].value              = 10  # mm
    #    GP['deg_alpha_rs'].value         = 0.975*GP['deg_alpha_rm'].value / SI['no_segmented_magnets']
        # GP['mm_d_rs'].value              = 0.20*GP['mm_d_rp'].value # d_pm > d_rp and d_pm > d_rs

        # Those are some obsolete variables that are convenient to have.
        # template.Radius_OuterStatorYoke = spec_geometry_dict['Radius_OuterStatorYoke'] = 1e3*0.5*stator_outer_diameter_Dse # mm
        # template.Radius_OuterRotor      = spec_geometry_dict['Radius_OuterRotor'] = 1e3*rotor_outer_radius_r_or # mm

        # Those variables are for PMSM convenience
        # template.rotor_steel_outer_radius = spec_geometry_dict['rotor_steel_outer_radius'] = template.Radius_OuterRotor

        # required_torque = SI['mec_power']/(2*np.pi*speed_rpm)*60
        # rotor_volume_Vr = required_torque/(2*SI['TangentialStress'])
        # template.required_torque = required_torque

    def get_template_neighbor_bounds(self):
        ''' The bounds are determined around the template design.
        '''

        Q = self.SI['Qs']
        p = self.SI['p']
        # s = self.SI['no_segmented_magnets']

        GP = self.d['GP']

        original_template_neighbor_bounds = {
            # STATOR
            "deg_alpha_st": [ 0.35*360/Q, 0.9*360/Q],
            "mm_d_sto":      [  0.5,   5],
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.1*GP['mm_d_st'].value], # if mm_d_st is too large, the derived stator yoke can be negative
            # "mm_r_so":      [1.0*GP['mm_r_so'].value, 1.2*GP['mm_r_so'].value],
            "mm_d_sy":      [1.0*GP['mm_d_sy'].value, 1.2*GP['mm_d_sy'].value],
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],
            # ROTOR
            "mm_d_sleeve":  [0.5,   2],
            # "split_ratio":  [0.4, 0.6], # Binder-2020-MLMS-0953@Fig.7
            "split_ratio":  [0.35, 0.5], # Q12p4优化的时候，轭部经常不够用，所以就把split_ratio减小——Exception: ('Error: Negative derived parameter', "acmop_parameter(type='derived', name='stator_yoke_depth', value=-1.362043443071423, bounds=[None, None], calc=<function template_machine_as_numbers.__init__.<locals>.<lambda> at 0x00000237CC403D30>)")
            "mm_d_pm":      [2.5, 7],
            "mm_d_ri":      [0.8*GP['mm_d_ri'].value,  1.2*GP['mm_d_ri'].value],
            # SPMSM specific
            "deg_alpha_rm": [0.6*360/(2*p),          1.0*360/(2*p)],
            "mm_d_rp":      [2.5,   6],
            # "deg_alpha_rs": [0.8*360/(2*p)/s,        0.975*360/(2*p)/s],
            # "mm_d_rs":      [2.5,   6]
        }
        return original_template_neighbor_bounds

    """ Obsolete feature """
    def build_design_parameters_list(self):
        GP = self.d['GP']
        SI = self.SI
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

class bearingless_consequentPole_design_variant(inner_rotor_motor.variant_machine_as_objects):

    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        if x_denorm is None:
            raise

        # 初始化父类
        super(bearingless_consequentPole_design_variant, self).__init__(template, x_denorm, counter, counter_loop)

        # Give it a name
        self.name = f'ind{counter}'
        self.name += f'-redo{counter_loop}' if counter_loop > 1 else ''

        # Get geometric parameters and spec input
        GP = self.template.d['GP']
        SI = self.template.SI

        # 检查几何变量之间是否有冲突
        self.check_invalid_design(GP, SI)

        # Parts
        self.rotorCore = CrossSectInnerConsequentPoleRotor.CrossSectConsequentPoleRotor(
                            name = 'ConsequentPoleRotor',
                            mm_d_pm      = GP['mm_d_pm'].value,
                            deg_alpha_rm = GP['deg_alpha_rm'].value, # angular span of the pole: class type DimAngular
                            # deg_alpha_rs = GP['deg_alpha_rs'].value, # segment span: class type DimAngular
                            mm_d_ri      = GP['mm_d_ri'].value, # rotor iron thickness: class type DimLinear
                            mm_r_ri      = GP['mm_r_ri'].value, # inner radius of rotor: class type DimLinear
                            mm_d_rp      = GP['mm_d_rp'].value, # interpolar iron thickness: class type DimLinear
                            # mm_d_rs      = GP['mm_d_rs'].value, # inter segment iron thickness: class type DimLinear
                            p = template.SI['p'], # Set pole-pairs to 2
                            s = template.SI['no_segmented_magnets'], # Set magnet segments/pole to 4
                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0))

        self.shaft = CrossSectInnerConsequentPoleRotor.CrossSectConsequentPoleShaft(name = 'Shaft',
                                                      ConsequentPole_rotor = self.rotorCore
                                                    )

        self.rotorMagnet = CrossSectInnerConsequentPoleRotor.CrossSectConsequentPoleMagnet( name = 'ConsequentPole',
                                                      ConsequentPole_rotor = self.rotorCore
                                                    )

        self.stator_core = CrossSectStator.CrossSectInnerRotorStator( name = 'StatorCore',
                                            deg_alpha_st = GP['deg_alpha_st'].value, #40,
                                            deg_alpha_sto = GP['deg_alpha_sto'].value, #20,
                                            mm_r_si = GP['mm_r_si'].value,
                                            mm_d_sto = GP['mm_d_sto'].value,
                                            mm_d_stt = GP['mm_d_stt'].value,
                                            mm_d_st = GP['mm_d_st'].value,
                                            mm_d_sy = GP['mm_d_sy'].value,
                                            mm_w_st = GP['mm_w_st'].value,
                                            mm_r_st = 0.0, # =0
                                            mm_r_sf = 0.0, # =0
                                            mm_r_sb = 0.0, # =0
                                            Q = template.SI['Qs'],
                                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                            )

        self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                    stator_core = self.stator_core)

        self.sleeve = CrossSectInnerNotchedRotor.CrossSectSleeve(
                            name = 'Sleeve',
                            notched_magnet = self.rotorMagnet,
                            d_sleeve = GP['mm_d_sleeve'].value
                            )

        # p = SI['p']
        # s = SI['no_segmented_magnets']
        # Qs = SI['Qs']
        #                            #   轴 转子 永磁体  护套 定子 绕组
        # self.number_of_parts_in_JMAG = 1 + 1 + p*2*s + 1 + 1 + Qs*2

        # 暂时有点乱 先这样
        # self.Rotation_Axis = -1
        #03 Mechanical Parameters
        self.update_mechanical_parameters()



        if True:
            ''' This was moved from JMAG.preProcess()
            '''
            # Implementation of id=0 control:
            #   After rotate the rotor by half the inter-pole notch span, The d-axis initial position is at pole pitch angle divided by 2.
            #   The U-phase current is sin(omega_syn*t) = 0 at t=0 and requires the d-axis to be at the winding phase axis (to obtain id=0 control)
            deg_pole_span = 180/SI['p']
            wily = self.template.d['EX']['wily']
            #                                                              inter-pole notch (0.5 for half)         rotate to x-axis    winding placing bias (half adjacent slot angle)      reverse north and south pole to make torque positive.
            print('[bearingless_cppm_design.py] [CPPM JMAG] InitialRotationAngle :',(deg_pole_span-GP['deg_alpha_rm'].value)*0.5, - deg_pole_span*0.5, + wily.deg_winding_U_phase_phase_axis_angle,  + deg_pole_span)
            print('[bearingless_cppm_design.py] [CPPM JMAG] InitialRotationAngle =',(deg_pole_span-GP['deg_alpha_rm'].value)*0.5  - deg_pole_span*0.5  + wily.deg_winding_U_phase_phase_axis_angle   + deg_pole_span, 'deg')
            self.InitialRotationAngle = (deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + wily.deg_winding_U_phase_phase_axis_angle     + deg_pole_span

        self.boolCustomizedCircuit = False

    def check_invalid_design(self, GP, SI):

        # 不合理的变量选择（mm_d_rp）会导致：一个变量的取值范围是受到另一个变量的取值影响的。

        if GP['mm_d_rp'].value > GP['mm_d_pm'].value:
            GP['mm_d_rp'].value            = GP['mm_d_pm'].value
            # free_variables[11] = free_variables[6]

            msg = '[Warning from bearingless_spmsm_design.py]: Inter-pole notch depth mm_d_rp cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rp is set to mm_d_pm.'
            print(msg)
            logger = logging.getLogger(__name__).warn(msg)

        # 不合理的变量选择（mm_d_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        # if GP['mm_d_rs'].value > GP['mm_d_pm'].value:
        #     GP['mm_d_rs'].value            = GP['mm_d_pm'].value
        #     # free_variables[12] = free_variables[6]

        #     msg = '[Warning from bearingless_spmsm_design.py]: Inter-segment notch depth mm_d_rs cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rs is set to mm_d_pm.'
        #     print(msg)
        #     logger = logging.getLogger(__name__).warn(msg)

        # # 不合理的变量选择（deg_alpha_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        # if not (GP['deg_alpha_rs'].value > GP['deg_alpha_rm'].value/SI['no_segmented_magnets']):
        #     GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value/SI['no_segmented_magnets']
        #     msg = '[Warning from bearingless_spmsm_design.py]: deg_alpha_rs cannot be larger than deg_alpha_rm/s. Note deg_alpha_rs is set to deg_alpha_rm/s.'
        #     print(msg)
        #     logger = logging.getLogger(__name__).warn(msg)

        # 如果没有永磁体分段，那么alpha_rs应该等于alpha_rm。
        # if SI['no_segmented_magnets'] == 1:
        #     GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value # raise Exception('Invalid alpha_rs. Check that it is equal to alpha_rm for s=1')
        #     GP['mm_d_rs'].value = 0 # raise Exception('Invalid d_rs. Check that it is equal to 0 for s =1')
        # This is all we need

