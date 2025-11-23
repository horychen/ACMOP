import winding_layout
import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np; import math
from pprint import pprint

import CrossSectInnerNotchedRotor
import CrossSectStator
import Location2D

# import pin

def derive_mm_r_ri(GP,SI):
    GP       ['mm_r_ri'].value = GP['mm_r_ro'].value - GP['mm_d_pm'].value - GP['mm_d_ri'].value
    if GP    ['mm_r_ri'].value<=0:
        logger = logging.getLogger(__name__)
        logger.error('mm_r_ri: %s, mm_r_ro: %s, mm_d_pm: %s, mm_d_ri: %s', GP['mm_r_ri'].value, GP['mm_r_ro'].value, GP['mm_d_pm'].value, GP['mm_d_ri'].value)
        logger.error('背铁太厚了 或 split_ration太小了！建议增大split_ratio是下限')
    return GP['mm_r_ri'].value 

def derive_mm_d_ri(GP,SI):
    GP       ['mm_d_ri'].value = GP['mm_r_ro'].value - GP['mm_d_pm'].value - GP['mm_r_ri'].value
    return GP['mm_d_ri'].value

class bearingless_spmsm_closedStator_template(inner_rotor_motor.template_machine_as_numbers):
    ''' This is a surface mounted PM motor but it might have saliency on q-axis if alpha_rm is less than 180/p.
        就是说，允许永磁体陷入转子铁芯，只是永磁体外没有铁包裹防止永磁体飞出，而是需要额外增加碳纤维套。
    '''
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(bearingless_spmsm_closedStator_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'ClosedSPMSM'
        self.name = '__ClosedSPMSM'

        # 初始化搜索空间
        GP = self.SI['GP'] # Geometry Parameter
        EX = self.SI['EX-user'] # EXcitations (was OP: Other Property)
        SI = self.SI      # Specification Input dictionary (was SD)
        childGP = OrderedDict({
            # SPMSM Peculiar
            "mm_d_pm"           : acmop_parameter("fixed",     "magnet_depth",                  None, [None, None], lambda GP,SI:None),
            "mm_d_ri"           : acmop_parameter("fixed",     "rotor_iron (back iron) depth",  None, [None, None], lambda GP,SI:derive_mm_d_ri(GP,SI)),
            "mm_r_ri"           : acmop_parameter("fixed",  "rotor_inner_radius",               None, [None, None], lambda GP,SI:derive_mm_r_ri(GP,SI)),
            "deg_alpha_rm"    : acmop_parameter("fixed",     "magnet_pole_span_angle",          180/SI['p'], [None, None], lambda GP,SI:None),
            # "mm_d_rp"         : acmop_parameter("fixed",     "inter_polar_iron_thickness",    None, [None, None], lambda GP,SI:None),
        })
        GP.update(childGP)

        # 把free的变量收集起来组成x_denorm_dict，它的变量的顺序是GP-user决定的
        self.x_denorm_dict, self.bounds_denorm_dict = self.set_gp_values_and_types_based_on_SI(spec_input_dict)
        logger = logging.getLogger(__name__); logger.info('Initial design gives x_denorm: %s', str(self.SI['x_denorm_dict']))
        logger = logging.getLogger(__name__); logger.info('Initial design gives x_denorm_bounds: %s', str(self.SI["bounds_denorm_dict"]))

        # Get Analytical Design
        self.PracticalInitialDesign(fea_config_dict, SI, GP, EX)
        # 定义搜索空间，determine bounds
        self.define_search_space(SI, GP)

        logger = logging.getLogger(__name__); logger.info('Initial design gives x_denorm_bounds: %s', str(self.SI["bounds_denorm_dict"]))

        # Template's Other Properties (Shared by the swarm)
        if True:
            EX = self.SI['EX'] = self.SI['EX-user']

            # WINDING Layout
            if 'Wrap_Around' in self.SI.keys():
                EX['wily'] = wily = winding_layout.winding_layout_v2(SI['DPNV_or_SEPA'], SI['Qs'], SI['p'], SI['ps'], SI['coil_pitch_y'], m=SI['m'], Wrap_Around=SI['Wrap_Around'])
            else:
                EX['wily'] = wily = winding_layout.winding_layout_v2(SI['DPNV_or_SEPA'], SI['Qs'], SI['p'], SI['ps'], SI['coil_pitch_y'], m=SI['m'])

            # THERMAL Properties
            EX['DriveW_zQ']         =            pyrhonen_procedure_as_function.get_zQ(SI, wily, GP['mm_r_si'].value*2*1e-3, GP['mm_r_ro'].value*2*1e-3, specified_mm_stack_length=EX['mm_stack_length']) # TODO:
            EX['DriveW_CurrentAmp'] = math.sqrt(2)*pyrhonen_procedure_as_function.get_stator_phase_current_rms(SI) # TODO:
            logger = logging.getLogger(__name__)
            logger.info('DriveW_CurrentAmp is initialized as: %s A (considering the specified voltage). This will be overwritten by Js-constraint later.', EX['DriveW_CurrentAmp'])
            EX['DriveW_Freq']       = EX['ExcitationFreqSimulated']
            EX['DriveW_Rs']         = 1.0 # TODO: Must be greater than zero to let JMAG work
            EX['DriveW_poles']      = SI['p']*2

        # BEARING Winding Excitation Properties
        if True:
            EX['BeariW_zQ']         = EX['DriveW_zQ']
            EX['BeariW_CurrentAmp'] = fea_config_dict['circuit.SUSPENSION_CURRENT_RATIO'] * (EX['DriveW_CurrentAmp'] / fea_config_dict['circuit.TORQUE_CURRENT_RATIO'])
            EX['BeariW_Freq']       = EX['DriveW_Freq']
            EX['BeariW_Rs']         = EX['DriveW_Rs'] * EX['BeariW_zQ'] / EX['DriveW_zQ']
            EX['BeariW_poles']      = SI['ps']*2
            EX['slot_current_utilizing_ratio'] = fea_config_dict['circuit.SUSPENSION_CURRENT_RATIO'] + fea_config_dict['circuit.TORQUE_CURRENT_RATIO'] # will be less than 1 for separate winding

    def get_x_denorm(self):
        return list(self.SI['x_denorm_dict'].values())

    def set_gp_values_and_types_based_on_SI(self, spec_input_dict):
        for key, property in spec_input_dict['GP-user'].items():

            if property['value'] is not None: self.SI['GP'][key].value = property['value']
            if property['bounds'] is not None: self.SI['GP'][key].bounds = property['bounds']

            self.SI['GP'][key].type = property['type']

            if property['type'] == 'free':
                self.SI['x_denorm_dict'][key] = property['value']
                self.SI['bounds_denorm_dict'][key] = property['bounds']

            if property['type'] == 'derived' and self.SI['GP'][key].calc is None:
                raise Exception(f'calc method is not provided for derived parameter {key}')

        return self.SI['x_denorm_dict'], self.SI['bounds_denorm_dict']

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
        stator_yoke_depth_d_sy  = Bg * np.pi * stator_inner_diameter_Dsi * get_alpha_rm_over_alpha_rp(p) / (2*Bsy * 2*p)

        # 定子齿部深度依赖于轭部深度
        stator_tooth_depth_d_st = (stator_outer_diameter_Dso - stator_inner_diameter_Dsi) *0.5 - stator_yoke_depth_d_sy
        stator_slot_depth_d_ss = stator_tooth_depth_d_st

        # 单个永磁体极面下，气隙中的磁通全部进入Qs个定子齿部（的宽度）
        stator_tooth_width_w_st = Bg * np.pi * stator_inner_diameter_Dsi / (Bst* Q)

        # 计算槽面积
        def get_stator_slot_area(Q, stator_outer_diameter_Dso, stator_yoke_depth_d_sy, stator_inner_diameter_Dsi, stator_tooth_width_w_st, stator_tooth_depth_d_st):
            return np.pi/(4*Q) * ((stator_outer_diameter_Dso - 2*stator_yoke_depth_d_sy)**2 - stator_inner_diameter_Dsi**2) - stator_tooth_width_w_st * stator_tooth_depth_d_st
        EX['stator_slot_area'] = get_stator_slot_area(Q, stator_outer_diameter_Dso, stator_yoke_depth_d_sy, stator_inner_diameter_Dsi, stator_tooth_width_w_st, stator_tooth_depth_d_st)

        # 计算端部绕组长度
        slot_pitch_pps = np.pi * (stator_inner_diameter_Dsi + stator_slot_depth_d_ss) / Q
        kov = EX['end_winding_length_factor_kov']
        EX['end_winding_length_Lew'] = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_w_st) + slot_pitch_pps*kov * (SI['coil_pitch_y'] - 1)

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
        GP['mm_d_ri'].value              = 1e3*stator_inner_radius_r_is - equivalent_air_gap_length - GP['mm_d_pm'].value - GP['mm_r_ri'].value

    def define_search_space(self, SI, GP):

        Q = SI['Qs']
        p = SI['p']

        if GP['split_ratio'].bounds is None:
            if p < 10:
                GP['split_ratio'].bounds = [0.3, 0.5]
                # "split_ratio":  [0.4, 0.6], # Binder-2020-MLMS-0953@Fig.7
                # "split_ratio":  [0.35, 0.5], # Q12p4优化的时候，轭部经常不够用，所以就把split_ratio减小——Exception: ('Error: Negative derived parameter', "acmop_parameter(type='derived', name='stator_yoke_depth', value=-1.362043443071423, bounds=[None, None], calc=<function template_machine_as_numbers.__init__.<locals>.<lambda> at 0x00000237CC403D30>)")
            else:
                GP['split_ratio'].bounds = [0.15, 0.35]

        # yoke_split_ratio = [0.2, 0.45]
        GP['mm_d_sy'].bounds = [el * (GP['mm_r_so'].value - GP['mm_r_si'].value) for el in [0.2, 0.45]]
        # tooth_split_ratio_at_middle_slot = [0.25, 0.50]
        GP['mm_w_st'].bounds = [el / Q * np.pi * (GP['mm_r_so'].value + GP['mm_r_si'].value) for el in [0.25, 0.50]]

        # attention: the bounds are determined around the template design, which means any change of the template design will lead to a change of the order the bounds.
        # original_template_neighbor_bounds = {
        #     "deg_alpha_st": [ 0.35*360/Q, 1.0*360/Q],
        #     "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],
        #     "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.1*GP['mm_d_st'].value], # if mm_d_st is too large, the derived stator yoke can be negative
        #     "mm_d_sto":     [0.2,                                         0.8], # this will influence split_ratio
        #     "mm_r_so":      [1.0*GP['mm_r_so'].value, 1.2*GP['mm_r_so'].value],
        #     "mm_d_pm":      [1, 3],
        #     "mm_d_ri":      [0.8*GP['mm_d_ri'].value,  1.2*GP['mm_d_ri'].value],
        #     # SPMSM specific
        #     "deg_alpha_rm": [0.6*360/(2*p),          1.0*360/(2*p)],
        #     "mm_d_rp":      [GP['mm_d_pm'].value/2,   GP['mm_d_pm'].value],
        # }


    """ Obsolete feature """
    def build_design_parameters_list(self):
        GP = self.SI['GP']
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

class bearingless_spmsm_closedSlot_variant(inner_rotor_motor.variant_machine_as_objects):
    ''' A variant of bearingless_spmsm_design_variant with closed slots.
    '''
    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        # 初始化父类
        super(bearingless_spmsm_closedSlot_variant, self).__init__(template, x_denorm, counter, counter_loop)
        self.x_denorm = x_denorm # for visualization only

        # Give it a name
        self.name = f'ind{counter}'
        self.name += f'-redo{counter_loop}' if counter_loop > 1 else ''

        # Get geometric parameters and spec input
        GP = self.template.SI['GP']
        SI = self.template.SI

        # 检查几何变量之间是否有冲突
        self.check_invalid_design(GP, SI)


        # 修改定子截面为闭口槽
        self.statorCore = CrossSectStator.CrossSectInnerRotorClosedSlotStator( name = 'StatorCore',
                                            mm_d_stt = SI['GP']['mm_d_stt'].value,
                                            mm_r_si = SI['GP']['mm_r_si'].value,
                                            mm_d_st = SI['GP']['mm_d_st'].value,
                                            mm_d_sy = SI['GP']['mm_d_sy'].value,
                                            mm_w_st = SI['GP']['mm_w_st'].value,
                                            Q = SI['Qs'],
                                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                            )

        self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                    stator_core = self.statorCore)

        # Parts
        print('[bearingless_spmsm_closedSlot_design.py] Building parts for variant:', self.name)
        deg_alpha_rs         = GP['deg_alpha_rm'].value if 'deg_alpha_rs' not in GP.keys() else GP['deg_alpha_rs'].value
        mm_d_rs              = 0.0                      if 'mm_d_rs'      not in GP.keys() else GP['mm_d_rs'].value
        mm_d_rp              = GP['mm_d_pm'].value      if 'mm_d_rp'      not in GP.keys() else GP['mm_d_rp'].value
        no_segmented_magnets = 1                        if 'no_segmented_magnets' not in GP.keys() else template.SI['no_segmented_magnets']
        self.rotorCore = CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                            name = 'NotchedRotor',
                            mm_d_pm      = GP['mm_d_pm'].value,
                            deg_alpha_rm = GP['deg_alpha_rm'].value, # angular span of the pole: class type DimAngular
                            deg_alpha_rs = deg_alpha_rs, # segment span: class type DimAngular
                            mm_d_ri      = GP['mm_d_ri'].value, # rotor iron thickness: class type DimLinear
                            mm_r_ri      = GP['mm_r_ri'].value, # inner radius of rotor: class type DimLinear
                            mm_d_rp      = mm_d_rp, # interpolar iron thickness: class type DimLinear
                            mm_d_rs      = mm_d_rs, # inter segment iron thickness: class type DimLinear
                            p = template.SI['p'], # Set pole-pairs to 2
                            s = no_segmented_magnets, # Set magnet segments/pole to 4
                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0))

        self.shaft = CrossSectInnerNotchedRotor.CrossSectShaft(name = 'Shaft',
                                                      notched_rotor = self.rotorCore
                                                    )

        self.rotorMagnet = CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
                                                      notched_rotor = self.rotorCore
                                                    )

        # self.statorCore = CrossSectStator.CrossSectInnerRotorStator( name = 'StatorCore',
        #                                     deg_alpha_st = GP['deg_alpha_st'].value, #40,
        #                                     deg_alpha_sto = GP['deg_alpha_sto'].value, #20,
        #                                     mm_r_si = GP['mm_r_si'].value,
        #                                     mm_d_sto = GP['mm_d_sto'].value,
        #                                     mm_d_stt = GP['mm_d_stt'].value,
        #                                     mm_d_st = GP['mm_d_st'].value,
        #                                     mm_d_sy = GP['mm_d_sy'].value,
        #                                     mm_w_st = GP['mm_w_st'].value,
        #                                     mm_r_st = 0.0, # =0
        #                                     mm_r_sf = 0.0, # =0
        #                                     mm_r_sb = 0.0, # =0
        #                                     Q = template.SI['Qs'],
        #                                     location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
        #                                     )

        # self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
        #                                             stator_core = self.statorCore)

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


        #03 Mechanical Parameters
        self.update_mechanical_parameters()



        if True:
            ''' This was moved from JMAG.preProcess()
            '''
            # Implementation of id=0 control:
            #   After rotate the rotor by half the inter-pole notch span, The d-axis initial position is at pole pitch angle divided by 2.
            #   The U-phase current is sin(omega_syn*t) = 0 at t=0 and requires the d-axis to be at the winding phase axis (to obtain id=0 control)
            deg_pole_span = 180/SI['p']
            wily = self.template.SI['EX']['wily']
            #                                                              inter-pole notch (0.5 for half)         rotate to x-axis    winding placing bias (half adjacent slot angle)      reverse north and south pole to make torque positive.
            logger = logging.getLogger(__name__)
            logger.info('[PMSM JMAG] InitialRotationAngle : %s %s %s %s', (deg_pole_span-GP['deg_alpha_rm'].value)*0.5, - deg_pole_span*0.5, + wily.deg_winding_U_phase_phase_axis_angle,  + deg_pole_span)
            logger.info('[PMSM JMAG] InitialRotationAngle = %s deg', (deg_pole_span-GP['deg_alpha_rm'].value)*0.5  - deg_pole_span*0.5  + wily.deg_winding_U_phase_phase_axis_angle   + deg_pole_span)
            self.InitialRotationAngle = (deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + wily.deg_winding_U_phase_phase_axis_angle     + deg_pole_span
            logger.info('[PMSM JMAG] InitialRotationAngle = %s deg', self.InitialRotationAngle)

        self.boolCustomizedCircuit = False

    def check_invalid_design(self, GP, SI):
        pass

        # 不合理的变量选择（mm_d_rp）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        # if GP['mm_d_rp'].value > GP['mm_d_pm'].value:
        #     GP['mm_d_rp'].value            = GP['mm_d_pm'].value
        #     # free_variables[11] = free_variables[6]
        #     msg = '[Warning from bearingless_spmsm_design.py]: Inter-pole notch depth mm_d_rp cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rp is set to mm_d_pm.'
        #     logger = logging.getLogger(__name__)
        #     logger.warning(msg)

        # 不合理的变量选择（mm_d_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        # if GP['mm_d_rs'].value > GP['mm_d_pm'].value:
        #     GP['mm_d_rs'].value            = GP['mm_d_pm'].value
        #     # free_variables[12] = free_variables[6]
        #     msg = '[Warning from bearingless_spmsm_design.py]: Inter-segment notch depth mm_d_rs cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rs is set to mm_d_pm.'
        #     logger = logging.getLogger(__name__)
        #     logger.warning(msg)

        # 不合理的变量选择（deg_alpha_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        # if not (GP['deg_alpha_rs'].value > GP['deg_alpha_rm'].value/SI['no_segmented_magnets']):
        #     GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value/SI['no_segmented_magnets']
        #     msg = '[Warning from bearingless_spmsm_design.py]: deg_alpha_rs cannot be larger than deg_alpha_rm/s. Note deg_alpha_rs is set to deg_alpha_rm/s.'
        #     logger = logging.getLogger(__name__)
        #     logger.warning(msg)

        # 如果没有永磁体分段，那么alpha_rs应该等于alpha_rm。
        # if SI['no_segmented_magnets'] == 1:
        #     GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value # raise Exception('Invalid alpha_rs. Check that it is equal to alpha_rm for s=1')
        #     GP['mm_d_rs'].value = 0 # raise Exception('Invalid d_rs. Check that it is equal to 0 for s =1')

        # This is all we need

# class bearingless_spmsm_closedSlot_variant(inner_rotor_motor.variant_machine_as_objects):
#     ''' A variant of bearingless_spmsm_design_variant with closed slots.
#     '''
#     def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
#         # 初始化父类
#         super(bearingless_spmsm_closedSlot_variant, self).__init__(template, x_denorm, counter, counter_loop)

#         SI = self.template.SI

#         # 修改定子截面为闭口槽
#         self.statorCore = CrossSectStator.CrossSectInnerRotorClosedSlotStator( name = 'StatorCore',
#                                             # deg_alpha_st = SI['GP']['deg_alpha_st'].value, #40,
#                                             # deg_alpha_sto = SI['GP']['deg_alpha_sto'].value, #20,
#                                             mm_r_si = SI['GP']['mm_r_si'].value,
#                                             # mm_d_sto = SI['GP']['mm_d_sto'].value,
#                                             # mm_d_stt = SI['GP']['mm_d_stt'].value,
#                                             mm_d_st = SI['GP']['mm_d_st'].value,
#                                             mm_d_sy = SI['GP']['mm_d_sy'].value,
#                                             mm_w_st = SI['GP']['mm_w_st'].value,
#                                             Q = SI['Qs'],
#                                             location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
#                                             )

#         self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
#                                                     stator_core = self.statorCore)

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



