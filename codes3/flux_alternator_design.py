import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
from pprint import pprint

import CrossSectInnerSalientPoleRotor, CrossSectInnerNotchedRotor, CrossSectStator, Location2D

# import pint

class flux_alternator_template(inner_rotor_motor.template_machine_as_numbers):
    ''' This is a surface mounted PM motor but it might have saliency on q-axis if alpha_rm is less than 180/p.
        就是说，允许永磁体陷入转子铁芯，只是永磁体外没有铁包裹防止永磁体飞出，而是需要额外增加碳纤维套。
    '''
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(flux_alternator_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'Flux_Alternator'
        self.name = '__Flux_Alternator'

        # 初始化搜索空间
        GP = self.d['GP'] # Geometry Parameter
        EX = self.d['EX'] # EXcitations (was OP: Other Property)
        SI = self.SI      # Specification Input dictionary (was SD)
        childGP = OrderedDict({
            # Flux_Alternator Peculiar                     Type       Name                         Value  Bounds       Calc
            "split_ratio_rotor_salient" : acmop_parameter("free",     "rotor_salient_split_ratio",  None, [None, None], lambda GP,SI:None),
            "deg_alpha_rsp"             : acmop_parameter("free",     "rotor_salient_pole_angle",   None, [None, None], lambda GP,SI:None),
            "mm_d_pm"                   : acmop_parameter("free",     "permanent_magnet_depth",   None, [None, None], lambda GP,SI:None),
            "pm_yoke_ratio"             : acmop_parameter("free",     "permanent_magnet_to_yoke_ratio",   None, [None, None], lambda GP,SI:None),
        })
        GP.update(childGP)

        # Get Analytical Design
        self.Bianchi2006(fea_config_dict, SI, GP, EX)
        GP['split_ratio_rotor_salient'].value = 0.25
        GP['deg_alpha_rsp'].value = 360 / (SI['Qs']*4) * (SI['pm']*2)
        GP['mm_d_pm'].value = 5
        GP['pm_yoke_ratio'].value = 0.9



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
        air_gap_flux_density_B = SI['guess_air_gap_flux_density_B'] # 0.9 T
        stator_tooth_flux_density_B_ds = SI['guess_stator_tooth_flux_density_B_ds'] # 1.5 T
        stator_yoke_flux_density_Bys = SI['guess_stator_yoke_flux_density_Bys']

        if SI['p'] >= 2:
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.75
            alpha_rm_over_alpha_rp = 1.0
            # stator_yoke_flux_density_Bys = 1.2
        else:
            # penalty for p=1 motor, i.e., large yoke height
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.5
            alpha_rm_over_alpha_rp = 0.75
            # stator_yoke_flux_density_Bys = 1.5

        stator_outer_diameter_Dse = 0.15 # this is related to the stator current density and should be determined by Js and power.
        sleeve_length = 3

        speed_rpm = SI['ExcitationFreqSimulated'] * 60 / SI['p'] # rpm

        rotor_outer_radius_r_or = pyrhonen_procedure_as_function.eric_specify_tip_speed_get_radius(SI['tip_speed'], speed_rpm)
        rotor_outer_diameter_Dr = rotor_outer_radius_r_or*2
        stator_inner_radius_r_is  = rotor_outer_radius_r_or + (sleeve_length+SI['minimum_mechanical_air_gap_length_mm'])*1e-3 # m (sleeve 3 mm, air gap 0.75 mm)
        stator_inner_diameter_Dis = stator_inner_radius_r_is*2
        split_ratio = stator_inner_diameter_Dis / stator_outer_diameter_Dse

        stator_yoke_height_h_ys = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis * alpha_rm_over_alpha_rp / (2*stator_yoke_flux_density_Bys * 2*SI['p'])
        # print(stator_outer_diameter_Dse, stator_inner_diameter_Dis, stator_yoke_height_h_ys)
        stator_tooth_height_h_ds = (stator_outer_diameter_Dse - stator_inner_diameter_Dis) / 2 - stator_yoke_height_h_ys
        stator_slot_height_h_ss = stator_tooth_height_h_ds
        stator_tooth_width_b_ds = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis / (stator_tooth_flux_density_B_ds* SI['Qs'])

        EX['stator_slot_area'] = stator_slot_area = np.pi/(4*SI['Qs']) * ((stator_outer_diameter_Dse - 2*stator_yoke_height_h_ys)**2 - stator_inner_diameter_Dis**2) - stator_tooth_width_b_ds * stator_tooth_height_h_ds

        slot_pitch_pps = np.pi * (stator_inner_diameter_Dis + stator_slot_height_h_ss) / SI['Qs']
        kov = 1.8 # \in [1.6, 2.0]
        EX['end_winding_length_Lew'] = end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_b_ds) + slot_pitch_pps*kov * (SI['coil_pitch_y'] - 1)

        Q = SI['Qs']
        p = SI['p']
        # STATOR
        GP['deg_alpha_st'].value         = 360/Q/2 - 2 # deg
        GP['deg_alpha_sto'].value         = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
        GP['mm_r_si'].value              = 1e3*stator_inner_radius_r_is # mm
        GP['mm_r_so'].value              = 1e3*stator_outer_diameter_Dse/2 # mm
        GP['mm_d_sto'].value              = 1 # mm
        GP['mm_d_stt'].value              = 1.5*GP['mm_d_sto'].value
        GP['mm_d_st'].value              = 1e3*(0.5*stator_outer_diameter_Dse - stator_yoke_height_h_ys) - GP['mm_r_si'].value - GP['mm_d_stt'].value  # mm
        GP['mm_d_sy'].value              = 1e3*stator_yoke_height_h_ys # mm
        GP['mm_w_st'].value              = 1e3*stator_tooth_width_b_ds # mm
        # ROTOR
        GP['mm_d_sleeve'].value          = sleeve_length
        GP['mm_d_mech_air_gap'].value    = SI['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_r_ro'].value              = 1e3*rotor_outer_radius_r_or

    def get_template_neighbor_bounds(self):
        ''' The bounds are determined around the template design.
        '''

        Q = self.SI['Qs']
        p = self.SI['p']
        # s = self.SI['no_segmented_magnets']

        GP = self.d['GP']

        deg_alpha_rsp_ini = 360 / (self.SI['Qs']*2) * (self.SI['pm']*2)

        original_template_neighbor_bounds = {
            # STATOR
            "deg_alpha_st": [ 0.35*360/Q, 0.9*360/Q],
            "mm_d_sto":     [  0.5,   5],
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.1*GP['mm_d_st'].value], # if mm_d_st is too large, the derived stator yoke can be negative
            "mm_d_sy":      [1.0*GP['mm_d_sy'].value, 1.2*GP['mm_d_sy'].value],
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],
            "split_ratio":  [0.4, 0.6], # Binder-2020-MLMS-0953@Fig.7
            "mm_d_pm":      [3,6],
            "pm_yoke_ratio":[0.9,0.98],
            # ROTOR
            "mm_d_sleeve":  [1,   2],
            "split_ratio_rotor_salient": [0.2,0.3],
            "deg_alpha_rsp": [deg_alpha_rsp_ini*0.9, deg_alpha_rsp_ini*1.1]
        }
        return original_template_neighbor_bounds

class flux_alternator_design_variant(inner_rotor_motor.variant_machine_as_objects):

    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        if x_denorm is None:
            raise

        # 初始化父类
        super(flux_alternator_design_variant, self).__init__(template, x_denorm, counter, counter_loop)

        # Give it a name
        self.name = f'ind{counter}'
        self.name += f'-redo{counter_loop}' if counter_loop > 1 else ''

        # Get geometric parameters and spec input
        GP = self.template.d['GP']
        SI = self.template.SI

        # 检查几何变量之间是否有冲突
        self.check_invalid_design(GP, SI)

        # Parts
        self.rotorCore = CrossSectInnerSalientPoleRotor.CrossSectInnerSalientPoleRotor(
                            mm_r_ro = GP['mm_r_ro'].value,
                            mm_d_sleeve = GP['mm_d_sleeve'].value,
                            split_ratio_rotor_salient = GP['split_ratio_rotor_salient'].value,
                            deg_alpha_rsp = GP['deg_alpha_rsp'].value,
                            pm = SI['pm'],
                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0))

        self.shaft = CrossSectInnerNotchedRotor.CrossSectShaft(name = 'Shaft',
                                                      notched_rotor = self.rotorCore
                                                    )

        self.stator_core = CrossSectStator.CrossSectInnerRotorStator_PMAtYoke( name = 'StatorCore',
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
                                            mm_d_pm       = GP['mm_d_pm'].value,
                                            pm_yoke_ratio = GP['pm_yoke_ratio'].value,
                                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                            )

        self.statorMagnet = CrossSectStator.CrossSectStatorMagnetAtYoke(name = 'StatorPM',
                                                                stator_core = self.stator_core
                                                                )

        # self.coils = CrossSectStator.CrossSectToroidalWiniding(name = 'Coils',
        #                                                        stator_core = self.stator_core)

        #03 Mechanical Parameters
        self.update_mechanical_parameters()

    def check_invalid_design(self, GP, SI):
        pass
