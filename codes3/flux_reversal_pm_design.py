import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
from pprint import pprint

import CrossSectInnerSalientPoleRotor, CrossSectInnerNotchedRotor, CrossSectStator, Location2D

def derive_mm_d_pm(GP,SI):
    GP       ['mm_d_pm'].value = GP['mm_r_si'].value * np.sin(0.5*GP['deg_alpha_pm_at_airgap'].value/180*np.pi) * 2
    return GP['mm_d_pm'].value

def derive_deg_alpha_pm_at_airgap(GP,SI):
    GP       ['deg_alpha_pm_at_airgap'].value = np.arcsin(0.5*GP['mm_d_pm'].value/GP['mm_r_si'].value) / np.pi*180 * 2 # option 1
    # GP       ['deg_alpha_pm_at_airgap'].value = GP['mm_d_pm'].value / (2*np.pi*GP['mm_r_si'].value) * 360 # option 2 (aproximation line to arc)
    return GP['deg_alpha_pm_at_airgap'].value

def derive_mm_d_ri(GP,SI):
    GP       ['mm_d_ri'].value = (1-GP['split_ratio_rotor_salient'].value) * GP['mm_r_ro'].value - GP['mm_r_ri'].value
    return GP['mm_d_ri'].value

def derive_mm_d_sy_overwrite(GP,SI):
    k2 = (1 / GP['PM_to_tooth_width_ratio'].value - 1 ) * 0.5
    print('[flux_switching.py] k2 =', k2)
    GP       ['mm_d_sy'].value = GP['FRPM_k3'].value * k2*GP['mm_d_pm'].value
    return GP['mm_d_sy'].value

class FRPM_template(inner_rotor_motor.template_machine_as_numbers):
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(FRPM_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'FRPM'
        self.name = '__FRPM'

        # 初始化搜索空间
        GP = self.d['GP'] # Geometry Parameter
        EX = self.d['EX'] # EXcitations (was OP: Other Property)
        SI = self.SI      # Specification Input dictionary (was SD)
        childGP = OrderedDict({
            # FRPM Peculiar                     Type       Name                         Value  Bounds       Calc
            "split_ratio_rotor_salient" : acmop_parameter("free",     "rotor_salient_split_ratio",  None, [None, None], lambda GP,SI:None),
            "deg_alpha_rsp"             : acmop_parameter("free",     "rotor_salient_pole_angle",   None, [None, None], lambda GP,SI:None),
            "PM_to_tooth_width_ratio"   : acmop_parameter("free",     "PM_to_tooth_width_ratio",    None, [None, None], lambda GP,SI:None),
            "mm_d_pm"                   : acmop_parameter("fixed",    "permanent_magnet_depth",     None, [None, None], lambda GP,SI:None), #derive_mm_d_pm(GP,SI)),
            "mm_d_air_pm"               : acmop_parameter("fixed",    "permanent_magnet_air_depth", None, [None, None], lambda GP,SI:None),
            "deg_alpha_pm_at_airgap"    : acmop_parameter("derived",  "permanent_magnet_span_angle_at_airgap",   None, [None, None], lambda GP,SI:derive_deg_alpha_pm_at_airgap(GP,SI)),
            "mm_difference_pm_yoke"     : acmop_parameter("fixed",    "permanent_magnet_to_yoke_ratio",   None, [None, None], lambda GP,SI:None),
            "FRPM_k3"                   : acmop_parameter("fixed",    "FRPM k3 (yoke to U core iron leg ratio)",   None, [None, None], lambda GP,SI:None),
            "mm_r_ri"                   : acmop_parameter("fixed",    "rotor_inner_radius (shaft_radius)", None, [None, None], lambda GP,SI:None),
            "mm_d_ri"                   : acmop_parameter("derived",  "rotor_iron (back iron) depth",  None, [None, None], lambda GP,SI:derive_mm_d_ri(GP,SI)),
            "mm_d_sy"                   : acmop_parameter("derived",  "stator_yoke_depth (OVERWRITE)", None, [None, None], lambda GP,SI:derive_mm_d_sy_overwrite(GP,SI)),
        })
        GP.update(childGP)

        # Get Analytical Design
        # self.MyCustomizedInitialDesignScript(fea_config_dict, SI, GP, EX)
        self.MyCustomizedInitialDesignScript_Version2(fea_config_dict, SI, GP, EX)
        GP['split_ratio_rotor_salient'].value = 0.15
        GP['mm_difference_pm_yoke'].value = 0.5 # mm
        GP['mm_r_ri'].value = SI['mm_radius_shaft']
        GP['mm_d_air_pm'].value = 0

        # 定义搜索空间，determine bounds
        original_template_neighbor_bounds = self.get_template_neighbor_bounds()
        self.bounds_denorm = self.define_search_space(GP, original_template_neighbor_bounds)

        # Template's Other Properties (Shared by the swarm)
        EX = self.get_other_properties_after_geometric_parameters_are_initialized(GP, SI, SI['mm_stack_length'])
        # EX['RotorPoleNumber'] = SI['pm']*2
        # BEARING Winding Excitation Properties
        if True:
            EX['BeariW_zQ']         = EX['DriveW_zQ']
            # EX['BeariW_CurrentAmp'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] * (EX['DriveW_CurrentAmp'] / fea_config_dict['TORQUE_CURRENT_RATIO'])
            EX['BeariW_Freq']       = EX['DriveW_Freq']
            EX['BeariW_Rs']         = EX['DriveW_Rs'] * EX['BeariW_zQ'] / EX['DriveW_zQ']
            EX['BeariW_poles']      = SI['ps']*2
            EX['slot_current_utilizing_ratio'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] + fea_config_dict['TORQUE_CURRENT_RATIO'] # will be less than 1 for separate winding

    def MyCustomizedInitialDesignScript_Version2(self, fea_config_dict, SI, GP, EX):
        mm_stack_length = SI['mm_stack_length']
        Pa_TangentialStress = SI['TangentialStress']
        m  = SI['m']
        Qs = SI['Qs']
        pm = SI['pm']
        pa = SI['pa']
        mec_power = SI['mec_power']
        speed_rpm = SI['ExcitationFreqSimulated'] / pm * 60

        required_torque = mec_power/(2*np.pi*speed_rpm/60)

        # (6.2)
        RotorActiveVolumn_Vr = required_torque / (2*Pa_TangentialStress)
        RotorActiveCrossSectArea_Sr = RotorActiveVolumn_Vr / (mm_stack_length*1e-3)
        RotorOuterRadius_r_or = np.sqrt(RotorActiveCrossSectArea_Sr/np.pi)
        mm_r_ro = RotorOuterRadius_r_or*1e3

        # print(mm_r_ro)
        # quit()

        aspect_ratio__rotor_axial_to_diameter_ratio = 2*mm_r_ro/mm_stack_length

        # mm_air_gap_length = 2 # Gruber Habil: [3,4] mm
        mm_air_gap_length = 3 # Gruber Habil: [3,4] mm
        mm_r_si     = mm_r_ro + mm_air_gap_length
        mm_r_airgap = mm_r_ro + mm_air_gap_length*0.5

        # note this is only valid for 12s/10pp motor
        if pm == 10 or pm == 20:
            rotor_to_stator_tooth_width_ratio = 1.4 # 2010 诸自强 1点4倍的来源 Analysis of Electromagnetic Performance of Flux
        elif pm == 14 or pm == 28:
            rotor_to_stator_tooth_width_ratio = 1.0
        else:
            rotor_to_stator_tooth_width_ratio = 1.4
            # raise Exception('Not supported pm value:', pm)

        if Qs == 6:
            rotor_to_stator_tooth_width_ratio = 1.0

        if pm == 28 and Qs == 12:
            rotor_to_stator_tooth_width_ratio = 1.0
        if pm == 20 and Qs == 12:
            rotor_to_stator_tooth_width_ratio = 1.0
        k4 = rotor_to_stator_tooth_width_ratio

        # mm_stator_tooth_width_w_UCoreWidth = 2*np.pi*mm_r_si / Qs / 4
        if Qs < 18:
            mm_d_pm = 5 #mm_stator_tooth_width_w_UCoreWidth
        else:
            mm_d_pm = 3

        k2 = 1.2
        PM_to_tooth_width_ratio = 1 / (1+2*k2)
        mm_stator_tooth_width_w_UCoreWidth = k2*mm_d_pm
        mm_w_st = mm_d_pm + 2*mm_stator_tooth_width_w_UCoreWidth
        mm_w_rt = k4 * mm_stator_tooth_width_w_UCoreWidth # Empirical: slightly wider (5.18) Habil-Gruber

        k1 = split_ratio = 0.7
        mm_r_so = mm_r_ro / split_ratio

        k3 = 1.2
        mm_d_sy = k3*mm_stator_tooth_width_w_UCoreWidth
        mm_d_st = mm_r_so - (mm_r_si + mm_d_sy)

        弧长 = mm_w_rt
        deg_alpha_rsp = 弧长 / mm_r_ro /np.pi*180

        弧长 = mm_w_st
        deg_alpha_st = 弧长 / mm_r_si /np.pi*180 + 2 # 2 deg of tooth tip

        # STATOR
        GP['deg_alpha_st'].value         = deg_alpha_st
        GP['deg_alpha_sto'].value        = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
        GP['mm_d_sto'].value             = 2 # mm
        GP['mm_d_stt'].value             = 3 * GP['mm_d_sto'].value # 1.5* # reduce tooth tip saturation
        GP['mm_r_si'].value              = mm_r_si
        GP['mm_r_so'].value              = mm_r_so
        GP['mm_d_st'].value              = mm_d_st - GP['mm_d_stt'].value
        GP['mm_d_sy'].value              = mm_d_sy
        # print(mm_r_so, mm_d_sy)
        # print(mm_r_so, mm_d_sy)
        # print(mm_r_so, mm_d_sy)
        GP['FRPM_k3'].value              = k3
        GP['PM_to_tooth_width_ratio'].value = PM_to_tooth_width_ratio
        GP['mm_w_st'].value              = mm_w_st
        GP['mm_d_pm'].value              = mm_d_pm
        GP['deg_alpha_pm_at_airgap'].value = mm_d_pm / (2*np.pi*mm_r_si) * 360
        # ROTOR
        GP['mm_d_sleeve'].value          = mm_air_gap_length - SI['minimum_mechanical_air_gap_length_mm']
        GP['mm_d_mech_air_gap'].value    = SI['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_r_ro'].value              = mm_r_ro
        GP['deg_alpha_rsp'].value        = deg_alpha_rsp

    def MyCustomizedInitialDesignScript(self, fea_config_dict, SI, GP, EX):

        # inputs: air gap flux density
        guess_air_gap_flux_density_Bg = SI['guess_air_gap_flux_density_Bg']
        guess_stator_yoke_flux_density_Bys = SI['guess_stator_yoke_flux_density_Bsy']
        mm_stack_length = SI['mm_stack_length']
        Pa_TangentialStress = SI['TangentialStress']
        m  = SI['m']
        Qs = SI['Qs']
        pm = SI['pm']
        pa = SI['pa']
        mec_power = SI['mec_power']
        speed_rpm = SI['ExcitationFreqSimulated'] / pm * 60

        required_torque = mec_power/(2*np.pi*speed_rpm/60)
        guess_linear_current_density_A = Pa_TangentialStress*2/guess_air_gap_flux_density_Bg

        # (6.2)
        RotorActiveVolumn_Vr = required_torque / (2*Pa_TangentialStress)
        RotorActiveCrossSectArea_Sr = RotorActiveVolumn_Vr / (mm_stack_length*1e-3)
        RotorOuterRadius_r_or = np.sqrt(RotorActiveCrossSectArea_Sr/np.pi)
        mm_r_ro = RotorOuterRadius_r_or*1e3

        aspect_ratio__rotor_axial_to_diameter_ratio = 2*mm_r_ro/mm_stack_length

        mm_air_gap_length = 2.5 # Gruber Habil: [3,4] mm
        mm_r_si     = mm_r_ro + mm_air_gap_length
        mm_r_airgap = mm_r_ro + mm_air_gap_length*0.5

        # note this is only valid for 12s/10pp motor
        if pm == 10 or pm == 20:
            rotor_to_stator_tooth_width_ratio = 1.4 # 2010 诸自强 1点4倍的来源 Analysis of Electromagnetic Performance of Flux
        elif pm == 14 or pm == 28:
            rotor_to_stator_tooth_width_ratio = 1.0
        else:
            rotor_to_stator_tooth_width_ratio = 1.0
            # raise Exception('Not supported pm value:', pm)

        if Qs == 6:
            rotor_to_stator_tooth_width_ratio = 0.4

        if pm == 28 and Qs == 12:
            rotor_to_stator_tooth_width_ratio = 0.5
        if pm == 20 and Qs == 12:
            rotor_to_stator_tooth_width_ratio = 0.5

        mm_stator_tooth_width_w_UCoreWidth = 2*np.pi*mm_r_si / Qs / 4
        mm_d_pm = mm_stator_tooth_width_w_UCoreWidth
        if mm_d_pm > 6:
            mm_d_pm = 6

        if Qs > 18: 
            if mm_d_pm > 3:
                mm_d_pm = 3

        if Qs == 6: 
            mm_d_pm = 18

        mm_w_st = mm_stator_tooth_width_w_UCoreWidth * 3
        # print('Stator tooth', mm_stator_tooth_width_w_UCoreWidth, mm_w_st)
        mm_w_rt = mm_rotor_tooth_width_w_rt = rotor_to_stator_tooth_width_ratio * mm_stator_tooth_width_w_UCoreWidth # Empirical: slightly wider (5.18) Habil-Gruber

        mmf_per_phase = guess_linear_current_density_A * 2*np.pi*mm_r_airgap*1e-3 / m
        mmf_per_slot  = guess_linear_current_density_A * 2*np.pi*mm_r_airgap*1e-3 / Qs
        # stator_current_per_phase = mmf_per_phase /  2 * N * I

        stator_slot_area = mmf_per_slot / SI['Js'] / SI['WindingFill'] 

        # stator slot depth
        stator_inner_radius_r_is_eff = mm_r_si*1e-3
        temp = (2*np.pi*stator_inner_radius_r_is_eff - Qs*mm_stator_tooth_width_w_UCoreWidth*3*1e-3) # mm_stator_tooth_width_w_UCoreWidth * 3 = stator non-slot part width
        stator_tooth_depth_d_st = ( np.sqrt(temp**2 + 4*np.pi*stator_slot_area*Qs) - temp ) / (2*np.pi)
        mm_d_st = stator_tooth_depth_d_st * 1e3

        mm_d_sy = guess_air_gap_flux_density_Bg * mm_stator_tooth_width_w_UCoreWidth / guess_stator_yoke_flux_density_Bys
        mm_r_so = mm_r_si + mm_d_st + mm_d_sy
        split_ratio = mm_r_ro / mm_r_so

        # print('[flux_sw.py]', mm_d_st, mm_d_sy, split_ratio)

        弧长 = mm_w_rt
        deg_alpha_rsp = 弧长 / mm_r_ro /np.pi*180

        # slot_pitch_pps = np.pi * (stator_inner_diameter_Dis + stator_slot_height_h_ss) / SI['Qs']
        # kov = 1.8 # \in [1.6, 2.0]
        # EX['end_winding_length_Lew'] = end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_b_ds) + slot_pitch_pps*kov * (SI['coil_pitch_y'] - 1)

        # STATOR
        if pm == 10 or pm == 20:
            GP['deg_alpha_st'].value         = 360/Qs - 360/Qs/4/2 # deg
        elif pm == 14 or pm == 28:
            GP['deg_alpha_st'].value         = 360/Qs * 0.55 # deg # this shows lower torque density for 12s/10pp motor! But is good for 12s/14pp motor
        else:
            GP['deg_alpha_st'].value         = 360/Qs - 360/Qs/4/2 # deg
            # raise
        GP['deg_alpha_sto'].value        = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
        GP['mm_d_sto'].value             = 2 # mm
        GP['mm_d_stt'].value             = 3 * GP['mm_d_sto'].value # 1.5* # reduce tooth tip saturation
        GP['mm_r_si'].value              = mm_r_si
        GP['mm_r_so'].value              = mm_r_so
        GP['mm_d_st'].value              = mm_d_st
        GP['mm_d_sy'].value              = mm_d_sy
        # print(mm_r_so, mm_d_sy)
        # print(mm_r_so, mm_d_sy)
        # print(mm_r_so, mm_d_sy)
        GP['mm_w_st'].value              = mm_w_st
        GP['mm_d_pm'].value              = mm_d_pm
        GP['deg_alpha_pm_at_airgap'].value = mm_d_pm / (2*np.pi*mm_r_si) * 360
        # ROTOR
        GP['mm_d_sleeve'].value          = mm_air_gap_length - SI['minimum_mechanical_air_gap_length_mm']
        GP['mm_d_mech_air_gap'].value    = SI['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_r_ro'].value              = mm_r_ro
        GP['deg_alpha_rsp'].value        = deg_alpha_rsp

    def get_template_neighbor_bounds(self):
        ''' The bounds are determined around the template design.
        '''

        Q = self.SI['Qs']
        p = self.SI['p']
        # s = self.SI['no_segmented_magnets']

        GP = self.d['GP']

        # deg_alpha_rsp_ini = 360 / (self.SI['Qs']*2) * (self.SI['pm']*2)

        original_template_neighbor_bounds = {
            # STATOR
            "deg_alpha_st": [ 0.35*360/Q, 0.9*360/Q],
            "mm_d_sto":     [  0.5, 5],
            "mm_d_stt":     [  3,   8],
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.1*GP['mm_d_st'].value], # if mm_d_st is too large, the derived stator yoke can be negative
            "mm_d_sy":      [1.0*GP['mm_d_sy'].value, 1.2*GP['mm_d_sy'].value],
            "PM_to_tooth_width_ratio": [1/4, 0.3], # = 1/(1+2*k2), k2=1.5 (max), k2>=1
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],
            "split_ratio":  [0.5, 0.7], # Binder-2020-MLMS-0953@Fig.7
            "mm_d_pm":      [3,8],
            "mm_d_air_pm":  [1,10],
            "deg_alpha_pm_at_airgap": [0.8*GP['deg_alpha_pm_at_airgap'].value, 1.2*GP['deg_alpha_pm_at_airgap'].value],
            # ROTOR
            "mm_d_sleeve":  [1,   2],
            "split_ratio_rotor_salient": [0.1,0.25],
            "deg_alpha_rsp": [GP['deg_alpha_rsp'].value*0.8, GP['deg_alpha_rsp'].value*1.2]
        }
        return original_template_neighbor_bounds

    def build_design_parameters_list(self):
        return []

class FRPM_design_variant(inner_rotor_motor.variant_machine_as_objects):

    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        if x_denorm is None:
            raise

        # 初始化父类
        super(FRPM_design_variant, self).__init__(template, x_denorm, counter, counter_loop)

        # Give it a name
        self.name = f'ind{counter}'
        self.name += f'-redo{counter_loop}' if counter_loop > 1 else ''

        # Get geometric parameters and spec input
        GP = self.template.d['GP']
        SI = self.template.SI
        EX = self.template.d['EX']

        # 检查几何变量之间是否有冲突
        self.check_invalid_design(GP, SI)

        # Parts
        self.rotorCore = CrossSectInnerSalientPoleRotor.CrossSectInnerSalientPoleRotorV2(
                            mm_r_ro = GP['mm_r_ro'].value,
                            mm_d_sleeve = GP['mm_d_sleeve'].value,
                            split_ratio_rotor_salient = GP['split_ratio_rotor_salient'].value,
                            deg_alpha_rsp = GP['deg_alpha_rsp'].value,
                            pm = SI['pm'],
                            mm_r_ri=SI['mm_radius_shaft'],
                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0))

        self.shaft = CrossSectInnerNotchedRotor.CrossSectShaft(name = 'Shaft',
                                                      notched_rotor = self.rotorCore
                                                    )
        self.stator_core = CrossSectStator.CrossSectInnerRotorStator_PMAtToothTipSurface( name = 'StatorCore',
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
                                            mm_d_pm  = GP['mm_d_pm'].value,
                                            deg_alpha_pm_at_airgap  = GP['deg_alpha_pm_at_airgap'].value,
                                            mm_difference_pm_yoke = GP['mm_difference_pm_yoke'].value,
                                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                            )

        self.statorMagnet = CrossSectStator.CrossSectStatorMagnetAtToothTipSurface(name = 'StatorPM',
                                                                stator_core = self.stator_core,
                                                                mm_d_air_pm = GP['mm_d_air_pm'].value
                                                                )

        self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                               stator_core = self.stator_core)


        #03 Mechanical Parameters
        self.update_mechanical_parameters()


        self.InitialRotationAngle = EX['wily'].deg_winding_U_phase_phase_axis_angle + 360 / SI['pm'] / 4 # 假设四等分的Cell结构
        # self.InitialRotationAngle = EX['wily'].deg_winding_U_phase_phase_axis_angle - 360 / SI['pm'] / 4 # 先定位到q轴上然后再转90度电角度？
        print(f'[FRPM_design.py] self.InitialRotationAngle set to {self.InitialRotationAngle} = {EX["wily"].deg_winding_U_phase_phase_axis_angle} + {360 / SI["pm"] / 4}')
        print(f'[FRPM_design.py] self.InitialRotationAngle set to {self.InitialRotationAngle} = {EX["wily"].deg_winding_U_phase_phase_axis_angle} + {360 / SI["pm"] / 4}')
        print(f'[FRPM_design.py] self.InitialRotationAngle set to {self.InitialRotationAngle} = {EX["wily"].deg_winding_U_phase_phase_axis_angle} + {360 / SI["pm"] / 4}')
        # raise

        self.boolCustomizedCircuit = False

    def check_invalid_design(self, GP, SI):
        pass

