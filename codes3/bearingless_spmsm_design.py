import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
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

class bearingless_spmsm_template(inner_rotor_motor.template_machine_as_numbers):
    ''' This is a surface mounted PM motor but it might have saliency on q-axis if alpha_rm is less than 180/p.
        就是说，允许永磁体陷入转子铁芯，只是永磁体外没有铁包裹防止永磁体飞出，而是需要额外增加碳纤维套。
    '''
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(bearingless_spmsm_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'SPMSM'
        self.name = '__SPMSM'

        # 初始化搜索空间
        GP = self.SI['GP'] # Geometry Parameter
        EX = self.SI['EX'] # EXcitations (was OP: Other Property)
        SI = self.SI      # Specification Input dictionary (was SD)
        childGP = OrderedDict({
            # SPMSM Peculiar
            "mm_d_pm"           : acmop_parameter("fixed",     "magnet_depth",                  None, [None, None], lambda GP,SI:None),
            "mm_d_ri"           : acmop_parameter("derived",     "rotor_iron (back iron) depth",  None, [None, None], lambda GP,SI:derive_mm_d_ri(GP,SI)),
            "deg_alpha_rm"      : acmop_parameter("fixed",     "magnet_pole_span_angle",        None, [None, None], lambda GP,SI:None),
            "mm_d_rp"           : acmop_parameter("fixed",     "inter_polar_iron_thickness",    None, [None, None], lambda GP,SI:None),
            "deg_alpha_rs"      : acmop_parameter("fixed" if SI['no_segmented_magnets']!=1 else "fixed",   "magnet_segment_span_angle",     None, [None, None], lambda GP,SI:None),
            "mm_d_rs"           : acmop_parameter("fixed" if SI['no_segmented_magnets']!=1 else "fixed",   "inter_segment_iron_thickness",  None, [None, None], lambda GP,SI:None),
            "mm_r_ri"           : acmop_parameter("derived",  "rotor_inner_radius",            None, [None, None], lambda GP,SI:derive_mm_r_ri(GP,SI)),
        })
        GP.update(childGP)

        self.set_gp_values_and_types_based_on_fea_config(spec_input_dict)
        GP = self.update_geometric_parameters_using_x_denorm_dict(self.SI['x_denorm_dict'])

        # Get Analytical Design
        self.PracticalInitialDesign(fea_config_dict, SI, GP, EX)

        print(self.SI['x_denorm_dict'])

        # 定义搜索空间，determine bounds
        self.original_template_neighbor_bounds = self.get_template_neighbor_bounds()
        self.bounds_denorm = self.SIefine_search_space(GP, self.original_template_neighbor_bounds)

        # Template's Other Properties (Shared by the swarm)
        EX = self.get_other_properties_after_geometric_parameters_are_initialized(GP, SI)
        # BEARING Winding Excitation Properties
        if True:
            EX['BeariW_zQ']         = EX['DriveW_zQ']
            EX['BeariW_CurrentAmp'] = fea_config_dict['circuit.SUSPENSION_CURRENT_RATIO'] * (EX['DriveW_CurrentAmp'] / fea_config_dict['circuit.TORQUE_CURRENT_RATIO'])
            EX['BeariW_Freq']       = EX['DriveW_Freq']
            EX['BeariW_Rs']         = EX['DriveW_Rs'] * EX['BeariW_zQ'] / EX['DriveW_zQ']
            EX['BeariW_poles']      = SI['ps']*2
            EX['slot_current_utilizing_ratio'] = fea_config_dict['circuit.SUSPENSION_CURRENT_RATIO'] + fea_config_dict['circuit.TORQUE_CURRENT_RATIO'] # will be less than 1 for separate winding

    # def Bianchi2006(self, fea_config_dict, SI, GP, EX):

    #     # inputs: air gap flux density
    #     air_gap_flux_density_Bg = SI['guess_air_gap_flux_density_Bg'] # 0.9 T
    #     stator_tooth_flux_density_Bst = SI['guess_stator_tooth_flux_density_Bst'] # 1.5 T
    #     stator_yoke_flux_density_Bsy = SI['guess_stator_yoke_flux_density_Bsy']

    #     if SI['p'] >= 2:
    #         ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.75
    #         alpha_rm_over_alpha_rp = 1.0
    #         # stator_yoke_flux_density_Bsy = 1.2
    #     else:
    #         # penalty for p=1 motor, i.e., large yoke height
    #         ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.5
    #         alpha_rm_over_alpha_rp = 0.75
    #         # stator_yoke_flux_density_Bsy = 1.5

    #     # ureg = pint.UnitRegistry()  # 0.225* ureg.meter
    #     stator_outer_diameter_Dse = SI['mm_stator_outer_diameter'] * 1e-3 # this is related to the stator current density and should be determined by Js and power.
    #     sleeve_length = SI['mm_sleeve_length'] * 1e-3 # mm

    #     speed_rpm = SI['ExcitationFreqSimulated'] * 60 / SI['p'] # rpm
    #     rotor_outer_radius_r_or = pyrhonen_procedure_as_function.eric_specify_tip_speed_get_radius(SI['tip_speed'], speed_rpm)

    #     rotor_outer_diameter_Dr = rotor_outer_radius_r_or*2
    #     stator_inner_radius_r_is  = rotor_outer_radius_r_or + (sleeve_length+SI['minimum_mechanical_air_gap_length_mm'])*1e-3 # m (sleeve 3 mm, air gap 0.75 mm)
    #     stator_inner_diameter_Dis = stator_inner_radius_r_is*2
    #     split_ratio = stator_inner_diameter_Dis / stator_outer_diameter_Dse

    #     stator_yoke_height_h_ys = air_gap_flux_density_Bg * np.pi * stator_inner_diameter_Dis * alpha_rm_over_alpha_rp / (2*stator_yoke_flux_density_Bsy * 2*SI['p'])
    #     # print(stator_outer_diameter_Dse, stator_inner_diameter_Dis, stator_yoke_height_h_ys)
    #     stator_tooth_height_h_ds = (stator_outer_diameter_Dse - stator_inner_diameter_Dis) / 2 - stator_yoke_height_h_ys
    #     stator_slot_height_h_ss = stator_tooth_height_h_ds
    #     stator_tooth_width_b_ds = air_gap_flux_density_Bg * np.pi * stator_inner_diameter_Dis / (stator_tooth_flux_density_Bst* SI['Qs'])

    #     EX['stator_slot_area'] = stator_slot_area = np.pi/(4*SI['Qs']) * ((stator_outer_diameter_Dse - 2*stator_yoke_height_h_ys)**2 - stator_inner_diameter_Dis**2) - stator_tooth_width_b_ds * stator_tooth_height_h_ds

    #     slot_pitch_pps = np.pi * (stator_inner_diameter_Dis + stator_slot_height_h_ss) / SI['Qs']
    #     kov = 1.8 # \in [1.6, 2.0]
    #     EX['end_winding_length_Lew'] = end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_b_ds) + slot_pitch_pps*kov * (SI['coil_pitch_y'] - 1)

    #     Q = SI['Qs']
    #     p = SI['p']
    #     # STATOR
    #     GP['deg_alpha_st'].value         = 360/Q - 2 # deg
    #     GP['deg_alpha_sto'].value         = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
    #     GP['mm_r_si'].value              = 1e3*stator_inner_radius_r_is # mm
    #     GP['mm_r_so'].value              = 1e3*stator_outer_diameter_Dse/2 # mm
    #     GP['mm_d_sto'].value              = 5 # mm
    #     GP['mm_d_stt'].value              = 1.5*GP['mm_d_sto'].value
    #     GP['mm_d_st'].value              = 1e3*(0.5*stator_outer_diameter_Dse - stator_yoke_height_h_ys) - GP['mm_r_si'].value - GP['mm_d_stt'].value  # mm
    #     # print(GP['mm_d_st'].value)
    #     # print (1e3*stator_outer_diameter_Dse)
    #     # print(1e3*stator_yoke_height_h_ys)
    #     # print(GP['mm_r_si'].value)
    #     # print (GP['mm_d_stt'].value)
    #     # quit()
    #     GP['mm_d_sy'].value              = 1e3*stator_yoke_height_h_ys # mm
    #     GP['mm_w_st'].value              = 1e3*stator_tooth_width_b_ds # mm
    #     # ROTOR
    #     GP['mm_d_sleeve'].value          = sleeve_length
    #     GP['mm_d_mech_air_gap'].value    = SI['minimum_mechanical_air_gap_length_mm']
    #     GP['split_ratio'].value          = split_ratio
    #     GP['mm_d_pm'].value              = 4  # mm
    #     # GP['mm_d_ri'].value              = 1e3*ROTOR_STATOR_YOKE_HEIGHT_RATIO*stator_yoke_height_h_ys # TODO：This ratio (0.75) is epirically specified
    #     GP['mm_r_ro'].value              = 1e3*rotor_outer_radius_r_or
    #     # GP['mm_r_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_d_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value
    #     GP['mm_r_ri'].value              = SI['mm_radius_shaft']
    #     GP['mm_d_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_r_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value
    #     # SPMSM specific
    #     GP['deg_alpha_rm'].value         = 1.0*360/(2*p) # deg
    #     GP['mm_d_rp'].value              = 3  # mm
    #     GP['deg_alpha_rs'].value         = GP['deg_alpha_rm'].value / SI['no_segmented_magnets']
    #     GP['mm_d_rs'].value              = 0.20*GP['mm_d_rp'].value # d_pm > d_rp and d_pm > d_rs

    def PracticalInitialDesign(self, fea_config_dict, SI, GP, EX):

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

        stator_outer_diameter_Dse = GP['mm_r_so'].value*2 * 1e-3 # this is related to the stator current density and should be determined by Js and power.
        sleeve_length = GP['mm_d_sleeve'].value * 1e-3 # mm

        rotor_outer_radius_r_or = 1e-3* (GP['mm_r_ri'].value + GP['mm_d_pm'].value + GP['mm_d_ri'].value)
        stator_inner_radius_r_is  = rotor_outer_radius_r_or + (sleeve_length+SI['minimum_mechanical_air_gap_length_mm'])*1e-3 # [m]
        stator_inner_diameter_Dis = stator_inner_radius_r_is*2
        split_ratio = stator_inner_diameter_Dis / stator_outer_diameter_Dse

        stator_yoke_height_h_ys = air_gap_flux_density_Bg * np.pi * stator_inner_diameter_Dis * alpha_rm_over_alpha_rp / (2*stator_yoke_flux_density_Bsy * 2*SI['p'])
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
        GP['mm_d_sto'].value              = stator_yoke_height_h_ys*1e3*0.1 # mm
        GP['mm_d_stt'].value              = 1.5*GP['mm_d_sto'].value
        GP['mm_d_st'].value              = 1e3*(0.5*stator_outer_diameter_Dse - stator_yoke_height_h_ys) - GP['mm_r_si'].value - GP['mm_d_stt'].value  # mm
        if GP['mm_d_st'].value<=0:
            print('背铁太厚了')
            print(GP['mm_d_st'].value)
            print (1e3*stator_outer_diameter_Dse)
            print(1e3*stator_yoke_height_h_ys)
            print(GP['mm_r_si'].value)
            print (GP['mm_d_stt'].value)
            quit()
        GP['mm_d_sy'].value              = 1e3*stator_yoke_height_h_ys # mm
        GP['mm_w_st'].value              = 1e3*stator_tooth_width_b_ds # mm
        # ROTOR
        GP['mm_d_sleeve'].value          = sleeve_length*1e3
        GP['mm_d_mech_air_gap'].value    = SI['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_d_pm'].value              = SI['mm_d_pm']  # mm
        # GP['mm_d_ri'].value              = 1e3*ROTOR_STATOR_YOKE_HEIGHT_RATIO*stator_yoke_height_h_ys # TODO：This ratio (0.75) is epirically specified
        GP['mm_r_ro'].value              = 1e3*rotor_outer_radius_r_or
        # GP['mm_r_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_d_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value
        GP['mm_r_ri'].value              = SI['mm_radius_shaft']
        GP['mm_d_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_r_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value
        # SPMSM specific
        GP['deg_alpha_rm'].value         = 1.0*360/(2*p) # deg
        GP['mm_d_rp'].value              = SI['mm_d_pm']  # mm
        GP['deg_alpha_rs'].value         = GP['deg_alpha_rm'].value / SI['no_segmented_magnets']
        GP['mm_d_rs'].value              = 0.20*GP['mm_d_rp'].value # d_pm > d_rp and d_pm > d_rs

    def get_template_neighbor_bounds(self):
        ''' The bounds are determined around the template design.
        '''

        Q = self.SI['Qs']
        p = self.SI['p']
        s = self.SI['no_segmented_magnets']

        GP = self.SI['GP']

        ######################    get bounds have a misalignment    ######################
        # attention: the bounds are determined around the template design, which means any change of the template design will lead to a change of the order the bounds.
        original_template_neighbor_bounds = {
            # Sleeve
            "mm_d_sleeve":  [3,   6], 
            # ROTOR
            # "split_ratio":  [0.4, 0.6], # Binder-2020-MLMS-0953@Fig.7
            "split_ratio":  [0.35, 0.5], # Q12p4优化的时候，轭部经常不够用，所以就把split_ratio减小——Exception: ('Error: Negative derived parameter', "acmop_parameter(type='derived', name='stator_yoke_depth', value=-1.362043443071423, bounds=[None, None], calc=<function template_machine_as_numbers.__init__.<locals>.<lambda> at 0x00000237CC403D30>)")
            # STATOR
            "deg_alpha_st": [ 0.35*360/Q, 1.0*360/Q],
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.1*GP['mm_d_st'].value], # if mm_d_st is too large, the derived stator yoke can be negative
            "mm_d_sto":     [0.2,                                         0.8], # this will influence split_ratio
            "mm_r_so":      [1.0*GP['mm_r_so'].value, 1.2*GP['mm_r_so'].value],
            "mm_d_pm":      [1, 3],
            "mm_d_ri":      [0.8*GP['mm_d_ri'].value,  1.2*GP['mm_d_ri'].value],
            # SPMSM specific
            "deg_alpha_rm": [0.6*360/(2*p),          1.0*360/(2*p)],
            "mm_d_rp":      [GP['mm_d_pm'].value/2,   GP['mm_d_pm'].value],
            # Rest parameters haven't been determined yet (need to be determined by the user defined fixed or free variables)
            "deg_alpha_rs": [0.8*360/(2*p)/s,        GP['deg_alpha_rm'].value/s],
            "mm_d_rs":      [2.5,   6],
            "mm_d_sy":      [1.0*GP['mm_d_sy'].value, 1.2*GP['mm_d_sy'].value]
        }
        # print('原始约束空间为：')
        # for k,v in original_template_neighbor_bounds.items(): print('\t', k,v)
        return original_template_neighbor_bounds

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



