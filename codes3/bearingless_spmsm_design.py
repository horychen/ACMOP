import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
from pprint import pprint

import CrossSectInnerNotchedRotor
import CrossSectStator
import Location2D

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
        GP = self.d['GP']
        OP = self.d['OP']
        SD = self.SD
        childGP = OrderedDict({
            # SPMSM Peculiar
            "deg_alpha_rm"      : acmop_parameter("free",      "magnet_pole_span_angle",        None, [None, None], lambda GP,SD:None),
            "mm_d_rp"           : acmop_parameter("free",      "inter_polar_iron_thickness",    None, [None, None], lambda GP,SD:None),
            "deg_alpha_rs"      : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "magnet_segment_span_angle",     None, [None, None], lambda GP,SD:None),
            "mm_d_rs"           : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "inter_segment_iron_thickness",  None, [None, None], lambda GP,SD:None),
        })
        GP.update(childGP)

       # Get Analytical Design
        self.Bianchi2006(fea_config_dict, SD, GP, OP)

        # 定义搜索空间，determine bounds
        original_template_neighbor_bounds = self.get_template_neighbor_bounds(GP, SD)        
        self.bounds_denorm = self.define_search_space(GP, original_template_neighbor_bounds)

        # Template's Other Properties (Shared by the swarm)
        OP = self.get_other_properties_after_geometric_parameters_are_initialized(GP, SD)
        # BEARING Winding Excitation Properties
        if True:
            OP['BeariW_zQ']         = OP['DriveW_zQ']
            OP['BeariW_CurrentAmp'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] * (OP['DriveW_CurrentAmp'] / fea_config_dict['TORQUE_CURRENT_RATIO'])
            OP['BeariW_Freq']       = OP['DriveW_Freq']
            OP['BeariW_Rs']         = OP['DriveW_Rs'] * OP['BeariW_zQ'] / OP['DriveW_zQ']
            OP['BeariW_poles']      = SD['ps']*2
            OP['slot_current_utilizing_ratio'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] + fea_config_dict['TORQUE_CURRENT_RATIO'] # will be less than 1 for separate winding

    def Bianchi2006(self, fea_config_dict, SD, GP, OP):

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
        GP['mm_d_sleeve'].value          = sleeve_length
        GP['mm_d_fixed_air_gap'].value   = SD['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_d_pm'].value              = 4  # mm
        GP['mm_d_ri'].value              = 1e3*ROTOR_STATOR_YOKE_HEIGHT_RATIO*stator_yoke_height_h_ys # TODO：This ratio (0.75) is randomly specified
        GP['mm_r_or'].value              = 1e3*rotor_outer_radius_r_or
        GP['mm_r_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_d_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_fixed_air_gap'].value
        # SPMSM specific
        GP['deg_alpha_rm'].value         = 0.95*360/(2*p) # deg
        GP['mm_d_rp'].value              = 3  # mm
        GP['deg_alpha_rs'].value         = 0.975*GP['deg_alpha_rm'].value / SD['no_segmented_magnets']
        GP['mm_d_rs'].value              = 0.20*GP['mm_d_rp'].value # d_pm > d_rp and d_pm > d_rs

        # Those are some obsolete variables that are convenient to have.
        # template.Radius_OuterStatorYoke = spec_geometry_dict['Radius_OuterStatorYoke'] = 1e3*0.5*stator_outer_diameter_Dse # mm
        # template.Radius_OuterRotor      = spec_geometry_dict['Radius_OuterRotor'] = 1e3*rotor_outer_radius_r_or # mm

        # Those variables are for PMSM convenience
        # template.rotor_steel_outer_radius = spec_geometry_dict['rotor_steel_outer_radius'] = template.Radius_OuterRotor

        # required_torque = SD['mec_power']/(2*np.pi*speed_rpm)*60
        # rotor_volume_Vr = required_torque/(2*SD['TangentialStress'])
        # template.required_torque = required_torque

    def get_template_neighbor_bounds(self, GP, SD):
        ''' The bounds are determined around the template design.
        '''

        Q = self.SD['Qs']
        p = self.SD['p']
        s = self.SD['no_segmented_magnets']

        GP = self.d['GP']

        original_template_neighbor_bounds = {
            # STATOR
            "deg_alpha_st": [ 0.35*360/Q, 0.9*360/Q],
            "mm_d_so":      [  0.5,   5],                                                       
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.2*GP['mm_d_st'].value],                
            "mm_r_os":      [1.0*GP['mm_r_os'].value, 1.2*GP['mm_r_os'].value], 
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],                
            # ROTOR
            "mm_d_sleeve":  [3,   6],
            "split_ratio":  [0.4, 0.6], # Binder-2020-MLMS-0953@Fig.7
            "mm_d_pm":      [2.5, 7],
            "mm_d_ri":      [0.8*GP['mm_d_ri'].value,  1.2*GP['mm_d_ri'].value],                              
            # SPMSM specific
            "deg_alpha_rm": [0.6*360/(2*p),          1.0*360/(2*p)],                                     
            "mm_d_rp":      [2.5,   6],                                                         
            "deg_alpha_rs": [0.8*360/(2*p)/s,        0.975*360/(2*p)/s],                               
            "mm_d_rs":      [2.5,   6]
        }
        return original_template_neighbor_bounds

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

class bearingless_spmsm_design_variant(inner_rotor_motor.variant_machine_as_objects):

    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        if x_denorm is None:
            raise

        # 初始化父类
        super(bearingless_spmsm_design_variant, self).__init__(template, x_denorm, counter, counter_loop)

        # 检查几何变量之间是否有冲突
        GP = self.template.d['GP']
        SD = self.template.SD
        self.check_invalid_design(GP, SD)

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
                            p = template.SD['p'], # Set pole-pairs to 2
                            s = template.SD['no_segmented_magnets'], # Set magnet segments/pole to 4
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
                                            Q = template.SD['Qs'],
                                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                            )

        self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                    stator_core = self.stator_core)

        self.sleeve = CrossSectInnerNotchedRotor.CrossSectSleeve(
                            name = 'Sleeve',
                            notched_magnet = self.rotorMagnet,
                            d_sleeve = GP['mm_d_sleeve'].value
                            )

        #03 Mechanical Parameters
        self.update_mechanical_parameters()

    def check_invalid_design(self, GP, SD):

        # 不合理的变量选择（mm_d_rp）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        if GP['mm_d_rp'].value > GP['mm_d_pm'].value:
            GP['mm_d_rp'].value            = GP['mm_d_pm'].value
            # free_variables[11] = free_variables[6]

            msg = '[Warning from bearingless_spmsm_design.py]: Inter-pole notch depth mm_d_rp cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rp is set to mm_d_pm.'
            print(msg)
            logger = logging.getLogger(__name__).warn(msg)

        # 不合理的变量选择（mm_d_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        if GP['mm_d_rs'].value > GP['mm_d_pm'].value:
            GP['mm_d_rs'].value            = GP['mm_d_pm'].value
            # free_variables[12] = free_variables[6]

            msg = '[Warning from bearingless_spmsm_design.py]: Inter-segment notch depth mm_d_rs cannot be larger than mm_d_pm or else the sleeve cannot really hold or even touch the PM. So mm_d_rs is set to mm_d_pm.'
            print(msg)
            logger = logging.getLogger(__name__).warn(msg)

        # 不合理的变量选择（deg_alpha_rs）会导致：一个变量的取值范围是受到另一个变量的取值影响的。
        if not (GP['deg_alpha_rs'].value > GP['deg_alpha_rm'].value/SD['no_segmented_magnets']):
            GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value/SD['no_segmented_magnets']
            msg = '[Warning from bearingless_spmsm_design.py]: deg_alpha_rs cannot be larger than deg_alpha_rm/s. Note deg_alpha_rs is set to deg_alpha_rm/s.'
            print(msg)
            logger = logging.getLogger(__name__).warn(msg)

        # 如果没有永磁体分段，那么alpha_rs应该等于alpha_rm。
        if SD['no_segmented_magnets'] == 1:
            GP['deg_alpha_rs'].value = GP['deg_alpha_rm'].value # raise Exception('Invalid alpha_rs. Check that it is equal to alpha_rm for s=1')
            GP['mm_d_rs'].value = 0 # raise Exception('Invalid d_rs. Check that it is equal to 0 for s =1')
        # This is all we need



























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



