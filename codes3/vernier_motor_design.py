import inner_rotor_motor, pyrhonen_procedure_as_function
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
import logging

import CrossSectVShapeConsequentPoleRotor
import CrossSectStator
import Location2D

def derive_mm_w_pm(GP,SD):
    GP["mm_w_pm"].value = ( GP['mm_r_os'].value - GP["mm_d_bg_air"].value - (GP['mm_r_ri'].value + GP['mm_d_ri'].value) ) / np.cos(['deg_alpha_vspm']) - GP['mm_d_pm'].value * np.tan(['deg_alpha_vspm'])
    return GP["mm_w_pm"].value

class vernier_motor_VShapePM_template(inner_rotor_motor.template_machine_as_numbers):
    def __init__(self, fea_config_dict=None, spec_input_dict=None):
        # 基本信息
        self.machine_type = 'PMVM'
        self.name = '__PMVM'

        # 初始化搜索空间
        GP = self.d['GP']
        OP = self.d['OP']
        SD = self.SD
        childGP = OrderedDict({
            # Vernier specific
            "deg_alpha_vspm"    : acmop_parameter("free",     "v-shape_magnet_tilt_angle",     None, [None, None], lambda GP,SD:None),
            "mm_d_bg_air"       : acmop_parameter("free",     "magnet_bridge_depth_to_air",    None, [None, None], lambda GP,SD:None),
            "mm_d_bg_magnet"    : acmop_parameter("free",     "magnet_bridge_depth_to_magnet", None, [None, None], lambda GP,SD:None),
            "mm_w_pm"           : acmop_parameter("derived",  "magnet_width",                  None, [None, None], lambda GP,SD:derive_mm_w_pm(GP,SD)),
        })
        self.d.update( {"GP": childGP} )
        GP.update(childGP)

        # Get Analytical Design
        self.ModifiedBianchi2006(fea_config_dict, SD, GP, OP)

        # 定义搜索空间，determine bounds
        self.bounds_denorm = self.define_search_space(GP)

        # Template's Other Properties (Shared by the swarm)
        OP = self.get_other_properties_after_geometric_parameters_are_initialized(GP, SD)

    def ModifiedBianchi2006(self, fea_config_dict, SD, GP, OP):

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
        GP['mm_d_sleeve'].value          = sleeve_length
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
        GP["mm_w_pm"].value              = ( GP['mm_r_os'].value - GP["mm_d_bg_air"].value - (GP['mm_r_ri'].value + GP['mm_d_ri'].value) ) / np.cos(['deg_alpha_vspm']) - GP['mm_d_pm'].value * np.tan(['deg_alpha_vspm'])

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
            "deg_alpha_vspm"    : [10, 40], # TODO Should be related to pr
            "mm_d_bg_air"       : [1, 5],
            "mm_d_bg_magnet"    : [1, 5],
        }
        return original_template_neighbor_bounds

class vernier_motor_VShapePM_design_variant(inner_rotor_motor.variant_machine_as_objects):

    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        # 初始化父类
        super(vernier_motor_VShapePM_design_variant, self).__init__(template, x_denorm, counter, counter_loop)

        # 检查几何变量之间是否有冲突
        GP = self.template.d['GP']
        SD = self.template.SD
        self.check_invalid_design(GP, SD)

        quit()

        # Parts
        self.rotorCore = CrossSectVShapeConsequentPoleRotor.CrossSectVShapeConsequentPoleRotor(
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

        self.shaft = CrossSectVShapeConsequentPoleRotor.CrossSectShaft(name = 'Shaft',
                                                      notched_rotor = self.rotorCore
                                                    )

        self.rotorMagnet = CrossSectVShapeConsequentPoleRotor.CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
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

        #03 Mechanical Parameters
        self.update_mechanical_parameters()

    def check_invalid_design(self, GP, SD):
        pass

