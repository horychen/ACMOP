import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
from pprint import pprint

import CrossSectInnerNotchedRotor
import CrossSectStator
import Location2D

class bearingless_induction_template(inner_rotor_motor.template_machine_as_numbers):
    ''' This is bearingless induction motor template
    '''
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(bearingless_induction_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'BLIM'
        self.name = '__BLIM'

        # 初始化搜索空间
        GP = self.d['GP']
        OP = self.d['OP']
        SD = self.SD
        childGP = OrderedDict({
            # IM Peculiar
            "mm_d_ro"       : acmop_parameter("free",   "rotor_slot_open_depth",    None, [None, None], lambda GP,SD:None),
            "mm_w_rt"       : acmop_parameter("free",   "rotor_tooth_width",        None, [None, None], lambda GP,SD:None),
            "deg_alpha_rt"  : acmop_parameter("free",   "rotor_tooth_span_angle",   None, [None, None], lambda GP,SD:None),
        })
        GP.update(childGP)

        if True:
            # Get Analytical Design

            # [2.1] Attain spec_derive_dict
            self.spec = spec = pyrhonen_procedure_as_function.desgin_specification(**spec_input_dict)
            print(spec.build_name())
            spec.bool_bad_specifications = spec.pyrhonen_procedure()

            # Update GP
            if True:
                sgd = spec.spec_geometry_dict
                stator_outer_diameter_Dse = sgd['Radius_OuterStatorYoke'] * 2
                stator_inner_diameter_Dis = sgd['Radius_InnerStatorYoke'] * 2
                stator_yoke_height_h_ys   = sgd['Radius_OuterStatorYoke'] - sgd['Radius_InnerStatorYoke']
                stator_tooth_width_b_ds   = sgd['Width_StatorTeethBody']
                stator_tooth_height_h_ds  = sgd['Width_StatorTeethNeck']
                stator_slot_height_h_ss   = stator_tooth_height_h_ds
                stator_inner_radius_r_is  = 0.5*stator_inner_diameter_Dis
                rotor_outer_radius_r_or   = sgd['Radius_OuterRotor']
                rotor_outer_diameter_Dr   = 2.0*rotor_outer_radius_r_or
                split_ratio               = stator_inner_diameter_Dis / stator_outer_diameter_Dse

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
                GP['mm_d_sleeve'].value          = 0.0
                GP['mm_d_fixed_air_gap'].value   = sgd['Length_AirGap']
                GP['split_ratio'].value          = split_ratio
                GP['mm_r_or'].value              = 1e3*rotor_outer_radius_r_or
                GP['mm_r_ri'].value              = sgd['Radius_Shaft']
                # BLIM specific
                GP['mm_d_ro'].value              = sgd['Length_HeadNeckRotorSlot']
                GP['mm_w_rt'].value              = sgd['rotor_tooth_width_b_dr']
                GP['deg_alpha_rt'].value         = 360/Q - sgd['Width_RotorSlotOpen']/(2*np.pi*rotor_outer_radius_r_or)*360

                GP = None
                print('!!![bearingless_induction_motor_design.py] GP is set to None for now.')
                print('!!![bearingless_induction_motor_design.py] DEBUG pending: GP has different order from x_denorm parsed in population.bearingless_induction_motor_design.local_design_variant().')

            import population
            # load initial design using the obsolete class bearingless_induction_motor_design
            self.obsolete_template = population.bearingless_induction_motor_design(spec_input_dict, spec.spec_derive_dict, spec.spec_geometry_dict, fea_config_dict)
            # self.spec_geometry_dict = spec.spec_geometry_dict


        # 定义搜索空间，determine bounds
        self.bounds_denorm = spec.get_im_classic_bounds(which_filter=fea_config_dict['which_filter'])
        bound_filter       = spec.bound_filter
        otnb               = spec.original_template_neighbor_bounds

        # Template's Other Properties (Shared by the swarm)
        self.obsolete_template

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

    def build_x_denorm(self):
        ''' 覆盖父类的同名方法
        '''
        return self.obsolete_template.build_x_denorm()

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

        im_variant = population.bearingless_induction_motor_design.local_design_variant(template, 0, counter, x_denorm)


        # # Parts
        # self.rotorCore = CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
        #                     name = 'NotchedRotor',
        #                     mm_d_pm      = GP['mm_d_pm'].value,
        #                     deg_alpha_rm = GP['deg_alpha_rm'].value, # angular span of the pole: class type DimAngular
        #                     deg_alpha_rs = GP['deg_alpha_rs'].value, # segment span: class type DimAngular
        #                     mm_d_ri      = GP['mm_d_ri'].value, # rotor iron thickness: class type DimLinear
        #                     mm_r_ri      = GP['mm_r_ri'].value, # inner radius of rotor: class type DimLinear
        #                     mm_d_rp      = GP['mm_d_rp'].value, # interpolar iron thickness: class type DimLinear
        #                     mm_d_rs      = GP['mm_d_rs'].value, # inter segment iron thickness: class type DimLinear
        #                     p = template.SD['p'], # Set pole-pairs to 2
        #                     s = template.SD['no_segmented_magnets'], # Set magnet segments/pole to 4
        #                     location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0))

        # self.shaft = CrossSectInnerNotchedRotor.CrossSectShaft(name = 'Shaft',
        #                                               notched_rotor = self.rotorCore
        #                                             )

        # self.rotorMagnet = CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
        #                                               notched_rotor = self.rotorCore
        #                                             )

        # self.stator_core = CrossSectStator.CrossSectInnerRotorStator( name = 'StatorCore',
        #                                     deg_alpha_st = GP['deg_alpha_st'].value, #40,
        #                                     deg_alpha_so = GP['deg_alpha_so'].value, #20,
        #                                     mm_r_si = GP['mm_r_si'].value,
        #                                     mm_d_so = GP['mm_d_so'].value,
        #                                     mm_d_sp = GP['mm_d_sp'].value,
        #                                     mm_d_st = GP['mm_d_st'].value,
        #                                     mm_d_sy = GP['mm_d_sy'].value,
        #                                     mm_w_st = GP['mm_w_st'].value,
        #                                     mm_r_st = 0.0, # =0
        #                                     mm_r_sf = 0.0, # =0
        #                                     mm_r_sb = 0.0, # =0
        #                                     Q = template.SD['Qs'],
        #                                     location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
        #                                     )

        # self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
        #                                             stator_core = self.stator_core)

        # self.sleeve = CrossSectInnerNotchedRotor.CrossSectSleeve(
        #                     name = 'Sleeve',
        #                     notched_magnet = self.rotorMagnet,
        #                     d_sleeve = GP['mm_d_sleeve'].value
        #                     )

        # #03 Mechanical Parameters
        # self.update_mechanical_parameters()

    def check_invalid_design(self, GP, SD):
        pass

