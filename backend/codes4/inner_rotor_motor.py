from dataclasses import dataclass
import pyrhonen_procedure_as_function, winding_layout
import numpy as np
import logging
import math
import utility
from utility import acmop_parameter
from time import time as clock_time
from collections import OrderedDict, namedtuple

# Abbreviations:
# GP = Geometric Parameters
# EX = Other Properties
# SI = Specification Dictionary

def derive_mm_r_si(GP,SI):
    # (option 1) depends on split_ratio and r_os
    GP       ['mm_r_si'].value = GP['mm_r_so'].value * GP['split_ratio'].value
    return GP['mm_r_si'].value

    # (option 2) depends on d_sy (which is bad, as d_sy is also derived) and r_os
    # GP['mm_r_si'].value = GP['mm_r_so'].value - GP['mm_d_sy'].value - GP['mm_d_st'].value - GP['mm_d_sts'].value
    # return GP['mm_r_si'].value

    # (option 3) depends on r_or and air gap length
    # GP       ['mm_r_si'].value = GP['mm_r_ro'].value + GP['mm_d_sleeve'].value + GP['mm_d_mech_air_gap'].value
    # return GP['mm_r_si'].value

def derive_mm_d_sy(GP,SI):
    GP       ['mm_d_sy'].value = GP['mm_r_so'].value - GP['mm_r_si'].value - GP['mm_d_st'].value - GP['mm_d_sts'].value
    return GP['mm_d_sy'].value

def derive_mm_r_ro(GP,SI):
    GP       ['mm_r_ro'].value = GP['mm_r_si'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value
    return GP['mm_r_ro'].value

def derive_split_ratio(GP,SI):
    GP       ['split_ratio'].value          = GP['mm_r_si'].value / GP['mm_r_so'].value
    return GP['split_ratio'].value

def derive_mm_d_st(GP,SI):
    if 'mm_d_sts' not in GP.keys():
        GP['mm_d_st'].value = GP['mm_r_so'].value - GP['mm_r_si'].value - GP['mm_d_sy'].value
    else:
        GP['mm_d_st'].value = GP['mm_r_so'].value - GP['mm_r_si'].value - GP['mm_d_sy'].value - GP['mm_d_sts'].value
    return GP['mm_d_st'].value

class template_machine_as_numbers(object):
    ''' # template 有点类似analytical的电机（由几何尺寸组成）
    '''
    def __init__(self, fea_config_dict=None, spec_input_dict=None):

        # 仿真输入
        self.FEA = self.fea_config_dict = fea_config_dict
        self.SI = self.spec_input_dict = spec_input_dict

        # 初始化搜索空间
        geometric_parameters = OrderedDict({
            # ROTOR                                  Type       Name                          Value  Bounds       Calc
            "mm_r_ro"           : acmop_parameter("fixed",    "rotor_outer_radius",            None, [None, None], lambda GP,SI:None), #derive_mm_r_ro(GP,SI)),
            "mm_d_mech_air_gap" : acmop_parameter("fixed",    "mechanical_air_gap_length",     None, [None, None], lambda GP,SI:None),
            "mm_d_sleeve"       : acmop_parameter("fixed",    "sleeve_length",                 None, [None, None], lambda GP,SI:None),
            "split_ratio"       : acmop_parameter("fixed",    "split_ratio_r_is_slash_r_os",   None, [None, None], lambda GP,SI:None), #derive_split_ratio(GP,SI)),
            # STATOR                           Type       Name                          Value  Bounds       Calc
            "mm_w_st"       : acmop_parameter("fixed",    "stator_tooth_width"         , None, [None, None], lambda GP,SI:None),
            "mm_r_si"       : acmop_parameter("derived", "stator_inner_radius"        , None, [None, None], lambda GP,SI:derive_mm_r_si(GP,SI)),
            "mm_r_so"       : acmop_parameter("fixed", "stator_outer_radius"        , None, [None, None], lambda GP,SI:None),
            "mm_d_sy"       : acmop_parameter("fixed",   "stator_yoke_depth"          , None, [None, None], lambda GP,SI:derive_mm_d_sy(GP,SI)),
            "mm_d_st"       : acmop_parameter("fixed",    "stator_tooth_depth"         , None, [None, None], lambda GP,SI:derive_mm_d_st(GP,SI)),
            "mm_d_sts"      : acmop_parameter("fixed",    "stator_tooth_shoe_depth"        , None, [None, None], lambda GP,SI:None),
        })

        # assert derived variables
        for k, v in geometric_parameters.items():
            if v.type == 'derived' and v.calc is None:
                raise Exception('calc method is not defined for the derived acmop_parameter:', v)

        # all in one place
        self.SI["GP"] = geometric_parameters
        self.SI["x_denorm_dict"] = OrderedDict()
        self.SI["bounds_denorm_dict"] = OrderedDict()

        # self.set_gp_values_and_types_based_on_SI(spec_input_dict)

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
        stator_yoke_depth_d_sy  = Bg * math.pi * stator_inner_diameter_Dsi * get_alpha_rm_over_alpha_rp(p) / (2*Bsy * 2*p)

        # 定子齿部深度依赖于轭部深度
        stator_tooth_depth_d_st = (stator_outer_diameter_Dso - stator_inner_diameter_Dsi) *0.5 - stator_yoke_depth_d_sy
        stator_slot_depth_d_ss = stator_tooth_depth_d_st

        # 单个永磁体极面下，气隙中的磁通全部进入Qs个定子齿部（的宽度）
        stator_tooth_width_w_st = Bg * math.pi * stator_inner_diameter_Dsi / (Bst* Q)

        # 计算槽面积
        def get_stator_slot_area(Q, stator_outer_diameter_Dso, stator_yoke_depth_d_sy, stator_inner_diameter_Dsi, stator_tooth_width_w_st, stator_tooth_depth_d_st):
            return math.pi/(4*Q) * ((stator_outer_diameter_Dso - 2*stator_yoke_depth_d_sy)**2 - stator_inner_diameter_Dsi**2) - stator_tooth_width_w_st * stator_tooth_depth_d_st
        EX['stator_slot_area'] = get_stator_slot_area(Q, stator_outer_diameter_Dso, stator_yoke_depth_d_sy, stator_inner_diameter_Dsi, stator_tooth_width_w_st, stator_tooth_depth_d_st)

        # 计算端部绕组长度
        slot_pitch_pps = math.pi * (stator_inner_diameter_Dsi + stator_slot_depth_d_ss) / Q
        kov = EX['end_winding_length_factor_kov']
        EX['end_winding_length_Lew'] = math.pi*0.5 * (slot_pitch_pps + stator_tooth_width_w_st) + slot_pitch_pps*kov * (SI['coil_pitch_y'] - 1)

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
        GP['mm_d_ri'].value              = 1e3*(stator_inner_radius_r_is - equivalent_air_gap_length) - GP['mm_d_pm'].value - GP['mm_r_ri'].value


    ''' 实用
    '''
    def get_rotor_volume(self, stack_length=None):
        if stack_length is None:
            return math.pi*(self.SI['GP']['mm_r_ro'].value*1e-3)**2 * (self.SI['EX']['mm_stack_length']*1e-3)
        else:
            return math.pi*(self.SI['GP']['mm_r_ro'].value*1e-3)**2 * (stack_length*1e-3)
    def get_rotor_weight(self, gravity=9.8, stack_length=None):
        material_density_rho = pyrhonen_procedure_as_function.get_material_data()[0]
        if stack_length is None:
            return gravity * self.get_rotor_volume() * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg
        else:
            return gravity * self.get_rotor_volume(stack_length=stack_length) * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg

    # ''' 玩弄几何变量
    # '''
    # def build_x_denorm(self):
    #     """ This is core function """
    #     # this is used in part_evaluation
    #     GP = self.SI['GP']
    #     self.x_denorm_dict = self.get_x_denorm_dict_from_geometric_parameters(GP)
    #     x_denorm = [val for key, val in self.x_denorm_dict.items()]
    #     if False:
    #         from pprint import pprint
    #         pprint(x_denorm_dict)
    #         pprint(x_denorm)
    #         pprint(GP)
    #         quit()
    #     return x_denorm
    # def get_x_denorm_dict_from_geometric_parameters(self, GP):
    #     x_denorm_dict = OrderedDict()
    #     for key, parameter in GP.items():
    #         if parameter.type == 'free':
    #             x_denorm_dict[key] = parameter.value
    #     return x_denorm_dict
    # def get_x_denorm_dict_from_x_denorm_list(self, x_denorm):
    #     # 先拿个模板来，但是几何尺寸的变量值是旧的
    #     x_denorm_dict = self.get_x_denorm_dict_from_geometric_parameters(self.SI['GP'])

    #     # print('[inner_rotor_motor.py] DEBUG x_denorm_dict:', x_denorm_dict.keys())

    #     # 对模板进行遍历，挨个把新的几何尺寸的值从x_denorm中读取出来并更新x_denorm_dict
    #     for key, new_val in zip(x_denorm_dict.keys(), x_denorm):
    #         x_denorm_dict[key] = new_val
    #     return x_denorm_dict
    # def update_geometric_parameters_using_x_denorm_dict(self, x_denorm_dict):
    #     # Update Free Parameters (a.k.a. x_denorm)
    #     for key, val in x_denorm_dict.items():
    #         self.SI['GP'][key].value = val

    #     # Update Derived Parameters
    #     for key, parameter in self.SI['GP'].items():
    #         if parameter.type == 'derived':
    #             # print(parameter)
    #             parameter.value = None # 先全部清空，防止编写derive时搞错依赖项关系，更新顺序是有先后的，后面的derived parameter 可以利用前面的derived parameter的值。
    #     # print()
    #     count_TypeError = 0
    #     list_parameter_to_derive = []
    #     for key, parameter in self.SI['GP'].items():
    #         if parameter.type == 'derived':
    #             try:
    #                 parameter.value = parameter.calc(self.SI['GP'], self.SI)
    #                 list_parameter_to_derive.append(key)
    #             except TypeError as e: # TypeError: unsupported operand type(s) for -: 'NoneType' and 'NoneType' 用来计算的变量还未被赋值
    #                 count_TypeError += 1
    #                 logger = logging.getLogger(__name__)
    #                 logger.warning('%s TypeError: None is used for derivation of %s', f'{count_TypeError=}', f'{parameter=}')
    #                 print('%s TypeError: None is used for derivation of %s', f'{count_TypeError=}', f'{parameter=}')
    #             else: # no exception
    #                 if parameter.value<=0:
    #                     logger = logging.getLogger(__name__)
    #                     logger.error('Negative geometric parameter:')
    #                     for k,v in self.SI['GP'].items():
    #                         logger.error('\t %s: %s', k, v)
    #                     raise Exception('Error: Negative derived parameter', str(parameter))
    #     for key in list_parameter_to_derive:
    #         parameter =  self.SI['GP'][key]
    #         try:
    #             parameter.value = parameter.calc(self.SI['GP'], self.SI)
    #         except TypeError as e:
    #             logger = logging.getLogger(__name__)
    #             logger.warning('TypeError fixed: %s', f'{parameter=}')
    #             pass
    #         count_TypeError -= 1
    #         logger = logging.getLogger(__name__)
    #         logger.info('TypeError fixed: %s', f'{count_TypeError=}')
    #     # 【太蠢啦】针对“用来计算的变量还未被赋值”的变量，再次调用它的calc方法。
    #     # while count_TypeError>0:
    #     #     for key, parameter in self.SI['GP'].items():
    #     #         if parameter.type == 'derived' and parameter.value is None:
    #     #             try:
    #     #                 parameter.value = parameter.calc(self.SI['GP'], self.SI)
    #     #             except TypeError as e: # TypeError: unsupported operand type(s) for -: 'NoneType' and 'NoneType' 用来计算的变量还未被赋值
    #     #                 count_TypeError += 1
    #     #                 print('[inner_rotor_motor.py] [Re] TypeError: None is used for derivation.')
    #     #                 pass
    #     #             finally:
    #     #                 count_TypeError -= 1
    #     #                 print(f'[inner_rotor_motor.py] {count_TypeError=} has derived: {parameter=}')
    #     return self.SI['GP']

class variant_machine_as_objects(object):
    ''' # variant则有点像是具体的电机实现类（由各个局部类，比如转子、定子等组成）
    '''
    def __init__(self, spmsm_template=None, x_denorm=None, counter=None, counter_loop=None, 
                verbose=True):

        self.template = spmsm_template

        #00 Settings
        # self.template.fea_config_dict = spmsm_template.fea_config_dict
        # self.template.spec_input_dict = spmsm_template.spec_input_dict
        # self.spec_geometry_dict = spmsm_template.spec_geometry_dict

        SI = self.template.spec_input_dict

        #01 Model ID
        # self.model_name_prefix
        self.counter = counter
        self.counter_loop = counter_loop
        if counter is not None:
            if counter_loop == 1:
                self.name = f"p{SI['p']}ps{SI['ps']}-Q{SI['Qs']}y{SI['coil_pitch_y']}-{counter}"
            else:
                self.name = f"p{SI['p']}ps{SI['ps']}-Q{SI['Qs']}y{SI['coil_pitch_y']}-{counter}-redo{counter_loop}"
        else:
            self.name = 'SPMSM_InitialDesign'
        self.ID = 'ID'

        #02 Geometry Data
        if x_denorm is None:
            # template as variant
            GP = self.template.d['GP'] # do nothing, use template's GP

    def reproduce_wily(self):
        ''' This method is only used for reproducing design from jsonpickle'''
        self.template.d['EX']['wily'] = winding_layout.winding_layout_v2(self.template.SI['DPNV_or_SEPA'], self.template.SI['Qs'], self.template.SI['p'], self.template.SI['ps'], self.template.SI['coil_pitch_y'])

    def update_mechanical_parameters(self, syn_freq=None):
        EX = self.template.SI['EX']
        SI = self.template.SI
        if syn_freq is None:
            if 'number_of_rotor_pole_pairs' not in SI.keys():
                number_of_rotor_pole_pairs = SI['p']
            else:
                number_of_rotor_pole_pairs = SI['number_of_rotor_pole_pairs']
            EX['the_speed'] = EX['DriveW_Freq']*60. / number_of_rotor_pole_pairs # rpm
            EX['Omega']     = EX['the_speed'] / 60. * 2*math.pi
            # self.omega = None # This variable name is devil! you can't tell its electrical or mechanical! #+ self.SIriveW_Freq * (1-self.the_slip) * 2*pi
        else:
            raise Exception('Not implemented.')

    def get_individual_name(self):
        if self.template.fea_config_dict['flag_optimization'] == True:
            return "ID%s" % (self.ID)
        else:
            return "%s_ID%s" % (self.name, self.ID)

