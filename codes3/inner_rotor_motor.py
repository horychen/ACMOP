from dataclasses import dataclass
import pyrhonen_procedure_as_function, winding_layout
import numpy as np
import logging
import utility
from utility import acmop_parameter, EPS
from time import time as clock_time
from collections import OrderedDict, namedtuple

# Abbreviations:
# GP = Geometric Parameters
# EX = Other Properties
# SI = Specification Dictionary

def derive_mm_r_si(GP,SI):
    # (option 1) depends on split_ratio and r_os
    GP       ['mm_r_si'].value = GP['mm_r_os'].value * GP['split_ratio'].value
    return GP['mm_r_si'].value
        # (option 2) depends on d_sy (which is bad, as d_sy is also derived) and r_os
        # GP['mm_r_si'].value = GP['mm_r_os'].value - GP['mm_d_sy'].value - GP['mm_d_st'].value - GP['mm_d_sp'].value
        # return GP['mm_r_si'].value

def derive_deg_alpha_so(GP,SI):
    # this allows a square tooth tip, or else the tooth tip might be pointy
    GP       ['deg_alpha_so'].value = GP['deg_alpha_st'].value/2
    return GP['deg_alpha_so'].value

def derive_mm_d_sp(GP,SI):
    GP       ['mm_d_sp'].value = 1.5*GP['mm_d_so'].value
    return GP['mm_d_sp'].value

def derive_mm_d_sy(GP,SI):
    GP       ['mm_d_sy'].value = GP['mm_r_os'].value - GP['mm_r_si'].value - GP['mm_d_st'].value - GP['mm_d_sp'].value
    return GP['mm_d_sy'].value

def derive_mm_r_or(GP,SI):
    GP       ['mm_r_or'].value = GP['mm_r_si'].value - GP['mm_d_sleeve'].value - GP['mm_d_fixed_air_gap'].value
    return GP['mm_r_or'].value

def derive_mm_r_ri(GP,SI):
    GP       ['mm_r_ri'].value = GP['mm_r_or'].value - GP['mm_d_pm'].value - GP['mm_d_ri'].value
    if GP    ['mm_r_ri'].value<=0:
        print()
        print(GP       ['mm_r_ri'].value, GP['mm_r_or'].value , GP['mm_d_pm'].value , GP['mm_d_ri'].value)
        print('背铁太厚了 或 split_ration太小了！建议增大split_ratio是下限')
        print()
    return GP['mm_r_ri'].value 

class template_machine_as_numbers(object):
    ''' # template 有点类似analytical的电机（由几何尺寸组成）
    '''
    def __init__(self, fea_config_dict=None, spec_input_dict=None):

        # 仿真输入
        self.fea_config_dict = fea_config_dict
        self.spec_input_dict = spec_input_dict
        self.SI = SI = spec_input_dict # short name

        # 初始化搜索空间
        geometric_parameters = OrderedDict({
            # STATOR                          Type       Name                          Value  Bounds       Calc
            "deg_alpha_st" : acmop_parameter("free",    "stator_tooth_span_angle"    , None, [None, None], lambda GP,SI:None),
            "mm_d_so"      : acmop_parameter("free",    "stator_tooth_open_depth"    , None, [None, None], lambda GP,SI:None),
            "mm_d_st"      : acmop_parameter("free",    "stator_tooth_depth"         , None, [None, None], lambda GP,SI:None),
            "mm_r_os"      : acmop_parameter("free",    "outer_stator_radius"        , None, [None, None], lambda GP,SI:None),
            "mm_w_st"      : acmop_parameter("free",    "stator_tooth_width"         , None, [None, None], lambda GP,SI:None),
            "mm_r_si"      : acmop_parameter("derived", "inner_stator_radius"        , None, [None, None], lambda GP,SI:derive_mm_r_si     (GP,SI)),
            "deg_alpha_so" : acmop_parameter("derived", "stator_tooth_open_angle"    , None, [None, None], lambda GP,SI:derive_deg_alpha_so(GP,SI)),
            "mm_d_sp"      : acmop_parameter("derived", "stator_tooth_tip_depth"     , None, [None, None], lambda GP,SI:derive_mm_d_sp     (GP,SI)),
            "mm_d_sy"      : acmop_parameter("derived", "stator_yoke_depth"          , None, [None, None], lambda GP,SI:derive_mm_d_sy     (GP,SI)),
            # ROTOR
            "mm_d_fixed_air_gap": acmop_parameter("fixed",    "mechanical_air_gap_length",     None, [None, None], lambda GP,SI:None),
            "mm_d_sleeve"       : acmop_parameter("fixed",    "sleeve_length",                 None, [None, None], lambda GP,SI:None),
            "split_ratio"       : acmop_parameter("free",     "split_ratio_r_is_slash_r_os",   None, [None, None], lambda GP,SI:None),
            "mm_r_or"           : acmop_parameter("derived",  "outer_rotor_radius",            None, [None, None], lambda GP,SI:derive_mm_r_or(GP,SI)),
            "mm_r_ri"           : acmop_parameter("derived",  "inner_rotor_radius",            None, [None, None], lambda GP,SI:derive_mm_r_ri(GP,SI)),
        })
        # all in one place
        self.d = {
            "which_filter": fea_config_dict['which_filter'],
            "GP": geometric_parameters,
            "EX": OrderedDict(), # "Other Properties" is now renamed to "EXcitations"
            "x_denorm_dict": OrderedDict(),
            "bounds_denorm": [],
        }

        if 'FixedSleeveLength' in self.d['which_filter']:
            self.d['GP']['mm_d_sleeve'].type = "fixed"
        elif 'VariableSleeveLength' in self.d['which_filter']:
            self.d['GP']['mm_d_sleeve'].type = "free"
        elif 'VariableStatorSlotDepth_VariableStatorYokeDepth' in self.d['which_filter']:
            # IM
            self.d['GP']['mm_d_fixed_air_gap'].type = "free"
        else:
            raise Exception(f"Not defined: {self.d['which_filter']}")

        # Get Analytical Design
        # self.ModifiedBianchi2006(fea_config_dict, SI)


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

    ''' 初始化，定义优化
    '''
    def define_search_space(self, GP, original_template_neighbor_bounds):
        # 定义搜索空间，determine bounds
        self.bounds_denorm = []
        for key, val in original_template_neighbor_bounds.items():
            parameter = GP[key]
            parameter.bounds = val
            if parameter.type == 'free':
                self.bounds_denorm.append(val)
        print(f'[inner_rotor_motor.py] template BOUNDS_denorm in R^{len(self.bounds_denorm)}:', self.bounds_denorm)
        return self.bounds_denorm
    def get_other_properties_after_geometric_parameters_are_initialized(self, GP, SI):

        # Template's Other Properties (Shared by the swarm)
        EX = self.d['EX']
        if True:
            # WINDING Layout
            EX['wily'] = wily = winding_layout.winding_layout_v2(SI['DPNV_or_SEPA'], SI['Qs'], SI['p'], SI['ps'], SI['coil_pitch_y'])
            # STACK LENGTH
            EX['mm_template_stack_length'] = pyrhonen_procedure_as_function.get_mm_template_stack_length(SI, GP['mm_r_or'].value*1e-3) # mm TODO:
            EX['mm_mechanical_air_gap_length'] = SI['minimum_mechanical_air_gap_length_mm']
            # THERMAL Properties
            EX['Js']                = SI['Js'] # Arms/mm^2 im_OP['Js'] 
            EX['fill_factor']       = SI['space_factor_kCu'] # im_OP['fill_factor'] 
            # MOTOR Winding Excitation Properties
            EX['DriveW_zQ']         =            pyrhonen_procedure_as_function.get_zQ(SI, wily, GP['mm_r_si'].value*2*1e-3, GP['mm_r_or'].value*2*1e-3) # TODO:
            EX['DriveW_CurrentAmp'] = np.sqrt(2)*pyrhonen_procedure_as_function.get_stator_phase_current_rms(SI) # TODO:
            EX['DriveW_Freq']       = SI['ExcitationFreqSimulated']
            EX['DriveW_Rs']         = 1.0 # TODO: Must be greater than zero to let JMAG work
            EX['DriveW_poles']      = SI['p']*2
        # self.d.update( {"EX": EX} )
        return EX

    ''' 实用
    '''
    def get_rotor_volume(self, stack_length=None):
        if stack_length is None:
            return np.pi*(self.d['GP']['mm_r_or'].value*1e-3)**2 * (self.d['EX']['mm_template_stack_length']*1e-3)
        else:
            return np.pi*(self.d['GP']['mm_r_or'].value*1e-3)**2 * (stack_length*1e-3)
    def get_rotor_weight(self, gravity=9.8, stack_length=None):
        material_density_rho = pyrhonen_procedure_as_function.get_material_data()[0]
        if stack_length is None:
            return gravity * self.get_rotor_volume() * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg
        else:
            return gravity * self.get_rotor_volume(stack_length=stack_length) * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg

    ''' 玩弄几何变量
    '''
    def build_x_denorm(self):
        """ This is core function """
        # this is used in part_evaluation
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
    def get_x_denorm_dict_from_geometric_parameters(self, GP):
        x_denorm_dict = OrderedDict()
        for key, parameter in GP.items():
            if parameter.type == 'free':
                x_denorm_dict[key] = parameter.value
        return x_denorm_dict
    def get_x_denorm_dict_from_x_denorm_list(self, x_denorm):
        # 先拿个模板来，但是几何尺寸的变量值是旧的
        x_denorm_dict = self.get_x_denorm_dict_from_geometric_parameters(self.d['GP'])

        print('[inner_rotor_motor.py] DEBUG x_denorm_dict:', x_denorm_dict.keys())

        # 对模板进行遍历，挨个把新的几何尺寸的值从x_denorm中读取出来并更新x_denorm_dict
        for key, new_val in zip(x_denorm_dict.keys(), x_denorm):
            x_denorm_dict[key] = new_val
        return x_denorm_dict
    def update_geometric_parametes_using_x_denorm_dict(self, x_denorm_dict):
        # Update Free Parameters (a.k.a. x_denorm)
        for key, val in x_denorm_dict.items():
            self.d['GP'][key].value = val

        # Update Derived Parameters
        for key, parameter in self.d['GP'].items():
            if parameter.type == 'derived':
                # print(parameter)
                parameter.value = None # 先全部清空，防止编写derive时搞错依赖项关系，更新顺序是有先后的，后面的derived parameter 可以利用前面的derived parameter的值。
        # print()
        for key, parameter in self.d['GP'].items():
            if parameter.type == 'derived':
                parameter.value = parameter.calc(self.d['GP'], self.SI)
                # print(parameter)
                if parameter.value<=0:
                    raise Exception('Error: Negative derived parameter', str(parameter))
        return self.d['GP']

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
        self.ID = 'Q%dp%ds%d'%(self.template.SI['Qs'], self.template.SI['p'],self.template.SI['no_segmented_magnets'])

        #02 Geometry Data
        if x_denorm is None:
            # template as variant
            GP = self.template.d['GP'] # do nothing, use template's GP
        else:
            x_denorm_dict = self.template.get_x_denorm_dict_from_x_denorm_list(x_denorm)
            if verbose:
                for k,v in x_denorm_dict.items():
                    print('\t', k,v)
            GP = self.template.update_geometric_parametes_using_x_denorm_dict(x_denorm_dict)

        print('[inner_rotor_motor.py] DEBUG x_denorm length is', len(x_denorm), x_denorm)

        #test: yes, it is by reference!
        # print(self.template.d['GP'])
        # GP['aaaaaa'] = 1239129
        # print(self.template.d['GP'])
        # quit()

        #03 Inherit properties
        # self.template.d['EX']
            # self.DriveW_CurrentAmp = None ########### will be assisned when drawing the coils
            # self.DriveW_poles      = spmsm_template.DriveW_poles

            # self.Js                = spmsm_template.Js 
            # self.fill_factor       = spmsm_template.fill_factor 

            # self.template.d['EX']['mm_template_stack_length']      = spmsm_template.stack_length 


        # #04 Material Condutivity Properties
        # if self.template.fea_config_dict is not None:
        #     self.End_Ring_Resistance = fea_config_dict['End_Ring_Resistance']
        #     self.Bar_Conductivity = fea_config_dict['Bar_Conductivity']
        # self.Copper_Loss = self.DriveW_CurrentAmp**2 / 2 * self.DriveW_Rs * 3
        # # self.Resistance_per_Turn = 0.01 # TODO

        #05 Winidng Excitation
        # self.CurrentAmp_per_phase = None # will be used in copper loss calculation

        #06 Meshing & Solver Properties
        # self.max_nonlinear_iteration = 50 # 30 for transient solve
        # self.meshSize_Magnet = 2 # mm

    def reproduce_wily(self):
        ''' This method is only used for reproducing design from jsonpickle'''
        self.template.d['EX']['wily'] = winding_layout.winding_layout_v2(self.template.SI['DPNV_or_SEPA'], self.template.SI['Qs'], self.template.SI['p'], self.template.SI['ps'], self.template.SI['coil_pitch_y'])

    def update_mechanical_parameters(self, syn_freq=None):
        EX = self.template.d['EX']
        if syn_freq is None:
            EX['the_speed'] = EX['DriveW_Freq']*60. / (0.5*EX['DriveW_poles']) # rpm
            EX['Omega']     = EX['the_speed'] / 60. * 2*np.pi
            # self.omega = None # This variable name is devil! you can't tell its electrical or mechanical! #+ self.DriveW_Freq * (1-self.the_slip) * 2*pi
        else:
            raise Exception('Not implemented.')

    ''' BELOW is for JMAG Designer
    '''
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
            EX = self.template.d['EX']
            CurrentAmp_in_the_slot = self.coils.mm2_slot_area * EX['fill_factor'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
            CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                # try:
                #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                # except AttributeError:
                #     # print(EX['wily'])
                #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily']['number_parallel_branch']
                #     print("[inner_rotor_motor.py] Reproduce design using jsonpickle will encounter error here: 'dict' object has no attribute 'number_parallel_branch', implying that the object wily has become a dict after jsonpickle.")
                #     # quit() 

            # Maybe there is a bug here... regarding the excitation for suspension winding...
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
            EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
            EX['DriveW_CurrentAmp'] = self.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            EX['BeariW_CurrentAmp'] = self.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            # self.spec_geometry_dict['DriveW_CurrentAmp'] = self.DriveW_CurrentAmp

            slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
            print('[inner_rotor_motor.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

            # print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
            # print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
            # print('---self.DriveW_CurrentAmp =', self.DriveW_CurrentAmp)
            # print('---self.BeariW_CurrentAmp =', self.BeariW_CurrentAmp)
            # print('---TORQUE_CURRENT_RATIO:', self.template.fea_config_dict['TORQUE_CURRENT_RATIO'])
            # print('---SUSPENSION_CURRENT_RATIO:', self.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            # Import Model into Designer
            toolJd.save(self.name, self.show(toString=True))

        return True

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
            print(', \n\t'.join("%s = %s" % item for item in variant_tuple_list))
            print(', \n\t'.join("%s = %s" % item for item in template_tuple_list))
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

