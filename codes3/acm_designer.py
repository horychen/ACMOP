from time import time as clock_time
from pylab import np # plt, mpl
import os, logging, win32com.client, json
import utility, utility_json
import population, FEMM_Solver, pyrhonen_procedure_as_function

# BPMSM codes
import JMAG, FEMM_SlidingMesh
import bearingless_spmsm_design, vernier_motor_design, flux_alternator_design, flux_switching_pm_design

class Swarm_Data_Analyzer(object):
    def __init__(self, fname, desired_x_denorm_dict):
        if not os.path.exists(fname):
            self.number_of_chromosome = 0
            self.swarm_data_xf = None
        else:
            ''' 1. Load json file
            '''
            with open(fname, 'r') as f:
                buf = f.read()
                swarm_data_as_dict = json.loads('{'+buf[1:]+'}')
                del buf
                if 'Test' in swarm_data_as_dict.keys():
                    del swarm_data_as_dict['Test']
            self.number_of_chromosome = len(swarm_data_as_dict)
            self.swarm_data_as_dict = swarm_data_as_dict

            ''' 2. Get swarm_data_xf
            '''
            # for k1, v1 in swarm_data_as_dict.items():
            #     for k2, v2 in v1.items():
            #         print(k1, k2, v2)

            # DEBUG
            # print('[acm_designer.py] Archive order:', list(list(swarm_data_as_dict.values())[0].values())[0]['x_denorm_dict'].keys() )
            # print('\tThe desired order is:', desired_x_denorm_dict.keys())
            # print('---'*30)
            def sort_as_desired(x_denorm_dict, desired_x_denorm_dict=None):
                if desired_x_denorm_dict is None:
                    return list(x_denorm_dict.values())
                else:
                    try:
                        return [x_denorm_dict[key] for key in desired_x_denorm_dict.keys()]
                    except KeyError as e:
                        print('Error Hint: some geometric parameters are renamed so the old json archive file now has different name from the new name.')
                        raise e

                ''' 目前只支持重新跑优化的时候减少自由的几何参数，如果要增加自由的几何参数，则要从GP里面拿出来对应的参数的取值。其实也很简单啦。'''
                x_denorm = []
                for key in desired_x_denorm_dict.keys():
                    if key not in x_denorm_dict:
                        x_denorm.append(self.decode(v)['Geometric parameters'][key]) # pseudo code for showing the concept, this will not work.
                    else:
                        x_denorm.append(x_denorm_dict[key])
            self.swarm_data_xf = [
                                    sort_as_desired(self.decode(v)['x_denorm_dict'], desired_x_denorm_dict) + [ self.decode(v)['Performance']['f1'], self.decode(v)['Performance']['f2'], self.decode(v)['Performance']['f3'] ]
                                    for v in swarm_data_as_dict.values() # v = {name:data}
                                    ]
            self.number_of_free_variables = len(self.swarm_data_xf[0]) - 3

            # DEBUG
            # print('[acm_designer.py]')
            # for ind, xf in enumerate(self.swarm_data_xf):
            #     print(f'{ind:04d}', ',\t'.join([f'{el:.2f}' for el in xf]))

            ''' 3. Get the list of other attribute by individuals (not needed for optimization)
            '''
                # self.swarm_data_project_names = [ self.decode(v)['Performance']['project_name'] for v in swarm_data_as_dict.values() ]
                # self.prepare_data_for_post_processing(swarm_data_as_dict)
            self.swarm_data_project_names = self.get_metric_of_the_whole_swarm('project_name')

    @staticmethod
    def decode(d):
        return list(d.values())[0]

    def get_metric_of_the_whole_swarm(self, metric):
        return [ self.decode(v)['Performance'][metric] for v in self.swarm_data_as_dict.values() ]
    def prepare_data_for_post_processing(self):

        ''' 3. Get the list of other attribute by individuals (not needed for optimization)
        '''
        # self.machine_data.append([float(x) for x in raw[3].split(',')])
        # self.rated_data.append(  [float(x) for x in raw[4].split(',')])
        self.FRW = self.get_metric_of_the_whole_swarm('FRW')
        self.Em = self.get_metric_of_the_whole_swarm('normalized_force_error_magnitude')
        self.Ea = self.get_metric_of_the_whole_swarm('force_error_angle')
        self.Tripple = self.get_metric_of_the_whole_swarm('normalized_torque_ripple')
        # self.RatedVol = self.get_metric_of_the_whole_swarm('rated_rotor_volume')
        # self.RatedWeight = self.get_metric_of_the_whole_swarm('rated_rotor_weight')
        self.RatedStkLen = self.get_metric_of_the_whole_swarm('rated_stack_length_mm')

        # self.f1 = [raw[-3] for raw in self.swarm_data_xf]
        # self.f2 = [raw[-2] for raw in self.swarm_data_xf]
        # self.f3 = [raw[-1] for raw in self.swarm_data_xf]

        # [power_factor, efficiency, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
        self.PowerFactor            = self.get_metric_of_the_whole_swarm('power_factor')
        try:
            self.Cost                = self.get_metric_of_the_whole_swarm('Cost')
        except:
            self.Cost                = np.array(self.get_metric_of_the_whole_swarm('f1'))
        try:
            self.TRV                = self.get_metric_of_the_whole_swarm('TRV')
        except:
            self.TRV                = - np.array(self.get_metric_of_the_whole_swarm('f1'))
        try:
            self.RatedEfficiency    = self.get_metric_of_the_whole_swarm('RatedEfficiency')
        except:
            self.RatedEfficiency    = - np.array(self.get_metric_of_the_whole_swarm('f2'))
        self.TorqueRipple           = self.get_metric_of_the_whole_swarm('normalized_torque_ripple')
        # self.torque_average                   = [raw[2] for raw in self.machine_data]
        # self.ss_avg_force_magnitude           = [raw[4] for raw in self.machine_data]
        # self.normalized_force_error_magnitude = [raw[5] for raw in self.machine_data]
        # self.force_error_angle                = [raw[6] for raw in self.machine_data]

        # self.l_rated_shaft_power                    = [raw[0] for raw in self.rated_data]
        # self.l_rated_efficiency                     = [raw[1] for raw in self.rated_data]
        self.l_rated_total_loss                     = self.get_metric_of_the_whole_swarm('rated_total_loss') # [raw[2] for raw in self.rated_data]
        self.l_rated_stator_copper_loss_along_stack = self.get_metric_of_the_whole_swarm('rated_stator_copper_loss_along_stack') # [raw[3] for raw in self.rated_data]
        self.l_rated_rotor_copper_loss_along_stack  = self.get_metric_of_the_whole_swarm('rated_rotor_copper_loss_along_stack') # [raw[4] for raw in self.rated_data]
        self.l_stator_copper_loss_in_end_turn       = self.get_metric_of_the_whole_swarm('stator_copper_loss_in_end_turn') # [raw[5] for raw in self.rated_data]
        self.l_rotor_copper_loss_in_end_turn        = self.get_metric_of_the_whole_swarm('rotor_copper_loss_in_end_turn') # [raw[6] for raw in self.rated_data]
        self.l_rated_iron_loss                      = self.get_metric_of_the_whole_swarm('rated_iron_loss') # [raw[7] for raw in self.rated_data]
        self.l_rated_windage_loss                   = self.get_metric_of_the_whole_swarm('rated_windage_loss') # [raw[8] for raw in self.rated_data]
        # self.l_rated_magnet_Joule_loss              = self.get_metric_of_the_whole_swarm('rated_magnet_Joule_loss')
        # self.l_rated_rotor_volume                   = self.get_metric_of_the_whole_swarm('rated_rotor_volume') # [raw[9] for raw in self.rated_data]
        # self.l_rated_rotor_weight                   = self.get_metric_of_the_whole_swarm('rated_rotor_weight') # [(V*8050*9.8) for V in self.l_rated_rotor_volume] # density of rotor is estimated to be that of steraw of 8050 g/cm^3
        self.l_rated_stack_length                   = self.get_metric_of_the_whole_swarm('rated_stack_length_mm') # [raw[10] for raw in self.rated_data] # new!
        # self.l_original_stack_length                = [raw[11] for raw in self.rated_data] # new!
        # self.l_original_rotor_weight                = [weight/rated*ori for weight, rated, ori in zip(self.l_rated_rotor_weight, self.l_rated_stack_length, self.l_original_stack_length)]

        # TODO: change to EX['mec_power'] and EX['the_speed']
        # required_torque = 50e3 / (30000/60*2*np.pi)         # TODO: should use rated stack length and torque average to compute this
        # self.l_TRV = [required_torque/raw for raw in self.l_rated_rotor_volume]
        # self.l_FRW = [F/W for W, F in zip(self.l_original_rotor_weight, self.l_ss_avg_force_magnitude)] # FRW
        pass

class swarm_data_container(object):
    def __init__(self, swarm_data_raw, fea_config_dict, swarm_data_json=None):

        self.swarm_data_raw = swarm_data_raw
        self.fea_config_dict = fea_config_dict

        # x, f(x)
        self.swarm_data_xf = []
        self.project_names = []
        self.machine_data = []
        self.rated_data = []
        self.Trip = []
        self.FRW = []
        self.Em = []
        self.Ea = []
        self.RatedVol = []
        self.RatedWeight = []
        self.RatedStkLen = []
        #IM 
            # if len(bound_filter) == 9: # This is induction motor
            #     for raw in swarm_data_raw:

            #         design_parameters_denorm = [float(x) for x in raw[5].split(',')]
            #         # print(design_parameters_denorm, len(design_parameters_denorm))
            #         # quit()

            #         loc1 = raw[2].find('f1')
            #         loc2 = raw[2].find('f2')
            #         loc3 = raw[2].find('f3')
            #         f1 = float(raw[2][loc1+3:loc2-1])
            #         f2 = float(raw[2][loc2+3:loc3-1])
            #         f3 = float(raw[2][loc3+3:])

            #         x_denorm = self.get_x_denorm_from_design_parameters(design_parameters_denorm, bound_filter)
            #         self.swarm_data_xf.append(x_denorm + [f1, f2, f3])
            #         # print(self.swarm_data_xf)
            #         # quit()

            #         self.project_names.append(raw[1][:-1])
            #         self.machine_data.append([float(x) for x in raw[3].split(',')])
            #         self.rated_data.append(  [float(x) for x in raw[4].split(',')])

            #         individual_Trip = [float(x) for x in raw[3].split(',')][3]
            #         self.Trip.append(individual_Trip)

            #         # Get FRW
            #         individual_ss_avg_force_magnitude = [float(x) for x in raw[3].split(',')][4]
            #         individual_Em                     = [float(x) for x in raw[3].split(',')][5]
            #         individual_Ea                     = [float(x) for x in raw[3].split(',')][6]
            #         individual_rated_rotor_volume     = [float(x) for x in raw[4].split(',')][9]
            #         individual_rated_rotor_weight     = (individual_rated_rotor_volume*8050*9.8)
            #         individual_rated_stack_length     = [float(x) for x in raw[4].split(',')][10]
            #         individual_original_stack_length  = [float(x) for x in raw[4].split(',')][11]
            #         individual_original_rotor_weight  = individual_rated_rotor_weight/individual_rated_stack_length*individual_original_stack_length
            #         individual_FRW = individual_ss_avg_force_magnitude/individual_original_rotor_weight
            #         self.FRW.append(individual_FRW)
            #         self.Em.append(individual_Em)
            #         self.Ea.append(individual_Ea)
            #         self.RatedVol.append(individual_rated_rotor_volume)
            #         self.RatedWeight.append(individual_rated_rotor_weight)
            #         self.RatedStkLen.append(individual_rated_stack_length)
        # else: # This is PM motor
        if True:
            # self.swarm_data_raw = swarm_data_raw
            # self.fea_config_dict = fea_config_dict

            # x, f(x)
            # self.swarm_data_xf = []
            # self.project_names = []

            # self.machine_data = []
            # self.rated_data = []
            # self.Trip = []
            # self.FRW = []
            # self.Em = []
            # self.Ea = []
            # self.RatedVol = []
            # self.RatedWeight = []
            # self.RatedStkLen = []
            self.deg_alpha_st = []
            self.mm_w_st = []
            self.mm_r_si = []

            # TODO: use swarm_data_json over raw
            if swarm_data_json is not None:
                for key in swarm_data_json.keys():
                    # print(swarm_data_json[key])
                    print('[acm_designer.py] DEBUG (swarm_data_json)', list(swarm_data_json[key].keys()))
                    # print('DEBUG', list(swarm_data_json[key].values()))

                    the_variant_dict = list(swarm_data_json[key].values())
                    x_denorm = list( the_variant_dict[0]['x_denorm_dict'].values() )
                    # x_denorm = [val for val in list(swarm_data_json[key].values())['x_denorm_dict'].items()]
                    # print(x_denorm)
                    # quit()
            if True:
                for raw in self.swarm_data_raw:

                    # spmsm_template.design_parameters = [
                    #                                   0 spmsm_template.deg_alpha_st 
                    #                                   1 spmsm_template.deg_alpha_sto 
                    #                                   2 spmsm_template.mm_r_si      
                    #                                   3 spmsm_template.mm_d_sto      
                    #                                   4 spmsm_template.mm_d_stt      
                    #                                   5 spmsm_template.mm_d_st      
                    #                                   6 spmsm_template.mm_d_sy      
                    #                                   7 spmsm_template.mm_w_st      
                    #                                   8 spmsm_template.mm_r_st      
                    #                                   9 spmsm_template.mm_r_sf      
                    #                                  10 spmsm_template.mm_r_sb      
                    #                                  11 spmsm_template.Q            
                    #                                  12 spmsm_template.sleeve_length
                    #                                  13 spmsm_template.fixed_air_gap_length
                    #                                  14 spmsm_template.mm_d_pm      
                    #                                  15 spmsm_template.deg_alpha_rm 
                    #                                  16 spmsm_template.deg_alpha_rs 
                    #                                  17 spmsm_template.mm_d_ri      
                    #                                  18 spmsm_template.mm_r_ri      
                    #                                  19 spmsm_template.mm_d_rp      
                    #                                  20 spmsm_template.mm_d_rs      
                    #                                  21 spmsm_template.p
                    #                                  22 spmsm_template.s
                    #                                 ]
                    design_parameters_denorm = [float(x) for x in raw[5].split(',')]
                    self.deg_alpha_st.append(design_parameters_denorm[0] )
                    self.mm_w_st.append(     design_parameters_denorm[7] )
                    self.mm_r_si.append(   design_parameters_denorm[2])

                    loc1 = raw[2].find('f1')
                    loc2 = raw[2].find('f2')
                    loc3 = raw[2].find('f3')
                    f1 = float(raw[2][loc1+3:loc2-1])
                    f2 = float(raw[2][loc2+3:loc3-1])
                    f3 = float(raw[2][loc3+3:])

                    if len(design_parameters_denorm) > 20:
                        ''' 永磁电机 复古 '''

                        # 在 acmop 里，我们已经放弃了使用 bound_filter 的概念。
                        # x_denorm = self.get_x_denorm_from_design_parameters(design_parameters_denorm, bound_filter)

                        """ This is consistent with bopt-python """
                        # x_denorm = [None]*11
                        # x_denorm[0]  = design_parameters_denorm[0] # spmsm_template.deg_alpha_st 
                        # x_denorm[1]  = design_parameters_denorm[3] # spmsm_template.mm_d_sto         
                        # x_denorm[2]  = design_parameters_denorm[5] # spmsm_template.mm_d_st
                        # x_denorm[3]  = design_parameters_denorm[7] # spmsm_template.mm_w_st         
                        # x_denorm[4]  = design_parameters_denorm[12] # spmsm_template.sleeve_length   
                        # x_denorm[5]  = design_parameters_denorm[14] # spmsm_template.mm_d_pm         
                        # x_denorm[6]  = design_parameters_denorm[15] # spmsm_template.deg_alpha_rm    
                        # x_denorm[7]  = design_parameters_denorm[16] # spmsm_template.deg_alpha_rs    
                        # x_denorm[8]  = design_parameters_denorm[17] # spmsm_template.mm_d_ri         
                        # x_denorm[9]  = design_parameters_denorm[19] # spmsm_template.mm_d_rp         
                        # x_denorm[10] = design_parameters_denorm[20] # spmsm_template.mm_d_rs         

                        """ This is consistent with ACMOP """
                        x_denorm = [None]*11
                        x_denorm[0]  = design_parameters_denorm[0] # spmsm_template.deg_alpha_st 
                        x_denorm[1]  = design_parameters_denorm[3] # spmsm_template.mm_d_sto         
                        x_denorm[2]  = design_parameters_denorm[5] # spmsm_template.mm_d_st
                        x_denorm[3]  = sum([design_parameters_denorm[i] for i in (2,4,5,6)]) # outer_stator_radius mm_r_so
                        x_denorm[4]  = design_parameters_denorm[7] # spmsm_template.mm_w_st   
                        x_denorm[5]  = design_parameters_denorm[12] #            mm_d_sleeve
                        r_si = design_parameters_denorm[2] # 2 spmsm_template.mm_r_si      
                        try:
                            x_denorm[6]  = r_si /  x_denorm[3] # split_ratio     r_is_slash_r_os 
                        except ZeroDivisionError as e:
                            print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                            print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                            print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                            raise e
                        x_denorm[7]  = design_parameters_denorm[14] # spmsm_template.mm_d_pm      
                        x_denorm[8]  = design_parameters_denorm[17] # spmsm_template.mm_d_ri         
                        # childGP
                        x_denorm[9]  = design_parameters_denorm[15] # spmsm_template.deg_alpha_rm    
                        x_denorm[10]  = design_parameters_denorm[19] # spmsm_template.mm_d_rp         
                        # x_denorm[11]  = design_parameters_denorm[16] # spmsm_template.deg_alpha_rs    
                        # x_denorm[12] = design_parameters_denorm[20] # spmsm_template.mm_d_rs         

                        # DEBUG
                        # odict_keys(['deg_alpha_st', 'mm_d_sto', 'mm_d_st', 'mm_r_so', 'mm_w_st', 'mm_d_sleeve', 'split_ratio', 'mm_d_pm', 'mm_d_ri', 'deg_alpha_rm', 'mm_d_rp'])
                        # deg_alpha_st 11.1183
                        # mm_d_sto 1.50079
                        # mm_d_st 42.9701
                        # mm_r_so 16.099
                        # mm_w_st 5.89091
                        # mm_d_sleeve 5.19948
                        # split_ratio 44.9638
                        # mm_d_pm 44.9638
                        # mm_d_ri 3.67901
                        # deg_alpha_rm 5.19948
                        # mm_d_rp 0.0

                    else:
                        '''感应电机 复古 '''
                        raise Exception('not implemented')

                    # print(design_parameters_denorm, f1, f2, f3)
                    # THERE IS A BUT HERE: slot_tip_open_ratio is less than 0.2---Not possible
                        # free_variables[0]  = design_parameters[0] # spmsm_template.deg_alpha_st 
                        # free_variables[4]  = design_parameters[7] # spmsm_template.mm_w_st         
                        # free_variables[10] = sum([design_parameters[i] for i in (18,17,19)]) # spmsm_template.mm_r_ri + spmsm_template.mm_d_ri + spmsm_template.mm_d_rp
                        # self.deg_alpha_st.append(x_denorm[0] ) 
                        # self.mm_w_st.append(x_denorm[4] ) 
                        # self.mm_radius.append(x_denorm[10]) 

                    self.project_names.append(raw[1][:-1])
                    self.machine_data.append([float(x) for x in raw[3].split(',')])
                    self.rated_data.append(  [float(x) for x in raw[4].split(',')])

                    individual_Trip = [float(x) for x in raw[3].split(',')][3]
                    self.Trip.append(individual_Trip)

                    # Get FRW
                    individual_ss_avg_force_magnitude = [float(x) for x in raw[3].split(',')][4]
                    individual_Em                     = [float(x) for x in raw[3].split(',')][5]
                    individual_Ea                     = [float(x) for x in raw[3].split(',')][6]
                    individual_rated_rotor_volume     = [float(x) for x in raw[4].split(',')][9]
                    individual_rated_rotor_weight     = (individual_rated_rotor_volume*8050*9.8)
                    individual_rated_stack_length     = [float(x) for x in raw[4].split(',')][10]
                    individual_original_stack_length  = [float(x) for x in raw[4].split(',')][11]
                    individual_original_rotor_weight  = individual_rated_rotor_weight/individual_rated_stack_length*individual_original_stack_length
                    individual_FRW = individual_ss_avg_force_magnitude/individual_original_rotor_weight
                    self.FRW.append(individual_FRW)
                    self.Em.append(individual_Em)
                    self.Ea.append(individual_Ea)
                    self.RatedVol.append(individual_rated_rotor_volume)
                    self.RatedWeight.append(individual_rated_rotor_weight)
                    self.RatedStkLen.append(individual_rated_stack_length)

                    # Add FRW to xf (This will cause re-starting error)
                    # self.swarm_data_xf.append(x_denorm + [individual_FRW, f1, f2, f3])

                    self.swarm_data_xf.append(x_denorm + [f1, f2, f3])

        self.number_of_free_variables = len(x_denorm)
        print('\tCount of individuals:', len(self.swarm_data_raw))

        self.l_OA = [raw[-3] for raw in self.swarm_data_xf]
        self.l_OB = [raw[-2] for raw in self.swarm_data_xf]
        self.l_OC = [raw[-1] for raw in self.swarm_data_xf]
        self.l_design_parameters = [raw[:-3] for raw in self.swarm_data_xf]

        # [power_factor, efficiency, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
        self.l_power_factor                     = [raw[0] for raw in self.machine_data]
        self.l_efficiency                       = [raw[1] for raw in self.machine_data]
        self.l_torque_average                   = [raw[2] for raw in self.machine_data]
        self.l_normalized_torque_ripple         = [raw[3] for raw in self.machine_data]
        self.l_ss_avg_force_magnitude           = [raw[4] for raw in self.machine_data]
        self.l_normalized_force_error_magnitude = [raw[5] for raw in self.machine_data]
        self.l_force_error_angle                = [raw[6] for raw in self.machine_data]

        self.l_rated_shaft_power                    = [raw[0] for raw in self.rated_data]
        self.l_rated_efficiency                     = [raw[1] for raw in self.rated_data]
        self.l_rated_total_loss                     = [raw[2] for raw in self.rated_data]
        self.l_rated_stator_copper_loss_along_stack = [raw[3] for raw in self.rated_data]
        self.l_rated_rotor_copper_loss_along_stack  = [raw[4] for raw in self.rated_data]
        self.l_stator_copper_loss_in_end_turn       = [raw[5] for raw in self.rated_data]
        self.l_rotor_copper_loss_in_end_turn        = [raw[6] for raw in self.rated_data]
        self.l_rated_iron_loss                      = [raw[7] for raw in self.rated_data]
        self.l_rated_windage_loss                   = [raw[8] for raw in self.rated_data]
        self.l_rated_rotor_volume                   = [raw[9] for raw in self.rated_data]
        self.l_rated_rotor_weight                   = [(V*8050*9.8) for V in self.l_rated_rotor_volume] # density of rotor is estimated to be that of steraw of 8050 g/cm^3
        self.l_rated_stack_length                   = [raw[10] for raw in self.rated_data] # new!
        self.l_original_stack_length                = [raw[11] for raw in self.rated_data] # new!
        self.l_original_rotor_weight                = [weight/rated*ori for weight, rated, ori in zip(self.l_rated_rotor_weight, self.l_rated_stack_length, self.l_original_stack_length)]

        # TODO: change to EX['mec_power'] and EX['the_speed']
        required_torque = 50e3 / (30000/60*2*np.pi)         # TODO: should use rated stack length and torque average to compute this
        self.l_TRV = [required_torque/raw for raw in self.l_rated_rotor_volume]
        self.l_FRW = [F/W for W, F in zip(self.l_original_rotor_weight, self.l_ss_avg_force_magnitude)] # FRW

    def get_list_y_data(self):

        list_y_data = [ self.l_TRV, ##self.l_rated_stack_length,
                        [100*raw for raw in self.l_OB], 
                        self.l_force_error_angle,
                        ]
        return list_y_data

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Utility
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

    def get_x_denorm_from_design_parameters(self, design_parameters, bound_filter=None):
        if bound_filter is None:
            x_denorm = design_parameters
            return x_denorm

        if len(bound_filter) == 13:
            # step 1: get free_variables from design_parameters
            free_variables = [None]*13
            free_variables[0]  = design_parameters[0] # spmsm_template.deg_alpha_st 
            free_variables[1]  = design_parameters[3] # spmsm_template.mm_d_sto         
            free_variables[2]  = design_parameters[5] # spmsm_template.mm_d_st
            free_variables[3]  = sum([design_parameters[i] for i in (2,4,5,6)]) # spmsm_template.mm_r_si + spmsm_template.mm_d_stt + spmsm_template.mm_d_st + spmsm_template.mm_d_sy # stator outer radius
            free_variables[4]  = design_parameters[7] # spmsm_template.mm_w_st         
            free_variables[5]  = design_parameters[12] # spmsm_template.sleeve_length   
            free_variables[6]  = design_parameters[14] # spmsm_template.mm_d_pm         
            free_variables[7]  = design_parameters[15] # spmsm_template.deg_alpha_rm    
            free_variables[8]  = design_parameters[16] # spmsm_template.deg_alpha_rs    
            free_variables[9]  = design_parameters[17] # spmsm_template.mm_d_ri         
            free_variables[10] = sum([design_parameters[i] for i in (18,17,19)]) # spmsm_template.mm_r_ri + spmsm_template.mm_d_ri + spmsm_template.mm_d_rp -> rotor_outer_steel_radius
            free_variables[11] = design_parameters[19] # spmsm_template.mm_d_rp         
            free_variables[12] = design_parameters[20] # spmsm_template.mm_d_rs         
        elif len(bound_filter) == 9:
            free_variables = design_parameters # For IM, free_variables are design_parameters (even always having the same length)
            # print(free_variables)

        # step 2: get x_denorm from free_variables
        x_denorm = []
        for idx, boo in enumerate(bound_filter):
            if boo == 1:
                # print(idx)
                x_denorm.append( free_variables[idx] )
        return x_denorm

    def sensitivity_bar_charts(self):
        number_of_variant = self.fea_config_dict['local_sensitivity_analysis_number_of_variants'] + 1
        number_of_free_variables = self.number_of_free_variables

        from pylab import subplots, mpl, plt
        mpl.style.use('classic')
        mpl.rcParams['legend.fontsize'] = 12
        # mpl.rcParams['legend.family'] = 'Times New Roman'
        mpl.rcParams['font.family'] = ['Times New Roman']
        # mpl.rcParams['font.size'] = 15.0
        font = {'family' : 'Times New Roman', #'serif',
                'color' : 'darkblue',
                'weight' : 'normal',
                'size' : 14,}
        textfont = {'family' : 'Times New Roman', #'serif',
                    'color' : 'darkblue',
                    'weight' : 'normal',
                    'size' : 11.5,}

        fig, axeses = subplots(4, 2, sharex=True, dpi=150, figsize=(16*0.75, 8*0.75), facecolor='w', edgecolor='k', constrained_layout=True)
        ax_list = []
        for i in range(4):
            ax_list.extend(axeses[i].tolist())
        # O2_prototype_ax.plot(O2_prototype_data[1], 'o-', lw=0.75, alpha=0.5, label=r'$\delta$'         )
        # O2_prototype_ax.plot(O2_prototype_data[0], 'v-', lw=0.75, alpha=0.5, label=r'$b_{\rm tooth,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[3], 's-', lw=0.75, alpha=0.5, label=r'$b_{\rm tooth,r}$')
        # O2_prototype_ax.plot(O2_prototype_data[5], '^-', lw=0.75, alpha=0.5, label=r'$w_{\rm open,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[2], 'd-', lw=0.75, alpha=0.5, label=r'$w_{\rm open,r}$')
        # O2_prototype_ax.plot(O2_prototype_data[6], '*-', lw=0.75, alpha=0.5, label=r'$h_{\rm head,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[4], 'X-', lw=0.75, alpha=0.5, label=r'$h_{\rm head,r}$')

        # Extract data
        free_param_list = [
        r'$L_g$',        
        r'$w_{st}$',     
        r'$w_{rt}$',     
        r'$\theta_{so}$',
        r'$w_{ro}$',
        r'$d_{so}$',
        r'$d_{ro}$']
        y_label_list = ['PF', r'$\eta$ [100%]', r'$T_{em} [N]$', r'$T_{rip}$ [100%]', r'$|F|$ [N]', r'$E_m$ [100%]', r'$E_a$ [deg]', 
                        r'$P_{Cu,s,JMAG}$', r'$P_{Cu,r,JMAG}$', r'$P_{Fe}$ [W]', r'$P_{eddy}$', r'$P_{hyst}$', r'$P_{Cu,s,FEMM}$', r'$P_{Cu,r,FEMM}$', 
                        r'Windage loss', r'Total loss']

        list_y_label = [r'$O_A$ [$\rm kNm/m^3$]', 
                         '$O_C$ [1]', 
                         'FRW [p.u.]',
                         '$E_a$ [deg]', 
                         '$O_B$ [%]', 
                         '$E_m$ [%]', 
                         r'$P_{\rm loss}$ [W]',
                         r'$T_{\rm rip}$ [%]',
                         # 'Rotor Weight [N]', #'Power Factor [1]',
                         ]
        list_y_data_max = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        list_y_data_min = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

        list_y_data = [ [el/1e3 for el in self.l_OA], 
                        self.l_OC,
                        [F/W for W, F in zip(self.l_original_rotor_weight, self.l_ss_avg_force_magnitude)], # FRW
                        self.l_force_error_angle,
                        [100*el for el in self.l_OB], 
                        [100*el for el in self.l_normalized_force_error_magnitude],
                        self.l_rated_total_loss,
                        [100*el for el in self.l_normalized_torque_ripple],
                        # self.l_rated_rotor_weight,
                        ]
        for i in range(len(list_y_label)):
            ax = ax_list[i]
            y_data = list_y_data[i]
            y_value_reference = y_data[0]
            ax.plot(y_value_reference*np.ones(len(y_data)), '-k', alpha=1, zorder=10)
            y_data = y_data[1:]
            number_of_points_per_geometry_variable = len(y_data)/number_of_free_variables
            for idx, part_of_y_data in enumerate([y_data[int(number_of_points_per_geometry_variable*_)\
                                                        :int(number_of_points_per_geometry_variable*(_+1))]\
                                                        for _ in range(number_of_free_variables)]):
                if idx%2 == 0:
                    line_style = '--bo'
                else:
                    line_style = '--ro'
                ax.plot(list(range(len(y_data)))[int(number_of_points_per_geometry_variable*idx)\
                                                :int(number_of_points_per_geometry_variable*(idx+1))], 
                                                part_of_y_data, line_style, alpha=0.33)

            low, high = ax.get_ylim()
            # ax.legend()
            ax.grid()
            ax.set_ylabel(list_y_label[i], **font)
            ax.set_xlim([0,140])
            for j in range(number_of_free_variables):
                if j%2==0:
                    alpha = 0.05
                else:
                    alpha = 0.15
                ax.axvspan(j*number_of_variant-0.5, (j+1)*number_of_variant-0.5, facecolor='k', alpha=alpha)
                ax.text(0.33*number_of_free_variables+j*number_of_variant, high-(high-low)*0.125, free_param_list[j])

            if i == 0:
                ax.set_yticks([-24, -22, -20, -18, -16])
            list_y_data_max[i].append(max(y_data))
            list_y_data_min[i].append(min(y_data))

        ax_list[-2].set_xlabel('Count of design variant', **font)
        ax_list[-1].set_xlabel('Count of design variant', **font)
        fig.savefig(r'C:\Users\horyc\Desktop/'+ 'LSA_curves.png', dpi=300)
        # plt.show()
        return













        INDEX_TOTAL_LOSS = 15 + 4 # index of total loss in the machine_data list

        # ------------------------------------ Sensitivity Analysis Bar Chart Scripts

        # print next(self.get_list_objective_function())
        data_max = []
        data_min = []
        eta_at_50kW_max = []
        eta_at_50kW_min = []
        O1_max   = []
        O1_min   = []
        for ind, i in enumerate(list(range(7))+[INDEX_TOTAL_LOSS]):
            print('\n-----------', y_label_list[i])
            l = list(self.get_certain_objective_function(i))
            y = l
            print('ind=', ind, 'i=', i, 'len(y)=', len(y))

            data_max.append([])
            data_min.append([])

            for j in range(int(len(y)/number_of_variant)): # iterate design parameters
                y_vs_design_parameter = y[j*number_of_variant:(j+1)*number_of_variant]

                try:
                    # if j == 6:
                    ax_list[ind].plot(y_vs_design_parameter, 'o-', lw=0.75, label=str(j)+' '+param_list[j], alpha=0.5)
                except IndexError as e:
                    print('Check the length of y should be 7*(%d+1)=%d, or else you should remove the redundant results in swarm_data.txt (they are produced because of the interrupted/resumed script run.)'%(number_of_variant, 7*number_of_variant))
                    raise e
                print('\tj=', j, param_list[j], '\t\t Max-Min:', max(y_vs_design_parameter) - min(y_vs_design_parameter))

                data_max[ind].append(max(y_vs_design_parameter))
                data_min[ind].append(min(y_vs_design_parameter))            

            if i==1:
                ax_list[ind].legend(prop={'family':'Times New Roman'})
            ax_list[ind].grid()
            ax_list[ind].set_ylabel(y_label_list[i], **font)

        print('\nObjectives vs. geometry variables:')
        for ind, el in enumerate(data_max):
            print(ind, 'Max', el)
        print('\nObjectives vs. geometry variables:')
        for ind, el in enumerate(data_min):
            print(ind, 'Min', el)

        if self.reference_design is not None:
            print('\n-------------------- Here goes the reference design:')
            for el in self.reference_design[1:]:
                print(el, end=' ')
            self.reference_data = [float(el) for el in self.reference_design[3].split(',')]
            O2_ref = fobj_scalar(self.reference_data[2],
                                 self.reference_data[4],
                                 self.reference_data[3],
                                 self.reference_data[5],
                                 self.reference_data[6],
                                 self.reference_data[INDEX_TOTAL_LOSS],
                                 weights=use_weights('O2'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
            O1_ref = fobj_scalar(self.reference_data[2],
                                 self.reference_data[4],
                                 self.reference_data[3],
                                 self.reference_data[5],
                                 self.reference_data[6],
                                 self.reference_data[INDEX_TOTAL_LOSS],
                                 weights=use_weights('O1'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
        else:
            raise Exception('self.reference_design is None.')

        print('Objective function 1')
        O1 = fobj_list( list(self.get_certain_objective_function(2)), 
                        list(self.get_certain_objective_function(4)), 
                        list(self.get_certain_objective_function(3)), 
                        list(self.get_certain_objective_function(5)), 
                        list(self.get_certain_objective_function(6)), 
                        np.array(list(self.get_certain_objective_function(9))) + np.array(list(self.get_certain_objective_function(12))) + np.array(list(self.get_certain_objective_function(13))),
                        weights=use_weights('O1'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
        O1_max = []
        O1_min = []
        from pylab import figure
        O1_ax  = figure().gca()
        O2_prototype_data = []
        results_for_refining_bounds = {}
        results_for_refining_bounds['O1'] = []
        for j in range(int(len(O1)/number_of_variant)): # iterate design parameters
            O1_vs_design_parameter = O1[j*number_of_variant:(j+1)*number_of_variant]
            O2_prototype_data.append(O1_vs_design_parameter)

            O1_ax.plot(O1_vs_design_parameter, label=str(j)+' '+param_list[j], alpha=0.5)
            print('\t', j, param_list[j], '\t\t max O1 - min O1:', max(O1_vs_design_parameter) - min(O1_vs_design_parameter), '\t\t', end=' ')

            # narrow bounds (refine bounds)
            results_for_refining_bounds['O1'].append( [ind for ind, el in enumerate(O1_vs_design_parameter) if el < O1_ref*1.0] )
            print(results_for_refining_bounds['O1'][j]) #'<- to derive new original_bounds.'

            O1_max.append(max(O1_vs_design_parameter))
            O1_min.append(min(O1_vs_design_parameter))            
        O1_ax.legend()
        O1_ax.grid()
        O1_ax.set_ylabel('O1 [1]', **font)
        O1_ax.set_xlabel('Count of design variants', **font)

        # fig_prototype = figure(500, figsize=(10, 5), facecolor='w', edgecolor='k')
        # O2_prototype_ax = fig_prototype.gca()
        # O2_prototype_ax.plot(list(range(-1, 22)), O1_ref*np.ones(23), 'k--', label='Reference design')
        # O2_prototype_ax.plot(O2_prototype_data[1], 'o-', lw=0.75, alpha=0.5, label=r'$L_g$')
        # O2_prototype_ax.plot(O2_prototype_data[0], 'v-', lw=0.75, alpha=0.5, label=r'$w_{st}$')
        # O2_prototype_ax.plot(O2_prototype_data[3], 's-', lw=0.75, alpha=0.5, label=r'$w_{rt}$')
        # O2_prototype_ax.plot(O2_prototype_data[5], '^-', lw=0.75, alpha=0.5, label=r'$\theta_{so}$')
        # O2_prototype_ax.plot(O2_prototype_data[2], 'd-', lw=0.75, alpha=0.5, label=r'$w_{ro}$')
        # O2_prototype_ax.plot(O2_prototype_data[6], '*-', lw=0.75, alpha=0.5, label=r'$d_{so}$')
        # O2_prototype_ax.plot(O2_prototype_data[4], 'X-', lw=0.75, alpha=0.5, label=r'$d_{ro}$')
        # O2_prototype_ax.legend()
        # O2_prototype_ax.set_ylabel('$O_2(x)$ [1]', **font)

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # O2
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        print('Objective function 2')
        O2 = fobj_list( list(self.get_certain_objective_function(2)), 
                        list(self.get_certain_objective_function(4)), 
                        list(self.get_certain_objective_function(3)), 
                        list(self.get_certain_objective_function(5)), 
                        list(self.get_certain_objective_function(6)), 
                        np.array(list(self.get_certain_objective_function(9))) + np.array(list(self.get_certain_objective_function(12))) + np.array(list(self.get_certain_objective_function(13))),
                        weights=use_weights('O2'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight )
        O2_max = []
        O2_min = []
        O2_ax  = figure().gca()
        O2_ecce_data = []
        results_for_refining_bounds['O2'] = []
        for j in range(int(len(O2)/number_of_variant)): # iterate design parameters: range(7)
            O2_vs_design_parameter = O2[j*number_of_variant:(j+1)*number_of_variant]
            O2_ecce_data.append(O2_vs_design_parameter)

            # narrow bounds (refine bounds)
            O2_ax.plot(O2_vs_design_parameter, 'o-', label=str(j)+' '+param_list[j], alpha=0.5)
            print('\t', j, param_list[j], '\t\t max O2 - min O2:', max(O2_vs_design_parameter) - min(O2_vs_design_parameter), '\t\t', end=' ')
            results_for_refining_bounds['O2'].append( [ind for ind, el in enumerate(O2_vs_design_parameter) if el < O2_ref*1.0] )
            print(results_for_refining_bounds['O2'][j]) #'<- to derive new original_bounds.'

            O2_max.append(max(O2_vs_design_parameter))
            O2_min.append(min(O2_vs_design_parameter))
        O2_ax.legend()
        O2_ax.grid()
        O2_ax.set_ylabel('O2 [1]', **font)
        O2_ax.set_xlabel('Count of design variants', **font)

        # for ecce digest
        fig_ecce = figure(figsize=(10, 5), facecolor='w', edgecolor='k')
        O2_ecce_ax = fig_ecce.gca()
        O2_ecce_ax.plot(list(range(-1, 22)), O2_ref*np.ones(23), 'k--', label='Reference design')
        O2_ecce_ax.plot(O2_ecce_data[1], 'o-', lw=0.75, alpha=0.5,      label=r'$L_g$')
        O2_ecce_ax.plot(O2_ecce_data[0], 'v-', lw=0.75, alpha=0.5,      label=r'$w_{st}$')
        O2_ecce_ax.plot(O2_ecce_data[3], 's-', lw=0.75, alpha=0.5,      label=r'$w_{rt}$')
        O2_ecce_ax.plot(O2_ecce_data[5], '^-', lw=0.75, alpha=0.5,      label=r'$\theta_{so}$')
        O2_ecce_ax.plot(O2_ecce_data[2], 'd-', lw=0.75, alpha=0.5,      label=r'$w_{ro}$')
        O2_ecce_ax.plot(O2_ecce_data[6], '*-', lw=0.75, alpha=0.5,      label=r'$d_{so}$')
        O2_ecce_ax.plot(O2_ecce_data[4], 'X-', lw=0.75, alpha=0.5,      label=r'$d_{ro}$')

        myfontsize = 12.5
        from pylab import plt
        plt.rcParams.update({'font.size': myfontsize})


        # Reference candidate design
        ref = np.zeros(8)
            # ref[0] = 0.635489                                   # PF
            # ref[1] = 0.963698                                   # eta
            # ref[1] = efficiency_at_50kW(1817.22+216.216+224.706)# eta@50kW

        if self.reference_design is not None:
            list_plotting_weights = [8, 3, self.required_torque, 0.1, self.rotor_weight, 0.2, 10, 2500]
            ref[0] = O2_ref                  / list_plotting_weights[0] 
            ref[1] = O1_ref                  / list_plotting_weights[1] 
            ref[2] = self.reference_data[2]  / list_plotting_weights[2]  # 100%
            ref[3] = self.reference_data[3]  / list_plotting_weights[3]  # 100%
            ref[4] = self.reference_data[4]  / list_plotting_weights[4]  # 100% = FRW
            ref[5] = self.reference_data[5]  / list_plotting_weights[5]  # 100%
            ref[6] = self.reference_data[6]  / list_plotting_weights[6]  # deg
            ref[7] = self.reference_data[INDEX_TOTAL_LOSS] / list_plotting_weights[7]  # W

        O1_ax.plot(list(range(-1, 22)), O1_ref*np.ones(23), 'k--')
        O2_ax.plot(list(range(-1, 22)), O2_ref*np.ones(23), 'k--')
        O2_ecce_ax.legend()
        O2_ecce_ax.grid()
        O2_ecce_ax.set_xticks(list(range(21)))
        O2_ecce_ax.annotate('Lower bound', xytext=(0.5, 5.5), xy=(0, 4), xycoords='data', arrowprops=dict(arrowstyle="->"))
        O2_ecce_ax.annotate('Upper bound', xytext=(18.0, 5.5),  xy=(20, 4), xycoords='data', arrowprops=dict(arrowstyle="->"))
        O2_ecce_ax.set_xlim((-0.5,20.5))
        O2_ecce_ax.set_ylim((0,14)) # 4,14
        O2_ecce_ax.set_xlabel(r'Number of design variant', **font)
        O2_ecce_ax.set_ylabel(r'$O_2(x)$ [1]', **font)
        fig_ecce.tight_layout()
        # fig_ecce.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction_full_paper\images\O2_vs_params.png', dpi=150)
        # plt.show()
        # quit() ###################################


        # Maximum
        data_max = np.array(data_max)
        O1_max   = np.array(O1_max)
        O2_max   = np.array(O2_max)
            # data_max[0] = (data_max[0])                   # PF
            # data_max[1] = (data_max[1])                   # eta
            # data_max[1] = efficiency_at_50kW(data_max[7]) # eta@50kW # should use data_min[7] because less loss, higher efficiency
        data_max[0] = O2_max       / list_plotting_weights[0]  
        data_max[1] = O1_max       / list_plotting_weights[1]  
        data_max[2] = (data_max[2])/ list_plotting_weights[2]  # 100%
        data_max[3] = (data_max[3])/ list_plotting_weights[3]  # 100%
        data_max[4] = (data_max[4])/ list_plotting_weights[4]  # 100% = FRW
        data_max[5] = (data_max[5])/ list_plotting_weights[5]  # 100%
        data_max[6] = (data_max[6])/ list_plotting_weights[6]  # deg
        data_max[7] = (data_max[7])/ list_plotting_weights[7]  # W
        y_max_vs_design_parameter_0 = [el[0] for el in data_max]
        y_max_vs_design_parameter_1 = [el[1] for el in data_max]
        y_max_vs_design_parameter_2 = [el[2] for el in data_max]
        y_max_vs_design_parameter_3 = [el[3] for el in data_max]
        y_max_vs_design_parameter_4 = [el[4] for el in data_max]
        y_max_vs_design_parameter_5 = [el[5] for el in data_max]
        y_max_vs_design_parameter_6 = [el[6] for el in data_max]

        # Minimum
        data_min = np.array(data_min)
        O1_min   = np.array(O1_min)
        O2_min   = np.array(O2_min)
            # data_min[0] = (data_min[0])                    # PF
            # data_min[1] = (data_min[1])                    # eta
            # data_min[1] = efficiency_at_50kW(data_min[7])  # eta@50kW
        data_min[0] = O2_min        / list_plotting_weights[0] 
        data_min[1] = O1_min        / list_plotting_weights[1] 
        data_min[2] = (data_min[2]) / list_plotting_weights[2] # 100%
        data_min[3] = (data_min[3]) / list_plotting_weights[3] # 100%
        data_min[4] = (data_min[4]) / list_plotting_weights[4] # 100% = FRW
        data_min[5] = (data_min[5]) / list_plotting_weights[5] # 100%
        data_min[6] = (data_min[6]) / list_plotting_weights[6] # deg
        data_min[7] = (data_min[7]) / list_plotting_weights[7] # W
        y_min_vs_design_parameter_0 = [el[0] for el in data_min]
        y_min_vs_design_parameter_1 = [el[1] for el in data_min]
        y_min_vs_design_parameter_2 = [el[2] for el in data_min]
        y_min_vs_design_parameter_3 = [el[3] for el in data_min]
        y_min_vs_design_parameter_4 = [el[4] for el in data_min]
        y_min_vs_design_parameter_5 = [el[5] for el in data_min]
        y_min_vs_design_parameter_6 = [el[6] for el in data_min]

        count = np.arange(len(y_max_vs_design_parameter_0))  # the x locations for the groups
        width = 1.0  # the width of the bar

        fig = figure(dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = fig.gca()
        # fig, ax = plt.subplots(dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')                                      #  #1034A
        rects1 = ax.bar(count - 3*width/8, y_min_vs_design_parameter_0, width/8, alpha=0.5, label=r'$L_g$, Air gap length', color='#6593F5')
        rects2 = ax.bar(count - 2*width/8, y_min_vs_design_parameter_1, width/8, alpha=0.5, label=r'$b_{st}$, Stator tooth width', color='#1D2951') # https://digitalsynopsis.com/design/beautiful-color-palettes-combinations-schemes/
        rects3 = ax.bar(count - 1*width/8, y_min_vs_design_parameter_2, width/8, alpha=0.5, label=r'$b_{rt}$, Rotor tooth width', color='#03396c')
        rects4 = ax.bar(count - 0*width/8, y_min_vs_design_parameter_3, width/8, alpha=0.5, label=r'$\theta_{so}$, Stator open width', color='#6497b1')
        rects5 = ax.bar(count + 1*width/8, y_min_vs_design_parameter_4, width/8, alpha=0.5, label=r'$w_{ro}$, Rotor open width',  color='#0E4D92')
        rects6 = ax.bar(count + 2*width/8, y_min_vs_design_parameter_5, width/8, alpha=0.5, label=r'$d_{so}$, Stator open depth', color='#005b96')
        rects7 = ax.bar(count + 3*width/8, y_min_vs_design_parameter_6, width/8, alpha=0.5, label=r'$d_{ro}$, Rotor open depth', color='#b3cde0') 
        print('ylim=', ax.get_ylim())
        autolabel(ax, rects1, bias=-0.10, textfont=textfont)
        autolabel(ax, rects2, bias=-0.10, textfont=textfont)
        autolabel(ax, rects3, bias=-0.10, textfont=textfont)
        autolabel(ax, rects4, bias=-0.10, textfont=textfont)
        autolabel(ax, rects5, bias=-0.10, textfont=textfont)
        autolabel(ax, rects6, bias=-0.10, textfont=textfont)
        autolabel(ax, rects7, bias=-0.10, textfont=textfont)
        one_one = np.array([1, 1])
        minus_one_one = np.array([-1, 1])
        ax.plot(rects4[0].get_x() + 0.5*width*minus_one_one, ref[0]*one_one, 'k--', lw=1.0, alpha=0.6, label='Reference design' )
        ax.plot(rects4[1].get_x() + 0.5*width*minus_one_one, ref[1]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[2].get_x() + 0.5*width*minus_one_one, ref[2]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[3].get_x() + 0.5*width*minus_one_one, ref[3]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[4].get_x() + 0.5*width*minus_one_one, ref[4]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[5].get_x() + 0.5*width*minus_one_one, ref[5]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[6].get_x() + 0.5*width*minus_one_one, ref[6]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[7].get_x() + 0.5*width*minus_one_one, ref[7]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.legend(loc='upper right', prop={'family':'Times New Roman'})
        # text for indicating reference values
        ax.text(rects4[0].get_x() - 3.5/8*width, ref[0]*1.01, '%.2f'%(ref[0]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[1].get_x() - 3.5/8*width, ref[1]*1.01, '%.2f'%(ref[1]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[2].get_x() - 3.5/8*width, ref[2]*1.01, '%.2f'%(ref[2]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[3].get_x() - 3.5/8*width, ref[3]*1.01, '%.2f'%(ref[3]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[4].get_x() - 3.5/8*width, ref[4]*1.01, '%.2f'%(ref[4]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[5].get_x() - 3.5/8*width, ref[5]*1.01, '%.2f'%(ref[5]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[6].get_x() - 3.5/8*width, ref[6]*1.01, '%.2f'%(ref[6]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[7].get_x() - 3.5/8*width, ref[7]*1.01, '%.2f'%(ref[7]), ha='center', va='bottom', rotation=90, **textfont)

        rects1 = ax.bar(count - 3*width/8, y_max_vs_design_parameter_0, width/8, alpha=0.5, label=r'$L_g$,         Air gap length', color='#6593F5')    # bottom=y_min_vs_design_parameter_0, 
        rects2 = ax.bar(count - 2*width/8, y_max_vs_design_parameter_1, width/8, alpha=0.5, label=r'$b_{st}$, Stator tooth width', color='#1D2951')     # bottom=y_min_vs_design_parameter_1, 
        rects3 = ax.bar(count - 1*width/8, y_max_vs_design_parameter_2, width/8, alpha=0.5, label=r'$b_{rt}$, Rotor tooth width', color='#03396c')      # bottom=y_min_vs_design_parameter_2, 
        rects4 = ax.bar(count - 0*width/8, y_max_vs_design_parameter_3, width/8, alpha=0.5, label=r'$\theta_{so}$, Stator open width', color='#6497b1') # bottom=y_min_vs_design_parameter_3, 
        rects5 = ax.bar(count + 1*width/8, y_max_vs_design_parameter_4, width/8, alpha=0.5, label=r'$w_{ro}$, Rotor open width',  color='#0E4D92')      # bottom=y_min_vs_design_parameter_4, 
        rects6 = ax.bar(count + 2*width/8, y_max_vs_design_parameter_5, width/8, alpha=0.5, label=r'$d_{so}$, Stator open depth', color='#005b96')      # bottom=y_min_vs_design_parameter_5, 
        rects7 = ax.bar(count + 3*width/8, y_max_vs_design_parameter_6, width/8, alpha=0.5, label=r'$d_{ro}$, Rotor open depth', color='#b3cde0')       # bottom=y_min_vs_design_parameter_6, 
        autolabel(ax, rects1, textfont=textfont)
        autolabel(ax, rects2, textfont=textfont)
        autolabel(ax, rects3, textfont=textfont)
        autolabel(ax, rects4, textfont=textfont)
        autolabel(ax, rects5, textfont=textfont)
        autolabel(ax, rects6, textfont=textfont)
        autolabel(ax, rects7, textfont=textfont)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Normalized Objective Functions', **font)
        ax.set_xticks(count)
        # ax.set_xticklabels(('Power Factor [100%]', r'$\eta$@$T_{em}$ [100%]', r'$T_{em}$ [15.9 N]', r'$T_{rip}$ [10%]', r'$|F|$ [51.2 N]', r'    $E_m$ [20%]', r'      $E_a$ [10 deg]', r'$P_{\rm Cu,Fe}$ [2.5 kW]')))
        # ax.set_xticklabels(('Power Factor [100%]', r'$O_1$ [3]', r'$T_{em}$ [15.9 N]', r'$T_{rip}$ [10%]', r'$|F|$ [51.2 N]', r'    $E_m$ [20%]', r'      $E_a$ [10 deg]', r'$P_{\rm Cu,Fe}$ [2.5 kW]'))
        ax.set_xticklabels(('$O_2$ [%g]'               %(list_plotting_weights[0]), 
                            '$O_1$ [%g]'               %(list_plotting_weights[1]), 
                            '$T_{em}$ [%g Nm]'         %(list_plotting_weights[2]), 
                            '$T_{rip}$ [%g%%]'         %(list_plotting_weights[3]*100), 
                            '$|F|$ [%g N]'             %(list_plotting_weights[4]), 
                            '    $E_m$ [%g%%]'         %(list_plotting_weights[5]*100), 
                            '      $E_a$ [%g deg]'     %(list_plotting_weights[6]), 
                            '$P_{\\rm Cu,Fe}$ [%g kW]' %(list_plotting_weights[7]*1e-3) ), **font)
        ax.grid()
        ax.set_ylim([0,4])
        # fig.tight_layout()
        # fig.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction\images\sensitivity_results.png', dpi=150)

        # plt.show()
        return results_for_refining_bounds

# class FEA_Solver:
#     def __init__(self, fea_config_dict, spec_input_dict, select_fea_config_dict, select_spec):


class acm_designer(object):
    def __init__(self, select_spec, spec_input_dict, select_fea_config_dict, fea_config_dict, acm_template=None):
        # spec and acm_template are objects
        self.acm_template = acm_template
        try:
            self.spec = acm_template.spec
            # print('[acm_designer.py] DEBUG: spec is', self.spec)
        except:
            # print('[acm_designer.py] DEBUG: spec is None.')
            self.spec = None

        self.fea_config_dict        = fea_config_dict
        self.spec_input_dict        = spec_input_dict
        self.select_spec            = select_spec
        self.select_fea_config_dict = select_fea_config_dict

        # to be used with PYGMO
        self.flag_do_not_evaluate_when_init_pop = False
        self.counter_fitness_called = 0
        self.counter_fitness_return = 0

        ###################
        # Initialize FEA solver and results
        ###################
        # self.solver = FEA_Solver(fea_config_dict, spec_input_dict, select_fea_config_dict, select_spec) # The two json inputs are used here
        self.app = None

        # 如果文件夹不存在，那就造起来！
        output_dir = self.fea_config_dict['output_dir']
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        self.dir_csv_output_folder = output_dir + 'csv/'
        if not os.path.isdir(self.dir_csv_output_folder):
            os.makedirs(self.dir_csv_output_folder)
        self.dir_jsonpickle_folder = output_dir + 'jsonpickle/'
        if not os.path.isdir(self.dir_jsonpickle_folder):
            os.makedirs(self.dir_jsonpickle_folder)

        # if fea_config_dict['bool_post_processing'] == False:
        #     self.fig_main, self.axeses = plt.subplots(2, 2, sharex=True, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        #     utility.pyplot_clear(self.axeses)
        # else:
        #     self.fig_main, self.axeses = None, None

        self.folder_to_be_deleted = None

        # if os.path.exists(output_dir+'swarm_MOO_log.txt'):
        #     os.rename(output_dir+'swarm_MOO_log.txt', output_dir+'swarm_MOO_log_backup.txt')
        open(output_dir+'swarm_MOO_log.txt', 'a').close()

        self.acm_template.build_x_denorm()
        # print(self.acm_template.x_denorm_dict)
        swarm_data_file = self.read_swarm_data_json(self.select_spec, self.acm_template.x_denorm_dict)


    def init_logger(self, prefix='pygmo_'):
        # self.logger = utility.myLogger(self.fea_config_dict['output_dir']+'../', prefix=prefix+self.fea_config_dict['run_folder'][:-1])
        self.logger = utility.myLogger(self.fea_config_dict['output_dir']+'../', prefix=prefix)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Automatic Performance Evaluation (This is just a wraper)
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def evaluate_design(self, acm_template, x_denorm, counter, counter_loop=1):
        if 'PM' in acm_template.name:
            return self.fea_wrapper(acm_template, x_denorm, counter, counter_loop)

        elif 'IM' in acm_template.name:
            return self.fea_bearingless_induction(acm_template.obsolete_template, x_denorm, counter, counter_loop)

        elif 'Flux_Alternator' in acm_template.name:
            return self.fea_wrapper(acm_template, x_denorm, counter, counter_loop)

        else:
            raise Exception(f'{acm_template.name} is not regognized as a valid machine type.')

    def save_to_disk(self, acm_variant, spec_performance_dict, GP, EX):

        # wily is not json serilizable, so is recordtype type object: acmop_parameters
        wily = EX['wily']
        EX['wily'] = None
        # EX['wily'] = {
        #     # 'layer_X_phases': wily.layer_X_phases,
        #     # 'layer_X_signs': wily.layer_X_signs,
        #     'coil_pitch_y': wily.coil_pitch_y,
        #     # 'layer_Y_phases': wily.layer_Y_phases,
        #     # 'layer_Y_signs': wily.layer_Y_signs,
        #     #
        #     'grouping_AC': wily.grouping_AC,
        #     'number_parallel_branch': wily.number_parallel_branch,
        #     'number_winding_layer': wily.number_winding_layer,
        #     #
        #     'bool_3PhaseCurrentSource': wily.bool_3PhaseCurrentSource,
        #     'CommutatingSequenceD': wily.CommutatingSequenceD,
        #     'CommutatingSequenceB': wily.CommutatingSequenceB,
        #     #
        #     'deg_winding_U_phase_phase_axis_angle': wily.deg_winding_U_phase_phase_axis_angle,
        #     #
        #     'Qs': wily.Qs,
        #     'p': wily.p,
        #     'ps': wily.ps,
        #     'pr': wily.pr,
        #     'SPP': wily.SPP,
        #     #
        #     'dict_coil_connection': wily.dict_coil_connection,
        # }
        list_of_GP_as_dict = [{key: val._asdict()} for key, val in GP.items()] # see _asdict in https://www.python.org/dev/peps/pep-0557/
        for parameter_key_val_pair in list_of_GP_as_dict:
            # print('DEBUG', parameter_key_val_pair)
            for key, val in parameter_key_val_pair.items():
                val['calc'] = None # function .calc cannot be serialized 

        number_current_generation = spec_performance_dict['number_current_generation'] #= int(acm_variant.counter//popsize), 
        individual_index = spec_performance_dict['individual_index'] #= acm_variant.counter

        big_dict = dict()
        with open(self.fea_config_dict['output_dir'] + self.select_spec + '.json', 'a') as f:
            big_dict[self.select_spec+f'-gen{number_current_generation}-ind{individual_index}'] = {
                'Spec inputs'  :         acm_variant.template.spec_input_dict,
                'x_denorm_dict':         acm_variant.template.get_x_denorm_dict_from_geometric_parameters(GP),
                'Geometric parameters':  list_of_GP_as_dict,
                'Excitations'  :         EX,
                'Performance'  :         spec_performance_dict
                # 'Derived'     :        self.spec.spec_derive_dict,
                # 'Geometry'    :        acm_variant.spec_geometry_dict,
            }
            f.write(f',\n"{acm_variant.counter}":')
            try:
                json.dump(big_dict, f, indent=4)
            except Exception as e:
                print(f'[acm_designer.py] [Warning] You might need to clean up .json data file yourself at {self.fea_config_dict["output_dir"]}\n'*3)
                raise e

        EX['wily'] = wily

    def evaluate_design_json_wrapper(self, acm_template, x_denorm=None, counter='Test', counter_loop=1):
        # This is a wrapper for the wrapper, in order to build up a json profile for the design variant

        # 这里应该返回新获得的设计，然后可以获得geometry_dict，然后包括x_denorm的信息方便重构设计。
        acm_variant = self.evaluate_design(acm_template, x_denorm, counter, counter_loop)

        if 'FEMM' in self.select_fea_config_dict:
            acm_variant.results_for_optimization = acm_variant.analyzer.build_results_for_optimization()

            # Save spec_performance_dict and others to disk
            GP = acm_variant.template.d['GP']
            EX = acm_variant.template.d['EX']
            self.save_to_disk(acm_variant, acm_variant.analyzer.spec_performance_dict, GP, EX)

            # Save also the object (acm_variant) to disk, but this takes a lot of disk space!
            if self.fea_config_dict['moo.save_acm_variant_object_as_jsonpickle'] == True:
                utility_json.to_json_recursively(acm_variant, acm_variant.name, save_here=self.fea_config_dict['output_dir']+'jsonpickle/')

            # Save time domain data to disk
            acm_variant.analyzer.save_time_domain_data(self.fea_config_dict['output_dir']+self.select_spec+f'-ind{counter}.pkl') # counter could be string

            return acm_variant

        elif 'JMAG' in self.select_fea_config_dict:

            cost_function, f1, f2, f3, FRW, \
            normalized_torque_ripple, \
            normalized_force_error_magnitude, \
            force_error_angle, \
            project_name, individual_name, \
            number_current_generation, individual_index,\
            power_factor, \
            rated_ratio, \
            rated_stack_length_mm, \
            rated_total_loss, \
            rated_stator_copper_loss_along_stack, \
            rated_magnet_Joule_loss, \
            rated_rotor_copper_loss_along_stack, \
            stator_copper_loss_in_end_turn, \
            rotor_copper_loss_in_end_turn, \
            rated_iron_loss, \
            rated_windage_loss, \
            str_results, \
            mm2_slot_area, \
            coil_flux_linkage_peak2peak_value, \
            TRV, Cost, Cost_Fe, Cost_Cu, Cost_PM, \
            ss_avg_force_magnitude, rotor_weight, torque_average = acm_variant.results_to_be_unpacked

            # acm_variant.spec_geometry_dict['x_denorm'] = list(x_denorm)

            spec_performance_dict = dict()
            spec_performance_dict['project_name'] = project_name
            spec_performance_dict['individual_name'] = individual_name
            spec_performance_dict['number_current_generation'] = number_current_generation
            spec_performance_dict['individual_index'] = individual_index
            # spec_performance_dict['cost_function'] = cost_function
            spec_performance_dict['f1'] = f1
            spec_performance_dict['f2'] = f2
            spec_performance_dict['f3'] = f3
            spec_performance_dict['TRV'] = TRV
            spec_performance_dict['FRW'] = FRW
            spec_performance_dict['torque_average'] = torque_average
            spec_performance_dict['ss_avg_force_magnitude'] = ss_avg_force_magnitude
            spec_performance_dict['rotor_weight'] = rotor_weight
            spec_performance_dict['normalized_torque_ripple'] = normalized_torque_ripple
            spec_performance_dict['normalized_force_error_magnitude'] = normalized_force_error_magnitude
            spec_performance_dict['force_error_angle'] = force_error_angle
            spec_performance_dict['coil_flux_linkage_peak2peak_value'] = coil_flux_linkage_peak2peak_value
            spec_performance_dict['mm2_slot_area'] = mm2_slot_area
            spec_performance_dict['Cost'] = Cost
            spec_performance_dict['Cost_Fe'] = Cost_Fe
            spec_performance_dict['Cost_Cu'] = Cost_Cu
            spec_performance_dict['Cost_PM'] = Cost_PM
            spec_performance_dict['power_factor'] = power_factor
            spec_performance_dict['rated_ratio'] = rated_ratio
            spec_performance_dict['rated_stack_length_mm'] = rated_stack_length_mm
            spec_performance_dict['rated_total_loss'] = rated_total_loss
            spec_performance_dict['rated_stator_copper_loss_along_stack'] = rated_stator_copper_loss_along_stack
            spec_performance_dict['rated_rotor_copper_loss_along_stack'] = rated_rotor_copper_loss_along_stack
            spec_performance_dict['rated_magnet_Joule_loss'] = rated_magnet_Joule_loss
            spec_performance_dict['stator_copper_loss_in_end_turn'] = stator_copper_loss_in_end_turn
            spec_performance_dict['rotor_copper_loss_in_end_turn'] = rotor_copper_loss_in_end_turn
            spec_performance_dict['rated_iron_loss'] = rated_iron_loss
            spec_performance_dict['rated_windage_loss'] = rated_windage_loss
            # spec_performance_dict['str_results'] = str_results
            spec_performance_dict['select_fea_config_dict'] = self.select_fea_config_dict
            spec_performance_dict['moo.fitness_OA'] = self.fea_config_dict['moo.fitness_OA']
            spec_performance_dict['moo.fitness_OB'] = self.fea_config_dict['moo.fitness_OB']
            spec_performance_dict['moo.fitness_OC'] = self.fea_config_dict['moo.fitness_OC']

            GP = acm_variant.template.d['GP']
            EX = acm_variant.template.d['EX']

            # Save to disk
            self.save_to_disk(acm_variant, spec_performance_dict, GP, EX)

            # save object (acm_variant) to disk
            utility_json.to_json_recursively(acm_variant, acm_variant.name, save_here=self.fea_config_dict['output_dir']+'jsonpickle/')

            # this is for optimization
            acm_variant.results_for_optimization = (cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle)

            # this is for comparison to FEMM
            def compare_with_FEMM(acm_variant):
                EX = acm_variant.template.d['EX']
                acm_variant.analyzer = FEMM_SlidingMesh.Individual_Analyzer_FEMM_Edition(p=EX['wily'].p)
                basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list, \
                    DisplacementAngle_list, \
                    circuit_current_GroupACU, \
                    circuit_current_GroupACV, \
                    circuit_current_GroupACW, \
                    circuit_current_GroupBDU, \
                    circuit_current_GroupBDV, \
                    circuit_current_GroupBDW, \
                    terminal_voltage_GroupACU, \
                    terminal_voltage_GroupACV, \
                    terminal_voltage_GroupACW, \
                    terminal_voltage_GroupBDU, \
                    terminal_voltage_GroupBDV, \
                    terminal_voltage_GroupBDW, \
                    coil_fluxLinkage_GroupACU, \
                    coil_fluxLinkage_GroupACV, \
                    coil_fluxLinkage_GroupACW, \
                    coil_fluxLinkage_GroupBDU, \
                    coil_fluxLinkage_GroupBDV, \
                    coil_fluxLinkage_GroupBDW = self.toolJd.dm.unpack(bool_more_info=True)
                electrical_period = acm_variant.template.fea_config_dict['designer.number_cycles_in_2ndTSS']/EX['DriveW_Freq']
                number_of_steps   = acm_variant.template.fea_config_dict['designer.number_of_steps_2ndTSS']
                step_size_sec = electrical_period / number_of_steps
                step_size_mech_deg = EX['Omega'] * step_size_sec / np.pi * 180

                for index in range(-self.toolJd.dm.number_of_steps_at_steady_state, 0):
                    time                         = float(time_list[index])
                    RotorAngle_MechanicalDegrees = float(DisplacementAngle_list[index])
                    torque = float(TorCon_list[index])
                    forces = ( float(ForConX_list[index]), float(ForConY_list[index]) )
                    energy = 0.0
                    circuitProperties = ( [ circuit_current_GroupACU[index], terminal_voltage_GroupACU[index], coil_fluxLinkage_GroupACU[index] ],
                                        [ circuit_current_GroupACV[index], terminal_voltage_GroupACV[index], coil_fluxLinkage_GroupACV[index] ],
                                        [ circuit_current_GroupACW[index], terminal_voltage_GroupACW[index], coil_fluxLinkage_GroupACW[index] ],
                                        [ circuit_current_GroupBDU[index], terminal_voltage_GroupBDU[index], coil_fluxLinkage_GroupBDU[index] ],
                                        [ circuit_current_GroupBDV[index], terminal_voltage_GroupBDV[index], coil_fluxLinkage_GroupBDV[index] ],
                                        [ circuit_current_GroupBDW[index], terminal_voltage_GroupBDW[index], coil_fluxLinkage_GroupBDW[index] ] )
                    acm_variant.analyzer.add(time, RotorAngle_MechanicalDegrees, torque, forces, energy, circuitProperties)
                acm_variant.analyzer.get_ss_data()
            compare_with_FEMM(acm_variant)
            acm_variant.analyzer.save_time_domain_data(counter) # TODO
            return acm_variant


    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # FEA solver and FEA results parser
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def read_swarm_survivor(self, popsize):
        if not os.path.exists(self.fea_config_dict['output_dir'] + 'swarm_survivor.txt'):
            return None

        with open(self.fea_config_dict['output_dir'] + 'swarm_survivor.txt', 'r') as f:
            buf = f.readlines()
            survivor_data_raw = buf[-popsize:]
            survivor_data = [[float(s) for s in line.split(',')] for line in survivor_data_raw]

        self.survivor_title_number = float(buf[-popsize-1][len('---------'):])

        # for el in survivor_data:
        #     print('\t', el)
        # quit()
        return survivor_data

    def write_swarm_survivor(self, pop, counter_fitness_return):
        with open(self.fea_config_dict['output_dir'] + 'swarm_survivor.txt', 'a') as f:
            f.write('\n---------%d\n'%(counter_fitness_return) \
                    + '\n'.join(','.join('%.16f'%(x) for x in el[0].tolist() + el[1].tolist() ) for el in zip(pop.get_x(), pop.get_f()) )) # convert 2d array to string

    def read_swarm_data_old(self, bound_filter=None, read_from_here=None):
        if read_from_here is not None:
            self.fea_config_dict['output_dir'] = read_from_here
        print('\t', self.fea_config_dict['output_dir'] + 'swarm_data.txt')
        if not os.path.exists(self.fea_config_dict['output_dir'] + 'swarm_data.txt'):
            msg = '\tNo file @ ' + self.fea_config_dict['output_dir'] + 'swarm_data.txt'
            print(msg)
            # raise Exception(msg)
            return None

        print('\t\tRead in', self.fea_config_dict['output_dir'] + 'swarm_data.txt')
        with open(self.fea_config_dict['output_dir'] + 'swarm_data.txt', 'r') as f:
            buf = f.readlines()
            buf = buf[1:]
            length_buf = len(buf) 

            if length_buf % 21 == 0:
                pass
            else:
                raise Exception('Invalid swarm_data.txt!')

            number_of_chromosome = length_buf / 21
            if number_of_chromosome == 0:
                return None

            self.swarm_data_raw = [buf[i:i+21] for i in range(0, len(buf), 21)]
            self.swarm_data_container = swarm_data_container(self.swarm_data_raw, self.fea_config_dict, bound_filter)
            self.swarm_data = self.swarm_data_container.swarm_data_xf
            # for el in self.swarm_data:
            #     print(el)
            # quit()
            return int(number_of_chromosome)

    def read_swarm_data(self, select_spec, read_from_here=None):
        if read_from_here is not None:
            self.fea_config_dict['output_dir'] = read_from_here
        print('[acm_designer.py]', self.fea_config_dict['output_dir'] + 'swarm_data.txt')
        if not os.path.exists(self.fea_config_dict['output_dir'] + 'swarm_data.txt'):
            msg = '\tNo file @ ' + self.fea_config_dict['output_dir'] + 'swarm_data.txt'
            print(msg)
            # raise Exception(msg)
            return None

        print('[acm_designer.py]', '\tRead in', self.fea_config_dict['output_dir'] + 'swarm_data.txt')
        with open(self.fea_config_dict['output_dir'] + 'swarm_data.txt', 'r') as f:
            buf = f.readlines()
            buf = buf[1:]
            length_buf = len(buf) 

            if length_buf % 21 == 0:
                pass
            else:
                raise Exception('Invalid swarm_data.txt!')

            number_of_chromosome = length_buf / 21
            if number_of_chromosome == 0:
                return None

            self.swarm_data_raw = [buf[i:i+21] for i in range(0, len(buf), 21)]

        print('[acm_designer.py]', self.fea_config_dict['output_dir'], select_spec)

        if os.path.exists(self.fea_config_dict['output_dir'] + select_spec+'.json'):
            with open(self.fea_config_dict['output_dir'] + select_spec+'.json', 'r') as f:
                # skip the first line
                for _ in range(1): next(f)
                buf = f.read()
                self.swarm_data_json = json.loads('{\n' + buf + '\n}')
                print('[acm_designer.py] ACMOP: read in swarm_data_json...', self.fea_config_dict['output_dir'] + select_spec+'.json')
                # print(buf) # debug
        else:
            self.swarm_data_json = None
            print('[acm_designer.py] These are bopt-python results. No swarm_data_json available.')

        self.swarm_data_container = swarm_data_container(self.swarm_data_raw, self.fea_config_dict, self.swarm_data_json)
        self.swarm_data = self.swarm_data_container.swarm_data_xf
        # for el in self.swarm_data:
        #     print(el)
        # quit()
        return int(number_of_chromosome)

    def read_swarm_data_json(self, select_spec, desired_x_denorm_dict=None):
        ''' In case of desired_x_denorm_dict being None, your swarm_data_xf must be of the same size.
        '''
        if select_spec is None: select_spec = self.select_spec
        self.analyzer = Swarm_Data_Analyzer(self.fea_config_dict['output_dir'] + select_spec + '.json', desired_x_denorm_dict)
        self.swarm_data = self.analyzer.swarm_data_xf
        self.swarm_data_file = self.fea_config_dict['output_dir'] + select_spec + '.json'
        return self.swarm_data_file

    def build_acm_variant(self, template, x_denorm, counter, counter_loop=1):
        if 'SPMSM' in template.machine_type:
            acm_variant = bearingless_spmsm_design.bearingless_spmsm_design_variant(template=template, x_denorm=x_denorm, counter=counter, counter_loop=counter_loop)
        elif 'PMVM' in template.machine_type:
            acm_variant = vernier_motor_design.vernier_motor_VShapePM_design_variant(template=template, x_denorm=x_denorm, counter=counter, counter_loop=counter_loop)
        elif 'Flux_Alternator' in template.machine_type:
            acm_variant = flux_alternator_design.flux_alternator_design_variant(template=template, x_denorm=x_denorm, counter=counter, counter_loop=counter_loop)
        elif 'FSPM' in template.machine_type:
            acm_variant = flux_switching_pm_design.FSPM_design_variant(template=template, x_denorm=x_denorm, counter=counter, counter_loop=counter_loop)
        else:
            raise Exception('Not supported machine_type:', template.machine_type)
        return acm_variant

    def fea_wrapper(self, template, x_denorm, counter, counter_loop, bool_re_evaluate=False):
        logger = logging.getLogger(__name__)
        msg = f'Run FEA for individual #{counter}'
        logger.info(msg)
        # print(msg)

        # get local design variant
        acm_variant = self.build_acm_variant(template, x_denorm, counter, counter_loop)

        if 'JMAG' in self.select_fea_config_dict:
            # project name
            self.project_name = acm_variant.name
            self.expected_project_file = self.fea_config_dict['output_dir'] + "%s.jproj"%(self.project_name)

            # study name
            study_name = acm_variant.name + "-Transient" # Change here and there 

            # project meta data
            project_meta_data = {
                "expected_project_file": self.expected_project_file,
                "project_name": self.project_name,
                "study_name": study_name,
                "dir_csv_output_folder": self.dir_csv_output_folder,
                "output_dir": self.fea_config_dict['output_dir']
            }

            # Leave the solving task to JMAG
            self.toolJd = self.build_jmag_project(acm_variant, project_meta_data, bool_re_evaluate=bool_re_evaluate)

            ################################################################
            # Load data for cost function evaluation
            ################################################################
            acm_variant.results_to_be_unpacked = results_to_be_unpacked = self.toolJd.build_str_results(acm_variant, self.project_name, study_name, self.dir_csv_output_folder, self.fea_config_dict, femm_solver=None)
            if results_to_be_unpacked is not None:
                if self.toolJd.fig_main is not None:
                    try:
                        if False:
                            self.toolJd.fig_main.savefig(self.fea_config_dict['output_dir'] + acm_variant.name + 'results.png', dpi=150)
                    except Exception as e:
                        print(e)
                        print('\n\n\nIgnore error and continue.')
                    finally:
                        utility.pyplot_clear(self.toolJd.axeses)
                # show()
                return acm_variant 
            else:
                raise Exception('[acm_designer] results_to_be_unpacked is None.')

        elif 'FEMM' in self.select_fea_config_dict:
            self.toolFEMM = self.build_femm_project(acm_variant)
            # acm_variant.results_to_be_unpacked = results_to_be_unpacked = toolFEMM.build_str_results(self.axeses, acm_variant, self.project_name, study_name, self.dir_csv_output_folder, self.fea_config_dict, femm_solver=None)
            return acm_variant
        else:
            raise Exception('[acm_designer.py] Wrong string of select_fea_config_dict:', self.select_fea_config_dict)

    ''' Produce JMAG Project
    '''
    def build_jmag_project(self, acm_variant, project_meta_data=None, bool_re_evaluate=False):
        if project_meta_data is None:
            # this object is reloaded, so retrieve saved meta data
            project_meta_data = self.project_meta_data
        else:
            # save for reload this object from json
            self.project_meta_data = project_meta_data
        expected_project_file = project_meta_data["expected_project_file"]
        study_name            = project_meta_data["study_name"]
        dir_csv_output_folder = project_meta_data["dir_csv_output_folder"]
        output_dir            = project_meta_data["output_dir"]

        # use alias
        # acm_variant = self

        # Leave the solving task to JMAG
        if bool_re_evaluate==False:

            # def draw_jmag_bpmsm():
            toolJd = JMAG.JMAG(self.fea_config_dict, self.spec_input_dict)

            toolJd.open(expected_project_file)
            if 'PMSM' in acm_variant.template.name:
                DRAW_SUCCESS = toolJd.draw_spmsm(acm_variant)
            elif 'Flux_Alternator' in acm_variant.template.name:
                DRAW_SUCCESS = toolJd.draw_doublySalient(acm_variant, bool_draw_whole_model=True)
            elif 'FSPM' in acm_variant.template.name:
                DRAW_SUCCESS = toolJd.draw_FSPM(acm_variant, bool_draw_whole_model=False)
            elif 'CPPM' in acm_variant.template.name:
                DRAW_SUCCESS = toolJd.draw_CPPM(acm_variant)
            else:
                raise Exception ('Add a new motor type')
            if DRAW_SUCCESS != 1:
                raise Exception('Drawer failed.')

            app = toolJd.app

            # JMAG
            if app.NumModels()>=1:
                model = app.GetModel(acm_variant.name)
            else:
                logger = logging.getLogger(__name__)
                logger.error('there is no model yet for %s'%(acm_variant.name))
                raise Exception('why is there no model yet? %s'%(acm_variant.name))

            if 'PMSM' in acm_variant.template.name:
                toolJd.pre_process_PMSM(app, model, acm_variant)
            elif 'FSPM' in acm_variant.template.name:
                toolJd.pre_process_FSPM(app, model, acm_variant)
            elif 'Flux_Alternator' in acm_variant.template.name:
                toolJd.pre_process_fluxAlternator(app, model, acm_variant)
            elif 'CPPM' in acm_variant.template.name:
                toolJd.pre_process_CPPM(app, model, acm_variant)

            study = toolJd.add_magnetic_transient_study(app, model, dir_csv_output_folder, study_name, acm_variant)
            toolJd.mesh_study(acm_variant, app, model, study, output_dir=output_dir)
            # raise KeyboardInterrupt
            toolJd.run_study(acm_variant, app, study, acm_variant.template.fea_config_dict, clock_time())


            # export Voltage if field data exists.
            if acm_variant.template.fea_config_dict['delete_results_after_calculation'] == False:
                # Export Circuit Voltage
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(dir_csv_output_folder + study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
        else:
            ''' This is a special debugging feature that allows to recover FEA results (swarm_data.txt) from JMAG Deisnger's csv results.
            '''
            THE_mm2_magnet_area = acm_variant.rotorMagnet.draw(None, bool_re_evaluate=True)
            THE_mm2_slot_area = acm_variant.coils.draw(None, bool_re_evaluate=True)

            CurrentAmp_in_the_slot = THE_mm2_slot_area * acm_variant.fill_factor * acm_variant.Js*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / acm_variant.template.d['EX']['DriveW_zQ']
            CurrentAmp_per_phase = CurrentAmp_per_conductor * acm_variant.template.d['EX']['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            acm_variant.CurrentAmp_per_phase = CurrentAmp_per_phase
            acm_variant.DriveW_CurrentAmp = acm_variant.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            acm_variant.BeariW_CurrentAmp = acm_variant.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            slot_area_utilizing_ratio = (acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp) / acm_variant.CurrentAmp_per_phase
            print('[acm_designer.py]---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)

            print('[acm_designer.py]---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
            print('[acm_designer.py]---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
            print('[acm_designer.py]---acm_variant.DriveW_CurrentAmp =', acm_variant.DriveW_CurrentAmp)
            print('[acm_designer.py]---acm_variant.BeariW_CurrentAmp =', acm_variant.BeariW_CurrentAmp)
            print('[acm_designer.py]---TORQUE_CURRENT_RATIO:', acm_variant.fea_config_dict['TORQUE_CURRENT_RATIO'])
            print('[acm_designer.py]---SUSPENSION_CURRENT_RATIO:', acm_variant.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

        return toolJd

    ''' Create FEMM Project
    '''
    def build_femm_project(self, acm_variant):
        # Leave the solving task to FEMM

        # 1. Drawer
        toolFEMM = FEMM_SlidingMesh.FEMM_SlidingMesh(acm_variant)
        toolFEMM.open()
        toolFEMM.probdef(stack_length=acm_variant.template.d['EX']['mm_template_stack_length'], time_harmonic_study_frequency=0)
        if 'PMSM' in acm_variant.template.name:
            if toolFEMM.draw_spmsm(acm_variant) != 1: raise Exception('Drawer failed.')
        elif 'Flux_Alternator' in acm_variant.template.name:
            if toolFEMM.draw_doublySalient(acm_variant) != 1: raise Exception('Drawer failed.')
        elif 'FSPM' in acm_variant.template.name:
            if toolFEMM.draw_doublySalient(acm_variant) != 1: raise Exception('Drawer failed.')
        else:
            raise

        # 2. Preprocess
        toolFEMM.pre_process(project_file_name=acm_variant.template.fea_config_dict['output_dir'] + acm_variant.name + '.fem')

        # 3. Run transient FEA and collect results
        toolFEMM.run_transient_study()

        # 4. Compute objective
        acm_variant.analyzer.get_ss_data()
        acm_variant.analyzer._compute_objectives(acm_variant, self.select_fea_config_dict, toolFEMM) # generate spec_performance_dict

        # 5. Save results to dist (this is done in evaluate_json_wrapper)

        self.acm_variant = acm_variant
        return toolFEMM

    ''' BELOW ARE BOPT-PYTHON CODES for BLIM ONLY
    '''
    def fea_bearingless_induction(self, im_template, x_denorm, counter, counter_loop):
        logger = logging.getLogger(__name__)
        print('Run FEA for individual #%d'%(counter))

        # get local design variant
        im_variant = population.bearingless_induction_motor_design.local_design_variant(im_template, 0, counter, x_denorm)

        # print('::', im_template.Radius_OuterRotor, im_template.Width_RotorSlotOpen)
        # print('::', im_variant.Radius_OuterRotor, im_variant.Width_RotorSlotOpen)
        # quit()

        # TODO: Change indivudal name to be more useful
        # if counter_loop == 1:
        #     im_variant.name = 'ind%d'%(counter)
        # else:
        #     im_variant.name = 'ind%d-redo%d'%(counter, counter_loop)
        # im_variant.spec = im_template.spec
        self.im_variant = im_variant
        self.femm_solver = FEMM_Solver.FEMM_Solver(self.im_variant, flag_read_from_jmag=False, freq=50) # eddy+static
        im = None



        if counter_loop == 1:
            self.project_name          = 'proj%d'%(counter)
        else:
            self.project_name          = 'proj%d-redo%d'%(counter, counter_loop)
        self.expected_project_file = self.fea_config_dict['output_dir'] + "%s.jproj"%(self.project_name)

        original_study_name = im_variant.name + "Freq"
        tran2tss_study_name = im_variant.name + 'Tran2TSS'

        self.dir_femm_temp         = self.fea_config_dict['output_dir'] + 'femm_temp/'
        self.femm_output_file_path = self.dir_femm_temp + original_study_name + '.csv'

        # self.jmag_control_state = False

        # local scripts
        def open_jmag(expected_project_file_path):
            if self.app is None:
                # app = win32com.client.Dispatch('designer.Application.181')
                app = win32com.client.Dispatch('designer.Application.171')
                # app = win32com.client.gencache.EnsureDispatch('designer.Application.171') # https://stackoverflow.com/questions/50127959/win32-dispatch-vs-win32-gencache-in-python-what-are-the-pros-and-cons

                if self.fea_config_dict['designer.Show'] == True:
                    app.Show()
                else:
                    app.Hide()
                # app.Quit()
                self.app = app # means that the JMAG Designer is turned ON now.
                self.bool_run_in_JMAG_Script_Editor = False

                def add_steel(self):
                    print('[First run on this computer detected]', im_template.spec_input_dict['Steel'], 'is added to jmag material library.')
                    import population
                    if 'M15' in im_template.spec_input_dict['Steel']:
                        population.add_M1xSteel(self.app, self.fea_config_dict['dir.parent'], steel_name="M-15 Steel")
                    elif 'M19' in im_template.spec_input_dict['Steel']:
                        population.add_M1xSteel(self.app, self.fea_config_dict['dir.parent'])
                    elif 'Arnon5' == im_template.spec_input_dict['Steel']:
                        population.add_Arnon5(self.app, self.fea_config_dict['dir.parent'])

                # too avoid tons of the same material in JAMG's material library
                fname = self.fea_config_dict['dir.parent'] + '.jmag_state.txt'
                if not os.path.exists(fname):
                    with open(fname, 'w') as f:
                        f.write(self.fea_config_dict['pc_name'] + '/' + im_template.spec_input_dict['Steel'] + '\n')
                    add_steel(self)
                else:
                    with open(fname, 'r') as f:
                        for line in f.readlines():
                            if self.fea_config_dict['pc_name'] + '/' + im_template.spec_input_dict['Steel'] not in line:
                                add_steel(self)

            else:
                app = self.app

            print('[acm_designer.py]', expected_project_file_path)
            if os.path.exists(expected_project_file_path):
                print('[acm_designer.py] JMAG project exists already. I learned my lessions. I will NOT delete it but create a new one with a different name instead.')
                # os.remove(expected_project_file_path)
                attempts = 2
                temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)
                while os.path.exists(temp_path):
                    attempts += 1
                    temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)

                expected_project_file_path = temp_path

            app.NewProject("Untitled")
            app.SaveAs(expected_project_file_path)
            logger.debug('Create JMAG project file: %s'%(expected_project_file_path))

            return app

        def draw_jmag(app):
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Draw the model in JMAG Designer
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            DRAW_SUCCESS = self.draw_jmag_induction(app,
                                                counter, 
                                                im_variant,
                                                im_variant.name)
            if DRAW_SUCCESS == 0:
                # TODO: skip this model and its evaluation
                cost_function = 99999 # penalty
                logging.getLogger(__name__).warn('Draw Failed for %s-%s\nCost function penalty = %g.%s', self.project_name, im_variant.name, cost_function, self.im_variant.show(toString=True))
                raise Exception('Draw Failed: Are you working on the PC? Sometime you by mistake operate in the JMAG Geometry Editor, then it fails to draw.')
                return None
            elif DRAW_SUCCESS == -1:
                raise Exception(' DRAW_SUCCESS == -1:')

            # JMAG
            if app.NumModels()>=1:
                model = app.GetModel(im_variant.name)
            else:
                logger.error('there is no model yet for %s'%(im_variant.name))
                raise Exception('why is there no model yet? %s'%(im_variant.name))
            return model

        def rotating_static_FEA():

            # wait for femm to finish, and get your slip of breakdown
            new_fname = self.dir_femm_temp + original_study_name + '.csv'
            with open(new_fname, 'r') as f:
                data = f.readlines()
                freq = float(data[0][:-1])
                torque = float(data[1][:-1])
            slip_freq_breakdown_torque, breakdown_torque, breakdown_force = freq, torque, None
            # Now we have the slip, set it up!
            im_variant.update_mechanical_parameters(slip_freq_breakdown_torque) # do this for records only

            # Must run this script after slip_freq_breakdown_torque is known to get new results, but if results are present, it is okay to use this script for post-processing.
            if not self.femm_solver.has_results():
                print('run_rotating_static_FEA')
                # utility.blockPrint()
                self.femm_solver.run_rotating_static_FEA()
                self.femm_solver.parallel_solve()
                # utility.enablePrint()

            # collecting parasolve with post-process
            # wait for .ans files
            # data_femm_solver = self.femm_solver.show_results_static(bool_plot=False) # this will wait as well?
            while not self.femm_solver.has_results():
                print(clock_time())
                sleep(3)
            results_dict = {}
            for f in [f for f in os.listdir(self.femm_solver.dir_run) if 'static_results' in f]:
                data = np.loadtxt(self.femm_solver.dir_run + f, unpack=True, usecols=(0,1,2,3))
                for i in range(len(data[0])):
                    results_dict[data[0][i]] = (data[1][i], data[2][i], data[3][i]) 
            keys_without_duplicates = [key for key, item in results_dict.items()]
            keys_without_duplicates.sort()
            with open(self.femm_solver.dir_run + "no_duplicates.txt", 'w') as fw:
                for key in keys_without_duplicates:
                    fw.writelines('%g %g %g %g\n' % (key, results_dict[key][0], results_dict[key][1], results_dict[key][2]))
            data_femm_solver = np.array([ keys_without_duplicates, 
                                             [results_dict[key][0] for key in keys_without_duplicates], 
                                             [results_dict[key][1] for key in keys_without_duplicates], 
                                             [results_dict[key][2] for key in keys_without_duplicates]])
            # print(data_femm_solver)
            # from pylab import plt
            # plt.figure(); plt.plot(data_femm_solver[0])
            # plt.figure(); plt.plot(data_femm_solver[1])
            # plt.figure(); plt.plot(data_femm_solver[2])
            # plt.figure(); plt.plot(data_femm_solver[3])
            # plt.show()
            # quit()
            return data_femm_solver

        def rotating_eddy_current_FEA(im, app, model):

            # Freq Sweeping for break-down Torque Slip
            # remember to export the B data using subroutine 
            # and check export table results only
            study = im.add_study(app, model, self.dir_csv_output_folder, choose_study_type='frequency')

            # Freq Study: you can choose to not use JMAG to find the breakdown slip.
            # Option 1: you can set im.slip_freq_breakdown_torque by FEMM Solver
            # Option 2: Use JMAG to sweeping the frequency
            # Does study has results?
            if study.AnyCaseHasResult():
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.check_csv_results(study.GetName())
            else:
                # mesh
                im.add_mesh(study, model)

                # Export Image
                    # for i in range(app.NumModels()):
                    #     app.SetCurrentModel(i)
                    #     model = app.GetCurrentModel()
                    #     app.ExportImage(r'D:\Users\horyc\OneDrive - UW-Madison\pop\run#10/' + model.GetName() + '.png')
                app.View().ShowAllAirRegions()
                # app.View().ShowMeshGeometry() # 2nd btn
                app.View().ShowMesh() # 3rn btn
                app.View().Zoom(3)
                app.View().Pan(-im.Radius_OuterRotor, 0)
                app.ExportImageWithSize('./' + model.GetName() + '.png', 2000, 2000)
                app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted because only ouput table results are selected.

                # run
                study.RunAllCases()
                app.Save()

                def check_csv_results(dir_csv_output_folder, study_name, returnBoolean=False, file_suffix='_torque.csv'): # '_iron_loss_loss.csv'
                    # print self.dir_csv_output_folder + study_name + '_torque.csv'
                    if not os.path.exists(dir_csv_output_folder + study_name + file_suffix):
                        if returnBoolean == False:
                            print('Nothing is found when looking into:', dir_csv_output_folder + study_name + file_suffix)
                            return None
                        else:
                            return False
                    else:
                        if returnBoolean == True:
                            return True

                    try:
                        # check csv results 
                        l_slip_freq = []
                        l_TorCon    = []
                        l_ForCon_X  = []
                        l_ForCon_Y  = []

                        with open(dir_csv_output_folder + study_name + '_torque.csv', 'r') as f: 
                            for ind, row in enumerate(utility.csv_row_reader(f)):
                                if ind >= 5:
                                    try:
                                        float(row[0])
                                    except:
                                        continue
                                    l_slip_freq.append(float(row[0]))
                                    l_TorCon.append(float(row[1]))

                        with open(dir_csv_output_folder + study_name + '_force.csv', 'r') as f: 
                            for ind, row in enumerate(utility.csv_row_reader(f)):
                                if ind >= 5:
                                    try:
                                        float(row[0])
                                    except:
                                        continue
                                    # l_slip_freq.append(float(row[0]))
                                    l_ForCon_X.append(float(row[1]))
                                    l_ForCon_Y.append(float(row[2]))

                        breakdown_force = max(np.sqrt(np.array(l_ForCon_X)**2 + np.array(l_ForCon_Y)**2))

                        index, breakdown_torque = utility.get_index_and_max(l_TorCon)
                        slip_freq_breakdown_torque = l_slip_freq[index]
                        return slip_freq_breakdown_torque, breakdown_torque, breakdown_force
                    except NameError as e:
                        logger = logging.getLogger(__name__)
                        logger.error('No CSV File Found.', exc_info=True)
                        raise e

                # evaluation based on the csv results
                print(':::ZZZ', self.dir_csv_output_folder)
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = check_csv_results(self.dir_csv_output_folder, study.GetName())

            # this will be used for other duplicated studies
            original_study_name = study.GetName()
            im.csv_previous_solve = self.dir_csv_output_folder + original_study_name + '_circuit_current.csv'
            im.update_mechanical_parameters(slip_freq_breakdown_torque, syn_freq=im.DriveW_Freq)


            # EC Rotate: Rotate the rotor to find the ripples in force and torque # 不关掉这些云图，跑第二个study的时候，JMAG就挂了：app.View().SetVectorView(False); app.View().SetFluxLineView(False); app.View().SetContourView(False)
            ecrot_study_name = original_study_name + "-FFVRC"
            casearray = [0 for i in range(1)]
            casearray[0] = 1
            model.DuplicateStudyWithCases(original_study_name, ecrot_study_name, casearray)

            app.SetCurrentStudy(ecrot_study_name)
            study = app.GetCurrentStudy()

            divisions_per_slot_pitch = 24 # self.fea_config_dict['ec_rotate_divisions_per_slot_pitch']  # 24
            study.GetStep().SetValue("Step", divisions_per_slot_pitch) 
            study.GetStep().SetValue("StepType", 0)
            study.GetStep().SetValue("FrequencyStep", 0)
            study.GetStep().SetValue("Initialfrequency", slip_freq_breakdown_torque)

                # study.GetCondition(u"RotCon").SetValue(u"MotionGroupType", 1)
            study.GetCondition("RotCon").SetValue("Displacement", + 360.0/im.Qr/divisions_per_slot_pitch)

            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

            # model.RestoreCadLink()
            study.Run()
            app.Save()
            # model.CloseCadLink()

        def transient_FEA_as_reference(im, slip_freq_breakdown_torque):
            # Transient Reference
            tranRef_study_name = "TranRef"
            if model.NumStudies()<4:
                model.DuplicateStudyWithType(tran2tss_study_name, "Transient2D", tranRef_study_name)
                app.SetCurrentStudy(tranRef_study_name)
                study = app.GetCurrentStudy()

                # 将一个滑差周期和十个同步周期，分成 400 * end_point / (1.0/im.DriveW_Freq) 份。
                end_point = 0.5/slip_freq_breakdown_torque + 10.0/im.DriveW_Freq
                # Pavel Ponomarev 推荐每个电周期400~600个点来捕捉槽效应。
                division = self.fea_config_dict['designer.TranRef-StepPerCycle'] * end_point / (1.0/im.DriveW_Freq)  # int(end_point * 1e4)
                                                                        # end_point = division * 1e-4
                study.GetStep().SetValue("Step", division + 1) 
                study.GetStep().SetValue("StepType", 1) # regular inverval
                study.GetStep().SetValue("StepDivision", division)
                study.GetStep().SetValue("EndPoint", end_point)

                # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
                study.GetStudyProperties().SetValue("DirectSolverType", 1)

                study.RunAllCases()
                app.Save()

        # debug for tia-iemdc-ecce-2019
        # data_femm_solver = rotating_static_FEA()
        # from show_results_iemdc19 import show_results_iemdc19
        # show_results_iemdc19(   self.dir_csv_output_folder, 
        #                         im_variant, 
        #                         femm_solver_data=data_femm_solver, 
        #                         femm_rotor_current_function=self.femm_solver.get_rotor_current_function()
        #                     )
        # quit()

        if self.fea_config_dict['bool_re_evaluate']==False:

            # this should be summoned even before initializing femm, and it will decide whether the femm results are reliable
            app = open_jmag(self.expected_project_file) # will set self.jmag_control_state to True

            ################################################################
            # Begin from where left: Frequency Study
            ################################################################
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Eddy Current Solver for Breakdown Torque and Slip
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # check for existing results
            if os.path.exists(self.femm_output_file_path):
                # for file in os.listdir(self.dir_femm_temp):
                #     if original_study_name in file:
                #         print('----------', original_study_name, file)
                print('Remove legacy femm output files @ %s'%(self.femm_output_file_path))
                os.remove(self.femm_output_file_path)
                os.remove(self.femm_output_file_path[:-4]+'.fem')
                # quit()
                # quit()

            # delete me
            # model = draw_jmag(app)

            # At this point, no results exist from femm.
            print('Run greedy_search_for_breakdown_slip...')
            femm_tic = clock_time()
            # self.femm_solver.__init__(im_variant, flag_read_from_jmag=False, freq=50.0)
            if im_variant.DriveW_poles == 4 and self.fea_config_dict['femm.use_fraction'] == True:
                print('FEMM model only solves for a fraction of 2.\n'*3)
                self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                    bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=2)
            else:
                # p >= 3 is not tested so do not use fraction for now
                self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                    bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=1) # 转子导条必须形成通路

            ################################################################
            # Begin from where left: Transient Study
            ################################################################
            model = draw_jmag(app)

            # EC-Rotate
            # rotating_eddy_current_FEA(im_variant, app, model)

            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # TranFEAwi2TSS for ripples and iron loss
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # add or duplicate study for transient FEA denpending on jmag_run_list
            # FEMM+JMAG (注意，这里我们用50Hz作为滑差频率先设置起来，等拿到breakdown slip freq的时候，再更新变量slip和study properties的时间。)
            study = im_variant.add_TranFEAwi2TSS_study( 50.0, app, model, self.dir_csv_output_folder, tran2tss_study_name, logger)
            app.SetCurrentStudy(tran2tss_study_name)
            study = app.GetCurrentStudy()
            self.mesh_study(im_variant, app, model, study)

            # wait for femm to finish, and get your slip of breakdown
            slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.femm_solver.wait_greedy_search(femm_tic)

            # Now we have the slip, set it up!
            im_variant.update_mechanical_parameters(slip_freq_breakdown_torque) # do this for records only
            if im_variant.the_slip != slip_freq_breakdown_torque / im_variant.DriveW_Freq:
                raise Exception('Check update_mechanical_parameters().')
            study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im_variant.the_slip))
            if True:
                number_of_steps_2ndTSS = self.fea_config_dict['designer.number_of_steps_2ndTSS'] 
                DM = app.GetDataManager()
                DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
                refarray = [[0 for i in range(3)] for j in range(3)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = 0.5/slip_freq_breakdown_torque #0.5 for 17.1.03l # 1 for 17.1.02y
                refarray[1][1] =    number_of_steps_2ndTSS                          # 16 for 17.1.03l #32 for 17.1.02y
                refarray[1][2] =        50
                refarray[2][0] = refarray[1][0] + 0.5/im_variant.DriveW_Freq #0.5 for 17.1.03l 
                refarray[2][1] =    number_of_steps_2ndTSS  # also modify range_ss! # don't forget to modify below!
                refarray[2][2] =        50
                DM.GetDataSet("SectionStepTable").SetTable(refarray)
                number_of_total_steps = 1 + 2 * number_of_steps_2ndTSS # [Double Check] don't forget to modify here!
                study.GetStep().SetValue("Step", number_of_total_steps)
                study.GetStep().SetValue("StepType", 3)
                study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

            # static FEA solver with FEMM (need eddy current FEA results)
            # print('::', im_variant.Omega, im_variant.Omega/2/np.pi)
            # print('::', femm_solver.im.Omega, femm_solver.im.Omega/2/np.pi)
            # quit()
            # rotating_static_FEA()

            # debug JMAG circuit
            # app.Save()
            # quit()

            # Run JMAG study
            self.run_study(im_variant, app, study, clock_time())

            # export Voltage if field data exists.
            if self.fea_config_dict['delete_results_after_calculation'] == False:
                # Export Circuit Voltage
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                # app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(self.dir_csv_output_folder + im_variant.name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
                print('[acm_designer.py] WriteTable to:', self.dir_csv_output_folder + tran2tss_study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
                print('[acm_designer.py] WriteTable to (converter):', os.path.abspath(self.dir_csv_output_folder + tran2tss_study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv"))
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(os.path.abspath(self.dir_csv_output_folder + tran2tss_study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")) # must be absolute path to JMAG

            # TranRef
            # transient_FEA_as_reference(im_variant, slip_freq_breakdown_torque)

        else:
            # FEMM的转子电流，哎，是个麻烦事儿。。。

            self.femm_solver.vals_results_rotor_current = []

            new_fname = self.dir_femm_temp + original_study_name + '.csv'
            try:
                with open(new_fname, 'r') as f:
                    buf = f.readlines()
                    for idx, el in enumerate(buf):
                        if idx == 0:
                            slip_freq_breakdown_torque = float(el)
                            continue
                        if idx == 1:
                            breakdown_torque = float(el)
                            continue
                        if idx == 2:
                            self.femm_solver.stator_slot_area = float(el)
                            continue
                        if idx == 3:
                            self.femm_solver.rotor_slot_area = float(el)
                            continue
                        # print(el)
                        # print(el.split(','))
                        temp = el.split(',')
                        self.femm_solver.vals_results_rotor_current.append( float(temp[0])+ 1j*float(temp[1]) )
                        # print(self.femm_solver.vals_results_rotor_current)
                self.dirty_backup_stator_slot_area = self.femm_solver.stator_slot_area
                self.dirty_backup_rotor_slot_area = self.femm_solver.rotor_slot_area
                self.dirty_backup_vals_results_rotor_current = self.femm_solver.vals_results_rotor_current
            except FileNotFoundError as error:
                print(error)
                print('Use dirty_backup to continue...') # 有些时候，不知道为什么femm的结果文件（.csv）没了，这时候曲线救国，凑活一下吧
                self.femm_solver.stator_slot_area = self.dirty_backup_stator_slot_area           
                self.femm_solver.rotor_slot_area = self.dirty_backup_rotor_slot_area            
                self.femm_solver.vals_results_rotor_current = self.dirty_backup_vals_results_rotor_current 




            # 电机的电流值取决于槽的面积。。。。
            THE_mm2_slot_area = self.femm_solver.stator_slot_area*1e6

            if 'VariableStatorSlotDepth' in self.fea_config_dict['which_filter']:
                # set DriveW_CurrentAmp using the calculated stator slot area.
                print('[A]: DriveW_CurrentAmp is updated.')

                # 槽深变化，电密不变，所以电流也会变化。
                CurrentAmp_in_the_slot = THE_mm2_slot_area * im_variant.fill_factor * im_variant.Js*1e-6 * np.sqrt(2)
                CurrentAmp_per_conductor = CurrentAmp_in_the_slot / im_variant.DriveW_zQ
                CurrentAmp_per_phase = CurrentAmp_per_conductor * im_variant.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor

                im_variant.CurrentAmp_per_phase = CurrentAmp_per_phase

                im_variant.DriveW_CurrentAmp = self.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                im_variant.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

                slot_area_utilizing_ratio = (im_variant.DriveW_CurrentAmp + im_variant.BeariW_CurrentAmp) / im_variant.CurrentAmp_per_phase
                print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)
                
                print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                print('---im_variant.DriveW_CurrentAmp =', im_variant.DriveW_CurrentAmp)
                print('---im_variant.BeariW_CurrentAmp =', im_variant.BeariW_CurrentAmp)
                print('---TORQUE_CURRENT_RATIO:', self.fea_config_dict['TORQUE_CURRENT_RATIO'])
                print('---SUSPENSION_CURRENT_RATIO:', self.fea_config_dict['SUSPENSION_CURRENT_RATIO'])





        ################################################################
        # Load data for cost function evaluation
        ################################################################
        im_variant.results_to_be_unpacked = results_to_be_unpacked = utility.build_str_results(self.axeses, im_variant, self.project_name, tran2tss_study_name, self.dir_csv_output_folder, self.fea_config_dict, self.femm_solver)
        if results_to_be_unpacked is not None:
            if self.fig_main is not None:
                try:
                    # self.fig_main.savefig(self.fea_config_dict['output_dir'] + im_variant.name + 'results.png', dpi=150)
                    pass
                except Exception as e:
                    print('Directory exists?', self.fea_config_dict['output_dir'] + im_variant.name + 'results.png')
                    print(e)
                    print('\n\n\nIgnore error and continue.')
                finally:
                    utility.pyplot_clear(self.axeses)
            # show()
            return im_variant
        else:
            raise Exception('results_to_be_unpacked is None.')
        # winding analysis? 之前的python代码利用起来啊
        # 希望的效果是：设定好一个设计，马上进行运行求解，把我要看的数据都以latex报告的形式呈现出来。
        # OP_PS_Qr36_M19Gauge29_DPNV_NoEndRing.jproj

    def draw_jmag_induction(self, app, individual_index, im_variant, model_name, bool_trimDrawer_or_vanGogh=True, doNotRotateCopy=False):

        if individual_index == -1: # 后处理是-1
            print('Draw model for post-processing')
            if individual_index+1 + 1 <= app.NumModels():
                logger = logging.getLogger(__name__)
                logger.debug('The model already exists for individual with index=%d. Skip it.' % individual_index)
                return -1 # the model is already drawn

        elif individual_index+1 <= app.NumModels(): # 一般是从零起步
            logger = logging.getLogger(__name__)
            logger.debug('The model already exists for individual with index=%d. Skip it.' % individual_index)
            return -1 # the model is already drawn

        # open JMAG Geometry Editor
        app.LaunchGeometryEditor()
        geomApp = app.CreateGeometryEditor()
        # geomApp.Show()
        geomApp.NewDocument()
        doc = geomApp.GetDocument()
        ass = doc.GetAssembly()

        # draw parts
        try:
            if bool_trimDrawer_or_vanGogh:
                d = population.TrimDrawer(im_variant) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.plot_shaft("Shaft")

                d.plot_rotorCore("Rotor Core")
                d.plot_cage("Cage")

                d.plot_statorCore("Stator Core")

                d.plot_coil("Coil")
                # d.plot_airWithinRotorSlots(u"Air Within Rotor Slots")

                if 'VariableStatorSlotDepth' in self.fea_config_dict['which_filter']:
                    # set DriveW_CurrentAmp using the calculated stator slot area.
                    print('[A]: DriveW_CurrentAmp is updated.')

                    # 槽深变化，电密不变，所以电流也会变化。
                    CurrentAmp_in_the_slot = d.mm2_slot_area * im_variant.fill_factor * im_variant.Js*1e-6 * np.sqrt(2)
                    CurrentAmp_per_conductor = CurrentAmp_in_the_slot / im_variant.DriveW_zQ
                    CurrentAmp_per_phase = CurrentAmp_per_conductor * im_variant.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                    variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor

                    im_variant.CurrentAmp_per_phase = CurrentAmp_per_phase

                    im_variant.DriveW_CurrentAmp = self.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                    im_variant.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

                    slot_area_utilizing_ratio = (im_variant.DriveW_CurrentAmp + im_variant.BeariW_CurrentAmp) / im_variant.CurrentAmp_per_phase
                    print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)

                    print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                    print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                    print('---im_variant.DriveW_CurrentAmp =', im_variant.DriveW_CurrentAmp)
                    print('---im_variant.BeariW_CurrentAmp =', im_variant.BeariW_CurrentAmp)
                    print('---TORQUE_CURRENT_RATIO:', self.fea_config_dict['TORQUE_CURRENT_RATIO'])
                    print('---SUSPENSION_CURRENT_RATIO:', self.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            else:
                d = VanGogh_JMAG(im_variant, doNotRotateCopy=doNotRotateCopy) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.draw_model()
            self.d = d
        except Exception as e:
            print('See log file to plotting error.')
            logger = logging.getLogger(__name__)
            logger.error('The drawing is terminated. Please check whether the specified bounds are proper.', exc_info=True)

            raise e

            # print 'Draw Failed'
            # if self.pc_name == 'Y730':
            #     # and send the email to hory chen
            #     raise e

            # or you can skip this model and continue the optimization!
            return False # indicating the model cannot be drawn with the script.

        # Import Model into Designer
        doc.SaveModel(True) # True=on : Project is also saved. 
        model = app.GetCurrentModel() # model = app.GetModel(u"IM_DEMO_1")
        model.SetName(model_name)
        model.SetDescription(im_variant.model_name_prefix + '\n' + im_variant.show(toString=True))

        if doNotRotateCopy:
            im_variant.pre_process_structural(app, d.listKeyPoints)
        else:
            im_variant.pre_process(app)

        model.CloseCadLink() # this is essential if you want to create a series of models
        return True

    # def run_study(self, im_variant, app, study, toc):
    #     logger = logging.getLogger(__name__)
    #     if self.fea_config_dict['designer.JMAG_Scheduler'] == False:
    #         print('[acm_designer.py] Run jam.exe...')
    #         # if run_list[1] == True:
    #         try:
    #             study.RunAllCases()
    #         except Exception as error:
    #             raise error
    #         msg = '[acm_designer.py] Time spent on %s is %g s.'%(study.GetName() , clock_time() - toc)
    #         logger.debug(msg)
    #         print(msg)
    #     else:
    #         print('[acm_designer.py] Submit to JMAG_Scheduler...')
    #         job = study.CreateJob()
    #         job.SetValue("Title", study.GetName())
    #         job.SetValue("Queued", True)
    #         job.Submit(False) # Fallse:CurrentCase, True:AllCases
    #         logger.debug('Submit %s to queue (Tran2TSS).'%(im_variant.individual_name))
    #         # wait and check
    #         # study.CheckForCaseResults()
    #     app.Save()

    #         # if the jcf file already exists, it pops a msg window
    #         # study.WriteAllSolidJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Solid', True) # True : Outputs cases that do not have results 
    #         # study.WriteAllMeshJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Mesh', True)

    #         # # run
    #         # if self.fea_config_dict['JMAG_Scheduler'] == False:
    #         #     study.RunAllCases()
    #         #     app.Save()
    #         # else:
    #         #     job = study.CreateJob()
    #         #     job.SetValue(u"Title", study.GetName())
    #         #     job.SetValue(u"Queued", True)
    #         #     job.Submit(True)
    #         #     logger.debug('Submit %s to queue (Freq).'%(im_variant.individual_name))
    #         #     # wait and check
    #         #     # study.CheckForCaseResults()

    # def mesh_study(self, im_variant, app, model, study):

    #     # this `if' judgment is effective only if JMAG-DeleteResultFiles is False 
    #     # if not study.AnyCaseHasResult(): 
    #     # mesh
    #     im_variant.add_mesh(study, model)

    #     # Export Image
    #     app.View().ShowAllAirRegions()
    #     # app.View().ShowMeshGeometry() # 2nd btn
    #     app.View().ShowMesh() # 3rn btn
    #     app.View().Zoom(3)
    #     app.View().Pan(-im_variant.Radius_OuterRotor, 0)
    #     app.ExportImageWithSize(self.fea_config_dict['output_dir'] + model.GetName() + '.png', 2000, 2000)
    #     app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.
  


    ''' BELOW IS OBSOLETE
    '''

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Post-processing
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def check_results_of_local_sensitivity_analysis(self):
        if self.fea_config_dict['local_sensitivity_analysis'] == True:
            run_folder = self.fea_config_dict['run_folder'][:-1] + 'lsa/'
        else:
            run_folder = self.fea_config_dict['run_folder']

        # Sensitivity Bar Charts
        return os.path.exists(fea_config_dict['dir.parent'] + 'pop/' + run_folder + 'swarm_data.txt')

    def collect_results_of_local_sensitivity_analysis(self):
        print('Start to collect results for local sensitivity analysis...')
        try:
            results_for_refining_bounds = utility.build_sensitivity_bar_charts(self.spec, self.sw)
            # quit()
        except Exception as e:
            raise e
            os.remove(fea_config_dict['dir.parent'] + 'pop/' + run_folder + 'swarm_data.txt')
            print('[acm_designer.py] Remove ' + fea_config_dict['dir.parent'] + 'pop/' + run_folder + 'swarm_data.txt')
            print('[acm_designer.py] Continue for sensitivity analysis...')

        self.results_for_refining_bounds = results_for_refining_bounds
        # return results_for_refining_bounds

    def build_refined_bounds(self, the_bounds):

        de_config_dict = self.de_config_dict
        number_of_variants= self.fea_config_dict['local_sensitivity_analysis_number_of_variants'] # 故意不加1的

        print('-'*20, 'results_for_refining_bounds (a.k.a. results_for_refining_bounds):')
        str_to_file = 'When you are done, replace this line with "Done", save file and close all figures.\n'
        for key, val in self.results_for_refining_bounds.items():
            print('---', key)
            for el in val:
                print('\t', el)
            str_to_file += key + '\n' + '\n'.join( ','.join(f'{x}' for x in y) for y in val ) + '\n'

        self.logger.debug('The default refining factors will be for: %s'%(self.fea_config_dict['use_weights']))

        with open('./refining_factors.txt', 'w') as f:
            f.write(str_to_file)
        file_backup_user_input = './refining_factors_%s.txt'%(self.fea_config_dict['run_folder'][:-1])
        if os.path.exists(file_backup_user_input):
            os.system('start %s'%(file_backup_user_input))
        os.system('start ./refining_factors.txt')
        from pylab import show
        show()
        from time import sleep
        while True:
            with open('./refining_factors.txt') as f:
                buf = f.read()
                if buf[:4] == 'Done' or buf[:4] == 'done':
                    print('Done')
                    break
                else:
                    print('.', end='')
                    sleep(1)
        if os.path.exists(file_backup_user_input):
            os.remove(file_backup_user_input)
        os.rename('./refining_factors.txt', file_backup_user_input)

        buf_list = buf.split('\n')
        # for el in buf_list:
        #     print(el)

        print('fea_config_dict - use_weights changed from %s to %s'%(self.fea_config_dict['use_weights'], buf_list[1]))
        self.fea_config_dict['use_weights'] = buf_list[1]

        user_input_for_refining_bounds = [[float(x) for x in y.split(',')] for y in buf_list[2:2+7]]
        print('user_input_for_refining_bounds')
        for el in user_input_for_refining_bounds:
            print('\t', el)

        for ind, bound in enumerate(user_input_for_refining_bounds):
            下界 = bound[0]
            if 下界 != 0:
                下界 -= 1
            上界 = bound[-1]
            if 上界 != number_of_variants:
                上界 += 1
            self.de_config_dict['narrow_bounds_normalized'][ind].append(下界/number_of_variants)
            self.de_config_dict['narrow_bounds_normalized'][ind].append(上界/number_of_variants)

        self.de_config_dict['bounds'] = []
        for bnd1, bnd2 in zip(self.de_config_dict['original_bounds'], self.de_config_dict['narrow_bounds_normalized']):
            diff = bnd1[1] - bnd1[0]
            self.de_config_dict['bounds'].append( [ bnd1[0]+diff*bnd2[0] , bnd1[0]+diff*bnd2[1] ]) # 注意，都是乘以original_bounds的上限哦！

        print('-'*40)
        print('narrow_bounds_normalized:')
        for el in self.de_config_dict['narrow_bounds_normalized']:
            print('\t', el)

        print('original_bounds:')
        for el in self.de_config_dict['original_bounds']:
            print('\t', el)

        print('refined bounds:')
        for el in self.de_config_dict['bounds']:
            print('\t', el)

        return de_config_dict['bounds']

    def build_local_bounds_from_best_design(self, best_design):
        raise Exception('build_local_bounds_from_best_design')
        return local_bounds

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Automatic Report Generation
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def build_oneReport(self):
        if '730' in self.fea_config_dict['pc_name']:
            os.system('cd /d "'+ self.fea_config_dict['dir.parent'] + 'release/OneReport/OneReport_TEX" && z_nul"')

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Talk to Database
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def talk_to_mysql_database(self):
        if self.spec.bool_bad_specifications:
            print('\nThe specifiaction can not be fulfilled. Read script log or OneReport.pdf for information and revise the specifiaction for $J_r$ or else your design name is wrong.')
        else:
            print('\nThe specifiaction is meet. Now check the database of blimuw if on Y730.')
            if '730' in self.fea_config_dict['pc_name']:
                utility.communicate_database(self.spec)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 2. Initilize Swarm and Initial Pyrhonen's Design (Run this part in JMAG) and femm solver (if required by run_list)
    #    Bounds: 1e-1也还是太小了（第三次报错），至少0.5mm长吧 
    #    # 1e-1 is the least geometry value. 
    #    A 1e-2 will leads to：转子闭口槽极限，会导致edge过小，从而报错：small arc entity exists.png
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def init_swarm(self):
        self.sw = population.swarm(self.fea_config_dict, de_config_dict=self.de_config_dict)
        # sw.show(which='all')
        # print sw.im.show(toString=True)
        # quit()

        if self.fea_config_dict['jmag_run_list'][0] == 0:
            self.sw.init_femm_solver() # define app.sw.femm_solver        

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 3. Run DE Optimization
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def run_de(self):
        sw = self.sw
        logger = self.logger

        count_abort = 0
        logger.debug('-------------------------count_abort=%d' % (count_abort))

        # if optimization_flat == True:
        # generate the initial generation
        sw.generate_pop(specified_initial_design_denorm=None)
        # quit()

        # add initial_design of Pyrhonen09 to the initial generation
        if sw.fea_config_dict['local_sensitivity_analysis'] == False:
            if count_abort == 0:
                utility.add_Pyrhonen_design_to_first_generation(sw, self.de_config_dict, logger)

        # write FEA config to disk
        sw.write_to_file_fea_config_dict()

        if True:
            try:
                de_generator = sw.de()
                # run
                # result = list(de_generator)
                for result in de_generator:
                    print(result)
            except Exception as e:
                print('See log file for the error msg.')
                logger.error('Optimization aborted.', exc_info=True)

                raise e
                # 避免死循环
                count_abort+1
                if count_abort > 10:
                    quit()
            else:
                logger.info('Done.')
                utility.send_notification('Done.')

    def get_de_config(self):

        self.de_config_dict = { 'original_bounds': self.get_original_bounds(),
                                'mut':        0.8,
                                'crossp':     0.7,
                                'popsize':    35, # 5~10 \times number of geometry parameters --JAC223
                                'iterations': 70,
                                'narrow_bounds_normalized':[[],
                                                            [],
                                                            [],
                                                            [],
                                                            [],
                                                            [],
                                                            [] ], # != []*7 （完全是两回事）
                                'bounds':None}
        return self.de_config_dict

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 4. Post-processing
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def best_design_by_weights(self, use_weights):
        if use_weights != self.fea_config_dict['use_weights']:
            raise Exception('Not implemented')

        # import utility
        swda = utility.build_Pareto_plot(self.spec, self.sw)
        # sw.fobj(999, individual_denorm)

        # Final LaTeX Report
        print('According to Pyrhonen09, 300 MPa is the typical yield stress of iron core.')
        initial_design = utility.Pyrhonen_design(self.sw.im, self.de_config_dict['bounds'])
        print (initial_design.design_parameters_denorm)
        print (swda.best_design_denorm)
        print (self.de_config_dict['bounds'])

        best_report_dir_prefix = '../release/OneReport/BestReport_TEX/contents/'
        file_name = 'final_report'
        file_suffix = '.tex'
        fname = open(best_report_dir_prefix+file_name+'_s01'+file_suffix, 'w', encoding='utf-8')
        def combine_lists_alternating_2(list1, list2):
            if abs(len(list1) - len(list2))<=1:
                result = [None]*(len(list1)+len(list2))
                result[::2] = list1
                result[1::2] = list2
                return result
            else:
                raise Exception('Try this (not tested).') # https://stackoverflow.com/questions/3678869/pythonic-way-to-combine-two-lists-in-an-alternating-fashion
                import itertools 
                return [x for x in itertools.chain.from_iterable(itertools.izip_longest(list1,list2)) if x]
        def combine_lists_alternating_4(list1, list2, list3, list4):
            result = [None]*(len(list1)+len(list2)+len(list3)+len(list4))
            result[::4] = list1
            result[1::4] = list2
            result[2::4] = list3
            result[3::4] = list4
            return result
        lower_bounds = [el[0] for el in self.de_config_dict['bounds']]
        upper_bounds = [el[1] for el in self.de_config_dict['bounds']]
        design_data = combine_lists_alternating_4(  initial_design.design_parameters_denorm, 
                                                    swda.best_design_denorm,
                                                    lower_bounds,
                                                    upper_bounds,)
        latex_table = r'''
            \begin{table*}[!t]
              \caption{Comparison of the key geometry parameters between best design and initial design}
              \centering
                \begin{tabular}{ccccc}
                    \hline
                    \hline
                    \thead{Geometry parameters} &
                    \thead{Initial\\\relax design} &
                    \thead{The best\\\relax design} &
                    \thead{Lower\\\relax bounds} &
                    \thead{Upper\\\relax bounds} \\
                    \hline
                    Air gap length $L_g$                    & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Stator tooth width $w_{st}$             & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Rotor tooth width $w_{rt}$              & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Stator slot open angle $\theta_{so}$    & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Rotor slot open width $w_{ro}$          & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Stator slot open depth $d_{so}$         & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Rotor slot open depth $d_{ro}$          & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    \hline
                    \vspace{-2.5ex}
                    \\
                    \multicolumn{3}{l}{*Note: Blah.}
                \end{tabular}
              \label{tab:001}
              \vspace{-3ex}
            \end{table*}
            ''' % tuple(design_data)
        print('''\\subsection{Best Design of Objective Function %s}\n\n
                    ''' % (self.fea_config_dict['use_weights']) + latex_table, file=fname)
        print(swda.str_best_design_details, file=fname)
        fname.close()

        os.system('cd /d '+ r'"D:\OneDrive - UW-Madison\c\release\OneReport\BestReport_TEX" && z_nul"') # 必须先关闭文件！否则编译不起来的
        # import subprocess
        # subprocess.call(r"D:\OneDrive - UW-Madison\c\release\OneReport\OneReport_TEX\z_nul", shell=True)

    def run_local_sensitivity_analysis(self, the_bounds, design_denorm=None):
        # if design_denorm not in the_bounds: then raise 
        if design_denorm is not None:
            for ind, el in enumerate(design_denorm):
                if el < the_bounds[ind][0] or el > the_bounds[ind][1]:
                    raise Exception('给的设计不在边界的内部')

        self.logger.debug('---------\nBegin Local Sensitivity Analysis')

        # de_config_dict['bounds'] 还没有被赋值
        self.de_config_dict['bounds'] = the_bounds

        self.init_swarm() # define app.sw
        self.sw.generate_pop(specified_initial_design_denorm=design_denorm)
        de_generator = self.sw.de()
        for result in de_generator:
            print(result)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 5. Check mechanical strength for the best design
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def run_static_structural_fea(self, design_denorm):
        sw = self.sw
        im_best = population.bearingless_induction_motor_design.local_design_variant(sw.im, \
                     -1, -1, design_denorm)

        # initialize JMAG Designer
        sw.designer_init()
        sw.app.Show()
        project_name = fea_config_dict['run_folder'][:-1]+'_Best' # spec.build_name() is too long...
        expected_project_file = sw.dir_project_files + "%s.jproj"%(project_name)
        print(expected_project_file)
        if not os.path.exists(expected_project_file):
            sw.app.NewProject("Untitled")
            sw.app.SaveAs(expected_project_file)
            logger.debug('Create JMAG project file: %s'%(expected_project_file))
        else:
            sw.app.Load(expected_project_file)
            logger.debug('Load JMAG project file: %s'%(expected_project_file))
            logger.debug('Existing models of %d are found in %s' % (sw.app.NumModels(), sw.app.GetDefaultModelFolderPath()))

        # draw the model in JMAG Designer
        DRAW_SUCCESS = sw.draw_jmag_induction( -1, 
                                             im_best,
                                             'Best %s %s'%(fea_config_dict['use_weights'], fea_config_dict['run_folder'][:-1]),
                                             bool_trimDrawer_or_vanGogh=False,
                                             doNotRotateCopy=True)
        if DRAW_SUCCESS == 0:
            raise Exception('Drawing failed')
        elif DRAW_SUCCESS == -1:
            print('Model Already Exists')
        # quit()

        model = sw.app.GetCurrentModel()
        if model.NumStudies() == 0:
            expected_csv_output_dir = sw.dir_csv_output_folder+'structural/'
            if not os.path.isdir(expected_csv_output_dir):
                os.makedirs(expected_csv_output_dir)
            study = im_best.add_structural_study(sw.app, model, expected_csv_output_dir) # 文件夹名应该与jproj同名
        else:
            study = model.GetStudy(0)


        if study.AnyCaseHasResult():
            pass
        else:
            # Add cases 
            study.GetDesignTable().AddParameterVariableName(u"Centrifugal_Force (CentrifugalForce2D): AngularVelocity")
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(1, 0, 45000) # r/min
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(2, 0, 15000) # r/min
            study.GetDesignTable().AddParameterVariableName(u"MeshSizeControl (ElementSizeOnPart): Size")
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(3, 1, 0.5) # mm
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(4, 1, 0.05) # mm

            # run (mesh is included in add_study_structural)
            study.RunAllCases()
            sw.app.Save()

            # Results-Graphs-Calculations-Add Part Calculation
            study.CreateCalculationDefinition(u"VonMisesStress")
            study.GetCalculationDefinition(u"VonMisesStress").SetResultType(u"MisesStress", u"")
            study.GetCalculationDefinition(u"VonMisesStress").SetResultCoordinate(u"Global Rectangular")
            study.GetCalculationDefinition(u"VonMisesStress").SetCalculationType(u"max")
            study.GetCalculationDefinition(u"VonMisesStress").ClearParts()
            study.GetCalculationDefinition(u"VonMisesStress").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)

            # show graph, Tab File|Edit|Calculation, click on Calculation - Response Graph Data to register response data
            parameter = sw.app.CreateResponseDataParameter(u"VonMisesStress")
            parameter.SetCalculationType(u"SingleValue")
            parameter.SetStartValue(u"-1")
            parameter.SetEndValue(u"-1")
            parameter.SetUnit(u"s")
            parameter.SetVariable(u"VonMisesStress")
            parameter.SetAllLine(False)
            parameter.SetCaseRangeType(1)
            parameter.SetLine(u"Maximum Value")
            sw.app.GetDataManager().CreateParametricDataWithParameter(study.GetDataSet(u"VonMisesStress", 4), parameter)

            # Contour results
            study.CreateScaling(u"Scale100")
            study.GetScaling(u"Scale100").SetScalingFactor(100)
            # sw.app.View().SetOriginalModelView(True)
            sw.app.View().ShowMeshGeometry()
            sw.app.View().SetScaledDisplacementView(True)
            study.CreateContour(u"MisesElement")
            study.GetContour(u"MisesElement").SetResultType(u"MisesStress", u"")
            study.GetContour(u"MisesElement").SetResultCoordinate(u"Global Rectangular")
            study.GetContour(u"MisesElement").SetContourType(2)
            study.GetContour(u"MisesElement").SetDigitsNotationType(2)

            study.GetContour(u"MisesElement").SetLogScale(True)
            study.GetContour(u"MisesElement").SetNumLabels(u"11")
            study.GetContour(u"MisesElement").SetPrecision(u"1")
            study.GetContour(u"MisesElement").SetGradient(u"PurpleRed", u"11", False)

            sw.app.View().SetContourView(True)



            # study.CreateContour(u"PrincipleStressElement")
            # study.GetContour(u"PrincipleStressElement").SetResultType(u"PrincipalStress", u"")
            # study.GetContour(u"PrincipleStressElement").SetComponent(u"I")
            # study.GetContour(u"PrincipleStressElement").SetContourType(2)
            # study.GetContour(u"PrincipleStressElement").SetDigitsNotationType(2)

            sw.app.Save()


        from pylab import show
        show()

def get_bad_fintess_values(machine_type='IM', ref=False):
    if ref == False:
        if 'IM' in machine_type:
            return 0, 0, 99
        elif 'PMSM' in machine_type:
            return 9999, 0, 999
    else:
        if 'IM' in machine_type:
            return 1,     10, 100
        elif 'PMSM' in machine_type:
            return 10000, 10, 1000





