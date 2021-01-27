# https://scipy-lectures.org/packages/3d_plotting/index.html#figure-management

# coding:u8
import shutil
from utility import my_execfile
from utility_moo import *
from win32com.client import pywintypes
bool_post_processing = True # solve or post-processing
bool_re_evaluate = False # re-evaluate the designs using csv (without calling FEA softwares)

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# 0. FEA Setting / General Information & Packages Loading
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
my_execfile('./default_setting.py', g=globals(), l=locals())
fea_config_dict
fea_config_dict['local_sensitivity_analysis'] = False
fea_config_dict['bool_refined_bounds'] = False
fea_config_dict['use_weights'] = 'O2' # this is not working
if 'Y730' in fea_config_dict['pc_name']:
    ################################################################
    # Y730
    ################################################################
    # Combined winding IM
    fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.975
    fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.025
    fea_config_dict['which_filter'] = 'VariableStatorSlotDepth'
    run_folder = r'run#550/'

    # # Separate winding IM
    # fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.60
    # fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.025
    # fea_config_dict['which_filter'] = 'VariableStatorSlotDepth'
    # run_folder = r'run#551/'

elif 'Severson01' in fea_config_dict['pc_name']:
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Severson01
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Combined winding IM
    fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.975
    fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.025
    fea_config_dict['which_filter'] = 'VariableStatorSlotDepth'
    run_folder = r'run#550010/'

elif 'Severson02' in fea_config_dict['pc_name']:
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Severson02
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Separate winding IM
    fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.60
    fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.025
    fea_config_dict['which_filter'] = 'VariableStatorSlotDepth'
    run_folder = r'run#550020/'

else: #T440p 
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # T440p
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Separate winding IM
    fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.60
    fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.025
    fea_config_dict['which_filter'] = 'VariableStatorSlotDepth'
    run_folder = r'run#550040/'

fea_config_dict['run_folder'] = run_folder

# spec's
my_execfile('./spec_TIA_ITEC_.py', g=globals(), l=locals())
spec
fea_config_dict['Active_Qr'] = spec.Qr
fea_config_dict['use_drop_shape_rotor_bar'] = spec.use_drop_shape_rotor_bar
build_model_name_prefix(fea_config_dict) # rebuild the model name for fea_config_dict
spec.build_im_template(fea_config_dict)

# select motor type ehere
print('Build ACM template...')
spec.acm_template = spec.im_template


import acm_designer
global ad
ad = acm_designer.acm_designer(fea_config_dict, spec)
# if 'Y730' in fea_config_dict['pc_name']:
#     ad.build_oneReport() # require LaTeX
#     ad.talk_to_mysql_database() # require MySQL
#     quit()
ad.init_logger()
ad.bool_re_evaluate = bool_re_evaluate


ad.bounds_denorm = spec.get_im_classic_bounds(which_filter=fea_config_dict['which_filter'])
ad.bound_filter  = spec.bound_filter
print('---------------------\nBounds:')
idx_ad = 0
for idx, f in enumerate(ad.bound_filter):
    if f == True:
        print(idx, f, '[%g,%g]'%tuple(spec.original_template_neighbor_bounds[idx]), '[%g,%g]'%tuple(ad.bounds_denorm[idx_ad]))
        idx_ad += 1
    else:
        print(idx, f, '[%g,%g]'%tuple(spec.original_template_neighbor_bounds[idx]))
# quit()

if ad.bool_re_evaluate:
    ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + ad.solver.fea_config_dict['run_folder']
    number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
    print('Count of chromosomes:', len(ad.solver.swarm_data))

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Optimization
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
import pygmo as pg
global counter_fitness_called, counter_fitness_return
from acm_designer import get_bad_fintess_values

class Problem_BearinglessInductionDesign(object):

    # Define objectives
    def fitness(self, x):
        global ad, counter_fitness_called, counter_fitness_return

        if ad.flag_do_not_evaluate_when_init_pop == True:
            return [0, 0, 0]

        ad, counter_fitness_called, counter_fitness_return
        if counter_fitness_called == counter_fitness_return:
            counter_fitness_called += 1
        else:
            # This is not reachable
            raise Exception('counter_fitness_called')
        print('Call fitness: %d, %d'%(counter_fitness_called, counter_fitness_return))

        # 不要标幺化了！统一用真的bounds，见get_bounds()
        x_denorm = x

        # evaluate x_denorm via FEA tools
        counter_loop = 0
        stuck_at = 0
        while True:
            if ad.bool_re_evaluate:
                if counter_fitness_return >= len(ad.solver.swarm_data):
                    quit()
                x_denorm = ad.solver.swarm_data[counter_fitness_return][:-3]
                print(ad.solver.swarm_data[counter_fitness_return])

            if stuck_at < counter_fitness_called:
                stuck_at = counter_fitness_called
                counter_loop = 0 # reset
            if stuck_at == counter_fitness_called:
                counter_loop += 1
                if counter_loop > 3:
                    raise Exception('Abort the optimization. Three attemps to evaluate the design have all failed for individual #%d'%(counter_fitness_called))

            try:
                cost_function, f1, f2, f3, FRW, \
                normalized_torque_ripple, \
                normalized_force_error_magnitude, \
                force_error_angle = \
                    ad.evaluate_design(ad.spec.im_template, x_denorm, counter_fitness_called, counter_loop=counter_loop)

                # remove folder .jfiles to save space (we have to generate it first in JMAG Designer to have field data and voltage profiles)
                if ad.solver.folder_to_be_deleted is not None and os.path.isdir(ad.solver.folder_to_be_deleted):
                    try:
                        shutil.rmtree(ad.solver.folder_to_be_deleted) # .jfiles directory
                    except PermissionError as error:
                        print(error)
                        print('Skip deleting this folder...')

                # update to be deleted when JMAG releases the use
                ad.solver.folder_to_be_deleted = ad.solver.expected_project_file[:-5]+'jfiles'

            except utility.ExceptionBadNumberOfParts as error:
                print(str(error)) 
                print("Detail: {}".format(error.payload))
                f1, f2, f3 = get_bad_fintess_values()
                utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nExceptionBadNumberOfParts:' + str(error) + '\n'*3 + "Detail: {}".format(error.payload))
                break

            except (utility.ExceptionReTry, pywintypes.com_error) as error:
                print(error)

                msg = 'FEA tool failed for individual #%d: attemp #%d.'%(counter_fitness_called, counter_loop)
                logger = logging.getLogger(__name__)
                logger.error(msg)
                print(msg)

                # if False:
                    # msg = 'Removing all files for individual #%d and try again...'%(counter_fitness_called)
                    # logger.error(msg)
                    # print(msg)
                    # try:
                    #         # turn off JMAG Designer
                    #         # try:
                    #         #     ad.solver.app.Quit()
                    #         # except:
                    #         #     print('I think there is no need to Quit the app')
                    #     ad.solver.app = None

                    #     # JMAG files
                    #     # os.remove(ad.solver.expected_project_file) # .jproj
                    #     # shutil.rmtree(ad.solver.expected_project_file[:-5]+'jfiles') # .jfiles directory # .jplot file in this folder will be used by JSOL softwares even JMAG Designer is closed.

                    #     # FEMM files
                    #     if os.path.exists(ad.solver.femm_output_file_path):
                    #         os.remove(ad.solver.femm_output_file_path) # .csv
                    #     if os.path.exists(ad.solver.femm_output_file_path[:-3]+'fem'):
                    #         os.remove(ad.solver.femm_output_file_path[:-3]+'fem') # .fem
                    #     for file in os.listdir(ad.solver.dir_femm_temp):
                    #         if 'femm_temp_' in file or 'femm_found' in file:
                    #             os.remove(ad.solver.dir_femm_temp + file)

                    # except Exception as e2:
                    #     utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nException 1:' + str(error) + '\n'*3 + 'Exception 2:' + str(e2))
                    #     raise e2
                from time import sleep
                print('\n\n\nSleep for 3 sec and continue.')
                sleep(3)

                continue

            except AttributeError as error:
                print(str(error)) 
                print("Detail: {}".format(error.payload))

                msg = 'FEA tool failed for individual #%d: attemp #%d.'%(counter_fitness_called, counter_loop)
                logger = logging.getLogger(__name__)
                logger.error(msg)
                print(msg)

                if 'designer.Application' in str(error):
                    from time import sleep
                    print('\n\n\nSleep for 3 sec and continue.')
                    sleep(3)

                    continue
                else:
                    raise error

            except Exception as e: # raise and need human inspection

                print('-'*40 + 'Unexpected error is caught.')
                print(str(e)) 
                utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nUnexpected expection:' + str(e))
                raise e

            else:
                break

        # - Torque per Rotor Volume
        f1 #= - ad.spec.required_torque / rated_rotor_volume
        # - Efficiency @ Rated Power
        f2 #= - rated_efficiency
        # Ripple Performance (Weighted Sum)
        f3 #= sum(list_weighted_ripples)
        print('f1,f2,f3:',f1,f2,f3)

        # Constraints (Em<0.2 and Ea<10 deg):
        # if abs(normalized_torque_ripple)>=0.2 or abs(normalized_force_error_magnitude) >= 0.2 or abs(force_error_angle) > 10:
        if abs(normalized_torque_ripple)>=0.3 or abs(normalized_force_error_magnitude) >= 0.3 or abs(force_error_angle) > 10 or FRW < 0.75:
            print('Constraints are violated:')
            if abs(normalized_torque_ripple)>=0.3:
                print('\tabs(normalized_torque_ripple)>=0.3')
            if abs(normalized_force_error_magnitude) >= 0.3:
                print('\tabs(normalized_force_error_magnitude) >= 0.3')
            if abs(force_error_angle) > 10:
                print('\tabs(force_error_angle) > 10')
            if FRW < 0.75:
                print('\tFRW < 0.75')
            f1, f2, f3 = get_bad_fintess_values()
        print('f1,f2,f3:',f1,f2,f3)

        counter_fitness_return += 1
        print('Fitness: %d, %d\n----------------'%(counter_fitness_called, counter_fitness_return))
        return [f1, f2, f3]

    # Return number of objectives
    def get_nobj(self):
        return 3

    # Return bounds of decision variables (a.k.a. chromosome)
    def get_bounds(self):
        global ad

        # denormalize the normalized chromosome x to x_denorm
        min_b, max_b = np.asarray(ad.bounds_denorm).T 
        # diff = np.fabs(min_b - max_b)
        # x_denorm = min_b + x * diff

        # print(min_b.tolist(), max_b.tolist())
        # print(([0]*7, [1]*7))
        # quit()
        return ( min_b.tolist(), max_b.tolist() )
        # return ([0]*7, [1]*7)

    # Return function name
    def get_name(self):
        return "Bearingless Induction Motor Design"



if bool_post_processing == True:
    # Combine all data 

    # Select optimal design by user-defined criteria
    if r'run#540' in fea_config_dict['run_folder']:
        ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + r'run#540011/' # severson01
        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
        swarm_data_severson01 = ad.solver.swarm_data

        def selection_criteria(swarm_data_):
            global best_idx, best_chromosome
            # HIGH TORQUE DENSITY
            # Severson01
            #       L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.
            # 1624 [1.15021, 8.97302, 8.33786, 3.22996, 0.759612, 2.81857, 1.11651, -22668.7, -0.953807, 4.79617]
            # Y730
            #       L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.
            # 794 [1.16163, 9.00566, 8.34039, 2.91474, 0.786231, 2.76114, 1.2485, -22666.7, -0.953681, 4.89779]
            # ----------------------------------------
            # HIGH EFFICIENCY
            # Y730
            #       L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.
            # 186 [1.25979, 6.80075, 5.64435, 5.8548, 1.59461, 2.11656, 2.58401, -17633.3, -0.958828, 5.53104]
            # 615 [1.2725, 5.6206, 4.60947, 3.56502, 2.27635, 0.506179, 2.78758, -17888.9, -0.958846, 8.56211]
            # ----------------------------------------
            # LOW RIPPLE PERFORMANCE
            # Severson02
            #       L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.
            # 1043 [1.38278, 8.91078, 7.43896, 2.66259, 0.611812, 1.50521, 1.51402, -19125.4, -0.953987, 2.91096]
            # 1129 [1.38878, 8.68378, 7.97301, 2.82904, 0.586374, 1.97867, 1.45825, -19169.0, -0.954226, 2.99944]
            # 1178 [1.36258, 8.9625, 7.49621, 2.5878, 0.503512, 0.678909, 1.74283, -19134.3, -0.952751, 2.90795]

            for idx, chromosome in enumerate(swarm_data_):
                # if chromosome[-1] < 5 and chromosome[-2] < -0.95 and chromosome[-3] < -22500: # best Y730     #1625, 0.000702091 * 8050 * 9.8 = 55.38795899 N.  FRW = 223.257 / 55.38795899 = 4.0
                # if chromosome[-1] < 10 and chromosome[-2] < -0.9585 and chromosome[-3] < -17500: # best Y730  #187, 0.000902584 * 8050 * 9.8 = 71.204851760 N. FRW = 151.246 / 71.204851760 = 2.124
                if chromosome[-1] < 3 and chromosome[-2] < -0.95 and chromosome[-3] < -19000: # best severson02 #1130, 0.000830274 * 8050 * 9.8 = 65.50031586 N.  FRW = 177.418 / 65.5 = 2.7
                    print(idx, chromosome)

                    def pyx_script():
                        # Plot cross section view
                        import population
                        im_best = population.bearingless_induction_motor_design.local_design_variant(ad.spec.im_template, 99, 999, best_chromosome[:-3])
                        im_best.ID = str(best_idx)
                        pyx_draw_model(im_best)
                        quit()


                    # # Take high torque density design for LSA
                    # if idx == 1625 - 1:
                    #     best_idx = idx
                    #     best_chromosome = chromosome
                    #     pyx_script()

                    # # Take high efficiency design for LSA
                    # if idx == 187 - 1:
                    #     best_idx = idx
                    #     best_chromosome = chromosome
                    #     pyx_script()

                    # Take low ripple performance design for LSA
                    if idx == 1130 - 1:
                        best_idx = idx
                        best_chromosome = chromosome
                        # pyx_script()


        print('-'*40+'\nSeverson01' + '\n      L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.')
        selection_criteria(swarm_data_severson01)

        ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + r'run#540021/' # severson02
        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
        swarm_data_severson02 = ad.solver.swarm_data

        print('-'*40+'\nSeverson02' + '\n      L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.')
        selection_criteria(swarm_data_severson02)

        # swarm_data_severson01 = []
        # swarm_data_severson02 = []
        ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + r'run#540/' # ad.solver.fea_config_dict['run_folder'] 
        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
        swarm_data_Y730 = ad.solver.swarm_data

        print('-'*40+'\nY730' + '\n      L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.')
        selection_criteria(swarm_data_Y730)

        # Set the output_dir back!
        ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + ad.solver.fea_config_dict['run_folder']
        # quit()

    elif fea_config_dict['run_folder'] == r'run#538/':
        ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + r'run#538011/' # severson01
        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
        swarm_data_severson01 = ad.solver.swarm_data

        ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + r'run#538021/' # severson02
        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
        swarm_data_severson02 = ad.solver.swarm_data

        ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + r'run#538/' # ad.solver.fea_config_dict['run_folder'] 
        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
        swarm_data_Y730 = ad.solver.swarm_data

    print('Sizes of the 3 populations (in order):', len(swarm_data_severson01), len(swarm_data_severson02), len(swarm_data_Y730))
    ad.solver.swarm_data = swarm_data_severson01 + swarm_data_severson02 + swarm_data_Y730 # list add

    udp = Problem_BearinglessInductionDesign()
    ad.flag_do_not_evaluate_when_init_pop = True
    counter_fitness_called, counter_fitness_return = 0, 0
    prob = pg.problem(udp)

    # LSA
    if fea_config_dict['local_sensitivity_analysis'] == True:

        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)

        if number_of_chromosome is not None:
            ad.solver.swarm_data

            # Learn Pareto front rank and plot
            for el in ad.solver.swarm_data:
                print('\t', el)
            print('count:', len(ad.solver.swarm_data))
            swarm_data_on_pareto_front = learn_about_the_archive(prob, ad.solver.swarm_data, len(ad.solver.swarm_data), fea_config_dict)

            # plot LSA
            ad.solver.swarm_data_container.sensitivity_bar_charts()

            plt.show()
            quit()


        else:
            def local_sensitivity_analysis(reference_design_denorm, percent=0.2):
                # 敏感性检查：以基本设计为准，检查不同的参数取极值时的电机性能变化！这是最简单有效的办法。七个设计参数，那么就有14种极值设计。

                if False:
                    # Within the original bounds
                    min_b, max_b = udp.get_bounds()
                    min_b, max_b = np.array(min_b), np.array(max_b)
                    diff = np.fabs(min_b - max_b)
                else:
                    # near the reference design
                    min_b = [el*(1.0-percent) for el in reference_design_denorm]
                    max_b = [el*(1.0+percent) for el in reference_design_denorm]
                    min_b, max_b = np.array(min_b), np.array(max_b)
                    diff = np.fabs(min_b - max_b)

                reference_design = (np.array(reference_design_denorm) - min_b) / diff
                print('reference_design_denorm:', reference_design_denorm)
                print('reference_design:\t\t', reference_design.tolist())
                base_design = reference_design.tolist()
                # quit()
                number_of_variants = fea_config_dict['local_sensitivity_analysis_number_of_variants']
                lsa_swarm = [base_design] # include reference design!
                for i in range(len(base_design)): # 7 design parameters
                    for j in range(number_of_variants+1): # 21 variants interval
                        # copy list
                        design_variant = base_design[::]
                        design_variant[i] = j * 1./number_of_variants
                        lsa_swarm.append(design_variant)

                lsa_swarm_denorm = min_b + lsa_swarm * diff 
                print(lsa_swarm)
                print(lsa_swarm_denorm)
                return lsa_swarm, lsa_swarm_denorm

            print('Best index', best_idx, '#%d'%(best_idx+1), 'Best chromosome', best_chromosome)
            _, lsa_swarm_denorm = local_sensitivity_analysis( reference_design_denorm=best_chromosome[:-3],
                                                              percent=fea_config_dict['local_sensitivity_analysis_percent'] )
            lsa_popsize = len(lsa_swarm_denorm)

            # quit()
            lsa_pop = pg.population(prob, size=lsa_popsize)
            print('Set ad.flag_do_not_evaluate_when_init_pop to False...')
            ad.flag_do_not_evaluate_when_init_pop = False
            for i, design_denorm in enumerate(lsa_swarm_denorm):
                print('Evaluate', i)
                lsa_pop.set_x(i, design_denorm)
            print('LSA is done for a pop size of %d'%(lsa_popsize))
        quit()

    # plot pareto plot for three objectives...
    else:
        print('------------for tutorial\n'*3)
        print(ad.solver.swarm_data)
        swarm_data_on_pareto_front = learn_about_the_archive(prob, ad.solver.swarm_data, len(ad.solver.swarm_data), len_s01=len(swarm_data_severson01), len_s02=len(swarm_data_severson02))
        plt.show()    
        quit()

        # # Reproduce a design 
        # cost_function, f1, f2, f3, \
        # normalized_torque_ripple, \
        # normalized_force_error_magnitude, \
        # force_error_angle = \
        #     ad.evaluate_design(ad.spec.im_template, best_chromosome[:-3], 1130)
        quit()
        run_static_structural_fea(swda.best_design_denorm)


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Multi-Objective Optimization
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
if True:
################################################################
# MOO Step 1:
#   Create UserDefinedProblem and create population
#   The magic method __init__ cannot be fined for UDP class
################################################################
    udp = Problem_BearinglessInductionDesign()
    counter_fitness_called, counter_fitness_return = 0, 0
    prob = pg.problem(udp)

    popsize = 78
        # Traceback (most recent call last):
        #   File "D:\OneDrive - UW-Madison\c\codes3\one_script.py", line 1189, in <module>
        #     pop = algo.evolve(pop)
        # ValueError: 
        # function: decomposition_weights
        # where: C:\bld\pygmo_1557474762576\_h_env\Library\include\pagmo/utils/multi_objective.hpp, 642
        # what: Population size of 72 is detected, but not supported by the 'grid' weight generation method selected. A size of 66 or 78 is possible.
    print('-'*40 + '\nPop size is', popsize)

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Add Restarting Feature
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

    # 检查swarm_data.txt，如果有至少一个数据，返回就不是None。
    print('Check swarm_data.txt...')
    number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
    if number_of_chromosome is not None:
        # 禁止在初始化pop时运行有限元
        ad.flag_do_not_evaluate_when_init_pop               = True

        number_of_finished_iterations                       = number_of_chromosome // popsize
        number_of_finished_chromosome_in_current_generation = number_of_chromosome % popsize

        # 如果刚好整除，把余数0改为popsize
        if number_of_finished_chromosome_in_current_generation == 0:
            number_of_finished_chromosome_in_current_generation = popsize
            print(f'\tThere are {number_of_chromosome} chromosomes found in swarm_data.txt.')
            print('\tWhat is the odds! The script just stopped when the evaluation of the whole pop is finished.')
            print('\tSet number_of_finished_chromosome_in_current_generation to popsize %d'%(number_of_finished_chromosome_in_current_generation))

        print('This is a restart of '+ fea_config_dict['run_folder'][:-1])
        print('\tNumber of finished iterations is %d'%(number_of_finished_iterations))
        # print('This means the initialization of the population class is interrupted. So the pop in swarm_data.txt is used as the survivor.')

        # 继续从swarm_survivor.txt中读取数据，注意，survivor总是是完整的一代的，除非popsize被修改了。
        print('\tCheck swarm_survivor.txt...', end='')
        ad.solver.survivor = ad.solver.read_swarm_survivor(popsize)

        # 如果发现ad.solver.survivor是None，那就说明是初始化pop的时候被中断了，此时就用swarm_data来生成pop。
        if ad.solver.survivor is not None:
            print('Found survivor!\nRestart the optimization based on the swarm_survivor.txt.')

            if len(ad.solver.survivor) != popsize:
                print('popsize is reduced') # 如果popsize增大了，read_swarm_survivor(popsize)就会报错了，因为-----不能被split后转为float
                raise Exception('This is a feature not tested. However, you can cheat to change popsize by manually modify swarm_data.txt or swarm_survivor.txt.')
        else:
            print('Gotta make do with swarm_data to generate survivor.')

        # 这些计数器的值永远都是评估过的chromosome的个数。
        counter_fitness_called = counter_fitness_return = number_of_chromosome
        print('counter_fitness_called = counter_fitness_return = number_of_chromosome = %d'%(number_of_chromosome))

    else:
        print('Nothing exists in swarm_data.txt.\nThis is a whole new run.')
        ad.flag_do_not_evaluate_when_init_pop = False
        number_of_finished_chromosome_in_current_generation = None
        number_of_finished_iterations = 0 # 实际上跑起来它不是零，而是一，因为我们认为初始化的一代也是一代。或者，我们定义number_of_finished_iterations = number_of_chromosome // popsize

    if bool_re_evaluate:
        counter_fitness_called = counter_fitness_return = 0        

    # 初始化population，如果ad.flag_do_not_evaluate_when_init_pop是False，那么就说明是 new run，否则，整代个体的fitness都是[0,0,0]。
    pop = pg.population(prob, size=popsize) 

    # 如果整代个体的fitness都是[0,0,0]，那就需要调用set_xf，把txt文件中的数据写入pop。如果发现数据的个数不够，那就调用set_x()来产生数据，形成初代个体。
    if ad.flag_do_not_evaluate_when_init_pop == True:
        pop_array = pop.get_x()
        if number_of_chromosome <= popsize:
            for i in range(popsize):
                if i < number_of_chromosome: #number_of_finished_chromosome_in_current_generation:
                    pop.set_xf(i, ad.solver.swarm_data[i][:-3], ad.solver.swarm_data[i][-3:])
                else:
                    print('Set ad.flag_do_not_evaluate_when_init_pop to False...')
                    ad.flag_do_not_evaluate_when_init_pop = False
                    print('Calling pop.set_x()---this is a restart for individual#%d during pop initialization.'%(i))
                    print(i, 'get_fevals:', prob.get_fevals())
                    pop.set_x(i, pop_array[i]) # evaluate this guy
        else:
            # 新办法，直接从swarm_data.txt（相当于archive）中判断出当前最棒的群体
            swarm_data_on_pareto_front = learn_about_the_archive(prob, ad.solver.swarm_data, popsize, fea_config_dict)
            # print(swarm_data_on_pareto_front)
            for i in range(popsize):
                pop.set_xf(i, swarm_data_on_pareto_front[i][:-3], swarm_data_on_pareto_front[i][-3:])

            # if False:
                # # 老办法（蠢办法），依赖于swarm_survivor.txt
                # for i in range(popsize):
                #     pop.set_xf(i, ad.solver.survivor[i][:-3], ad.solver.survivor[i][-3:])

                # # 小测试，该式子应该永远成立，如果不成立，说明分析完一代以后，write survivor没有被正常调用。
                # if ad.solver.survivor is not None:
                #     if ad.solver.survivor_title_number // popsize == number_of_finished_iterations:
                #         print('survivor_title_number is', ad.solver.survivor_title_number, 'number_of_chromosome is', number_of_chromosome)
                #         if ad.solver.survivor_title_number == number_of_chromosome:
                #             # 刚好swarm_data也是完整的一代
                #             print('\t刚刚好swarm_data也是完整的一代！！！')
                #     else:
                #         raise

                # # 手动把当前pop和swarm_data中的最新个体进行dominance比较
                #     # 经常出现的情况，survivor和swarm_data都存在，前者总是完整的一代，也就是说，
                #     # 搞到非初始化的某一代的中间的某个个体的时候断了，PyGMO不支持从中间搞起，那么只能我自己来根据swarm_data.txt中最后的数据来判断是否产生了值得留在种群中的个体了。
                # if ad.solver.survivor is not None and number_of_chromosome % popsize != 0: # number_of_finished_chromosome_in_current_generation < popsize: # recall if number_of_finished_chromosome_in_current_generation == 0, it is set to popsize.
                #     pop_array = pop.get_x()
                #     fits_array = pop.get_f()

                #     # number_of_finished_iterations = number_of_chromosome // popsize
                #     base = number_of_finished_iterations*popsize
                #     for i in range(number_of_chromosome % popsize):
                #         obj_challenger = ad.solver.swarm_data[i+base][-3:]
                #         if pg.pareto_dominance(obj_challenger, fits_array[i]):
                #             pop.set_xf(i, ad.solver.swarm_data[i+base][:-3] , obj_challenger)
                #             print(i+base, '\t', obj_challenger, '\n\t', fits_array[i].tolist())
                #             print('\t', pop.get_x()[i], pop.get_f()[i].tolist())
                #         else:
                #             print(i+base)

        # 必须放到这个if的最后，因为在 learn_about_the_archive 中是有初始化一个 pop_archive 的，会调用fitness方法。
        ad.flag_do_not_evaluate_when_init_pop = False

    print('-'*40, '\nPop is initialized:\n', pop)
    hv = pg.hypervolume(pop)
    quality_measure = hv.compute(ref_point=get_bad_fintess_values(ref=True)) # ref_point must be dominated by the pop's pareto front
    print('quality_measure: %g'%(quality_measure))
    # raise KeyboardInterrupt

    # 初始化以后，pop.problem.get_fevals()就是popsize，但是如果大于popsize，说明“pop.set_x(i, pop_array[i]) # evaluate this guy”被调用了，说明还没输出过 survivors 数据，那么就写一下。
    if pop.problem.get_fevals() > popsize:
        print('Write survivors.')
        ad.solver.write_swarm_survivor(pop, counter_fitness_return)

################################################################
# MOO Step 2:
#   Select algorithm (another option is pg.nsga2())
################################################################
    # Don't forget to change neighbours to be below popsize (default is 20) decomposition="bi"
    algo = pg.algorithm(pg.moead(gen=1, weight_generation="grid", decomposition="tchebycheff", 
                                 neighbours=20, 
                                 CR=1, F=0.5, eta_m=20, 
                                 realb=0.9, 
                                 limit=2, preserve_diversity=True)) # https://esa.github.io/pagmo2/docs/python/algorithms/py_algorithms.html#pygmo.moead
    print('-'*40, '\n', algo)
    # quit()

################################################################
# MOO Step 3:
#   Begin optimization
################################################################
    number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
    number_of_finished_iterations = number_of_chromosome // popsize
    number_of_iterations = 50
    logger = logging.getLogger(__name__)
    # try:
    if True:
        for _ in range(number_of_finished_iterations, number_of_iterations):
            msg = 'This is iteration #%d. '%(_)
            print(msg)
            logger.info(msg)
            pop = algo.evolve(pop)

            msg += 'Write survivors to file. '
            ad.solver.write_swarm_survivor(pop, counter_fitness_return)

            hv = pg.hypervolume(pop)
            quality_measure = hv.compute(ref_point=get_bad_fintess_values(ref=True)) # ref_point must be dominated by the pop's pareto front
            msg += 'Quality measure by hyper-volume: %g'% (quality_measure)
            print(msg)
            logger.info(msg)
            
            my_print(ad, pop, _)
            # my_plot(fits, vectors, ndf)
    # except Exception as e:
    #     print(pop.get_x())
    #     print(pop.get_f().tolist())
    #     raise e

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Weighted objective function optimmization
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
else:
    # 微分进化的配置
    ad.get_de_config()
    print('Run: ' + run_folder + '\nThe auto generated bounds are:', ad.de_config_dict['original_bounds'])
    # quit()

    # 如果需要局部敏感性分析，那就先跑了再说
    if fea_config_dict['local_sensitivity_analysis'] == True:
        if not ad.check_results_of_local_sensitivity_analysis():
            ad.run_local_sensitivity_analysis(ad.de_config_dict['original_bounds'], design_denorm=None)
        else:
            ad.de_config_dict['bounds'] = ad.de_config_dict['original_bounds']
            ad.init_swarm() # define ad.sw
            ad.sw.generate_pop(specified_initial_design_denorm=None)
        ad.collect_results_of_local_sensitivity_analysis() # need sw to work
        fea_config_dict['local_sensitivity_analysis'] = False # turn off lsa mode

    # Build the final bounds
    if fea_config_dict['bool_refined_bounds'] == -1:
        ad.build_local_bounds_from_best_design(None)
    elif fea_config_dict['bool_refined_bounds'] == True:
        ad.build_refined_bounds(ad.de_config_dict['original_bounds'])
    elif fea_config_dict['bool_refined_bounds'] == False:
        ad.de_config_dict['bounds'] = ad.de_config_dict['original_bounds']
    else:
        raise Exception('bool_refined_bounds')
    print('The final bounds are:')
    for el in ad.de_config_dict['bounds']:
        print('\t', el)

    if False == bool_post_processing:
        if fea_config_dict['flag_optimization'] == True:
            ad.init_swarm() # define ad.sw
            ad.run_de()
        else:
            print('Do something.')

    elif True == bool_post_processing:
        ad.init_swarm()
        swda = ad.best_design_by_weights(fea_config_dict['use_weights'])
        from pylab import show
        show()
        quit()
        run_static_structural_fea(swda.best_design_denorm)


