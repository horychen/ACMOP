import builtins
if hasattr(builtins, 'ad'):
    print('Global variables are shared between modules...')
else:
    raise Exception('Please add global variable (address) "ad" to module __builtins__.')
print(builtins.ad)
print(ad)
print(ad.counter_fitness_called)
print(ad.counter_fitness_return)

import logging, os, shutil
import numpy as np
from acm_designer import get_bad_fintess_values
import utility
import pywintypes
class Problem_BearinglessInductionDesign(object):

    # Define objectives
    def fitness(self, x):
        global ad

        if ad.flag_do_not_evaluate_when_init_pop == True:
            return [0, 0, 0]

        if ad.counter_fitness_called == ad.counter_fitness_return:
            ad.counter_fitness_called += 1
        else:
            # This is not reachable
            raise Exception('ad.counter_fitness_called')
        print('Call fitness: %d, %d'%(ad.counter_fitness_called, ad.counter_fitness_return))

        # 不要标幺化了！统一用真的bounds，见get_bounds()
        x_denorm = x

        # evaluate x_denorm via FEA tools
        counter_loop = 0
        stuck_at = 0
        while True:
            if ad.fea_config_dict['bool_re_evaluate']:
                if ad.counter_fitness_return >= len(ad.solver.swarm_data):
                    quit()
                x_denorm = ad.solver.swarm_data[ad.counter_fitness_return][:-3]
                print(ad.solver.swarm_data[ad.counter_fitness_return])

            if stuck_at < ad.counter_fitness_called:
                stuck_at = ad.counter_fitness_called
                counter_loop = 0 # reset
            if stuck_at == ad.counter_fitness_called:
                counter_loop += 1
                if counter_loop > 3:
                    raise Exception('Abort the optimization. Three attemps to evaluate the design have all failed for individual #%d'%(ad.counter_fitness_called))

            # try:
            if True:
                acm_variant = ad.evaluate_design_json_wrapper(ad.spec.acm_template, x_denorm, ad.counter_fitness_called, counter_loop=counter_loop)

                cost_function, f1, f2, f3, FRW, \
                normalized_torque_ripple, \
                normalized_force_error_magnitude, \
                force_error_angle = acm_variant.results_for_optimization

                # remove folder .jfiles to save space (we have to generate it first in JMAG Designer to have field data and voltage profiles)
                if ad.solver.folder_to_be_deleted is not None and os.path.isdir(ad.solver.folder_to_be_deleted):
                    try:
                        shutil.rmtree(ad.solver.folder_to_be_deleted) # .jfiles directory
                    except PermissionError as error:
                        print(error)
                        print('Skip deleting this folder...')

                # update to be deleted when JMAG releases the use
                ad.solver.folder_to_be_deleted = ad.solver.expected_project_file[:-5]+'jfiles'

            # except utility.ExceptionBadNumberOfParts as error:
            #     print(type(error), str(error)) 
            #     print("Detail: {}".format(error.payload))
            #     f1, f2, f3 = get_bad_fintess_values()
            #     utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nExceptionBadNumberOfParts:' + str(error) + '\n'*3 + "Detail: {}".format(error.payload))
            #     break

            # except (utility.ExceptionReTry, pywintypes.com_error) as error:
            #     print(type(error), error)

            #     msg = 'FEA tool failed for individual #%d: attemp #%d.'%(ad.counter_fitness_called, counter_loop)
            #     logger = logging.getLogger(__name__)
            #     logger.error(msg)
            #     print(msg)

            #     # if False:
            #         # msg = 'Removing all files for individual #%d and try again...'%(ad.counter_fitness_called)
            #         # logger.error(msg)
            #         # print(msg)
            #         # try:
            #         #         # turn off JMAG Designer
            #         #         # try:
            #         #         #     ad.solver.app.Quit()
            #         #         # except:
            #         #         #     print('I think there is no need to Quit the app')
            #         #     ad.solver.app = None

            #         #     # JMAG files
            #         #     # os.remove(ad.solver.expected_project_file) # .jproj
            #         #     # shutil.rmtree(ad.solver.expected_project_file[:-5]+'jfiles') # .jfiles directory # .jplot file in this folder will be used by JSOL softwares even JMAG Designer is closed.

            #         #     # FEMM files
            #         #     if os.path.exists(ad.solver.femm_output_file_path):
            #         #         os.remove(ad.solver.femm_output_file_path) # .csv
            #         #     if os.path.exists(ad.solver.femm_output_file_path[:-3]+'fem'):
            #         #         os.remove(ad.solver.femm_output_file_path[:-3]+'fem') # .fem
            #         #     for file in os.listdir(ad.solver.dir_femm_temp):
            #         #         if 'femm_temp_' in file or 'femm_found' in file:
            #         #             os.remove(ad.solver.dir_femm_temp + file)

            #         # except Exception as e2:
            #         #     utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nException 1:' + str(error) + '\n'*3 + 'Exception 2:' + str(e2))
            #         #     raise e2
            #     from time import sleep
            #     print('\n\n\nSleep for 3 sec and continue.')
            #     sleep(3)

            #     continue

            # except AttributeError as error:
            #     print(type(error), str(error)) 
            #     print("Detail: {}".format(error.payload))

            #     msg = 'FEA tool failed for individual #%d: attemp #%d.'%(ad.counter_fitness_called, counter_loop)
            #     logger = logging.getLogger(__name__)
            #     logger.error(msg)
            #     print(msg)

            #     if 'designer.Application' in str(error):
            #         from time import sleep
            #         print('\n\n\nSleep for 3 sec and continue.')
            #         sleep(3)

            #         continue
            #     else:
            #         raise error

            # except Exception as e: # raise and need human inspection

            #     print('-'*40 + 'Unexpected error is caught.')
            #     print(type(error), str(e)) 
            #     utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nUnexpected expection:' + str(e))
            #     raise e

            # else:
            #     break
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

        ad.counter_fitness_return += 1
        print('Fitness: %d, %d\n----------------'%(ad.counter_fitness_called, ad.counter_fitness_return))
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

import pygmo as pg
def get_prob_and_popsize():
    udp = Problem_BearinglessInductionDesign()
    prob = pg.problem(udp)
    popsize = 78
    return udp, prob, popsize

print('Module Problem_BearinglessInductionDesign is imported...')

