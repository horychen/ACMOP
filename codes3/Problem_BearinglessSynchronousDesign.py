import builtins
if hasattr(builtins, 'ad'):
    print('[Problem_BlessSyn] Global variables are shared between modules...')
else:
    raise Exception('Please add global variable (address) "ad" to module __builtins__.')
print('[Problem_BlessSyn]', builtins.ad)
print('[Problem_BlessSyn]', ad)
print('[Problem_BlessSyn]', ad.counter_fitness_called)
print('[Problem_BlessSyn]', ad.counter_fitness_return)

import logging, os, shutil
import numpy as np
from acm_designer import get_bad_fintess_values
import utility
class Problem_BearinglessSynchronousDesign(object):

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
        print('[Problem_BlessSyn] Call fitness: %d, %d'%(ad.counter_fitness_called, ad.counter_fitness_return))

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
                print('[Problem_BlessSyn]', ad.solver.swarm_data[ad.counter_fitness_return])

            if stuck_at < ad.counter_fitness_called:
                stuck_at = ad.counter_fitness_called
                counter_loop = 0 # reset
            if stuck_at == ad.counter_fitness_called:
                counter_loop += 1

            # if True:
            try:
                cost_function, f1, f2, f3, FRW, \
                normalized_torque_ripple, \
                normalized_force_error_magnitude, \
                force_error_angle = \
                    ad.evaluate_design_json_wrapper(ad.acm_template, x_denorm, ad.counter_fitness_called, counter_loop=counter_loop)

                # remove folder .jfiles to save space (we have to generate it first in JMAG Designer to have field data and voltage profiles)
                if ad.solver.folder_to_be_deleted is not None and os.path.isdir(ad.solver.folder_to_be_deleted):
                    try:
                        shutil.rmtree(ad.solver.folder_to_be_deleted) # .jfiles directory
                    except PermissionError as error:
                        print(error)
                        print('Skip deleting this folder...')
                # update to be deleted when JMAG releases the use
                ad.solver.folder_to_be_deleted = ad.solver.expected_project_file[:-5]+'jfiles'

            except KeyboardInterrupt as error:
                raise error

            # except utility.ExceptionReTry as error: # The copy region target is not found
            #     print(str(error))
            #     print('CJH: "the ind***TranPMSM_torque.csv is not found" means the mesher or the solver has failed. For now, simply consider it to be bad design.')
            #     f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')
            #     logger = logging.getLogger(__name__)
            #     logger.error(str(error))
            #     break

            except utility.ExceptionBadNumberOfParts as error:
                print('ExceptionBadNumberOfParts captured:', str(error)) 
                # print("Detail: {}".format(error.payload))
                f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')
                # utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nExceptionBadNumberOfParts:' + str(error) + '\n'*3)
                break

            except Exception as error:
                # if ad.bool_re_evaluate == True:
                #     print('bool_re_evaluate is True...')
                #     raise error

                raise error
                f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')
                print(str(error))
                logger = logging.getLogger(__name__)
                logger.error(str(error))
                break
                # except FileNotFoundError as error: # The copy region target is not found
                #     print(str(error))
                #     print('CJH: "the ind***TranPMSM_torque.csv is not found" means the mesher or the solver has failed. For now, simply consider it to be bad design.')
                #     f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')

                # except utility.ExceptionBadNumberOfParts as error:
                #     print(str(error)) 
                #     # print("Detail: {}".format(error.payload))
                #     f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')
                #     utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nExceptionBadNumberOfParts:' + str(error) + '\n'*3)
                #     break

                # except (utility.ExceptionReTry, pywintypes.com_error) as error:
                #     print(error)

                #     msg = 'FEA tool failed for individual #%d: attemp #%d.'%(ad.counter_fitness_called, counter_loop)
                #     logger = logging.getLogger(__name__)
                #     logger.error(msg)
                #     print(msg)

                #     if counter_loop > 1: # > 1 = two attemps; > 2 = three attemps
                #         print(error)
                #         raise Exception('Abort the optimization. Two attemps to evaluate the design have all failed for individual #%d'%(ad.counter_fitness_called))
                #     else:
                #         from time import sleep
                #         print('\n\n\nSleep for 3 sec and continue.')
                #         sleep(3)
                #         continue

                # except AttributeError as error:
                #     print(str(error)) 
                #     # print("Detail: {}".format(error.payload))

                #     msg = 'FEA tool failed for individual #%d: attemp #%d.'%(ad.counter_fitness_called, counter_loop)
                #     logger = logging.getLogger(__name__)
                #     logger.error(msg)
                #     print(msg)

                #     if 'designer.Application' in str(error):
                #         if counter_loop > 1: 
                #             print(error)
                #             raise Exception('Abort the optimization. Two attemps to evaluate the design have all failed for individual #%d'%(ad.counter_fitness_called))
                #         else:
                #             from time import sleep
                #             print('\n\n\nSleep for 3 sec and continue.')
                #             sleep(3)                        
                #             continue
                #     else:
                #         raise error

                # except Exception as e: # raise and need human inspection

                #     # raise e
                #     print('-'*40 + 'Unexpected error is caught.')
                #     print(str(e)) 
                #     utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nUnexpected expection:' + str(e))
                #     raise e
            else:
                # - Price
                f1 
                # - Efficiency @ Rated Power
                f2 
                # Ripple Performance (Weighted Sum)
                f3 
                print('f1,f2,f3:',f1,f2,f3)

                try:
                    # Constraints (Em<0.2 and Ea<10 deg):
                    # if abs(normalized_torque_ripple)>=0.2 or abs(normalized_force_error_magnitude) >= 0.2 or abs(force_error_angle) > 10 or SafetyFactor < 1.5:
                    # if abs(normalized_torque_ripple)>=0.2 or abs(normalized_force_error_magnitude) >= 0.2 or abs(force_error_angle) > 10 or FRW < 1:
                    # if abs(normalized_torque_ripple)>=0.2 or abs(normalized_force_error_magnitude) >= 0.2 or abs(force_error_angle) > 10:
                    if abs(normalized_torque_ripple)>=0.3 or abs(normalized_force_error_magnitude) >= 0.35 or abs(force_error_angle) > 20 or FRW < 0.5:
                        print('[Problem_BlessSyn] Constraints are violated:')
                        if abs(normalized_torque_ripple)>=0.3:
                            print('\tabs(normalized_torque_ripple)>=0.3 | (=%f)' % (normalized_torque_ripple))
                        if abs(normalized_force_error_magnitude) >= 0.35:
                            print('\tabs(normalized_force_error_magnitude) >= 0.35 | (=%f)' % (normalized_force_error_magnitude))
                        if abs(force_error_angle) > 20:
                            print('\tabs(force_error_angle) > 20 | (=%f)' % (force_error_angle))
                        if FRW < 0.5:
                            print('\tFRW < 0.5 | (=%f)' % (FRW))
                        f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')
                    print('f1,f2,f3:',f1,f2,f3)
                except:
                    msg = 'This design causes an error in JMAG and hence is discarded..'
                    print(msg)
                    logger = logging.getLogger(__name__)
                    logger.warn(msg)

                break

        ad.counter_fitness_return += 1
        print('[Problem_BlessSyn] Fitness: %d, %d\n----------------'%(ad.counter_fitness_called, ad.counter_fitness_return))
        # raise KeyboardInterrupt
        return [f1, f2, f3]

    # Return number of objectives
    def get_nobj(self):
        return 3

    # Return bounds of decision variables (a.k.a. chromosome)
    def get_bounds(self):
        global ad
        print('[Problem_BlessSyn] Problem_BearinglessSynchronousDesign.get_bounds:', ad.acm_template.bounds_denorm)
        min_b, max_b = np.asarray(ad.acm_template.bounds_denorm).T 
        return ( min_b.tolist(), max_b.tolist() )

    # Return function name
    def get_name(self):
        return "Bearingless PMSM Design"

import pygmo as pg
def get_prob_and_popsize():
    udp = Problem_BearinglessSynchronousDesign()
    prob = pg.problem(udp)
    popsize = 78
    return udp, prob, popsize
