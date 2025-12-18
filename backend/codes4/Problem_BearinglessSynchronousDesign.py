import numpy as np
# from acm_designer import get_bad_fintess_values
import utility
import logging, os, shutil
import pywintypes
import builtins
logger = logging.getLogger(__name__)
if hasattr(builtins, 'ad'):
    logger.info('Global variable ad is shared between modules as we cannot pass new argument to udp class.')
else:
    raise Exception('[Problem_BlessSyn] Please add global variable (address) "ad" to module __builtins__.')

# print('[Problem_BlessSyn]', builtins.ad)

# print('[Problem_BlessSyn]', ad)
# print('[Problem_BlessSyn]', ad.counter_fitness_called)
# print('[Problem_BlessSyn]', ad.counter_fitness_return)

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
            raise Exception(f'ad.counter_fitness_called = {ad.counter_fitness_called} != ad.counter_fitness_return = {ad.counter_fitness_return}!!!')
        logger.info('-'*40)
        logger.debug('Call fitness: %d, %d', ad.counter_fitness_called, ad.counter_fitness_return)

        # 不要标幺化了！统一用真的bounds，见get_bounds()
        x_denorm = x

        # evaluate x_denorm via FEA tools
        counter_loop = 0
        stuck_at = 0
        while True:
            if stuck_at < ad.counter_fitness_called:
                stuck_at = ad.counter_fitness_called
                counter_loop = 0 # reset
            if stuck_at == ad.counter_fitness_called:
                counter_loop += 1

            # if True:
            try:
                ad.name = ad.machine_class + f'gen-{ad.generation}-ind-{ad.counter_fitness_called}'
                cost_function, f1, f2, f3, FRW, \
                normalized_torque_ripple, \
                normalized_force_error_magnitude, \
                force_error_angle = ad.evaluate_design_json_wrapper(x_denorm, ad.counter_fitness_called, counter_loop=counter_loop)

                # For JMAG, remove folder ".jfiles" to save space (we have to generate it first in JMAG Designer to have field data and voltage profiles)
                # if 'JMAG' in ad.select_fea_config_dict:
                #     if ad.folder_to_be_deleted is not None and os.path.isdir(ad.folder_to_be_deleted):
                #         try:
                #             shutil.rmtree(ad.   folder_to_be_deleted) # .jfiles directory
                #         except PermissionError as error:
                #             logger.warning('PermissionError: %s', error)
                #             logger.warning('Skip deleting this folder...')
                #     # update to be deleted when JMAG releases the use
                #     ad.folder_to_be_deleted = ad.expected_project_file[:-5]+'jfiles'

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
                logger.error('ExceptionBadNumberOfParts captured: %s', str(error)) 
                # print("Detail: {}".format(error.payload))
                f1, f2, f3 = ad.get_bad_fintess_values(machine_type=ad.name)
                # utility.send_notification(ad.solver.fea_config_dict['pc_name'] + '\n\nExceptionBadNumberOfParts:' + str(error) + '\n'*3)
                break

            except pywintypes.com_error as error:
                logger.error('pywintypes.com_error: %s', error)
                logger.error('The call to JMAG has failed. Restart?')
                raise error

            except Exception as error:
                # if ad.bool_re_evaluate == True:
                #     print('bool_re_evaluate is True...')
                #     raise error

                raise error
                f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')
                logger.error('Exception: %s', str(error))
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
                pass
            else:
                # - Price
                f1 
                # - Efficiency @ Rated Power
                f2 
                # Ripple Performance (Weighted Sum)
                f3 
                logger.debug('f1,f2,f3: %s, %s, %s', f1, f2, f3)

                if ad.fea_config_dict['moo.apply_constraints']==True:
                    # Constraints (Em<0.2 and Ea<10 deg):
                    # if abs(normalized_torque_ripple)>=0.2 or abs(normalized_force_error_magnitude) >= 0.2 or abs(force_error_angle) > 10 or SafetyFactor < 1.5:
                    # if abs(normalized_torque_ripple)>=0.2 or abs(normalized_force_error_magnitude) >= 0.2 or abs(force_error_angle) > 10 or FRW < 1:
                    # if abs(normalized_torque_ripple)>=0.2 or abs(normalized_force_error_magnitude) >= 0.2 or abs(force_error_angle) > 10:
                    if abs(normalized_torque_ripple)>=0.3 or abs(normalized_force_error_magnitude) >= 0.35 or abs(force_error_angle) > 20 or FRW < 0.5:
                        logger.warning('Constraints are violated:')
                        if abs(normalized_torque_ripple)>=0.3:
                            logger.warning('\tabs(normalized_torque_ripple)>=0.3 | (=%f)', normalized_torque_ripple)
                        if abs(normalized_force_error_magnitude) >= 0.35:
                            logger.warning('\tabs(normalized_force_error_magnitude) >= 0.35 | (=%f)', normalized_force_error_magnitude)
                        if abs(force_error_angle) > 20:
                            logger.warning('\tabs(force_error_angle) > 20 | (=%f)', force_error_angle)
                        if FRW < 0.5:
                            logger.warning('\tFRW < 0.5 | (=%f)', FRW)
                        # f1, f2, f3 = get_bad_fintess_values(machine_type='PMSM')
                        # f1, f2, f3 = get_bad_fintess_values(machine_type='CPPM')
                        f1, f2, f3 = ad.get_bad_fintess_values(machine_type='CSPPM')
                    logger.debug('f1,f2,f3: %s, %s, %s', f1, f2, f3)

                break

        ad.counter_fitness_return += 1
        logger.debug('Fitness: %d, %d', ad.counter_fitness_called, ad.counter_fitness_return)
        return [f1, f2, f3]

    # Return number of objectives
    def get_nobj(self):
        global ad
        return ad.nobj

    # Return bounds of decision variables (a.k.a. chromosome)
    def get_bounds(self):
        global ad
        bounds_denorm = list(ad.get_free_variable_bounds_dict().values())
        min_b, max_b = np.asarray(bounds_denorm).T 
        return ( min_b.tolist(), max_b.tolist() )

    # Return function name
    def get_name(self):
        return "Bearingless PMSM Design"

import pygmo as pg
# algorithm = pg.algorithm(pg.sade(gen=100))
# pop = pg.population(prob, size=50)
# pop = algorithm.evolve(pop)
# best_fitness = pop.champion_f
# best_solution = pop.champion_x
# print("Best Fitness:", best_fitness)
# print("Best Solution:", best_solution)
def get_prob():
    udp = Problem_BearinglessSynchronousDesign()
    prob = pg.problem(udp)
    return udp, prob #, popsize
