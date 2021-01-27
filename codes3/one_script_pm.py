# coding:u8
import utility
from utility import my_execfile
import utility_moo
from win32com.client import pywintypes
bool_post_processing = False # solve or post-processing
bool_re_evaluate = False

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# 0. FEA Setting / General Information & Packages Loading
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# FEA setting
my_execfile('./default_setting.py', g=globals(), l=locals())
fea_config_dict
fea_config_dict['bool_post_processing'] = bool_post_processing
if True:
    # PEMD2020

    # Combined winding PMSM
    fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.95
    fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.05
        # run_folder = r'run#610/' # FRW constraint is added and sleeve_length is 3 (not varying). Excitation ratio is 95%:5% between Torque and Suspension windings.
        # run_folder = r'run#611/' # zero Rs is not allowed
        # run_folder = r'run#61495/' # spec_PEMD_BPMSM_Q12p2, 99 zQ is fixed to 10 | 98 zQ is derived | 97 sleeve length is reduced to 1 mm | 96 Jingwei's layout | 95 alpha_rm is fixed to be 360/2/p | 94 full alpha_rm bug is fixed | )

    # run_folder = r'run#62399/' # spec_ECCE_PMSM_ (Q6p2)
    # run_folder = r'run#62499/' # spec_PEMD_BPMSM_Q12p2
    # run_folder = r'run#62599/' # spec_PEMD_BPMSM_Q6p1)
    # run_folder = r'run#62699/' # spec_PEMD_BPMSM_Q12p4)
    # run_folder = r'run#62799/' # spec_PEMD_BPMSM_Q24p1
    run_folder = r'run#63899/' # spec_PEMD_BPMSM_Q12p1 (run by Ashad)

    # spec's
    # my_execfile('./spec_ECCE_PMSM_Q6p2.py', g=globals(), l=locals()) # Q=6, p=2
    # my_execfile('./spec_PEMD_BPMSM_Q12p2.py', g=globals(), l=locals()) # Q=12, p=2
    # my_execfile('./spec_PEMD_BPMSM_Q6p1.py', g=globals(), l=locals()) # Q=6, p=1
    # my_execfile('./spec_PEMD_BPMSM_Q12p4.py', g=globals(), l=locals()) # Q=12, p=4, ps=5
    # my_execfile('./spec_PEMD_BPMSM_Q24p1.py', g=globals(), l=locals()) # Q=24, p=1, ps=2
    my_execfile('./spec_PEMD_BPMSM_Q12p1.py', g=globals(), l=locals()) # Q=12, p=1, ps=2

    fea_config_dict['run_folder'] = run_folder

    # Adopt Bianchi 2006 for a SPM motor template
    spec.build_pmsm_template(fea_config_dict, im_template=None)
else:
    # TIA_ITEC2018_PAPER

    if 'Y730' in fea_config_dict['pc_name']:
        ################################################################
        # Y730
        ################################################################
        # run_folder = r'run#600/' # FRW constraint is removed and sleeve_length is 3 (not varying)
        # run_folder = r'run#601/' # FRW constraint is removed and sleeve_length is 2.5 (not varying)

        # Combined winding PMSM
        fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.95
        fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.05
        run_folder = r'run#603/' # FRW constraint is added and sleeve_length is 3 (not varying). Excitation ratio is 95%:5% between Torque and Suspension windings.

        # # Separate winding PMSM
        # fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.60
        # fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.05
        # run_folder = r'run#604/'
        # raise
    elif 'Severson01' in fea_config_dict['pc_name']:
        ################################################################
        # Severson01
        ################################################################
        print('Severson01')
        # Separate winding PMSM
        fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.60
        fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.05
        run_folder = r'run#604010/'
    elif 'Severson02' in fea_config_dict['pc_name']:
        ################################################################
        # Severson02
        ################################################################
        print('Severson02')
        # Combined winding PMSM
        fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.95
        fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.05
        run_folder = r'run#603020/'
    else:
        ################################################################
        # T440p
        ################################################################
        print('T440p')
        # Combined winding PMSM
        fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.95
        fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.05
        run_folder = r'run#603010/' # continued run from severson01

    fea_config_dict['run_folder'] = run_folder

    # spec's
    my_execfile('./spec_TIA_ITEC_.py', g=globals(), l=locals()) # Q=24, p=1

    # Case Q=24 can use IM's stator for PMSM's
    spec.build_im_template(fea_config_dict)
    spec.build_pmsm_template(fea_config_dict, im_template=spec.im_template)

# select motor type here
spec.acm_template = spec.pmsm_template
print('Build ACM template...')

import acm_designer
global ad
ad = acm_designer.acm_designer(fea_config_dict, spec)
ad.init_logger(prefix='bpmsm')
ad.bool_re_evaluate = bool_re_evaluate

ad.bounds_denorm = spec.acm_template.get_classic_bounds(which_filter='VariableSleeveLength')
# ad.bounds_denorm = spec.acm_template.get_classic_bounds(which_filter='FixedSleeveLength') # ad.get_classic_bounds() <- obsolete
ad.bound_filter  = spec.acm_template.bound_filter
# print(ad.bounds_denorm)
# quit()
print('---------------------\nBounds:')
idx_ad = 0
for idx, f in enumerate(ad.bound_filter):
    if f == True:
        print(idx, f, '[%g,%g]'%tuple(spec.acm_template.original_template_neighbor_bounds[idx]), '[%g,%g]'%tuple(ad.bounds_denorm[idx_ad]))
        idx_ad += 1
    else:
        print(idx, f, '[%g,%g]'%tuple(spec.acm_template.original_template_neighbor_bounds[idx]))
quit()

if ad.bool_re_evaluate:
    ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + ad.solver.fea_config_dict['run_folder']
    number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
    print('Count of chromosomes:', len(ad.solver.swarm_data))

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Optimization
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
from acm_designer import get_bad_fintess_values
import pygmo as pg
ad.counter_fitness_called = 0
ad.counter_fitness_return = 0
__builtins__.ad = ad # share global variable between modules # https://stackoverflow.com/questions/142545/how-to-make-a-cross-module-variable
print(__builtins__.ad)
from Problem_BearinglessSynchronousDesign import Problem_BearinglessSynchronousDesign

if bool_post_processing == True:
    import one_script_pm_post_processing 
    one_script_pm_post_processing.post_processing(ad, fea_config_dict)
    quit()

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Multi-Objective Optimization
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
if True:
################################################################
# MOO Step 1:
#   Create UserDefinedProblem and create population
#   The magic method __init__ cannot be fined for UDP class
################################################################
    udp = Problem_BearinglessSynchronousDesign()
    prob = pg.problem(udp)

    popsize = 78
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
        ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome
        print('ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome = %d'%(number_of_chromosome))

    else:
        print('Nothing exists in swarm_data.txt.\nThis is a whole new run.')
        ad.flag_do_not_evaluate_when_init_pop = False
        number_of_finished_chromosome_in_current_generation = None
        number_of_finished_iterations = 0 # 实际上跑起来它不是零，而是一，因为我们认为初始化的一代也是一代。或者，我们定义number_of_finished_iterations = number_of_chromosome // popsize

    if bool_re_evaluate:
        ad.counter_fitness_called = ad.counter_fitness_return = 0        

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
                    print('Set "ad.flag_do_not_evaluate_when_init_pop" to False...')
                    ad.flag_do_not_evaluate_when_init_pop = False
                    print('Calling pop.set_x()---this is a restart for individual#%d during pop initialization.'%(i))
                    print(i, 'get_fevals:', prob.get_fevals())
                    pop.set_x(i, pop_array[i]) # evaluate this guy
        else:
            # 新办法，直接从swarm_data.txt（相当于archive）中判断出当前最棒的群体
            swarm_data_on_pareto_front = utility_moo.learn_about_the_archive(prob, ad.solver.swarm_data, popsize, fea_config_dict)
            # print(swarm_data_on_pareto_front)
            for i in range(popsize):
                pop.set_xf(i, swarm_data_on_pareto_front[i][:-3], swarm_data_on_pareto_front[i][-3:])

        # 必须放到这个if的最后，因为在 learn_about_the_archive 中是有初始化一个 pop_archive 的，会调用fitness方法。
        ad.flag_do_not_evaluate_when_init_pop = False

    print('-'*40, '\nPop is initialized:\n', pop)
    hv = pg.hypervolume(pop)
    quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
    print('quality_measure: %g'%(quality_measure))
    # raise KeyboardInterrupt

    # 初始化以后，pop.problem.get_fevals()就是popsize，但是如果大于popsize，说明“pop.set_x(i, pop_array[i]) # evaluate this guy”被调用了，说明还没输出过 survivors 数据，那么就写一下。
    if pop.problem.get_fevals() > popsize:
        print('Write survivors.')
        ad.solver.write_swarm_survivor(pop, ad.counter_fitness_return)

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
    number_of_iterations = 250
    logger = logging.getLogger(__name__)
    # try:
    if True:
        for _ in range(number_of_finished_iterations, number_of_iterations):
            msg = 'This is iteration #%d. '%(_)
            print(msg)
            logger.info(msg)
            pop = algo.evolve(pop)

            msg += 'Write survivors to file. '
            ad.solver.write_swarm_survivor(pop, ad.counter_fitness_return)

            hv = pg.hypervolume(pop)
            quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
            msg += 'Quality measure by hyper-volume: %g'% (quality_measure)
            print(msg)
            logger.info(msg)
            
            utility_moo.my_print(ad, pop, _)
            # my_plot(fits, vectors, ndf)
    # except Exception as e:
    #     print(pop.get_x())
    #     print(pop.get_f().tolist())
    #     raise e


