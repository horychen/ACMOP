def main(bool_post_processing=False):
    # Vernier Machine
    # mop = AC_Machine_Optiomization_Wrapper(
    #         select_fea_config_dict = "#03 JMAG Non-Nearingless Motor Evaluation Setting",
    #         select_spec            = "PMVM p2pr10-Q12y3 Wenbo",
    #         project_loc            = r'D:/DrH/acmop/_WenboVShapeVernier/'
    #     )

    # TEC-ISMB-2021
        # select_spec = "IM Q24p1y9 A"
        # select_spec = "IM Q24p1y9 Qr32"
        # select_spec = "IM Q24p2y6 Qr32"
        # select_spec = "IM Q24p2y6 Qr16"
        # select_spec = "IM Q36p3y5ps2 Qr20-FSW Round Bar"
        # select_spec = "IM Q36p3y5ps2 Qr24-ISW Round Bar"
        # select_spec = "IM Q36p3y5ps2 Qr20-FSW Round Bar Separate Winding"
        # select_spec = "IM Q36p3y5ps2 Qr24-ISW Round Bar Separate Winding"
        # select_spec = "IM Q24p1y9 Qr14 Round Bar"
        # select_spec = "IM Q24p1y9 Qr16 Round Bar"
        # select_spec = "IM p2ps3Qs18y4 Qr30-FSW Round Bar EquivDoubleLayer"
        # select_spec = "IM p2ps3Qs24y5 Qr18 Round Bar EquivDoubleLayer"
    # TIA-IEMDC-ECCE-2020
        # select_spec = "IM Q24p1y9 A"
        # select_spec = "IM Q24p1y9 Qr32"
        # select_spec = "IM Q24p2y6 Qr32"
        # select_spec = "IM Q24p2y6 Qr16"
        # select_spec = "IM Q24p1y9 Qr32 Round Bar" # RevisionInNov: rev2 reviewer 1
        # select_spec = "IM Q24p1y9 Qr16 Round Bar" # Prototype
    # mop = AC_Machine_Optiomization_Wrapper(
    #         select_fea_config_dict = "#019 JMAG IM Nine Variables",
    #         select_spec            = 'IM p2ps3Qs18y4 Qr30-FSW Round Bar EquivDoubleLayer',
    #         project_loc            = r'D:/DrH/acmop/_default/'
    #     )

    # PEMD 2020 paper 
        # select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting'
        # select_fea_config_dict = '#021   PMSM Re-evaluation wo/ CSV Setting'
        # select_fea_config_dict = "#0210 JMAG PMSM Re-evaluation wo/ CSV Setting Q12p4ps5 Sub-hamonics"
        # select_spec =  "PMSM Q06p1y2 A"            # [ './spec_PEMD_BPMSM_Q6p1.py',
        # select_spec =  "PMSM Q06p2y1 A"            #   './spec_ECCE_PMSM_Q6p2.py',
        # select_spec =  "PMSM Q12p1y5 A"            #   './spec_PEMD_BPMSM_Q12p1.py',
        # select_spec =  "PMSM Q12p2y3 A"            #   './spec_PEMD_BPMSM_Q12p2.py',
        # select_spec =  "PMSM Q12p4y1 A"            #   './spec_PEMD_BPMSM_Q12p4.py',
        # select_spec =  "PMSM Q24p1y9 A"            #   './spec_PEMD_BPMSM_Q24p1.py'],
    mop = AC_Machine_Optiomization_Wrapper(
            select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting',
            select_spec            = 'PMSM Q12p2y3 A',
            project_loc            = r'D:/DrH/acmop/_default/'
        )


    #########################
    # Call the five modules
    #########################

    # mop.part_winding()
    acm_template = mop.part_initialDesign()

    if not bool_post_processing:
        # mop.part_evaluation()
        mop.part_optimization(acm_template)
        # mop.part_reportWithStreamlit()
    else:
        # Recover a design from its jsonpickle-object file
        import utility_json
        try:
            motor_design_variant = utility_json.from_json_recursively('p2ps1-Q12y3-0999', load_here=mop.ad.output_dir+'jsonpickle/')
            motor_design_variant.build_jmag_project(motor_design_variant.project_meta_data)
        except FileNotFoundError as e:
            print(e)
            print("你有没有搞错啊？json文件找不到啊，忘记把bool_post_processing改回来了？")

from dataclasses import dataclass
import sys; sys.path.insert(0, './codes3/')
import main_utility

# part_winding
# import winding_layout, PyX_Utility, math

# part_initialDesign
import pyrhonen_procedure_as_function, acm_designer, bearingless_spmsm_design
global ad

@dataclass
class AC_Machine_Optiomization_Wrapper(object):
    ''' Inputs
    '''
    # A. select FEA setting
    select_fea_config_dict: str
    # B. select design specification
    select_spec: str
    # C. decide output directory 20210127
    project_loc: str

    ''' Derived
    '''
    output_dir: str = None
    spec_input_dict: dict = None
    fea_config_dict: dict = None

    def __post_init__(self):
        Help = r'''[Steps for adding a new slot pole combination for IM]
        1. Update machine_specifications.json
        2. Run winding_layout_derivation_ismb2020.py to get a new stator winding layout and paste the code into winding_layout.py
        3. Run Pole-specific_winding_with_neutral_plate_the_design_table_generator.py to get a new rotor winding layout and paste the code into winding_layout.py
        4. Update this file with new "select_spec".
        '''
        self.output_dir, self.spec_input_dict, self.fea_config_dict = main_utility.load_settings(self.select_spec, self.select_fea_config_dict, self.project_loc)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[1] Winding Part (Can be skipped)'
    # Use winding_layout_derivation.py to derive windings defined in class winding_layout_v2 during choosing winding phase.
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_winding(self):
        # wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=False, Qs=24, p=2, ps=1)
        # wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=2, ps=1, coil_pitch_y=6)
        # wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=1, ps=2, coil_pitch_y=9)
        wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=self.spec_input_dict['DPNV_or_SEPA'], 
                                                Qs=self.spec_input_dict['Qs'], 
                                                p=self.spec_input_dict['p'], 
                                                ps=self.spec_input_dict['ps'], 
                                                coil_pitch_y=self.spec_input_dict['coil_pitch_y'])

        '[1.1] Winding in the slot'
        if False:
            def draw_winding_in_the_slot(u, Qs, list_layer_phases, list_layer_signs, text=''):

                for i in range(Qs):
                    radius_slot = 30
                    LRIF = layer_radius_incremental_factor = 0.1
                    angular_loc = 2*math.pi/Qs*i
                    x_slot = radius_slot*math.cos(angular_loc)
                    y_slot = radius_slot*math.sin(angular_loc)
                    u.pyx_text(   [ x_slot*(1.0+LRIF), 
                                    y_slot*(1.0+LRIF)],
                                  str(i+1) )
                    u.pyx_marker( [ x_slot*(1.0+2*LRIF), 
                                    y_slot*(1.0+2*LRIF)], size=0.05)

                    radius_tooth = radius_slot + 5 
                    x_tooth = radius_tooth*math.cos(angular_loc + math.pi/Qs)
                    y_tooth = radius_tooth*math.sin(angular_loc + math.pi/Qs)
                    radius_airgap = radius_slot - 5
                    x_toothtip = radius_airgap*math.cos(angular_loc+math.pi/Qs)
                    y_toothtip = radius_airgap*math.sin(angular_loc+math.pi/Qs)
                    u.pyx_line([x_toothtip, y_toothtip], [x_tooth, y_tooth])

                    for ind, phases in enumerate(list_layer_phases):
                        signs = list_layer_signs[ind]
                        u.pyx_text(   [ x_slot*(1.0-ind*LRIF), 
                                        y_slot*(1.0-ind*LRIF)],
                                      '$' + phases[i].lower() + '^' + signs[i]
                                      + '$' )

                u.pyx_text([0,0], (r'DPNV Winding' if wily.DPNV_or_SEPA else r'Separate Winding') + text)

            u = PyX_Utility.PyX_Utility()
            draw_winding_in_the_slot(u, wily.Qs, wily.list_layer_motor_phases, wily.list_layer_motor_signs, text=' Motor Mode' )
            u.cvs.writePDFfile(mop.output_dir + 'pyx_output_M')

            u = PyX_Utility.PyX_Utility()
            draw_winding_in_the_slot(u, wily.Qs, wily.list_layer_suspension_phases, wily.list_layer_suspension_signs, text=' Suspension Mode' )
            u.cvs.writePDFfile(mop.output_dir + 'pyx_output_S')
            # u.cvs.writeSVGfile(r'C:\Users\horyc\Desktop\pyx_output')
            # u.cvs.writeEPSfile(r'C:\Users\horyc\Desktop\pyx_output')
            # quit()

        '[1.2] Winding function / Current Linkage waveform'
        if False:
            from pylab import plt, np
            zQ = 100 # number of conductors/turns per slot (Assume to be 100 for now)
            turns_per_layer = zQ / wily.number_winding_layer
            U_phase = winding_layout.PhaseWinding(wily.Qs, wily.m, turns_per_layer, wily.ox_distribution_phase_U)
            U_phase.plotFuncObj(U_phase.winding_func)
            U_phase.fig_plotFuncObj.savefig(mop.output_dir + 'winding_function.png')
            U_phase.plot2piFft(U_phase.winding_func, Fs=1/(2*np.pi/3600), L=32000*2**4) # 在2pi的周期内取360个点
            U_phase.fig_plot2piFft.savefig(mop.output_dir + 'winding_function_DFT.png')
            plt.show()

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[2] Initial Design Part'
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_initialDesign(self):
        if 'PM' in self.select_spec:
            function = bearingless_spmsm_design.bearingless_spmsm_template
        elif 'VM' in self.select_spec:
            function = vernier_motor_design.vernier_motor_VShapePM_template
        acm_template = function(self.fea_config_dict, self.spec_input_dict)

        self.ad = acm_designer.acm_designer(
                    self.fea_config_dict, 
                    self.spec_input_dict, 
                    self.output_dir, 
                    self.select_spec, 
                    self.select_fea_config_dict,
                    acm_template=acm_template
                )

        if False:
            if 'Y730' in self.fea_config_dict['pc_name']:
                self.ad.build_oneReport() # require LaTeX
                # ad.talk_to_mysql_database() # require MySQL

        return acm_template

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[3] Evaluation Part (Can be skipped)'
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_evaluation(self):
        # build x_denorm for the template design
        x_denorm = self.ad.acm_template.build_x_denorm()

        # evaluate design (with json output)
        cost_function, f1, f2, f3, FRW, \
            normalized_torque_ripple, \
            normalized_force_error_magnitude, \
            force_error_angle = self.ad.evaluate_design_json_wrapper(self.ad.acm_template, x_denorm)

        print('[part_evaluation]:', cost_function, f1, f2, f3, FRW, \
        normalized_torque_ripple, \
        normalized_force_error_magnitude, \
        force_error_angle)

        print('Check several things: 1. the winding initial excitation angle; 2. the rotor d-axis initial position should be orthoganal to winding excitation field.')

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[4] Optimization Part' Multi-Objective Optimization
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_optimization(self, acm_template):
        ad = self.ad
        ad.init_logger(prefix='acmdm')

        # [4.1] Get bounds
        # ad.bounds_denorm = acm_template.bounds_denorm # get_classic_bounds(which_filter=self.fea_config_dict['which_filter'])
        # ad.bound_filter  = acm_template.bound_filter
        # otnb = acm_template.original_template_neighbor_bounds
        # print('---------------------\nBounds: (if there are two bounds within one line, they should be the same)')
        # idx_ad = 0
        # for idx, f in enumerate(ad.bound_filter):
        #     if f == True:
        #         print(idx, f, '[%g,%g]'%tuple(otnb[idx]), '[%g,%g]'%tuple(ad.bounds_denorm[idx_ad]))
        #         idx_ad += 1
        #     else:
        #         print(idx, f, '[%g,%g]'%tuple(otnb[idx]))

        # if self.fea_config_dict['bool_post_processing'] == True: # use the new script file instead: main_post_processing_pm.py
        #     import one_script_pm_post_processing 
        #     one_script_pm_post_processing.post_processing(ad, self.fea_config_dict)
        #     quit()

        # [4.3] MOO
        from acm_designer import get_bad_fintess_values
        import logging
        import utility_moo
        import pygmo as pg
        ad.counter_fitness_called = 0
        ad.counter_fitness_return = 0
        __builtins__.ad = ad # share global variable between modules # https://stackoverflow.com/questions/142545/how-to-make-a-cross-module-variable
        import codes3.Problem_BearinglessSynchronousDesign # must import this after __builtins__.ad = ad
        print('[acmop.py]', __builtins__.ad)

        ################################################################
        # MOO Step 1:
        #   Create UserDefinedProblem and create population
        #   The magic method __init__ cannot be fined for UDP class
        ################################################################
        # [4.3.1] Basic setup
        _, prob, popsize = codes3.Problem_BearinglessSynchronousDesign.get_prob_and_popsize()

        print('[acmop.py]', '-'*40 + '\n[acmop.py] Pop size is', popsize)

        # [4.3.2] Generate the pop
        if False:
            pop = pg.population(prob, size=popsize) 
        # Add Restarting Feature when generating pop
        else:
            from main_utility import get_sorted_swarm_data_from_the_archive
            # def get_sorted_swarm_data_from_the_archive(path_to_archive):
            #     output_dir_backup = ad.solver.output_dir
            #     ad.solver.output_dir = ad.solver.fea_config_dict['dir.parent'] + path_to_archive
            #     number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
            #     ad.solver.output_dir = output_dir_backup

            #     ad.flag_do_not_evaluate_when_init_pop = True
            #     pop = pg.population(prob, size=popsize)
            #     swarm_data_on_pareto_front = utility_moo.learn_about_the_archive(prob, ad.solver.swarm_data, popsize, ad.solver.fea_config_dict, bool_plot_and_show=False)
            #     ad.flag_do_not_evaluate_when_init_pop = False
            #     return swarm_data_on_pareto_front

            # 检查swarm_data.txt，如果有至少一个数据，返回就不是None。
            print('[acmop.py] Check swarm_data.txt...')
            number_of_chromosome = ad.solver.read_swarm_data()

            # case 1: swarm_data.txt exists
            if number_of_chromosome is not None:

                number_of_finished_iterations                       = number_of_chromosome // popsize
                number_of_finished_chromosome_in_current_generation = number_of_chromosome % popsize

                # 如果刚好整除，把余数0改为popsize
                if number_of_finished_chromosome_in_current_generation == 0:
                    number_of_finished_chromosome_in_current_generation = popsize
                    print(f'\tThere are {number_of_chromosome} chromosomes found in swarm_data.txt.')
                    print('\tWhat is the odds! The script just stopped when the evaluation of the whole pop is finished.')
                    print('\tSet number_of_finished_chromosome_in_current_generation to popsize %d'%(number_of_finished_chromosome_in_current_generation))

                print('[acmop.py] This is a restart of '+ self.fea_config_dict['run_folder'][:-1])
                print('\tNumber of finished iterations is %d'%(number_of_finished_iterations))
                # print('This means the initialization of the population class is interrupted. So the pop in swarm_data.txt is used as the survivor.')

                # swarm_survivor feature. Not sure if this is needed anymore...
                if True:
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
                print('[acmop.py] ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome = %d'%(number_of_chromosome))

                # case 1-A: swarm_data.txt exists and this is a re-evaluation run using the existing csv files (比如我们修改了计算铜损的代码，那就必须借助已有的有限元结果重新生成swarm_data.txt)
                if self.fea_config_dict['bool_re_evaluate']:
                    ad.counter_fitness_called = ad.counter_fitness_return = 0        

                # 禁止在初始化pop时运行有限元
                ad.flag_do_not_evaluate_when_init_pop = True
                # 初始化population，如果ad.flag_do_not_evaluate_when_init_pop是False，那么就说明是 new run，否则，整代个体的fitness都是[0,0,0]。
                pop = pg.population(prob, size=popsize)
                if self.fea_config_dict['bool_re_evaluate_wo_csv']:
                    swarm_data_backup = ad.solver.swarm_data[::] # This is going to be over-written in next line
                    swarm_data_on_pareto_front, _ = get_sorted_swarm_data_from_the_archive(prob, popsize, path_to_archive)
                    ad.flag_do_not_evaluate_when_init_pop = True # When you call function get_sorted_swarm_data_from_the_archive, flag_do_not_evaluate_when_init_pop is set to False at the end. Sometimes we do not want this, for example, restarting restart re-evaluation without csv.
                    ad.solver.swarm_data = swarm_data_backup
                    for i in range(popsize):
                        # print(path_to_archive, ':', swarm_data_on_pareto_front[i][::-1])
                        pop.set_xf(i, swarm_data_on_pareto_front[i][:-3], swarm_data_on_pareto_front[i][-3:])
                    print('Old pop:')
                    print(pop)

                # Restarting feature related codes
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

            # case 2: swarm_data.txt does not exist
            else:
                number_of_finished_chromosome_in_current_generation = None
                number_of_finished_iterations = 0 # 实际上跑起来它不是零，而是一，因为我们认为初始化的一代也是一代。或者，我们定义number_of_finished_iterations = number_of_chromosome // popsize

                # case 2-A: swarm_data.txt does not exist and this is a whole new run.
                if not self.fea_config_dict['bool_re_evaluate_wo_csv']:
                    print('Nothing exists in swarm_data.txt.\nThis is a whole new run.')
                    ad.flag_do_not_evaluate_when_init_pop = False
                    pop = pg.population(prob, size=popsize)

                # case 2-B: swarm_data.txt does not exist and this is a re-evalation run (without csv)
                else:
                    print('Nothing exists in swarm_data.txt.\nRe-start from %s'%(path_to_archive))
                    ad.flag_do_not_evaluate_when_init_pop = True
                    pop = pg.population(prob, size=popsize)
                    # read in swarm data from another older run's archive and start from it!
                    swarm_data_on_pareto_front, _ = get_sorted_swarm_data_from_the_archive(prob, popsize, path_to_archive)
                    ad.flag_do_not_evaluate_when_init_pop = False
                    for i in range(popsize):
                        print(path_to_archive, ':', swarm_data_on_pareto_front[i][::-1])
                        pop.set_x(i, swarm_data_on_pareto_front[i][:-3]) # re-evaluate this guy

            # this flag must be false to move on
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
        # [4.3.3] Selecting algorithm
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
        # [4.3.4] Begin optimization
        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
        number_of_finished_iterations = number_of_chromosome // popsize
        number_of_iterations = 20
        logger = logging.getLogger(__name__)
        # try:
        if True:
            for _ in range(number_of_finished_iterations, number_of_iterations):
                ad.number_of
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

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[5] Report Part'
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_reportWithStreamlit(self):
        # bool_post_processing
        # best desigh?

        # Status report: Generation, individuals, geometry as input, fea tools, performance as output (based on JSON files)
        # Do not show every step, but only those key steps showing how this population is built 

        # 从json重构JMAG模型（x_denorm）
        pass

if __name__ == '__main__':
    main()
