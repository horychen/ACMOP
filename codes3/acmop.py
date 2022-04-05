# Please use shortcut "ctrl+k,ctrl+1" to fold the code for better navigation
# Please use shortcut "ctrl+k,ctrl+2" to fold the code for better navigation

import os, json, acm_designer, bearingless_spmsm_design, vernier_motor_design, bearingless_induction_design, flux_alternator_design, flux_switching_pm_design

from soupsieve import select
# from codes3.population import VanGogh_JMAG

import VanGogh_Cairo

from matplotlib import projections # for part_initialDesign

from dataclasses import dataclass
@dataclass
class AC_Machine_Optiomization_Wrapper(object):
    ''' Inputs
    '''
    # A. select FEA setting
    select_fea_config_dict: str
    # B. select design specification
    select_spec: str
    # C. decide output directory (initialize either one)
    project_loc: str = None
    path2SwarmData: str = None
    # D. this is up to you
    bool_show_GUI: bool = False

    ''' Derived
    '''
    spec_input_dict: dict = None
    fea_config_dict: dict = None

    def __post_init__(self):
        self.Help = r'''[Steps for adding a new slot pole combination for IM]
        1. Update machine_specifications.json
        2. Run winding_layout_derivation_ismb2020.py to get a new stator winding layout and paste the code into winding_layout.py
        3. Run Pole-specific_winding_with_neutral_plate_the_design_table_generator.py to get a new rotor winding layout and paste the code into winding_layout.py
        4. Update this file with new "select_spec".
        '''
        self.spec_input_dict, self.fea_config_dict = self.load_settings( 
                                            self.select_spec, 
                                            self.select_fea_config_dict, 
                                            project_loc=self.project_loc, 
                                            path2SwarmData=self.path2SwarmData)
        self.fea_config_dict['designer.Show'] = self.bool_show_GUI

        print('[acmop.py] project_loc (user-input):', self.project_loc)
        if self.path2SwarmData is None:
            self.path2SwarmData = self.project_loc + self.select_spec.replace(' ', '_') + '/'
        if self.project_loc is None:
            self.project_loc = os.path.abspath(os.path.join(self.path2SwarmData, '..',))

        # Convert to abs path (JMAG requires absolute path)
        self.project_loc                   = os.path.abspath(self.project_loc) + '/'
        self.path2SwarmData                = os.path.abspath(self.path2SwarmData) + '/'
        self.fea_config_dict['output_dir'] = os.path.abspath(self.fea_config_dict['output_dir']) + '/'
        print('[acmop.py] project_loc (converted) :', self.project_loc)
        print('[acmop.py] path2SwarmData          :', self.path2SwarmData)
        print('[acmop.py] output_dir              :', self.fea_config_dict['output_dir'])

        self.acm_template = self.part_initialDesign() # Module 2 (mop.ad is available now)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[1] Winding Part (Can be skipped)'
    # Use winding_layout_derivation.py to derive windings defined in class winding_layout_v2 during choosing winding phase.
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_winding(self):
        import winding_layout, PyX_Utility, math # for part_winding

        # wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=False, Qs=24, p=2, ps=1)
        # wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=2, ps=1, coil_pitch_y=6)
        # wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=1, ps=2, coil_pitch_y=9)
        wily = winding_layout.winding_layout_v2(DPNV_or_SEPA=self.spec_input_dict['DPNV_or_SEPA'], 
                                                Qs=self.spec_input_dict['Qs'], 
                                                p=self.spec_input_dict['p'], 
                                                ps=self.spec_input_dict['ps'], 
                                                coil_pitch_y=self.spec_input_dict['coil_pitch_y'])

        '[1.1] Winding in the slot'
        if True:
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
        if True:
            from pylab import plt, np
            zQ = 100 # number of conductors/turns per slot (Assume to be 100 for now)
            turns_per_layer = zQ / wily.number_winding_layer
            U_phase = winding_layout.PhaseWinding(wily.Qs, wily.m, turns_per_layer, wily.ox_distribution_phase_U)
            U_phase.plotFuncObj(U_phase.winding_func)
            U_phase.fig_plotFuncObj.savefig(mop.output_dir + 'winding_function.png')
            U_phase.plot2piFft(U_phase.winding_func, Fs=1/(2*np.pi/3600), L=32000*2**4) # 在2pi的周期内取360个点
            U_phase.fig_plot2piFft.savefig(mop.output_dir + 'winding_function_·.png')
            plt.show()

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[2] Initial Design Part'
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_initialDesign(self):
        if 'PMSM' in self.select_spec:
            function = bearingless_spmsm_design.bearingless_spmsm_template
        elif 'PMVM' in self.select_spec:
            function = vernier_motor_design.vernier_motor_VShapePM_template
        elif 'IM' in self.select_spec:
            function = bearingless_induction_design.bearingless_induction_template
        elif 'Flux Alternator' in self.select_spec: 
            function = flux_alternator_design.flux_alternator_template
        elif 'FSPM' in self.select_spec: 
            function = flux_switching_pm_design.FSPM_template
        acm_template = function(self.fea_config_dict, self.spec_input_dict)

        self.ad = acm_designer.acm_designer(
                    self.select_spec, 
                    self.spec_input_dict, 
                    self.select_fea_config_dict,
                    self.fea_config_dict, 
                    acm_template=acm_template,
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
        print('[acmop.py] x_denorm:',  x_denorm)
        print('[acmop.py] x_denorm_dict:', self.ad.acm_template.x_denorm_dict)

        if True:
            ''' Default transient FEA
            '''
            motor_design_variant = self.ad.evaluate_design_json_wrapper(self.ad.acm_template, x_denorm, counter='Initial')

            from pylab import plt, np
            fig, axes = plt.subplots(4)
            axes[0].plot(motor_design_variant.analyzer.motor_current_U)
            axes[0].plot(motor_design_variant.analyzer.motor_current_V)
            axes[0].plot(motor_design_variant.analyzer.motor_current_W)
            axes[1].plot(motor_design_variant.analyzer.bearing_current_U)
            axes[1].plot(motor_design_variant.analyzer.bearing_current_V)
            axes[1].plot(motor_design_variant.analyzer.bearing_current_W)
            axes[2].plot(motor_design_variant.analyzer.femm_motor_currents_d)
            axes[2].plot(motor_design_variant.analyzer.femm_motor_currents_q)
            axes[3].plot(motor_design_variant.analyzer.femm_motor_fluxLinkage_d)
            axes[3].plot(motor_design_variant.analyzer.femm_motor_fluxLinkage_q)

            fig, axes = plt.subplots(8)
            axes[0].plot(motor_design_variant.analyzer.femm_time, motor_design_variant.analyzer.femm_torque)
            axes[1].plot(motor_design_variant.analyzer.femm_time, motor_design_variant.analyzer.sfv.force_abs)
            axes[2].plot(motor_design_variant.analyzer.femm_time, motor_design_variant.analyzer.sfv.force_x)
            axes[3].plot(motor_design_variant.analyzer.femm_time, motor_design_variant.analyzer.sfv.force_y)
            axes[4].plot(motor_design_variant.analyzer.femm_time, motor_design_variant.analyzer.femm_energy)
            axes[5].plot(motor_design_variant.analyzer.femm_time, list(map(lambda el: el[0], motor_design_variant.analyzer.femm_circuit_currents)))
            axes[6].plot(motor_design_variant.analyzer.femm_time, list(map(lambda el: el[0], motor_design_variant.analyzer.femm_circuit_voltages)))
            axes[7].plot(motor_design_variant.analyzer.femm_time, list(map(lambda el: el[0], motor_design_variant.analyzer.femm_circuit_fluxLinkages)))
            plt.show()


            print('[acmop.py] Listing analyzer.spec_performance_dict:')
            for k,v in motor_design_variant.analyzer.spec_performance_dict.items():
                print('\t', k, v)

        else:

            ''' An example showing how to change the initial angle between between rotor d-axis and current vector.
                At t=0, current vector is aligned with beta-axis, so we need to align rotor d-axis with the alpha-axis, such that id=0 control is implemented.
                However, salient pole motor can produce more torque if we apply some id. 
                In this case, we can sweep variable self.ad.acm_template.fea_config_dict['femm.MechDeg_IdEqualToNonZeroAngle'] to reach maximum torque.
            '''
            from pylab import plt, np
            from utility import suspension_force_vector
            fig, axes = plt.subplots(5)
            for angle in np.arange(-15, 15.1, 5):
                self.ad.acm_template.fea_config_dict['femm.MechDeg_IdEqualToNonZeroAngle'] = angle
                print('User shifts the initial rotor position angle by', self.ad.acm_template.fea_config_dict['femm.MechDeg_IdEqualToNonZeroAngle'], 'deg')
                motor_design_variant = self.ad.evaluate_design_json_wrapper(self.ad.acm_template, x_denorm)

                force_x = list(map(lambda el: el[0], motor_design_variant.analyzer.femm_forces))
                force_y = list(map(lambda el: el[1], motor_design_variant.analyzer.femm_forces))
                sfv = suspension_force_vector(force_x, force_y)

                axes[0].plot(motor_design_variant.analyzer.femm_time, motor_design_variant.analyzer.femm_torque, label=str(angle))
                axes[1].plot(motor_design_variant.analyzer.femm_time, sfv.force_abs, label=str(angle))
                axes[2].plot(motor_design_variant.analyzer.femm_time, force_x, label=str(angle))
                axes[3].plot(motor_design_variant.analyzer.femm_time, force_y, label=str(angle))
                axes[4].plot(motor_design_variant.analyzer.femm_time, motor_design_variant.analyzer.femm_energy, label=str(angle))
            for ax in axes:
                ax.legend()
            plt.show()

        print('[acmop.py] Check several things: 1. the winding initial excitation angle; 2. the rotor d-axis initial position should be orthoganal to winding excitation field.')

    def part_evaluation_geometry(self, xf=[], counter='Cairo'):
        print('\n ------- part_evaluation_geometry -------')
        x_denorm = self.ad.acm_template.build_x_denorm()
        print('--------xf is ', xf)
        if xf != []:
            x_denorm = xf[:len(x_denorm)]
        acm_variant = self.ad.build_acm_variant(self.ad.acm_template, x_denorm, counter=counter)
        toolCairo = VanGogh_Cairo.VanGogh_Cairo(acm_variant, width_in_points=acm_variant.template.d['GP']['mm_r_so'].value*2.1, 
                                                            height_in_points=acm_variant.template.d['GP']['mm_r_so'].value*2.1 )
        if 'PMSM' in acm_variant.template.name:
            toolCairo.draw_spmsm(acm_variant)
        elif 'Alternator' in acm_variant.template.name:
            toolCairo.draw_doubly_salient(acm_variant)
        elif 'FSPM' in acm_variant.template.name:
            toolCairo.draw_doubly_salient(acm_variant, bool_draw_whole_model=False)
        else:
            raise

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[4] Optimization Part' Multi-Objective Optimization
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def part_optimization(self):
        ad = self.ad
        ad.init_logger(prefix='acmdm')

        # [4.1] Get bounds

        # [4.3] MOO (need to share global variables to the Problem class)
        from acm_designer import get_bad_fintess_values
        import logging, builtins, utility_moo
        import pygmo as pg
        ad.counter_fitness_called = 0
        ad.counter_fitness_return = 0
        builtins.ad = ad # share global variable between modules # https://stackoverflow.com/questions/142545/how-to-make-a-cross-module-variable
        import Problem_BearinglessSynchronousDesign # must import this after __builtins__.ad = ad
        # print('[acmop.py]', builtins.ad)

        ################################################################
        # MOO Step 1:
        #   Create UserDefinedProblem and create population
        #   The magic method __init__ cannot be fined for UDP class
        ################################################################
        # [4.3.1] Basic setup
        _, prob = Problem_BearinglessSynchronousDesign.get_prob()
        popsize = self.fea_config_dict["moo.popsize"]

        print('[acmop.py]', '-'*40 + '\n[acmop.py] Pop size is', popsize)

        # [4.3.2] Generate the pop
        if False:
            pop = pg.population(prob, size=popsize) 
        # Add Restarting Feature when generating pop
        else:

            # 检查swarm_data.txt，如果有至少一个数据，返回就不是None。
            print(f'[acmop.py] Check for swarm data from: {self.select_spec}.json ...')
            self.ad.acm_template.build_x_denorm()
            # swarm_data_file = ad.   read_swarm_data_json(self.select_spec, self.ad.acm_template.x_denorm_dict)
            number_of_chromosome = ad.analyzer.number_of_chromosome

            # case 1: swarm_data.txt exists # Restarting feature related codes
            if number_of_chromosome != 0:

                number_of_finished_iterations                       = number_of_chromosome // popsize
                number_of_finished_chromosome_in_current_generation = number_of_chromosome % popsize

                # 如果刚好整除，把余数0改为popsize
                if number_of_finished_chromosome_in_current_generation == 0:
                    number_of_finished_chromosome_in_current_generation = popsize
                    print(f'\tThere are {number_of_chromosome} chromosomes found in {ad.swarm_data_file}.')
                    print('\tWhat is the odds! The script just stopped when the evaluation of the whole pop is finished.')
                    print('\tSet number_of_finished_chromosome_in_current_generation to popsize %d'%(number_of_finished_chromosome_in_current_generation))

                print('[acmop.py] This is a restart of '+ self.path2SwarmData)
                print('\tNumber of finished iterations is %d'%(number_of_finished_iterations))
                # print('This means the initialization of the population class is interrupted. So the pop in swarm_data.txt is used as the survivor.')

                # 这些计数器的值永远都是评估过的chromosome的个数。
                ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome
                print('[acmop.py] ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome = %d'%(number_of_chromosome))

                # 禁止在初始化pop时运行有限元
                ad.flag_do_not_evaluate_when_init_pop = True

                # 初始化population，如果ad.flag_do_not_evaluate_when_init_pop是False，那么就说明是 new run，否则，整代个体的fitness都是[0,0,0]。
                pop = pg.population(prob, size=popsize)

                # 如果整代个体的fitness都是[0,0,0]，那就需要调用set_xf，把txt文件中的数据写入pop。如果发现数据的个数不够，那就调用set_x()来产生数据，形成初代个体。
                if ad.flag_do_not_evaluate_when_init_pop == True:
                    pop_array = pop.get_x()
                    if number_of_chromosome <= popsize: # 个体数不够一代的情况
                        for i in range(popsize):
                            if i < number_of_chromosome: #number_of_finished_chromosome_in_current_generation:
                                pop.set_xf(i, ad.   swarm_data[i][:-3], ad.   swarm_data[i][-3:])
                            else:
                                print('[acmop.py] Set "ad.flag_do_not_evaluate_when_init_pop" to False...')
                                ad.flag_do_not_evaluate_when_init_pop = False
                                print('[acmop.py] Calling pop.set_x()---this is a restart for individual#%d during pop initialization.'%(i))
                                print('[acmop.py]', i, 'get_fevals:', prob.get_fevals())
                                pop.set_x(i, pop_array[i]) # evaluate this guy
                    else:
                        # 新办法，直接从swarm_data.txt（相当于archive）中判断出当前最棒的群体
                        swarm_data_on_pareto_front = utility_moo.learn_about_the_archive(prob, ad.   swarm_data, popsize, self.fea_config_dict)
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
                print('[acmop.py] Nothing exists in swarm_data.txt.\nThis is a whole new run.')
                ad.flag_do_not_evaluate_when_init_pop = False
                pop = pg.population(prob, size=popsize)

            # this flag must be false before moving on
            ad.flag_do_not_evaluate_when_init_pop = False

        print('[acmop.py]', '-'*40, '\nPop is initialized:\n', pop)
        hv = pg.hypervolume(pop)
        quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
        print('[acmop.py] quality_measure: %g'%(quality_measure))
        # raise KeyboardInterrupt

        # 初始化以后，pop.problem.get_fevals()就是popsize，但是如果大于popsize，说明“pop.set_x(i, pop_array[i]) # evaluate this guy”被调用了，说明还没输出过 survivors 数据，那么就写一下。
        if pop.problem.get_fevals() > popsize:
            print('[acmop.py] Write survivors.')
            ad.   write_swarm_survivor(pop, ad.counter_fitness_return)


        ################################################################
        # MOO Step 2:
        #   Select algorithm (another option is pg.nsga2())
        ################################################################
        # [4.3.3] Selecting algorithm
        # Don't forget to change neighbours to be below popsize (default is 20) decomposition="bi"
        algo = pg.algorithm(pg.moead(gen=1, weight_generation="grid", decomposition="tchebycheff", 
                                     neighbours=int(popsize/4), 
                                     CR=1, F=0.5, eta_m=20, 
                                     realb=0.9, 
                                     limit=2, preserve_diversity=True)) # https://esa.github.io/pagmo2/docs/python/algorithms/py_algorithms.html#pygmo.moead
        print('[acmop.py]', '-'*40, '\n', algo)
        print('\tThe neighbourhood size is set to 1/4 of the popsize', int(popsize/4))
        # quit()

        ################################################################
        # MOO Step 3:
        #   Begin optimization
        ################################################################
        # [4.3.4] Begin optimization
        # number_of_chromosome = ad.   read_swarm_data(self.select_spec)
        # swarm_data_file = ad.   read_swarm_data_json(self.select_spec, self.ad.acm_template.x_denorm_dict)
        number_of_chromosome = ad.analyzer.number_of_chromosome
        number_of_finished_iterations = number_of_chromosome // popsize
        number_of_iterations = 50
        logger = logging.getLogger(__name__)

        for _ in range(number_of_finished_iterations, number_of_iterations):
            msg = '[acmop.py] This is iteration #%d. '%(_)
            print(msg)
            logger.info(msg)
            pop = algo.evolve(pop)

            msg += 'Write survivors to file. '
            ad.   write_swarm_survivor(pop, ad.counter_fitness_return)

            hv = pg.hypervolume(pop)
            quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
            msg += 'Quality measure by hyper-volume: %g'% (quality_measure)
            print('[acmop.py]', msg)
            logger.info(msg)

            utility_moo.my_print(ad, pop, _)
            # my_plot(fits, vectors, ndf)
        pass

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[5] Report Part'
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def reproduce_design_from_jsonpickle(self, fname, bool_evaluate=False):
        # 从json重构JMAG模型（x_denorm）
        # Recover a design from its jsonpickle-object file
        import utility_json
        try:
            reproduced_design_variant = utility_json.from_json_recursively(fname, load_here=self.ad.fea_config_dict['output_dir']+'jsonpickle/')
            reproduced_design_variant.reproduce_wily()
            if bool_evaluate:
                reproduced_design_variant.build_jmag_project(reproduced_design_variant.project_meta_data)
        except FileNotFoundError as e:
            raise e
        self.reproduced_design_variant = reproduced_design_variant
        return reproduced_design_variant

    def part_post_optimization_analysis(self, project_name):
        # Status report: Generation, individuals, geometry as input, fea tools, performance as output (based on JSON files)
        # Do not show every step, but only those key steps showing how this population is built 
        # refer to D:\DrH\Codes\visualize

        ## Do `streamlit run visualize.py` and find an optimal design first
        number_of_chromosome = self.ad.read_swarm_data_json(self.select_spec) # ad.swarm_data 在此处被赋值

        _swarm_data          = self.ad.swarm_data
        _swarm_project_names = self.ad.swarm_data_container.project_names

        best_index = _swarm_project_names.index(project_name)
        best_chromosome = _swarm_data[best_index]

        print('---Module 5')
        print(project_name, best_index, best_chromosome)

        self.reproduce_design_from_x_denorm_and_acm_template(best_chromosome)

    def reproduce_design_from_x_denorm_and_acm_template(self, best_chromosome):
        # re-build the jmag project
        import numpy as np
        x_denorm = np.array(best_chromosome[:-3]) 

        # evaluate design (with json output)
        _ = self.ad.evaluate_design_json_wrapper(self.acm_template, x_denorm, counter=self.ad.counter_fitness_called)

    @staticmethod
    def load_settings(select_spec, select_fea_config_dict, project_loc=None, path2SwarmData=None, bool_post_processing=False):
        __file__dirname_as_in_python39 = os.path.dirname(os.path.abspath(__file__))

        with open((__file__dirname_as_in_python39)+'/machine_specifications.json', 'r') as f:
            raw_specs = json.load(f)
        with open((__file__dirname_as_in_python39)+'/machine_simulation.json', 'r') as f:
            raw_fea_config_dicts = json.load(f)

        spec_input_dict = raw_specs[select_spec]['Inputs']
        fea_config_dict = raw_fea_config_dicts[select_fea_config_dict]
        fea_config_dict['bool_post_processing'] = bool_post_processing

        # import where_am_i
        # where_am_i.where_am_i_v2(fea_config_dict, bool_post_processing)
        def get_pc_name():
            import platform
            import socket
            n1 = platform.node()
            n2 = socket.gethostname()
            n3 = os.environ["COMPUTERNAME"]
            if n1 == n2 == n3:
                return n1
            elif n1 == n2:
                return n1
            elif n1 == n3:
                return n1
            elif n2 == n3:
                return n2
            else:
                raise Exception("Computer names are not equal to each other.")

        dir_parent = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/'
        dir_codes  = os.path.abspath(os.path.dirname(__file__)) + '/'
        pc_name = get_pc_name()
        os.chdir(dir_codes)
        print('[acmop.py] CD to:', dir_codes)
        fea_config_dict['dir.parent'] = dir_parent
        fea_config_dict['pc_name']    = pc_name

        if path2SwarmData is None:
            path2SwarmData = project_loc + select_spec.replace(' ', '_')+'/'
        if project_loc is None:
            project_loc = os.path.abspath(os.path.join(path2SwarmData, '..',))

        output_dir = fea_config_dict['output_dir'] = path2SwarmData

        # create output folder only when not post-processing? No, sometimes in post-processing we run FEA simulation.
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        with open(output_dir+'acmop-settings.txt', 'w') as f:
            f.write(select_spec + ' | ' + select_fea_config_dict)
        # print(spec_input_dict)
        # quit()

        return spec_input_dict, fea_config_dict

def main(number_which_part):
    mop = AC_Machine_Optiomization_Wrapper(
        # select_spec='IM Q24p1y9 Qr32 Round Bar',
        # select_fea_config_dict = '#019 JMAG IM Nine Variables',

        # select_spec            = 'PMSM Q12p4y1 PEMD-2020', #'PMSM Q24p1y9 PEMD'
        # select_spec            = 'Flux Alternator 1955',
        # select_spec              = "FSPM-12s14pp",
        select_spec              = "FSPM-12s10pp",
        select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting',
        # select_fea_config_dict = "#029 JMAG PMSM No-load EMF",
        # select_fea_config_dict = '#04 FEMM PMSM Evaluation Setting',

        project_loc            = fr'../_default/',
        bool_show_GUI          = False
        # TODO: make bool_show_GUI a property of class (see the codes in unit conversion)
    )

    #########################
    # Call the five modules
    #########################
    if number_which_part == 1:
        mop.part_winding() # Module 1
    elif number_which_part == 2:
        mop.acm_template   # Module 2 (the execution code has been moved to the end of __post_init__ of AC_Machine_Optiomization_Wrapper)
    elif number_which_part == 3:
        mop.part_evaluation() # Module 3
    elif number_which_part == 31:
        mop.part_evaluation_geometry()
    elif number_which_part == 4:
        mop.part_optimization() # Module 4
    elif number_which_part == 5:
        # motor_design_variant = mop.reproduce_design_from_jsonpickle('p4ps5-Q12y1-0999') # Module 5 - reproduction of any design variant object
        motor_design_variant = mop.reproduce_design_from_jsonpickle('__indInitial.json')
    elif number_which_part == 51:
        # mop.part_post_optimization_analysis(project_name='proj212-SPMSM_IDQ12p1s1') # Module 5
        mop.part_post_optimization_analysis(project_name='proj12-SPMSM_IDQ12p4s1') # Module 5 - visualize swarm data
    else:
        pass
    return mop

if __name__ == '__main__':
    # main(31)
    # main(3)
    main(4)

    ''' Interactive variable checking examples:

        Example 1:
            open cmd.exe
            cd codes3
            python 
            >>> import acmop
            >>> mop = acmop.main()
            >>> mop.reproduced_design_variant.analyzer.g

        Example 2:
            >>> GP = mop.ad.acm_variant.template.d['GP'] 
            >>> GP['mm_r_ri'].value + GP['mm_d_ri'].value + GP['mm_d_rp'].value 
            46.7464829275686
            >>> GP['mm_r_ro'] 
            acmop_parameter(type='derived', name='outer_rotor_radius', value=47.7464829275686, bounds=[None, None], calc=<function template_machine_as_numbers.__init__.<locals>.<lambda> at 0x00000130CF961790>)

        Example 3:
            >>> dir(mop.ad.acm_variant.analyzer)

        Example 4 (show element |B| in the air gap):
            >>> import utility, acmop
            >>> from pylab import np, plt
            >>> mop = acmop.main(5)
            >>> var = mop.reproduced_design_variant
            >>> z = var.analyzer.z 
            >>> z_in_question = var.template.d["GP"]['mm_r_ro'].value + 10e-2 + 1j * 0
            >>> index, _ = utility.get_index_and_min(np.abs(z - z_in_question))
            >>> plt.plot(np.abs(var.analyzer.b[:,index])); plt.show() 
    '''

def examples_from_the_publications(bool_post_processing=True):
    # Vernier Machine
    # mop = AC_Machine_Optiomization_Wrapper(
    #         select_fea_config_dict = "#03 JMAG Non-Bearingless Motor Evaluation Setting",
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
    mop = AC_Machine_Optiomization_Wrapper(
            select_fea_config_dict = "#019 JMAG IM Nine Variables",
            select_spec            = 'IM p2ps3Qs18y4 Qr30-FSW Round Bar EquivDoubleLayer',
            project_loc            = r'D:/DrH/acmop/_default/'
        )

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
    # mop = AC_Machine_Optiomization_Wrapper(
    #         select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting',
    #         select_spec            = 'PMSM Q12p1y5 A',
    #         project_loc = fr'{os.path.dirname(__file__)}/_PEMD_2020_swarm_data_collected\_Q12p1y5_restart_from_optimal_and_reevaluate_wo_csv/',
    #         bool_show_GUI         = True
    #     )
    # mop = AC_Machine_Optiomization_Wrapper(
    #         select_fea_config_dict = '#0211 JMAG PMSM Q12p4ps5 Sub-hamonics',
    #         select_spec            = 'PMSM Q12p4y1 A',
    #         # project_loc            = fr'{os.path.dirname(__file__)}/_default/',
    #         project_loc = fr'{os.path.dirname(__file__)}/_PEMD_2020_swarm_data_collected\_Q12p4y1_restart_from_optimal_and_reevaluate_wo_csv_Subharmonics/',
    #         bool_show_GUI         = True
    #     )


    #########################
    # Call the five modules
    #########################

    # mop.part_winding() # Module 1
    # acm_template = mop.part_initialDesign() # Module 2 (moved to __post_init__)

    if not bool_post_processing:
        mop.part_evaluation() # Module 3
        # mop.part_optimization(acm_template) # Module 4
    else:
        if False:
            mop.reproduce_design_from_jsonpickle('p2ps1-Q12y3-0999')
        else:
            # mop.part_post_optimization_analysis(project_name='proj212-SPMSM_IDQ12p1s1') # Module 5
            mop.part_post_optimization_analysis(project_name='proj12-SPMSM_IDQ12p4s1') # Module 5
