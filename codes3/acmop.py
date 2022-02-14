# Please use shortcut "ctrl+k,ctrl+1" to fold the code for better navigation
# Please use shortcut "ctrl+k,ctrl+2" to fold the code for better navigation

import os, main_utility, acm_designer, bearingless_spmsm_design, vernier_motor_design, bearingless_induction_design

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
        self.spec_input_dict, self.fea_config_dict = \
                main_utility.load_settings( self.select_spec, 
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

        r""" <ACMOP parent dir> = D:/DrH/Codes/acmop/
             <Data folder name> = <ACMOP parent dir>/_default/, /_WenboVShapeVernier/, or /_PEMD_2020_swarm_data_collected/_Q12p4y1_restart_from_optimal_and_reevaluate_wo_csv_Subharmonics/
             <Project location> = <Data folder name>
             <path2SwarmData>   = <Project location>/PMSM_Q12p4y1_A/(swarm_data.txt)
             <fea_config_dict['run_folder']> = <path2SwarmData>
             <output_dir> = <path2SwarmData>
        """

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
        if 'PMSM' in self.select_spec:
            function = bearingless_spmsm_design.bearingless_spmsm_template
        elif 'PMVM' in self.select_spec:
            function = vernier_motor_design.vernier_motor_VShapePM_template
        elif 'IM' in self.select_spec:
            function = bearingless_induction_design.bearingless_induction_template
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
            motor_design_variant = self.ad.evaluate_design_json_wrapper(self.ad.acm_template, x_denorm)
            print('[acmop.py] Listing analyzer.spec_performance_dict:')
            for k,v in motor_design_variant.analyzer.spec_performance_dict.items():
                print('\t', k, v)

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

        # [4.3] MOO (need to share global variables to the Problem class)
        from acm_designer import get_bad_fintess_values
        import logging, builtins
        import utility_moo
        import pygmo as pg
        ad.counter_fitness_called = 0
        ad.counter_fitness_return = 0
        builtins.ad = ad # share global variable between modules # https://stackoverflow.com/questions/142545/how-to-make-a-cross-module-variable
        import codes3.Problem_BearinglessSynchronousDesign # must import this after __builtins__.ad = ad
        # print('[acmop.py]', builtins.ad)

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
            number_of_chromosome = ad.solver.read_swarm_data(self.select_spec)

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
                    print('[acmop.py] Old pop:')
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
                                print('[acmop.py] Set "ad.flag_do_not_evaluate_when_init_pop" to False...')
                                ad.flag_do_not_evaluate_when_init_pop = False
                                print('[acmop.py] Calling pop.set_x()---this is a restart for individual#%d during pop initialization.'%(i))
                                print('[acmop.py]', i, 'get_fevals:', prob.get_fevals())
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
                    print('[acmop.py] Nothing exists in swarm_data.txt.\nThis is a whole new run.')
                    ad.flag_do_not_evaluate_when_init_pop = False
                    pop = pg.population(prob, size=popsize)

                # case 2-B: swarm_data.txt does not exist and this is a re-evalation run (without csv)
                else:
                    print('[acmop.py] Nothing exists in swarm_data.txt.\nRe-start from %s'%(path_to_archive))
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

        print('[acmop.py]', '-'*40, '\nPop is initialized:\n', pop)
        hv = pg.hypervolume(pop)
        quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
        print('[acmop.py] quality_measure: %g'%(quality_measure))
        # raise KeyboardInterrupt

        # 初始化以后，pop.problem.get_fevals()就是popsize，但是如果大于popsize，说明“pop.set_x(i, pop_array[i]) # evaluate this guy”被调用了，说明还没输出过 survivors 数据，那么就写一下。
        if pop.problem.get_fevals() > popsize:
            print('[acmop.py] Write survivors.')
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
        print('[acmop.py]', '-'*40, '\n', algo)
        # quit()

        ################################################################
        # MOO Step 3:
        #   Begin optimization
        ################################################################
        # [4.3.4] Begin optimization
        number_of_chromosome = ad.solver.read_swarm_data(self.select_spec)
        number_of_finished_iterations = number_of_chromosome // popsize
        number_of_iterations = 20
        logger = logging.getLogger(__name__)
        # try:
        if True:
            for _ in range(number_of_finished_iterations, number_of_iterations):
                ad.number_of
                msg = '[acmop.py] This is iteration #%d. '%(_)
                print(msg)
                logger.info(msg)
                pop = algo.evolve(pop)

                msg += 'Write survivors to file. '
                ad.solver.write_swarm_survivor(pop, ad.counter_fitness_return)

                hv = pg.hypervolume(pop)
                quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
                msg += 'Quality measure by hyper-volume: %g'% (quality_measure)
                print('[acmop.py]', msg)
                logger.info(msg)

                utility_moo.my_print(ad, pop, _)
                # my_plot(fits, vectors, ndf)
        # except Exception as e:
        #     print(pop.get_x())
        #     print(pop.get_f().tolist())
        #     raise e        
        pass

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # '[5] Report Part'
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def reproduce_design_from_jsonpickle(self, fname):
        # 从json重构JMAG模型（x_denorm）
        # Recover a design from its jsonpickle-object file
        import utility_json
        try:
            motor_design_variant = utility_json.from_json_recursively(fname, load_here=self.ad.output_dir+'jsonpickle/')
            motor_design_variant.reproduce_wily()
            motor_design_variant.build_jmag_project(motor_design_variant.project_meta_data)
        except FileNotFoundError as e:
            print(e)
            print("你有没有搞错啊？json文件找不到啊，忘记把bool_post_processing改回来了？")
        return motor_design_variant

    def part_post_optimization_analysis(self, project_name):
        # Status report: Generation, individuals, geometry as input, fea tools, performance as output (based on JSON files)
        # Do not show every step, but only those key steps showing how this population is built 
        # refer to D:\DrH\Codes\visualize

        ## Do `streamlit run visualize.py` and find an optimal design first
        self.ad.solver.output_dir = self.path2SwarmData
        number_of_chromosome = self.ad.solver.read_swarm_data(self.select_spec) # ad.solver.swarm_data 在此处被赋值

        _swarm_data          = self.ad.solver.swarm_data
        _swarm_project_names = self.ad.solver.swarm_data_container.project_names

        best_index = _swarm_project_names.index(project_name)
        best_chromosome = _swarm_data[best_index]

        print('---Module 5')
        print(project_name, best_index, best_chromosome)

        self.reproduce_design_from_x_denorm_and_acm_template(best_chromosome)

    def reproduce_design_from_x_denorm_and_acm_template(self, best_chromosome):
        # re-build the jmag project
        x_denorm = np.array(best_chromosome[:-3]) 

        # evaluate design (with json output)
        cost_function, f1, f2, f3, FRW, \
            normalized_torque_ripple, \
            normalized_force_error_magnitude, \
            force_error_angle = self.ad.evaluate_design_json_wrapper(self.acm_template, x_denorm, counter=self.ad.counter_fitness_called)

def main():
    mop = AC_Machine_Optiomization_Wrapper(
        # select_spec='IM Q24p1y9 Qr32 Round Bar',
        # select_fea_config_dict = '#019 JMAG IM Nine Variables',

        select_spec            = 'PMSM Q12p4y1 PEMD-2020', #'PMSM Q18p4y2 Beijing ShiDaiChaoQun',
        # select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting',
        select_fea_config_dict = '#04 FEMM PMSM Evaluation Setting',

        project_loc            = fr'../_default/',
        bool_show_GUI          = True
    )

    #########################
    # Call the five modules
    #########################
    # mop.part_winding() # Module 1
    # mop.acm_template   # Module 2 (the execution code has been moved to the end of __post_init__ of AC_Machine_Optiomization_Wrapper)
    if True:
        mop.part_evaluation() # Module 3
        # mop.part_optimization(acm_template) # Module 4
    else:
        if True:
            motor_design_variant = mop.reproduce_design_from_jsonpickle('p4ps5-Q12y1-0999') # Module 5 - reproduction of any design variant object
        else:
            # mop.part_post_optimization_analysis(project_name='proj212-SPMSM_IDQ12p1s1') # Module 5
            mop.part_post_optimization_analysis(project_name='proj12-SPMSM_IDQ12p4s1') # Module 5 - visualize swarm data

if __name__ == '__main__':
    from pylab import np
    main()


    ''' Interactive variable checking examples:

        Example 1:
        >>> GP = mop.ad.acm_variant.template.d['GP'] 
        >>> GP['mm_r_ri'].value + GP['mm_d_ri'].value + GP['mm_d_rp'].value 
        46.7464829275686
        >>> GP['mm_r_or'] 
        acmop_parameter(type='derived', name='outer_rotor_radius', value=47.7464829275686, bounds=[None, None], calc=<function template_machine_as_numbers.__init__.<locals>.<lambda> at 0x00000130CF961790>)

        Example 2:
        >>> dir(mop.ad.acm_variant.analyzer)
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
