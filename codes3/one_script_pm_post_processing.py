from pylab import mpl, np
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

# mpl.style.use('classic')
mpl.rcParams['legend.fontsize'] = 12.5
# mpl.rcParams['legend.family'] = 'Times New Roman'
mpl.rcParams['font.family'] = ['Times New Roman']
mpl.rcParams['font.size'] = 14.0
font = {'family' : 'Times New Roman', #'serif',
        'color' : 'darkblue',
        'weight' : 'normal',
        'size' : 14,}
textfont = {'family' : 'Times New Roman', #'serif',
            'color' : 'darkblue',
            'weight' : 'normal',
            'size' : 11.5,}

def pyx_script(ad, best_chromosome, best_idx, Q, p, proj_name):
    if proj_name is None:
        name = 'Q%dp%didx%d'%(Q,p,best_idx)
    else:
        name = 'Q%dp%didx%d%s'%(Q,p,best_idx, proj_name)

    # Plot cross section view
    import bearingless_spmsm_design
    spmsm_best = bearingless_spmsm_design.bearingless_spmsm_design(
                                        spmsm_template=ad.spec.acm_template,
                                        x_denorm=best_chromosome[:-3],
                                        counter=999,
                                        counter_loop=1
                                        )
    spmsm_best.ID = name
    import VanGogh
    tool_tikz = VanGogh.VanGogh_TikZPlotter()
    spmsm_best.draw_spmsm(tool_tikz, bool_pyx=True) # collecting track_path list for tool_tikz

    def redraw_cross_section_with_pyx(tikz, no_repeat_stator, no_repeat_rotor, mm_rotor_outer_radius, mm_air_gap_length, mm_rotor_outer_steel_radius, mm_rotor_inner_radius):
        # PyX
        import pyx
        tikz.c = pyx.canvas.canvas() # clear the canvas because we want to redraw 90 deg with the data tikz.track_path
        from copy import deepcopy
        def pyx_draw_path(path, sign=1, bool_exclude_path=False):
            if bool_exclude_path == False:
                if len(path) == 4: # line
                    tikz.draw_line(path[:2], path[2:4], untrack=True)
                else: # == 6 for arc
                    tikz.draw_arc(path[:2], path[2:4], path[4:6], relangle=sign*path[6], untrack=True)
        def rotate(_, x, y):
            return np.cos(_)*x + np.sin(_)*y, -np.sin(_)*x + np.cos(_)*y
        def is_at_stator(path):
            return np.sqrt(path[0]**2 + path[1]**2) > mm_rotor_outer_radius + 0.5*mm_air_gap_length

        print('Index   | Path data')
        for index, path in enumerate(tikz.track_path): # track_path is passed by reference and is changed by mirror

                # Failed to fill the closed path, because there is no arc-like path available.
                # p = pyx.path.line(4, 0, 5, 0) << pyx.path.line(5, 0, 5, 1) << pyx.path.line(5, 1, 4, 1)
                # p.append(path.closepath())
                # tikz.c.stroke(p)
                # tikz.c.stroke(path.rect(0, 0, 1, 1), [pyx.style.linewidth.Thick,
                #                      pyx.color.rgb.red,
                #                      pyx.deco.filled([pyx.color.rgb.green])])

            path_mirror = deepcopy(path)
            # for mirror copy (along x-axis)
            path_mirror[1] = path[1]*-1
            path_mirror[3] = path[3]*-1
            # for mirror copy (along y-axis)
            # path_mirror[0] = path[0]*-1
            # path_mirror[2] = path[2]*-1

            bool_exclude_path = False

            # rotate path and plot
            if is_at_stator(path):
                Q = no_repeat_stator
            else:
                Q = no_repeat_rotor*2



            EPS = 1e-6
            if is_at_stator(path):
                # 按照Eric的要求，把不必要的线给删了。
                if abs(path[1] + path[3]) < EPS: # 镜像对称线
                    bool_exclude_path = True
                if abs(path[0] - path[2]) + np.cos(2*np.pi/Q/2) < EPS: # 旋转对称线（特别情况，tan(90°) = ∞
                    bool_exclude_path = True                
                else:
                    if abs( abs((path[1] - path[3])/(path[0] - path[2])) - abs(np.tan(2*np.pi/Q/2)) ) < EPS: # 旋转对称线
                        bool_exclude_path = True

            if not is_at_stator(path):
                # 按照Eric的要求，把不必要的线给删了。
                if  (abs(np.sqrt(path[0]**2+path[1]**2) - mm_rotor_inner_radius)<EPS or abs(np.sqrt(path[2]**2+path[3]**2) - mm_rotor_inner_radius)<EPS) \
                    and (len(path)==4): # 转子铁芯内径到外径的线
                    bool_exclude_path = True

            #     # 特别的是，画永磁体的时候，边界要闭合哦。
            #     if abs(np.sqrt(path[0]**2+path[1]**2) - mm_rotor_outer_steel_radius) < EPS or abs(np.sqrt(path[2]**2+path[3]**2) - mm_rotor_outer_steel_radius) < EPS:
            #         bool_exclude_path = False

            # Make sure models with different outer diameters have the same scale.
            # tikz.draw_arc([125,0], [-125,0], relangle=sign*180, untrack=True)
            tikz.c.fill(pyx.path.circle(0, 0, 125), [pyx.color.transparency(1)]) # use this if THICK is used.


            _ = 2*np.pi/Q
 
            if True: # full model
                for counter in range(Q):

                    # 转子：旋转复制
                    if not is_at_stator(path):
                        path[0], path[1] = rotate(_, path[0], path[1])
                        path[2], path[3] = rotate(_, path[2], path[3])
                        pyx_draw_path(path, sign=1, bool_exclude_path=bool_exclude_path)

                    # 定子：镜像+旋转复制
                    if is_at_stator(path):

                        path[0], path[1] = rotate(_, path[0], path[1])
                        path[2], path[3] = rotate(_, path[2], path[3])
                        pyx_draw_path(path, sign=1, bool_exclude_path=bool_exclude_path)
                        # print(index, '\t|', ',\t'.join(['%g'%(el) for el in path]))

                        path_mirror[0], path_mirror[1] = rotate(_, path_mirror[0], path_mirror[1])
                        path_mirror[2], path_mirror[3] = rotate(_, path_mirror[2], path_mirror[3])
                        pyx_draw_path(path_mirror, sign=-1, bool_exclude_path=bool_exclude_path)

                    # break

            else: # backup

                # 转子：旋转复制
                if not is_at_stator(path):
                    path[0], path[1] = rotate(0.5*np.pi - 0.5*0.5*_, path[0], path[1])
                    path[2], path[3] = rotate(0.5*np.pi - 0.5*0.5*_, path[2], path[3])
                    pyx_draw_path(path, sign=1, bool_exclude_path=bool_exclude_path)
                    # print(index, '\t|', ',\t'.join(['%g'%(el) for el in path]))

                    # path[0], path[1] = rotate(0.5*np.pi - 0*0.5*_, path[0], path[1])
                    # path[2], path[3] = rotate(0.5*np.pi - 0*0.5*_, path[2], path[3])
                    # pyx_draw_path(path, sign=1, bool_exclude_path=bool_exclude_path)

                # 定子：镜像+旋转复制
                if is_at_stator(path):

                    path[0], path[1] = rotate(0.5*np.pi - 0.5*_, path[0], path[1])
                    path[2], path[3] = rotate(0.5*np.pi - 0.5*_, path[2], path[3])
                    pyx_draw_path(path, sign=1, bool_exclude_path=bool_exclude_path)
                    # print(index, '\t|', ',\t'.join(['%g'%(el) for el in path]))

                    path_mirror[0], path_mirror[1] = rotate(0.5*np.pi - 0.5*_, path_mirror[0], path_mirror[1])
                    path_mirror[2], path_mirror[3] = rotate(0.5*np.pi - 0.5*_, path_mirror[2], path_mirror[3])
                    pyx_draw_path(path_mirror, sign=-1, bool_exclude_path=bool_exclude_path)

                    # 注意，所有 tack_path 中的 path 都已经转动了90度了！
                    # for mirror copy (along y-axis)
                    path[0] *= -1
                    path[2] *= -1
                    pyx_draw_path(path, sign=-1, bool_exclude_path=bool_exclude_path)

                    path_mirror[0] *= -1
                    path_mirror[2] *= -1
                    pyx_draw_path(path_mirror, sign=1, bool_exclude_path=bool_exclude_path)

    redraw_cross_section_with_pyx(tool_tikz, spmsm_best.Q, spmsm_best.p, spmsm_best.Radius_OuterRotor, spmsm_best.Length_AirGap, spmsm_best.Radius_OuterRotorSteel, spmsm_best.Radius_InnerRotor)
    tool_tikz.c.writePDFfile("Figure_selected_optimal_design_%s"%(name))
    tool_tikz.c.writeEPSfile("Figure_selected_optimal_design_%s"%(name))
    tool_tikz.c.writeSVGfile("Figure_selected_optimal_design_%s"%(name))
    print('Write to pdf file: Figure_selected_optimal_design_%s.pdf.'%(name))
    # os.system('start %s'%("selected_optimal_design%s.pdf"%(spmsm_best.name)))
    # quit()


    # from pdf2image import convert_from_path
    # import PIL
    # PIL.Image.MAX_IMAGE_PIXELS = 933120000
    # pages = convert_from_path("selected_optimal_design_%s.pdf"%(name), 300)
    # for page in pages:
    #     page.save("selected_optimal_design_%s.jpg"%(name), 'JPEG')

def selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[20, -0.92, 200], best_idx=None, Q=None, p=None, proj_name=None):

    for idx, chromosome in enumerate(_swarm_data):
        OC = upper_bound_objectives[0] # ripple performance
        OB = upper_bound_objectives[1] # -efficiency
        OA = upper_bound_objectives[2] # cost
        if chromosome[-1] < OC and chromosome[-2] < OB and chromosome[-3] < OA:
            print(idx, _swarm_project_names[idx], chromosome[::-1])

            if idx == best_idx:
                best_chromosome = chromosome
                pyx_script(ad, best_chromosome, best_idx, Q, p, proj_name)

                # re-build the jmag project
                if True:
                    x_denorm = best_chromosome[:-3]
                    cost_function, f1, f2, f3, FRW, \
                    normalized_torque_ripple, \
                    normalized_force_error_magnitude, \
                    force_error_angle = \
                        ad.evaluate_design(ad.spec.acm_template, x_denorm, ad.counter_fitness_called)

                return best_chromosome

def post_processing(ad, fea_config_dict):

    ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + ad.solver.fea_config_dict['run_folder'] 
    number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
    print('number_of_chromosome =', number_of_chromosome)

    _swarm_data          = ad.solver.swarm_data
    _swarm_project_names = ad.solver.swarm_data_container.project_names

    # print('-------for tutorial\n'*3)
    # print(_swarm_data)
    # print('-------for tutorial ends\n'*3)
    # quit()

    # print('-'*40+'\nY730' + '\n      L_g,    w_st,   w_rt,   theta_so,   w_ro,    d_so,    d_ro,    -TRV,    -eta,    OC.')
    print('-'*40+'\n%s | %s'%(fea_config_dict['pc_name'], fea_config_dict['run_folder']) + '\n      _____________________________________________________________________________________')

    # Select optimal design by user-defined criteria
    if r'run#62399' in fea_config_dict['run_folder']:
        selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[22, -0.92, 200], best_idx=3228, proj_name='proj3280', Q=ad.spec.acm_template.Q, p=ad.spec.acm_template.p)
    elif r'run#62499' in fea_config_dict['run_folder']:
        selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[22, -0.94, 200], best_idx=4235, proj_name='proj4240', Q=ad.spec.acm_template.Q, p=ad.spec.acm_template.p)
    elif r'run#62599' in fea_config_dict['run_folder']:
        selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[13, -0.94, 200], best_idx=3538, proj_name='proj3543', Q=ad.spec.acm_template.Q, p=ad.spec.acm_template.p)
    elif r'run#62699' in fea_config_dict['run_folder']:
        selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[3, -0.95, 160], best_idx=2623, proj_name='proj2761', Q=ad.spec.acm_template.Q, p=ad.spec.acm_template.p)
    elif r'run#62799' in fea_config_dict['run_folder']:
        selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[6, -0.94, 200], best_idx=17, proj_name='proj18', Q=ad.spec.acm_template.Q, p=ad.spec.acm_template.p)
    elif r'run#62899' in fea_config_dict['run_folder']:
        selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[14, -0.94, 200], best_idx=3680, proj_name='proj3681', Q=ad.spec.acm_template.Q, p=ad.spec.acm_template.p)
    else:
        raise Exception('Not implmented in post_processing.py')

    # Set the output_dir back! Obsolete?
    ad.solver.output_dir = ad.solver.fea_config_dict['dir_parent'] + ad.solver.fea_config_dict['run_folder']



if __name__ == '__main__':
    # PEMD 2020

    from utility import my_execfile
    bool_post_processing = True # solve or post-processing

    my_execfile('./default_setting.py', g=globals(), l=locals())
    fea_config_dict
    fea_config_dict['bool_post_processing'] = bool_post_processing

    # Combined winding PMSM
    fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.95
    fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.05

    # Collect all data 
    list_of_swarm_data_on_pareto_front = []
    fig, ax = plt.subplots(constrained_layout=False)
    fig.set_rasterized(True) # https://stackoverflow.com/questions/19638773/matplotlib-plots-lose-transparency-when-saving-as-ps-eps
    plt.subplots_adjust(left=None, bottom=None, right=0.85, top=None, wspace=None, hspace=None)
    index = 0
    for run_folder, spec_file, marker, label in zip(
                                                      [ r'run#62599/',  # spec_PEMD_BPMSM_Q6p1)
                                                        r'run#62399/',  # spec_ECCE_PMSM_ (Q6p2)
                                                        r'run#62899/',  # spec_PEMD_BPMSM_Q12p1
                                                        r'run#62499/',  # spec_PEMD_BPMSM_Q12p2
                                                        r'run#62699/',  # spec_PEMD_BPMSM_Q12p4)
                                                        r'run#62799/'], # spec_PEMD_BPMSM_Q24p1
                                                      [ './spec_PEMD_BPMSM_Q6p1.py',
                                                        './spec_ECCE_PMSM_Q6p2.py',
                                                        './spec_PEMD_BPMSM_Q12p1.py',
                                                        './spec_PEMD_BPMSM_Q12p2.py',
                                                        './spec_PEMD_BPMSM_Q12p4.py',
                                                        './spec_PEMD_BPMSM_Q24p1.py'],
                                                      [ '$1$', '$2$', '$3$', '$4$', '$5$', '$6$' ], # [ ',', '+', '.', 'o', '*' ]
                                                      [ 'Q6p1', 'Q6p2', 'Q12p1', 'Q12p2', 'Q12p4', 'Q24p1' ] 
                                                    ):  
        # # debug
        # index += 1
        # if index < 5:
        #     continue
        # if '623' not in run_folder and '627' not in run_folder:
        # if '623' not in run_folder:
        # if '627' not in run_folder:
        if '626' not in run_folder:
            continue

        fea_config_dict['run_folder'] = run_folder

        my_execfile(spec_file, g=globals(), l=locals()) # Q=24, p=1, ps=2
        spec

        # Adopt Bianchi 2006 for a SPM motor template
        spec.build_pmsm_template(fea_config_dict, im_template=None)

        # select motor type here
        spec.acm_template = spec.pmsm_template
        print('Build ACM template...')

        import acm_designer
        ad = acm_designer.acm_designer(fea_config_dict, spec)
        ad.bool_re_evaluate = False
        ad.flag_do_not_evaluate_when_init_pop = True

        ad.bounds_denorm = spec.acm_template.get_classic_bounds(which_filter='VariableSleeveLength')
        ad.bound_filter  = spec.acm_template.bound_filter

        ad.counter_fitness_called = 99999
        ad.counter_fitness_return = 99999

        __builtins__.ad = ad # share global variable between modules # https://stackoverflow.com/questions/142545/how-to-make-a-cross-module-variable

        # # re-build the jmag project for the template
        # if True:
        #     x_denorm = spec.acm_template.build_x_denorm()
        #     cost_function, f1, f2, f3, FRW, \
        #     normalized_torque_ripple, \
        #     normalized_force_error_magnitude, \
        #     force_error_angle = \
        #         ad.evaluate_design(ad.spec.acm_template, x_denorm, ad.counter_fitness_called)
        #     continue

        # select the optimal design and draw its cross section sketch
        if bool_post_processing == True:
            post_processing(ad, fea_config_dict)

        # plot the Pareto front for the archive
        if bool_post_processing == True:
            import utility_moo

            from Problem_BearinglessSynchronousDesign import Problem_BearinglessSynchronousDesign
            udp = Problem_BearinglessSynchronousDesign()
            import pygmo as pg
            prob = pg.problem(udp)
            popsize = 78

            swarm_data_on_pareto_front = utility_moo.learn_about_the_archive(prob, ad.solver.swarm_data, popsize, fea_config_dict, bool_plot_and_show=False)
            print(len(swarm_data_on_pareto_front), len(swarm_data_on_pareto_front[0]))

            # list_of_swarm_data_on_pareto_front.append(swarm_data_on_pareto_front)

            fits = [el[-3:] for el in swarm_data_on_pareto_front]
            list_alpha_st = [el[0] for el in swarm_data_on_pareto_front]
            list_ripple_sum = [el[-1] for el in swarm_data_on_pareto_front]
            scatter_handle = utility_moo.my_2p5d_plot_non_dominated_fronts(fits, comp=[0,1], marker=marker, up_to_rank_no=1, ax=ax, fig=fig, no_colorbar=True, z_filter=20, label=label)

            # ripple_ax = plt.figure().gca()
            # ripple_ax.plot(list_alpha_st, list_ripple_sum, 'ko')
            # ripple_ax.set_xlabel('alpha_st')
            # ripple_ax.set_ylabel('ripple sum')

        # plot other performance values other than the 3 objectives, such as Ea and FRW.
        if bool_post_processing == False:
            # plt.figure()
            # plt.plot(ad.solver.swarm_data_container.FRW, '--')
            # plt.plot(ad.solver.swarm_data_container.l_FRW, 'o')
            
            # Is wide slot open design existing on all Pareto front? 
            # Plot slot open vs. force ripple plot using the archive. 
            slot_tip_open_ratio = np.array(ad.solver.swarm_data_container.deg_alpha_st)/180*np.pi * np.array(ad.solver.swarm_data_container.mm_r_si)*1e-3 / (np.array(ad.solver.swarm_data_container.mm_w_st)*1e-3)
            plt.figure()

            # plt.plot(slot_tip_open_ratio, ad.solver.swarm_data_container.Ea, 'o')
            # plt.ylim([0,10])
            # plt.ylabel('$E_a$ [deg]')

            plt.plot(slot_tip_open_ratio, ad.solver.swarm_data_container.Em, 'o')
            plt.ylim([0,0.30])
            plt.xlabel('$\\alpha_{st}r_{si}/w_{st}$ [1]')
            plt.ylabel('$E_m$ [%]')


    # color bar
    cbar_ax = fig.add_axes([0.875, 0.15, 0.02, 0.7])
    cbar_ax.get_yaxis().labelpad = 10
    clb = fig.colorbar(scatter_handle, cax=cbar_ax)
    clb.ax.set_ylabel('Ripple Performance Sum [1]', rotation=270, labelpad=14)
    ax.legend()
    ax.plot(np.arange(50,250), np.ones(200)*-97, '--k')
    ax.set_xlim([50, 250])
    ax.set_ylim([-99, -75])
    ax.text( 31, -97.3, '$-97$', fontdict=font)
    # ax.set_yticks([-75, -85, -90, -95, -97], [-75, -85, -90, -95, -97])

    # fig.savefig(r'./Figure_Combined_Pareto_Front.eps', format='eps', dpi=1000)
    # fig.savefig(r'./Figure_Combined_Pareto_Front.png', format='png', dpi=300)
    fig.savefig(r'./Figure_Combined_Pareto_Front.svg', format='svg', dpi=1000) # fast
    plt.show()




