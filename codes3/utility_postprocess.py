
#禁止在cache时打印 
#also added to  D:\DrH\bopt-python\codes3\pyrhonen_procedure_as_function.py
#also added to  D:\DrH\bopt-python\codes3\winding_layout.py
#also added to  D:\DrH\bopt-python\codes3\population.py
# def print(*arg, **kwarg):
#     pass

# from visualize_config import * # markdown strings and global variables

# Consider to replace this with CjhStylePlot
from pylab import mpl, plt, np
def get_plot():
    # mpl.style.use('classic')
    mpl.rcParams['mathtext.fontset'] = 'stix'
    # mpl.rcParams['font.family'] = 'STIXGeneral'
    mpl.rcParams['font.family'] = 'sans-serif' # (Streamlit) 2021-03-16 20:37:40.259 font.family must be one of (serif, sans-serif, cursive, monospace) when text.usetex is True. serif will be used by default. 
    mpl.rcParams['legend.fontsize'] = 12.5
    mpl.rcParams['font.size'] = 14.0
    font = {'family' : 'Times New Roman', #'serif',
            'color' : 'darkblue',
            'weight' : 'normal',
            'size' : 14,}
    textfont = {'family' : 'Times New Roman', #'serif',
                'color' : 'darkblue',
                'weight' : 'normal',
                'size' : 11.5,}
    # plt.rc('text', usetex=True) # https://github.com/matplotlib/matplotlib/issues/4495/
    # plt.rc('pgf', texsystem='pdflatex')
    fig, ax = plt.subplots(figsize=(8,5), constrained_layout=False)
    fig.set_rasterized(True) # https://stackoverflow.com/questions/19638773/matplotlib-plots-lose-transparency-when-saving-as-ps-eps
    plt.subplots_adjust(left=None, bottom=None, right=0.85, top=None, wspace=None, hspace=None)

    return fig, ax, font

def selection_criteria(ad, _swarm_data, _swarm_project_names, upper_bound_objectives=[20, -0.92, 200], best_idx=None, Q=None, p=None, proj_name=None):

    OC = upper_bound_objectives[0] # ripple performance
    OB = upper_bound_objectives[1] # -efficiency
    OA = upper_bound_objectives[2] # -TRV
    return_tuples = []
    for idx, chromosome in enumerate(_swarm_data):
        if chromosome[-1] < OC and chromosome[-2] < OB and chromosome[-3] < OA:
            # print('\n', idx, _swarm_project_names[idx], chromosome[::-1])
            return_tuples.append((idx, _swarm_project_names[idx], chromosome[::-1]))

            # Do stuff to the best design specified by proj_name
            if proj_name is not None and proj_name in _swarm_project_names[idx]:
                best_chromosome = chromosome

                print('|||||||||||||||| Validate:', best_chromosome[::-1])
                # plot the cross-sectional view
                if False:
                    PyX_Utility.pyx_script_im(ad, best_chromosome, best_idx, Q, p, proj_name=None)
                    quit()

                # re-build the jmag project
                if True:
                    print('- Now re-build the JMAG project...')
                    x_denorm = np.array(best_chromosome[:-3]) 
                    print(x_denorm)
                    # quit()

                    # evaluate design (with json output)
                    cost_function, f1, f2, f3, FRW, \
                        normalized_torque_ripple, \
                        normalized_force_error_magnitude, \
                        force_error_angle = ad.evaluate_design_json_wrapper(ad.spec.acm_template, x_denorm, counter=ad.counter_fitness_called)

                # Sensitivity analysis with respect to alpha_st
                if False:
                    # debug
                    print('需要转换开口角度转换为alpha_st然后做敏感性分析。')
                    for name, value in zip(ad.spec.acm_template.x_denorm_names, best_chromosome):
                        print(name, '\t', value)

                    # Detect wide slot open design?
                    if True:
                        def StatorSlotOpenAngleSo_to_AlphaSt(theta, Q):
                            return 360.0/Q - theta
                        deg_alpha_st = StatorSlotOpenAngleSo_to_AlphaSt(best_chromosome[3], Q)
                        mm_w_st = best_chromosome[1]
                        mm_r_si = best_chromosome[0] + ad.spec.acm_template.Radius_OuterRotor

                        def straight_slot_alpha_st(deg_alpha_st, mm_w_st):
                            deg_alpha_st_at_w_st = ( np.pi - 2 * np.arccos(0.5*mm_w_st / mm_r_si) ) /np.pi*180
                            return deg_alpha_st_at_w_st
                        deg_alpha_st_at_w_st = straight_slot_alpha_st(deg_alpha_st, mm_w_st)
                        直槽所对应开口角度 = (2*np.pi - deg_alpha_st_at_w_st * Q) / Q

                        print('||||||||||||||||||||||||||||||||| Wide slot open design?', proj_name, Q, p, deg_alpha_st_at_w_st, deg_alpha_st, end='|||')
                        if deg_alpha_st <= deg_alpha_st_at_w_st:
                            print('Yes.')
                        else:
                            print('No.')

                    counter_bias = 700000

                    # Change to False if only plotting is needed.
                    if False:
                        x_denorm = np.array(best_chromosome[:-3])
                        deg_alpha_st = StatorSlotOpenAngleSo_to_AlphaSt(x_denorm[3], Q)
                        percent_ori = int(100*(deg_alpha_st / deg_alpha_st_at_w_st)) # (跟直槽比)

                        print('percent_ori =', percent_ori)
                        # quit()

                        # evaluate design (with json output)
                        cost_function, f1, f2, f3, FRW, \
                            normalized_torque_ripple, \
                            normalized_force_error_magnitude, \
                            force_error_angle = ad.evaluate_design_json_wrapper(ad.spec.acm_template, x_denorm, counter=counter_bias+percent_ori)

                        percent_per_step = 10
                        for step in range(-8, 5):
                            # change deg_alpha_st
                            new_deg_alpha_st = deg_alpha_st_at_w_st * (100 + step*percent_per_step)/100
                            def AlphaSt_to_StatorSlotOpenAngleSo(alpha, Q):
                                return 360.0/Q - alpha
                            x_denorm[3] = AlphaSt_to_StatorSlotOpenAngleSo(new_deg_alpha_st, Q)

                            # evaluate design (with json output)
                            cost_function, f1, f2, f3, FRW, \
                                normalized_torque_ripple, \
                                normalized_force_error_magnitude, \
                                force_error_angle = ad.evaluate_design_json_wrapper(ad.spec.acm_template, x_denorm, counter=counter_bias+100+step*percent_per_step)
                    else:
                        number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
                        _swarm_data          = ad.solver.swarm_data
                        _swarm_project_names = ad.solver.swarm_data_container.project_names

                        import re
                        list_percent = [int(re.split('-|j', name)[1])-counter_bias for name in _swarm_project_names]
                        print('|||', _swarm_project_names)
                        print('|||', list_percent)
                        print('|||', ad.solver.swarm_data_container.Ea)
                        fig = plt.figure()
                        plt.subplot(221)
                        plt.subplots_adjust(left=0.2, bottom=None, right=0.99, top=None, wspace=None, hspace=None)
                        plt.plot(list_percent[1:], [individual[-3]*1e-3 for individual in _swarm_data][1:], '-o', label=r'$O_A$ [$\rm kNm/m^3$]')
                        plt.ylabel(r'$O_A$ [$\rm kNm/m^3$]')
                        # plt.xlabel(r'$\alpha_{\rm st}/\alpha_{\rm st}^*$ [\%]'))
                        plt.title(r'$p=%d$, $Q_r=%d$'%(p, ad.spec.acm_template.Qr))
                        plt.grid()
                        plt.xticks(np.arange(20, 141, step=20))
                        plt.yticks(np.arange(-35, -19, step=5))
                        # plt.yticks(np.arange(140, 165, step=5))
                        plt.subplot(223)
                        plt.plot(list_percent[1:], [individual[-2]*100 for individual in _swarm_data][1:], '-o', label=r'$O_B$ [\%]')
                        plt.ylabel(r'$O_B$ [\%]')
                        plt.xlabel(r'$\alpha_{\rm st}/\alpha_{\rm st}^*$ [\%]')
                        plt.grid()
                        plt.xticks(np.arange(20, 141, step=20))
                        plt.yticks(np.arange(-97, -96.4, step=0.1))
                        # plt.yticks(np.arange(-0.945, -0.96, step=-0.005))
                        # plt.yticks(np.arange(-0.967, -0.97, step=-0.001))
                        plt.subplot(122)
                        plt.plot(list_percent[1:], [individual[-1] for individual in _swarm_data][1:], '-o', label=r'$O_C$ [1]')
                        plt.plot(list_percent[1:], [100*el for el in ad.solver.swarm_data_container.Trip][1:], '--x', label=r'$T_{\rm rip}$ [\%]')
                        plt.plot(list_percent[1:], [100*el for el in ad.solver.swarm_data_container.Em][1:], '--v', label=r'$E_m$ [\%]')
                        plt.plot(list_percent[1:], ad.solver.swarm_data_container.Ea[1:], '--^', label=r'$E_a$ [deg]')
                        plt.xlabel(r'$\alpha_{\rm st}/\alpha_{\rm st}^*$ [\%]')
                        plt.title(r'$p=%d$, $Q_r=%d$'%(p, ad.spec.acm_template.Qr))
                        plt.grid()
                        plt.xticks(np.arange(20, 141, step=20))
                        plt.yticks(np.arange(1, 10, step=1))
                        # plt.ylabel(r'')
                        plt.legend()
                        fname = 'p%dQr%dQs%d.pdf'%(0.5*ad.spec.acm_template.DriveW_poles, ad.spec.acm_template.Qr, ad.spec.acm_template.Qs)
                        fig.savefig(f'{os.path.dirname(__file__)}/{fname}', format='pdf', dpi=600, bbox_inches='tight', pad_inches=0.0, transparent=True)
                        # plt.show()
                        # quit()

                    # scatter_handle = pareto_front_plot_script(_swarm_data)
                    # pareto_front_plot_color_bar_etc(scatter_handle, bool_no_limit=True)
                    # plt.show()
                    # quit()

                # Sensitivity analysis with respect to yoke height
                if False:
                    counter_bias = 800000
                    x_denorm = np.array(best_chromosome[:-3])
                    mm_stator_outer_radius_ori = x_denorm[3]
                    percent_ori = 100 # (跟自己比)

                    percent_per_step = 5
                    # for step in range(-5, 3):
                    #     # change mm_stator_outer_radius
                    #     x_denorm[3] = mm_stator_outer_radius_ori * (100 + step*percent_per_step)/100

                    #     # evaluate design (with json output)
                    #     cost_function, f1, f2, f3, FRW, \
                    #         normalized_torque_ripple, \
                    #         normalized_force_error_magnitude, \
                    #         force_error_angle = ad.evaluate_design_json_wrapper(ad.spec.acm_template, x_denorm, counter=counter_bias+100+step*percent_per_step)

                    print(ad.solver.output_dir)
                    number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
                    _swarm_data          = ad.solver.swarm_data
                    _swarm_project_names = ad.solver.swarm_data_container.project_names

                    import re
                    list_percent = [int(re.split('-|j', name)[1])-counter_bias for name in _swarm_project_names]
                    print('|||', _swarm_project_names)
                    print('|||', list_percent)
                    print('|||', ad.solver.swarm_data_container.Ea)
                    plt.figure()
                    plt.subplot(221)
                    plt.plot(list_percent[:], [individual[-3] for individual in _swarm_data][:], '*', label=r'$O_A$ [USD]')
                    plt.ylabel(r'$O_A$ [USD]')
                    plt.xlabel(r'$\Delta d_{\rm sy}$ [%]')
                    plt.title('Q=%d, p=%d'%(Q,p))
                    plt.subplot(223)
                    plt.plot(list_percent[:], [individual[-2] for individual in _swarm_data][:], '*', label=r'$O_B$ [%]')
                    plt.ylabel(r'$O_B$ [%]')
                    plt.xlabel(r'$\Delta d_{\rm sy}$ [%]')
                    plt.subplot(122)
                    plt.plot(list_percent[:], [100*el for el in ad.solver.swarm_data_container.Trip][:], '--*', label=r'$T_{\rm rip}$ [%]')
                    plt.plot(list_percent[:], ad.solver.swarm_data_container.Ea[:], '--*', label=r'$E_a$ [deg]')
                    plt.plot(list_percent[:], [100*el for el in ad.solver.swarm_data_container.Em][:], '--*', label=r'$E_m$ [%]')
                    plt.plot(list_percent[:], [individual[-1] for individual in _swarm_data][:], '--*', label=r'$O_C$ [1]')
                    plt.xlabel(r'$\Delta d_{\rm sy}$ [%]')
                    plt.title('Q=%d, p=%d'%(Q,p))
                    plt.grid()
                    # plt.ylabel(r'')
                    plt.legend()
                    plt.show()
                    # quit()

                # Show those new evaluated designs (sensitivity analysis) in Pareto front style
                if False:
                    scatter_handle = pareto_front_plot_script(_swarm_data)
                    pareto_front_plot_color_bar_etc(scatter_handle, bool_no_limit=True)
                    plt.show()
                    quit()

    return return_tuples

import acmop, os, sys, builtins
import pandas as pd

class SwarmAnalyzer(object):
    def __init__(self, path2acmop):
        self.path2acmop = path2acmop if path2acmop[-1] == '/' else path2acmop+'/'

        self.get_folders_of_collections()

    def get_folders_of_collections(self):
        folders_of_collections = []
        for el in os.listdir(self.path2acmop):
            if 'swarm_data_collected' in el:
                "Found swarm_data_collected:", (el)
                folders_of_collections.append(el)
        self.folders_of_collections = folders_of_collections
        return folders_of_collections

    def get_swarm_group(self, folder_of_collection):
        path = self.path2acmop + folder_of_collection
        print('DEBUG: Look into', path)
        self.dict_path2SwarmDataOfTheSpecification = dict()
        self.dict_settingsOfTheSpecification = dict()
        list_specifications = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if '.json' in file:
                    # print(root)
                    normpath = os.path.normpath(root)
                    specification = normpath.split(os.sep)[-1].replace('_', ' ') # convert folder name to spec name
                    list_specifications.append(specification)
                    self.dict_path2SwarmDataOfTheSpecification[specification] = root + '/'
                    # print(specification)
                if 'settings.txt' in file:
                    with open(root+'/'+file, 'r') as f:
                        buf = f.read()
                        lst = buf.split('|')
                        specification = lst[0].strip()
                        select_fea_config_dict = lst[1].strip()
                        self.dict_settingsOfTheSpecification[specification] = select_fea_config_dict
        return list_specifications, self.dict_path2SwarmDataOfTheSpecification, self.dict_settingsOfTheSpecification

    @staticmethod
    def call_selection_criteria(ad, upper_bound_objectives, best_idx=None, proj_name=None):
        return selection_criteria(ad, 
                                ad.analyzer.swarm_data_xf, 
                                ad.analyzer.swarm_data_project_names, 
                                upper_bound_objectives=upper_bound_objectives, 
                                    best_idx=best_idx, proj_name=proj_name,
                                    # best_idx=-1, proj_name='NONE', 
                                    # best_idx=1867, proj_name='proj1868-PS-variant0-1868_IDBLIM-PS-variant0-1868', 
                                    Q=ad.spec_input_dict['Qs'], p=ad.spec_input_dict['p'])

    # @st.cache # (hash_funcs={0:get_ad})
    def get_ad_list(self, selected_specifications):
        print('[get_ad_list] for cache')
        def get_ad(specification):
            """ 醉翁之意不在酒，要的不是ad，而是ad.solver.swarm_data """
            mop = acmop.AC_Machine_Optiomization_Wrapper(
                select_fea_config_dict = self.dict_settingsOfTheSpecification[specification], 
                select_spec    = specification,
                path2SwarmData = self.dict_path2SwarmDataOfTheSpecification[specification],
                bool_show_GUI = True
            )
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Collect all swarm data from optimization results
            #~*~*~*~~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            swarm_data_file = mop.ad.read_swarm_data_json(mop.select_spec)
            number_of_chromosome = mop.ad.analyzer.number_of_chromosome
            # mop.ad.solver.output_dir = mop.path2SwarmData # tell solver where to load swarm_data.txt
            # number_of_chromosome = mop.ad.solver.read_swarm_data(specification) # ad.solver.swarm_data 在此处被赋值
            return mop.ad
        ad_list = []
        for ind, specification in enumerate(selected_specifications):

            # 读取 Swarm data
            ad = get_ad(specification)
            ad_list.append(ad)
            print(f'\t get_ad({specification})')

        return ad_list

    def inspect_swarm_and_show_table_plus_Pareto_front(self, ad_list):
        fig, ax, font = get_plot()
        df_dict = dict()

        # 调整图例legends的顺序
        if len(builtins.order)!=len(ad_list):
            print('[utility_postprocess.py] Skip re-ordering.')
            ad_list_sorted = ad_list
        else:
            ad_list_sorted = [ad_list[排名-1] for 排名 in builtins.order ]

        for ind, ad in enumerate(ad_list_sorted):
            builtins.ad = ad # 
            print(f'\n[{ind+1}]', ad.select_spec, len(ad.analyzer.swarm_data_xf),  '\n')

            Qs = ad.spec_input_dict['Qs']
            ps = ad.spec_input_dict['ps']
            p  = ad.spec_input_dict['p']
            if 'Qr' in ad.spec_input_dict:
                Qr = ad.spec_input_dict['Qr']
                label = f'p{p}Qr{Qr}Qs{Qs}ps{ps}'
            else:
                label = f'p{p}ps{ps}Qs{Qs}'
            marker = f'${ind+1}$'


            # 绘制 Pareto front
            sys.stdout = open(os.devnull, 'w')
            scatter_handle, more_info, auto_optimal_designs = self.pareto_front_plot_script(ad.analyzer.swarm_data_xf, fig, ax, marker, label, fea_config_dict=ad.fea_config_dict, z_filter=16, bool_return_more_details=True) # z_filter=20 filtered individual that has OC larger than 20
            sys.stdout = sys.__stdout__

            df_dict[ad.select_spec] = [label, len(ad.analyzer.swarm_data_xf), more_info[0][1], auto_optimal_designs[0], auto_optimal_designs[1], auto_optimal_designs[2]] # tier 1 size
            # print(more_info)
            # print(auto_optimal_designs)
            # break

        # 列出基本信息表，包括自动选择的最优个体参数，
        df = pd.DataFrame(data=df_dict, index=['label', 'archive size', 'Rank 1 PF size', 'Low Cost Design', 'High Efficiency Design', 'Low Ripple Design']).T
        # st.table(df) #  df.style.format("{:.2%}")

        # 绘制 Pareto front 的Z轴（云图色彩）
        fig = self.pareto_front_plot_color_bar_etc(scatter_handle, fig, ax, font, settings=None)
        fname = f'{os.path.dirname(__file__)}/ParetoFrontOverlapped-{builtins.folder_of_collection}.pdf'
        print('[utility_postprocess.py] save to ', fname)
        fig.savefig(fname, format='pdf', dpi=400, transparent=True) # 不能用bbox_inches='tight'，否则colorbar 会偏
        return df, fig

    def performance_table_plus_donut_chart(self, df_dict, ad, select_spec, _best_index, _best_individual_data):
        swarm_data_container = ad.swarm_data_container

        # print('\tBest:', _best_index, _best_individual_data)

        list_of_table_column = ['%.3f'%(el) for el in
                        (swarm_data_container.l_TRV[_best_index]/1000,
                         swarm_data_container.FRW[_best_index],
                         swarm_data_container.Trip[_best_index]*100,
                         swarm_data_container.Em[_best_index]*100,
                         swarm_data_container.Ea[_best_index],
                         -_best_individual_data[-2]*100, # +effciency
                         -_best_individual_data[-3],     # TRV or -cost
                          swarm_data_container.l_power_factor[0]) ]
        df_dict[select_spec] = list_of_table_column

        print('-'*len(select_spec) + '\tTRV, FRW, Trip, Em, Ea, eta, TRV, PF')
        print('\t', select_spec, ', '.join(list_of_table_column))
        print('\trated_total_loss                    :', swarm_data_container.l_rated_total_loss                    [_best_index], 'W')
        print('\trated_stator_copper_loss_along_stack:', swarm_data_container.l_rated_stator_copper_loss_along_stack[_best_index], 'W')
        print('\trated_rotor_copper_loss_along_stack :', swarm_data_container.l_rated_rotor_copper_loss_along_stack [_best_index], 'W')
        print('\tstator_copper_loss_in_end_turn      :', swarm_data_container.l_stator_copper_loss_in_end_turn      [_best_index], 'W')
        print('\trotor_copper_loss_in_end_turn       :', swarm_data_container.l_rotor_copper_loss_in_end_turn       [_best_index], 'W')
        print('\trated_iron_loss                     :', swarm_data_container.l_rated_iron_loss                     [_best_index], 'W')
        print('\trated_windage_loss                  :', swarm_data_container.l_rated_windage_loss                  [_best_index], 'W')
        print('\trated_rotor_volume                  :', swarm_data_container.RatedVol    [_best_index], 'm3')
        print('\trated_rotor_weight                  :', swarm_data_container.RatedWeight [_best_index], 'N')
        print('\trated_stack_length                  :', swarm_data_container.RatedStkLen [_best_index], 'mm')
        total_loss = swarm_data_container.l_rated_total_loss[_best_index]
        sizes = np.array( [ swarm_data_container.l_rated_iron_loss[_best_index],
                            swarm_data_container.l_rated_rotor_copper_loss_along_stack [_best_index] + swarm_data_container.l_rotor_copper_loss_in_end_turn[_best_index], # <- l_rotor_copper_loss_in_end_turn only exists for IM, there is no l_rotor_copper_loss_in_end_turn for PM motor.
                            swarm_data_container.l_rated_stator_copper_loss_along_stack[_best_index] + swarm_data_container.l_stator_copper_loss_in_end_turn[_best_index], 
                            swarm_data_container.l_rated_windage_loss[_best_index]
                          ]
                        ) / total_loss

        # print('\n'.join([el for el in dir(ad) if not el.startswith('__')]))
        # print()
        # print('\n'.join([el for el in dir(ad) if not el.startswith('__')]))
        fig = self.donut_chart(total_loss, sizes, 
                    ad.spec_input_dict['Qs'], 
                    ad.spec_input_dict['p'],
                    ad.spec_input_dict['ps'],
                    ad.spec_input_dict['Qr'] if 'Qr' in ad.spec_input_dict.keys() else None,
                    )
        return list_of_table_column, fig

    def donut_chart(self, total_loss, sizes, Qs, p, ps, Qr=None):
        print('\tValidate:', sum(sizes)*total_loss, '=', total_loss)

        #colors # https://www.schemecolor.com/color/green
        # colors = ['#C766A1','#F49762','#FFEC8A','#A1D47B', '#32D081', '#CCEFAB', '#F66867', '#F7DD7D', '#5BC5EA', '#3D93DD']
        colors = ['#C766A1','#F49762','#FFEC8A','#A1D47B']
        labels = ['Iron', 'Magnet', 'Copper', 'Windage']

        #explsion
        explode = tuple([0.025 for i in range(len(labels))]) # (0.025,0.025,0.025,0.025,0.025,0.025)

        import matplotlib as mpl
        mpl.rcParams['font.size'] = 15.0
        mpl.rcParams['font.family'] = ['Times New Roman']
        mpl.rcParams['font.family'] = 'sans-serif'
        font = {'family' : 'Times New Roman', #'serif',
            'color' : 'darkblue',
            'weight' : 'normal',
            'size' : 14,}
        fig = plt.figure()
        ax1 = plt.gca()
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5, edgecolor='white')
        ax1.text(-0.30, 0, r'Total: $%.0f$ W'%(total_loss), color='tomato', bbox=props, fontsize=20)
        patches, texts, autotexts = ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.0f%%', startangle=0, pctdistance=0.50, explode = explode, normalize=True) # MatplotlibDeprecationWarning: normalize=None does not normalize if the sum is less than 1 but this behavior is deprecated since 3.3 until two minor releases later. After the deprecation period the default value will be normalize=True. To prevent normalization pass normalize=False
        [ _.set_fontsize(20) for _ in texts]
        [ _.set_fontsize(20) for _ in autotexts]
        #draw circle
        centre_circle = plt.Circle((0,0),0.70,fc='white')
        ax1.add_artist(centre_circle)
        # Equal aspect ratio ensures that pie is drawn as a circle
        ax1.axis('equal')  

        # plt.tight_layout()
        # fig.tight_layout()
        if Qr is not None:
            figname = f'{os.path.dirname(__file__)}/Figure_loss_donut_chart_Qs{Qs}p{p}ps{ps}Qr{Qr}-{builtins.folder_of_collection}.pdf'
        else:
            figname = f'{os.path.dirname(__file__)}/Figure_loss_donut_chart_Qs{Qs}p{p}ps{ps}-{builtins.folder_of_collection}.pdf'
        fig.savefig(figname, format='pdf', dpi=600, bbox_inches='tight', pad_inches=0.0, transparent=True)
        print('\tSave donuts loss chart to', figname, '\n')
        # https://medium.com/@kvnamipara/a-better-visualisation-of-pie-charts-by-matplotlib-935b7667d77f

        # 如果转速可以变化，可以画成Stack Plots
        # https://www.youtube.com/watch?v=xN-Supd4H38
        # ax1.legend(loc=(0.05, 0.07))

        return fig

    def select_optimal_designs_manually(self, st, session_state, global_recovered_values, 
                                            ad_list, selected_specifications):

        # 再循环一次做别的事
        dict_of_list_of_table_column = dict()
        df_dict = dict()
        number_of_one_optimal_design_selected = 0
        for ind, specification in enumerate(selected_specifications):

            st.write(f'\t### [{ind}] {specification}')
            text = rf"Input upper bounds of objectives as [$O_C$, $O_B$, $O_A$] for filtering {specification}:"
            # print('The text:', text)

            if text in session_state:
                user_input_upper_bounds_4filter = st.text_input(\
                    text,  # label
                    session_state[text] # st.empty() # global_recovered_values[text] # default value 
                )
                print('\t[default             ]', session_state[text])
            else:
                user_input_upper_bounds_4filter = st.text_input(\
                    text,  # label
                    '[0,0,0]' # st.empty() # global_recovered_values[text] # default value 
                )
                print('\t[default             ]', '[0,0,0]')
            print('\t[user manually select]', user_input_upper_bounds_4filter)

            # manually save dict instead of using wrapper
            # global_recovered_values[text] = user_input_upper_bounds_4filter

            # 读取 Swarm data
            ad = ad_list[ind]
            # ad.swarm_data

            # 手动选择最优个体
            _best_index, _best_individual_data = None, None
            for ind, el in enumerate(self.call_selection_criteria(ad, eval(user_input_upper_bounds_4filter))):
                _best_index, _proj_name, _best_individual_data_reversed = el
                _best_individual_data = _best_individual_data_reversed[::-1]

                st.write('\t', el)

            if ind == 0 and _best_index is not None:
                st.write('There is only one individual left, do you want to re-produce it?')
                number_of_one_optimal_design_selected += 1

                # print(dir(ad.swarm_data_container))
                list_of_table_column, fig = self.performance_table_plus_donut_chart(df_dict, ad, specification, _best_index, _best_individual_data)
                dict_of_list_of_table_column[specification] = list_of_table_column

                # st.pyplot(fig)

        if number_of_one_optimal_design_selected == len(builtins.order ): # for writing paper (table results)

            ## 打印成latex文档直接可以用的表格形式
            print('\n\n[Performance table] ready to be copied:')

            if True:
                ## 自动顺序
                # 初始化字符串列表，作为表格的行，待添加分隔符“&”
                list_of_strings = []
                for _ in range(len(dict_of_list_of_table_column[specification])):
                    list_of_strings.append('')

                # 添加数据和分隔符（性能纵列）
                index = 0
                for specification, list_of_table_column in dict_of_list_of_table_column.items():

                    print('\t', index, specification); index += 1
                    for ind, entry in enumerate(list_of_table_column):
                        value = float(entry)
                        list_of_strings[ind] += f'{value:.1f}' + ' & '

            else:
                ## 手动顺序并添加第一列性能符号
                list_of_strings = [
                r'$\rm TRV$~[$\rm \frac{kNm}{m^3}$]  &',
                r'$\rm FRW$~[1]                      &',
                r'$T_{\rm rip}$~[\%]                 &',
                r'$E_m$~[\%]                         &',
                r'$E_a$~[deg]                        &',
                r'$\eta$~[\%]                        &',
                r'$TRV$~[$\rm USD$]                  &',
                r'Power factor [1]                   &',
                ]


                # 手动修改列表顺序
                ordered_specification_list = [
                'IM Q24p1y9 Qr16 Round Bar',
                'IM Q24p1y9 Qr14 Round Bar',
                'IM p2ps3Qs18y4 Qr30-FSW Round Bar EquivDoubleLayer',
                'IM p2ps3Qs24y5 Qr18 Round Bar EquivDoubleLayer',
                'IM Q36p3y5ps2 Qr24-ISW Round Bar',
                'IM Q36p3y5ps2 Qr20-FSW Round Bar',
                ]


                print('\t#（性能横列）打出来看看')
                print('\t', ['TRV', 'FRW', '$T_\\mathrm{rip}$', '$E_m$', '$E_a$', '$\\eta$', 'TRV', 'PF'])
                for ind, specification in enumerate(ordered_specification_list):
                    list_of_table_column = dict_of_list_of_table_column[specification]
                    print('\t', specification + ' & ' + ' & '.join([f'{float(el):.1f}' for el in list_of_table_column]))

                # 添加数据和分隔符（性能纵列）
                index = 0
                for specification in ordered_specification_list:
                    list_of_table_column = dict_of_list_of_table_column[specification]

                    print('\t', index, specification); index += 1
                    for ind, entry in enumerate(list_of_table_column):
                        value = float(entry)
                        list_of_strings[ind] += f'{value:.1f}' + ' & '

            print('\t# （性能横列）打出来看看')
            for s in list_of_strings:
                print('\t', s)

            print('\t# （性能纵列）打出来看看')
            for specification, list_of_table_column in dict_of_list_of_table_column.items():
                print(specification, end='')
                for performance in list_of_table_column:
                    print(performance, end=r' & ')
                print()

        if 'PMSM' in selected_specifications[0]:
            df = pd.DataFrame(data=df_dict, index=['TRV', 'FRW', '$T_\\mathrm{rip}$', '$E_m$', '$E_a$', '$\\eta$', 'Cost', 'disp.PF']).T
        elif 'IM' in selected_specifications[0]:
            df = pd.DataFrame(data=df_dict, index=['TRV', 'FRW', '$T_\\mathrm{rip}$', '$E_m$', '$E_a$', '$\\eta$', 'TRV', 'disp.PF']).T
        elif 'PMVM' in selected_specifications[0]:
            raise Exception('not implemented')
        return df





    ''' Below are moved from main_utility.py
    '''
    @staticmethod
    def get_sorted_swarm_data_from_the_archive(prob, popsize, path_to_archive, bool_absolute_path=False):
        output_dir_backup = ad.output_dir
        if not bool_absolute_path:
            ad.output_dir = ad.fea_config_dict['dir.parent'] + path_to_archive
        else:
            ad.output_dir = path_to_archive
        number_of_chromosome = ad.read_swarm_data(ad.bound_filter)
        if number_of_chromosome is None:
            raise Exception('Correct path?', ad.output_dir)
        ad.output_dir = output_dir_backup

        ad.flag_do_not_evaluate_when_init_pop = True
        pop = pg.population(prob, size=popsize)
        swarm_data_on_pareto_front, more_info = utility_moo.learn_about_the_archive(prob, ad.swarm_data, popsize, ad.solver.fea_config_dict, bool_plot_and_show=False, bool_more_info=True)
        ad.flag_do_not_evaluate_when_init_pop = False
        return swarm_data_on_pareto_front, more_info

    @staticmethod
    def pareto_front_plot_script(_swarm_data, fig, ax, marker, label, fea_config_dict=None, z_filter=None,
                                bool_return_more_details=False):
        import utility_moo

        from Problem_BearinglessSynchronousDesign import Problem_BearinglessSynchronousDesign
        ad.flag_do_not_evaluate_when_init_pop = True # this is very important, or else Problem_BearinglessSynchronousDesign.fitness() will invoke JMAG Designer as unexpected.
        udp = Problem_BearinglessSynchronousDesign()
        import pygmo as pg
        prob = pg.problem(udp)
        popsize = fea_config_dict["moo.popsize"]

        swarm_data_on_pareto_front, more_info = utility_moo.learn_about_the_archive(prob, _swarm_data, popsize, fea_config_dict, bool_plot_and_show=False, bool_more_info=True)
        print(len(swarm_data_on_pareto_front), len(swarm_data_on_pareto_front[0]))

        # list_of_swarm_data_on_pareto_front.append(swarm_data_on_pareto_front)

        fits = [el[-3:] for el in swarm_data_on_pareto_front]
        list_alpha_st = [el[0] for el in swarm_data_on_pareto_front]
        list_ripple_sum = [el[-1] for el in swarm_data_on_pareto_front]
        scatter_handle, auto_optimal_designs = utility_moo.my_2p5d_plot_non_dominated_fronts(fits, comp=[0,1], marker=marker, up_to_rank_no=1, ax=ax, fig=fig, no_colorbar=True, z_filter=z_filter, label=label, bool_return_auto_optimal_design=True)

        if bool_return_more_details:
            return scatter_handle, more_info, auto_optimal_designs
        else:
            return scatter_handle
        # ripple_ax = plt.figure().gca()
        # ripple_ax.plot(list_alpha_st, list_ripple_sum, 'ko')
        # ripple_ax.set_xlabel('alpha_st')
        # ripple_ax.set_ylabel('ripple sum')

    @staticmethod
    def pareto_front_plot_color_bar_etc(scatter_handle, fig, ax, font, bool_no_limit=False, settings=1):
        # color bar
        cbar_ax = fig.add_axes([0.875, 0.15, 0.02, 0.7])
        cbar_ax.get_yaxis().labelpad = 10
        clb = fig.colorbar(scatter_handle, cax=cbar_ax)
        clb.ax.set_ylabel('Ripple Performance Sum [1]', rotation=270, labelpad=14)
        if not bool_no_limit:
            if settings == 1:
                ax.legend()
                ax.set_xlim([50, 250])
                ax.set_ylim([-99, -75])

                # refernce line
                ax.plot(np.arange(50,250), np.ones(200)*-97, '--k')
                ax.text( 31, -97.3, '$-97$', fontdict=font)
                # ax.set_yticks([-75, -85, -90, -95, -97], [-75, -85, -90, -95, -97])

            elif settings == 2:
                ax.legend()
                ax.set_xlim([115, 195])
                ax.set_ylim([-98.5, -90])

                # refernce line
                ax.plot(np.arange(50,250), np.ones(200)*-97, '--k')
                ax.text( 31, -97.3, '$-97$', fontdict=font)
                # ax.set_yticks([-75, -85, -90, -95, -97], [-75, -85, -90, -95, -97])

            elif settings == 3:
                ax.legend()
                ax.set_xlim([120, 250])
                ax.set_ylim([-98.5, -90])

                ax.plot(np.arange(120,250), np.ones(130)*-97, '--k')
                # ax.text( 31, -97.3, '$-97$', fontdict=font)

            elif settings == 4:
                ax.legend(loc='upper left')
                ax.set_xlim([-40000, -15000])
                ax.set_ylim([-97.2, -95])
                ax.grid()
                ax.set_xlabel(r'$\rm {-TRV}$ [$\rm Nm/m^3$]')

            elif settings == 5:
                ax.legend(loc='best')
                # ax.set_xlim([-30000, -8000])
                # ax.set_ylim([-94, -89])
                ax.grid(True)
                ax.set_xlabel(r'$\rm {-TRV}$ [$\rm Nm/m^3$]')

            else:
                ax.legend(loc='best')
                ax.grid(True)
                ax.set_xlabel(r'$\rm {-TRV}$ [$\rm Nm/m^3$]')


        # Eric asked about non-transparent legend
        from pylab import mpl
        mpl.rcParams["legend.framealpha"] = None # default is 0.8 # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
        mpl.rcParams["legend.shadow"] = True
        ax.legend().set_zorder(555)

        # from pylab import plt
            # fig.tight_layout() # not compatable to use with color bar
            # fig.savefig(r'./Figure_Combined_Pareto_Front.eps', format='eps', dpi=1000)
            # fig.savefig(r'./Figure_Combined_Pareto_Front.png', format='png', dpi=600)
            # fig.savefig(r'./Figure_Combined_Pareto_Front.svg', format='svg', dpi=1000) # fast
        # fig.savefig(r'./Figure_Combined_Pareto_Front.pdf', format='pdf', dpi=600)
        return fig


