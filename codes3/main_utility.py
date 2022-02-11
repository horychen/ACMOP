from pylab import np
import pygmo as pg
import utility_moo
import json, os
def load_settings(select_spec, select_fea_config_dict, project_loc=None, path2SwarmData=None, bool_post_processing=False):
    # print('-'*40, '[main_utility.py] load_settings()')
    __file__dirname_as_in_python39 = os.path.dirname(os.path.abspath(__file__))
    # print(__file__dirname_as_in_python39)

    with open((__file__dirname_as_in_python39)+'/machine_specifications.json', 'r') as f:
        raw_specs = json.load(f)
    with open((__file__dirname_as_in_python39)+'/machine_simulation.json', 'r') as f:
        raw_fea_config_dicts = json.load(f)

    # def decode_raw_specs(raw_specs, select_spec=None):
    #     for key, val in raw_specs.items():
    #         print('\n', key)
    #         for ke, va in val.items():
    #             print('\t', ke)
    #             for k, v in va.items():
    #                 print('\t\t', k + ':', v)
    # decode_raw_specs(raw_specs, select_spec)
    # def decode_raw_fea_configs(raw_fea_config_dicts):
    #     for key, val in raw_fea_config_dicts.items():
    #         print('\n', key)
    #         for ke, va in val.items():
    #             print('\t', ke+':', va)
    # decode_raw_fea_configs(raw_fea_config_dicts)

    spec_input_dict = raw_specs[select_spec]['Inputs']
    fea_config_dict = raw_fea_config_dicts[select_fea_config_dict]
    fea_config_dict['bool_post_processing'] = bool_post_processing

    import where_am_i
    where_am_i.where_am_i_v2(fea_config_dict, bool_post_processing)
    if path2SwarmData is None:
        path2SwarmData = project_loc + select_spec.replace(' ', '_')+'/'
    if project_loc is None:
        project_loc = os.path.abspath(os.path.join(path2SwarmData, '..',))

    fea_config_dict['output_dir'] = path2SwarmData

    output_dir = fea_config_dict['output_dir'] #[:-1] + r'_json_files/'

    # create output folder only when not post-processing? No, sometimes in post-processing we run FEA simulation.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # print('[main_utility.py] output_dir:', output_dir)
    with open(output_dir+'acmop-settings.txt', 'w') as f:
        f.write(select_spec + ' | ' + select_fea_config_dict)
    # print(spec_input_dict)
    # quit()

    return spec_input_dict, fea_config_dict

def get_sorted_swarm_data_from_the_archive(prob, popsize, path_to_archive, bool_absolute_path=False):
    output_dir_backup = ad.solver.output_dir
    if not bool_absolute_path:
        ad.solver.output_dir = ad.solver.fea_config_dict['dir.parent'] + path_to_archive
    else:
        ad.solver.output_dir = path_to_archive
    number_of_chromosome = ad.solver.read_swarm_data(ad.bound_filter)
    if number_of_chromosome is None:
        raise Exception('Correct path?', ad.solver.output_dir)
    ad.solver.output_dir = output_dir_backup

    ad.flag_do_not_evaluate_when_init_pop = True
    pop = pg.population(prob, size=popsize)
    swarm_data_on_pareto_front, more_info = utility_moo.learn_about_the_archive(prob, ad.solver.swarm_data, popsize, ad.solver.fea_config_dict, bool_plot_and_show=False, bool_more_info=True)
    ad.flag_do_not_evaluate_when_init_pop = False
    return swarm_data_on_pareto_front, more_info


def pareto_front_plot_script(_swarm_data, fig, ax, marker, label, fea_config_dict=None, z_filter=None,
                            bool_return_more_details=False):
    import utility_moo

    from Problem_BearinglessSynchronousDesign import Problem_BearinglessSynchronousDesign
    ad.flag_do_not_evaluate_when_init_pop = True # this is very important, or else Problem_BearinglessSynchronousDesign.fitness() will invoke JMAG Designer as unexpected.
    udp = Problem_BearinglessSynchronousDesign()
    import pygmo as pg
    prob = pg.problem(udp)
    popsize = 78

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


