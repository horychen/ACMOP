import sys
sys.path.insert(0, './codes3/')
import pyrhonen_procedure_as_function # main_utility, 
import utility_postprocess

from pylab import mpl, np, plt
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
plt.rc('text', usetex=True) # https://github.com/matplotlib/matplotlib/issues/4495/
plt.rc('pgf', texsystem='pdflatex')
fig, ax = plt.subplots(figsize=(8,5), constrained_layout=False)
fig.set_rasterized(True) # https://stackoverflow.com/questions/19638773/matplotlib-plots-lose-transparency-when-saving-as-ps-eps
plt.subplots_adjust(left=None, bottom=None, right=0.85, top=None, wspace=None, hspace=None)

def load_settings(select_spec, select_fea_config_dict, path2swarmData, bool_post_processing=False):
    import json, os
    # print(__file__[:-len('main_utility.py')]+'./machine_specifications.json', 'r')
    with open(os.path.dirname(__file__)+'/machine_specifications.json', 'r') as f:
        raw_specs = json.load(f)
    with open(os.path.dirname(__file__)+'/machine_simulation.json', 'r') as f:
        raw_fea_config_dicts = json.load(f)
    # quit()

    def decode_raw_specs(raw_specs, select_spec=None):
        for key, val in raw_specs.items():
            print('\n', key)
            for ke, va in val.items():
                print('\t', ke)
                for k, v in va.items():
                    print('\t\t', k + ':', v)
    # decode_raw_specs(raw_specs, select_spec)
    def decode_raw_fea_configs(raw_fea_config_dicts):
        for key, val in raw_fea_config_dicts.items():
            print('\n', key)
            for ke, va in val.items():
                print('\t', ke+':', va)
    # decode_raw_fea_configs(raw_fea_config_dicts)

    spec_input_dict = raw_specs[select_spec]['Inputs']
    fea_config_dict = raw_fea_config_dicts[select_fea_config_dict]
    fea_config_dict['bool_post_processing'] = bool_post_processing

    import where_am_i
    where_am_i.where_am_i_v2(fea_config_dict, bool_post_processing)
    fea_config_dict['run_folder'] = path2swarmData
    fea_config_dict['output_dir'] = path2swarmData

    # output_dir = fea_config_dict['dir.parent'] + fea_config_dict['run_folder'][:-1] + r'_json_files/'
    output_dir = fea_config_dict['run_folder'][:-1] + r'_json_files/'

    # create output folder only when not post-processing? No, sometimes in post-processing we run FEA simulation.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print('\t[main_utility.py]', output_dir)
    with open(output_dir+'settings.txt', 'w') as f:
        f.write(select_spec + ' | ' + select_fea_config_dict)
    # print(spec_input_dict)
    # quit()

    return output_dir, spec_input_dict, fea_config_dict

if True:
    # PEMD 2020
    specs = [   "PMSM Q06p1y2 A",
                "PMSM Q06p2y1 A",
                "PMSM Q12p1y5 A",
                "PMSM Q12p2y3 A",
                "PMSM Q12p4y1 A",
                "PMSM Q12p4y1 A",
                "PMSM Q24p1y9 A" ]
    labels = [ 'Q6p1', 'Q6p2', 'Q12p1', 'Q12p2', 'Q12p4', 'Q12p4', 'Q24p1' ] 
    if True:
        # 
        folder_of_collection = 'C:\lab\_chenjh2\ACMOP\_PEMD_2020_swarm_data_collected/'
        data_folder_names    = ['_Q06p1y2_restart_from_optimal_and_reevaluate_wo_csv/', 
                                '_Q06p2y1_restart_from_optimal_and_reevaluate_wo_csv/', 
                                '_Q12p1y5_restart_from_optimal_and_reevaluate_wo_csv/', 
                                '_Q12p2y3_restart_from_optimal_and_reevaluate_wo_csv/',
                                '_Q12p4y1_restart_from_optimal_and_reevaluate_wo_csv/',
                                '_Q12p4y1_restart_from_optimal_and_reevaluate_wo_csv_Subharmonics/',
                                '_Q24p1y9_restart_from_optimal_and_reevaluate_wo_csv/']
        markers = [ '$1$', '$2$', '$3$', '$4$', '$5$', '$8$', '$6$' ] # [ ',', '+', '.', 'o', '*' ]
        plotting_setting = 2
        select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting'
    else:
        # re-optimize since Eric found a bug in steel cost calculation
        folder_of_collection = '_re0_swarm_data_collected/'
        data_folder_names    = ['_re0_Q06p1y2/', 
                                '_re0_Q06p2y1/', 
                                '_re0_Q12p1y5/', 
                                '_re0_Q12p2y3/',
                                '_re0_Q12p4y1/',
                                '_re0_Q12p4y1_Subharmonics/',
                                '_re0_Q24p1y9/' ]
        markers = [ '$1$', '$2$', '$3$', '$4$', None, '$5$', '$6$' ]
        plotting_setting = 3
        select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting'
elif False:
    # TIA-IEMDC-ECCE 2019
    specs  = [ "IM Q24p1y9 A", "IM Q24p1y9 Qr32", "IM Q24p2y6 Qr16", "IM Q24p2y6 Qr32"] 
    labels = [ '$p=1,Q_r=16$',    '$p=1,Q_r=32$',    '$p=2,Q_r=16$',    '$p=2,Q_r=32$'] 
    # optimize 2 pole induction motor for tia-iemdc-ecce paper
    folder_of_collection = '_TIA_IEMDC_swarm_data_collected/'
    data_folder_names    = ['Q24p1y9/', 
                            'Q24p1y9/', 
                            'Q24p2y6/',
                            'Q24p2y6/']
    markers = [ '$1$', '$2$', '$3$', '$4$' ]
    plotting_setting = 4
    select_fea_config_dict = '#011 JMAG IM Re-evaluation wo/ CSV Setting'
else:
    # TIA-ISMB 2020
    folder_of_collection = '_TIA_ISMB_swarm_data_collected/'
    specs  = [  "IM Q24p1y9 Qr14 Round Bar", 
                "IM Q24p1y9 Qr16 Round Bar", 
                "IM Q36p3y5ps2 Qr20-FSW Round Bar Separate Winding", 
                "IM Q36p3y5ps2 Qr24-ISW Round Bar Separate Winding",
                "IM Q36p3y5ps2 Qr20-FSW Round Bar", 
                "IM Q36p3y5ps2 Qr24-ISW Round Bar"]
    labels = [  'Qs24ps2p1Qr14',
                'Qs24ps2p1Qr16',
                'Qs36ps2p3Qr20(Sepa)',
                'Qs36ps2p3Qr24(Sepa)',
                'Qs36ps2p3Qr20(Comb)',
                'Qs36ps2p3Qr24(Comb)']
    data_folder_names    = ['DataFolder_p1Qs24ps2Qr14/', 
                            'DataFolder_p1Qs24ps2Qr16/', 
                            'Single-Layer-FSW-Separate-Winding/', 
                            'Single-Layer-ISW-Separate-Winding/',
                            'Single-Layer-FSW/', 
                            'Single-Layer-ISW/']
    markers = [ '$1$', '$2$', '$3$', '$4$', '$5$', '$6$']
    plotting_setting = None
    select_fea_config_dict = "#019 JMAG IM Nine Variables"


for select_spec, data_folder_name, marker, label in zip(specs,
                                                        data_folder_names,
                                                        markers,
                                                        labels
                                                        ):

    if marker is None:
        continue

    path_to_archive = folder_of_collection + data_folder_name + select_spec.replace(' ', '_') + '/'

    output_dir, spec_input_dict, fea_config_dict = load_settings(select_spec, select_fea_config_dict, folder_of_collection+data_folder_name, bool_post_processing = True)

    spec = pyrhonen_procedure_as_function.desgin_specification(**spec_input_dict)
    if 'SM' in select_spec:
        spec.acm_template = spec.build_pmsm_template(fea_config_dict, spec_input_dict, im_template=None)
    elif 'IM' in select_spec:
        # This template generation code is copied from main.py, and it needs some revision to reduce the codes.
        print(spec.build_name())
        spec.bool_bad_specifications = spec.pyrhonen_procedure()
        for k,v in spec.spec_geometry_dict.items():
            print(k+':', v)
        print(spec.build_name()) # rebuild for new guess air gap flux density # TODO：自动修正转子电流密度的设置值？
        import population
        spec.acm_template = population.bearingless_induction_motor_design(spec_input_dict, spec.spec_derive_dict, spec.spec_geometry_dict, fea_config_dict)
    else:
        raise

    import acm_designer
    global ad
    ad = acm_designer.acm_designer(select_spec, spec_input_dict, select_fea_config_dict, fea_config_dict, acm_template=spec.acm_template)

    # dummy population
    ad.flag_do_not_evaluate_when_init_pop = True

    if 'IM' in select_spec:
        ad.bounds_denorm = spec.get_im_classic_bounds(which_filter=fea_config_dict['which_filter'])
        ad.bound_filter  = spec.bound_filter
        otnb = spec.original_template_neighbor_bounds
    elif 'PMSM' in select_spec:
        ad.bounds_denorm = spec.acm_template.bounds_denorm
        # ad.bound_filter  = spec.bound_filter
        otnb = spec.acm_template.original_template_neighbor_bounds

    ad.counter_fitness_called = 99999
    ad.counter_fitness_return = 99999

    __builtins__.ad = ad # share global variable between modules # https://stackoverflow.com/questions/142545/how-to-make-a-cross-module-variable


    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Collect all data 
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    path_to_target = fea_config_dict['dir.parent'] + folder_of_collection + data_folder_name+select_spec.replace(' ', '_') + '/'  # '_Q24p1y9_restart_from_optimal_and_reevaluate_wo_csv/PMSM_Q24p1y9_A/'

    ad.analyzer.output_dir = path_to_target
    number_of_chromosome = ad.read_swarm_data(select_spec)
    print('number_of_chromosome =', number_of_chromosome)

    # Set the output_dir back right away
    # ad.analyzer.output_dir = ad.analyzer.fea_config_dict['dir.parent'] + ad.analyzer.fea_config_dict['run_folder']

    _swarm_data          = ad.swarm_data
    _swarm_project_names = ad.swarm_data_container.project_names

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # plot the Pareto front for the archive
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    if True:
        scatter_handle = utility_postprocess.pareto_front_plot_script(ad.swarm_data, fig, ax, marker, label, 
                                                                      fea_config_dict=fea_config_dict, z_filter=60) # z_filter=20 filtered individual that has OC larger than 20
    # break

# plt.show()
# quit()

fig = utility_postprocess.pareto_front_plot_color_bar_etc(scatter_handle, fig, ax, font, settings=plotting_setting)

# plt.gca().set_axis_off()
# plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
#             hspace = 0, wspace = 0)
# plt.margins(0,0)

fig.savefig(fea_config_dict['dir.parent'] + folder_of_collection + r'Figure_Combined_Pareto_Front.pdf', format='pdf', dpi=600)
plt.show()


