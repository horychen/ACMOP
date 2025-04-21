from pylab import np, plt
import pandas as pd
import os, json, builtins, datetime, re
import streamlit as st
import utility_postprocess, acmop, utility, base64
# import acm_designer
from io import BytesIO
import pygmo as pg
# import Problem_BearinglessSynchronousDesign

def basic_information_about_optimization():
    # dimension of x and f? 
    # pop size?
    # algorithm?
    pass

def displayPDF(file, width=1800, height=500):
    # Opening file from file path
    with open(file, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')

    # Embedding PDF in HTML
    # pdf_display = F'<embed src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf">'
    pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}#view=FitH" width="{width}" height="{height}" type="application/pdf"></iframe>' # https://discuss.streamlit.io/t/rendering-pdf-on-ui/13505

    # Displaying File
    st.markdown(pdf_display, unsafe_allow_html=True)

def displayPDF_side_by_side(list_files, width=400, height=450):
    list_base64_pdf = []
    for file in list_files:
        with open(file, "rb") as f:
            base64_pdf = base64.b64encode(f.read()).decode('utf-8')
            list_base64_pdf.append(base64_pdf)

    # Embedding PDF in HTML
    # 添加 #view=FitH 可以最大化显示pdf，参考：https://stackoverflow.com/questions/20562543/zoom-to-fit-pdf-embedded-in-html
    pdf_display = ''
    for base64_pdf in list_base64_pdf:
        pdf_display += F'''
        <div class="box">
            <iframe src="data:application/pdf;base64,{base64_pdf}#view=FitH"
                        frameborder="0" 
                        scrolling="no" 
                        width="{width}"
                        height="{height}"
                        align="left"
                        type="application/pdf"> 
            </iframe> 
        </div>
        '''
        # <div class="box"><iframe src="data:application/pdf;base64,{base64_pdf}"
        #                         frameborder="0" 
        #                         scrolling="no" 
        #                         width={width}
        #                         height={height}
        #                         align="left">
        #                         type="application/pdf"
        #                         </iframe>
        # </div>
        # '''
    # Displaying File
    st.markdown(pdf_display, unsafe_allow_html=True)

def pyplot_width(fig):
    # with st.echo():
    buf = BytesIO()
    fig.savefig(buf, format="png")
    image_width = st.number_input("Width of the loss breakdown donut plot", 1, 2000, 700)
    use_column_width = st.checkbox("Use column width?")
    st.image(buf, width=int(image_width), use_column_width=use_column_width)


# from emy-c
def get_user_config():
    history = {}
    fname_session_state = f'{os.path.dirname(__file__)}/visualize_streamlit_user_session_data.json'
    if not os.path.exists(fname_session_state):
        with open(fname_session_state, 'w') as f:
            f.write('{\n}')
    with open(fname_session_state, 'r') as f:
        d = json.load(f)
        for k, v in d.items():
            history[k] = v
    return history


# from emy-c
if __name__ == '__main__':

    with st.sidebar:
        st.markdown(
        """
        <style>
        [data-testid="stSidebar"][aria-expanded="true"]{
            min-width: 450px;
            max-width: 1450px;
        }
        """,
        unsafe_allow_html=True,
    )

    if 'look' not in st.session_state:
        st.session_state.look = []

    ## Session State as hisotry    
    history = get_user_config()
    
    print(f'\n\n\n{history=}')
    
    ## 标题
    st.title(f'ACMOP Visualization {datetime.date.today()}')
    
    with st.sidebar:
        st.title('User Inputs')
        st.sidebar.header('User Inputs')


    tab1, tab2 = st.tabs(["Optimization", "Select"])

    with tab1:
        st.subheader('Optimal Design')
        # st.pyplot(fig)

        ## 选择项目路径
        path2acmop = os.path.abspath(os.path.dirname(__file__) + '\\..')
        if not os.path.exists(path2acmop + '/_default'): os.mkdir(path2acmop + '/_default')
        value = st.session_state['1.path2project'] if '1.path2project' in st.session_state.keys() else path2acmop + '/_default'
        path2project = st.text_input(label='[User] Input path2project:', value=value, on_change=None, key='1.path2project')
        if path2project[-1]!='/' or path2project[-1]!='\\': path2project += '/'

        ## 选择电机规格
        _, list_specifications, _ = next(os.walk(path2project))
        selected_specifications = st.multiselect(label="[User] Select folder(s):", options=list_specifications, default=None, key='2.selected_specifications')

        ## 按照所选的电机规格，显示用户输入信息
        swarm_dict = {}
        if selected_specifications == []:
            st.error("Please select at least one specification.")
        else:
            for folder in selected_specifications:
                with open(path2project+folder+'/acmop-settings.txt', 'r') as f:
                    buf = f.read()
                    lst = buf.split('|')
                    select_spec = lst[0].strip()
                    select_fea_config_dict = lst[1].strip()
                utility.blockPrint()
                swarm_dict[folder] = mop = acmop.AC_Machine_Optiomization_Wrapper(select_fea_config_dict, select_spec, project_loc=path2project)
                utility.enablePrint()

            ## 侧边栏 Sidebar
            user_selected_folder = st.sidebar.selectbox('Select a folder to show its inputs', selected_specifications, key='3.user_selected_folder')

            ## The current user selected MOP object
            print(f'{user_selected_folder=}')
            mop = swarm_dict[user_selected_folder]
            
            # Show user selected mop's inputs
            st.sidebar.table(pd.DataFrame(data=list(mop.spec_input_dict.values()), index=list(mop.spec_input_dict.keys()), dtype="string", columns=['Value',]))
            st.sidebar.table(pd.DataFrame(data=list(mop.fea_config_dict.values()), index=list(mop.fea_config_dict.keys()), dtype="string", columns=['Value',]))

            ## 集群基本信息
            st.write('# 1. Swarms Information')
            ## 是否显示自动最优个体表格和帕累托前沿？
            optimal_xf_dict = None
            if st.checkbox("Show Swarm Table and Pareto Front?"):
                df_swarm, fig_Pareto, optimal_fitness_dict, optimal_xf_dict = utility_postprocess.inspect_swarm_and_show_table_plus_Pareto_front(swarm_dict, output_dir=path2project)
                st.table(df_swarm)
                st.pyplot(fig_Pareto)

            st.write('# 2. Template/Initial Design Information')
            wily_fname = 'wily_p%dps%dQ%dy%d'%(mop.spec_input_dict['p'], mop.spec_input_dict['ps'], mop.spec_input_dict['Qs'], mop.spec_input_dict['coil_pitch_y'])
            st.write('## 2.1. Winding Information of', wily_fname)
            try:
                displayPDF(f'{path2acmop}_wily/{wily_fname}.pdf')
            except FileNotFoundError:
                st.write('The winding derivation file is absent')

            st.write('## 2.2. Cross Section (Initial and Optimal Designs)')
            cairo_fname           = mop.fea_config_dict['output_dir'] + 'indCairo.pdf'
            cairo_fname_optimal_1 = mop.fea_config_dict['output_dir'] + 'indCairoOptimal1.pdf'
            cairo_fname_optimal_2 = mop.fea_config_dict['output_dir'] + 'indCairoOptimal2.pdf'
            cairo_fname_optimal_3 = mop.fea_config_dict['output_dir'] + 'indCairoOptimal3.pdf'
            if not os.path.exists(cairo_fname):     mop.part_evaluation_geometry()

            # Show user selected mop's auto optimal designs
            if optimal_xf_dict is not None:
                # auto_optimal_designs_fitnesses = optimal_fitness_dict[mop.ad.select_spec] # obsolete
                auto_optimal_designs_xf        = optimal_xf_dict[mop.ad.select_spec]
                if auto_optimal_designs_xf[0] !=[]: mop.part_evaluation_geometry(auto_optimal_designs_xf[0], counter='CairoOptimal1') # and not os.path.exists(cairo_fname_optimal_1)
                if auto_optimal_designs_xf[1] !=[]: mop.part_evaluation_geometry(auto_optimal_designs_xf[1], counter='CairoOptimal2') # and not os.path.exists(cairo_fname_optimal_2)
                if auto_optimal_designs_xf[2] !=[]: mop.part_evaluation_geometry(auto_optimal_designs_xf[2], counter='CairoOptimal3') # and not os.path.exists(cairo_fname_optimal_3)
                st.write('### 2.2.1 Initial Design:')
                # displayPDF(cairo_fname, width=500, height=500)

                st.write('### 2.2.2 Auto Optimal Design minimum OA:')
                st.write(str(auto_optimal_designs_xf[0]))
                # displayPDF(cairo_fname_optimal_1, width=500, height=500)

                st.write('### 2.2.3 Auto Optimal Design minimum OB:')
                st.write(str(auto_optimal_designs_xf[1]))
                # displayPDF(cairo_fname_optimal_2, width=500, height=500)

                st.write('### 2.2.4 Auto Optimal Design minimum OC:')
                st.write(str(auto_optimal_designs_xf[2]))
                # displayPDF(cairo_fname_optimal_3, width=500, height=500)

                displayPDF_side_by_side([cairo_fname, cairo_fname_optimal_1, cairo_fname_optimal_2, cairo_fname_optimal_3])


            ## 根据用户在 text_input 的输入来筛选符合条件的最优个体
            if True:
                print(f'{selected_specifications=}')
                def select_optimal_designs_manually(selected_specifications):
                    # 再遍历selected_specifications一次，做别的事
                    dict_of_list_of_table_column = dict()
                    number_of_one_optimal_design_selected = 0
                    for ind, folder in enumerate(selected_specifications):

                        st.write(f'\t### [{ind}] {folder}')
                        user_input_upper_bounds_4filter = st.text_input(label=rf"Input upper bounds of objectives as [$O_C$, $O_B$, $O_A$] for filtering {folder}:", 
                            value='[20, -0.8, 200]', 
                            key=f'4.user_input_upper_bounds_4filter:{folder}'
                        )

                        # 获取ad
                        mop = swarm_dict[folder]
                        ad = mop.ad

                        # 手动选择最优个体
                        _best_index, _best_individual_data, = None, None
                        for ind, el in enumerate(utility_postprocess.call_selection_criteria(ad, eval(user_input_upper_bounds_4filter))):
                            if el is None:
                                raise Exception(str(el))
                            else:
                                _best_index, _proj_name, _best_individual_data_reversed = el
                                _best_individual_data = _best_individual_data_reversed[::-1]

                            st.write(F'\t{el[0]}, {el[1]}, f3={el[-1][0]:.1f}, f2={el[-1][1]:.4f}, f1={el[-1][2]:.1f}, ' + ', '.join(F'{x:.2f}' for x in el[-1][3:]))

                        if ind == 0 and _best_index is not None:
                            if st.checkbox('There is only one individual left, do you want to re-produce it?'):
                                mop.part_evaluation_geometry(_best_individual_data, counter='UserSelectedOptimal')
                                displayPDF_side_by_side([mop.fea_config_dict['output_dir'] + 'indUserSelectedOptimal.pdf', ])
                            number_of_one_optimal_design_selected += 1

                            # print(dir(ad.swarm_data_container))
                            list_of_table_column, fig_donut = utility_postprocess.performance_table_plus_donut_chart(ad, folder, _best_index, _best_individual_data, output_dir=path2project)
                            dict_of_list_of_table_column[folder] = list_of_table_column

                            # for writing paper (table results)

                            ## 打印成latex文档直接可以用的表格形式
                            # print('\n\n[Performance table] ready to be copied:')

                            if True:
                                ## 自动顺序
                                # 初始化字符串列表，作为表格的行，待添加分隔符“&”
                                list_of_strings = []
                                for _ in range(len(dict_of_list_of_table_column[folder])):
                                    list_of_strings.append('')

                                # 添加数据和分隔符（性能纵列）
                                index = 0
                                for _folder, list_of_table_column in dict_of_list_of_table_column.items():

                                    # print('\t', index, _folder); index += 1
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

                            # print('\t# （性能横列）打出来看看')
                            # for s in list_of_strings:
                            #     print('\t', s)

                            # print('\t# （性能纵列）打出来看看')
                            # for specification, list_of_table_column in dict_of_list_of_table_column.items():
                            #     print(specification, end='')
                            #     for performance in list_of_table_column:
                            #         print(performance, end=r' & ')
                            #     print()

                            df_performancce = pd.DataFrame(data=dict_of_list_of_table_column, index=['TRV', 'FRW', '$T_\\mathrm{rip}$', '$E_m$', '$E_a$', '$\\eta$', 'Cost', 'disp.PF']).T
                            return df_performancce, fig_donut

                # [1.1, -0.93, 170]
                _ = select_optimal_designs_manually(selected_specifications)
                if _ is not None:
                    df_performancce, fig_donut = _
                    st.table(df_performancce)
                    pyplot_width(fig_donut)
                    # mop.ad.acm_variant.analyzer.load_time_domain_data(_best_index) # TODO

                # save user input filters as json file
                with open(f'{os.path.dirname(__file__)}/streamlit_user_session_data.json', 'w') as f:
                    json.dump(dict(st.session_state), f, ensure_ascii=False, indent=4)

        ## 结束 (Streamlit widgets automatically run the script from top to bottom. Since this button is not connected to any other logic, it just causes a plain rerun.)
        # st.button("Re-run")


    with tab2:
        st.subheader('Sensitivity Analysis')
        # st.pyplot(fig)

        def local_sensitivity_analysis(specify_x_denorm=None, x_denorm_bounds = None, number_of_steps = 10):
                    # 敏感性检查：以基本设计为准，检查不同的参数取极值时的电机性能变化！这是最简单有效的办法。七个设计参数，那么就有14种极值设计。
            if specify_x_denorm is None:
                # build x_denorm for the template design
                raise Exception('specify_x_denorm is None')
            else:
                x_denorm = specify_x_denorm
                x_denorm_bounds = x_denorm_bounds
            print('[acmop.py] x_denorm:',  x_denorm)
            print('[acmop.py] x_denorm_bounds:',  x_denorm_bounds)
            # print('[acmop.py] x_denorm_dict:', x_denorm_dict)
            # quit()
            pop_x_values = []
            bounds_min = []
            bounds_max = []
            diff = []
            if True:
                for i, bounds in enumerate(x_denorm_bounds):
                    bounds_min.append(bounds[0])
                    bounds_max.append(bounds[1])
                    diff.append(bounds[1] - bounds[0])

                print(f'[acmop.py] bounds_min: {bounds_min}, bounds_max: {bounds_max}') 
                print(f'[acmop.py] diff: {diff}')
                # quit()
                # if specify_x_denorm is None:
                    # initial_design_denorm = np.array(utility.Pyrhonen_design(pmsm).design_parameters_denorm)
                # else:
                    # initial_design_denorm = specified_initial_design_denorm
                initial_design = x_denorm
                number_of_steps = number_of_steps
                pop_x_values = [initial_design] # include initial design!
                # print(len(initial_design))
                # quit()
                for i in range(len(initial_design)): # 10 design parameters
                    for j in range(number_of_steps+1): # 21 variants interval
                        # copy list
                        design_variant = initial_design[::]
                        design_variant[i] = j * diff[i] / number_of_steps + bounds_min[i]
                        pop_x_values.append(design_variant)
                # for ind, el in enumerate(pop_x_values):
                    # print(ind)
                    # print(el)
            return pop_x_values
            
        ## 选择项目路径
        path2acmop = os.path.abspath(os.path.dirname(__file__) + '\\..')
        if not os.path.exists(path2acmop + '/_default'): os.mkdir(path2acmop + '/_default')
        value = st.session_state['1.path2project'] if '1.path2project' in st.session_state.keys() else path2acmop + '/_default'
        path2project = st.text_input(label='[User] Input path2project:', value=value, on_change=None, key='2.path2project')
        if path2project[-1]!='/' or path2project[-1]!='\\': path2project += '/'

        ## 选择电机规格
        _, list_specifications, _ = next(os.walk(path2project))
        # print(list_specifications)
        # quit()
        selected_specifications_2 = st.multiselect(label="[User] Select folder(s):", options=list_specifications, default=None, key='2.selected_specifications_2')

        ## create swarm_dict
        swarm_dict = {}
        if selected_specifications_2 == []:
            st.error("Please select at least one specification.")
        else:
            for folder in selected_specifications_2:
                # print(folder)
                # quit()
                with open(path2project+folder+'/acmop-settings.txt', 'r') as f:
                    buf = f.read()
                    lst = buf.split('|')
                    select_spec = lst[0].strip()
                    select_fea_config_dict = lst[1].strip()
                # quit()
                json_file_name = folder.replace('_', ' ') + '.json'
                json_file_path = path2project + folder + '/' + json_file_name
                # print(json_file_path)
                # quit()
                try:
                    # 读取 JSON 文件
                    with open(json_file_path, 'r') as jf:
                        json_content = jf.read()  # 加载 JSON 数据
                        start_index = json_content.find('{')
                        if start_index != -1:
                            pseudo_json_data = json_content[start_index:]
                        # print(pseudo_json_data)
                        # quit()
                except FileNotFoundError:
                    st.error(f"JSON file '{json_file_name}' not found in folder '{folder}'.")
                    continue
                except json.JSONDecodeError:
                    st.error(f"Error decoding JSON file '{json_file_name}' in folder '{folder}'.")
                    continue
                utility.blockPrint()
                swarm_dict[folder] = mop = acmop.AC_Machine_Optiomization_Wrapper(select_fea_config_dict, select_spec, project_loc=path2project)
                utility.enablePrint()

                json_data = json.loads(pseudo_json_data)
                for key in json_data:
                    if isinstance(json_data[key], dict) and 'Geometric parameters' in json_data[key]: # 直接找最顶级的键值
                        geometric_parameters = json_data[key]['Geometric parameters']
                        break

                # print(geometric_parameters)
                # print(GP)
                # quit()

                design_variable_name = []
                param_type = []
                param_name = []
                param_value = []
                param_bounds = []
                param_calc = []

                for param in geometric_parameters:
                    for key, value in param.items():
                        if value.get("type") == "free":
                            bounds = value.get("bounds")
                            if bounds != [None, None]:
                                design_variable_name.append(key)
                                param_type.append(value.get("type", "N/A"))
                                param_name.append(value.get("name", "N/A"))
                                param_value.append(value.get("value", "N/A"))
                                param_bounds.append(value.get("bounds", "N/A"))
                                param_calc.append(value.get("calc", "N/A"))
                number_of_variants = len(design_variable_name)
                print(f"There are {number_of_variants} free parameters:")
                print(f"Name: {param_name}, Type: {param_type}, Value: {param_value}, Bounds: {param_bounds}, Calc: {param_calc}")
                # quit()
                
            
            # print(f'mop')
            # quit()
        
        ## 之前选择的电机规格中 选择一个较优的设计 然后进行灵敏度分析（还需要一个选择
        ## 读json 拿到 x f 准备re-evaluation            done
        ## re-evaluation 准备一个新的文件夹吧要不
        ## 一个参数一个文件夹 然后跑20个变化
        ## 画图
        

        # ## 选择分析变量
        selected_sensitivity_variables = st.multiselect(
                label='Select Sensitivity Variables',
                options = design_variable_name # TODO: 需要读swarm_dict 找到free的变量 拉出来表格选择 目前可以先不管
            )
        
        pop_x_values = local_sensitivity_analysis(param_value, param_bounds, number_of_steps = 20)
        print(pop_x_values)
        # quit()
        
        # pop.set_x(pop_x_values)


        # ## 执行灵敏度分析并展示图像
        # if st.button("Run Sensitivity Analysis"):
        #     # fig, ax = plt.subplots()
                        
            
            
        #     for ind, el in enumerate(sensitivity_denorm):
        #         # sensitivity_denorm = local_sensitivity_analysis(specified_initial_design_denorm=specify_x_denorm)
        #         acmop.part_evaluation(specify_counter=None, specify_x_denorm=sensitivity_denorm)

        #         # TODO: plot
        #         # TODO: 1 number of x_denorm is 10 while optimization is 5
        #         # TODO: 2 The call of JMAG need to be fixed # done

            # sw = population.swarm(fea_config_dict, de_config_dict=de_config_dict)
            # sw, spec = utility.load_data(path2project) # TODO: initialize sw and spec
            # analyzer = utility.SwarmDataAnalyzer(sw, spec, dir_run=path2project, run_integer=1)
            # analyzer.sensitivity_bar_charts(ax)

            # st.pyplot(fig)


# if st.session_state.user_selected_motor in history and 'd_user_input_motor_dict' in history[st.session_state.user_selected_motor]:
#     st.session_state.d_user_input_motor_dict = history[st.session_state.user_selected_motor]['d_user_input_motor_dict']
# else:
#     st.session_state.d_user_input_motor_dict = {
#         # Timing
#         'CL_TS': 1e-4,
