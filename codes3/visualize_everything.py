from pylab import np, plt
import pandas as pd
import os, json, builtins, datetime
import streamlit as st
import utility_postprocess, acmop, utility, base64, acm_designer, bearingless_spmsm_design, vernier_motor_design, bearingless_induction_design, flux_alternator_design, flux_switching_pm_design, bearingless_consequentPole_design, bearingless_VShapeconsequentPole_design, bearingless_consequentsinglePole_design
from io import BytesIO

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
        ## 选择项目路径
        path2acmop = os.path.abspath(os.path.dirname(__file__) + '\\..')
        if not os.path.exists(path2acmop + '/_default'): os.mkdir(path2acmop + '/_default')
        value = st.session_state['1.path2project'] if '1.path2project' in st.session_state.keys() else path2acmop + '/_default'
        path2project = st.text_input(label='[User] Input path2project:', value=value, on_change=None, key='2.path2project')
        if path2project[-1]!='/' or path2project[-1]!='\\': path2project += '/'

        ## 选择分析变量
        selected_sensitivity_variables = st.multiselect(
                label='Select Sensitivity Variables',
                options=[]
            )
        ## 执行灵敏度分析并展示图像
        if st.button("Run Sensitivity Analysis"):
            # fig, ax = plt.subplots()
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
                    elif 'CPPM' in self.select_spec:
                        function = bearingless_consequentPole_design.bearingless_consequentPole_template
                    elif 'VCPPM' in self.select_spec:
                        function = bearingless_VShapeconsequentPole_design.bearingless_VconsequentPole_template
                    elif 'CSPPM' in self.select_spec:
                        function = bearingless_consequentsinglePole_design.bearingless_consequentsinglePole_template
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

                def local_sensitivity_analysis(self, specify_x_denorm=None):
                            # 敏感性检查：以基本设计为准，检查不同的参数取极值时的电机性能变化！这是最简单有效的办法。七个设计参数，那么就有14种极值设计。
                    if specify_x_denorm is None:
                        # build x_denorm for the template design
                        x_denorm = self.ad.acm_template.build_x_denorm()
                    else:
                        x_denorm = specify_x_denorm
                    print('[acmop.py] x_denorm:',  x_denorm)
                    print('[acmop.py] x_denorm_dict:', self.ad.acm_template.x_denorm_dict)
                    # quit()
                    self.init_pop = []
                    if False:
                        diff = np.array(self.fea_config_dict['local_sensitivity_analysis_diff_bounds'])
                        min_b = np.array(self.fea_config_dict['local_sensitivity_analysis_min_bounds'])
                        if specify_x_denorm is None:
                            initial_design_denorm = np.array(utility.Pyrhonen_design(self.im).design_parameters_denorm )
                        else:
                            initial_design_denorm = specified_initial_design_denorm
                        initial_design = (initial_design_denorm - min_b) / diff
                        print(initial_design_denorm.tolist())
                        print(initial_design.tolist())
                        base_design = initial_design.tolist()
                        print('base_design:', base_design, '\n-------------')
                        # quit()
                        number_of_variants = self.fea_config_dict['local_sensitivity_analysis_number_of_variants']
                        self.init_pop = [initial_design] # include initial design!
                        for i in range(len(base_design)): # 10 design parameters
                            for j in range(number_of_variants+1): # 21 variants interval
                                # copy list
                                design_variant = base_design[::]
                                design_variant[i] = j * 1./number_of_variants
                                self.init_pop.append(design_variant)
                        for ind, el in enumerate(self.init_pop):
                            print(ind)
                            print(el)
                    return self.init_pop
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

            sensitivity_denorm = AC_Machine_Optiomization_Wrapper(
                select_spec = "CPPM-24s8pp-ps1-RippleRedunction", # 补充sleeve的部分以改变转矩密度过低
                select_fea_config_dict = "#02 JMAG PMSM Evaluation Setting",
                project_loc            = fr'../_ICEMS_new/',
                bool_show_GUI          = True
            ).local_sensitivity_analysis()
            for ind, el in enumerate(sensitivity_denorm):
                # sensitivity_denorm = local_sensitivity_analysis(specified_initial_design_denorm=specify_x_denorm)
                acmop.part_evaluation(specify_counter=None, specify_x_denorm=sensitivity_denorm)

                # TODO: plot
                # TODO: 1 number of x_denorm is 10 while optimization is 5
                # TODO: 2 The call of JMAG need to be fixed

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