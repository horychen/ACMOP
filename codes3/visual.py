from pylab import np, plt
import pandas as pd
import os
import json
import builtins
import datetime
import streamlit as st
import utility_postprocess
import acmop
import utility
import base64
import acm_designer
import population
# import acm_designer
from io import BytesIO


def displayPDF(file, width=1800, height=500):
    with open(file, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
    # https://discuss.streamlit.io/t/rendering-pdf-on-ui/13505
    pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}#view=FitH" width="{width}" height="{height}" type="application/pdf"></iframe>'
    st.markdown(pdf_display, unsafe_allow_html=True)


def displayPDF_side_by_side(list_files, width=400, height=450):
    list_base64_pdf = []
    for file in list_files:
        with open(file, "rb") as f:
            base64_pdf = base64.b64encode(f.read()).decode('utf-8')
            list_base64_pdf.append(base64_pdf)
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
    st.markdown(pdf_display, unsafe_allow_html=True)


def pyplot_width(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png")
    image_width = st.number_input(
        "Width of the loss breakdown donut plot", 1, 2000, 700)
    use_column_width = st.checkbox("Use column width?")
    st.image(buf, width=int(image_width), use_column_width=use_column_width)



def OptimizationSetupContent(user_selected_folder):
    mop = swarm_dict[user_selected_folder]
    # mop


    col1, col2 = st.columns(2)

    # 左栏内容
    with col1:
        st.header("左栏")
        st.subheader(f'Optimization Setup of {user_selected_folder}')
        mop.ad.acm_template.d['GP']

        st.write('#### Decision Variables:')
        convert_to_dict = {key: (val.type, val.value, val.bounds[0], val.bounds[1]) for key, val in mop.ad.acm_template.d['GP'].items() if 'free' in val.type}
        convert_to_df = pd.DataFrame(convert_to_dict)
        st.table(convert_to_df.T)
        st.write('#### Derived Variables:')
        convert_to_dict = {key: (val.type, val.value) for key, val in mop.ad.acm_template.d['GP'].items() if 'derived' in val.type}
        convert_to_df = pd.DataFrame(convert_to_dict)
        st.table(convert_to_df.T)
        st.write('#### Fixed Variables:')
        convert_to_dict = {key: (val.type, val.value) for key, val in mop.ad.acm_template.d['GP'].items() if 'fixed' in val.type}
        convert_to_df = pd.DataFrame(convert_to_dict)
        st.table(convert_to_df.T)

        st.write('#### Objectives (Non-dominated Sorting)')
        st.write(f'''
                    {mop.fea_config_dict['moo.fitness_OA']=}
                    {mop.fea_config_dict['moo.fitness_OB']=}
                    {mop.fea_config_dict['moo.fitness_OC']=}
        ''')

    # 右栏内容
    with col2:
        st.header("右栏")
        mop.part_evaluation_geometry()
        displayPDF_side_by_side(
            [mop.fea_config_dict['output_dir'] + 'indCairo.pdf', mop.fea_config_dict['output_dir'] + 'indCairoOptimal1.pdf', mop.fea_config_dict['output_dir'] + 'indCairoOptimal2.pdf', mop.fea_config_dict['output_dir'] + 'indCairoOptimal3.pdf'],
            width=300,  height=300)


def PopulationContent():
    # 集群基本信息
    st.write('# 1. Swarms Information')
    # 是否显示自动最优个体表格和帕累托前沿？
    optimal_xf_dict = None
    if st.checkbox("Show Swarm Table and Pareto Front?"):
        # , key='Input z_filter for Pareto front'))
        z_filter = float(st.text_input(
            label='Input z-filter for Pareto front:', value='20'))
        df_swarm, fig_Pareto, optimal_fitness_dict, optimal_xf_dict = utility_postprocess.inspect_swarm_and_show_table_plus_Pareto_front(
            swarm_dict, z_filter, output_dir=path2project)
        st.pyplot(fig_Pareto)
        st.table(df_swarm)

    st.write('# 2. Template/Initial Design Information')
    wily_fname = 'wily_p%dps%dQ%dy%d' % (
        mop.spec_input_dict['p'], mop.spec_input_dict['ps'], mop.spec_input_dict['Qs'], mop.spec_input_dict['coil_pitch_y'])
    st.write('## 2.1. Winding Information of', wily_fname)
    try:
        displayPDF(f'{path2acmop}_wily/{wily_fname}.pdf')
    except FileNotFoundError:
        st.write('The winding derivation file is absent')

    st.write('## 2.2. Cross Section (Initial and Optimal Designs)')
    cairo_fname = mop.fea_config_dict['output_dir'] + 'indCairo.pdf'
    cairo_fname_optimal_1 = mop.fea_config_dict['output_dir'] + \
        'indCairoOptimal1.pdf'
    cairo_fname_optimal_2 = mop.fea_config_dict['output_dir'] + \
        'indCairoOptimal2.pdf'
    cairo_fname_optimal_3 = mop.fea_config_dict['output_dir'] + \
        'indCairoOptimal3.pdf'
    if not os.path.exists(cairo_fname):
        mop.part_evaluation_geometry()

    # Show user selected mop's auto optimal designs
    if optimal_xf_dict is not None:
        # auto_optimal_designs_fitnesses = optimal_fitness_dict[mop.ad.select_spec] # obsolete

        auto_optimal_designs_xf = optimal_xf_dict[mop.ad.select_spec]
        if auto_optimal_designs_xf[0] != []:
            # and not os.path.exists(cairo_fname_optimal_1)
            mop.part_evaluation_geometry(
                auto_optimal_designs_xf[0], counter='CairoOptimal1')
        if auto_optimal_designs_xf[1] != []:
            # and not os.path.exists(cairo_fname_optimal_2)
            mop.part_evaluation_geometry(
                auto_optimal_designs_xf[1], counter='CairoOptimal2')
        if auto_optimal_designs_xf[2] != []:
            # and not os.path.exists(cairo_fname_optimal_3)
            mop.part_evaluation_geometry(
                auto_optimal_designs_xf[2], counter='CairoOptimal3')
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

        displayPDF_side_by_side(
            [cairo_fname, cairo_fname_optimal_1, cairo_fname_optimal_2, cairo_fname_optimal_3])


def SelectIndividualContent():
    st.subheader('Select Individual')

    # 根据用户在 text_input 的输入来筛选符合条件的最优个体
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
                # print(dir(acmop.AC_Machine_Optiomization_Wrapper.part_initialDesign))
                # quit()
                mop.part_initialDesign()
                ad = mop.ad
                # print(ad.acm_template)
                # quit()
                # 手动选择最优个体
                _best_index, _best_individual_data, = None, None
                for ind, el in enumerate(utility_postprocess.call_selection_criteria(ad, eval(user_input_upper_bounds_4filter))):
                    if el is None:
                        raise Exception(str(el))
                    else:
                        _best_index, _proj_name, _best_individual_data_reversed = el
                        _best_individual_data = _best_individual_data_reversed[::-1]

                    st.write(F'\t{el[0]}, {el[1]}, f3={el[-1][0]:.1f}, f2={el[-1][1]:.4f}, f1={el[-1][2]:.1f}, ' +
                             ', '.join(F'{x:.2f}' for x in el[-1][3:]))

                if ind == 0 and _best_index is not None:
                    if st.checkbox('There is only one individual left, do you want to re-produce it?'):
                        mop.part_evaluation_geometry(
                            _best_individual_data, counter='UserSelectedOptimal')
                        displayPDF_side_by_side(
                            [mop.fea_config_dict['output_dir'] + 'indUserSelectedOptimal.pdf', ])
                    number_of_one_optimal_design_selected += 1

                    # print(dir(ad.swarm_data_container))
                    list_of_table_column, fig_donut = utility_postprocess.performance_table_plus_donut_chart(
                        ad, folder, _best_index, _best_individual_data, output_dir=path2project)
                    dict_of_list_of_table_column[folder] = list_of_table_column

                    # for writing paper (table results)

                    # 打印成latex文档直接可以用的表格形式
                    # print('\n\n[Performance table] ready to be copied:')

                    if True:
                        # 自动顺序
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
                        # 手动顺序并添加第一列性能符号
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
                        print('\t', [
                              'TRV', 'FRW', '$T_\\mathrm{rip}$', '$E_m$', '$E_a$', '$\\eta$', 'TRV', 'PF'])
                        for ind, specification in enumerate(ordered_specification_list):
                            list_of_table_column = dict_of_list_of_table_column[specification]
                            print('\t', specification + ' & ' + ' & '.join(
                                [f'{float(el):.1f}' for el in list_of_table_column]))

                        # 添加数据和分隔符（性能纵列）
                        index = 0
                        for specification in ordered_specification_list:
                            list_of_table_column = dict_of_list_of_table_column[specification]

                            print('\t', index, specification)
                            index += 1
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

                    df_performancce = pd.DataFrame(data=dict_of_list_of_table_column, index=[
                                                   'TRV', 'FRW', '$T_\\mathrm{rip}$', '$E_m$', '$E_a$', '$\\eta$', 'Cost', 'disp.PF']).T
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


def SensitivityAnalysisContent():
    st.subheader('Sensitivity Analysis')

    # 之前选择的电机规格中 选择一个较优的设计 然后进行灵敏度分析（还需要一个选择）

    # 选择分析变量
    selected_sensitivity_variables = st.multiselect(
        label='Select Sensitivity Variables',
        options=[]  # TODO: 需要读swarm_dict 找到free的变量 拉出来表格选择 目前可以先不管
        # 暂时不用了 我们直接读GP里面的变量
    )

    # 执行灵敏度分析并展示图像
    if st.button("Run Sensitivity Analysis"):
        mop = swarm_dict[folder]
        ad = mop.acm_template.d['GP']
        print(ad)


def main():

    # 标签页
    OptimizationSetupTab, PopulationTab, SelectIndividualTab, SensitivityAnalysisTab = st.tabs(
        ["Optimization Setup", "Population", "Select Individual", "Sensitivity Analysis"])

    with OptimizationSetupTab:  # Optimization Setup
        OptimizationSetupContent(user_selected_folder)

    with PopulationTab:  # Population
        PopulationContent()

    with SelectIndividualTab:
        SelectIndividualContent()

    with SensitivityAnalysisTab:
        SensitivityAnalysisContent()


if __name__ == '__main__':

    # """ DO NOT MODIFY BEGINS """
    # """ DO NOT MODIFY BEGINS """
    # """ DO NOT MODIFY BEGINS """
    print(f"======================{datetime.date.today()}======================")
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

    # 标题
    with st.sidebar:
        st.title(f'ACMOP Visualization {datetime.date.today()}')

        # 选择项目路径
        path2acmop = os.path.abspath(os.path.dirname(__file__) + '\\..')
        if not os.path.exists(path2acmop + '/_default'):
            os.mkdir(path2acmop + '/_default')
        value = st.session_state['1.path2project'] if '1.path2project' in st.session_state.keys(
        ) else path2acmop + '/_default'
        path2project = st.text_input(
            label='[User] Input path2project:', value=value, on_change=None, key='1.path2project')
        if path2project[-1] != '/' or path2project[-1] != '\\':
            path2project += '/'

        # 选择电机规格
        _, list_specifications, _ = next(os.walk(path2project))
        selected_specifications = st.multiselect(
            label="[User] Select folder(s):", options=list_specifications, default=None, key='2.selected_specifications')

    # 读入硬盘里的电机优化结果的数据
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

            print('Blocking the print statements from AC_Machine_Optiomization_Wrapper...')
            utility.blockPrint()
            swarm_dict[folder] = mop = acmop.AC_Machine_Optiomization_Wrapper(select_fea_config_dict, select_spec, project_loc=path2project)
            utility.enablePrint()
            # 侧边栏 Sidebar
        # Show user selected mop's inputs
        with st.sidebar:
            st.title('User Configurations')

            # 按照所选的电机规格，显示用户输入信息
            user_selected_folder = st.sidebar.selectbox(
                'Select a folder to show its user configurations', selected_specifications, key='3.user_selected_folder')
            mop = swarm_dict[user_selected_folder]

            print(f'{user_selected_folder=}')
            show_user_configurations = st.checkbox('Show User Configurations')
            if show_user_configurations:
                st.sidebar.header('Specifications')
                st.sidebar.table(pd.DataFrame(data=list(mop.spec_input_dict.values()), index=list(
                    mop.spec_input_dict.keys()), dtype="string", columns=['Value',]))
                st.sidebar.header('Simulation Settings')
                st.sidebar.table(pd.DataFrame(data=list(mop.fea_config_dict.values()), index=list(
                    mop.fea_config_dict.keys()), dtype="string", columns=['Value',]))
        main()
    # """ DO NOT MODIFY ENDS """
    # """ DO NOT MODIFY ENDS """
    # """ DO NOT MODIFY ENDS """
