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
from itertools import product

import copy

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

def BasicInformationContent():

    wily_fname = 'wily_p%dps%dQ%dy%d' % (
        mop.spec_input_dict['p'], mop.spec_input_dict['ps'], mop.spec_input_dict['Qs'], mop.spec_input_dict['coil_pitch_y'])
    st.write('## 2.1. Winding Information of', wily_fname)
    try:
        f'{path2acmop}/_wily/{wily_fname}.pdf'
        # displayPDF(f'{path2acmop}/_wily/{wily_fname}.pdf')
        displayPDF_side_by_side([
                f'{path2acmop}/_wily/{wily_fname}_T1.pdf',
                f'{path2acmop}/_wily/{wily_fname}_T2.pdf',
                f'{path2acmop}/_wily/{wily_fname}_T3abc.pdf',
            ], width=600, height=400
        )
        displayPDF(f'{path2acmop}/_wily/{wily_fname}_T4.pdf', width=600, height=230)
        displayPDF_side_by_side([
                f'{path2acmop}/_wily/{wily_fname}_T4a.pdf',
                f'{path2acmop}/_wily/{wily_fname}_T4b.pdf',
                f'{path2acmop}/_wily/{wily_fname}_T4c.pdf',
            ], width=600, height=300
        )
    except FileNotFoundError:
        st.write('The winding derivation file is absent', f'{path2acmop}/_wily/{wily_fname}.pdf')


def InitialDesignContent(user_selected_folder):
    mop = swarm_dict[user_selected_folder]
    # mop

    col1, col2 = st.columns(2)

    # 左栏内容
    with col1:
        st.header(f'Optimization Setup of {user_selected_folder}')

        st.write('#### Decision Variables:')
        gp_items = mop.ad.acm_template.d['GP'].items()
        decision_vars = [(key, val) for key, val in gp_items if 'free' in val.type]
        
        if decision_vars:
            # 初始化 session_state 存储用户修改的值
            if 'decision_var_values' not in st.session_state:
                st.session_state.decision_var_values = {
                    var_name: var_param.value for var_name, var_param in decision_vars
                }
            
            # 初始化原始值（用于检测修改）
            if 'decision_var_original_values' not in st.session_state:
                st.session_state.decision_var_original_values = {
                    var_name: var_param.value for var_name, var_param in decision_vars
                }
            
            # 使用列布局创建可编辑表格
            # 显示表格标题
            cols_header = st.columns([2, 2, 2, 3])
            with cols_header[0]:
                st.write("**Variable**")
            with cols_header[1]:
                st.write("**Value**")
            with cols_header[2]:
                st.write("**Original**")
            with cols_header[3]:
                st.write("**Bounds**")
            
            # 为每行创建可编辑的输入
            edited_values = {}
            for idx, (var_name, var_param) in enumerate(decision_vars):
                original_value = st.session_state.decision_var_original_values.get(var_name, var_param.value)
                current_value = st.session_state.decision_var_values.get(var_name, var_param.value)
                
                cols = st.columns([2, 2, 2, 3])
                
                with cols[1]:
                    new_value = st.number_input(
                        f"Value for {var_name}",
                        value=float(current_value),
                        step=0.001,
                        format="%.6f",
                        key=f"decision_var_{var_name}",
                        label_visibility="collapsed"
                    )
                    edited_values[var_name] = new_value
                    st.session_state.decision_var_values[var_name] = new_value
                
                # 在更新值后重新检查是否被修改
                is_modified = abs(new_value - original_value) > 1e-10
                
                with cols[0]:
                    if is_modified:
                        st.markdown(f'<div style="background-color: #ffeb3b; padding: 5px; border-radius: 3px;">{var_name}</div>', unsafe_allow_html=True)
                    else:
                        st.write(var_name)
                
                with cols[2]:
                    st.write(f"{original_value:.6f}")
                
                with cols[3]:
                    if hasattr(var_param, 'bounds') and var_param.bounds:
                        st.write(f"[{var_param.bounds[0]:.6f}, {var_param.bounds[1]:.6f}]")
                    else:
                        st.write("无 bounds")
        else:
            st.write("No decision variables.")
        
        
        st.write('#### Derived Variables:')
        gp_items = mop.ad.acm_template.d['GP'].items()
        derived_vars = [(key, val) for key, val in gp_items if 'derived' in val.type]
        
        if derived_vars:
            # 初始化 session_state 存储 Derived Variables 的原始值和修改后的值
            if 'derived_var_original_values' not in st.session_state:
                st.session_state.derived_var_original_values = {
                    var_name: var_param.value for var_name, var_param in derived_vars
                }
            if 'derived_var_values' not in st.session_state:
                st.session_state.derived_var_values = {
                    var_name: var_param.value for var_name, var_param in derived_vars
                }
            
            # 显示表格标题
            cols_header = st.columns([2, 2, 2])
            with cols_header[0]:
                st.write("**Variable**")
            with cols_header[1]:
                st.write("**Value**")
            with cols_header[2]:
                st.write("**Original**")
            
            # 为每行创建可编辑的输入
            for var_name, var_param in derived_vars:
                original_value = st.session_state.derived_var_original_values.get(var_name, var_param.value)
                current_value = st.session_state.derived_var_values.get(var_name, var_param.value)
                
                cols = st.columns([2, 2, 2])
                
                with cols[1]:
                    new_value = st.number_input(
                        f"Value for {var_name}",
                        value=float(current_value),
                        step=0.001,
                        format="%.6f",
                        key=f"derived_var_{var_name}",
                        label_visibility="collapsed"
                    )
                    st.session_state.derived_var_values[var_name] = new_value
                
                # 在更新值后重新检查是否被修改
                is_modified = abs(new_value - original_value) > 1e-10
                
                with cols[0]:
                    if is_modified:
                        st.markdown(f'<div style="background-color: #ffeb3b; padding: 5px; border-radius: 3px;">{var_name}</div>', unsafe_allow_html=True)
                    else:
                        st.write(var_name)
                
                with cols[2]:
                    st.write(f"{original_value:.6f}")
        else:
            st.write("No derived variables.")

        st.write('#### Fixed Variables:')
        gp_items = mop.ad.acm_template.d['GP'].items()
        fixed_vars = [(key, val) for key, val in gp_items if 'fixed' in val.type]
        
        if fixed_vars:
            # 初始化 session_state 存储 Fixed Variables 的原始值和修改后的值
            if 'fixed_var_original_values' not in st.session_state:
                st.session_state.fixed_var_original_values = {
                    var_name: var_param.value for var_name, var_param in fixed_vars
                }
            if 'fixed_var_values' not in st.session_state:
                st.session_state.fixed_var_values = {
                    var_name: var_param.value for var_name, var_param in fixed_vars
                }
            
            # 显示表格标题
            cols_header = st.columns([2, 2, 2])
            with cols_header[0]:
                st.write("**Variable**")
            with cols_header[1]:
                st.write("**Value**")
            with cols_header[2]:
                st.write("**Original**")
            
            # 为每行创建可编辑的输入
            for var_name, var_param in fixed_vars:
                original_value = st.session_state.fixed_var_original_values.get(var_name, var_param.value)
                current_value = st.session_state.fixed_var_values.get(var_name, var_param.value)
                
                cols = st.columns([2, 2, 2])
                
                with cols[1]:
                    new_value = st.number_input(
                        f"Value for {var_name}",
                        value=float(current_value),
                        step=0.001,
                        format="%.6f",
                        key=f"fixed_var_{var_name}",
                        label_visibility="collapsed"
                    )
                    st.session_state.fixed_var_values[var_name] = new_value
                
                # 在更新值后重新检查是否被修改
                is_modified = abs(new_value - original_value) > 1e-10
                
                with cols[0]:
                    if is_modified:
                        st.markdown(f'<div style="background-color: #ffeb3b; padding: 5px; border-radius: 3px;">{var_name}</div>', unsafe_allow_html=True)
                    else:
                        st.write(var_name)
                
                with cols[2]:
                    st.write(f"{original_value:.6f}")
        else:
            st.write("No fixed variables.")

        st.write('#### Objectives (Non-dominated Sorting)')
        st.write(f'''
                    {mop.fea_config_dict['moo.fitness_OA']=}
                    {mop.fea_config_dict['moo.fitness_OB']=}
                    {mop.fea_config_dict['moo.fitness_OC']=}
        ''')

        # mop.ad.acm_template.d['GP']

    # 右栏内容
    with col2:


        # 获取所有 Decision Variables
        gp_items = mop.ad.acm_template.d['GP'].items()
        decision_vars = [(key, val) for key, val in gp_items if 'free' in val.type]
        
        if decision_vars and 'decision_var_values' in st.session_state:
            # 根据用户修改的值构建 x_denorm
            # 需要按照 x_denorm_dict 的顺序来构建
            edited_values = st.session_state.decision_var_values
            x_denorm_dict = copy.deepcopy(mop.ad.acm_template.x_denorm_dict)
            
            # 更新 x_denorm_dict 中用户修改的值
            for var_name in x_denorm_dict.keys():
                if var_name in edited_values:
                    x_denorm_dict[var_name] = edited_values[var_name]
            
            # 更新 GP 中的 Fixed Variables（如果用户修改了）
            if 'fixed_var_values' in st.session_state:
                for var_name, new_value in st.session_state.fixed_var_values.items():
                    if var_name in mop.ad.acm_template.d['GP']:
                        mop.ad.acm_template.d['GP'][var_name].value = new_value
            
            # 更新 GP 中的 Derived Variables（如果用户修改了）
            if 'derived_var_values' in st.session_state:
                for var_name, new_value in st.session_state.derived_var_values.items():
                    if var_name in mop.ad.acm_template.d['GP']:
                        mop.ad.acm_template.d['GP'][var_name].value = new_value
            
            # 按照 x_denorm_dict 的顺序构建 x_denorm 列表
            specify_x_denorm = list(x_denorm_dict.values())
            
            st.write("#### 当前参数值:")
            # 显示修改的变量（高亮）
            modified_vars = []
            if 'decision_var_original_values' in st.session_state:
                for var_name in x_denorm_dict.keys():
                    if var_name in st.session_state.decision_var_original_values:
                        orig_val = st.session_state.decision_var_original_values[var_name]
                        curr_val = x_denorm_dict[var_name]
                        if abs(curr_val - orig_val) > 1e-10:
                            modified_vars.append(var_name)
            
            param_display_parts = []
            for name, val in x_denorm_dict.items():
                if name in modified_vars:
                    param_display_parts.append(f"**{name}={val:.6f}**")
                else:
                    param_display_parts.append(f"{name}={val:.6f}")
            st.write(" | ".join(param_display_parts))
            
            # 显示修改的 Fixed Variables
            if 'fixed_var_values' in st.session_state and 'fixed_var_original_values' in st.session_state:
                modified_fixed = []
                for var_name, curr_val in st.session_state.fixed_var_values.items():
                    orig_val = st.session_state.fixed_var_original_values.get(var_name, curr_val)
                    if abs(curr_val - orig_val) > 1e-10:
                        modified_fixed.append((var_name, orig_val, curr_val))
                
                if modified_fixed:
                    st.write("#### 修改的 Fixed Variables:")
                    fixed_display = " | ".join([f"**{name}**: {orig:.6f} → {curr:.6f}" for name, orig, curr in modified_fixed])
                    st.write(fixed_display)
            
            # 显示修改的 Derived Variables
            if 'derived_var_values' in st.session_state and 'derived_var_original_values' in st.session_state:
                modified_derived = []
                for var_name, curr_val in st.session_state.derived_var_values.items():
                    orig_val = st.session_state.derived_var_original_values.get(var_name, curr_val)
                    if abs(curr_val - orig_val) > 1e-10:
                        modified_derived.append((var_name, orig_val, curr_val))
                
                if modified_derived:
                    st.write("#### 修改的 Derived Variables:")
                    derived_display = " | ".join([f"**{name}**: {orig:.6f} → {curr:.6f}" for name, orig, curr in modified_derived])
                    st.write(derived_display)
            
            # 自动生成图形（不需要按钮）
            with st.spinner("正在生成图形..."):
                try:
                    # 生成 PDF
                    saved_filename = mop.part_evaluation_geometry(specify_x_denorm=specify_x_denorm, counter='CurrentDesign')
                    
                    # 尝试多个可能的路径
                    pdf_paths_to_try = []
                    if saved_filename:
                        pdf_paths_to_try.append(saved_filename)
                    pdf_paths_to_try.append(mop.fea_config_dict['output_dir'] + 'indCurrentDesign.pdf')
                    if saved_filename and isinstance(saved_filename, str):
                        # 如果返回的是相对路径，尝试构建完整路径
                        if not os.path.isabs(saved_filename):
                            pdf_paths_to_try.append(os.path.join(mop.fea_config_dict['output_dir'], saved_filename))
                    
                    pdf_found = None
                    for pdf_path in pdf_paths_to_try:
                        if pdf_path and os.path.exists(pdf_path):
                            pdf_found = pdf_path
                            break
                    
                    if pdf_found:
                        st.write("#### 生成的图形:")
                        displayPDF_side_by_side([pdf_found], width='100%', height=600)
                    else:
                        st.warning(f"PDF 文件未找到。尝试的路径: {pdf_paths_to_try}")
                except Exception as e:
                    st.error(f"生成图形时出错: {str(e)}")
                    import traceback
                    st.code(traceback.format_exc())
        else:
            st.info("请先在左栏修改 Decision Variables 的值")


def SearchSpaceContent(user_selected_folder):
    mop = swarm_dict[user_selected_folder]
    
    st.header("参数扫描功能")
    
    # 获取所有 Decision Variables
    gp_items = mop.ad.acm_template.d['GP'].items()
    decision_vars = [(key, val) for key, val in gp_items if 'free' in val.type]
    
    if len(decision_vars) == 0:
        st.warning("没有找到 Decision Variables")
    else:
        # 为每个变量计算大中小三个值
        var_values = []
        var_names = []
        for var_name, var_param in decision_vars:
            if hasattr(var_param, 'bounds') and var_param.bounds:
                lower = var_param.bounds[0]
                upper = var_param.bounds[1]
                middle = (lower + upper) / 2
                var_values.append([lower, middle, upper])
                var_names.append(var_name)
            else:
                st.error(f"变量 {var_name} 没有 bounds")
                return
        
        # 生成所有组合（3^n 种，n 是 Decision Variables 的数量）
        combinations = list(product(*var_values))
        num_vars = len(var_names)
        total_combinations = 3 ** num_vars
        
        st.write(f"#### 扫描参数空间: {len(combinations)} 种组合 (3^{num_vars} = {total_combinations})")
        st.write(f"变量数量: {num_vars}")
        st.write(f"变量: {', '.join(var_names)}")
        
        # 添加按钮，只有点击后才开始扫描和画图
        if st.button("开始扫描和画图", type="primary"):
            # 进度条
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # 生成所有 PDF
            pdf_files = []
            pdf_params = []  # 保存每个PDF对应的参数信息
            
            for idx, combo in enumerate(combinations):
                status_text.text(f"处理中: {idx+1}/{len(combinations)}")
                progress_bar.progress((idx + 1) / len(combinations))
                
                # 构建 x_denorm：使用所有变量的组合值
                # combo 的长度应该等于 decision_vars 的数量
                x_denorm = list(combo)
                
                # 生成唯一的 counter 名称
                counter_name = f'Scan_{idx:03d}'
                
                # 运行 part_evaluation_geometry
                try:
                    mop.part_evaluation_geometry(specify_x_denorm=x_denorm, counter=counter_name)
                    pdf_path = mop.fea_config_dict['output_dir'] + f'ind{counter_name}.pdf'
                    if os.path.exists(pdf_path):
                        pdf_files.append(pdf_path)
                        # 保存参数信息：变量名和对应的值
                        param_info = {var_names[j]: combo[j] for j in range(num_vars)}
                        pdf_params.append(param_info)
                except Exception as e:
                    st.warning(f"组合 {idx+1} 处理失败: {str(e)}")
            
            status_text.text(f"完成！共生成 {len(pdf_files)} 个 PDF 文件")
            
            # 显示所有 PDF，每行显示三个
            if pdf_files:
                st.write("#### 生成的 PDF 文件:")
                # 每行显示三个PDF
                for i in range(0, len(pdf_files), 3):
                    # 创建三列
                    cols = st.columns(3)
                    
                    # 在每列中显示 PDF 和参数信息
                    for col_idx in range(3):
                        idx = i + col_idx
                        if idx < len(pdf_files):
                            pdf_path = pdf_files[idx]
                            param_info = pdf_params[idx]
                            
                            with cols[col_idx]:
                                # 显示参数信息
                                param_text = " | ".join([f"{name}={val:.3f}" for name, val in param_info.items()])
                                st.write(f"**[{idx+1}]** {param_text}")
                                
                                # 单独显示每个PDF
                                displayPDF_side_by_side([pdf_path], width=300, height=300)


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
    # wily_fname = 'wily_p%dps%dQ%dy%d' % (
    #     mop.spec_input_dict['p'], mop.spec_input_dict['ps'], mop.spec_input_dict['Qs'], mop.spec_input_dict['coil_pitch_y'])
    # st.write('## 2.1. Winding Information of', wily_fname)
    # try:
    #     displayPDF(f'{path2acmop}_wily/{wily_fname}.pdf')
    # except FileNotFoundError:
    #     st.write('The winding derivation file is absent')

    st.write('## 2.2. Cross Section (Initial and Optimal Designs)')
    cairo_fname = mop.fea_config_dict['output_dir'] + 'indCairo.pdf'
    cairo_fname_optimal_1 = mop.fea_config_dict['output_dir'] + 'indCairoOptimal1.pdf'
    cairo_fname_optimal_2 = mop.fea_config_dict['output_dir'] + 'indCairoOptimal2.pdf'
    cairo_fname_optimal_3 = mop.fea_config_dict['output_dir'] + 'indCairoOptimal3.pdf'
    if not os.path.exists(cairo_fname):
        mop.part_evaluation_geometry()

    # Show user selected mop's auto optimal designs
    if optimal_xf_dict is not None:
        # auto_optimal_designs_fitnesses = optimal_fitness_dict[mop.ad.select_spec] # obsolete

        auto_optimal_designs_xf = optimal_xf_dict[mop.ad.select_spec]
        if auto_optimal_designs_xf[0] != []:
            # and not os.path.exists(cairo_fname_optimal_1)
            mop.part_evaluation_geometry(auto_optimal_designs_xf[0], counter='CairoOptimal1')
        if auto_optimal_designs_xf[1] != []:
            # and not os.path.exists(cairo_fname_optimal_2)
            mop.part_evaluation_geometry(auto_optimal_designs_xf[1], counter='CairoOptimal2')
        if auto_optimal_designs_xf[2] != []:
            # and not os.path.exists(cairo_fname_optimal_3)
            mop.part_evaluation_geometry(auto_optimal_designs_xf[2], counter='CairoOptimal3')

        pop_col1, pop_col2 = st.columns(2)

        # 左栏内容
        with pop_col1:

            st.write('### 2.2.1 Initial Design:')
            # displayPDF(cairo_fname, width=500, height=500)
            st.write(mop.acm_variant.x_denorm)

            st.write('### 2.2.2 Minimum OA Design xf')
            auto_optimal_designs_xf[0]
            # displayPDF(cairo_fname_optimal_1, width=500, height=500)

            st.write('### 2.2.3 Minimum OB Design xf')
            auto_optimal_designs_xf[1]
            # displayPDF(cairo_fname_optimal_2, width=500, height=500)

            st.write('### 2.2.4 Minimum OC Design xf')
            auto_optimal_designs_xf[2]
            # displayPDF(cairo_fname_optimal_3, width=500, height=500)

        with pop_col2:
            displayPDF_side_by_side(
                [cairo_fname, cairo_fname_optimal_1, cairo_fname_optimal_2, cairo_fname_optimal_3])


def SelectIndividualContent():
    st.subheader('Select Individual')

    # 根据用户在 text_input 的输入来筛选符合条件的最优个体
    if True:
        print(f'\n\n\n--------------------------{selected_specifications=}')

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
        mop
        mop.acm_template.d['GP']


def main():

    # 标签页
    BasicInformationTab, InitialDesignTab, SearchSpaceTab, PopulationTab, SelectIndividualTab, SensitivityAnalysisTab = st.tabs(
        ["BasicInformation", "Initial Design", "Search Space", "Population", "Select Individual", "Sensitivity Analysis"])

    with BasicInformationTab:
        BasicInformationContent()

    with InitialDesignTab:  # Initial Design
        InitialDesignContent(user_selected_folder)

    with SearchSpaceTab:  # Search Space
        SearchSpaceContent(user_selected_folder)

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
    st.set_page_config(layout="wide")
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
            show_user_configurations = st.checkbox('Show User Configurations', value=True)
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
