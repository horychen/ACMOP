import SessionState
import utility_streamlit
session_state = SessionState.get(
                    default_value_global_recovered_values=None,
                    global_recovered_values=None,
                    ad_list=None,
                    folder_of_collection=None,
                )
if session_state.default_value_global_recovered_values is None:
    print('First Session!')
    session_state.global_recovered_values = global_recovered_values = utility_streamlit.load_global_recovered_values()
    from copy import deepcopy
    session_state.default_value_global_recovered_values = deepcopy(global_recovered_values)
    # 如果不把默认值保存下来而且永远不变，你就会发现在网页端修改 st.text_input 的输入是隔一次有效的，因为无效的那一次更新了 st.text_input 的 default value。
    # 而 st.text_input 检测到 default value 变化的话，就会显示更新为 default value 的值，而不是用户这一次输入的值，所以在用户看来，就是第一次输入的无效的，要输入第二次才有效。
else:
    global_recovered_values = session_state.global_recovered_values

_=""" Regular Python Packages """
from pylab import mpl, np, plt
import pandas as pd
import os, json
import builtins

_=""" Streamlit Packages """
import streamlit as st
from bokeh.plotting import figure
from bokeh.palettes import Dark2_5 as palette
import itertools; colors = itertools.cycle(palette) # bokeh colors has a list of colors which can be used in plots 

_=""" Should not load ACMOP/codes3 Modules """
import acmop
import utility_postprocess

if __name__ == '__main__':
    mop = acmop.AC_Machine_Optiomization_Wrapper(
            select_fea_config_dict = '#0211 JMAG PMSM Q12p4ps5 Sub-hamonics',
            select_spec            = 'PMSM Q12p4y1 A',
            project_loc            = r'D:/DrH/acmop/_default/',
            bool_show_jmag         = True
        )
    mop.part_reportWithStreamlit()

    ## 标题
    builtins.title = 'Swarm Analyzer and Visualization'
    builtins.order  = [4,3,5,6,2,1,7] # [1,2,3,4,5,6] # 只影响帕累托图里的marker顺序
    st.title(builtins.title)
    st.write(r'''### 2021-03-16''')

    ## 选择优化集合

    SA = utility_postprocess.SwarmAnalyzer(os.path.dirname(__file__))
    SA.get_folders_of_collections()

    folder_of_collection = st.selectbox(
        "Choose a collection:", SA.folders_of_collections, 0
    )

    ## 多项选择电机规格
    list_specifications, dict_path2SwarmDataOfTheSpecification, dict_settingsOfTheSpecification \
        = SA.get_swarm_group(folder_of_collection=folder_of_collection)
    # print(list_specifications)
    selected_specifications = st.multiselect(
        "Choose a specification", list_specifications, 
        [el for el in list_specifications if 'Separate Winding' not in el]
    )

    ## 侧边栏测试
    st.sidebar.header('HEADER')
    st.sidebar.selectbox('Test', ('1', '2', '3'))


# if __name__ == '!__main__':


    ## 按照所选的电机规格，显示信息
    if not selected_specifications:
        st.error("Please select at least one trace.")
    else:

        ## 读入 ad_list 并保存到 session_state
        if session_state.ad_list is None or session_state.folder_of_collection != folder_of_collection:
            session_state.ad_list = ad_list = SA.get_ad_list(selected_specifications)
            session_state.folder_of_collection = folder_of_collection
        else:
            ad_list = session_state.ad_list

        ## 是否显示自动最优个体表格和帕累托前沿？
        if st.checkbox("Show Table and Pareto Front."):
            # st.checkbox("Great", value = True)
            df, fig = SA.inspect_swarm_and_show_table_plus_Pareto_front(ad_list)
            st.table(df) #  df.style.format("{:.2%}")
            st.pyplot(fig)

        ## 根据用户在 text_input 的输入来筛选符合条件的最优个体
        if True:
            df = SA.select_optimal_designs_manually(st, session_state, global_recovered_values, 
                                                    ad_list, selected_specifications)
            st.table(df)

        ## 以Json格式保存用户输入到文件
        st.write("Recorded values: ", global_recovered_values)
        with open(f'{os.path.dirname(__file__)}/global_recovered_values.json', 'w') as f:
            json.dump(global_recovered_values, f, ensure_ascii=False, indent=4)








    ## 结束
    # Streamlit widgets automatically run the script from top to bottom. Since this button is not connected to any other logic, it just causes a plain rerun.
    st.button("Re-run")
