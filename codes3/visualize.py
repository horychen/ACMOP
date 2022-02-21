_=""" Regular Python Packages """
from pylab import mpl, np, plt
import pandas as pd
import os, json, builtins

_=""" Streamlit Packages """
import streamlit as st
# from bokeh.plotting import figure
# from bokeh.palettes import Dark2_5 as palette
# import itertools; colors = itertools.cycle(palette) # bokeh colors has a list of colors which can be used in plots 

_=""" My Visualization Utilities """
import utility_postprocess

if __name__ == '__main__':

    # Init session_state
    "st.session_state object:", st.session_state
    if 'ad_list' not in st.session_state:
        st.session_state.ad_list = None
        st.session_state.folder_of_collection = None

    ## 侧边栏测试
    st.sidebar.header('HEADER')
    st.sidebar.selectbox('Test', ('1', '2', '3'))

    ## 标题
    builtins.order = [4,3,5,6,2,1,7] # [1,2,3,4,5,6] # 只影响帕累托图里的marker顺序
    st.title('Swarm Analyzer and Its Visualization')
    st.write(r'''### 2021-03-16''')

    ## 选择优化集合
    path2acmop = os.path.abspath(os.path.dirname(__file__) + '/../')
    SA = utility_postprocess.SwarmAnalyzer(path2acmop)
    last_folder_of_collection = st.session_state.folder_of_collection
    st.session_state.folder_of_collection = st.selectbox("Choose a collection:", SA.folders_of_collections, 0)
    

    if len(SA.folders_of_collections) > 0:
        ## 多项选择电机规格
        list_specifications, dict_path2SwarmDataOfTheSpecification, dict_settingsOfTheSpecification = SA.get_swarm_group(folder_of_collection=st.session_state.folder_of_collection)
        # print(list_specifications)
        selected_specifications = st.multiselect(
            "Choose a specification", list_specifications, 
            [el for el in list_specifications if 'Separate Winding' not in el]
        )

# if __name__ == '!__main__':

        ## 按照所选的电机规格，显示信息
        if not selected_specifications:
            st.error("Please select at least one trace.")

        else:
            ## 读入 ad_list 并保存到 session_state
            # 优点是后面选最优设计的时候快，但是缺点是一旦ad_list或ad代码有变动，由于不是用的st.cache，不会自动刷新，需要重启Streamlit，有时候会忘记，debug不方便。
            if st.session_state.ad_list is None or st.session_state.folder_of_collection != last_folder_of_collection:
                st.session_state.ad_list = ad_list = SA.get_ad_list(selected_specifications)
                # st.session_state.folder_of_collection = last_folder_of_collection
            else:
                ad_list = st.session_state.ad_list

            for ad in ad_list:
                ad.read_swarm_data_json()

            ## 是否显示自动最优个体表格和帕累托前沿？
            if st.checkbox("Show Table and Pareto Front."):
                # st.checkbox("Great", value = True)
                df, fig = SA.inspect_swarm_and_show_table_plus_Pareto_front(ad_list)
                st.table(df) #  df.style.format("{:.2%}")
                st.pyplot(fig)

            ## 根据用户在 text_input 的输入来筛选符合条件的最优个体
            if True:
                df = SA.select_optimal_designs_manually(st, st.session_state, None, ad_list, selected_specifications)
                st.table(df)

            ## 以Json格式保存用户输入到文件
                # st.write("Recorded values: ", st.session_state)
                # with open(f'{os.path.dirname(__file__)}/global_recovered_values_new.json', 'w') as f:
                #     json.dump(st.session_state, f, ensure_ascii=False, indent=4)








        ## 结束
        # Streamlit widgets automatically run the script from top to bottom. Since this button is not connected to any other logic, it just causes a plain rerun.
        st.button("Re-run")
