from pylab import np, plt
import pandas as pd
import os, json, builtins, datetime
import streamlit as st
import utility_postprocess, acmop, utility

# Init session_state
"st.session_state object:", st.session_state

## 标题
st.title('ACMOP Visualization')
st.write(f'''### {datetime.date.today()}''')

## 选择电机规格
path2acmop   = os.path.abspath(os.path.dirname(__file__) + '/../')
path2project = path2acmop + '/_default/'
path2project = st.text_input(label='[User] Input path2project:', value=path2project, on_change=None, key='1.path2project')
if path2project[-1]!='/' and path2project[-1]!='\\': path2project += '/'
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

    ## 侧边栏
    st.sidebar.header('User Inputs')
    user_select_input_dict = st.sidebar.selectbox('Select a folder to show its inputs', selected_specifications, key='3.user_select_input_dict')
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
            width: 500px;
        }
        [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
            width: 500px;
            margin-left: -500px;
        }
        </style>
        """,
        unsafe_allow_html=True,
    ) # make sidebar wider!

    mop = swarm_dict[user_select_input_dict]
    st.sidebar.table(pd.DataFrame(data=list(mop.spec_input_dict.values()), index=list(mop.spec_input_dict.keys()), dtype="string", columns=['Value',]))
    st.sidebar.table(pd.DataFrame(data=list(mop.fea_config_dict.values()), index=list(mop.fea_config_dict.keys()), dtype="string", columns=['Value',]))

    ## 集群基本信息
    ## 是否显示自动最优个体表格和帕累托前沿？
    if st.checkbox("Show Swarm Table and Pareto Front?"):
        df, fig = utility_postprocess.inspect_swarm_and_show_table_plus_Pareto_front(swarm_dict, output_dir=path2project)
        st.table(df)
        st.pyplot(fig)

# ## 根据用户在 text_input 的输入来筛选符合条件的最优个体
# if True:
#     df = SA.select_optimal_designs_manually(st, st.session_state, None, ad_list, selected_specifications)
#     st.table(df)

## 结束 (Streamlit widgets automatically run the script from top to bottom. Since this button is not connected to any other logic, it just causes a plain rerun.)
st.button("Re-run")
