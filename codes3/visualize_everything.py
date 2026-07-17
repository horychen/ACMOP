from pylab import np, plt
import pandas as pd
import os, json, builtins, datetime
import streamlit as st
import pygmo as pg
import utility_postprocess, acmop, utility, base64, acm_designer, population, postprocess_data
# import acm_designer
from io import BytesIO

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


def _load_metrics_for_selected_folders(path2project, swarm_dict):
    selected_specs = [mop.select_spec for mop in swarm_dict.values()]
    df = postprocess_data.load_acmop_metrics(path2project)
    if selected_specs and not df.empty:
        df = df[df["select_spec"].isin(selected_specs)].copy()
    return df


def _metric_columns_available(df):
    preferred_columns = [
        "select_spec",
        "archive_record_no",
        "design_key",
        "project_name",
        "Cost",
        "rated_efficiency_pct",
        "TRV_kNm_per_m3",
        "FRW",
        "normalized_torque_ripple_pct",
        "normalized_force_error_magnitude_pct",
        "force_error_angle",
        "power_factor",
        "rated_total_loss",
        "rated_stack_length_mm",
        "torque_average",
        "ss_avg_force_magnitude",
        "Cost_Fe",
        "Cost_Cu",
        "Cost_PM",
        "f1",
        "f2",
        "f3",
        "valid_metrics",
    ]
    return [column for column in preferred_columns if column in df.columns]


def _show_selected_design_metrics(row):
    metric_specs = [
        ("Cost", "Cost", "{:.3g}"),
        ("Efficiency", "rated_efficiency_pct", "{:.2f}%"),
        ("TRV", "TRV_kNm_per_m3", "{:.2f} kNm/m3"),
        ("FRW", "FRW", "{:.3f}"),
        ("Torque Ripple", "normalized_torque_ripple_pct", "{:.2f}%"),
        ("Force Error", "normalized_force_error_magnitude_pct", "{:.2f}%"),
        ("Error Angle", "force_error_angle", "{:.2f} deg"),
        ("PF", "power_factor", "{:.3f}"),
    ]
    cols = st.columns(4)
    for index, (label, column, fmt) in enumerate(metric_specs):
        value = row.get(column, np.nan)
        if pd.isna(value):
            display_value = "N/A"
        else:
            display_value = fmt.format(float(value))
        cols[index % 4].metric(label, display_value)


def _plot_metric_scatter(df):
    plot_df = df.dropna(subset=["Cost", "rated_efficiency_pct", "f3"])
    if plot_df.empty:
        st.info("No Cost / Efficiency / f3 data available for scatter plot.")
        return
    fig, ax = plt.subplots(figsize=(8, 5), dpi=140)
    scatter = ax.scatter(
        plot_df["Cost"],
        plot_df["rated_efficiency_pct"],
        c=plot_df["f3"],
        cmap="viridis",
        s=28,
        alpha=0.78,
        edgecolors="none",
    )
    ax.set_xlabel("Cost")
    ax.set_ylabel("Rated efficiency [%]")
    ax.grid(True, alpha=0.25)
    ax.set_title("Cost vs Efficiency, colored by ripple objective")
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("f3 ripple objective")
    st.pyplot(fig)


def _plot_ripple_scatter(df):
    plot_df = df.dropna(subset=["normalized_torque_ripple_pct", "normalized_force_error_magnitude_pct", "force_error_angle"])
    if plot_df.empty:
        st.info("No ripple/error data available for scatter plot.")
        return
    fig, ax = plt.subplots(figsize=(8, 5), dpi=140)
    scatter = ax.scatter(
        plot_df["normalized_torque_ripple_pct"],
        plot_df["normalized_force_error_magnitude_pct"],
        c=plot_df["force_error_angle"],
        cmap="plasma",
        s=28,
        alpha=0.78,
        edgecolors="none",
    )
    ax.set_xlabel("Torque ripple [%]")
    ax.set_ylabel("Force error magnitude [%]")
    ax.grid(True, alpha=0.25)
    ax.set_title("Ripple and suspension force error")
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Force error angle [deg]")
    st.pyplot(fig)


def _plot_loss_breakdown(row):
    loss_columns = {
        "Stator copper": "rated_stator_copper_loss_along_stack",
        "Magnet Joule": "rated_magnet_Joule_loss",
        "Rotor copper": "rated_rotor_copper_loss_along_stack",
        "End turn stator": "stator_copper_loss_in_end_turn",
        "End turn rotor": "rotor_copper_loss_in_end_turn",
        "Iron": "rated_iron_loss",
        "Windage": "rated_windage_loss",
    }
    values = []
    labels = []
    for label, column in loss_columns.items():
        value = row.get(column, np.nan)
        if pd.notna(value) and float(value) > 0:
            labels.append(label)
            values.append(float(value))
    if not values:
        st.info("No positive loss breakdown data available for the selected design.")
        return
    fig, ax = plt.subplots(figsize=(6, 5), dpi=140)
    ax.pie(values, labels=labels, autopct="%1.1f%%", startangle=90)
    ax.axis("equal")
    ax.set_title("Rated loss breakdown")
    st.pyplot(fig)


PARETO_HIGHLIGHT_METRICS = {
    "TRV": {
        "label": "Max TRV",
        "short_label": "TRV max",
        "sources": [("TRV", 1.0), ("f1", -1.0)],
        "direction": "max",
        "scale": 0.001,
        "unit": "kNm/m3",
        "color": "#1f77b4",
        "marker": "o",
    },
    "FRW": {
        "label": "Max FRW",
        "short_label": "FRW max",
        "sources": [("FRW", 1.0)],
        "direction": "max",
        "scale": 1.0,
        "unit": "1",
        "color": "#e07a1f",
        "marker": "^",
    },
    "normalized_torque_ripple": {
        "label": "Min torque ripple",
        "short_label": r"$T_{ripple}$ min",
        "sources": [("normalized_torque_ripple", 1.0)],
        "direction": "min",
        "scale": 100.0,
        "unit": "%",
        "color": "#c43c35",
        "marker": "*",
    },
    "normalized_force_error_magnitude": {
        "label": r"Min $E_m$",
        "short_label": r"$E_m$ min",
        "sources": [("normalized_force_error_magnitude", 1.0)],
        "direction": "min",
        "scale": 100.0,
        "unit": "%",
        "color": "#168a6b",
        "marker": "D",
    },
    "force_error_angle": {
        "label": r"Min $E_a$",
        "short_label": r"$E_a$ min",
        "sources": [("force_error_angle", 1.0)],
        "direction": "min",
        "scale": 1.0,
        "unit": "deg",
        "color": "#2f6fb0",
        "marker": "P",
    },
    "rated_efficiency": {
        "label": r"Max $\eta$",
        "short_label": r"$\eta$ max",
        "sources": [("RatedEfficiency", 1.0), ("rated_efficiency", 1.0), ("f2", -1.0)],
        "direction": "max",
        "scale": 100.0,
        "unit": "%",
        "color": "#7a5195",
        "marker": "h",
    },
    "Cost": {
        "label": "Min Cost",
        "short_label": "Cost min",
        "sources": [("Cost", 1.0), ("f1", 1.0)],
        "direction": "min",
        "scale": 1.0,
        "unit": "USD",
        "color": "#6b6258",
        "marker": "X",
    },
    "power_factor": {
        "label": "Max PF",
        "short_label": "PF max",
        "sources": [("power_factor", 1.0)],
        "direction": "max",
        "scale": 1.0,
        "unit": "1",
        "color": "#d95f8d",
        "marker": "v",
    },
}


def _rank_one_indices(fitnesses):
    if len(fitnesses) == 0:
        return np.asarray([], dtype=int)
    if len(fitnesses) == 1:
        return np.asarray([0], dtype=int)
    fronts, _, _, _ = pg.fast_non_dominated_sorting(fitnesses)
    return np.asarray(fronts[0], dtype=int)


def _get_analyzer_metric_values(analyzer, metric_style):
    for source_key, multiplier in metric_style["sources"]:
        try:
            values = np.asarray(
                analyzer.get_metric_of_the_whole_swarm(source_key),
                dtype=float,
            )
        except (KeyError, TypeError, ValueError):
            continue
        return values * multiplier
    return None


def _find_pareto_metric_optima(mop, z_filter, metric_keys):
    analyzer = mop.ad.analyzer
    swarm_data_xf = analyzer.swarm_data_xf or []
    if not swarm_data_xf:
        return []

    fitnesses = np.asarray([individual[-3:] for individual in swarm_data_xf], dtype=float)
    rank_one = _rank_one_indices(fitnesses)
    candidate_mask = (
        np.isfinite(fitnesses[rank_one]).all(axis=1)
        & (fitnesses[rank_one, 2] < z_filter)
    )
    candidate_indices = rank_one[candidate_mask]
    if len(candidate_indices) == 0:
        return []

    project_names = getattr(analyzer, "swarm_data_project_names", [])
    selections = []
    for metric_key in metric_keys:
        style = PARETO_HIGHLIGHT_METRICS[metric_key]
        metric_values = _get_analyzer_metric_values(analyzer, style)
        if metric_values is None:
            continue
        if len(metric_values) != len(swarm_data_xf):
            continue

        finite_candidates = candidate_indices[np.isfinite(metric_values[candidate_indices])]
        if len(finite_candidates) == 0:
            continue
        direction_multiplier = 1.0 if style["direction"] == "min" else -1.0
        best_index = min(
            finite_candidates.tolist(),
            key=lambda index: (
                direction_multiplier * metric_values[index],
                fitnesses[index, 2],
            ),
        )
        selections.append(
            {
                "metric_key": metric_key,
                "archive_index": int(best_index),
                "project_name": project_names[best_index] if best_index < len(project_names) else "",
                "metric_value": float(metric_values[best_index]),
                "OA": float(fitnesses[best_index, 0]),
                "OB": float(fitnesses[best_index, 1]),
                "OC": float(fitnesses[best_index, 2]),
            }
        )
    return selections


def _highlight_pareto_metric_optima(fig, swarm_dict, z_filter, metric_keys):
    if not metric_keys:
        return pd.DataFrame()

    ax = fig.axes[0]
    rows = []
    legend_labels = set()
    annotation_groups = {}
    for spec_number, (_, mop) in enumerate(swarm_dict.items(), start=1):
        for selection in _find_pareto_metric_optima(mop, z_filter, metric_keys):
            style = PARETO_HIGHLIGHT_METRICS[selection["metric_key"]]
            legend_label = style["short_label"]
            scatter_label = legend_label if legend_label not in legend_labels else None
            legend_labels.add(legend_label)

            x_coord = selection["OA"]
            y_coord = selection["OB"] * 100.0
            ax.scatter(
                [x_coord],
                [y_coord],
                s=190,
                marker=style["marker"],
                c=style["color"],
                edgecolors="black",
                linewidths=1.0,
                zorder=300,
                label=scatter_label,
            )
            annotation_key = (spec_number, round(x_coord, 12), round(y_coord, 12))
            annotation_group = annotation_groups.setdefault(
                annotation_key,
                {"x": x_coord, "y": y_coord, "labels": []},
            )
            annotation_group["labels"].append(style["short_label"])
            rows.append(
                {
                    "Specification": mop.ad.select_spec,
                    "Metric": style["label"],
                    "Individual": selection["project_name"] or f'archive #{selection["archive_index"] + 1}',
                    "Value": selection["metric_value"] * style["scale"],
                    "Unit": style["unit"],
                    "OA": selection["OA"],
                    "OB": selection["OB"],
                    "OC": selection["OC"],
                }
            )

    if rows:
        x_limits = ax.get_xlim()
        y_limits = ax.get_ylim()
        for (spec_number, _, _), annotation_group in annotation_groups.items():
            x_coord = annotation_group["x"]
            y_coord = annotation_group["y"]
            x_is_right = x_coord >= sum(x_limits) / 2.0
            y_is_top = y_coord >= sum(y_limits) / 2.0
            ax.annotate(
                f'{", ".join(annotation_group["labels"])} ({spec_number})',
                (x_coord, y_coord),
                xytext=(-8 if x_is_right else 8, -8 if y_is_top else 8),
                textcoords="offset points",
                fontsize=8.5,
                color="#202124",
                horizontalalignment="right" if x_is_right else "left",
                verticalalignment="top" if y_is_top else "bottom",
                zorder=301,
            )
        legend = ax.legend(loc="best", ncol=2, fontsize=9)
        legend.set_zorder(555)
    return pd.DataFrame(rows)


# from emy-c (maintained by zjl)
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


if __name__ == '__main__':

    # """ DO NOT MODIFY BEGINS """
    # """ DO NOT MODIFY BEGINS """
    # """ DO NOT MODIFY BEGINS """
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
    history = get_user_config(); print(f'\n\n\n{history=}')
    ## 标题
    st.title(f'ACMOP Visualization {datetime.date.today()}')

    ## 选择项目路径
    path2acmop = os.path.abspath(os.path.dirname(__file__) + '\\..')
    if not os.path.exists(path2acmop + '/_default'): os.mkdir(path2acmop + '/_default')
    value = st.session_state['1.path2project'] if '1.path2project' in st.session_state.keys() else path2acmop + '/_default'
    path2project = st.text_input(label='[User] Input path2project:', value=value, on_change=None, key='1.path2project')
    if path2project[-1]!='/' or path2project[-1]!='\\': path2project += '/'

    ## 选择电机规格
    _, list_specifications, _ = next(os.walk(path2project))
    selected_specifications = st.multiselect(label="[User] Select folder(s):", options=list_specifications, default=None, key='2.selected_specifications')

    ## 读入硬盘里的电机优化结果的数据
    swarm_dict = {}
    if selected_specifications == []:
        st.error("Please select at least one specification.")
        st.stop()
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

    bool_filter_pareto_front = False
    if bool_filter_pareto_front == True:
        st.info(swarm_dict["PMSM_Q12p4y1_PEMD-2020"])
        st.info(swarm_dict["PMSM_Q12p4y1_PEMD-2020"].ad.swarm_data_file)
        with open(swarm_dict["PMSM_Q12p4y1_PEMD-2020"].ad.swarm_data_file, 'r', encoding='utf-8') as f:
            buf = f.read()
            data = json.loads('{'+buf[1:]+'}')
            del buf
        st.info(data.keys() )
        key = 'split_ratio'
        for item in data['2']['PMSM Q12p4y1 PEMD-2020-gen0-ind2']['Geometric parameters']:
            if key in item.keys():
                st.info(item)
                st.info(item[key]['value'])
        st.info(data['2']['PMSM Q12p4y1 PEMD-2020-gen0-ind2']['Geometric parameters'])


    ## 侧边栏 Sidebar
    # Show user selected mop's inputs
    with st.sidebar:
        st.title('User Configurations')

        ## 按照所选的电机规格，显示用户输入信息
        user_selected_folder = st.sidebar.selectbox('Select a folder to show its user configurations', selected_specifications, key='3.user_selected_folder')
        mop = swarm_dict[user_selected_folder]

        print(f'{user_selected_folder=}')
        st.sidebar.header('Specifications')
        st.sidebar.table(pd.DataFrame(data=list(mop.spec_input_dict.values()), index=list(mop.spec_input_dict.keys()), dtype="string", columns=['Value',]))
        st.sidebar.header('Simulation Settings')
        st.sidebar.table(pd.DataFrame(data=list(mop.fea_config_dict.values()), index=list(mop.fea_config_dict.keys()), dtype="string", columns=['Value',]))
    # """ DO NOT MODIFY ENDS """
    # """ DO NOT MODIFY ENDS """
    # """ DO NOT MODIFY ENDS """

    # 标签页
    tab_metrics, tab0, tab1, tab2, tab3 = st.tabs(["Metrics Dashboard", "Optimization Setup", "Population", "Select Individual", "Sensitivity Analysis"])

    with tab_metrics:
        st.subheader("Post-process Metrics Dashboard")
        metrics_df = _load_metrics_for_selected_folders(path2project, swarm_dict)
        if metrics_df.empty:
            st.warning("No ACMOP post-process metric records were found for the selected folders.")
        else:
            valid_only = st.checkbox("Use only valid metric rows", value=True, key="metrics.valid_only")
            if valid_only and "valid_metrics" in metrics_df.columns:
                metrics_df = metrics_df[metrics_df["valid_metrics"]].copy()

            available_specs = sorted(metrics_df["select_spec"].dropna().unique().tolist())
            selected_metric_specs = st.multiselect(
                "Filter specifications in dashboard",
                options=available_specs,
                default=available_specs,
                key="metrics.selected_specs",
            )
            if selected_metric_specs:
                metrics_df = metrics_df[metrics_df["select_spec"].isin(selected_metric_specs)].copy()

            st.caption(f"{len(metrics_df)} metric records loaded from existing ACMOP/JMAG result archives.")
            summary_df = postprocess_data.summarize_metric_archive(metrics_df)
            if not summary_df.empty:
                st.write("#### Archive Summary")
                st.dataframe(summary_df, use_container_width=True)

            table_columns = _metric_columns_available(metrics_df)
            st.download_button(
                "Download filtered metrics CSV",
                data=metrics_df[table_columns].to_csv(index=False).encode("utf-8-sig"),
                file_name="acmop_metrics_filtered.csv",
                mime="text/csv",
            )

            st.write("#### All Metrics")
            st.dataframe(metrics_df[table_columns], use_container_width=True, height=360)

            st.write("#### Selected Design")
            label_series = metrics_df.apply(
                lambda row: f"{row.get('select_spec', '')} | #{row.get('archive_record_no', '')} | {row.get('project_name', row.get('design_key', ''))}",
                axis=1,
            )
            selected_label = st.selectbox("Choose one design record", options=label_series.tolist(), key="metrics.selected_design")
            selected_row = metrics_df.loc[label_series[label_series == selected_label].index[0]]
            st.code(selected_row.get("design_key", ""), language="text")
            _show_selected_design_metrics(selected_row)

            col_left, col_right = st.columns(2)
            with col_left:
                _plot_metric_scatter(metrics_df)
            with col_right:
                _plot_ripple_scatter(metrics_df)
            _plot_loss_breakdown(selected_row)

    with tab0: # Optimization Setup
        mop = swarm_dict[user_selected_folder]
        st.subheader(f'Optimization Setup of {user_selected_folder}')
        mop.ad.acm_template.d['GP']

        st.write('#### Decision Variables:')
        convert_to_dict = {key: (val.type, val.value) for key, val in mop.ad.acm_template.d['GP'].items() if 'free' in val.type}; convert_to_df = pd.DataFrame(convert_to_dict)
        st.table(convert_to_df.T)
        st.write('#### Derived Variables:')
        convert_to_dict = {key: (val.type, val.value) for key, val in mop.ad.acm_template.d['GP'].items() if 'derived' in val.type}; convert_to_df = pd.DataFrame(convert_to_dict)
        st.table(convert_to_df.T)
        st.write('#### Fixed Variables:')
        convert_to_dict = {key: (val.type, val.value) for key, val in mop.ad.acm_template.d['GP'].items() if 'fixed' in val.type}; convert_to_df = pd.DataFrame(convert_to_dict)
        st.table(convert_to_df.T)


        st.write('#### Objectives (Non-dominated Sorting)')
        st.write(f'''
                 {mop.fea_config_dict['moo.fitness_OA']=}
                 {mop.fea_config_dict['moo.fitness_OB']=}
                 {mop.fea_config_dict['moo.fitness_OC']=}
        ''')

    with tab1: # Population

        ## 集群基本信息
        st.write('# 1. Swarms Information')
        ## 是否显示自动最优个体表格和帕累托前沿？
        optimal_xf_dict = None
        if st.checkbox("Show Swarm Table and Pareto Front?"):
            oc_values = [
                float(individual[-1])
                for mop_item in swarm_dict.values()
                for individual in mop_item.ad.analyzer.swarm_data_xf
                if len(individual) >= 3 and np.isfinite(individual[-1])
            ]
            minimum_oc = min(oc_values) if oc_values else None
            default_z_filter = 20.0
            if minimum_oc is not None and minimum_oc >= default_z_filter:
                default_z_filter = float(np.floor(max(oc_values)) + 1.0)
            z_filter = st.number_input(
                label='Input z-filter for Pareto front:',
                min_value=0.0,
                value=default_z_filter,
            )
            df_swarm, fig_Pareto, optimal_fitness_dict, optimal_xf_dict = utility_postprocess.inspect_swarm_and_show_table_plus_Pareto_front(swarm_dict, z_filter, output_dir=path2project)

            st.write('#### Highlight metric optima on Rank-1 Pareto front')
            if "pareto.highlight_metrics" not in st.session_state:
                st.session_state["pareto.highlight_metrics"] = []
            metric_keys = list(PARETO_HIGHLIGHT_METRICS)
            button_rows = [st.columns(5), st.columns(5)]
            for index, metric_key in enumerate(metric_keys):
                column = button_rows[index // 5][index % 5]
                if column.button(
                    PARETO_HIGHLIGHT_METRICS[metric_key]["label"],
                    key=f"pareto.button.{metric_key}",
                    use_container_width=True,
                ):
                    st.session_state["pareto.highlight_metrics"] = [metric_key]
            if button_rows[1][3].button("Show all", key="pareto.button.all", use_container_width=True):
                st.session_state["pareto.highlight_metrics"] = metric_keys
            if button_rows[1][4].button("Clear", key="pareto.button.clear", use_container_width=True):
                st.session_state["pareto.highlight_metrics"] = []

            highlighted_df = _highlight_pareto_metric_optima(
                fig_Pareto,
                swarm_dict,
                z_filter,
                st.session_state["pareto.highlight_metrics"],
            )
            st.pyplot(fig_Pareto)
            st.table(df_swarm)
            if not highlighted_df.empty:
                st.dataframe(highlighted_df, use_container_width=True, hide_index=True)
            elif st.session_state["pareto.highlight_metrics"]:
                st.warning("No Rank-1 individual satisfies the current OC filter and selected metric.")
            missing_specs = [
                mop_item.ad.select_spec
                for mop_item in swarm_dict.values()
                if mop_item.ad.select_spec not in optimal_fitness_dict
            ]
            if missing_specs:
                minimum_text = f"；当前数据中的最小 OC 为 {minimum_oc:.3f}" if minimum_oc is not None else ""
                st.warning(
                    f"以下规格在 OC < {z_filter:g} 的条件下没有帕累托个体："
                    f"{', '.join(missing_specs)}{minimum_text}。"
                )

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
        if optimal_xf_dict is not None and mop.ad.select_spec in optimal_xf_dict:
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


    with tab2:
        st.subheader('Select Individual')

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



    with tab3:
        st.subheader('Sensitivity Analysis')

        ## 之前选择的电机规格中 选择一个较优的设计 然后进行灵敏度分析（还需要一个选择）

        ## 选择分析变量
        selected_sensitivity_variables = st.multiselect(
                label='Select Sensitivity Variables',
                options=[] # TODO: 需要读swarm_dict 找到free的变量 拉出来表格选择 目前可以先不管
                # 暂时不用了 我们直接读GP里面的变量
            )

        ## 执行灵敏度分析并展示图像
        if st.button("Run Sensitivity Analysis"):
            # fig, ax = plt.subplots()
            # pmsm = sw.pmsm # 在读swarm_dict之后得到和im结构类似的pmsm 事实上应该就是变化的几个参数
            # class InitialDesign(object):
                # def __init__(self, pmsm, bounds=None):
                    # unit: mm 
                    # self.air_gap_length_delta           = pmsm.template.d['GP']['mm_d_sleeve'].value
                    # self.stator_tooth_width_b_ds        = pmsm.template.d['GP']['mm_mm_w_st'].value*1e3
                    # self.rotor_tooth_width_b_dr         = ( 2*np.pi*(pmsm.template.d['GP']['mm_r_ro'].value - pmsm.Length_HeadNeckRotorSlot)  - pmsm.Radius_of_RotorSlot * (2*Qr+2*pi) ) / Qr
                    # self.Angle_StatorSlotOpen           = pmsm.Angle_StatorSlotOpen # deg
                    # self.b1                             = pmsm.Width_RotorSlotOpen
                    # self.Width_StatorTeethHeadThickness = pmsm.Width_StatorTeethHeadThickness
                    # self.Length_HeadNeckRotorSlot       = pmsm.Length_HeadNeckRotorSlot

                    # self.design_parameters_denorm = [   self.air_gap_length_delta,
                                                        # self.stator_tooth_width_b_ds,
                                                        # self.rotor_tooth_width_b_dr,
                                                        # self.Angle_StatorSlotOpen,
                                                        # self.b1,
                                                        # self.Width_StatorTeethHeadThickness,
                                                        # self.Length_HeadNeckRotorSlot ]

                    # if bounds is None:
                        # self.design_parameters_denorm
                    # else:
                        # self.show_norm(bounds, self.design_parameters_denorm)


                # def show_denorm(self, bounds, design_parameters_norm):
                #     pop = design_parameters_norm
                #     min_b, max_b = np.asarray(bounds).T 
                #     diff = np.fabs(min_b - max_b)
                #     pop_denorm = min_b + pop * diff
                #     print('[De-normalized]:', end=' ')
                #     print(pop_denorm.tolist())
                    
                # def show_norm(self, bounds, design_parameters_denorm):
                #     min_b, max_b = np.asarray(bounds).T 
                #     diff = np.fabs(min_b - max_b)
                #     print(design_parameters_denorm)
                #     print(min_b)
                #     print(bounds)
                #     self.design_parameters_norm = (design_parameters_denorm - min_b)/diff #= pop
                #     # print type(self.design_parameters_norm)
                #     print('[Normalized]:', end=' ')
                #     print(self.design_parameters_norm.tolist())
                        
            # def local_sensitivity_analysis(self, specify_x_denorm=None):
                        # 敏感性检查：以基本设计为准，检查不同的参数取极值时的电机性能变化！这是最简单有效的办法。七个设计参数，那么就有14种极值设计。
                # if specify_x_denorm is None:
                    # build x_denorm for the template design
                    # x_denorm = self.ad.acm_template.build_x_denorm()
                # else:
                    # x_denorm = specify_x_denorm
                # print('[acmop.py] x_denorm:',  x_denorm)
                # print('[acmop.py] x_denorm_dict:', self.ad.acm_template.x_denorm_dict)
                # quit()
                # self.init_pop = []
                # if True:
                    # diff = np.array(self.fea_config_dict['local_sensitivity_analysis_diff_bounds'])
                    # min_b = np.array(self.fea_config_dict['local_sensitivity_analysis_min_bounds'])
                    # if specify_x_denorm is None:
                        # initial_design_denorm = np.array(utility.Pyrhonen_design(self.pmsm).design_parameters_denorm)
                    # else:
                        # initial_design_denorm = specified_initial_design_denorm
                    # initial_design = (initial_design_denorm - min_b) / diff
                    # print(initial_design_denorm.tolist())
                    # print(initial_design.tolist())
                    # base_design = initial_design.tolist()
                    # print('base_design:', base_design, '\n-------------')
                    # quit()
                    # number_of_variants = self.fea_config_dict['local_sensitivity_analysis_number_of_variants']
                    # self.init_pop = [initial_design] # include initial design!
                    # for i in range(len(base_design)): # 10 design parameters
                        # for j in range(number_of_variants+1): # 21 variants interval
                            # copy list
                            # design_variant = base_design[::]
                            # design_variant[i] = j * 1./number_of_variants
                            # self.init_pop.append(design_variant)
                    # for ind, el in enumerate(self.init_pop):
                        # print(ind)
                        # print(el)
                # return self.init_pop
            mop = swarm_dict[folder]
            ad = mop.acm_template.d['GP']
            print(ad)
            # acmop.reproduce_design_from_design_parameters(ad.design_parameters_denorm)
            
            # for ind, el in enumerate(sensitivity_denorm):
                # sensitivity_denorm = local_sensitivity_analysis(specified_initial_design_denorm=specify_x_denorm)
                # acmop.reproduce_design_from_design_parameters(el)

            for ind,el in enumerate():
                pass
                # TODO: plot
                # TODO: 1 number of x_denorm is 10 while optimization is 5
                # TODO: 2 The call of JMAG need to be fixed # done

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
