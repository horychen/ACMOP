# Backend Dependency Report

This report serves to help with refactoring the codebase away from the global object pattern. It lists the state variables and the functions that depend on them, followed by an exhaustive list of functions and their implicit inputs/outputs.

## Variable Directory (Categorized by where they are needed)

This section lists global/state variables and the functions that read or write to them, making it easy to see the dependency graph of each variable.

### Variables: `DisplacementAngle_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `EX['DriveW_zQ']`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `EX['Js']`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `EX['WindingFill']`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `EX['wily']`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `EX['wily'].number_parallel_branch`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `ForConX_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
  - `utility.py::read_csv_results_4_comparison__transient`
  - `utility.py::read_csv_results_4_comparison_eddycurrent`

### Variables: `ForConY_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
  - `utility.py::read_csv_results_4_comparison__transient`
  - `utility.py::read_csv_results_4_comparison_eddycurrent`

### Variables: `O1_max.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `O1_min.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `O2_ecce_data.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `O2_max.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `O2_min.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `O2_prototype_data.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `TorCon_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
  - `utility.py::read_csv_results_4_comparison__transient`
  - `utility.py::read_csv_results_4_comparison_eddycurrent`

### Variables: `acm_variant.rotorMagnet.notched_rotor.p`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_spmsm`

### Variables: `acm_variant.wily`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `acm_variant.winding.wily`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `app`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.draw_jmag_model`
  - `JMAG.py::JMAG.mesh_study`
  - `JMAG.py::JMAG.open`
  - `JMAG.py::JMAG.pre_process_PMSM`
  - `JMAG.py::JMAG.run_study`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.open`

### Variables: `app.CreateGeometryEditor`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_jmag_model`

### Variables: `app.ExportImageWithSize`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.mesh_study`
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `app.FunctionFactory`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `app.GetCurrentModel`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_jmag_model`

### Variables: `app.GetDataManager`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `app.GetMaterialLibrary`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`

### Variables: `app.Hide`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`

### Variables: `app.LaunchGeometryEditor`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_jmag_model`

### Variables: `app.NewProject`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`

### Variables: `app.NumModels`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_jmag_model`

### Variables: `app.Save`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.run_study`

### Variables: `app.SaveAs`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`

### Variables: `app.SetCurrentStudy`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `app.Show`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`

### Variables: `app.ShowCircuitGrid`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `app.VersionString`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`

### Variables: `app.View`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.mesh_study`
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `basic_info.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
  - `utility.py::read_csv_results_4_comparison__transient`

### Variables: `coil_flux_linkage_peak2peak_value_results.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.build_str_results`

### Variables: `created_parts.append`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.add_radial_array`

### Variables: `data['wily']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.from_dict_full`

### Variables: `data_max.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `data_min.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `errors.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.validate_parameters`

### Variables: `filtered_O_list.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `filtered_filtered_O_list.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `filtered_index_list.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `filtered_loss_list.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `filtered_torque_list.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `filtered_x.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `filtered_y.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `fitness_mapping`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.build_str_results`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.build_str_results`

### Variables: `full_dict['wily']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `func_wrapper`
- **Read by (Depends on this)**:
  - `utility.py::blockPrinting`

### Variables: `im_variant.winding.EX['RatedSpeed']`
- **Read by (Depends on this)**:
  - `utility.py::get_windage_loss`

### Variables: `index_list.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `instance.winding.wily`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.from_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.from_dict_full`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.from_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.from_dict_full`

### Variables: `key_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
  - `utility.py::read_csv_results_4_comparison__transient`

### Variables: `kp_cjh_list.append`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw_per_phase`

### Variables: `kp_els_list.append`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw_per_phase`

### Variables: `l.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.list_cost_function`

### Variables: `l_ForCon_X.append`
- **Read by (Depends on this)**:
  - `utility.py::check_csv_results_4_general_purpose`

### Variables: `l_ForCon_Y.append`
- **Read by (Depends on this)**:
  - `utility.py::check_csv_results_4_general_purpose`

### Variables: `l_TorCon.append`
- **Read by (Depends on this)**:
  - `utility.py::check_csv_results_4_general_purpose`

### Variables: `l_cost_function.append`
- **Read by (Depends on this)**:
  - `utility.py::fobj_list`

### Variables: `l_rated_stack_length.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `l_rated_total_loss.append`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`

### Variables: `l_slip_freq.append`
- **Read by (Depends on this)**:
  - `utility.py::check_csv_results_4_general_purpose`

### Variables: `list_region_objects.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.prepareSection`

### Variables: `list_regions.append`
- **Read by (Depends on this)**:
  - `machine_geometry.py::draw_instruction_parser`
  - `machine_geometry_utils.py::draw_instruction_parser`

### Variables: `list_xy_magnets.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `lst_x.append`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`

### Variables: `lst_y.append`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`

### Variables: `max_idxs.append`
- **Read by (Depends on this)**:
  - `utility.py::max_indices_2`

### Variables: `more_info.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.learn_about_the_archive`

### Variables: `motor.rotor.OD`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `motor.rotor.airGap`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `motor.stator.ID`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `motor.stator.OD`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `motor.stator.liner`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `motor.stator.toothDepth`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `motor.stator.toothWidth`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `motor.stator.yoke`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `new_key_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `out.append`
- **Read by (Depends on this)**:
  - `utility.py::to_precision`

### Variables: `reformat_wily_info`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `result.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.get_metric_of_the_whole_swarm`

### Variables: `results.append`
- **Read by (Depends on this)**:
  - `BH/Inspect_BH_curve.py::collect_and_convert_bh_curves`

### Variables: `results_for_refining_bounds['O1'].append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `results_for_refining_bounds['O2'].append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `rotor.ID`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `rotor.OD`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `rotor.magnetDepth`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `rotor_Joule_loss_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `rows.append`
- **Read by (Depends on this)**:
  - `BH/Inspect_BH_curve.py::read_bh_data`

### Variables: `self.Angle_StatorSlotOpen`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.CommutatingSequenceB`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.CommutatingSequenceD`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.Cost`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.Cost_Cu`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.Cost_Fe`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.Cost_PM`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.Current_dict`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.Current_dict['Time(s)']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.DisplacementAngle_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.Ea`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.Ea.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.Em`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.Em.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.ExcitationFreqSimulated`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.FRW`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.FRW.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.FluxLinkage_dict`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.ForConAbs_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.ForConX_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.ForConY_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.GP`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`
  - `modern_machine_designer_utility.py::Geometry.__init__`
  - `modern_machine_designer_utility.py::Geometry.__repr__`
  - `modern_machine_designer_utility.py::Geometry.print_parameters`
  - `modern_machine_designer_utility.py::Geometry.update_from_GP`
- **Written by (Modifies this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`
  - `modern_machine_designer_utility.py::Geometry.__init__`

### Variables: `self.GP.items`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.__init__`
  - `modern_machine_designer_utility.py::Geometry.__repr__`
  - `modern_machine_designer_utility.py::Geometry.print_parameters`
  - `modern_machine_designer_utility.py::Geometry.update_from_GP`

### Variables: `self.HP`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`
- **Written by (Modifies this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['1']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['1']['0']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['1']['1']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['2']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['2']['0']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['2']['1']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['3']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['3']['0']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['3']['1']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['4']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['4']['0']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['4']['1']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['6']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['6']['0']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['6']['1']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['7']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['7']['0']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['7']['1']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['9']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['9']['0']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP['9']['1']`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.HP_mirror`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`
- **Written by (Modifies this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.JMAG_version_number`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.open`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.open`

### Variables: `self.JMAG_version_string`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.open`

### Variables: `self.La`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_bemf_constants`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_motor_losses`
- **Written by (Modifies this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.Length_HeadNeckRotorSlot`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.Omega`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.PowerFactor`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.Q`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw_per_phase`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.Qr`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.Qs`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `modern_machine_designer_utility.py::Winding.get_winding_factor`
  - `modern_machine_designer_utility.py::Winding.to_dict`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.RP`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`
- **Written by (Modifies this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.RatedEfficiency`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.RatedStkLen`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.RatedStkLen.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.RatedVol`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.RatedVol.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.RatedWeight`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.RatedWeight.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.SI`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_jmag_model`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.draw_jmag_model`

### Variables: `self.SIPNV_or_SEPA`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.SIesign_display_generator`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.SIesign_parameters_denorm`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.SIesign_parameters_generator`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.get_best_generation`
  - `utility.py::SwarmDataAnalyzer.my_population_distribution_plots`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.SIesign_parameters_norm`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.show_norm`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.show_norm`

### Variables: `self.SIesign_parameters_norm.tolist`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.show_norm`

### Variables: `self.SIir_run`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.SPP`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.TRV`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.TorCon_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.TorqueRipple`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.Trip`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.Trip.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.Tripple`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.Width_StatorTeethHeadThickness`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.__class__`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.__class__.__module__`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.__class__.__name__`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.__dict__`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.print_parameters`
  - `modern_machine_designer_utility.py::Geometry.to_dict`

### Variables: `self.__dict__.items`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.print_parameters`
  - `modern_machine_designer_utility.py::Geometry.to_dict`

### Variables: `self._extract_performance_lists`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self._get_parameter_logger`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`

### Variables: `self._initialize_empty`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self._load_from_json`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.__init__`

### Variables: `self._load_from_raw`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.__init__`

### Variables: `self._next_index`
- **Read by (Depends on this)**:
  - `machine_geometry.py::MachineGeometry.add_part`
  - `machine_geometry_utils.py::MachineGeometry.add_part`

### Variables: `self.accumSquaredData`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.add_circuit`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `self.add_material`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `self.add_part`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.add_radial_array`

### Variables: `self.add_plots`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.build_str_results`

### Variables: `self.air_gap`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.air_gap_length_delta`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.all_points`
- **Read by (Depends on this)**:
  - `machine_geometry.py::Coil.draw_instruction`
  - `machine_geometry.py::MachineGeometry.add_part`
  - `machine_geometry.py::MachineGeometry.draw_machine_using_CairoDrawer`
  - `machine_geometry.py::MachineGeometry.sync`
  - `machine_geometry.py::Magnet.draw_instruction`
  - `machine_geometry.py::RotorCore.draw_instruction`
  - `machine_geometry.py::StatorCore.draw_instruction`
  - `machine_geometry_utils.py::Coil.draw_instruction`
  - `machine_geometry_utils.py::MachineGeometry.add_part`
  - `machine_geometry_utils.py::Magnet.draw_instruction`
  - `machine_geometry_utils.py::RotorCore.draw_instruction`
  - `machine_geometry_utils.py::StatorCore.draw_instruction`
- **Written by (Modifies this)**:
  - `machine_geometry.py::MachineGeometry.sync`

### Variables: `self.all_points.pole_count`
- **Read by (Depends on this)**:
  - `machine_geometry.py::Magnet.draw_instruction`
  - `machine_geometry.py::RotorCore.draw_instruction`
  - `machine_geometry_utils.py::Magnet.draw_instruction`
  - `machine_geometry_utils.py::RotorCore.draw_instruction`

### Variables: `self.all_points.slot_count`
- **Read by (Depends on this)**:
  - `machine_geometry.py::Coil.draw_instruction`
  - `machine_geometry.py::StatorCore.draw_instruction`
  - `machine_geometry_utils.py::Coil.draw_instruction`
  - `machine_geometry_utils.py::StatorCore.draw_instruction`

### Variables: `self.ampl`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.app`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.checkGeomApp`
  - `JMAG.py::JMAG.close`
  - `JMAG.py::JMAG.open`
  - `JMAG.py::JMAG.save`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.open`

### Variables: `self.app.CreateGeometryEditor`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.checkGeomApp`

### Variables: `self.app.GetCurrentModel`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.save`

### Variables: `self.app.LaunchGeometryEditor`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.checkGeomApp`

### Variables: `self.app.Quit`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.close`

### Variables: `self.ass`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.getSketch`
  - `JMAG.py::JMAG.regionMirrorCopy`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.ass.CreateSketch`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.ass.GetItem`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.getSketch`
  - `JMAG.py::JMAG.regionMirrorCopy`

### Variables: `self.avg_val_of_turn_func`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.awg`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_motor_losses`
- **Written by (Modifies this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.b1`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.bFillRegion`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.bMirror`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.draw_spmsm`
  - `JMAG.py::JMAG.prepareSection`
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.draw_spmsm`
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.basic_info`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.best_design_denorm`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.best_design_denorm['0']`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.best_design_display`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.best_design_display.split`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.bool_3PhaseCurrentSource`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.bool_CustomizedCircuit`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.bool_DPNVorSEPA`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.bool_PermanentMagnet`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.bool_RotorNotched`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.bool_StatorSlotClosed`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.bool_distributed_or_concentrated`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.bool_double_layer_winding`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.bool_initialized`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.bool_jmagDeleteResultsAfterCalculation`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.bool_suppressShaft`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.bounds`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
  - `modern_machine_designer_utility.py::Parameter.to_dict`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.bounds['0']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.bounds['1']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.buf`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`
  - `utility.py::SwarmDataAnalyzer.design_display_generator`
  - `utility.py::SwarmDataAnalyzer.design_parameters_generator`
  - `utility.py::SwarmDataAnalyzer.find_individual`
  - `utility.py::SwarmDataAnalyzer.get_certain_objective_function`
  - `utility.py::SwarmDataAnalyzer.get_list_objective_function`
  - `utility.py::SwarmDataAnalyzer.get_windage_loss`
  - `utility.py::SwarmDataAnalyzer.list_cost_function`
  - `utility.py::SwarmDataAnalyzer.list_generations`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.buf_length`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.build_basic_info`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.calc`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.calc_bemf_constants`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`

### Variables: `self.calc_bounds`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.calc_motor_losses`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`

### Variables: `self.calculate_excitation_current`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_spmsm`

### Variables: `self.calculate_slot_area`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.analyze_thermal_performance`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`

### Variables: `self.checkGeomApp`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.circuit_current`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.coeff`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.coil_fluxLinkage`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.coil_pitch`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.coil_pitch_y`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `modern_machine_designer_utility.py::Winding.get_winding_factor`
  - `modern_machine_designer_utility.py::Winding.to_dict`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw_per_phase`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.color`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.__init__`
  - `modern_machine_designer_utility.py::Geometry.__repr__`
  - `modern_machine_designer_utility.py::Geometry.print_parameters`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Geometry.__init__`

### Variables: `self.components_make_region`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.draw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Geometry.draw`

### Variables: `self.connection_star_raw_dict`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.consts`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.cosine`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.count`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.counter`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.ctx`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.apply_stroke`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.ctx.arc`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.arc_negative`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.fill_preserve`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.line_to`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.move_to`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.new_path`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.paint`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.ctx.path_extents`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.restore`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.rotate`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.save`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.scale`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.set_line_cap`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.apply_stroke`

### Variables: `self.ctx.set_line_width`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.apply_stroke`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.set_source_rgb`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.ctx.set_source_rgba`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.apply_stroke`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.stroke`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.apply_stroke`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.ctx.transform`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.d_air_gap`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.split_ratio`

### Variables: `self.d_tooth`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.d_stator_yoke`

### Variables: `self.decode_py_reduce_ordered_dict`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.filter_data`

### Variables: `self.defaultUnit`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.deg_alpha_st`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.deg_alpha_st.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.deg_winding_U_phase_phase_axis_angle`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.degree_between_slots`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.derivation`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.to_dict`

### Variables: `self.derivation.__dict__`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.to_dict`

### Variables: `self.derivation.__dict__.items`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.to_dict`

### Variables: `self.dict_coil_connection`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.dict_suspension_kw_cjh`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.dict_suspension_kw_els`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.dict_suspension_kw_els['A_angle']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.dict_suspension_kw_els['B_angle']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.dict_suspension_kw_els['C_angle']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.dict_torque_kw_cjh`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.dict_torque_kw_els`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.distributed_or_concentrated`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.dm`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.build_str_results`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.build_str_results`

### Variables: `self.doc`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.addConstraintCocentricity`
  - `JMAG.py::JMAG.checkGeomApp`
  - `JMAG.py::JMAG.getSketch`
  - `JMAG.py::JMAG.pre_process_PMSM`
  - `JMAG.py::JMAG.prepareSection`
  - `JMAG.py::JMAG.regionCircularPattern360Origin`
  - `JMAG.py::JMAG.regionMirrorCopy`
  - `JMAG.py::JMAG.save`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.checkGeomApp`
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.doc.CreateReferenceFromItem`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.addConstraintCocentricity`
  - `JMAG.py::JMAG.getSketch`
  - `JMAG.py::JMAG.regionCircularPattern360Origin`
  - `JMAG.py::JMAG.regionMirrorCopy`

### Variables: `self.doc.GetAssembly`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.doc.GetSelection`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.pre_process_PMSM`
  - `JMAG.py::JMAG.prepareSection`

### Variables: `self.doc.SaveModel`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.save`

### Variables: `self.dpnv_grouping_dict_a`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.dpnv_grouping_dict_a['GAC']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`

### Variables: `self.dpnv_grouping_dict_a['GBD']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`

### Variables: `self.dpnv_grouping_dict_b`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.dpnv_grouping_dict_b['GAC']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`

### Variables: `self.dpnv_grouping_dict_b['GBD']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`

### Variables: `self.dpnv_grouping_dict_c`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.dpnv_grouping_dict_c['GAC']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`

### Variables: `self.dpnv_grouping_dict_c['GBD']`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`

### Variables: `self.draw_function`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.__init__`
  - `modern_machine_designer_utility.py::Geometry.draw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Geometry.__init__`

### Variables: `self.draw_machine_using_CairoDrawer`
- **Read by (Depends on this)**:
  - `machine_geometry.py::MachineGeometry.show_geometry_svg`

### Variables: `self.edge4Ref`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.prepareSection`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.edge4ref`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.prepareSection`

### Variables: `self.estimate_max_wires_in_slot`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.analyze_thermal_performance`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`

### Variables: `self.f1`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.f2`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.f3`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.fea_config_dict`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.build_str_results`
  - `JMAG.py::JMAG.open`
  - `JMAG.py::JMAG.pre_process_PMSM`
  - `modern_machine_designer_utility.py::swarm_data_container.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container.__init__`

### Variables: `self.fea_config_dict['designer.StepPerCycle_3rdTSS']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.build_str_results`

### Variables: `self.fea_config_dict['designer.max_nonlinear_iteration']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `self.fea_config_dict['designer.show']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `self.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`

### Variables: `self.fea_config_dict['pc_name']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.open`

### Variables: `self.femm_loss_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.fig_plot2piFft`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.plot2piFft`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.plot2piFft`

### Variables: `self.fig_plotFuncObj`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.plotFuncObj`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.plotFuncObj`

### Variables: `self.filename`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.filter_data`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`

### Variables: `self.flag_do_not_evaluate_when_init_pop`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`

### Variables: `self.flag_material_already_loaded`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.open`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.open`

### Variables: `self.force_abs`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_ang`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_ang.append`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_err_abs`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_err_ang`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_err_ang_new_way`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_err_ang_old_way`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_error_angle`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_x`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.force_y`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.geomApp`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.checkGeomApp`
  - `JMAG.py::JMAG.getSketch`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.checkGeomApp`
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.geomApp.GetDocument`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.geomApp.NewDocument`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.checkGeomApp`

### Variables: `self.geometry`
- **Read by (Depends on this)**:
  - `machine.py::Machine.__init__`
  - `machine.py::Machine.draw_jmag`
  - `machine.py::Machine.draw_svg`
  - `machine.py::Machine.sync`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`
- **Written by (Modifies this)**:
  - `machine.py::Machine.__init__`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`

### Variables: `self.geometry.add_part`
- **Read by (Depends on this)**:
  - `machine.py::Machine.sync`

### Variables: `self.geometry.all_points`
- **Read by (Depends on this)**:
  - `machine.py::Machine.draw_jmag`

### Variables: `self.geometry.machineGeometry`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`

### Variables: `self.geometry.machineGeometry.items`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`

### Variables: `self.geometry.machineGeometry.values`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`

### Variables: `self.geometry.parts`
- **Read by (Depends on this)**:
  - `machine.py::Machine.draw_jmag`

### Variables: `self.geometry.show_geometry_svg`
- **Read by (Depends on this)**:
  - `machine.py::Machine.draw_svg`

### Variables: `self.geometry.sync`
- **Read by (Depends on this)**:
  - `machine.py::Machine.sync`

### Variables: `self.getSketch`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.drawArc`
  - `JMAG.py::JMAG.drawCircle`
  - `JMAG.py::JMAG.drawLine`

### Variables: `self.get_certain_objective_function`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_torque_force`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.get_complex_number_kw`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.get_complex_number_kw_per_phase`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`

### Variables: `self.get_complex_number_winding_factor_of_coil_i`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw_per_phase`

### Variables: `self.get_free_variables`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_free_variable_bounds_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_free_variables_as_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`

### Variables: `self.get_metric_of_the_whole_swarm`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.get_parameter`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.set_free_variables_from_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.set_parameter_value`

### Variables: `self.get_parameter_dict_by_name`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`

### Variables: `self.get_parameter_fields`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_free_variables`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_parameter`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_parameter_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_parameter_dict_by_name`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_parameters_by_type`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_parameters_summary`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.list_parameters`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_derived_parameters`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.validate_parameters`

### Variables: `self.get_parameters_by_type`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.list_parameters`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`

### Variables: `self.get_parameters_summary`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__repr__`

### Variables: `self.get_rotor_volume`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_rotor_weight`

### Variables: `self.get_voltage_and_current`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.get_winding_factor`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.get_wire_properties`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.analyze_thermal_performance`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_motor_losses`

### Variables: `self.gp`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.d_air_gap`
  - `machine_geometry_utils.py::MachineGeometry.d_magnet`
  - `machine_geometry_utils.py::MachineGeometry.d_tooth`
  - `machine_geometry_utils.py::MachineGeometry.d_tooth_shoe`
  - `machine_geometry_utils.py::MachineGeometry.r_rotor_outer`
  - `machine_geometry_utils.py::MachineGeometry.r_shaft`
  - `machine_geometry_utils.py::MachineGeometry.r_stator_outer`
  - `machine_geometry_utils.py::MachineGeometry.sync`
  - `machine_geometry_utils.py::MachineGeometry.tooth_shape`
  - `machine_geometry_utils.py::MachineGeometry.w_tooth`
- **Written by (Modifies this)**:
  - `machine_geometry_utils.py::MachineGeometry.sync`

### Variables: `self.gp['d_air_gap']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.d_air_gap`

### Variables: `self.gp['d_magnet']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.d_magnet`

### Variables: `self.gp['d_tooth']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.d_tooth`

### Variables: `self.gp['d_tooth_shoe']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.d_tooth_shoe`

### Variables: `self.gp['r_rotor_outer']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.r_rotor_outer`

### Variables: `self.gp['r_shaft']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.r_shaft`

### Variables: `self.gp['r_stator_outer']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.r_stator_outer`

### Variables: `self.gp['tooth_shape']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.tooth_shape`

### Variables: `self.gp['w_tooth']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.w_tooth`

### Variables: `self.grouping_AC`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.hex_to_rgb`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`

### Variables: `self.horizontal_position`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`
- **Written by (Modifies this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.iRotateCopy`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.draw_spmsm`
  - `JMAG.py::JMAG.prepareSection`
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.draw_spmsm`
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.id`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.id_rotorCore`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.pre_process_PMSM`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `self.id_statorCore`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.pre_process_PMSM`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `self.im`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.DriveW_poles`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.Qr`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.design_parameters`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.design_parameters['2']`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.rotor_slot_height_h_sr`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.stack_length`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.template`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.template.SI`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.template.SI['GP']`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.template.SI['GP']['mm_r_ro']`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.im.template.SI['GP']['mm_r_ro'].value`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.imag`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.initial_excitation_bias_compensation_deg`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.initialized`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.jd`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.jmag_loss_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.k`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.kd1`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.kp1`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.kw1`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `modern_machine_designer_utility.py::Winding.to_dict`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.l21`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l22`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l41`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l42`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l_FRW`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_OA`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_OB`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.get_list_y_data`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_OC`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_TRV`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.get_list_y_data`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_design_parameters`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_efficiency`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_force_error_angle`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.get_list_y_data`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_leftlayer1`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l_leftlayer2`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l_normalized_force_error_magnitude`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_normalized_torque_ripple`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_original_rotor_weight`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_original_stack_length`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_power_factor`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_efficiency`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_iron_loss`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_magnet_Joule_loss`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.l_rated_rotor_copper_loss_along_stack`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_rotor_volume`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_rotor_weight`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_shaft_power`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_stack_length`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_stator_copper_loss_along_stack`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_total_loss`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rated_windage_loss`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_rightlayer1`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l_rightlayer2`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.l_rotor_copper_loss_in_end_turn`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_ss_avg_force_magnitude`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_stack`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.aspect_ratio`

### Variables: `self.l_stator_copper_loss_in_end_turn`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.l_torque_average`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.layer_A1`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.layer_A2`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.layer_B1`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.layer_B2`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.layer_X_phases`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.layer_X_signs`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.layer_Y_phases`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.layer_Y_signs`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.liner`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `WireSlot_v2.py::MotorThermalAnalyzer.estimate_max_wires_in_slot`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.estimate_max_wires_in_slot`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.list_cost_function`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.get_best_generation`

### Variables: `self.list_layer_motor_phases`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.list_layer_motor_signs`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.list_layer_suspension_phases`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.list_layer_suspension_signs`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.list_phase_u_slot_number`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.list_phase_v_slot_number`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.list_phase_w_slot_number`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.list_rotor_current_amp`
- **Read by (Depends on this)**:
  - `utility.py::get_copper_loss_Bolognani`

### Variables: `self.list_slot_number_of_phase`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.lst_x`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`

### Variables: `self.lst_y`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`

### Variables: `self.lst_y.index`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`

### Variables: `self.m`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `modern_machine_designer_utility.py::Winding.get_winding_factor`
  - `modern_machine_designer_utility.py::Winding.to_dict`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.machine_data`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.machine_data.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.mec_power`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.message`
- **Read by (Depends on this)**:
  - `utility.py::ExceptionBadNumberOfParts.__init__`
  - `utility.py::ExceptionBadNumberOfParts.__str__`
  - `utility.py::ExceptionReTry.__init__`
  - `utility.py::ExceptionReTry.__str__`
- **Written by (Modifies this)**:
  - `utility.py::ExceptionBadNumberOfParts.__init__`
  - `utility.py::ExceptionReTry.__init__`

### Variables: `self.mm_r_ro`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_rotor_volume`

### Variables: `self.mm_r_ro.value`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_rotor_volume`

### Variables: `self.mm_r_si`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.mm_r_si.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.mm_w_st`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.mm_w_st.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.model`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.motor`
- **Read by (Depends on this)**:
  - `machine.py::Machine.__init__`
  - `machine.py::Machine.sync`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
- **Written by (Modifies this)**:
  - `machine.py::Machine.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.motor.stator`
- **Read by (Depends on this)**:
  - `machine.py::Machine.sync`

### Variables: `self.motor.stator.tooth_shape`
- **Read by (Depends on this)**:
  - `machine.py::Machine.sync`

### Variables: `self.my_scatter_plot`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_torque_force`

### Variables: `self.mycurrent`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.mytime`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.myvoltage`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.name`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.__init__`
  - `modern_machine_designer_utility.py::Geometry.__repr__`
  - `modern_machine_designer_utility.py::Geometry.print_parameters`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility._get_parameter_logger`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
  - `modern_machine_designer_utility.py::Parameter.__init__`
  - `modern_machine_designer_utility.py::Parameter.__repr__`
  - `modern_machine_designer_utility.py::Parameter.to_dict`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Geometry.__init__`
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.no_winding_layer`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`

### Variables: `self.normalized_force_error_magnitude`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.number_of_chromosome`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`

### Variables: `self.number_of_designs`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`
  - `utility.py::SwarmDataAnalyzer.design_display_generator`
  - `utility.py::SwarmDataAnalyzer.design_parameters_generator`
  - `utility.py::SwarmDataAnalyzer.find_individual`
  - `utility.py::SwarmDataAnalyzer.get_certain_objective_function`
  - `utility.py::SwarmDataAnalyzer.get_list_objective_function`
  - `utility.py::SwarmDataAnalyzer.get_windage_loss`
  - `utility.py::SwarmDataAnalyzer.list_cost_function`
  - `utility.py::SwarmDataAnalyzer.list_generations`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.number_of_free_variables`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.number_of_parallel_branch`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `modern_machine_designer_utility.py::Winding.to_dict`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.number_of_winding_layer`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`

### Variables: `self.number_parallel_branch`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.number_winding_layer`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.options`
- **Read by (Depends on this)**:
  - `machine_geometry.py::RotorCore.draw_instruction`
  - `machine_geometry_utils.py::RotorCore.draw_instruction`
  - `machine_geometry_utils.py::StatorCore.draw_instruction`

### Variables: `self.overwritten`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.ox_distribution_phase_U`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.ox_distribution_three_phase`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.p`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `modern_machine_designer_utility.py::Winding.get_winding_factor`
  - `modern_machine_designer_utility.py::Winding.to_dict`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.pairs`
- **Read by (Depends on this)**:
  - `winding_layout.py::pole_specific_winding_with_neutral.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::pole_specific_winding_with_neutral.__init__`

### Variables: `self.parameter_dict`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.parameter_dict_by_name`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.update_geometric_parameters`

### Variables: `self.parameter_dict_by_name.get`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.apply_parameter_dict`

### Variables: `self.parts`
- **Read by (Depends on this)**:
  - `machine_geometry.py::MachineGeometry.add_part`
  - `machine_geometry.py::MachineGeometry.draw_machine_using_CairoDrawer`
  - `machine_geometry_utils.py::MachineGeometry.add_part`

### Variables: `self.parts.append`
- **Read by (Depends on this)**:
  - `machine_geometry.py::MachineGeometry.add_part`
  - `machine_geometry_utils.py::MachineGeometry.add_part`

### Variables: `self.path2SwarmData`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.draw_individual_from_swarm`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.write_swarm_survivor`
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`
  - `modern_machine_designer_utility.py::Winding.plot_winding_function`

### Variables: `self.payload`
- **Read by (Depends on this)**:
  - `utility.py::ExceptionBadNumberOfParts.__init__`
  - `utility.py::ExceptionReTry.__init__`
- **Written by (Modifies this)**:
  - `utility.py::ExceptionBadNumberOfParts.__init__`
  - `utility.py::ExceptionReTry.__init__`

### Variables: `self.phase`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.pole_count`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::MachineGeometry.deg_alpha_rm`
- **Written by (Modifies this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.poles`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_bemf_constants`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_motor_losses`
- **Written by (Modifies this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.pr`
- **Read by (Depends on this)**:
  - `winding_layout.py::winding_layout_v2.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::winding_layout_v2.__init__`

### Variables: `self.prepareSection`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_spmsm`

### Variables: `self.print_out_string`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`

### Variables: `self.projName`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.project_names`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.project_names.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.ps`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `modern_machine_designer_utility.py::Winding.get_winding_factor`
  - `modern_machine_designer_utility.py::Winding.to_dict`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Winding.__init__`
  - `winding_layout.py::winding_layout_v2.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.q`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.q2`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.qs`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.r_rotor_outer`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.d_stator_yoke`
  - `machine_geometry_utils.py::MachineGeometry.split_ratio`

### Variables: `self.r_stator_outer`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.aspect_ratio`
  - `machine_geometry_utils.py::MachineGeometry.d_stator_yoke`
  - `machine_geometry_utils.py::MachineGeometry.split_ratio`

### Variables: `self.radian_between_slots`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.range_ss`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.rated_data`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.rated_data.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.rated_speed`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`
- **Written by (Modifies this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.read_csv_results_4_general_purpose`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.build_str_results`

### Variables: `self.real`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.reference_data`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.reference_data['2']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.reference_data['3']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.reference_data['4']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.reference_data['5']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.reference_data['6']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.reference_design`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.__init__`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.reference_design['3']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.reference_design['3'].split`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.regionCircularPattern360Origin`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.prepareSection`

### Variables: `self.regionMirrorCopy`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.prepareSection`

### Variables: `self.required_torque`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.rotated_position`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`
- **Written by (Modifies this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.rotor_od`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.rotor_tooth_width_b_dr`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.rotor_volume`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_torque_force`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.rotor_weight`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container.sensitivity_bar_charts`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_torque_force`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.run_integer`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.save`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_spmsm`

### Variables: `self.scale`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.prepareSection`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.scalingFactor`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.select_FEA_tool`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.select_fea_config_dict`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.setSymPos`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.setTurnFuncObject`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.show`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.draw_spmsm`

### Variables: `self.show_geometry`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.draw_individual_from_swarm`

### Variables: `self.show_norm`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.sine`
- **Read by (Depends on this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Goertzel_Data_Struct.__init__`

### Variables: `self.sketch`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.addConstraintCocentricity`
  - `JMAG.py::JMAG.drawArc`
  - `JMAG.py::JMAG.drawCircle`
  - `JMAG.py::JMAG.drawLine`
  - `JMAG.py::JMAG.getSketch`
  - `JMAG.py::JMAG.prepareSection`
  - `JMAG.py::JMAG.regionCircularPattern360Origin`
  - `JMAG.py::JMAG.regionMirrorCopy`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.drawArc`
  - `JMAG.py::JMAG.drawCircle`
  - `JMAG.py::JMAG.drawLine`
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.sketch.CloseSketch`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.prepareSection`

### Variables: `self.sketch.CreateArc`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.drawArc`

### Variables: `self.sketch.CreateBiConstraint`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.addConstraintCocentricity`

### Variables: `self.sketch.CreateCircle`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.drawCircle`

### Variables: `self.sketch.CreateLine`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.drawLine`

### Variables: `self.sketch.CreateRegionCircularPattern`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.regionCircularPattern360Origin`

### Variables: `self.sketch.CreateRegionMirrorCopy`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.regionMirrorCopy`

### Variables: `self.sketch.CreateRegions`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.prepareSection`

### Variables: `self.sketch.GetItem`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.addConstraintCocentricity`
  - `JMAG.py::JMAG.prepareSection`
  - `JMAG.py::JMAG.regionMirrorCopy`

### Variables: `self.sketch.OpenSketch`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.drawArc`
  - `JMAG.py::JMAG.drawCircle`
  - `JMAG.py::JMAG.drawLine`
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.sketch.SetProperty`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.sketchNameList`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `JMAG.py::JMAG.getSketch`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.sketchNameList.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.getSketch`

### Variables: `self.sketch_color`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.convert_to_pdf`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.convert_to_pdf`

### Variables: `self.slot_count`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::MachineGeometry.alpha_stator_tooth_span`
- **Written by (Modifies this)**:
  - `machine_geometry.py::AllPoints.__post_init__`
  - `machine_geometry_utils.py::AllPoints.__post_init__`

### Variables: `self.slot_per_phase`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.slots`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `WireSlot_v2.py::MotorThermalAnalyzer.analyze_thermal_performance`
  - `WireSlot_v2.py::MotorThermalAnalyzer.calculate_slot_area`
  - `WireSlot_v2.py::MotorThermalAnalyzer.estimate_max_wires_in_slot`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_bemf_constants`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calculate_slot_area`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.estimate_max_wires_in_slot`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.spec`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.spec.Jr`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.spec.Js`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.spec.Steel`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.spec.VoltageRating`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.spec.stator_phase_current_rms`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.specs`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`

### Variables: `self.speed_rpm`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.ss_avg_force_angle`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_avg_force_magnitude`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_avg_force_vector`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_avg_force_vector['0']`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_avg_force_vector['1']`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_max_force_err_abs`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_max_force_err_abs['0']`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_max_force_err_abs['1']`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_max_force_err_ang`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`
- **Written by (Modifies this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_max_force_err_ang['0']`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.ss_max_force_err_ang['1']`
- **Read by (Depends on this)**:
  - `utility.py::suspension_force_vector.__init__`

### Variables: `self.stack_length`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.stack_length_max`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.stator`
- **Read by (Depends on this)**:
  - `machine_core.py::MotorParameters.rInner`
  - `machine_core.py::MotorParameters.rOuter`

### Variables: `self.stator.ID`
- **Read by (Depends on this)**:
  - `machine_core.py::MotorParameters.rInner`

### Variables: `self.stator.OD`
- **Read by (Depends on this)**:
  - `machine_core.py::MotorParameters.rOuter`

### Variables: `self.stator.yoke`
- **Read by (Depends on this)**:
  - `machine_core.py::MotorParameters.rOuter`

### Variables: `self.stator_id`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `WireSlot_v2.py::MotorThermalAnalyzer.analyze_thermal_performance`
  - `WireSlot_v2.py::MotorThermalAnalyzer.calculate_slot_area`
  - `WireSlot_v2.py::MotorThermalAnalyzer.estimate_max_wires_in_slot`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calculate_slot_area`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.estimate_max_wires_in_slot`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.stator_od`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.stator_tooth_width_b_ds`
- **Read by (Depends on this)**:
  - `utility.py::Pyrhonen_design.__init__`
- **Written by (Modifies this)**:
  - `utility.py::Pyrhonen_design.__init__`

### Variables: `self.str_best_design_details`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.study`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.study_name`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `self.surface`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.convert_to_pdf`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.surface.finish`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::CairoDrawer.convert_to_pdf`

### Variables: `self.suspen_kd_at_h`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.suspen_kp_at_h`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.suspen_kw_at_h`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.sw`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.__init__`

### Variables: `self.sw.fea_config_dict`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.sw.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.sensitivity_bar_charts`

### Variables: `self.sw.fea_config_dict['use_weights']`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.BeariW_poles`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.DriveW_poles`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.Radius_OuterStatorYoke`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.template`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.template.SI`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.template.SI['GP']`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.template.SI['GP']['mm_r_ro']`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.sw.im.template.SI['GP']['mm_r_ro'].value`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.swarm_data_as_dict`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.get_metric_of_the_whole_swarm`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`

### Variables: `self.swarm_data_as_dict.items`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.get_metric_of_the_whole_swarm`

### Variables: `self.swarm_data_project_names`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`

### Variables: `self.swarm_data_raw`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.__init__`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::swarm_data_container.__init__`

### Variables: `self.swarm_data_xf`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
  - `modern_machine_designer_utility.py::swarm_data_container._extract_performance_lists`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container._initialize_empty`

### Variables: `self.swarm_data_xf.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.swarm_data_xf['0']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.__init__`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_json`
  - `modern_machine_designer_utility.py::swarm_data_container._load_from_raw`

### Variables: `self.sym_begin_pos`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
  - `winding_layout.py::PhaseWinding.setSymPos`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`

### Variables: `self.sym_begin_pos_1`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`

### Variables: `self.sym_begin_pos_2`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.setSymPos`

### Variables: `self.sym_turn_func`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.sym_winding_func`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.t`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.target`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__repr__`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`

### Variables: `self.target.machine_class`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__repr__`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.target_j`
- **Read by (Depends on this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`
- **Written by (Modifies this)**:
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.template`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.template.SI`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.template.SI['GP']`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.template.SI['GP']['mm_r_ro']`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`

### Variables: `self.template.SI['GP']['mm_r_ro'].value`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.terminal_voltage`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.time_list`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.to_dict`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.save_to_file`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_json`

### Variables: `self.to_dict_full`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.save_to_file_full`

### Variables: `self.tooth_depth`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `WireSlot_v2.py::MotorThermalAnalyzer.calculate_slot_area`
  - `WireSlot_v2.py::MotorThermalAnalyzer.estimate_max_wires_in_slot`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calculate_slot_area`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.estimate_max_wires_in_slot`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.tooth_shape`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.alpha_stator_tooth_span`
  - `machine_geometry_utils.py::MachineGeometry.d_tooth_shoe`

### Variables: `self.tooth_width`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `WireSlot_v2.py::MotorThermalAnalyzer.calculate_slot_area`
  - `WireSlot_v2.py::MotorThermalAnalyzer.estimate_max_wires_in_slot`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_motor_losses`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calculate_slot_area`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.estimate_max_wires_in_slot`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `self.torque_average`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Swarm_Data_Analyzer.prepare_data_for_post_processing`

### Variables: `self.torque_kd_at_h`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.torque_kp_at_h`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.torque_kw_at_h`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.ts`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.turn_func`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`

### Variables: `self.turn_func_bias`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.turns_per_slot`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
  - `winding_layout.py::PhaseWinding.setTurnFuncObject`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.type`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
  - `modern_machine_designer_utility.py::Parameter.__repr__`
  - `modern_machine_designer_utility.py::Parameter.to_dict`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.ui_info`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `self.unit`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
  - `modern_machine_designer_utility.py::Parameter.__repr__`
  - `modern_machine_designer_utility.py::Parameter.to_dict`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.value`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`
  - `modern_machine_designer_utility.py::Parameter.__repr__`
  - `modern_machine_designer_utility.py::Parameter.to_dict`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Parameter.__init__`

### Variables: `self.verbose`
- **Read by (Depends on this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.format_print_out_string`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_kw_per_phase`
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.get_complex_number_winding_factor_of_coil_i`
- **Written by (Modifies this)**:
  - `winding_layout_derivation_ismb2021_asymetry_no_drawing.py::Winding_Derivation.__init__`

### Variables: `self.verbose_drawing`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.drawArc`
  - `modern_machine_designer_utility.py::CairoDrawer.drawLine`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`
  - `modern_machine_designer_utility.py::CairoDrawer.__init__`

### Variables: `self.view`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.visualization_points`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.draw`
  - `modern_machine_designer_utility.py::Geometry.print_parameters`
- **Written by (Modifies this)**:
  - `modern_machine_designer_utility.py::Geometry.draw`

### Variables: `self.visualization_points.items`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Geometry.print_parameters`

### Variables: `self.weights_name`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_torque_force`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.weights_used`
- **Read by (Depends on this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`
  - `utility.py::SwarmDataAnalyzer.my_scatter_plot`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth`
  - `utility.py::SwarmDataAnalyzer.pareto_plot_torque_force`
- **Written by (Modifies this)**:
  - `utility.py::SwarmDataAnalyzer.build_basic_info`

### Variables: `self.winding`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.l_stack`
  - `machine_geometry_utils.py::MachineGeometry.pole_count`
  - `machine_geometry_utils.py::MachineGeometry.slot_count`
  - `machine_geometry_utils.py::MachineGeometry.sync`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_rotor_volume`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`
- **Written by (Modifies this)**:
  - `machine_geometry_utils.py::MachineGeometry.sync`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.__init__`

### Variables: `self.winding.EX`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_rotor_volume`
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.winding.EX.copy`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.winding.EX['mm_stack_length_specified']`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.get_rotor_volume`

### Variables: `self.winding.wily`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.winding.wily.to_dict`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.to_dict_full`

### Variables: `self.winding['l_stack']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.l_stack`

### Variables: `self.winding['pole_count']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.pole_count`

### Variables: `self.winding['slot_count']`
- **Read by (Depends on this)**:
  - `machine_geometry_utils.py::MachineGeometry.slot_count`

### Variables: `self.winding_func`
- **Read by (Depends on this)**:
  - `winding_layout.py::PhaseWinding.__init__`
- **Written by (Modifies this)**:
  - `winding_layout.py::PhaseWinding.__init__`

### Variables: `self.workDir`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.__init__`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.__init__`

### Variables: `self.yoke_thickness`
- **Read by (Depends on this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `WireSlot_v2.py::MotorThermalAnalyzer.analyze_thermal_performance`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.analyze`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.calc_bemf_constants`
- **Written by (Modifies this)**:
  - `WireSlot_v2.py::MotorThermalAnalyzer.__init__`
  - `machine_analyzer.py::MotorPerformanceAnalyzer.__init__`

### Variables: `stator.ID`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `stator.OD`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `stator.toothDepth`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `stator.toothWidth`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `stator.yoke`
- **Read by (Depends on this)**:
  - `machine_geometry.py::AllPoints.__post_init__`

### Variables: `time_list.append`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
  - `utility.py::read_csv_results_4_comparison__transient`

### Variables: `warnings.append`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Modern_Machine_Designer_Utility.validate_parameters`

### Variables: `wily`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.pre_process_PMSM`
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`
  - `modern_machine_designer_utility.py::Winding.plot_winding_function`
- **Written by (Modifies this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.add_magnetic_transient_study`
  - `JMAG.py::JMAG.pre_process_PMSM`
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `wily.CommutatingSequenceB`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.CommutatingSequenceD`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.DPNV_or_SEPA`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`

### Variables: `wily.Qs`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`
  - `modern_machine_designer_utility.py::Winding.plot_winding_function`

### Variables: `wily.bool_3PhaseCurrentSource`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `wily.bool_CustomizedCircuit`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `wily.bool_DPNVorSEPA`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.add_magnetic_transient_study`

### Variables: `wily.bool_distributed_or_concentrated`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.coil_pitch_y`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `wily.dict_coil_connection`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.dict_coil_connection['layer X phases']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.dict_coil_connection['layer X signs']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.dict_coil_connection['layer Y phases']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.dict_coil_connection['layer Y signs']`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.grouping_AC`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.layer_X_phases`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `wily.layer_X_signs`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `wily.layer_Y_phases`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `wily.layer_Y_signs`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.pre_process_PMSM`

### Variables: `wily.list_layer_motor_phases`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`

### Variables: `wily.list_layer_motor_signs`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`

### Variables: `wily.list_layer_suspension_phases`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`

### Variables: `wily.list_layer_suspension_signs`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.draw_winding_in_the_slot`

### Variables: `wily.m`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.plot_winding_function`

### Variables: `wily.number_of_parallel_branch`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`
  - `JMAG.py::JMAG.read_csv_results_4_general_purpose`

### Variables: `wily.number_of_winding_layer`
- **Read by (Depends on this)**:
  - `JMAG.py::JMAG.add_circuit`

### Variables: `wily.number_winding_layer`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.plot_winding_function`

### Variables: `wily.ox_distribution_phase_U`
- **Read by (Depends on this)**:
  - `modern_machine_designer_utility.py::Winding.plot_winding_function`

## Function Input/Output List


## File: `angle_error_nick.py`

### `angle_error` (Line 2)
- **Arguments**: alpha_star, alpha_actual
- **State Inputs**: None detected
- **State Outputs**: None detected

### `compute_angle_error` (Line 63)
- **Arguments**: alpha_star, alpha_actual
- **State Inputs**: None detected
- **State Outputs**: None detected

### `main` (Line 90)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected


## File: `JMAG.py`

### `JMAG.__init__` (Line 31)
- **Arguments**: self, fea_config_dict
- **State Inputs (Reads)**:
  - `self.JMAG_version_number`
  - `self.app`
  - `self.ass`
  - `self.bMirror`
  - `self.bool_suppressShaft`
  - `self.consts`
  - `self.defaultUnit`
  - `self.doc`
  - `self.edge4Ref`
  - `self.fea_config_dict`
  - `self.flag_material_already_loaded`
  - `self.geomApp`
  - `self.iRotateCopy`
  - `self.jd`
  - `self.model`
  - `self.projName`
  - `self.sketch`
  - `self.sketchNameList`
  - `self.study`
  - `self.verbose_drawing`
  - `self.view`
  - `self.workDir`
- **State Outputs (Writes)**:
  - `self.JMAG_version_number`
  - `self.app`
  - `self.ass`
  - `self.bMirror`
  - `self.bool_suppressShaft`
  - `self.consts`
  - `self.defaultUnit`
  - `self.doc`
  - `self.edge4Ref`
  - `self.fea_config_dict`
  - `self.flag_material_already_loaded`
  - `self.geomApp`
  - `self.iRotateCopy`
  - `self.jd`
  - `self.model`
  - `self.projName`
  - `self.sketch`
  - `self.sketchNameList`
  - `self.study`
  - `self.verbose_drawing`
  - `self.view`
  - `self.workDir`

### `JMAG.open` (Line 59)
- **Arguments**: self, Steel_name, expected_project_file_path, pc_name, dir_parent, bool_jmagDesignerShow
- **State Inputs (Reads)**:
  - `app`
  - `app.GetMaterialLibrary`
  - `app.Hide`
  - `app.NewProject`
  - `app.SaveAs`
  - `app.Show`
  - `app.VersionString`
  - `self.JMAG_version_number`
  - `self.JMAG_version_string`
  - `self.app`
  - `self.fea_config_dict`
  - `self.fea_config_dict['pc_name']`
  - `self.flag_material_already_loaded`
- **State Outputs (Writes)**:
  - `app`
  - `self.JMAG_version_number`
  - `self.JMAG_version_string`
  - `self.app`
  - `self.flag_material_already_loaded`

### `JMAG.close` (Line 269)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.app`
  - `self.app.Quit`
- **State Outputs**: None detected

### `JMAG.save` (Line 272)
- **Arguments**: self, name, description
- **State Inputs (Reads)**:
  - `self.app`
  - `self.app.GetCurrentModel`
  - `self.doc`
  - `self.doc.SaveModel`
- **State Outputs**: None detected

### `JMAG.pre_process_PMSM` (Line 280)
- **Arguments**: self, app, model, acm_variant
- **State Inputs (Reads)**:
  - `acm_variant.winding.wily`
  - `app`
  - `app.ExportImageWithSize`
  - `app.View`
  - `list_xy_magnets.append`
  - `self.doc`
  - `self.doc.GetSelection`
  - `self.fea_config_dict`
  - `self.fea_config_dict['designer.show']`
  - `self.id_rotorCore`
  - `self.id_statorCore`
  - `wily`
  - `wily.layer_X_phases`
  - `wily.layer_X_signs`
  - `wily.layer_Y_phases`
  - `wily.layer_Y_signs`
- **State Outputs (Writes)**:
  - `self.id_rotorCore`
  - `self.id_statorCore`
  - `wily`

### `JMAG.add_magnetic_transient_study` (Line 456)
- **Arguments**: self, app, model, path2FEACsv, study_name, acm_variant
- **State Inputs (Reads)**:
  - `acm_variant.winding.wily`
  - `app`
  - `app.GetDataManager`
  - `app.SetCurrentStudy`
  - `app.ShowCircuitGrid`
  - `self.JMAG_version_number`
  - `self.add_circuit`
  - `self.add_material`
  - `self.fea_config_dict`
  - `self.fea_config_dict['designer.max_nonlinear_iteration']`
  - `self.id_rotorCore`
  - `self.id_statorCore`
  - `self.study_name`
  - `wily`
  - `wily.bool_3PhaseCurrentSource`
  - `wily.bool_CustomizedCircuit`
  - `wily.bool_DPNVorSEPA`
- **State Outputs (Writes)**:
  - `self.study_name`
  - `wily`

### `JMAG.add_structural_static_study` (Line 746)
- **Arguments**: self
- **State Inputs**: None detected
- **State Outputs**: None detected

### `JMAG.add_mesh` (Line 748)
- **Arguments**: self, study, model
- **State Inputs**: None detected
- **State Outputs**: None detected

### `JMAG.add_material` (Line 751)
- **Arguments**: self, study, acm_variant
- **State Inputs**: None detected
- **State Outputs**: None detected

### `JMAG.add_circuit` (Line 795)
- **Arguments**: self, app, model, study, acm_variant, bool_3PhaseCurrentSource
- **State Inputs (Reads)**:
  - `acm_variant.winding.wily`
  - `app`
  - `app.FunctionFactory`
  - `app.ShowCircuitGrid`
  - `wily`
  - `wily.CommutatingSequenceB`
  - `wily.CommutatingSequenceD`
  - `wily.bool_DPNVorSEPA`
  - `wily.bool_distributed_or_concentrated`
  - `wily.coil_pitch_y`
  - `wily.dict_coil_connection`
  - `wily.dict_coil_connection['layer X phases']`
  - `wily.dict_coil_connection['layer X signs']`
  - `wily.dict_coil_connection['layer Y phases']`
  - `wily.dict_coil_connection['layer Y signs']`
  - `wily.grouping_AC`
  - `wily.layer_X_phases`
  - `wily.layer_X_signs`
  - `wily.layer_Y_phases`
  - `wily.layer_Y_signs`
  - `wily.number_of_parallel_branch`
  - `wily.number_of_winding_layer`
- **State Outputs (Writes)**:
  - `wily`

### `JMAG.addConstraintCocentricity` (Line 1040)
- **Arguments**: self, vA, vB
- **State Inputs (Reads)**:
  - `self.doc`
  - `self.doc.CreateReferenceFromItem`
  - `self.sketch`
  - `self.sketch.CreateBiConstraint`
  - `self.sketch.GetItem`
- **State Outputs**: None detected

### `JMAG.drawLine` (Line 1060)
- **Arguments**: self, startxy, endxy, returnVertexName
- **State Inputs (Reads)**:
  - `self.getSketch`
  - `self.sketch`
  - `self.sketch.CreateLine`
  - `self.sketch.OpenSketch`
- **State Outputs (Writes)**:
  - `self.sketch`

### `JMAG.drawArc` (Line 1078)
- **Arguments**: self, centerxy, startxy, endxy, returnVertexName
- **State Inputs (Reads)**:
  - `self.getSketch`
  - `self.sketch`
  - `self.sketch.CreateArc`
  - `self.sketch.OpenSketch`
- **State Outputs (Writes)**:
  - `self.sketch`

### `JMAG.drawCircle` (Line 1095)
- **Arguments**: self, centerxy, radius, returnVertexName
- **State Inputs (Reads)**:
  - `self.getSketch`
  - `self.sketch`
  - `self.sketch.CreateCircle`
  - `self.sketch.OpenSketch`
- **State Outputs (Writes)**:
  - `self.sketch`

### `JMAG.checkGeomApp` (Line 1109)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.app`
  - `self.app.CreateGeometryEditor`
  - `self.app.LaunchGeometryEditor`
  - `self.doc`
  - `self.geomApp`
  - `self.geomApp.NewDocument`
- **State Outputs (Writes)**:
  - `self.doc`
  - `self.geomApp`

### `JMAG.getSketch` (Line 1117)
- **Arguments**: self, sketchName, color
- **State Inputs (Reads)**:
  - `self.ass`
  - `self.ass.CreateSketch`
  - `self.ass.GetItem`
  - `self.checkGeomApp`
  - `self.doc`
  - `self.doc.CreateReferenceFromItem`
  - `self.doc.GetAssembly`
  - `self.geomApp`
  - `self.geomApp.GetDocument`
  - `self.sketch`
  - `self.sketch.OpenSketch`
  - `self.sketch.SetProperty`
  - `self.sketchNameList`
  - `self.sketchNameList.append`
- **State Outputs (Writes)**:
  - `self.ass`
  - `self.doc`
  - `self.geomApp`
  - `self.sketch`

### `JMAG.prepareSection` (Line 1142)
- **Arguments**: self, token, bMirrorMerge, bRotateMerge
- **State Inputs (Reads)**:
  - `list_region_objects.append`
  - `self.bMirror`
  - `self.doc`
  - `self.doc.GetSelection`
  - `self.edge4Ref`
  - `self.edge4ref`
  - `self.iRotateCopy`
  - `self.regionCircularPattern360Origin`
  - `self.regionMirrorCopy`
  - `self.sketch`
  - `self.sketch.CloseSketch`
  - `self.sketch.CreateRegions`
  - `self.sketch.GetItem`
- **State Outputs**: None detected

### `JMAG.regionMirrorCopy` (Line 1197)
- **Arguments**: self, region, edge4Ref, symmetryType, bMerge
- **State Inputs (Reads)**:
  - `self.ass`
  - `self.ass.GetItem`
  - `self.doc`
  - `self.doc.CreateReferenceFromItem`
  - `self.sketch`
  - `self.sketch.CreateRegionMirrorCopy`
  - `self.sketch.GetItem`
- **State Outputs**: None detected

### `JMAG.regionCircularPattern360Origin` (Line 1219)
- **Arguments**: self, idx, region, Q_float, bMerge
- **State Inputs (Reads)**:
  - `self.doc`
  - `self.doc.CreateReferenceFromItem`
  - `self.sketch`
  - `self.sketch.CreateRegionCircularPattern`
- **State Outputs**: None detected

### `JMAG.draw_jmag_model` (Line 1265)
- **Arguments**: self, app, individual_index, im_variant, model_name, bool_trimDrawer_or_vanGogh, doNotRotateCopy
- **State Inputs (Reads)**:
  - `app`
  - `app.CreateGeometryEditor`
  - `app.GetCurrentModel`
  - `app.LaunchGeometryEditor`
  - `app.NumModels`
  - `self.SI`
- **State Outputs (Writes)**:
  - `self.SI`

### `JMAG.run_study` (Line 1335)
- **Arguments**: acm_variant, app, study, fea_config_dict, toc
- **State Inputs (Reads)**:
  - `app`
  - `app.Save`
- **State Outputs**: None detected

### `JMAG.mesh_study` (Line 1360)
- **Arguments**: self, acm_variant, app, model, study, output_dir
- **State Inputs (Reads)**:
  - `app`
  - `app.ExportImageWithSize`
  - `app.View`
- **State Outputs**: None detected

### `JMAG.draw_spmsm` (Line 1444)
- **Arguments**: self, acm_variant, bool_pyx
- **State Inputs (Reads)**:
  - `acm_variant.rotorMagnet.notched_rotor.p`
  - `self.bMirror`
  - `self.calculate_excitation_current`
  - `self.iRotateCopy`
  - `self.prepareSection`
  - `self.save`
  - `self.show`
- **State Outputs (Writes)**:
  - `self.bMirror`
  - `self.iRotateCopy`

### `JMAG.show` (Line 1585)
- **Arguments**: self, acm_variant, toString
- **State Inputs**: None detected
- **State Outputs**: None detected

### `JMAG.add_plots` (Line 1610)
- **Arguments**: axeses, dm, title, label, zorder, time_list, sfv, torque, range_ss, alpha
- **State Inputs**: None detected
- **State Outputs**: None detected

### `JMAG.read_csv_results_4_general_purpose` (Line 1669)
- **Arguments**: study_name, path_prefix, fea_config_dict, femm_solver, acm_variant
- **State Inputs (Reads)**:
  - `DisplacementAngle_list.append`
  - `ForConX_list.append`
  - `ForConY_list.append`
  - `TorCon_list.append`
  - `acm_variant.wily`
  - `basic_info.append`
  - `key_list.append`
  - `new_key_list.append`
  - `rotor_Joule_loss_list.append`
  - `self.Current_dict`
  - `self.Current_dict['Time(s)']`
  - `self.DisplacementAngle_list`
  - `self.FluxLinkage_dict`
  - `self.ForConAbs_list`
  - `self.ForConX_list`
  - `self.ForConY_list`
  - `self.TorCon_list`
  - `self.basic_info`
  - `self.circuit_current`
  - `self.coil_fluxLinkage`
  - `self.femm_loss_list`
  - `self.get_voltage_and_current`
  - `self.jmag_loss_list`
  - `self.mycurrent`
  - `self.mytime`
  - `self.myvoltage`
  - `self.terminal_voltage`
  - `self.time_list`
  - `self.ui_info`
  - `time_list.append`
  - `wily`
  - `wily.coil_pitch_y`
  - `wily.number_of_parallel_branch`
- **State Outputs (Writes)**:
  - `self.ForConAbs_list`
  - `self.ForConX_list`
  - `self.ForConY_list`
  - `self.TorCon_list`
  - `self.basic_info`
  - `self.femm_loss_list`
  - `self.jmag_loss_list`
  - `self.mycurrent`
  - `self.mytime`
  - `self.myvoltage`
  - `self.time_list`
  - `self.ui_info`
  - `wily`

### `JMAG.build_str_results` (Line 2097)
- **Arguments**: self, acm_variant, project_name, tran_study_name, path2FEACsv, fea_config_dict, femm_solver
- **State Inputs (Reads)**:
  - `coil_flux_linkage_peak2peak_value_results.append`
  - `fitness_mapping`
  - `self.add_plots`
  - `self.dm`
  - `self.fea_config_dict`
  - `self.fea_config_dict['designer.StepPerCycle_3rdTSS']`
  - `self.read_csv_results_4_general_purpose`
- **State Outputs (Writes)**:
  - `fitness_mapping`
  - `self.dm`


## File: `machine.py`

### `Machine.__init__` (Line 6)
- **Arguments**: self, motor
- **State Inputs (Reads)**:
  - `self.geometry`
  - `self.motor`
- **State Outputs (Writes)**:
  - `self.geometry`
  - `self.motor`

### `Machine.sync` (Line 10)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.geometry`
  - `self.geometry.add_part`
  - `self.geometry.sync`
  - `self.motor`
  - `self.motor.stator`
  - `self.motor.stator.tooth_shape`
- **State Outputs**: None detected

### `Machine.draw_svg` (Line 22)
- **Arguments**: self, filename
- **State Inputs (Reads)**:
  - `self.geometry`
  - `self.geometry.show_geometry_svg`
- **State Outputs**: None detected

### `Machine.draw_jmag` (Line 25)
- **Arguments**: self, filename
- **State Inputs (Reads)**:
  - `self.geometry`
  - `self.geometry.all_points`
  - `self.geometry.parts`
- **State Outputs**: None detected


## File: `machine_analyzer.py`

### `MotorPerformanceAnalyzer.__init__` (Line 16)
- **Arguments**: self, motor
- **State Inputs (Reads)**:
  - `motor.rotor.OD`
  - `motor.rotor.airGap`
  - `motor.stator.ID`
  - `motor.stator.OD`
  - `motor.stator.liner`
  - `motor.stator.toothDepth`
  - `motor.stator.toothWidth`
  - `motor.stator.yoke`
  - `self.La`
  - `self.air_gap`
  - `self.awg`
  - `self.liner`
  - `self.motor`
  - `self.poles`
  - `self.rated_speed`
  - `self.rotor_od`
  - `self.slots`
  - `self.stator_id`
  - `self.stator_od`
  - `self.target_j`
  - `self.tooth_depth`
  - `self.tooth_width`
  - `self.yoke_thickness`
- **State Outputs (Writes)**:
  - `self.La`
  - `self.air_gap`
  - `self.awg`
  - `self.liner`
  - `self.motor`
  - `self.poles`
  - `self.rated_speed`
  - `self.rotor_od`
  - `self.slots`
  - `self.stator_id`
  - `self.stator_od`
  - `self.target_j`
  - `self.tooth_depth`
  - `self.tooth_width`
  - `self.yoke_thickness`

### `MotorPerformanceAnalyzer.get_wire_properties` (Line 35)
- **Arguments**: self, awg_size
- **State Inputs**: None detected
- **State Outputs**: None detected

### `MotorPerformanceAnalyzer.calculate_slot_area` (Line 42)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.slots`
  - `self.stator_id`
  - `self.tooth_depth`
  - `self.tooth_width`
- **State Outputs**: None detected

### `MotorPerformanceAnalyzer.estimate_max_wires_in_slot` (Line 52)
- **Arguments**: self, d_od
- **State Inputs (Reads)**:
  - `self.liner`
  - `self.slots`
  - `self.stator_id`
  - `self.tooth_depth`
  - `self.tooth_width`
- **State Outputs**: None detected

### `MotorPerformanceAnalyzer.calc_bemf_constants` (Line 89)
- **Arguments**: self, z_slot, B_sat
- **State Inputs (Reads)**:
  - `self.La`
  - `self.poles`
  - `self.slots`
  - `self.yoke_thickness`
- **State Outputs**: None detected

### `MotorPerformanceAnalyzer.calc_motor_losses` (Line 113)
- **Arguments**: self, turns_per_phase, current, rpm
- **State Inputs (Reads)**:
  - `self.La`
  - `self.awg`
  - `self.get_wire_properties`
  - `self.poles`
  - `self.tooth_width`
- **State Outputs**: None detected

### `MotorPerformanceAnalyzer.analyze` (Line 141)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.awg`
  - `self.calc_bemf_constants`
  - `self.calc_motor_losses`
  - `self.calculate_slot_area`
  - `self.estimate_max_wires_in_slot`
  - `self.get_wire_properties`
  - `self.rated_speed`
  - `self.slots`
  - `self.stator_id`
  - `self.target_j`
  - `self.yoke_thickness`
- **State Outputs**: None detected


## File: `machine_core.py`

### `MotorParameters.rOuter` (Line 37)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.stator`
  - `self.stator.OD`
  - `self.stator.yoke`
- **State Outputs**: None detected

### `MotorParameters.rInner` (Line 41)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.stator`
  - `self.stator.ID`
- **State Outputs**: None detected


## File: `machine_geometry.py`

### `rotate_point` (Line 11)
- **Arguments**: p, deg
- **State Inputs**: None detected
- **State Outputs**: None detected

### `AllPoints.__post_init__` (Line 23)
- **Arguments**: self, motor
- **State Inputs (Reads)**:
  - `rotor.ID`
  - `rotor.OD`
  - `rotor.magnetDepth`
  - `self.HP`
  - `self.HP['4']`
  - `self.HP['4']['0']`
  - `self.HP['4']['1']`
  - `self.HP['6']`
  - `self.HP['6']['0']`
  - `self.HP['6']['1']`
  - `self.HP['7']`
  - `self.HP['7']['0']`
  - `self.HP['7']['1']`
  - `self.HP['9']`
  - `self.HP['9']['0']`
  - `self.HP['9']['1']`
  - `self.HP_mirror`
  - `self.RP`
  - `self.horizontal_position`
  - `self.pole_count`
  - `self.rotated_position`
  - `self.slot_count`
  - `stator.ID`
  - `stator.OD`
  - `stator.toothDepth`
  - `stator.toothWidth`
  - `stator.yoke`
- **State Outputs (Writes)**:
  - `self.HP`
  - `self.HP_mirror`
  - `self.RP`
  - `self.horizontal_position`
  - `self.pole_count`
  - `self.rotated_position`
  - `self.slot_count`

### `parse_point_name` (Line 87)
- **Arguments**: name, all_points, rotation_deg
- **State Inputs**: None detected
- **State Outputs**: None detected

### `draw_instruction_parser` (Line 119)
- **Arguments**: part, all_points, drawer
- **State Inputs (Reads)**:
  - `list_regions.append`
- **State Outputs**: None detected

### `StatorCore.draw_instruction` (Line 211)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.slot_count`
- **State Outputs**: None detected

### `RotorCore.draw_instruction` (Line 230)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.pole_count`
  - `self.options`
- **State Outputs**: None detected

### `Magnet.draw_instruction` (Line 248)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.pole_count`
- **State Outputs**: None detected

### `Coil.draw_instruction` (Line 268)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.slot_count`
- **State Outputs**: None detected

### `MachineGeometry.add_part` (Line 298)
- **Arguments**: self, part
- **State Inputs (Reads)**:
  - `self._next_index`
  - `self.all_points`
  - `self.parts`
  - `self.parts.append`
- **State Outputs**: None detected

### `MachineGeometry.sync` (Line 305)
- **Arguments**: self, motor
- **State Inputs (Reads)**:
  - `self.all_points`
- **State Outputs (Writes)**:
  - `self.all_points`

### `MachineGeometry.draw_machine_using_CairoDrawer` (Line 308)
- **Arguments**: self, drawer
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.parts`
- **State Outputs**: None detected

### `MachineGeometry.show_geometry_svg` (Line 324)
- **Arguments**: self, filename, scale
- **State Inputs (Reads)**:
  - `self.draw_machine_using_CairoDrawer`
- **State Outputs**: None detected


## File: `machine_geometry_utils.py`

### `rotate_point` (Line 9)
- **Arguments**: p, deg
- **State Inputs**: None detected
- **State Outputs**: None detected

### `AllPoints.__post_init__` (Line 21)
- **Arguments**: self, user_input
- **State Inputs (Reads)**:
  - `self.GP`
  - `self.HP`
  - `self.HP['1']`
  - `self.HP['1']['0']`
  - `self.HP['1']['1']`
  - `self.HP['2']`
  - `self.HP['2']['0']`
  - `self.HP['2']['1']`
  - `self.HP['3']`
  - `self.HP['3']['0']`
  - `self.HP['3']['1']`
  - `self.HP['4']`
  - `self.HP['4']['0']`
  - `self.HP['4']['1']`
  - `self.HP['6']`
  - `self.HP['6']['0']`
  - `self.HP['6']['1']`
  - `self.HP['7']`
  - `self.HP['7']['0']`
  - `self.HP['7']['1']`
  - `self.HP['9']`
  - `self.HP['9']['0']`
  - `self.HP['9']['1']`
  - `self.HP_mirror`
  - `self.RP`
  - `self.horizontal_position`
  - `self.pole_count`
  - `self.rotated_position`
  - `self.slot_count`
- **State Outputs (Writes)**:
  - `self.GP`
  - `self.HP`
  - `self.HP_mirror`
  - `self.RP`
  - `self.horizontal_position`
  - `self.pole_count`
  - `self.rotated_position`
  - `self.slot_count`

### `parse_point_name` (Line 103)
- **Arguments**: name, all_points, rotation_deg
- **State Inputs**: None detected
- **State Outputs**: None detected

### `draw_instruction_parser` (Line 138)
- **Arguments**: part, all_points, drawer
- **State Inputs (Reads)**:
  - `list_regions.append`
- **State Outputs**: None detected

### `RotorCore.draw_instruction` (Line 251)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.pole_count`
  - `self.options`
- **State Outputs**: None detected

### `StatorCore.draw_instruction` (Line 269)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.slot_count`
  - `self.options`
- **State Outputs**: None detected

### `Magnet.draw_instruction` (Line 305)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.pole_count`
- **State Outputs**: None detected

### `Coil.draw_instruction` (Line 325)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.all_points`
  - `self.all_points.slot_count`
- **State Outputs**: None detected

### `MachineGeometry.slot_count` (Line 360)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.winding`
  - `self.winding['slot_count']`
- **State Outputs**: None detected

### `MachineGeometry.pole_count` (Line 364)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.winding`
  - `self.winding['pole_count']`
- **State Outputs**: None detected

### `MachineGeometry.l_stack` (Line 368)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.winding`
  - `self.winding['l_stack']`
- **State Outputs**: None detected

### `MachineGeometry.d_air_gap` (Line 372)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['d_air_gap']`
- **State Outputs**: None detected

### `MachineGeometry.r_stator_outer` (Line 376)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['r_stator_outer']`
- **State Outputs**: None detected

### `MachineGeometry.r_rotor_outer` (Line 380)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['r_rotor_outer']`
- **State Outputs**: None detected

### `MachineGeometry.r_shaft` (Line 384)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['r_shaft']`
- **State Outputs**: None detected

### `MachineGeometry.w_tooth` (Line 388)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['w_tooth']`
- **State Outputs**: None detected

### `MachineGeometry.d_tooth` (Line 392)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['d_tooth']`
- **State Outputs**: None detected

### `MachineGeometry.d_magnet` (Line 396)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['d_magnet']`
- **State Outputs**: None detected

### `MachineGeometry.tooth_shape` (Line 400)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['tooth_shape']`
- **State Outputs**: None detected

### `MachineGeometry.deg_alpha_rm` (Line 404)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.pole_count`
- **State Outputs**: None detected

### `MachineGeometry.split_ratio` (Line 409)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.d_air_gap`
  - `self.r_rotor_outer`
  - `self.r_stator_outer`
- **State Outputs**: None detected

### `MachineGeometry.d_stator_yoke` (Line 413)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.d_tooth`
  - `self.r_rotor_outer`
  - `self.r_stator_outer`
- **State Outputs**: None detected

### `MachineGeometry.d_tooth_shoe` (Line 417)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.gp['d_tooth_shoe']`
  - `self.tooth_shape`
- **State Outputs**: None detected

### `MachineGeometry.alpha_stator_tooth_span` (Line 423)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.slot_count`
  - `self.tooth_shape`
- **State Outputs**: None detected

### `MachineGeometry.aspect_ratio` (Line 429)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.l_stack`
  - `self.r_stator_outer`
- **State Outputs**: None detected

### `MachineGeometry.__post_init__` (Line 432)
- **Arguments**: self
- **State Inputs**: None detected
- **State Outputs**: None detected

### `MachineGeometry.add_part` (Line 436)
- **Arguments**: self, part
- **State Inputs (Reads)**:
  - `self._next_index`
  - `self.all_points`
  - `self.parts`
  - `self.parts.append`
- **State Outputs**: None detected

### `MachineGeometry.add_radial_array` (Line 443)
- **Arguments**: self, part_type, base_name, count, options, color
- **State Inputs (Reads)**:
  - `created_parts.append`
  - `self.add_part`
- **State Outputs**: None detected

### `MachineGeometry.show_geometry_svg` (Line 453)
- **Arguments**: self, filename, scale
- **State Inputs**: None detected
- **State Outputs**: None detected

### `MachineGeometry.sync` (Line 459)
- **Arguments**: self, user_input
- **State Inputs (Reads)**:
  - `self.gp`
  - `self.winding`
- **State Outputs (Writes)**:
  - `self.gp`
  - `self.winding`


## File: `main.py`

### `get_motor_parameters` (Line 26)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected

### `update_motor_parameters` (Line 30)
- **Arguments**: params
- **State Inputs**: None detected
- **State Outputs**: None detected


## File: `modern_machine_designer_utility.py`

### `Parameter.__init__` (Line 6)
- **Arguments**: self, name, type, value, bounds, calc, calc_bounds, unit, parameter_dict
- **State Inputs (Reads)**:
  - `self.bounds`
  - `self.bounds['0']`
  - `self.bounds['1']`
  - `self.calc`
  - `self.calc_bounds`
  - `self.initialized`
  - `self.name`
  - `self.overwritten`
  - `self.parameter_dict`
  - `self.type`
  - `self.unit`
  - `self.value`
- **State Outputs (Writes)**:
  - `self.bounds`
  - `self.calc`
  - `self.calc_bounds`
  - `self.initialized`
  - `self.name`
  - `self.parameter_dict`
  - `self.type`
  - `self.unit`
  - `self.value`

### `Parameter.__repr__` (Line 43)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.name`
  - `self.type`
  - `self.unit`
  - `self.value`
- **State Outputs**: None detected

### `Parameter.sensitivity` (Line 46)
- **Arguments**: self, param_name
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Parameter.to_dict` (Line 49)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.bounds`
  - `self.name`
  - `self.type`
  - `self.unit`
  - `self.value`
- **State Outputs**: None detected

### `Parameter.from_dict` (Line 68)
- **Arguments**: cls, data
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Winding.__init__` (Line 92)
- **Arguments**: self, phase_number_m, stator_slot_number_Qs, pole_pair_number_p, suspension_pole_pair_number_ps, coil_pitch_y, bool_DPNVorSEPA, number_of_parallel_branch
- **State Inputs (Reads)**:
  - `self.CommutatingSequenceB`
  - `self.CommutatingSequenceD`
  - `self.Qs`
  - `self.SPP`
  - `self.bool_3PhaseCurrentSource`
  - `self.bool_CustomizedCircuit`
  - `self.bool_DPNVorSEPA`
  - `self.bool_distributed_or_concentrated`
  - `self.coil_pitch_y`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.dict_coil_connection`
  - `self.get_winding_factor`
  - `self.grouping_AC`
  - `self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y`
  - `self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y`
  - `self.kw1`
  - `self.layer_X_phases`
  - `self.layer_X_signs`
  - `self.layer_Y_phases`
  - `self.layer_Y_signs`
  - `self.m`
  - `self.number_of_parallel_branch`
  - `self.number_of_winding_layer`
  - `self.p`
  - `self.ps`
- **State Outputs (Writes)**:
  - `self.CommutatingSequenceB`
  - `self.CommutatingSequenceD`
  - `self.Qs`
  - `self.SPP`
  - `self.bool_3PhaseCurrentSource`
  - `self.bool_CustomizedCircuit`
  - `self.bool_DPNVorSEPA`
  - `self.coil_pitch_y`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.dict_coil_connection`
  - `self.grouping_AC`
  - `self.kw1`
  - `self.layer_X_phases`
  - `self.layer_X_signs`
  - `self.layer_Y_phases`
  - `self.layer_Y_signs`
  - `self.m`
  - `self.number_of_parallel_branch`
  - `self.number_of_winding_layer`
  - `self.p`
  - `self.ps`

### `Winding.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y` (Line 196)
- **Arguments**: self, layer_X_phases, coil_pitch
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Winding.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y` (Line 198)
- **Arguments**: self, layer_X_signs, coil_pitch
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Winding.get_winding_factor` (Line 202)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.Qs`
  - `self.coil_pitch_y`
  - `self.m`
  - `self.p`
  - `self.ps`
- **State Outputs**: None detected

### `Winding.draw_winding_in_the_slot` (Line 208)
- **Arguments**: u, Qs, list_layer_phases, list_layer_signs, text
- **State Inputs (Reads)**:
  - `self.path2SwarmData`
  - `wily`
  - `wily.DPNV_or_SEPA`
  - `wily.Qs`
  - `wily.list_layer_motor_phases`
  - `wily.list_layer_motor_signs`
  - `wily.list_layer_suspension_phases`
  - `wily.list_layer_suspension_signs`
- **State Outputs**: None detected

### `Winding.plot_winding_function` (Line 251)
- **Arguments**: wily
- **State Inputs (Reads)**:
  - `self.path2SwarmData`
  - `wily`
  - `wily.Qs`
  - `wily.m`
  - `wily.number_winding_layer`
  - `wily.ox_distribution_phase_U`
- **State Outputs**: None detected

### `Winding.to_dict` (Line 263)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.Qs`
  - `self.coil_pitch_y`
  - `self.derivation`
  - `self.derivation.__dict__`
  - `self.derivation.__dict__.items`
  - `self.kw1`
  - `self.m`
  - `self.number_of_parallel_branch`
  - `self.p`
  - `self.ps`
- **State Outputs**: None detected

### `Winding.from_dict` (Line 317)
- **Arguments**: cls, data
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Geometry.__init__` (Line 362)
- **Arguments**: self, name, GP, draw_function, color
- **State Inputs (Reads)**:
  - `self.GP`
  - `self.GP.items`
  - `self.color`
  - `self.draw_function`
  - `self.name`
- **State Outputs (Writes)**:
  - `self.GP`
  - `self.color`
  - `self.draw_function`
  - `self.name`

### `Geometry.update_from_GP` (Line 371)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.GP`
  - `self.GP.items`
- **State Outputs**: None detected

### `Geometry.__repr__` (Line 382)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.GP`
  - `self.GP.items`
  - `self.color`
  - `self.name`
- **State Outputs**: None detected

### `Geometry.print_parameters` (Line 393)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.GP`
  - `self.GP.items`
  - `self.__dict__`
  - `self.__dict__.items`
  - `self.color`
  - `self.name`
  - `self.visualization_points`
  - `self.visualization_points.items`
- **State Outputs**: None detected

### `Geometry.to_dict` (Line 439)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.__dict__`
  - `self.__dict__.items`
- **State Outputs**: None detected

### `Geometry.from_dict` (Line 484)
- **Arguments**: cls, data
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Geometry.draw` (Line 490)
- **Arguments**: self, drawer
- **State Inputs (Reads)**:
  - `self.components_make_region`
  - `self.draw_function`
  - `self.visualization_points`
- **State Outputs (Writes)**:
  - `self.components_make_region`
  - `self.visualization_points`

### `CairoDrawer.__init__` (Line 514)
- **Arguments**: self, width_in_points, height_in_points, filename, verbose_drawing, scale, bFillRegion
- **State Inputs (Reads)**:
  - `self.bFillRegion`
  - `self.bMirror`
  - `self.ctx`
  - `self.ctx.paint`
  - `self.ctx.restore`
  - `self.ctx.save`
  - `self.ctx.scale`
  - `self.ctx.set_source_rgb`
  - `self.ctx.transform`
  - `self.filename`
  - `self.iRotateCopy`
  - `self.scale`
  - `self.surface`
  - `self.verbose_drawing`
- **State Outputs (Writes)**:
  - `self.bFillRegion`
  - `self.bMirror`
  - `self.ctx`
  - `self.filename`
  - `self.iRotateCopy`
  - `self.scale`
  - `self.surface`
  - `self.verbose_drawing`

### `CairoDrawer.apply_stroke` (Line 537)
- **Arguments**: self, lw
- **State Inputs (Reads)**:
  - `self.ctx`
  - `self.ctx.set_line_cap`
  - `self.ctx.set_line_width`
  - `self.ctx.set_source_rgba`
  - `self.ctx.stroke`
- **State Outputs**: None detected

### `CairoDrawer.convert_to_pdf` (Line 544)
- **Arguments**: self, bool_open_pdf, filename
- **State Inputs (Reads)**:
  - `self.sketch_color`
  - `self.surface`
  - `self.surface.finish`
- **State Outputs (Writes)**:
  - `self.sketch_color`

### `CairoDrawer.hex_to_rgb` (Line 561)
- **Arguments**: self, hex_color
- **State Inputs**: None detected
- **State Outputs**: None detected

### `CairoDrawer.drawLine` (Line 577)
- **Arguments**: self, p1, p2
- **State Inputs (Reads)**:
  - `self.verbose_drawing`
- **State Outputs**: None detected

### `CairoDrawer.drawArc` (Line 582)
- **Arguments**: self, centerxy, startxy, endxy
- **State Inputs (Reads)**:
  - `self.verbose_drawing`
- **State Outputs**: None detected

### `CairoDrawer.getSketch` (Line 587)
- **Arguments**: self, name, color
- **State Inputs**: None detected
- **State Outputs**: None detected

### `CairoDrawer.prepareSection` (Line 591)
- **Arguments**: self, region_dict, color
- **State Inputs (Reads)**:
  - `self.bFillRegion`
  - `self.ctx`
  - `self.ctx.arc`
  - `self.ctx.arc_negative`
  - `self.ctx.fill_preserve`
  - `self.ctx.line_to`
  - `self.ctx.move_to`
  - `self.ctx.new_path`
  - `self.ctx.path_extents`
  - `self.ctx.restore`
  - `self.ctx.rotate`
  - `self.ctx.save`
  - `self.ctx.scale`
  - `self.ctx.set_line_width`
  - `self.ctx.set_source_rgba`
  - `self.ctx.stroke`
  - `self.hex_to_rgb`
  - `self.scale`
- **State Outputs**: None detected

### `CairoDrawer.finalize_part` (Line 675)
- **Arguments**: self
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.__init__` (Line 680)
- **Arguments**: self, specs
- **State Inputs (Reads)**:
  - `self.flag_do_not_evaluate_when_init_pop`
  - `self.geometry`
  - `self.specs`
  - `self.target`
  - `self.winding`
- **State Outputs (Writes)**:
  - `self.flag_do_not_evaluate_when_init_pop`
  - `self.geometry`
  - `self.specs`
  - `self.target`
  - `self.winding`

### `Modern_Machine_Designer_Utility._get_parameter_logger` (Line 690)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.name`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.apply_parameter_dict` (Line 711)
- **Arguments**: self, prev_params, key_map
- **State Inputs (Reads)**:
  - `self._get_parameter_logger`
  - `self.geometry`
  - `self.geometry.machineGeometry`
  - `self.geometry.machineGeometry.values`
  - `self.get_parameter_dict_by_name`
  - `self.get_parameters_by_type`
  - `self.parameter_dict_by_name`
  - `self.parameter_dict_by_name.get`
- **State Outputs (Writes)**:
  - `self.parameter_dict_by_name`

### `Modern_Machine_Designer_Utility.update_geometric_parameters` (Line 781)
- **Arguments**: self, x_denorm, x_denorm_dict
- **State Inputs (Reads)**:
  - `self.geometry`
  - `self.geometry.machineGeometry`
  - `self.geometry.machineGeometry.items`
  - `self.get_free_variables`
  - `self.get_parameter_dict_by_name`
  - `self.get_parameters_by_type`
  - `self.parameter_dict_by_name`
- **State Outputs (Writes)**:
  - `self.parameter_dict_by_name`

### `Modern_Machine_Designer_Utility.get_pc_name` (Line 809)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.learn_about_the_archive` (Line 822)
- **Arguments**: self, prob, swarm_data, popsize, bool_plot_and_show, bool_more_info
- **State Inputs (Reads)**:
  - `more_info.append`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.write_swarm_survivor` (Line 901)
- **Arguments**: self, pop, counter_fitness_return
- **State Inputs (Reads)**:
  - `self.path2SwarmData`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_bad_fintess_values` (Line 917)
- **Arguments**: self, machine_type, ref
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_rotor_volume` (Line 941)
- **Arguments**: self, stack_length
- **State Inputs (Reads)**:
  - `self.mm_r_ro`
  - `self.mm_r_ro.value`
  - `self.winding`
  - `self.winding.EX`
  - `self.winding.EX['mm_stack_length_specified']`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_rotor_weight` (Line 949)
- **Arguments**: self, gravity, stack_length
- **State Inputs (Reads)**:
  - `self.get_rotor_volume`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_free_variables` (Line 966)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_free_variables_as_dict` (Line 969)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_free_variables`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_free_variable_bounds_dict` (Line 980)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_free_variables`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.set_free_variables_from_dict` (Line 992)
- **Arguments**: self, free_variables_dict
- **State Inputs (Reads)**:
  - `self.get_parameter`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.update_derived_parameters` (Line 999)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_parameter_fields` (Line 1017)
- **Arguments**: self
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_parameters_by_type` (Line 1032)
- **Arguments**: self, param_type
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_parameter` (Line 1045)
- **Arguments**: self, name
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.set_parameter_value` (Line 1058)
- **Arguments**: self, name, value
- **State Inputs (Reads)**:
  - `self.get_parameter`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_parameter_dict` (Line 1075)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_parameter_dict_by_name` (Line 1084)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.list_parameters` (Line 1093)
- **Arguments**: self, param_type
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
  - `self.get_parameters_by_type`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.get_parameters_summary` (Line 1107)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_parameter_fields`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.validate_parameters` (Line 1136)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `errors.append`
  - `self.get_parameter_fields`
  - `warnings.append`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.__repr__` (Line 1168)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.get_parameters_summary`
  - `self.target`
  - `self.target.machine_class`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.to_dict` (Line 1177)
- **Arguments**: self
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.to_json` (Line 1219)
- **Arguments**: self, indent, ensure_ascii
- **State Inputs (Reads)**:
  - `self.to_dict`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.save_to_file` (Line 1232)
- **Arguments**: self, filepath, indent
- **State Inputs (Reads)**:
  - `self.to_dict`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.to_dict_full` (Line 1243)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `full_dict['wily']`
  - `self.__class__`
  - `self.__class__.__module__`
  - `self.__class__.__name__`
  - `self.bool_PermanentMagnet`
  - `self.bool_RotorNotched`
  - `self.bool_StatorSlotClosed`
  - `self.bool_jmagDeleteResultsAfterCalculation`
  - `self.counter`
  - `self.geometry`
  - `self.geometry.machineGeometry`
  - `self.geometry.machineGeometry.items`
  - `self.get_parameter_fields`
  - `self.name`
  - `self.select_FEA_tool`
  - `self.select_fea_config_dict`
  - `self.target`
  - `self.target.machine_class`
  - `self.winding`
  - `self.winding.EX`
  - `self.winding.EX.copy`
  - `self.winding.wily`
  - `self.winding.wily.to_dict`
- **State Outputs (Writes)**:
  - `full_dict['wily']`

### `Modern_Machine_Designer_Utility.save_to_file_full` (Line 1399)
- **Arguments**: self, filepath, indent
- **State Inputs (Reads)**:
  - `self.to_dict_full`
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.from_dict_full` (Line 1411)
- **Arguments**: cls, data
- **State Inputs (Reads)**:
  - `data['wily']`
  - `instance.winding.wily`
- **State Outputs (Writes)**:
  - `instance.winding.wily`

### `Modern_Machine_Designer_Utility.load_from_file_full` (Line 1552)
- **Arguments**: cls, filepath
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.from_dict` (Line 1567)
- **Arguments**: cls, data
- **State Inputs (Reads)**:
  - `instance.winding.wily`
- **State Outputs (Writes)**:
  - `instance.winding.wily`

### `Modern_Machine_Designer_Utility.from_json` (Line 1639)
- **Arguments**: cls, json_str
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.load_from_file` (Line 1653)
- **Arguments**: cls, filepath
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.remove_jfiles_folders` (Line 1669)
- **Arguments**: root_dir
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.myLogger` (Line 1691)
- **Arguments**: dir_log, prefix
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Modern_Machine_Designer_Utility.draw_individual_from_swarm` (Line 1720)
- **Arguments**: self, index
- **State Inputs (Reads)**:
  - `self.path2SwarmData`
  - `self.show_geometry`
- **State Outputs**: None detected

### `Swarm_Data_Analyzer.decode_py_reduce_ordered_dict` (Line 1764)
- **Arguments**: x_denorm_dict_raw
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Swarm_Data_Analyzer.__init__` (Line 1807)
- **Arguments**: self, fname, desired_x_denorm_dict, bool_filter_pareto_front
- **State Inputs (Reads)**:
  - `self.decode_py_reduce_ordered_dict`
  - `self.filter_data`
  - `self.get_metric_of_the_whole_swarm`
  - `self.number_of_chromosome`
  - `self.number_of_free_variables`
  - `self.swarm_data_as_dict`
  - `self.swarm_data_project_names`
  - `self.swarm_data_xf`
  - `self.swarm_data_xf.append`
  - `self.swarm_data_xf['0']`
- **State Outputs (Writes)**:
  - `self.number_of_chromosome`
  - `self.number_of_free_variables`
  - `self.swarm_data_as_dict`
  - `self.swarm_data_project_names`
  - `self.swarm_data_xf`

### `Swarm_Data_Analyzer.filter_data` (Line 1976)
- **Arguments**: self, data, param_type, filter_key, direction, filter_value
- **State Inputs (Reads)**:
  - `self.decode_py_reduce_ordered_dict`
- **State Outputs**: None detected

### `Swarm_Data_Analyzer.decode` (Line 1999)
- **Arguments**: d
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Swarm_Data_Analyzer.get_metric_of_the_whole_swarm` (Line 2005)
- **Arguments**: self, metric
- **State Inputs (Reads)**:
  - `result.append`
  - `self.swarm_data_as_dict`
  - `self.swarm_data_as_dict.items`
- **State Outputs**: None detected

### `Swarm_Data_Analyzer.prepare_data_for_post_processing` (Line 2025)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.Cost`
  - `self.Cost_Cu`
  - `self.Cost_Fe`
  - `self.Cost_PM`
  - `self.Ea`
  - `self.Em`
  - `self.FRW`
  - `self.PowerFactor`
  - `self.RatedEfficiency`
  - `self.RatedStkLen`
  - `self.TRV`
  - `self.TorqueRipple`
  - `self.Tripple`
  - `self.f1`
  - `self.f2`
  - `self.f3`
  - `self.force_error_angle`
  - `self.get_metric_of_the_whole_swarm`
  - `self.l_rated_iron_loss`
  - `self.l_rated_magnet_Joule_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.normalized_force_error_magnitude`
  - `self.rotor_weight`
  - `self.ss_avg_force_magnitude`
  - `self.swarm_data_xf`
  - `self.torque_average`
- **State Outputs (Writes)**:
  - `self.Cost`
  - `self.Cost_Cu`
  - `self.Cost_Fe`
  - `self.Cost_PM`
  - `self.Ea`
  - `self.Em`
  - `self.FRW`
  - `self.PowerFactor`
  - `self.RatedEfficiency`
  - `self.RatedStkLen`
  - `self.TRV`
  - `self.TorqueRipple`
  - `self.Tripple`
  - `self.f1`
  - `self.f2`
  - `self.f3`
  - `self.force_error_angle`
  - `self.l_rated_iron_loss`
  - `self.l_rated_magnet_Joule_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.normalized_force_error_magnitude`
  - `self.rotor_weight`
  - `self.ss_avg_force_magnitude`
  - `self.torque_average`

### `swarm_data_container.__init__` (Line 2111)
- **Arguments**: self, swarm_data_raw, fea_config_dict, swarm_data_json, swarm_data_json_file_path
- **State Inputs (Reads)**:
  - `self._initialize_empty`
  - `self._load_from_json`
  - `self._load_from_raw`
  - `self.fea_config_dict`
  - `self.swarm_data_raw`
- **State Outputs (Writes)**:
  - `self.fea_config_dict`
  - `self.swarm_data_raw`

### `swarm_data_container._initialize_empty` (Line 2141)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.Ea`
  - `self.Em`
  - `self.FRW`
  - `self.RatedStkLen`
  - `self.RatedVol`
  - `self.RatedWeight`
  - `self.Trip`
  - `self.deg_alpha_st`
  - `self.l_FRW`
  - `self.l_OA`
  - `self.l_OB`
  - `self.l_OC`
  - `self.l_TRV`
  - `self.l_design_parameters`
  - `self.l_efficiency`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_normalized_torque_ripple`
  - `self.l_original_rotor_weight`
  - `self.l_original_stack_length`
  - `self.l_power_factor`
  - `self.l_rated_efficiency`
  - `self.l_rated_iron_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_rotor_weight`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.machine_data`
  - `self.mm_r_si`
  - `self.mm_w_st`
  - `self.number_of_free_variables`
  - `self.project_names`
  - `self.rated_data`
  - `self.swarm_data_xf`
- **State Outputs (Writes)**:
  - `self.Ea`
  - `self.Em`
  - `self.FRW`
  - `self.RatedStkLen`
  - `self.RatedVol`
  - `self.RatedWeight`
  - `self.Trip`
  - `self.deg_alpha_st`
  - `self.l_FRW`
  - `self.l_OA`
  - `self.l_OB`
  - `self.l_OC`
  - `self.l_TRV`
  - `self.l_design_parameters`
  - `self.l_efficiency`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_normalized_torque_ripple`
  - `self.l_original_rotor_weight`
  - `self.l_original_stack_length`
  - `self.l_power_factor`
  - `self.l_rated_efficiency`
  - `self.l_rated_iron_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_rotor_weight`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.machine_data`
  - `self.mm_r_si`
  - `self.mm_w_st`
  - `self.number_of_free_variables`
  - `self.project_names`
  - `self.rated_data`
  - `self.swarm_data_xf`

### `swarm_data_container._load_from_json` (Line 2188)
- **Arguments**: self, swarm_data_json
- **State Inputs (Reads)**:
  - `self.Ea`
  - `self.Ea.append`
  - `self.Em`
  - `self.Em.append`
  - `self.FRW`
  - `self.FRW.append`
  - `self.RatedStkLen`
  - `self.RatedStkLen.append`
  - `self.RatedVol`
  - `self.RatedVol.append`
  - `self.RatedWeight`
  - `self.RatedWeight.append`
  - `self.Trip`
  - `self.Trip.append`
  - `self._extract_performance_lists`
  - `self._initialize_empty`
  - `self.machine_data`
  - `self.machine_data.append`
  - `self.number_of_free_variables`
  - `self.project_names`
  - `self.project_names.append`
  - `self.rated_data`
  - `self.rated_data.append`
  - `self.swarm_data_xf`
  - `self.swarm_data_xf.append`
  - `self.swarm_data_xf['0']`
- **State Outputs (Writes)**:
  - `self.number_of_free_variables`

### `swarm_data_container._extract_performance_lists` (Line 2302)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.l_FRW`
  - `self.l_OA`
  - `self.l_OB`
  - `self.l_OC`
  - `self.l_TRV`
  - `self.l_design_parameters`
  - `self.l_efficiency`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_normalized_torque_ripple`
  - `self.l_original_rotor_weight`
  - `self.l_original_stack_length`
  - `self.l_power_factor`
  - `self.l_rated_efficiency`
  - `self.l_rated_iron_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_rotor_weight`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.machine_data`
  - `self.rated_data`
  - `self.swarm_data_xf`
- **State Outputs (Writes)**:
  - `self.l_FRW`
  - `self.l_OA`
  - `self.l_OB`
  - `self.l_OC`
  - `self.l_TRV`
  - `self.l_design_parameters`
  - `self.l_efficiency`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_normalized_torque_ripple`
  - `self.l_original_rotor_weight`
  - `self.l_original_stack_length`
  - `self.l_power_factor`
  - `self.l_rated_efficiency`
  - `self.l_rated_iron_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_rotor_weight`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`

### `swarm_data_container._load_from_raw` (Line 2347)
- **Arguments**: self, swarm_data_raw
- **State Inputs (Reads)**:
  - `self.Ea`
  - `self.Ea.append`
  - `self.Em`
  - `self.Em.append`
  - `self.FRW`
  - `self.FRW.append`
  - `self.RatedStkLen`
  - `self.RatedStkLen.append`
  - `self.RatedVol`
  - `self.RatedVol.append`
  - `self.RatedWeight`
  - `self.RatedWeight.append`
  - `self.Trip`
  - `self.Trip.append`
  - `self._extract_performance_lists`
  - `self._initialize_empty`
  - `self.deg_alpha_st`
  - `self.deg_alpha_st.append`
  - `self.l_FRW`
  - `self.l_OA`
  - `self.l_OB`
  - `self.l_OC`
  - `self.l_TRV`
  - `self.l_design_parameters`
  - `self.l_efficiency`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_normalized_torque_ripple`
  - `self.l_original_rotor_weight`
  - `self.l_original_stack_length`
  - `self.l_power_factor`
  - `self.l_rated_efficiency`
  - `self.l_rated_iron_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_rotor_weight`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.machine_data`
  - `self.machine_data.append`
  - `self.mm_r_si`
  - `self.mm_r_si.append`
  - `self.mm_w_st`
  - `self.mm_w_st.append`
  - `self.number_of_free_variables`
  - `self.project_names`
  - `self.project_names.append`
  - `self.rated_data`
  - `self.rated_data.append`
  - `self.swarm_data_xf`
  - `self.swarm_data_xf.append`
  - `self.swarm_data_xf['0']`
- **State Outputs (Writes)**:
  - `self.deg_alpha_st`
  - `self.l_FRW`
  - `self.l_OA`
  - `self.l_OB`
  - `self.l_OC`
  - `self.l_TRV`
  - `self.l_design_parameters`
  - `self.l_efficiency`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_normalized_torque_ripple`
  - `self.l_original_rotor_weight`
  - `self.l_original_stack_length`
  - `self.l_power_factor`
  - `self.l_rated_efficiency`
  - `self.l_rated_iron_loss`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_rotor_weight`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stack_length`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.mm_r_si`
  - `self.mm_w_st`
  - `self.number_of_free_variables`

### `swarm_data_container.get_list_y_data` (Line 2580)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.l_OB`
  - `self.l_TRV`
  - `self.l_force_error_angle`
- **State Outputs**: None detected

### `swarm_data_container.sensitivity_bar_charts` (Line 2593)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `O1_max.append`
  - `O1_min.append`
  - `O2_ecce_data.append`
  - `O2_max.append`
  - `O2_min.append`
  - `O2_prototype_data.append`
  - `data_max.append`
  - `data_min.append`
  - `results_for_refining_bounds['O1'].append`
  - `results_for_refining_bounds['O2'].append`
  - `self.fea_config_dict`
  - `self.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
  - `self.get_certain_objective_function`
  - `self.l_OA`
  - `self.l_OB`
  - `self.l_OC`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_normalized_torque_ripple`
  - `self.l_original_rotor_weight`
  - `self.l_rated_total_loss`
  - `self.l_ss_avg_force_magnitude`
  - `self.number_of_free_variables`
  - `self.reference_data`
  - `self.reference_data['2']`
  - `self.reference_data['3']`
  - `self.reference_data['4']`
  - `self.reference_data['5']`
  - `self.reference_data['6']`
  - `self.reference_design`
  - `self.reference_design['3']`
  - `self.reference_design['3'].split`
  - `self.required_torque`
  - `self.rotor_volume`
  - `self.rotor_weight`
- **State Outputs (Writes)**:
  - `self.reference_data`


## File: `utility.py`

### `my_execfile` (Line 12)
- **Arguments**: filename, g, l
- **State Inputs**: None detected
- **State Outputs**: None detected

### `json_dump_ignoring_unserializable` (Line 16)
- **Arguments**: obj
- **State Inputs**: None detected
- **State Outputs**: None detected

### `ExceptionReTry.__init__` (Line 43)
- **Arguments**: self, message, payload
- **State Inputs (Reads)**:
  - `self.message`
  - `self.payload`
- **State Outputs (Writes)**:
  - `self.message`
  - `self.payload`

### `ExceptionReTry.__str__` (Line 46)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.message`
- **State Outputs**: None detected

### `ExceptionBadNumberOfParts.__init__` (Line 51)
- **Arguments**: self, message, payload
- **State Inputs (Reads)**:
  - `self.message`
  - `self.payload`
- **State Outputs (Writes)**:
  - `self.message`
  - `self.payload`

### `ExceptionBadNumberOfParts.__str__` (Line 54)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.message`
- **State Outputs**: None detected

### `communicate_database` (Line 57)
- **Arguments**: spec
- **State Inputs**: None detected
- **State Outputs**: None detected

### `get_index_and_max` (Line 149)
- **Arguments**: the_list
- **State Inputs**: None detected
- **State Outputs**: None detected

### `get_index_and_min` (Line 153)
- **Arguments**: the_list
- **State Inputs**: None detected
- **State Outputs**: None detected

### `gcd` (Line 157)
- **Arguments**: a, b
- **State Inputs**: None detected
- **State Outputs**: None detected

### `lcm` (Line 162)
- **Arguments**: a, b
- **State Inputs**: None detected
- **State Outputs**: None detected

### `myLogger` (Line 167)
- **Arguments**: dir_log, prefix
- **State Inputs**: None detected
- **State Outputs**: None detected

### `blockPrinting` (Line 218)
- **Arguments**: func
- **State Inputs (Reads)**:
  - `func_wrapper`
- **State Outputs**: None detected

### `blockPrint` (Line 229)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected

### `enablePrint` (Line 232)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected

### `to_precision` (Line 238)
- **Arguments**: x, p
- **State Inputs (Reads)**:
  - `out.append`
- **State Outputs**: None detected

### `singleSidedDFT` (Line 301)
- **Arguments**: signal, samp_freq
- **State Inputs**: None detected
- **State Outputs**: None detected

### `basefreqDFT` (Line 309)
- **Arguments**: signal, samp_freq, ax_time_domain, ax_freq_domain, base_freq
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Pyrhonen_design.__init__` (Line 334)
- **Arguments**: self, im, bounds
- **State Inputs (Reads)**:
  - `self.Angle_StatorSlotOpen`
  - `self.Length_HeadNeckRotorSlot`
  - `self.SIesign_parameters_denorm`
  - `self.Width_StatorTeethHeadThickness`
  - `self.air_gap_length_delta`
  - `self.b1`
  - `self.rotor_tooth_width_b_dr`
  - `self.show_norm`
  - `self.stator_tooth_width_b_ds`
- **State Outputs (Writes)**:
  - `self.Angle_StatorSlotOpen`
  - `self.Length_HeadNeckRotorSlot`
  - `self.SIesign_parameters_denorm`
  - `self.Width_StatorTeethHeadThickness`
  - `self.air_gap_length_delta`
  - `self.b1`
  - `self.rotor_tooth_width_b_dr`
  - `self.stator_tooth_width_b_ds`

### `Pyrhonen_design.show_denorm` (Line 375)
- **Arguments**: self, bounds, design_parameters_norm
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Pyrhonen_design.show_norm` (Line 383)
- **Arguments**: self, bounds, design_parameters_denorm
- **State Inputs (Reads)**:
  - `self.SIesign_parameters_norm`
  - `self.SIesign_parameters_norm.tolist`
- **State Outputs (Writes)**:
  - `self.SIesign_parameters_norm`

### `add_Pyrhonen_design_to_first_generation` (Line 402)
- **Arguments**: sw, de_config_dict, logger
- **State Inputs**: None detected
- **State Outputs**: None detected

### `send_notification` (Line 424)
- **Arguments**: text, subject
- **State Inputs**: None detected
- **State Outputs**: None detected

### `get_windage_loss` (Line 436)
- **Arguments**: im_variant, mm_stack_length, TEMPERATURE_OF_AIR
- **State Inputs (Reads)**:
  - `im_variant.winding.EX['RatedSpeed']`
- **State Outputs**: None detected

### `suspension_force_vector.__init__` (Line 515)
- **Arguments**: self, force_x, force_y, range_ss
- **State Inputs (Reads)**:
  - `self.force_abs`
  - `self.force_ang`
  - `self.force_ang.append`
  - `self.force_err_abs`
  - `self.force_err_ang`
  - `self.force_err_ang_new_way`
  - `self.force_err_ang_old_way`
  - `self.force_error_angle`
  - `self.force_x`
  - `self.force_y`
  - `self.normalized_force_error_magnitude`
  - `self.range_ss`
  - `self.ss_avg_force_angle`
  - `self.ss_avg_force_magnitude`
  - `self.ss_avg_force_vector`
  - `self.ss_avg_force_vector['0']`
  - `self.ss_avg_force_vector['1']`
  - `self.ss_max_force_err_abs`
  - `self.ss_max_force_err_abs['0']`
  - `self.ss_max_force_err_abs['1']`
  - `self.ss_max_force_err_ang`
  - `self.ss_max_force_err_ang['0']`
  - `self.ss_max_force_err_ang['1']`
- **State Outputs (Writes)**:
  - `self.force_abs`
  - `self.force_ang`
  - `self.force_err_abs`
  - `self.force_err_ang`
  - `self.force_err_ang_new_way`
  - `self.force_err_ang_old_way`
  - `self.force_error_angle`
  - `self.force_x`
  - `self.force_y`
  - `self.normalized_force_error_magnitude`
  - `self.range_ss`
  - `self.ss_avg_force_angle`
  - `self.ss_avg_force_magnitude`
  - `self.ss_avg_force_vector`
  - `self.ss_max_force_err_abs`
  - `self.ss_max_force_err_ang`

### `pyplot_clear` (Line 565)
- **Arguments**: axeses
- **State Inputs**: None detected
- **State Outputs**: None detected

### `read_csv_results_4_comparison__transient` (Line 580)
- **Arguments**: study_name, path_prefix
- **State Inputs (Reads)**:
  - `ForConX_list.append`
  - `ForConY_list.append`
  - `TorCon_list.append`
  - `basic_info.append`
  - `key_list.append`
  - `time_list.append`
- **State Outputs**: None detected

### `read_csv_results_4_comparison_eddycurrent` (Line 653)
- **Arguments**: study_name, path_prefix
- **State Inputs (Reads)**:
  - `ForConX_list.append`
  - `ForConY_list.append`
  - `TorCon_list.append`
- **State Outputs**: None detected

### `collect_jmag_Tran2TSSProlong_results` (Line 689)
- **Arguments**: im_variant, path_prefix, fea_config_dict, axeses, femm_solver_data
- **State Inputs**: None detected
- **State Outputs**: None detected

### `csv_row_reader` (Line 787)
- **Arguments**: handle
- **State Inputs**: None detected
- **State Outputs**: None detected

### `whole_row_reader` (Line 792)
- **Arguments**: reader
- **State Inputs**: None detected
- **State Outputs**: None detected

### `get_copper_loss_Bolognani` (Line 798)
- **Arguments**: stator_slot_area, rotor_slot_area, STATOR_SLOT_FILL_FACTOR, ROTOR_SLOT_FILL_FACTOR, TEMPERATURE_OF_COIL, copper_loss_parameters
- **State Inputs (Reads)**:
  - `EX['DriveW_zQ']`
  - `EX['Js']`
  - `EX['WindingFill']`
  - `EX['wily']`
  - `EX['wily'].number_parallel_branch`
  - `self.im`
  - `self.im.DriveW_poles`
  - `self.im.Qr`
  - `self.im.design_parameters`
  - `self.im.design_parameters['2']`
  - `self.im.rotor_slot_height_h_sr`
  - `self.im.stack_length`
  - `self.im.template`
  - `self.im.template.SI`
  - `self.im.template.SI['GP']`
  - `self.im.template.SI['GP']['mm_r_ro']`
  - `self.im.template.SI['GP']['mm_r_ro'].value`
  - `self.list_rotor_current_amp`
- **State Outputs**: None detected

### `check_csv_results_4_general_purpose` (Line 900)
- **Arguments**: study_name, path_prefix, returnBoolean
- **State Inputs (Reads)**:
  - `l_ForCon_X.append`
  - `l_ForCon_Y.append`
  - `l_TorCon.append`
  - `l_slip_freq.append`
- **State Outputs**: None detected

### `Goertzel_Data_Struct.__init__` (Line 961)
- **Arguments**: self, id
- **State Inputs (Reads)**:
  - `self.accumSquaredData`
  - `self.ampl`
  - `self.bool_initialized`
  - `self.coeff`
  - `self.cosine`
  - `self.count`
  - `self.id`
  - `self.imag`
  - `self.k`
  - `self.phase`
  - `self.q`
  - `self.q2`
  - `self.real`
  - `self.scalingFactor`
  - `self.sine`
- **State Outputs (Writes)**:
  - `self.accumSquaredData`
  - `self.ampl`
  - `self.bool_initialized`
  - `self.coeff`
  - `self.cosine`
  - `self.count`
  - `self.id`
  - `self.imag`
  - `self.k`
  - `self.phase`
  - `self.q`
  - `self.q2`
  - `self.real`
  - `self.scalingFactor`
  - `self.sine`

### `Goertzel_Data_Struct.goertzel_realtime` (Line 983)
- **Arguments**: gs, targetFreq, numSamples, samplingRate, data
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Goertzel_Data_Struct.goertzel_offline` (Line 1024)
- **Arguments**: gs, targetFreq, samplingRate, data_list
- **State Inputs**: None detected
- **State Outputs**: None detected

### `compute_power_factor_from_half_period` (Line 1065)
- **Arguments**: voltage, current, mytime, targetFreq, numPeriodicalExtension
- **State Inputs**: None detected
- **State Outputs**: None detected

### `compute_power_factor_from_full_period` (Line 1100)
- **Arguments**: voltage, current, mytime, targetFreq, numPeriodicalExtension
- **State Inputs**: None detected
- **State Outputs**: None detected

### `max_indices_2` (Line 1131)
- **Arguments**: arr, k
- **State Inputs (Reads)**:
  - `max_idxs.append`
- **State Outputs**: None detected

### `min_indices` (Line 1149)
- **Arguments**: arr, k
- **State Inputs**: None detected
- **State Outputs**: None detected

### `max_indices` (Line 1156)
- **Arguments**: arr, k
- **State Inputs**: None detected
- **State Outputs**: None detected

### `autolabel` (Line 1166)
- **Arguments**: ax, rects, xpos, bias, textfont
- **State Inputs**: None detected
- **State Outputs**: None detected

### `efficiency_at_50kW` (Line 1183)
- **Arguments**: total_loss
- **State Inputs**: None detected
- **State Outputs**: None detected

### `use_weights` (Line 1186)
- **Arguments**: which
- **State Inputs**: None detected
- **State Outputs**: None detected

### `compute_list_cost` (Line 1199)
- **Arguments**: weights, rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss
- **State Inputs**: None detected
- **State Outputs**: None detected

### `fobj_scalar` (Line 1225)
- **Arguments**: torque_average, ss_avg_force_magnitude, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle, total_loss, weights, rotor_volume, rotor_weight
- **State Inputs**: None detected
- **State Outputs**: None detected

### `fobj_list` (Line 1237)
- **Arguments**: l_torque_average, l_ss_avg_force_magnitude, l_normalized_torque_ripple, l_normalized_force_error_magnitude, l_force_error_angle, l_total_loss, weights, rotor_volume, rotor_weight
- **State Inputs (Reads)**:
  - `l_cost_function.append`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.__init__` (Line 1256)
- **Arguments**: self, sw, spec, dir_run, run_integer, bool_sensitivity_analysis
- **State Inputs (Reads)**:
  - `self.SIir_run`
  - `self.buf`
  - `self.buf_length`
  - `self.build_basic_info`
  - `self.number_of_designs`
  - `self.reference_design`
  - `self.run_integer`
  - `self.spec`
  - `self.sw`
- **State Outputs (Writes)**:
  - `self.SIir_run`
  - `self.buf`
  - `self.buf_length`
  - `self.number_of_designs`
  - `self.reference_design`
  - `self.run_integer`
  - `self.spec`
  - `self.sw`

### `SwarmDataAnalyzer.design_display_generator` (Line 1290)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.design_parameters_generator` (Line 1294)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.list_generations` (Line 1298)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.list_cost_function` (Line 1309)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `l.append`
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.find_individual` (Line 1315)
- **Arguments**: self, generation_index, individual_index
- **State Inputs (Reads)**:
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.get_best_generation` (Line 1326)
- **Arguments**: self, popsize, generator, returnMore
- **State Inputs (Reads)**:
  - `self.SIesign_parameters_generator`
  - `self.list_cost_function`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.get_list_objective_function` (Line 1353)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.get_certain_objective_function` (Line 1360)
- **Arguments**: self, which
- **State Inputs (Reads)**:
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.get_windage_loss` (Line 1369)
- **Arguments**: self, which
- **State Inputs (Reads)**:
  - `self.buf`
  - `self.number_of_designs`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.my_population_distribution_plots` (Line 1379)
- **Arguments**: self, de_config_dict
- **State Inputs (Reads)**:
  - `self.SIesign_parameters_generator`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.my_scatter_plot` (Line 1440)
- **Arguments**: self, x, y, O, fig, ax, s, marker, index_list
- **State Inputs (Reads)**:
  - `self.ExcitationFreqSimulated`
  - `self.Omega`
  - `self.Qr`
  - `self.Qs`
  - `self.SIesign_display_generator`
  - `self.SIesign_parameters_generator`
  - `self.best_design_denorm`
  - `self.best_design_denorm['0']`
  - `self.best_design_display`
  - `self.best_design_display.split`
  - `self.mec_power`
  - `self.required_torque`
  - `self.rotor_volume`
  - `self.rotor_weight`
  - `self.run_integer`
  - `self.spec`
  - `self.spec.Jr`
  - `self.spec.Js`
  - `self.spec.Steel`
  - `self.spec.VoltageRating`
  - `self.spec.stator_phase_current_rms`
  - `self.speed_rpm`
  - `self.stack_length`
  - `self.stack_length_max`
  - `self.str_best_design_details`
  - `self.sw`
  - `self.sw.fea_config_dict`
  - `self.sw.fea_config_dict['use_weights']`
  - `self.sw.im`
  - `self.sw.im.BeariW_poles`
  - `self.sw.im.DriveW_poles`
  - `self.sw.im.Radius_OuterStatorYoke`
  - `self.sw.im.template`
  - `self.sw.im.template.SI`
  - `self.sw.im.template.SI['GP']`
  - `self.sw.im.template.SI['GP']['mm_r_ro']`
  - `self.sw.im.template.SI['GP']['mm_r_ro'].value`
  - `self.template`
  - `self.template.SI`
  - `self.template.SI['GP']`
  - `self.template.SI['GP']['mm_r_ro']`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.weights_name`
  - `self.weights_used`
- **State Outputs (Writes)**:
  - `self.best_design_denorm`
  - `self.best_design_display`
  - `self.str_best_design_details`

### `SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth` (Line 1560)
- **Arguments**: self, fig, ax, marker, bool_filtered
- **State Inputs (Reads)**:
  - `filtered_O_list.append`
  - `filtered_filtered_O_list.append`
  - `filtered_index_list.append`
  - `filtered_loss_list.append`
  - `filtered_torque_list.append`
  - `filtered_x.append`
  - `filtered_y.append`
  - `index_list.append`
  - `l_rated_stack_length.append`
  - `l_rated_total_loss.append`
  - `self.get_certain_objective_function`
  - `self.mec_power`
  - `self.my_scatter_plot`
  - `self.required_torque`
  - `self.rotor_volume`
  - `self.rotor_weight`
  - `self.stack_length`
  - `self.stack_length_max`
  - `self.weights_used`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.pareto_plot_torque_force` (Line 1648)
- **Arguments**: self, fig2, axeses, marker
- **State Inputs (Reads)**:
  - `self.get_certain_objective_function`
  - `self.my_scatter_plot`
  - `self.rotor_volume`
  - `self.rotor_weight`
  - `self.weights_name`
  - `self.weights_used`
- **State Outputs**: None detected

### `SwarmDataAnalyzer.sensitivity_bar_charts` (Line 1765)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `O1_max.append`
  - `O1_min.append`
  - `O2_ecce_data.append`
  - `O2_max.append`
  - `O2_min.append`
  - `O2_prototype_data.append`
  - `data_max.append`
  - `data_min.append`
  - `results_for_refining_bounds['O1'].append`
  - `results_for_refining_bounds['O2'].append`
  - `self.get_certain_objective_function`
  - `self.reference_data`
  - `self.reference_data['2']`
  - `self.reference_data['3']`
  - `self.reference_data['4']`
  - `self.reference_data['5']`
  - `self.reference_data['6']`
  - `self.reference_design`
  - `self.reference_design['3']`
  - `self.reference_design['3'].split`
  - `self.required_torque`
  - `self.rotor_volume`
  - `self.rotor_weight`
  - `self.sw`
  - `self.sw.fea_config_dict`
  - `self.sw.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
- **State Outputs (Writes)**:
  - `self.reference_data`

### `SwarmDataAnalyzer.build_basic_info` (Line 2127)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.ExcitationFreqSimulated`
  - `self.Omega`
  - `self.Qr`
  - `self.Qs`
  - `self.mec_power`
  - `self.required_torque`
  - `self.rotor_volume`
  - `self.rotor_weight`
  - `self.spec`
  - `self.speed_rpm`
  - `self.stack_length`
  - `self.stack_length_max`
  - `self.sw`
  - `self.template`
  - `self.template.SI`
  - `self.template.SI['GP']`
  - `self.template.SI['GP']['mm_r_ro']`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.weights_name`
  - `self.weights_used`
- **State Outputs (Writes)**:
  - `self.ExcitationFreqSimulated`
  - `self.Omega`
  - `self.Qr`
  - `self.Qs`
  - `self.mec_power`
  - `self.required_torque`
  - `self.rotor_volume`
  - `self.rotor_weight`
  - `self.speed_rpm`
  - `self.stack_length`
  - `self.stack_length_max`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.weights_name`
  - `self.weights_used`

### `build_sensitivity_bar_charts` (Line 2159)
- **Arguments**: spec, sw
- **State Inputs**: None detected
- **State Outputs**: None detected

### `build_Pareto_plot` (Line 2164)
- **Arguments**: spec, sw
- **State Inputs**: None detected
- **State Outputs**: None detected


## File: `winding_layout.py`

### `infer_Y_layer_phases_from_X_layer_and_coil_pitch_y` (Line 9)
- **Arguments**: layer_X_phases, coil_pitch
- **State Inputs**: None detected
- **State Outputs**: None detected

### `infer_Y_layer_signs_from_X_layer_and_coil_pitch_y` (Line 11)
- **Arguments**: layer_X_signs, coil_pitch
- **State Inputs**: None detected
- **State Outputs**: None detected

### `infer_Y_layer_grpAC_from_X_layer_and_coil_pitch_y` (Line 14)
- **Arguments**: grouping_AC, coil_pitch
- **State Inputs**: None detected
- **State Outputs**: None detected

### `winding_layout_v2.__init__` (Line 21)
- **Arguments**: self, DPNV_or_SEPA, Qs, p, ps, coil_pitch_y, pr, m, Wrap_Around
- **State Inputs (Reads)**:
  - `self.CommutatingSequenceB`
  - `self.CommutatingSequenceD`
  - `self.Qs`
  - `self.SIPNV_or_SEPA`
  - `self.SPP`
  - `self.bool_3PhaseCurrentSource`
  - `self.coil_pitch_y`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.dict_coil_connection`
  - `self.distributed_or_concentrated`
  - `self.grouping_AC`
  - `self.kd1`
  - `self.kp1`
  - `self.layer_X_phases`
  - `self.layer_X_signs`
  - `self.layer_Y_phases`
  - `self.layer_Y_signs`
  - `self.list_layer_motor_phases`
  - `self.list_layer_motor_signs`
  - `self.list_layer_suspension_phases`
  - `self.list_layer_suspension_signs`
  - `self.m`
  - `self.number_parallel_branch`
  - `self.number_winding_layer`
  - `self.ox_distribution_phase_U`
  - `self.ox_distribution_three_phase`
  - `self.p`
  - `self.pr`
  - `self.ps`
- **State Outputs (Writes)**:
  - `self.CommutatingSequenceB`
  - `self.CommutatingSequenceD`
  - `self.Qs`
  - `self.SIPNV_or_SEPA`
  - `self.SPP`
  - `self.bool_3PhaseCurrentSource`
  - `self.coil_pitch_y`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.dict_coil_connection`
  - `self.distributed_or_concentrated`
  - `self.grouping_AC`
  - `self.kd1`
  - `self.kp1`
  - `self.layer_X_phases`
  - `self.layer_X_signs`
  - `self.layer_Y_phases`
  - `self.layer_Y_signs`
  - `self.list_layer_motor_phases`
  - `self.list_layer_motor_signs`
  - `self.list_layer_suspension_phases`
  - `self.list_layer_suspension_signs`
  - `self.m`
  - `self.number_parallel_branch`
  - `self.number_winding_layer`
  - `self.ox_distribution_phase_U`
  - `self.ox_distribution_three_phase`
  - `self.p`
  - `self.pr`
  - `self.ps`

### `pole_specific_winding_with_neutral.__init__` (Line 1187)
- **Arguments**: self, Qr, p, ps, coil_pitch_y
- **State Inputs (Reads)**:
  - `self.pairs`
- **State Outputs (Writes)**:
  - `self.pairs`

### `nextpow2` (Line 1277)
- **Arguments**: L
- **State Inputs**: None detected
- **State Outputs**: None detected

### `periodic2pi` (Line 1283)
- **Arguments**: x
- **State Inputs**: None detected
- **State Outputs**: None detected

### `segmented_func` (Line 1294)
- **Arguments**: x, lst_x, lst_y
- **State Inputs**: None detected
- **State Outputs**: None detected

### `PhaseWinding.__init__` (Line 1319)
- **Arguments**: self, Qs, m, turns_per_slot, ox_distribution_phase_U, desc_type
- **State Inputs (Reads)**:
  - `self.avg_val_of_turn_func`
  - `self.degree_between_slots`
  - `self.ox_distribution_phase_U`
  - `self.radian_between_slots`
  - `self.setSymPos`
  - `self.setTurnFuncObject`
  - `self.slot_per_phase`
  - `self.sym_begin_pos`
  - `self.sym_turn_func`
  - `self.sym_winding_func`
  - `self.turn_func`
  - `self.turns_per_slot`
  - `self.winding_func`
- **State Outputs (Writes)**:
  - `self.avg_val_of_turn_func`
  - `self.degree_between_slots`
  - `self.ox_distribution_phase_U`
  - `self.radian_between_slots`
  - `self.slot_per_phase`
  - `self.sym_turn_func`
  - `self.sym_winding_func`
  - `self.turns_per_slot`
  - `self.winding_func`

### `PhaseWinding.setTurnFuncObject` (Line 1346)
- **Arguments**: self, ox_distribution_phase_U
- **State Inputs (Reads)**:
  - `lst_x.append`
  - `lst_y.append`
  - `self.lst_x`
  - `self.lst_y`
  - `self.radian_between_slots`
  - `self.turn_func`
  - `self.turns_per_slot`
- **State Outputs (Writes)**:
  - `self.lst_x`
  - `self.lst_y`
  - `self.turn_func`

### `PhaseWinding.setSymPos` (Line 1378)
- **Arguments**: self, index
- **State Inputs (Reads)**:
  - `self.lst_x`
  - `self.lst_y`
  - `self.lst_y.index`
  - `self.sym_begin_pos`
  - `self.sym_begin_pos_1`
  - `self.sym_begin_pos_2`
- **State Outputs (Writes)**:
  - `self.sym_begin_pos`
  - `self.sym_begin_pos_1`
  - `self.sym_begin_pos_2`

### `PhaseWinding.plot2piFft` (Line 1403)
- **Arguments**: self, func, Fs, L
- **State Inputs (Reads)**:
  - `self.fig_plot2piFft`
- **State Outputs (Writes)**:
  - `self.fig_plot2piFft`

### `PhaseWinding.plotFuncObj` (Line 1457)
- **Arguments**: self, func
- **State Inputs (Reads)**:
  - `self.fig_plotFuncObj`
- **State Outputs (Writes)**:
  - `self.fig_plotFuncObj`

### `winding_layout.__init__` (Line 1521)
- **Arguments**: self, DPNV_or_SEPA, Qs, p, ps
- **State Inputs (Reads)**:
  - `self.CommutatingSequenceB`
  - `self.CommutatingSequenceD`
  - `self.Qs`
  - `self.bool_3PhaseCurrentSource`
  - `self.coil_pitch`
  - `self.distributed_or_concentrated`
  - `self.grouping_AC`
  - `self.initial_excitation_bias_compensation_deg`
  - `self.l21`
  - `self.l22`
  - `self.l41`
  - `self.l42`
  - `self.l_leftlayer1`
  - `self.l_leftlayer2`
  - `self.l_rightlayer1`
  - `self.l_rightlayer2`
  - `self.layer_A1`
  - `self.layer_A2`
  - `self.layer_B1`
  - `self.layer_B2`
  - `self.no_winding_layer`
  - `self.number_parallel_branch`
  - `self.p`
- **State Outputs (Writes)**:
  - `self.CommutatingSequenceB`
  - `self.CommutatingSequenceD`
  - `self.Qs`
  - `self.bool_3PhaseCurrentSource`
  - `self.coil_pitch`
  - `self.distributed_or_concentrated`
  - `self.grouping_AC`
  - `self.initial_excitation_bias_compensation_deg`
  - `self.l21`
  - `self.l22`
  - `self.l41`
  - `self.l42`
  - `self.l_leftlayer1`
  - `self.l_leftlayer2`
  - `self.l_rightlayer1`
  - `self.l_rightlayer2`
  - `self.layer_A1`
  - `self.layer_A2`
  - `self.layer_B1`
  - `self.layer_B2`
  - `self.no_winding_layer`
  - `self.number_parallel_branch`
  - `self.p`


## File: `winding_layout_derivation_ismb2021_asymetry_no_drawing.py`

### `_print` (Line 18)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected

### `set_verbose` (Line 24)
- **Arguments**: verbose
- **State Inputs**: None detected
- **State Outputs**: None detected

### `limit_to_360_deg` (Line 30)
- **Arguments**: PHI
- **State Inputs**: None detected
- **State Outputs**: None detected

### `belong_to_which_phase_belt` (Line 38)
- **Arguments**: PHI, phase_belt
- **State Inputs**: None detected
- **State Outputs**: None detected

### `belong_to_band` (Line 85)
- **Arguments**: LB, UB, PHI
- **State Inputs**: None detected
- **State Outputs**: None detected

### `phase_angle_of_slot_i_at_frequency_h` (Line 102)
- **Arguments**: slot_number, h, Q
- **State Inputs**: None detected
- **State Outputs**: None detected

### `compute_star_of_slots` (Line 108)
- **Arguments**: Q, p, m, verbose
- **State Inputs**: None detected
- **State Outputs**: None detected

### `compute_connection_star_at_another_frequency` (Line 146)
- **Arguments**: connection_star_raw_dict, frequency_ratio, which_phase, verbose
- **State Inputs**: None detected
- **State Outputs**: None detected

### `winding_distribution_factor` (Line 206)
- **Arguments**: Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict, Aa
- **State Inputs**: None detected
- **State Outputs**: None detected

### `winding_short_pitch_factor_v2` (Line 243)
- **Arguments**: h, coil_pitch_y, Q
- **State Inputs**: None detected
- **State Outputs**: None detected

### `Winding_Derivation.__init__` (Line 253)
- **Arguments**: self, slot_pole_comb, bool_double_layer_winding, verbose
- **State Inputs (Reads)**:
  - `self.Q`
  - `self.bool_double_layer_winding`
  - `self.coil_pitch_y`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.dpnv_grouping_dict_b`
  - `self.dpnv_grouping_dict_c`
  - `self.list_phase_u_slot_number`
  - `self.list_phase_v_slot_number`
  - `self.list_phase_w_slot_number`
  - `self.list_slot_number_of_phase`
  - `self.m`
  - `self.p`
  - `self.ps`
  - `self.q`
  - `self.qs`
  - `self.suspen_kd_at_h`
  - `self.suspen_kp_at_h`
  - `self.suspen_kw_at_h`
  - `self.t`
  - `self.torque_kd_at_h`
  - `self.torque_kp_at_h`
  - `self.torque_kw_at_h`
  - `self.ts`
  - `self.turn_func_bias`
  - `self.verbose`
- **State Outputs (Writes)**:
  - `self.Q`
  - `self.bool_double_layer_winding`
  - `self.coil_pitch_y`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.dpnv_grouping_dict_b`
  - `self.dpnv_grouping_dict_c`
  - `self.list_phase_u_slot_number`
  - `self.list_phase_v_slot_number`
  - `self.list_phase_w_slot_number`
  - `self.list_slot_number_of_phase`
  - `self.m`
  - `self.p`
  - `self.ps`
  - `self.q`
  - `self.qs`
  - `self.suspen_kd_at_h`
  - `self.suspen_kp_at_h`
  - `self.suspen_kw_at_h`
  - `self.t`
  - `self.torque_kd_at_h`
  - `self.torque_kp_at_h`
  - `self.torque_kw_at_h`
  - `self.ts`
  - `self.turn_func_bias`
  - `self.verbose`

### `Winding_Derivation.get_complex_number_winding_factor_of_coil_i` (Line 444)
- **Arguments**: self, i, coil_pitch_y, Q, v, p
- **State Inputs (Reads)**:
  - `self.verbose`
- **State Outputs**: None detected

### `Winding_Derivation.get_complex_number_kw_per_phase` (Line 458)
- **Arguments**: self, v, p, positive_connected_coils, negative_connected_coils
- **State Inputs (Reads)**:
  - `kp_cjh_list.append`
  - `kp_els_list.append`
  - `self.Q`
  - `self.coil_pitch_y`
  - `self.get_complex_number_winding_factor_of_coil_i`
  - `self.verbose`
- **State Outputs**: None detected

### `Winding_Derivation.get_complex_number_kw` (Line 486)
- **Arguments**: self, p_or_ps, v, bool_study_suspension_subharmonics
- **State Inputs (Reads)**:
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.dpnv_grouping_dict_a['GAC']`
  - `self.dpnv_grouping_dict_a['GBD']`
  - `self.dpnv_grouping_dict_b`
  - `self.dpnv_grouping_dict_b['GAC']`
  - `self.dpnv_grouping_dict_b['GBD']`
  - `self.dpnv_grouping_dict_c`
  - `self.dpnv_grouping_dict_c['GAC']`
  - `self.dpnv_grouping_dict_c['GBD']`
  - `self.get_complex_number_kw_per_phase`
  - `self.ps`
  - `self.verbose`
- **State Outputs**: None detected

### `Winding_Derivation.format_print_out_string` (Line 531)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `reformat_wily_info`
  - `self.Q`
  - `self.bool_double_layer_winding`
  - `self.coil_pitch_y`
  - `self.connection_star_raw_dict`
  - `self.dict_suspension_kw_cjh`
  - `self.dict_suspension_kw_els`
  - `self.dict_suspension_kw_els['A_angle']`
  - `self.dict_suspension_kw_els['B_angle']`
  - `self.dict_suspension_kw_els['C_angle']`
  - `self.dict_torque_kw_cjh`
  - `self.dict_torque_kw_els`
  - `self.dpnv_grouping_dict_a`
  - `self.dpnv_grouping_dict_b`
  - `self.dpnv_grouping_dict_c`
  - `self.get_complex_number_kw`
  - `self.grouping_AC`
  - `self.layer_X_phases`
  - `self.layer_X_signs`
  - `self.m`
  - `self.p`
  - `self.print_out_string`
  - `self.ps`
  - `self.verbose`
- **State Outputs (Writes)**:
  - `self.coil_pitch_y`
  - `self.grouping_AC`
  - `self.layer_X_phases`
  - `self.layer_X_signs`
  - `self.print_out_string`

### `main_derivation` (Line 655)
- **Arguments**: m, Qs, p, ps, coil_pitch_y, verbose
- **State Inputs**: None detected
- **State Outputs**: None detected


## File: `WireSlot_v1.py`

### `calc_stator_geometry` (Line 11)
- **Arguments**: OD, ID, tooth_width, tooth_depth, yoke, liner, slots
- **State Inputs**: None detected
- **State Outputs**: None detected

### `calc_winding_capacity` (Line 32)
- **Arguments**: net_area, gross_area, awg, orthocyclic_factor, fill_heuristic
- **State Inputs**: None detected
- **State Outputs**: None detected

### `calc_thermal_load` (Line 64)
- **Arguments**: ID, z_slot, a_bare, J, slots
- **State Inputs**: None detected
- **State Outputs**: None detected

### `calc_bemf_constants` (Line 84)
- **Arguments**: La, yoke, z_slot, B_sat, poles, slots
- **State Inputs**: None detected
- **State Outputs**: None detected

### `calc_motor_losses` (Line 120)
- **Arguments**: La, tooth_width, turns_per_phase, current, rpm, awg, poles
- **State Inputs**: None detected
- **State Outputs**: None detected


## File: `WireSlot_v2.py`

### `MotorThermalAnalyzer.__init__` (Line 9)
- **Arguments**: self, stator_od, rotor_od, air_gap, tooth_depth, tooth_width, slots, liner
- **State Inputs (Reads)**:
  - `self.air_gap`
  - `self.liner`
  - `self.rotor_od`
  - `self.slots`
  - `self.stator_id`
  - `self.stator_od`
  - `self.tooth_depth`
  - `self.tooth_width`
  - `self.yoke_thickness`
- **State Outputs (Writes)**:
  - `self.air_gap`
  - `self.liner`
  - `self.rotor_od`
  - `self.slots`
  - `self.stator_id`
  - `self.stator_od`
  - `self.tooth_depth`
  - `self.tooth_width`
  - `self.yoke_thickness`

### `MotorThermalAnalyzer.calculate_slot_area` (Line 33)
- **Arguments**: self
- **State Inputs (Reads)**:
  - `self.slots`
  - `self.stator_id`
  - `self.tooth_depth`
  - `self.tooth_width`
- **State Outputs**: None detected

### `MotorThermalAnalyzer.get_wire_properties` (Line 44)
- **Arguments**: self, awg_size
- **State Inputs**: None detected
- **State Outputs**: None detected

### `MotorThermalAnalyzer.estimate_max_wires_in_slot` (Line 61)
- **Arguments**: self, d_od
- **State Inputs (Reads)**:
  - `self.liner`
  - `self.slots`
  - `self.stator_id`
  - `self.tooth_depth`
  - `self.tooth_width`
- **State Outputs**: None detected

### `MotorThermalAnalyzer.analyze_thermal_performance` (Line 104)
- **Arguments**: self, awg_size, target_j
- **State Inputs (Reads)**:
  - `self.calculate_slot_area`
  - `self.estimate_max_wires_in_slot`
  - `self.get_wire_properties`
  - `self.slots`
  - `self.stator_id`
  - `self.yoke_thickness`
- **State Outputs**: None detected


## File: `BH/Inspect_BH_curve.py`

### `script_dir` (Line 21)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected

### `is_bh_curve_file` (Line 25)
- **Arguments**: path
- **State Inputs**: None detected
- **State Outputs**: None detected

### `read_bh_data` (Line 31)
- **Arguments**: path
- **State Inputs (Reads)**:
  - `rows.append`
- **State Outputs**: None detected

### `write_bh_file` (Line 60)
- **Arguments**: path, data
- **State Inputs**: None detected
- **State Outputs**: None detected

### `ensure_bh_file` (Line 68)
- **Arguments**: source
- **State Inputs**: None detected
- **State Outputs**: None detected

### `collect_and_convert_bh_curves` (Line 83)
- **Arguments**: dir_path
- **State Inputs (Reads)**:
  - `results.append`
- **State Outputs**: None detected

### `plot_bh_curves` (Line 111)
- **Arguments**: labels_and_data, out_path
- **State Inputs**: None detected
- **State Outputs**: None detected

### `main` (Line 141)
- **Arguments**: None
- **State Inputs**: None detected
- **State Outputs**: None detected

