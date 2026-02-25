# Variables & State Dependencies Report

This report extracts parameters read and modified by functions/methods across the codebase, focusing on state variables like `self.*`, `machine.*`, `EX[...]`, dictionary `.get()`, etc.

## Categorized Variables summary

### Variables Read (Inputs)
- `EX['DriveW_zQ']`
- `EX['Js']`
- `EX['WindingFill']`
- `EX['end_winding_length_Lew']`
- `EX['end_winding_length_factor_kov']`
- `EX['stator_slot_area']`
- `EX['wily']`
- `EX['wily'].number_parallel_branch`
- `acm_variant.rotorMagnet.notched_rotor.p`
- `im_variant.winding.EX['RatedSpeed']`
- `machine.all_points`
- `machine.parts`
- `my_machine.geometry`
- `my_machine.materials`
- `my_machine.target`
- `my_machine.winding`
- `self.Angle_StatorSlotOpen`
- `self.CommutatingSequenceB`
- `self.CommutatingSequenceD`
- `self.Cost`
- `self.Cost_Cu`
- `self.Cost_Fe`
- `self.Cost_PM`
- `self.Current_dict`
- `self.Current_dict['Time(s)']`
- `self.DisplacementAngle_list`
- `self.Ea`
- `self.Ea.append`
- `self.Em`
- `self.Em.append`
- `self.ExcitationFreqSimulated`
- `self.FRW`
- `self.FRW.append`
- `self.FluxLinkage_dict`
- `self.ForConAbs_list`
- `self.ForConX_list`
- `self.ForConY_list`
- `self.GP`
- `self.GP.items`
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
- `self.JMAG_version_number`
- `self.JMAG_version_string`
- `self.Length_HeadNeckRotorSlot`
- `self.Omega`
- `self.PowerFactor`
- `self.Q`
- `self.Q_prime`
- `self.Qr`
- `self.Qs`
- `self.RP`
- `self.RatedEfficiency`
- `self.RatedStkLen`
- `self.RatedStkLen.append`
- `self.RatedVol`
- `self.RatedVol.append`
- `self.RatedWeight`
- `self.RatedWeight.append`
- `self.SI`
- `self.SIPNV_or_SEPA`
- `self.SIesign_display_generator`
- `self.SIesign_parameters_denorm`
- `self.SIesign_parameters_generator`
- `self.SIesign_parameters_norm`
- `self.SIesign_parameters_norm.tolist`
- `self.SIir_run`
- `self.SPP`
- `self.TRV`
- `self.TorCon_list`
- `self.TorqueRipple`
- `self.Trip`
- `self.Trip.append`
- `self.Tripple`
- `self.Width_StatorTeethHeadThickness`
- `self.__class__`
- `self.__class__.__module__`
- `self.__class__.__name__`
- `self.__dict__`
- `self.__dict__.items`
- `self._extract_performance_lists`
- `self._get_parameter_logger`
- `self._initialize_empty`
- `self._load_from_json`
- `self._load_from_raw`
- `self._next_index`
- `self.accumSquaredData`
- `self.add_circuit`
- `self.add_material`
- `self.add_part`
- `self.add_plots`
- `self.air_gap_length_delta`
- `self.all_points`
- `self.all_points.pole_count`
- `self.all_points.slot_count`
- `self.ampl`
- `self.app`
- `self.app.CreateGeometryEditor`
- `self.app.GetCurrentModel`
- `self.app.LaunchGeometryEditor`
- `self.app.Quit`
- `self.ass`
- `self.ass.CreateSketch`
- `self.ass.GetItem`
- `self.avg_val_of_turn_func`
- `self.b1`
- `self.bFillRegion`
- `self.bMirror`
- `self.basic_info`
- `self.bearing_winding_conductors_per_slot`
- `self.bearing_winding_current`
- `self.best_design_denorm`
- `self.best_design_denorm['0']`
- `self.best_design_display`
- `self.best_design_display.split`
- `self.bool_3PhaseCurrentSource`
- `self.bool_CustomizedCircuit`
- `self.bool_DPNVorSEPA`
- `self.bool_PermanentMagnet`
- `self.bool_RotorNotched`
- `self.bool_StatorSlotClosed`
- `self.bool_distributed_or_concentrated`
- `self.bool_double_layer_winding`
- `self.bool_initialized`
- `self.bool_jmagDeleteResultsAfterCalculation`
- `self.bool_suppressShaft`
- `self.bounds`
- `self.bounds['0']`
- `self.bounds['1']`
- `self.buf`
- `self.buf_length`
- `self.build_basic_info`
- `self.calc`
- `self.calc_bounds`
- `self.calculate_excitation_current`
- `self.checkGeomApp`
- `self.circuit_current`
- `self.coeff`
- `self.coil_fluxLinkage`
- `self.coil_pitch`
- `self.coil_pitch_y`
- `self.color`
- `self.components_make_region`
- `self.conductor_current_amplitude`
- `self.conductors_per_slot`
- `self.connection_star_raw_dict`
- `self.connection_type`
- `self.consts`
- `self.cosine`
- `self.count`
- `self.counter`
- `self.ctx`
- `self.ctx.arc`
- `self.ctx.arc_negative`
- `self.ctx.fill_preserve`
- `self.ctx.line_to`
- `self.ctx.move_to`
- `self.ctx.new_path`
- `self.ctx.paint`
- `self.ctx.restore`
- `self.ctx.rotate`
- `self.ctx.save`
- `self.ctx.scale`
- `self.ctx.set_line_cap`
- `self.ctx.set_line_width`
- `self.ctx.set_source_rgb`
- `self.ctx.set_source_rgba`
- `self.ctx.stroke`
- `self.ctx.transform`
- `self.d_air_gap`
- `self.d_magnet`
- `self.d_tooth`
- `self.dc_bus_voltage`
- `self.decode_py_reduce_ordered_dict`
- `self.defaultUnit`
- `self.deg_alpha_st`
- `self.deg_alpha_st.append`
- `self.deg_winding_U_phase_phase_axis_angle`
- `self.degree_between_slots`
- `self.derivation`
- `self.derivation.__dict__`
- `self.derivation.__dict__.items`
- `self.dict_coil_connection`
- `self.dict_kw_cjh`
- `self.dict_kw_els`
- `self.dict_suspension_kw_cjh`
- `self.dict_suspension_kw_els`
- `self.dict_suspension_kw_els['A_angle']`
- `self.dict_suspension_kw_els['B_angle']`
- `self.dict_suspension_kw_els['C_angle']`
- `self.dict_torque_kw_cjh`
- `self.dict_torque_kw_els`
- `self.distributed_or_concentrated`
- `self.dl_grouping_AC`
- `self.dl_grouping_BD`
- `self.dl_leftlayer`
- `self.dl_leftlayer['U']`
- `self.dl_leftlayer['V']`
- `self.dl_leftlayer['W']`
- `self.dl_rightlayer`
- `self.dl_rightlayer['U']`
- `self.dl_rightlayer['V']`
- `self.dl_rightlayer['W']`
- `self.dm`
- `self.doc`
- `self.doc.CreateReferenceFromItem`
- `self.doc.GetAssembly`
- `self.doc.GetSelection`
- `self.doc.SaveModel`
- `self.dpnv_grouping_dict_a`
- `self.dpnv_grouping_dict_a['GAC']`
- `self.dpnv_grouping_dict_a['GBD']`
- `self.dpnv_grouping_dict_b`
- `self.dpnv_grouping_dict_b['GAC']`
- `self.dpnv_grouping_dict_b['GBD']`
- `self.dpnv_grouping_dict_c`
- `self.dpnv_grouping_dict_c['GAC']`
- `self.dpnv_grouping_dict_c['GBD']`
- `self.draw_function`
- `self.draw_jmag`
- `self.drawer_T1`
- `self.drawer_T2`
- `self.drawer_T3a`
- `self.drawer_T3b`
- `self.drawer_T3c`
- `self.drawer_T4`
- `self.drawer_T4a`
- `self.drawer_T4b`
- `self.drawer_T4c`
- `self.drawer_Text`
- `self.drive_winding_conductors_per_slot`
- `self.drive_winding_current`
- `self.edge4Ref`
- `self.edge4ref`
- `self.estimated_back_emf`
- `self.excitation_frequency`
- `self.excitation_frequency_simulated`
- `self.f1`
- `self.f2`
- `self.f3`
- `self.fea_config_dict`
- `self.fea_config_dict['designer.StepPerCycle_3rdTSS']`
- `self.fea_config_dict['designer.max_nonlinear_iteration']`
- `self.fea_config_dict['designer.show']`
- `self.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
- `self.fea_config_dict['pc_name']`
- `self.femm_loss_list`
- `self.fig_plot2piFft`
- `self.fig_plotFuncObj`
- `self.filename`
- `self.fill_factor`
- `self.filter_data`
- `self.flag_do_not_evaluate_when_init_pop`
- `self.flag_material_already_loaded`
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
- `self.geomApp`
- `self.geomApp.GetDocument`
- `self.geomApp.NewDocument`
- `self.geometry`
- `self.geometry.all_points`
- `self.geometry.l_stack`
- `self.geometry.machineGeometry`
- `self.geometry.machineGeometry.items`
- `self.geometry.machineGeometry.values`
- `self.geometry.parts`
- `self.geometry.r_rotor_outer`
- `self.geometry.r_rotor_outer.value`
- `self.geometry.sync`
- `self.getSketch`
- `self.get_InitialRotationAngle`
- `self.get_certain_objective_function`
- `self.get_complex_number_kw`
- `self.get_complex_number_kw_per_phase`
- `self.get_complex_number_winding_factor_of_coil_i`
- `self.get_free_variables`
- `self.get_metric_of_the_whole_swarm`
- `self.get_parameter`
- `self.get_parameter_dict_by_name`
- `self.get_parameter_fields`
- `self.get_parameters_by_type`
- `self.get_parameters_summary`
- `self.get_rotor_volume`
- `self.get_voltage_and_current`
- `self.get_winding_factor`
- `self.gp`
- `self.gp['d_air_gap']`
- `self.gp['d_magnet']`
- `self.gp['d_tooth']`
- `self.gp['d_tooth_shoe']`
- `self.gp['r_rotor_outer']`
- `self.gp['r_shaft']`
- `self.gp['r_stator_outer']`
- `self.gp['tooth_shape']`
- `self.gp['w_tooth']`
- `self.grouping_AC`
- `self.hex_to_rgb`
- `self.horizontal_position`
- `self.iRotateCopy`
- `self.id`
- `self.id_rotorCore`
- `self.id_statorCore`
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
- `self.imag`
- `self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y`
- `self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y`
- `self.initial_excitation_bias_compensation_deg`
- `self.initial_rotation_angle`
- `self.initialized`
- `self.is_wye_connection`
- `self.jd`
- `self.jmag_loss_list`
- `self.k`
- `self.kd1`
- `self.kp1`
- `self.kw1`
- `self.l21`
- `self.l22`
- `self.l41`
- `self.l42`
- `self.l_FRW`
- `self.l_OA`
- `self.l_OB`
- `self.l_OC`
- `self.l_TRV`
- `self.l_design_parameters`
- `self.l_efficiency`
- `self.l_force_error_angle`
- `self.l_lamination`
- `self.l_leftlayer1`
- `self.l_leftlayer2`
- `self.l_normalized_force_error_magnitude`
- `self.l_normalized_torque_ripple`
- `self.l_original_rotor_weight`
- `self.l_original_stack_length`
- `self.l_power_factor`
- `self.l_rated_efficiency`
- `self.l_rated_iron_loss`
- `self.l_rated_magnet_Joule_loss`
- `self.l_rated_rotor_copper_loss_along_stack`
- `self.l_rated_rotor_volume`
- `self.l_rated_rotor_weight`
- `self.l_rated_shaft_power`
- `self.l_rated_stack_length`
- `self.l_rated_stator_copper_loss_along_stack`
- `self.l_rated_total_loss`
- `self.l_rated_windage_loss`
- `self.l_rightlayer1`
- `self.l_rightlayer2`
- `self.l_rotor_copper_loss_in_end_turn`
- `self.l_ss_avg_force_magnitude`
- `self.l_stack`
- `self.l_stator_copper_loss_in_end_turn`
- `self.l_torque_average`
- `self.lamination_count`
- `self.lamination_factor`
- `self.layer_A1`
- `self.layer_A2`
- `self.layer_B1`
- `self.layer_B2`
- `self.layer_X_phases`
- `self.layer_X_signs`
- `self.layer_Y_phases`
- `self.layer_Y_signs`
- `self.list_cost_function`
- `self.list_layer_motor_phases`
- `self.list_layer_motor_signs`
- `self.list_layer_suspension_phases`
- `self.list_layer_suspension_signs`
- `self.list_phase_u_slot_number`
- `self.list_phase_v_slot_number`
- `self.list_phase_w_slot_number`
- `self.list_rotor_current_amp`
- `self.list_slot_number_of_phase`
- `self.logger`
- `self.lst_x`
- `self.lst_y`
- `self.lst_y.index`
- `self.m`
- `self.machine_data`
- `self.machine_data.append`
- `self.magnet_area`
- `self.magnet_grade`
- `self.magnet_name`
- `self.magnet_start_angle`
- `self.magnet_temperature`
- `self.materials`
- `self.materials.stator_steel`
- `self.mec_power`
- `self.message`
- `self.mm_r_ro`
- `self.mm_r_ro.value`
- `self.mm_r_si`
- `self.mm_r_si.append`
- `self.mm_w_st`
- `self.mm_w_st.append`
- `self.model`
- `self.my_scatter_plot`
- `self.mycurrent`
- `self.mytime`
- `self.myvoltage`
- `self.name`
- `self.no_winding_layer`
- `self.normalized_force_error_magnitude`
- `self.number_of_chromosome`
- `self.number_of_designs`
- `self.number_of_free_variables`
- `self.number_of_parallel_branch`
- `self.number_of_winding_layer`
- `self.number_parallel_branch`
- `self.number_winding_layer`
- `self.open_jmag`
- `self.options`
- `self.overwritten`
- `self.ox_distribution_phase_U`
- `self.ox_distribution_three_phase`
- `self.p`
- `self.pairs`
- `self.parallel_branch_count`
- `self.parameter_dict`
- `self.parameter_dict_by_name`
- `self.parameter_dict_by_name.get`
- `self.parts`
- `self.parts.append`
- `self.path2SwarmData`
- `self.path2boptPython`
- `self.payload`
- `self.phase`
- `self.phase_count`
- `self.phase_current_amplitude`
- `self.phase_resistance`
- `self.pole_count`
- `self.pr`
- `self.prepareSection`
- `self.print_out_string`
- `self.projName`
- `self.project_names`
- `self.project_names.append`
- `self.ps`
- `self.q`
- `self.q2`
- `self.qs`
- `self.r_rotor_outer`
- `self.r_stator_outer`
- `self.radian_between_slots`
- `self.range_ss`
- `self.rated_current_density`
- `self.rated_data`
- `self.rated_data.append`
- `self.rated_speed`
- `self.read_csv_results_4_general_purpose`
- `self.real`
- `self.reference_data`
- `self.reference_data['2']`
- `self.reference_data['3']`
- `self.reference_data['4']`
- `self.reference_data['5']`
- `self.reference_data['6']`
- `self.reference_design`
- `self.reference_design['3']`
- `self.reference_design['3'].split`
- `self.regionCircularPattern360Origin`
- `self.regionMirrorCopy`
- `self.regions`
- `self.regions.append`
- `self.required_torque`
- `self.rotated_position`
- `self.rotor_core_material`
- `self.rotor_tooth_width_b_dr`
- `self.rotor_volume`
- `self.rotor_weight`
- `self.run_integer`
- `self.save`
- `self.scale`
- `self.scalingFactor`
- `self.select_FEA_tool`
- `self.select_fea_config_dict`
- `self.series_turns`
- `self.setSymPos`
- `self.setTurnFuncObject`
- `self.show`
- `self.show_geometry`
- `self.show_norm`
- `self.sine`
- `self.sketch`
- `self.sketch.CloseSketch`
- `self.sketch.CreateArc`
- `self.sketch.CreateBiConstraint`
- `self.sketch.CreateCircle`
- `self.sketch.CreateLine`
- `self.sketch.CreateRegionCircularPattern`
- `self.sketch.CreateRegionMirrorCopy`
- `self.sketch.CreateRegions`
- `self.sketch.GetItem`
- `self.sketch.OpenSketch`
- `self.sketch.SetProperty`
- `self.sketchNameList`
- `self.sketchNameList.append`
- `self.sketch_color`
- `self.slot_area`
- `self.slot_count`
- `self.slot_current_amplitude`
- `self.slot_current_at`
- `self.slot_per_phase`
- `self.spec`
- `self.spec.Jr`
- `self.spec.Js`
- `self.spec.Steel`
- `self.spec.VoltageRating`
- `self.spec.stator_phase_current_rms`
- `self.specs`
- `self.speed_rpm`
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
- `self.stack_length`
- `self.stack_length_max`
- `self.stack_length_specified`
- `self.stator_core_material`
- `self.stator_steel`
- `self.stator_tooth_width_b_ds`
- `self.steel_stack_factor`
- `self.str_best_design_details`
- `self.study`
- `self.study_name`
- `self.surface`
- `self.surface.finish`
- `self.suspen_kd_at_h`
- `self.suspen_kp_at_h`
- `self.suspen_kw_at_h`
- `self.suspension_current_ratio`
- `self.sw`
- `self.sw.fea_config_dict`
- `self.sw.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
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
- `self.swarm_data_as_dict`
- `self.swarm_data_as_dict.items`
- `self.swarm_data_project_names`
- `self.swarm_data_raw`
- `self.swarm_data_xf`
- `self.swarm_data_xf.append`
- `self.swarm_data_xf['0']`
- `self.sym_begin_pos`
- `self.sym_begin_pos_1`
- `self.sym_begin_pos_2`
- `self.sym_turn_func`
- `self.sym_winding_func`
- `self.sync`
- `self.t`
- `self.target`
- `self.target.counter`
- `self.target.fea_config_dict`
- `self.target.fea_config_dict['pc_name']`
- `self.target.free_parameters`
- `self.target.free_parameters.copy`
- `self.target.machine_class`
- `self.target.path2Data`
- `self.target.path2SwarmData`
- `self.target.project_name`
- `self.target.results_for_optimization`
- `self.target.select_FEA_tool`
- `self.target.swarm_data_json_file_path`
- `self.target.update_free_parameters`
- `self.template`
- `self.template.SI`
- `self.template.SI['GP']`
- `self.template.SI['GP']['mm_r_ro']`
- `self.template.SI['GP']['mm_r_ro'].value`
- `self.terminal_voltage`
- `self.time_list`
- `self.to_dict`
- `self.to_dict_full`
- `self.tooth_shape`
- `self.torque_average`
- `self.torque_current_ratio`
- `self.torque_current_utilization_ratio`
- `self.torque_kd_at_h`
- `self.torque_kp_at_h`
- `self.torque_kw_at_h`
- `self.ts`
- `self.turn_func`
- `self.turn_func_bias`
- `self.turns_per_slot`
- `self.type`
- `self.ui_info`
- `self.unit`
- `self.update_excitations`
- `self.user_input`
- `self.user_input['geometry']`
- `self.user_input['winding']`
- `self.user_input['winding']['l_stack']`
- `self.user_input['winding']['slot_count']`
- `self.value`
- `self.verbose`
- `self.verbose_drawing`
- `self.view`
- `self.visualization_points`
- `self.visualization_points.items`
- `self.weights_name`
- `self.weights_used`
- `self.wily`
- `self.wily.bool_DPNVorSEPA`
- `self.wily.deg_winding_U_phase_phase_axis_angle`
- `self.wily.kw1`
- `self.wily.number_of_parallel_branch`
- `self.winding`
- `self.winding.EX`
- `self.winding.EX.copy`
- `self.winding.EX.get`
- `self.winding.EX['mm_stack_length_specified']`
- `self.winding.sync`
- `self.winding.wily`
- `self.winding.wily.to_dict`
- `self.winding['l_stack']`
- `self.winding['pole_count']`
- `self.winding['slot_count']`
- `self.winding_dict`
- `self.winding_dict['bearing_winding_resistance']`
- `self.winding_dict['coil_pitch']`
- `self.winding_dict['connection_type']`
- `self.winding_dict['dc_bus_voltage']`
- `self.winding_dict['drive_winding_resistance']`
- `self.winding_dict['fill_factor']`
- `self.winding_dict['phase_count']`
- `self.winding_dict['pole_count']`
- `self.winding_dict['rated_current_density']`
- `self.winding_dict['rated_power']`
- `self.winding_dict['rated_speed']`
- `self.winding_dict['slot_count']`
- `self.winding_dict['suspension_current_ratio']`
- `self.winding_dict['torque_current_ratio']`
- `self.winding_dict['wire_diameter']`
- `self.winding_dict['wire_diameter_with_insulation']`
- `self.winding_func`
- `self.wire_diameter`
- `self.wire_diameter_with_insulation`
- `self.workDir`

### Variables Written (Outputs)
- `EX['end_winding_length_Lew']`
- `EX['stator_slot_area']`
- `self.Angle_StatorSlotOpen`
- `self.CommutatingSequenceB`
- `self.CommutatingSequenceD`
- `self.Cost`
- `self.Cost_Cu`
- `self.Cost_Fe`
- `self.Cost_PM`
- `self.Ea`
- `self.Em`
- `self.ExcitationFreqSimulated`
- `self.FRW`
- `self.ForConAbs_list`
- `self.ForConX_list`
- `self.ForConY_list`
- `self.GP`
- `self.HP`
- `self.HP_mirror`
- `self.JMAG_version_number`
- `self.JMAG_version_string`
- `self.Length_HeadNeckRotorSlot`
- `self.Omega`
- `self.PowerFactor`
- `self.Q`
- `self.Q_prime`
- `self.Qr`
- `self.Qs`
- `self.RP`
- `self.RatedEfficiency`
- `self.RatedStkLen`
- `self.RatedVol`
- `self.RatedWeight`
- `self.SI`
- `self.SIPNV_or_SEPA`
- `self.SIesign_parameters_denorm`
- `self.SIesign_parameters_norm`
- `self.SIir_run`
- `self.SPP`
- `self.TRV`
- `self.TorCon_list`
- `self.TorqueRipple`
- `self.Trip`
- `self.Tripple`
- `self.Width_StatorTeethHeadThickness`
- `self.accumSquaredData`
- `self.air_gap_length_delta`
- `self.ampl`
- `self.app`
- `self.ass`
- `self.avg_val_of_turn_func`
- `self.b1`
- `self.bFillRegion`
- `self.bMirror`
- `self.basic_info`
- `self.bearing_winding_conductors_per_slot`
- `self.bearing_winding_current`
- `self.best_design_denorm`
- `self.best_design_display`
- `self.bool_3PhaseCurrentSource`
- `self.bool_CustomizedCircuit`
- `self.bool_DPNVorSEPA`
- `self.bool_double_layer_winding`
- `self.bool_initialized`
- `self.bool_suppressShaft`
- `self.bounds`
- `self.buf`
- `self.buf_length`
- `self.calc`
- `self.calc_bounds`
- `self.coeff`
- `self.coil_pitch`
- `self.coil_pitch_y`
- `self.color`
- `self.components_make_region`
- `self.conductor_current_amplitude`
- `self.conductors_per_slot`
- `self.connection_star_raw_dict`
- `self.consts`
- `self.cosine`
- `self.count`
- `self.ctx`
- `self.defaultUnit`
- `self.deg_alpha_st`
- `self.deg_winding_U_phase_phase_axis_angle`
- `self.degree_between_slots`
- `self.dict_coil_connection`
- `self.dict_kw_cjh`
- `self.dict_kw_els`
- `self.distributed_or_concentrated`
- `self.dl_grouping_AC`
- `self.dl_grouping_BD`
- `self.dl_leftlayer`
- `self.dl_rightlayer`
- `self.dm`
- `self.doc`
- `self.dpnv_grouping_dict_a`
- `self.dpnv_grouping_dict_b`
- `self.dpnv_grouping_dict_c`
- `self.draw_function`
- `self.drawer_T1`
- `self.drawer_T2`
- `self.drawer_T3a`
- `self.drawer_T3b`
- `self.drawer_T3c`
- `self.drawer_T4`
- `self.drawer_T4a`
- `self.drawer_T4b`
- `self.drawer_T4c`
- `self.drawer_Text`
- `self.drive_winding_conductors_per_slot`
- `self.drive_winding_current`
- `self.edge4Ref`
- `self.estimated_back_emf`
- `self.excitation_frequency`
- `self.excitation_frequency_simulated`
- `self.f1`
- `self.f2`
- `self.f3`
- `self.fea_config_dict`
- `self.femm_loss_list`
- `self.fig_plot2piFft`
- `self.fig_plotFuncObj`
- `self.filename`
- `self.flag_do_not_evaluate_when_init_pop`
- `self.flag_material_already_loaded`
- `self.force_abs`
- `self.force_ang`
- `self.force_err_abs`
- `self.force_err_ang`
- `self.force_err_ang_new_way`
- `self.force_err_ang_old_way`
- `self.force_error_angle`
- `self.force_x`
- `self.force_y`
- `self.geomApp`
- `self.geometry`
- `self.geometry.all_points`
- `self.gp`
- `self.grouping_AC`
- `self.horizontal_position`
- `self.iRotateCopy`
- `self.id`
- `self.id_rotorCore`
- `self.id_statorCore`
- `self.imag`
- `self.initial_excitation_bias_compensation_deg`
- `self.initial_rotation_angle`
- `self.initialized`
- `self.jd`
- `self.jmag_loss_list`
- `self.k`
- `self.kd1`
- `self.kp1`
- `self.kw1`
- `self.l21`
- `self.l22`
- `self.l41`
- `self.l42`
- `self.l_FRW`
- `self.l_OA`
- `self.l_OB`
- `self.l_OC`
- `self.l_TRV`
- `self.l_design_parameters`
- `self.l_efficiency`
- `self.l_force_error_angle`
- `self.l_leftlayer1`
- `self.l_leftlayer2`
- `self.l_normalized_force_error_magnitude`
- `self.l_normalized_torque_ripple`
- `self.l_original_rotor_weight`
- `self.l_original_stack_length`
- `self.l_power_factor`
- `self.l_rated_efficiency`
- `self.l_rated_iron_loss`
- `self.l_rated_magnet_Joule_loss`
- `self.l_rated_rotor_copper_loss_along_stack`
- `self.l_rated_rotor_volume`
- `self.l_rated_rotor_weight`
- `self.l_rated_shaft_power`
- `self.l_rated_stack_length`
- `self.l_rated_stator_copper_loss_along_stack`
- `self.l_rated_total_loss`
- `self.l_rated_windage_loss`
- `self.l_rightlayer1`
- `self.l_rightlayer2`
- `self.l_rotor_copper_loss_in_end_turn`
- `self.l_ss_avg_force_magnitude`
- `self.l_stator_copper_loss_in_end_turn`
- `self.l_torque_average`
- `self.lamination_count`
- `self.lamination_factor`
- `self.layer_A1`
- `self.layer_A2`
- `self.layer_B1`
- `self.layer_B2`
- `self.layer_X_phases`
- `self.layer_X_signs`
- `self.layer_Y_phases`
- `self.layer_Y_signs`
- `self.list_layer_motor_phases`
- `self.list_layer_motor_signs`
- `self.list_layer_suspension_phases`
- `self.list_layer_suspension_signs`
- `self.list_phase_u_slot_number`
- `self.list_phase_v_slot_number`
- `self.list_phase_w_slot_number`
- `self.list_slot_number_of_phase`
- `self.logger`
- `self.lst_x`
- `self.lst_y`
- `self.m`
- `self.machine_data`
- `self.magnet_area`
- `self.magnet_name`
- `self.magnet_start_angle`
- `self.magnet_temperature`
- `self.materials`
- `self.mec_power`
- `self.message`
- `self.mm_r_si`
- `self.mm_w_st`
- `self.model`
- `self.mycurrent`
- `self.mytime`
- `self.myvoltage`
- `self.name`
- `self.no_winding_layer`
- `self.normalized_force_error_magnitude`
- `self.number_of_chromosome`
- `self.number_of_designs`
- `self.number_of_free_variables`
- `self.number_of_parallel_branch`
- `self.number_of_winding_layer`
- `self.number_parallel_branch`
- `self.number_winding_layer`
- `self.ox_distribution_phase_U`
- `self.ox_distribution_three_phase`
- `self.p`
- `self.pairs`
- `self.parameter_dict`
- `self.parameter_dict_by_name`
- `self.payload`
- `self.phase`
- `self.phase_current_amplitude`
- `self.phase_resistance`
- `self.pole_count`
- `self.pr`
- `self.print_out_string`
- `self.projName`
- `self.project_names`
- `self.ps`
- `self.q`
- `self.q2`
- `self.qs`
- `self.radian_between_slots`
- `self.range_ss`
- `self.rated_data`
- `self.real`
- `self.reference_data`
- `self.reference_design`
- `self.regions`
- `self.required_torque`
- `self.rotated_position`
- `self.rotor_core_material`
- `self.rotor_tooth_width_b_dr`
- `self.rotor_volume`
- `self.rotor_weight`
- `self.run_integer`
- `self.scale`
- `self.scalingFactor`
- `self.series_turns`
- `self.sine`
- `self.sketch`
- `self.sketchNameList`
- `self.sketch_color`
- `self.slot_area`
- `self.slot_count`
- `self.slot_current_amplitude`
- `self.slot_current_at`
- `self.slot_per_phase`
- `self.spec`
- `self.specs`
- `self.speed_rpm`
- `self.ss_avg_force_angle`
- `self.ss_avg_force_magnitude`
- `self.ss_avg_force_vector`
- `self.ss_max_force_err_abs`
- `self.ss_max_force_err_ang`
- `self.stack_length`
- `self.stack_length_max`
- `self.stack_length_specified`
- `self.stator_core_material`
- `self.stator_tooth_width_b_ds`
- `self.str_best_design_details`
- `self.study`
- `self.study_name`
- `self.surface`
- `self.suspen_kd_at_h`
- `self.suspen_kp_at_h`
- `self.suspen_kw_at_h`
- `self.sw`
- `self.swarm_data_as_dict`
- `self.swarm_data_project_names`
- `self.swarm_data_raw`
- `self.swarm_data_xf`
- `self.sym_begin_pos`
- `self.sym_begin_pos_1`
- `self.sym_begin_pos_2`
- `self.sym_turn_func`
- `self.sym_winding_func`
- `self.t`
- `self.target`
- `self.target.counter`
- `self.target.fea_config_dict`
- `self.target.path2Data`
- `self.target.path2SwarmData`
- `self.target.project_name`
- `self.target.results_for_optimization`
- `self.target.swarm_data_json_file_path`
- `self.template.SI['GP']['mm_r_ro'].value`
- `self.time_list`
- `self.torque_average`
- `self.torque_current_utilization_ratio`
- `self.torque_kd_at_h`
- `self.torque_kp_at_h`
- `self.torque_kw_at_h`
- `self.ts`
- `self.turn_func`
- `self.turn_func_bias`
- `self.turns_per_slot`
- `self.type`
- `self.ui_info`
- `self.unit`
- `self.user_input`
- `self.value`
- `self.verbose`
- `self.verbose_drawing`
- `self.view`
- `self.visualization_points`
- `self.weights_name`
- `self.weights_used`
- `self.wily`
- `self.winding`
- `self.winding_func`
- `self.workDir`

## File: `clean_jfiles.py`

### `get_swarm_group` (Line 22)
- **Arguments**: self, folder_of_collection
- **Reads (State Dependencies)**:
  - `self.path2boptPython`
- **Writes**: None state variables detected

## File: `clean_large_files.py`

### `get_swarm_group` (Line 22)
- **Arguments**: self, folder_of_collection
- **Reads (State Dependencies)**:
  - `self.path2boptPython`
- **Writes**: None state variables detected

## File: `app\config_loader.py`

### `get_project_root` (Line 26)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `load_acmop_config` (Line 31)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_backend_virtual_env` (Line 56)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_frontend_backend_url` (Line 62)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `app\routers\acmop.py`

### `get_codes4_path` (Line 33)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_default_dir` (Line 40)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `convert_numpy_types` (Line 1076)
- **Arguments**: obj
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `simple_non_dominated_sorting` (Line 1101)
- **Arguments**: fits
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `app\routers\debug.py`

### `get_default_user_input` (Line 9)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `app\routers\design.py`

### `get_codes4_path` (Line 8)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `app\routers\machine_specs.py`

### `_agent_log` (Line 19)
- **Arguments**: payload
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `_robust_dict` (Line 32)
- **Arguments**: obj
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `_specs_to_api_response` (Line 67)
- **Arguments**: my_machine
- **Reads (State Dependencies)**:
  - `my_machine.winding`
  - `my_machine.target`
  - `my_machine.materials`
  - `my_machine.geometry`
- **Writes**: None state variables detected

## File: `app\routers\project.py`

### `load_specifications` (Line 10)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `app\services\optimization_service.py`

### `OptimizationService.__init__` (Line 15)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.logger`
- **Writes (State Changes)**:
  - `self.logger`

### `OptimizationService.run_optimization` (Line 18)
- **Arguments**: self, spec_name, fea_config, project_loc
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `OptimizationService.run_winding_part` (Line 34)
- **Arguments**: self, spec_name, fea_config, project_loc
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `BH\Inspect_BH_curve.py`

### `script_dir` (Line 21)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `is_bh_curve_file` (Line 25)
- **Arguments**: path
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `read_bh_data` (Line 31)
- **Arguments**: path
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `write_bh_file` (Line 60)
- **Arguments**: path, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `ensure_bh_file` (Line 68)
- **Arguments**: source
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `collect_and_convert_bh_curves` (Line 83)
- **Arguments**: dir_path
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `plot_bh_curves` (Line 111)
- **Arguments**: labels_and_data, out_path
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `main` (Line 141)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\angle_error_nick.py`

### `angle_error` (Line 2)
- **Arguments**: alpha_star, alpha_actual
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `compute_angle_error` (Line 63)
- **Arguments**: alpha_star, alpha_actual
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `main` (Line 90)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\fix_jmag.py`

### `repl` (Line 40)
- **Arguments**: m
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `repl_EX` (Line 46)
- **Arguments**: m
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\JMAG.py`

### `JMAG.__init__` (Line 31)
- **Arguments**: self, fea_config_dict
- **Reads (State Dependencies)**:
  - `self.app`
  - `self.fea_config_dict`
  - `self.consts`
  - `self.bMirror`
  - `self.ass`
  - `self.flag_material_already_loaded`
  - `self.sketchNameList`
  - `self.geomApp`
  - `self.model`
  - `self.projName`
  - `self.jd`
  - `self.bool_suppressShaft`
  - `self.workDir`
  - `self.view`
  - `self.doc`
  - `self.JMAG_version_number`
  - `self.verbose_drawing`
  - `self.sketch`
  - `self.study`
  - `self.defaultUnit`
  - `self.iRotateCopy`
  - `self.edge4Ref`
- **Writes (State Changes)**:
  - `self.app`
  - `self.fea_config_dict`
  - `self.consts`
  - `self.bMirror`
  - `self.ass`
  - `self.flag_material_already_loaded`
  - `self.sketchNameList`
  - `self.geomApp`
  - `self.model`
  - `self.projName`
  - `self.jd`
  - `self.bool_suppressShaft`
  - `self.workDir`
  - `self.view`
  - `self.doc`
  - `self.JMAG_version_number`
  - `self.verbose_drawing`
  - `self.sketch`
  - `self.study`
  - `self.defaultUnit`
  - `self.iRotateCopy`
  - `self.edge4Ref`

### `JMAG.open` (Line 59)
- **Arguments**: self, Steel_name, expected_project_file_path, pc_name, dir_parent, bool_jmagDesignerShow
- **Reads (State Dependencies)**:
  - `self.app`
  - `self.fea_config_dict`
  - `self.JMAG_version_string`
  - `self.fea_config_dict['pc_name']`
  - `self.JMAG_version_number`
  - `self.flag_material_already_loaded`
- **Writes (State Changes)**:
  - `self.app`
  - `self.flag_material_already_loaded`
  - `self.JMAG_version_string`
  - `self.JMAG_version_number`

### `JMAG.close` (Line 269)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.app`
  - `self.app.Quit`
- **Writes**: None state variables detected

### `JMAG.save` (Line 272)
- **Arguments**: self, name, description
- **Reads (State Dependencies)**:
  - `self.doc`
  - `self.app`
  - `self.doc.SaveModel`
  - `self.app.GetCurrentModel`
- **Writes**: None state variables detected

### `JMAG.pre_process_PMSM` (Line 280)
- **Arguments**: self, app, model, acm_variant
- **Reads (State Dependencies)**:
  - `self.doc.GetSelection`
  - `self.fea_config_dict`
  - `self.doc`
  - `self.id_rotorCore`
  - `self.fea_config_dict['designer.show']`
  - `self.id_statorCore`
- **Writes (State Changes)**:
  - `self.id_rotorCore`
  - `self.id_statorCore`

### `JMAG.add_magnetic_transient_study` (Line 456)
- **Arguments**: self, app, model, path2FEACsv, study_name, acm_variant
- **Reads (State Dependencies)**:
  - `self.fea_config_dict`
  - `self.add_circuit`
  - `self.study_name`
  - `self.add_material`
  - `self.JMAG_version_number`
  - `self.id_rotorCore`
  - `self.fea_config_dict['designer.max_nonlinear_iteration']`
  - `self.id_statorCore`
- **Writes (State Changes)**:
  - `self.study_name`

### `JMAG.add_structural_static_study` (Line 746)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.add_mesh` (Line 748)
- **Arguments**: self, study, model
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.add_material` (Line 751)
- **Arguments**: self, study, acm_variant
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.add_circuit` (Line 795)
- **Arguments**: self, app, model, study, acm_variant, bool_3PhaseCurrentSource
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.addConstraintCocentricity` (Line 1040)
- **Arguments**: self, vA, vB
- **Reads (State Dependencies)**:
  - `self.sketch.CreateBiConstraint`
  - `self.doc`
  - `self.doc.CreateReferenceFromItem`
  - `self.sketch.GetItem`
  - `self.sketch`
- **Writes**: None state variables detected

### `JMAG.drawLine` (Line 1060)
- **Arguments**: self, startxy, endxy, returnVertexName
- **Reads (State Dependencies)**:
  - `self.getSketch`
  - `self.sketch.OpenSketch`
  - `self.sketch.CreateLine`
  - `self.sketch`
- **Writes (State Changes)**:
  - `self.sketch`

### `JMAG.drawArc` (Line 1078)
- **Arguments**: self, centerxy, startxy, endxy, returnVertexName
- **Reads (State Dependencies)**:
  - `self.getSketch`
  - `self.sketch.CreateArc`
  - `self.sketch.OpenSketch`
  - `self.sketch`
- **Writes (State Changes)**:
  - `self.sketch`

### `JMAG.drawCircle` (Line 1095)
- **Arguments**: self, centerxy, radius, returnVertexName
- **Reads (State Dependencies)**:
  - `self.sketch.CreateCircle`
  - `self.getSketch`
  - `self.sketch.OpenSketch`
  - `self.sketch`
- **Writes (State Changes)**:
  - `self.sketch`

### `JMAG.checkGeomApp` (Line 1109)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.geomApp`
  - `self.app`
  - `self.app.CreateGeometryEditor`
  - `self.doc`
  - `self.geomApp.NewDocument`
  - `self.app.LaunchGeometryEditor`
- **Writes (State Changes)**:
  - `self.doc`
  - `self.geomApp`

### `JMAG.getSketch` (Line 1117)
- **Arguments**: self, sketchName, color
- **Reads (State Dependencies)**:
  - `self.geomApp`
  - `self.sketch.SetProperty`
  - `self.doc.GetAssembly`
  - `self.ass.CreateSketch`
  - `self.doc.CreateReferenceFromItem`
  - `self.doc`
  - `self.sketchNameList.append`
  - `self.sketch.OpenSketch`
  - `self.ass`
  - `self.sketch`
  - `self.sketchNameList`
  - `self.ass.GetItem`
  - `self.geomApp.GetDocument`
  - `self.checkGeomApp`
- **Writes (State Changes)**:
  - `self.doc`
  - `self.ass`
  - `self.geomApp`
  - `self.sketch`

### `JMAG.prepareSection` (Line 1142)
- **Arguments**: self, token, bMirrorMerge, bRotateMerge
- **Reads (State Dependencies)**:
  - `self.doc.GetSelection`
  - `self.edge4ref`
  - `self.regionCircularPattern360Origin`
  - `self.doc`
  - `self.bMirror`
  - `self.regionMirrorCopy`
  - `self.sketch.CloseSketch`
  - `self.iRotateCopy`
  - `self.sketch.CreateRegions`
  - `self.sketch.GetItem`
  - `self.edge4Ref`
  - `self.sketch`
- **Writes**: None state variables detected

### `JMAG.regionMirrorCopy` (Line 1197)
- **Arguments**: self, region, edge4Ref, symmetryType, bMerge
- **Reads (State Dependencies)**:
  - `self.sketch.CreateRegionMirrorCopy`
  - `self.doc.CreateReferenceFromItem`
  - `self.doc`
  - `self.sketch.GetItem`
  - `self.ass`
  - `self.sketch`
  - `self.ass.GetItem`
- **Writes**: None state variables detected

### `JMAG.regionCircularPattern360Origin` (Line 1219)
- **Arguments**: self, idx, region, Q_float, bMerge
- **Reads (State Dependencies)**:
  - `self.doc`
  - `self.sketch`
  - `self.sketch.CreateRegionCircularPattern`
  - `self.doc.CreateReferenceFromItem`
- **Writes**: None state variables detected

### `JMAG.draw_jmag_model` (Line 1265)
- **Arguments**: self, app, individual_index, im_variant, model_name, bool_trimDrawer_or_vanGogh, doNotRotateCopy
- **Reads (State Dependencies)**:
  - `self.SI`
- **Writes (State Changes)**:
  - `self.SI`

### `JMAG.run_study` (Line 1335)
- **Arguments**: acm_variant, app, study, fea_config_dict, toc
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.mesh_study` (Line 1360)
- **Arguments**: self, acm_variant, app, model, study, output_dir
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.draw_spmsm` (Line 1444)
- **Arguments**: self, acm_variant, bool_pyx
- **Reads (State Dependencies)**:
  - `self.bMirror`
  - `self.iRotateCopy`
  - `self.prepareSection`
  - `self.save`
  - `self.show`
  - `self.calculate_excitation_current`
  - `acm_variant.rotorMagnet.notched_rotor.p`
- **Writes (State Changes)**:
  - `self.bMirror`
  - `self.iRotateCopy`

### `JMAG.show` (Line 1585)
- **Arguments**: self, acm_variant, toString
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.add_plots` (Line 1610)
- **Arguments**: axeses, dm, title, label, zorder, time_list, sfv, torque, range_ss, alpha
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `JMAG.read_csv_results_4_general_purpose` (Line 1669)
- **Arguments**: study_name, path_prefix, fea_config_dict, femm_solver, acm_variant
- **Reads (State Dependencies)**:
  - `self.mycurrent`
  - `self.myvoltage`
  - `self.terminal_voltage`
  - `self.time_list`
  - `self.Current_dict`
  - `self.mytime`
  - `self.DisplacementAngle_list`
  - `self.get_voltage_and_current`
  - `self.ui_info`
  - `self.coil_fluxLinkage`
  - `self.jmag_loss_list`
  - `self.basic_info`
  - `self.ForConX_list`
  - `self.femm_loss_list`
  - `self.ForConY_list`
  - `self.TorCon_list`
  - `self.FluxLinkage_dict`
  - `self.ForConAbs_list`
  - `self.circuit_current`
  - `self.Current_dict['Time(s)']`
- **Writes (State Changes)**:
  - `self.mycurrent`
  - `self.mytime`
  - `self.myvoltage`
  - `self.time_list`
  - `self.ui_info`
  - `self.ForConAbs_list`
  - `self.ForConX_list`
  - `self.femm_loss_list`
  - `self.jmag_loss_list`
  - `self.basic_info`
  - `self.ForConY_list`
  - `self.TorCon_list`

### `JMAG.build_str_results` (Line 2097)
- **Arguments**: self, acm_variant, project_name, tran_study_name, path2FEACsv, fea_config_dict, femm_solver
- **Reads (State Dependencies)**:
  - `self.fea_config_dict['designer.StepPerCycle_3rdTSS']`
  - `self.read_csv_results_4_general_purpose`
  - `self.fea_config_dict`
  - `self.dm`
  - `self.add_plots`
- **Writes (State Changes)**:
  - `self.dm`

## File: `codes4\legacy_codes.py`

### `LegacyCodes.PracticalInitialDesign` (Line 6)
- **Arguments**: self, fea_config_dict, SI, GP, EX
- **Reads (State Dependencies)**:
  - `EX['end_winding_length_factor_kov']`
  - `EX['end_winding_length_Lew']`
  - `EX['stator_slot_area']`
- **Writes (State Changes)**:
  - `EX['end_winding_length_Lew']`
  - `EX['stator_slot_area']`

## File: `codes4\machine.py`

### `Machine.__init__` (Line 23)
- **Arguments**: self, user_input
- **Reads (State Dependencies)**:
  - `self.materials`
  - `self.geometry`
  - `self.target`
  - `self.user_input`
  - `self.winding`
- **Writes (State Changes)**:
  - `self.materials`
  - `self.geometry`
  - `self.target`
  - `self.user_input`
  - `self.winding`

### `Machine.get_rotor_volume` (Line 33)
- **Arguments**: self, stack_length
- **Reads (State Dependencies)**:
  - `self.winding.EX['mm_stack_length_specified']`
  - `self.geometry`
  - `self.winding.EX.get`
  - `self.winding.EX`
  - `self.geometry.r_rotor_outer`
  - `self.winding`
  - `self.geometry.r_rotor_outer.value`
- **Writes**: None state variables detected

### `Machine.get_rotor_weight` (Line 39)
- **Arguments**: self, gravity, stack_length
- **Reads (State Dependencies)**:
  - `self.get_rotor_volume`
- **Writes**: None state variables detected

### `Machine.get_free_variables_as_dict` (Line 43)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.target`
  - `self.target.free_parameters`
- **Writes**: None state variables detected

### `Machine.sync` (Line 46)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.user_input['geometry']`
  - `self.materials`
  - `self.geometry`
  - `self.geometry.parts`
  - `self.user_input['winding']`
  - `self.user_input['winding']['l_stack']`
  - `self.winding.sync`
  - `self.user_input`
  - `self.user_input['winding']['slot_count']`
  - `self.winding`
  - `self.geometry.sync`
  - `self.geometry.l_stack`
  - `self.geometry.all_points`
- **Writes (State Changes)**:
  - `self.geometry.all_points`

### `Machine.draw_machine_using_CairoDrawer` (Line 74)
- **Arguments**: self, drawer
- **Reads (State Dependencies)**:
  - `self.geometry.parts`
  - `self.geometry.all_points`
  - `self.geometry`
- **Writes**: None state variables detected

### `Machine.draw_machine_using_JMAG` (Line 88)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.draw_jmag`
  - `self.materials`
  - `self.materials.stator_steel`
  - `self.open_jmag`
- **Writes**: None state variables detected

### `Machine.open_jmag` (Line 99)
- **Arguments**: self, expected_project_file, Steel_name, bool_jmagDesignerShow
- **Reads (State Dependencies)**:
  - `self.target`
  - `self.target.fea_config_dict['pc_name']`
  - `self.target.fea_config_dict`
  - `self.target.project_name`
- **Writes (State Changes)**:
  - `self.target.fea_config_dict`
  - `self.target.project_name`

### `Machine.draw_jmag` (Line 111)
- **Arguments**: self, toolJd
- **Reads (State Dependencies)**:
  - `self.geometry`
  - `self.target.project_name`
  - `self.geometry.parts`
  - `self.target`
  - `self.geometry.all_points`
- **Writes**: None state variables detected

### `Machine.FEA_evaluate` (Line 137)
- **Arguments**: self, project_loc, bool_jmagDesignerShow, x_denorm, counter, counter_loop
- **Reads (State Dependencies)**:
  - `self.materials`
  - `self.target.swarm_data_json_file_path`
  - `self.target.project_name`
  - `self.sync`
  - `self.target`
  - `self.target.path2Data`
  - `self.target.path2SwarmData`
  - `self.target.update_free_parameters`
  - `self.open_jmag`
  - `self.target.select_FEA_tool`
  - `self.draw_jmag`
  - `self.materials.stator_steel`
  - `self.target.counter`
- **Writes (State Changes)**:
  - `self.target.swarm_data_json_file_path`
  - `self.target.project_name`
  - `self.target.path2Data`
  - `self.target.path2SwarmData`
  - `self.target.counter`

### `Machine.compile_results` (Line 193)
- **Arguments**: self, toolJd, study_name, path2FEACsv, counter
- **Reads (State Dependencies)**:
  - `self.target.swarm_data_json_file_path`
  - `self.target.project_name`
  - `self.target.results_for_optimization`
  - `self.target.free_parameters.copy`
  - `self.target`
  - `self.target.fea_config_dict`
  - `self.target.select_FEA_tool`
  - `self.target.free_parameters`
- **Writes (State Changes)**:
  - `self.target.results_for_optimization`

### `run_step_by_step` (Line 260)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\machine_geometry.py`

### `rotate_point` (Line 9)
- **Arguments**: p, deg
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `AllPoints.__post_init__` (Line 21)
- **Arguments**: self, user_input
- **Reads (State Dependencies)**:
  - `self.HP['9']`
  - `self.HP['4']['0']`
  - `self.HP['7']`
  - `self.HP['7']['0']`
  - `self.HP['6']['1']`
  - `self.HP['9']['1']`
  - `self.HP['7']['1']`
  - `self.HP`
  - `self.HP['9']['0']`
  - `self.rotated_position`
  - `self.HP['1']['1']`
  - `self.HP['4']`
  - `self.HP['3']['1']`
  - `self.slot_count`
  - `self.horizontal_position`
  - `self.HP['6']['0']`
  - `self.GP`
  - `self.HP['6']`
  - `self.HP['3']['0']`
  - `self.pole_count`
  - `self.HP_mirror`
  - `self.HP['2']['0']`
  - `self.HP['4']['1']`
  - `self.HP['1']['0']`
  - `self.HP['1']`
  - `self.RP`
  - `self.HP['2']['1']`
  - `self.HP['3']`
  - `self.HP['2']`
- **Writes (State Changes)**:
  - `self.HP`
  - `self.GP`
  - `self.rotated_position`
  - `self.RP`
  - `self.pole_count`
  - `self.HP_mirror`
  - `self.slot_count`
  - `self.horizontal_position`

### `parse_point_name` (Line 103)
- **Arguments**: name, all_points, rotation_deg
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `draw_instruction_parser` (Line 138)
- **Arguments**: part, all_points, drawer
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `RotorCore.draw_instruction` (Line 251)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.options`
  - `self.all_points.pole_count`
  - `self.all_points`
- **Writes**: None state variables detected

### `StatorCore.draw_instruction` (Line 269)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points.slot_count`
  - `self.options`
  - `self.all_points`
- **Writes**: None state variables detected

### `Magnet.draw_instruction` (Line 305)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points.pole_count`
  - `self.all_points`
- **Writes**: None state variables detected

### `Coil.draw_instruction` (Line 325)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points.slot_count`
  - `self.all_points`
- **Writes**: None state variables detected

### `MachineGeometry.slot_count` (Line 360)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding`
  - `self.winding['slot_count']`
- **Writes**: None state variables detected

### `MachineGeometry.pole_count` (Line 364)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding`
  - `self.winding['pole_count']`
- **Writes**: None state variables detected

### `MachineGeometry.l_stack` (Line 368)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding['l_stack']`
  - `self.winding`
- **Writes**: None state variables detected

### `MachineGeometry.d_air_gap` (Line 372)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['d_air_gap']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.r_stator_outer` (Line 376)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['r_stator_outer']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.r_rotor_outer` (Line 380)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['r_rotor_outer']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.r_shaft` (Line 384)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['r_shaft']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.w_tooth` (Line 388)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['w_tooth']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.d_tooth` (Line 392)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['d_tooth']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.d_magnet` (Line 396)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['d_magnet']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.tooth_shape` (Line 400)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp['tooth_shape']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.deg_alpha_rm` (Line 404)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.pole_count`
- **Writes**: None state variables detected

### `MachineGeometry.split_ratio` (Line 409)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.r_rotor_outer`
  - `self.d_air_gap`
  - `self.r_stator_outer`
- **Writes**: None state variables detected

### `MachineGeometry.d_stator_yoke` (Line 413)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.r_rotor_outer`
  - `self.r_stator_outer`
  - `self.d_tooth`
- **Writes**: None state variables detected

### `MachineGeometry.d_tooth_shoe` (Line 417)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.tooth_shape`
  - `self.gp['d_tooth_shoe']`
  - `self.gp`
- **Writes**: None state variables detected

### `MachineGeometry.alpha_stator_tooth_span` (Line 423)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.tooth_shape`
  - `self.slot_count`
- **Writes**: None state variables detected

### `MachineGeometry.aspect_ratio` (Line 429)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.l_stack`
  - `self.r_stator_outer`
- **Writes**: None state variables detected

### `MachineGeometry.__post_init__` (Line 432)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `MachineGeometry.add_part` (Line 436)
- **Arguments**: self, part
- **Reads (State Dependencies)**:
  - `self.parts.append`
  - `self._next_index`
  - `self.parts`
  - `self.all_points`
- **Writes**: None state variables detected

### `MachineGeometry.add_radial_array` (Line 443)
- **Arguments**: self, part_type, base_name, count, options, color
- **Reads (State Dependencies)**:
  - `self.add_part`
- **Writes**: None state variables detected

### `MachineGeometry.show_geometry_svg` (Line 453)
- **Arguments**: self, filename, scale
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `MachineGeometry.sync` (Line 459)
- **Arguments**: self, user_input
- **Reads (State Dependencies)**:
  - `self.winding`
  - `self.gp`
- **Writes (State Changes)**:
  - `self.winding`
  - `self.gp`

## File: `codes4\machine_materials.py`

### `MachineMaterial.__post_init__` (Line 20)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.lamination_count`
  - `self.l_lamination`
  - `self.steel_stack_factor`
  - `self.magnet_grade`
  - `self.l_stack`
  - `self.stator_steel`
  - `self.d_magnet`
- **Writes (State Changes)**:
  - `self.lamination_count`

## File: `codes4\machine_target.py`

### `MachineTarget.update_free_parameters` (Line 43)
- **Arguments**: self, user_input, x
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `MachineTarget.get_required_GP` (Line 50)
- **Arguments**: self, user_input
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\machine_winding.py`

### `MachineWinding.phase_count` (Line 13)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['phase_count']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.slot_count` (Line 15)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['slot_count']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.pole_count` (Line 17)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['pole_count']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.coil_pitch` (Line 19)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['coil_pitch']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.Qs` (Line 23)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.slot_count`
- **Writes**: None state variables detected

### `MachineWinding.p` (Line 25)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.pole_count`
- **Writes**: None state variables detected

### `MachineWinding.wire_diameter_with_insulation` (Line 31)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict`
  - `self.winding_dict['wire_diameter_with_insulation']`
- **Writes**: None state variables detected

### `MachineWinding.wire_diameter` (Line 33)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['wire_diameter']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.rated_current_density` (Line 35)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['rated_current_density']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.parallel_branch_count` (Line 37)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.wily`
  - `self.wily.number_of_parallel_branch`
- **Writes**: None state variables detected

### `MachineWinding.connection_type` (Line 41)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['connection_type']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.is_wye_connection` (Line 43)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.connection_type`
- **Writes**: None state variables detected

### `MachineWinding.rated_speed` (Line 45)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['rated_speed']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.rated_power` (Line 47)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict`
  - `self.winding_dict['rated_power']`
- **Writes**: None state variables detected

### `MachineWinding.dc_bus_voltage` (Line 49)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['dc_bus_voltage']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.fill_factor` (Line 51)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['fill_factor']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.torque_current_ratio` (Line 60)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['torque_current_ratio']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.suspension_current_ratio` (Line 62)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict`
  - `self.winding_dict['suspension_current_ratio']`
- **Writes**: None state variables detected

### `MachineWinding.drive_winding_resistance` (Line 64)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['drive_winding_resistance']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.bearing_winding_resistance` (Line 66)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.winding_dict['bearing_winding_resistance']`
  - `self.winding_dict`
- **Writes**: None state variables detected

### `MachineWinding.__post_init__` (Line 92)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.bool_DPNVorSEPA`
  - `self.slot_count`
  - `self.wily`
  - `self.pole_count`
  - `self.phase_count`
  - `self.ps`
  - `self.coil_pitch`
- **Writes (State Changes)**:
  - `self.wily`
  - `self.ps`

### `MachineWinding.get_InitialRotationAngle` (Line 114)
- **Arguments**: self, machineGeometry
- **Reads (State Dependencies)**:
  - `self.wily.deg_winding_U_phase_phase_axis_angle`
  - `self.p`
  - `self.wily`
  - `self.initial_rotation_angle`
- **Writes (State Changes)**:
  - `self.initial_rotation_angle`

### `MachineWinding.update_excitations` (Line 119)
- **Arguments**: self, machineGeometry, materials
- **Reads (State Dependencies)**:
  - `self.torque_current_ratio`
  - `self.phase_current_amplitude`
  - `self.drive_winding_conductors_per_slot`
  - `self.rated_current_density`
  - `self.wily`
  - `self.bearing_winding_conductors_per_slot`
  - `self.slot_current_amplitude`
  - `self.phase_count`
  - `self.slot_area`
  - `self.stack_length_specified`
  - `self.torque_current_utilization_ratio`
  - `self.wily.kw1`
  - `self.wily.number_of_parallel_branch`
  - `self.slot_count`
  - `self.magnet_area`
  - `self.excitation_frequency_simulated`
  - `self.dc_bus_voltage`
  - `self.get_InitialRotationAngle`
  - `self.is_wye_connection`
  - `self.fill_factor`
  - `self.wily.bool_DPNVorSEPA`
  - `self.series_turns`
  - `self.bearing_winding_current`
  - `self.rated_speed`
  - `self.suspension_current_ratio`
  - `self.p`
  - `self.initial_rotation_angle`
  - `self.conductors_per_slot`
  - `self.drive_winding_current`
  - `self.conductor_current_amplitude`
- **Writes (State Changes)**:
  - `self.bearing_winding_conductors_per_slot`
  - `self.slot_area`
  - `self.slot_current_amplitude`
  - `self.drive_winding_current`
  - `self.phase_current_amplitude`
  - `self.drive_winding_conductors_per_slot`
  - `self.initial_rotation_angle`
  - `self.torque_current_utilization_ratio`
  - `self.magnet_area`
  - `self.conductors_per_slot`
  - `self.series_turns`
  - `self.excitation_frequency_simulated`
  - `self.bearing_winding_current`
  - `self.conductor_current_amplitude`

### `MachineWinding.sync` (Line 179)
- **Arguments**: self, machineGeometry, materials, l_stack
- **Reads (State Dependencies)**:
  - `self.estimated_back_emf`
  - `self.phase_resistance`
  - `self.bool_DPNVorSEPA`
  - `self.stator_core_material`
  - `self.magnet_temperature`
  - `self.rated_current_density`
  - `self.wily`
  - `self.wire_diameter`
  - `self.stack_length_specified`
  - `self.phase_count`
  - `self.excitation_frequency`
  - `self.wily.number_of_parallel_branch`
  - `self.lamination_factor`
  - `self.magnet_start_angle`
  - `self.slot_count`
  - `self.coil_pitch`
  - `self.rotor_core_material`
  - `self.update_excitations`
  - `self.pole_count`
  - `self.ps`
  - `self.magnet_name`
  - `self.parallel_branch_count`
  - `self.rated_speed`
  - `self.p`
  - `self.conductors_per_slot`
  - `self.slot_current_at`
- **Writes (State Changes)**:
  - `self.magnet_name`
  - `self.estimated_back_emf`
  - `self.phase_resistance`
  - `self.stator_core_material`
  - `self.magnet_temperature`
  - `self.lamination_factor`
  - `self.rotor_core_material`
  - `self.wily`
  - `self.stack_length_specified`
  - `self.slot_current_at`
  - `self.magnet_start_angle`
  - `self.excitation_frequency`

### `MachineWinding.validate_fill_factor` (Line 242)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.wire_diameter_with_insulation`
  - `self.rated_current_density`
  - `self.wire_diameter`
  - `self.fill_factor`
- **Writes**: None state variables detected

## File: `codes4\modern_machine_designer_utility.py`

### `Parameter.__init__` (Line 6)
- **Arguments**: self, name, type, value, bounds, calc, calc_bounds, unit, parameter_dict
- **Reads (State Dependencies)**:
  - `self.type`
  - `self.name`
  - `self.parameter_dict`
  - `self.unit`
  - `self.bounds['1']`
  - `self.bounds['0']`
  - `self.calc`
  - `self.overwritten`
  - `self.value`
  - `self.calc_bounds`
  - `self.initialized`
  - `self.bounds`
- **Writes (State Changes)**:
  - `self.type`
  - `self.name`
  - `self.parameter_dict`
  - `self.unit`
  - `self.calc`
  - `self.value`
  - `self.initialized`
  - `self.calc_bounds`
  - `self.bounds`

### `Parameter.__repr__` (Line 43)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.type`
  - `self.unit`
  - `self.name`
  - `self.value`
- **Writes**: None state variables detected

### `Parameter.sensitivity` (Line 46)
- **Arguments**: self, param_name
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Parameter.to_dict` (Line 49)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.type`
  - `self.name`
  - `self.unit`
  - `self.value`
  - `self.bounds`
- **Writes**: None state variables detected

### `Parameter.from_dict` (Line 68)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Winding.__init__` (Line 92)
- **Arguments**: self, phase_number_m, stator_slot_number_Qs, pole_pair_number_p, suspension_pole_pair_number_ps, coil_pitch_y, bool_DPNVorSEPA, number_of_parallel_branch
- **Reads (State Dependencies)**:
  - `self.bool_distributed_or_concentrated`
  - `self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y`
  - `self.bool_DPNVorSEPA`
  - `self.bool_CustomizedCircuit`
  - `self.bool_3PhaseCurrentSource`
  - `self.layer_X_phases`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.Qs`
  - `self.CommutatingSequenceD`
  - `self.m`
  - `self.layer_X_signs`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.grouping_AC`
  - `self.layer_Y_phases`
  - `self.kw1`
  - `self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y`
  - `self.get_winding_factor`
  - `self.layer_Y_signs`
  - `self.ps`
  - `self.number_of_parallel_branch`
  - `self.number_of_winding_layer`
  - `self.p`
  - `self.dict_coil_connection`
  - `self.SPP`
- **Writes (State Changes)**:
  - `self.bool_DPNVorSEPA`
  - `self.bool_CustomizedCircuit`
  - `self.bool_3PhaseCurrentSource`
  - `self.layer_X_phases`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.Qs`
  - `self.CommutatingSequenceD`
  - `self.m`
  - `self.layer_X_signs`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.grouping_AC`
  - `self.layer_Y_phases`
  - `self.kw1`
  - `self.layer_Y_signs`
  - `self.ps`
  - `self.number_of_parallel_branch`
  - `self.number_of_winding_layer`
  - `self.p`
  - `self.dict_coil_connection`
  - `self.SPP`

### `Winding.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y` (Line 196)
- **Arguments**: self, layer_X_phases, coil_pitch
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Winding.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y` (Line 198)
- **Arguments**: self, layer_X_signs, coil_pitch
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Winding.get_winding_factor` (Line 202)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.p`
  - `self.m`
  - `self.ps`
  - `self.coil_pitch_y`
  - `self.Qs`
- **Writes**: None state variables detected

### `Winding.draw_winding_in_the_slot` (Line 208)
- **Arguments**: u, Qs, list_layer_phases, list_layer_signs, text
- **Reads (State Dependencies)**:
  - `self.path2SwarmData`
- **Writes**: None state variables detected

### `Winding.plot_winding_function` (Line 251)
- **Arguments**: wily
- **Reads (State Dependencies)**:
  - `self.path2SwarmData`
- **Writes**: None state variables detected

### `Winding.to_dict` (Line 263)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.number_of_parallel_branch`
  - `self.p`
  - `self.kw1`
  - `self.ps`
  - `self.derivation.__dict__`
  - `self.m`
  - `self.derivation.__dict__.items`
  - `self.coil_pitch_y`
  - `self.Qs`
  - `self.derivation`
- **Writes**: None state variables detected

### `Winding.from_dict` (Line 317)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Geometry.__init__` (Line 362)
- **Arguments**: self, name, GP, draw_function, color
- **Reads (State Dependencies)**:
  - `self.name`
  - `self.GP.items`
  - `self.GP`
  - `self.color`
  - `self.draw_function`
- **Writes (State Changes)**:
  - `self.name`
  - `self.draw_function`
  - `self.color`
  - `self.GP`

### `Geometry.update_from_GP` (Line 371)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.GP.items`
  - `self.GP`
- **Writes**: None state variables detected

### `Geometry.__repr__` (Line 382)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.name`
  - `self.GP.items`
  - `self.color`
  - `self.GP`
- **Writes**: None state variables detected

### `Geometry.print_parameters` (Line 393)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.visualization_points`
  - `self.GP.items`
  - `self.name`
  - `self.GP`
  - `self.__dict__`
  - `self.__dict__.items`
  - `self.visualization_points.items`
  - `self.color`
- **Writes**: None state variables detected

### `Geometry.to_dict` (Line 439)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.__dict__`
  - `self.__dict__.items`
- **Writes**: None state variables detected

### `Geometry.from_dict` (Line 484)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Geometry.draw` (Line 490)
- **Arguments**: self, drawer
- **Reads (State Dependencies)**:
  - `self.components_make_region`
  - `self.visualization_points`
  - `self.draw_function`
- **Writes (State Changes)**:
  - `self.components_make_region`
  - `self.visualization_points`

### `CairoDrawer.__init__` (Line 514)
- **Arguments**: self, width_in_points, height_in_points, filename, verbose_drawing, scale, bFillRegion
- **Reads (State Dependencies)**:
  - `self.ctx.set_source_rgb`
  - `self.ctx.restore`
  - `self.bFillRegion`
  - `self.ctx.scale`
  - `self.scale`
  - `self.surface`
  - `self.bMirror`
  - `self.iRotateCopy`
  - `self.ctx`
  - `self.ctx.transform`
  - `self.filename`
  - `self.verbose_drawing`
  - `self.ctx.save`
  - `self.ctx.paint`
- **Writes (State Changes)**:
  - `self.bFillRegion`
  - `self.scale`
  - `self.surface`
  - `self.bMirror`
  - `self.iRotateCopy`
  - `self.ctx`
  - `self.verbose_drawing`
  - `self.filename`

### `CairoDrawer.apply_stroke` (Line 537)
- **Arguments**: self, lw
- **Reads (State Dependencies)**:
  - `self.ctx.set_line_width`
  - `self.ctx.stroke`
  - `self.ctx.set_source_rgba`
  - `self.ctx`
  - `self.ctx.set_line_cap`
- **Writes**: None state variables detected

### `CairoDrawer.convert_to_pdf` (Line 544)
- **Arguments**: self, bool_open_pdf, filename
- **Reads (State Dependencies)**:
  - `self.sketch_color`
  - `self.surface.finish`
  - `self.surface`
- **Writes (State Changes)**:
  - `self.sketch_color`

### `CairoDrawer.hex_to_rgb` (Line 561)
- **Arguments**: self, hex_color
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `CairoDrawer.drawLine` (Line 577)
- **Arguments**: self, p1, p2
- **Reads (State Dependencies)**:
  - `self.verbose_drawing`
- **Writes**: None state variables detected

### `CairoDrawer.drawArc` (Line 582)
- **Arguments**: self, centerxy, startxy, endxy
- **Reads (State Dependencies)**:
  - `self.verbose_drawing`
- **Writes**: None state variables detected

### `CairoDrawer.getSketch` (Line 587)
- **Arguments**: self, name, color
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `CairoDrawer.prepareSection` (Line 591)
- **Arguments**: self, region_dict, color
- **Reads (State Dependencies)**:
  - `self.ctx.new_path`
  - `self.ctx.restore`
  - `self.bFillRegion`
  - `self.ctx.set_line_width`
  - `self.ctx.scale`
  - `self.ctx.save`
  - `self.ctx.fill_preserve`
  - `self.ctx.line_to`
  - `self.ctx.rotate`
  - `self.ctx.set_source_rgba`
  - `self.ctx.stroke`
  - `self.hex_to_rgb`
  - `self.ctx`
  - `self.ctx.arc_negative`
  - `self.ctx.move_to`
  - `self.ctx.arc`
- **Writes**: None state variables detected

### `CairoDrawer.finalize_part` (Line 674)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.__init__` (Line 679)
- **Arguments**: self, specs
- **Reads (State Dependencies)**:
  - `self.geometry`
  - `self.flag_do_not_evaluate_when_init_pop`
  - `self.specs`
  - `self.target`
  - `self.winding`
- **Writes (State Changes)**:
  - `self.geometry`
  - `self.flag_do_not_evaluate_when_init_pop`
  - `self.specs`
  - `self.target`
  - `self.winding`

### `Modern_Machine_Designer_Utility._get_parameter_logger` (Line 689)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.name`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.apply_parameter_dict` (Line 710)
- **Arguments**: self, prev_params, key_map
- **Reads (State Dependencies)**:
  - `self._get_parameter_logger`
  - `self.geometry`
  - `self.geometry.machineGeometry.values`
  - `self.geometry.machineGeometry`
  - `self.parameter_dict_by_name.get`
  - `self.parameter_dict_by_name`
  - `self.get_parameters_by_type`
  - `self.get_parameter_dict_by_name`
- **Writes (State Changes)**:
  - `self.parameter_dict_by_name`

### `Modern_Machine_Designer_Utility.update_geometric_parameters` (Line 780)
- **Arguments**: self, x_denorm, x_denorm_dict
- **Reads (State Dependencies)**:
  - `self.geometry`
  - `self.get_free_variables`
  - `self.geometry.machineGeometry`
  - `self.geometry.machineGeometry.items`
  - `self.parameter_dict_by_name`
  - `self.get_parameters_by_type`
  - `self.get_parameter_dict_by_name`
- **Writes (State Changes)**:
  - `self.parameter_dict_by_name`

### `Modern_Machine_Designer_Utility.get_pc_name` (Line 808)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.learn_about_the_archive` (Line 821)
- **Arguments**: self, prob, swarm_data, popsize, bool_plot_and_show, bool_more_info
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.write_swarm_survivor` (Line 900)
- **Arguments**: self, pop, counter_fitness_return
- **Reads (State Dependencies)**:
  - `self.path2SwarmData`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_bad_fintess_values` (Line 916)
- **Arguments**: self, machine_type, ref
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_rotor_volume` (Line 940)
- **Arguments**: self, stack_length
- **Reads (State Dependencies)**:
  - `self.winding.EX['mm_stack_length_specified']`
  - `self.winding.EX`
  - `self.winding`
  - `self.mm_r_ro.value`
  - `self.mm_r_ro`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_rotor_weight` (Line 948)
- **Arguments**: self, gravity, stack_length
- **Reads (State Dependencies)**:
  - `self.get_rotor_volume`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_free_variables` (Line 965)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_free_variables_as_dict` (Line 968)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_free_variables`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_free_variable_bounds_dict` (Line 979)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_free_variables`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.set_free_variables_from_dict` (Line 991)
- **Arguments**: self, free_variables_dict
- **Reads (State Dependencies)**:
  - `self.get_parameter`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.update_derived_parameters` (Line 998)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter_fields` (Line 1016)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameters_by_type` (Line 1031)
- **Arguments**: self, param_type
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter` (Line 1044)
- **Arguments**: self, name
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.set_parameter_value` (Line 1057)
- **Arguments**: self, name, value
- **Reads (State Dependencies)**:
  - `self.get_parameter`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter_dict` (Line 1074)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter_dict_by_name` (Line 1083)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.list_parameters` (Line 1092)
- **Arguments**: self, param_type
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
  - `self.get_parameters_by_type`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameters_summary` (Line 1106)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.validate_parameters` (Line 1135)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.__repr__` (Line 1167)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.target`
  - `self.target.machine_class`
  - `self.get_parameters_summary`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.to_dict` (Line 1176)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.to_json` (Line 1218)
- **Arguments**: self, indent, ensure_ascii
- **Reads (State Dependencies)**:
  - `self.to_dict`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.save_to_file` (Line 1231)
- **Arguments**: self, filepath, indent
- **Reads (State Dependencies)**:
  - `self.to_dict`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.to_dict_full` (Line 1242)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.bool_StatorSlotClosed`
  - `self.name`
  - `self.winding`
  - `self.winding.wily`
  - `self.geometry.machineGeometry`
  - `self.select_fea_config_dict`
  - `self.bool_PermanentMagnet`
  - `self.target.machine_class`
  - `self.geometry.machineGeometry.items`
  - `self.winding.wily.to_dict`
  - `self.select_FEA_tool`
  - `self.geometry`
  - `self.winding.EX.copy`
  - `self.__class__`
  - `self.get_parameter_fields`
  - `self.__class__.__name__`
  - `self.counter`
  - `self.__class__.__module__`
  - `self.target`
  - `self.winding.EX`
  - `self.bool_RotorNotched`
  - `self.bool_jmagDeleteResultsAfterCalculation`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.save_to_file_full` (Line 1398)
- **Arguments**: self, filepath, indent
- **Reads (State Dependencies)**:
  - `self.to_dict_full`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.from_dict_full` (Line 1410)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.load_from_file_full` (Line 1551)
- **Arguments**: cls, filepath
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.from_dict` (Line 1566)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.from_json` (Line 1638)
- **Arguments**: cls, json_str
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.load_from_file` (Line 1652)
- **Arguments**: cls, filepath
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.remove_jfiles_folders` (Line 1668)
- **Arguments**: root_dir
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.myLogger` (Line 1690)
- **Arguments**: dir_log, prefix
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.draw_individual_from_swarm` (Line 1719)
- **Arguments**: self, index
- **Reads (State Dependencies)**:
  - `self.show_geometry`
  - `self.path2SwarmData`
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.decode_py_reduce_ordered_dict` (Line 1763)
- **Arguments**: x_denorm_dict_raw
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.__init__` (Line 1806)
- **Arguments**: self, fname, desired_x_denorm_dict, bool_filter_pareto_front
- **Reads (State Dependencies)**:
  - `self.swarm_data_project_names`
  - `self.swarm_data_as_dict`
  - `self.get_metric_of_the_whole_swarm`
  - `self.swarm_data_xf['0']`
  - `self.number_of_free_variables`
  - `self.swarm_data_xf`
  - `self.decode_py_reduce_ordered_dict`
  - `self.swarm_data_xf.append`
  - `self.filter_data`
  - `self.number_of_chromosome`
- **Writes (State Changes)**:
  - `self.swarm_data_project_names`
  - `self.swarm_data_as_dict`
  - `self.number_of_free_variables`
  - `self.swarm_data_xf`
  - `self.number_of_chromosome`

### `Swarm_Data_Analyzer.filter_data` (Line 1975)
- **Arguments**: self, data, param_type, filter_key, direction, filter_value
- **Reads (State Dependencies)**:
  - `self.decode_py_reduce_ordered_dict`
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.decode` (Line 1998)
- **Arguments**: d
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.get_metric_of_the_whole_swarm` (Line 2004)
- **Arguments**: self, metric
- **Reads (State Dependencies)**:
  - `self.swarm_data_as_dict`
  - `self.swarm_data_as_dict.items`
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.prepare_data_for_post_processing` (Line 2024)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.torque_average`
  - `self.l_rated_windage_loss`
  - `self.f2`
  - `self.f3`
  - `self.Cost_PM`
  - `self.FRW`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self.l_rated_total_loss`
  - `self.force_error_angle`
  - `self.Cost_Cu`
  - `self.Ea`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.Cost_Fe`
  - `self.get_metric_of_the_whole_swarm`
  - `self.f1`
  - `self.rotor_weight`
  - `self.Tripple`
  - `self.Cost`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.Em`
  - `self.normalized_force_error_magnitude`
  - `self.PowerFactor`
  - `self.l_rated_magnet_Joule_loss`
  - `self.TorqueRipple`
  - `self.l_rated_iron_loss`
  - `self.swarm_data_xf`
  - `self.ss_avg_force_magnitude`
  - `self.TRV`
  - `self.RatedStkLen`
  - `self.RatedEfficiency`
- **Writes (State Changes)**:
  - `self.torque_average`
  - `self.l_rated_windage_loss`
  - `self.f2`
  - `self.f3`
  - `self.Cost_PM`
  - `self.FRW`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self.l_rated_total_loss`
  - `self.force_error_angle`
  - `self.Cost_Cu`
  - `self.Ea`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.Cost_Fe`
  - `self.f1`
  - `self.rotor_weight`
  - `self.Tripple`
  - `self.Cost`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.Em`
  - `self.normalized_force_error_magnitude`
  - `self.PowerFactor`
  - `self.l_rated_magnet_Joule_loss`
  - `self.TorqueRipple`
  - `self.l_rated_iron_loss`
  - `self.ss_avg_force_magnitude`
  - `self.TRV`
  - `self.RatedStkLen`
  - `self.RatedEfficiency`

### `swarm_data_container.__init__` (Line 2110)
- **Arguments**: self, swarm_data_raw, fea_config_dict, swarm_data_json, swarm_data_json_file_path
- **Reads (State Dependencies)**:
  - `self._load_from_json`
  - `self.fea_config_dict`
  - `self.swarm_data_raw`
  - `self._load_from_raw`
  - `self._initialize_empty`
- **Writes (State Changes)**:
  - `self.swarm_data_raw`
  - `self.fea_config_dict`

### `swarm_data_container._initialize_empty` (Line 2140)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.deg_alpha_st`
  - `self.l_force_error_angle`
  - `self.l_rated_rotor_volume`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_efficiency`
  - `self.mm_w_st`
  - `self.l_OC`
  - `self.l_rated_windage_loss`
  - `self.l_design_parameters`
  - `self.l_rated_rotor_weight`
  - `self.FRW`
  - `self.l_OB`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self.l_rated_total_loss`
  - `self.l_power_factor`
  - `self.machine_data`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.Ea`
  - `self.Trip`
  - `self.l_original_stack_length`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.l_original_rotor_weight`
  - `self.l_normalized_torque_ripple`
  - `self.l_OA`
  - `self.mm_r_si`
  - `self.RatedVol`
  - `self.l_FRW`
  - `self.project_names`
  - `self.Em`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_rated_iron_loss`
  - `self.number_of_free_variables`
  - `self.l_TRV`
  - `self.rated_data`
  - `self.swarm_data_xf`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_efficiency`
  - `self.RatedWeight`
  - `self.RatedStkLen`
- **Writes (State Changes)**:
  - `self.deg_alpha_st`
  - `self.l_force_error_angle`
  - `self.l_rated_rotor_volume`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_efficiency`
  - `self.mm_w_st`
  - `self.l_OC`
  - `self.l_rated_windage_loss`
  - `self.l_design_parameters`
  - `self.l_rated_rotor_weight`
  - `self.FRW`
  - `self.l_OB`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self.l_rated_total_loss`
  - `self.l_power_factor`
  - `self.machine_data`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.Ea`
  - `self.Trip`
  - `self.l_original_stack_length`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.l_original_rotor_weight`
  - `self.l_normalized_torque_ripple`
  - `self.l_OA`
  - `self.mm_r_si`
  - `self.RatedVol`
  - `self.l_FRW`
  - `self.project_names`
  - `self.Em`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_rated_iron_loss`
  - `self.number_of_free_variables`
  - `self.l_TRV`
  - `self.rated_data`
  - `self.swarm_data_xf`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_efficiency`
  - `self.RatedWeight`
  - `self.RatedStkLen`

### `swarm_data_container._load_from_json` (Line 2187)
- **Arguments**: self, swarm_data_json
- **Reads (State Dependencies)**:
  - `self.Ea.append`
  - `self.Em.append`
  - `self._initialize_empty`
  - `self.Trip.append`
  - `self.FRW`
  - `self.swarm_data_xf['0']`
  - `self._extract_performance_lists`
  - `self.RatedStkLen.append`
  - `self.machine_data`
  - `self.project_names.append`
  - `self.RatedWeight.append`
  - `self.Ea`
  - `self.Trip`
  - `self.RatedVol`
  - `self.machine_data.append`
  - `self.swarm_data_xf.append`
  - `self.project_names`
  - `self.RatedVol.append`
  - `self.Em`
  - `self.rated_data`
  - `self.number_of_free_variables`
  - `self.swarm_data_xf`
  - `self.rated_data.append`
  - `self.RatedWeight`
  - `self.RatedStkLen`
  - `self.FRW.append`
- **Writes (State Changes)**:
  - `self.number_of_free_variables`

### `swarm_data_container._extract_performance_lists` (Line 2301)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.l_force_error_angle`
  - `self.l_rated_rotor_volume`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_efficiency`
  - `self.l_OC`
  - `self.l_rated_windage_loss`
  - `self.l_design_parameters`
  - `self.l_rated_rotor_weight`
  - `self.l_OB`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self.l_rated_total_loss`
  - `self.l_power_factor`
  - `self.machine_data`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_original_stack_length`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.l_original_rotor_weight`
  - `self.l_normalized_torque_ripple`
  - `self.l_OA`
  - `self.l_FRW`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_rated_iron_loss`
  - `self.rated_data`
  - `self.l_TRV`
  - `self.swarm_data_xf`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_efficiency`
- **Writes (State Changes)**:
  - `self.l_force_error_angle`
  - `self.l_rated_rotor_volume`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_efficiency`
  - `self.l_OC`
  - `self.l_rated_windage_loss`
  - `self.l_design_parameters`
  - `self.l_rated_rotor_weight`
  - `self.l_OB`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self.l_rated_total_loss`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_original_stack_length`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.l_original_rotor_weight`
  - `self.l_normalized_torque_ripple`
  - `self.l_OA`
  - `self.l_FRW`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_rated_iron_loss`
  - `self.l_TRV`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_efficiency`

### `swarm_data_container._load_from_raw` (Line 2346)
- **Arguments**: self, swarm_data_raw
- **Reads (State Dependencies)**:
  - `self.deg_alpha_st`
  - `self.Ea.append`
  - `self.l_force_error_angle`
  - `self.l_rated_rotor_volume`
  - `self.l_normalized_force_error_magnitude`
  - `self.Em.append`
  - `self.l_rated_efficiency`
  - `self.mm_w_st`
  - `self.l_OC`
  - `self.l_rated_windage_loss`
  - `self.mm_r_si.append`
  - `self._initialize_empty`
  - `self.l_design_parameters`
  - `self.Trip.append`
  - `self.l_rated_rotor_weight`
  - `self.FRW`
  - `self.swarm_data_xf['0']`
  - `self.l_OB`
  - `self.deg_alpha_st.append`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self._extract_performance_lists`
  - `self.l_rated_total_loss`
  - `self.RatedStkLen.append`
  - `self.l_power_factor`
  - `self.machine_data`
  - `self.l_rated_shaft_power`
  - `self.RatedWeight.append`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.Ea`
  - `self.Trip`
  - `self.l_original_stack_length`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.project_names.append`
  - `self.l_original_rotor_weight`
  - `self.l_normalized_torque_ripple`
  - `self.l_OA`
  - `self.mm_r_si`
  - `self.RatedVol`
  - `self.machine_data.append`
  - `self.mm_w_st.append`
  - `self.swarm_data_xf.append`
  - `self.l_FRW`
  - `self.project_names`
  - `self.RatedVol.append`
  - `self.Em`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_rated_iron_loss`
  - `self.number_of_free_variables`
  - `self.l_TRV`
  - `self.rated_data`
  - `self.swarm_data_xf`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_efficiency`
  - `self.rated_data.append`
  - `self.RatedWeight`
  - `self.RatedStkLen`
  - `self.FRW.append`
- **Writes (State Changes)**:
  - `self.deg_alpha_st`
  - `self.l_force_error_angle`
  - `self.l_rated_rotor_volume`
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_efficiency`
  - `self.mm_w_st`
  - `self.l_OC`
  - `self.l_rated_windage_loss`
  - `self.l_design_parameters`
  - `self.l_rated_rotor_weight`
  - `self.l_OB`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_rated_stack_length`
  - `self.l_rated_total_loss`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_original_stack_length`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.l_original_rotor_weight`
  - `self.l_normalized_torque_ripple`
  - `self.l_OA`
  - `self.mm_r_si`
  - `self.l_FRW`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_rated_iron_loss`
  - `self.number_of_free_variables`
  - `self.l_TRV`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_efficiency`

### `swarm_data_container.get_list_y_data` (Line 2579)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.l_force_error_angle`
  - `self.l_OB`
  - `self.l_TRV`
- **Writes**: None state variables detected

### `swarm_data_container.sensitivity_bar_charts` (Line 2592)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.fea_config_dict`
  - `self.reference_data['2']`
  - `self.l_force_error_angle`
  - `self.l_normalized_force_error_magnitude`
  - `self.reference_design`
  - `self.l_OC`
  - `self.reference_data`
  - `self.reference_data['6']`
  - `self.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
  - `self.l_OB`
  - `self.l_rated_total_loss`
  - `self.reference_design['3'].split`
  - `self.get_certain_objective_function`
  - `self.l_original_rotor_weight`
  - `self.l_normalized_torque_ripple`
  - `self.required_torque`
  - `self.reference_data['4']`
  - `self.l_OA`
  - `self.rotor_weight`
  - `self.reference_design['3']`
  - `self.l_ss_avg_force_magnitude`
  - `self.rotor_volume`
  - `self.number_of_free_variables`
  - `self.reference_data['5']`
  - `self.reference_data['3']`
- **Writes (State Changes)**:
  - `self.reference_data`

## File: `codes4\Problem_BearinglessSynchronousDesign.py`

### `Problem_BearinglessSynchronousDesign.fitness` (Line 22)
- **Arguments**: self, x
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Problem_BearinglessSynchronousDesign.get_nobj` (Line 194)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Problem_BearinglessSynchronousDesign.get_bounds` (Line 199)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Problem_BearinglessSynchronousDesign.get_name` (Line 206)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_prob` (Line 217)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\ReactDrawer.py`

### `ReactDrawer.__init__` (Line 4)
- **Arguments**: self, verbose_drawing
- **Reads (State Dependencies)**:
  - `self.verbose_drawing`
  - `self.regions`
- **Writes (State Changes)**:
  - `self.verbose_drawing`
  - `self.regions`

### `ReactDrawer.hex_to_rgb` (Line 8)
- **Arguments**: self, hex_color
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `ReactDrawer.drawLine` (Line 24)
- **Arguments**: self, p1, p2
- **Reads (State Dependencies)**:
  - `self.verbose_drawing`
- **Writes**: None state variables detected

### `ReactDrawer.drawArc` (Line 29)
- **Arguments**: self, centerxy, startxy, endxy
- **Reads (State Dependencies)**:
  - `self.verbose_drawing`
- **Writes**: None state variables detected

### `ReactDrawer.getSketch` (Line 55)
- **Arguments**: self, name, color
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `ReactDrawer.prepareSection` (Line 58)
- **Arguments**: self, region_dict, color
- **Reads (State Dependencies)**:
  - `self.regions.append`
  - `self.regions`
- **Writes**: None state variables detected

### `ReactDrawer.draw_machine` (Line 63)
- **Arguments**: self, machine
- **Reads (State Dependencies)**:
  - `machine.parts`
  - `machine.all_points`
  - `self.regions`
  - `self.prepareSection`
- **Writes (State Changes)**:
  - `self.regions`

## File: `codes4\utility.py`

### `my_execfile` (Line 12)
- **Arguments**: filename, g, l
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `json_dump_ignoring_unserializable` (Line 16)
- **Arguments**: obj
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `ExceptionReTry.__init__` (Line 43)
- **Arguments**: self, message, payload
- **Reads (State Dependencies)**:
  - `self.payload`
  - `self.message`
- **Writes (State Changes)**:
  - `self.payload`
  - `self.message`

### `ExceptionReTry.__str__` (Line 46)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.message`
- **Writes**: None state variables detected

### `ExceptionBadNumberOfParts.__init__` (Line 51)
- **Arguments**: self, message, payload
- **Reads (State Dependencies)**:
  - `self.payload`
  - `self.message`
- **Writes (State Changes)**:
  - `self.payload`
  - `self.message`

### `ExceptionBadNumberOfParts.__str__` (Line 54)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.message`
- **Writes**: None state variables detected

### `communicate_database` (Line 57)
- **Arguments**: spec
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_index_and_max` (Line 149)
- **Arguments**: the_list
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_index_and_min` (Line 153)
- **Arguments**: the_list
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `gcd` (Line 157)
- **Arguments**: a, b
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `lcm` (Line 162)
- **Arguments**: a, b
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `myLogger` (Line 167)
- **Arguments**: dir_log, prefix
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `blockPrinting` (Line 218)
- **Arguments**: func
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `blockPrint` (Line 229)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `enablePrint` (Line 232)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `to_precision` (Line 238)
- **Arguments**: x, p
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `singleSidedDFT` (Line 301)
- **Arguments**: signal, samp_freq
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `basefreqDFT` (Line 309)
- **Arguments**: signal, samp_freq, ax_time_domain, ax_freq_domain, base_freq
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Pyrhonen_design.__init__` (Line 334)
- **Arguments**: self, im, bounds
- **Reads (State Dependencies)**:
  - `self.Width_StatorTeethHeadThickness`
  - `self.SIesign_parameters_denorm`
  - `self.rotor_tooth_width_b_dr`
  - `self.b1`
  - `self.Length_HeadNeckRotorSlot`
  - `self.show_norm`
  - `self.Angle_StatorSlotOpen`
  - `self.air_gap_length_delta`
  - `self.stator_tooth_width_b_ds`
- **Writes (State Changes)**:
  - `self.Width_StatorTeethHeadThickness`
  - `self.SIesign_parameters_denorm`
  - `self.rotor_tooth_width_b_dr`
  - `self.b1`
  - `self.Length_HeadNeckRotorSlot`
  - `self.Angle_StatorSlotOpen`
  - `self.air_gap_length_delta`
  - `self.stator_tooth_width_b_ds`

### `Pyrhonen_design.show_denorm` (Line 375)
- **Arguments**: self, bounds, design_parameters_norm
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Pyrhonen_design.show_norm` (Line 383)
- **Arguments**: self, bounds, design_parameters_denorm
- **Reads (State Dependencies)**:
  - `self.SIesign_parameters_norm`
  - `self.SIesign_parameters_norm.tolist`
- **Writes (State Changes)**:
  - `self.SIesign_parameters_norm`

### `add_Pyrhonen_design_to_first_generation` (Line 402)
- **Arguments**: sw, de_config_dict, logger
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `send_notification` (Line 424)
- **Arguments**: text, subject
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_windage_loss` (Line 436)
- **Arguments**: im_variant, mm_stack_length, TEMPERATURE_OF_AIR
- **Reads (State Dependencies)**:
  - `im_variant.winding.EX['RatedSpeed']`
- **Writes**: None state variables detected

### `suspension_force_vector.__init__` (Line 515)
- **Arguments**: self, force_x, force_y, range_ss
- **Reads (State Dependencies)**:
  - `self.force_err_abs`
  - `self.force_ang`
  - `self.force_x`
  - `self.ss_avg_force_vector['1']`
  - `self.force_abs`
  - `self.ss_max_force_err_ang['1']`
  - `self.ss_max_force_err_ang`
  - `self.ss_avg_force_vector['0']`
  - `self.range_ss`
  - `self.ss_avg_force_vector`
  - `self.force_ang.append`
  - `self.force_err_ang_old_way`
  - `self.force_error_angle`
  - `self.force_y`
  - `self.ss_max_force_err_abs['0']`
  - `self.ss_max_force_err_ang['0']`
  - `self.ss_max_force_err_abs`
  - `self.ss_avg_force_angle`
  - `self.force_err_ang_new_way`
  - `self.force_err_ang`
  - `self.normalized_force_error_magnitude`
  - `self.ss_avg_force_magnitude`
  - `self.ss_max_force_err_abs['1']`
- **Writes (State Changes)**:
  - `self.range_ss`
  - `self.force_err_abs`
  - `self.force_ang`
  - `self.force_err_ang`
  - `self.force_error_angle`
  - `self.force_y`
  - `self.normalized_force_error_magnitude`
  - `self.ss_avg_force_vector`
  - `self.force_err_ang_old_way`
  - `self.force_x`
  - `self.ss_max_force_err_abs`
  - `self.force_abs`
  - `self.ss_avg_force_magnitude`
  - `self.ss_max_force_err_ang`
  - `self.ss_avg_force_angle`
  - `self.force_err_ang_new_way`

### `pyplot_clear` (Line 565)
- **Arguments**: axeses
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `read_csv_results_4_comparison__transient` (Line 580)
- **Arguments**: study_name, path_prefix
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `read_csv_results_4_comparison_eddycurrent` (Line 653)
- **Arguments**: study_name, path_prefix
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `collect_jmag_Tran2TSSProlong_results` (Line 689)
- **Arguments**: im_variant, path_prefix, fea_config_dict, axeses, femm_solver_data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `csv_row_reader` (Line 787)
- **Arguments**: handle
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `whole_row_reader` (Line 792)
- **Arguments**: reader
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `get_copper_loss_Bolognani` (Line 798)
- **Arguments**: stator_slot_area, rotor_slot_area, STATOR_SLOT_FILL_FACTOR, ROTOR_SLOT_FILL_FACTOR, TEMPERATURE_OF_COIL, copper_loss_parameters
- **Reads (State Dependencies)**:
  - `EX['DriveW_zQ']`
  - `self.im.Qr`
  - `self.im.template.SI`
  - `self.im`
  - `self.im.DriveW_poles`
  - `self.im.template.SI['GP']`
  - `self.im.stack_length`
  - `self.list_rotor_current_amp`
  - `self.im.rotor_slot_height_h_sr`
  - `self.im.design_parameters['2']`
  - `self.im.template.SI['GP']['mm_r_ro']`
  - `EX['WindingFill']`
  - `self.im.template`
  - `self.im.design_parameters`
  - `EX['wily']`
  - `self.im.template.SI['GP']['mm_r_ro'].value`
  - `EX['Js']`
  - `EX['wily'].number_parallel_branch`
- **Writes**: None state variables detected

### `check_csv_results_4_general_purpose` (Line 900)
- **Arguments**: study_name, path_prefix, returnBoolean
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Goertzel_Data_Struct.__init__` (Line 961)
- **Arguments**: self, id
- **Reads (State Dependencies)**:
  - `self.ampl`
  - `self.bool_initialized`
  - `self.cosine`
  - `self.q`
  - `self.accumSquaredData`
  - `self.imag`
  - `self.q2`
  - `self.scalingFactor`
  - `self.sine`
  - `self.id`
  - `self.coeff`
  - `self.k`
  - `self.real`
  - `self.phase`
  - `self.count`
- **Writes (State Changes)**:
  - `self.ampl`
  - `self.bool_initialized`
  - `self.cosine`
  - `self.q`
  - `self.accumSquaredData`
  - `self.imag`
  - `self.q2`
  - `self.scalingFactor`
  - `self.sine`
  - `self.id`
  - `self.coeff`
  - `self.k`
  - `self.real`
  - `self.phase`
  - `self.count`

### `Goertzel_Data_Struct.goertzel_realtime` (Line 983)
- **Arguments**: gs, targetFreq, numSamples, samplingRate, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Goertzel_Data_Struct.goertzel_offline` (Line 1024)
- **Arguments**: gs, targetFreq, samplingRate, data_list
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `compute_power_factor_from_half_period` (Line 1065)
- **Arguments**: voltage, current, mytime, targetFreq, numPeriodicalExtension
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `compute_power_factor_from_full_period` (Line 1100)
- **Arguments**: voltage, current, mytime, targetFreq, numPeriodicalExtension
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `max_indices_2` (Line 1131)
- **Arguments**: arr, k
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `min_indices` (Line 1149)
- **Arguments**: arr, k
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `max_indices` (Line 1156)
- **Arguments**: arr, k
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `autolabel` (Line 1166)
- **Arguments**: ax, rects, xpos, bias, textfont
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `efficiency_at_50kW` (Line 1183)
- **Arguments**: total_loss
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `use_weights` (Line 1186)
- **Arguments**: which
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `compute_list_cost` (Line 1199)
- **Arguments**: weights, rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `fobj_scalar` (Line 1225)
- **Arguments**: torque_average, ss_avg_force_magnitude, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle, total_loss, weights, rotor_volume, rotor_weight
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `fobj_list` (Line 1237)
- **Arguments**: l_torque_average, l_ss_avg_force_magnitude, l_normalized_torque_ripple, l_normalized_force_error_magnitude, l_force_error_angle, l_total_loss, weights, rotor_volume, rotor_weight
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.__init__` (Line 1256)
- **Arguments**: self, sw, spec, dir_run, run_integer, bool_sensitivity_analysis
- **Reads (State Dependencies)**:
  - `self.SIir_run`
  - `self.number_of_designs`
  - `self.buf`
  - `self.buf_length`
  - `self.reference_design`
  - `self.build_basic_info`
  - `self.run_integer`
  - `self.spec`
  - `self.sw`
- **Writes (State Changes)**:
  - `self.SIir_run`
  - `self.number_of_designs`
  - `self.buf_length`
  - `self.reference_design`
  - `self.buf`
  - `self.run_integer`
  - `self.spec`
  - `self.sw`

### `SwarmDataAnalyzer.design_display_generator` (Line 1290)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.design_parameters_generator` (Line 1294)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.list_generations` (Line 1298)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.list_cost_function` (Line 1309)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.find_individual` (Line 1315)
- **Arguments**: self, generation_index, individual_index
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.get_best_generation` (Line 1326)
- **Arguments**: self, popsize, generator, returnMore
- **Reads (State Dependencies)**:
  - `self.SIesign_parameters_generator`
  - `self.list_cost_function`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.get_list_objective_function` (Line 1353)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.get_certain_objective_function` (Line 1360)
- **Arguments**: self, which
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.get_windage_loss` (Line 1369)
- **Arguments**: self, which
- **Reads (State Dependencies)**:
  - `self.number_of_designs`
  - `self.buf`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.my_population_distribution_plots` (Line 1379)
- **Arguments**: self, de_config_dict
- **Reads (State Dependencies)**:
  - `self.SIesign_parameters_generator`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.my_scatter_plot` (Line 1440)
- **Arguments**: self, x, y, O, fig, ax, s, marker, index_list
- **Reads (State Dependencies)**:
  - `self.sw.im.template`
  - `self.stack_length`
  - `self.spec.Js`
  - `self.speed_rpm`
  - `self.best_design_denorm['0']`
  - `self.sw.fea_config_dict['use_weights']`
  - `self.Qs`
  - `self.sw`
  - `self.template.SI['GP']['mm_r_ro']`
  - `self.template.SI`
  - `self.best_design_display`
  - `self.run_integer`
  - `self.sw.im.template.SI['GP']['mm_r_ro']`
  - `self.weights_name`
  - `self.best_design_display.split`
  - `self.sw.im.template.SI`
  - `self.sw.im.template.SI['GP']['mm_r_ro'].value`
  - `self.required_torque`
  - `self.sw.im`
  - `self.sw.im.DriveW_poles`
  - `self.mec_power`
  - `self.template`
  - `self.weights_used`
  - `self.sw.fea_config_dict`
  - `self.spec.Steel`
  - `self.rotor_weight`
  - `self.sw.im.template.SI['GP']`
  - `self.Omega`
  - `self.spec`
  - `self.template.SI['GP']`
  - `self.spec.VoltageRating`
  - `self.sw.im.BeariW_poles`
  - `self.stack_length_max`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.SIesign_display_generator`
  - `self.sw.im.Radius_OuterStatorYoke`
  - `self.Qr`
  - `self.ExcitationFreqSimulated`
  - `self.spec.stator_phase_current_rms`
  - `self.str_best_design_details`
  - `self.best_design_denorm`
  - `self.rotor_volume`
  - `self.spec.Jr`
  - `self.SIesign_parameters_generator`
- **Writes (State Changes)**:
  - `self.str_best_design_details`
  - `self.best_design_denorm`
  - `self.best_design_display`

### `SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth` (Line 1560)
- **Arguments**: self, fig, ax, marker, bool_filtered
- **Reads (State Dependencies)**:
  - `self.stack_length_max`
  - `self.get_certain_objective_function`
  - `self.required_torque`
  - `self.stack_length`
  - `self.rotor_volume`
  - `self.mec_power`
  - `self.weights_used`
  - `self.rotor_weight`
  - `self.my_scatter_plot`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.pareto_plot_torque_force` (Line 1648)
- **Arguments**: self, fig2, axeses, marker
- **Reads (State Dependencies)**:
  - `self.get_certain_objective_function`
  - `self.rotor_volume`
  - `self.weights_used`
  - `self.rotor_weight`
  - `self.my_scatter_plot`
  - `self.weights_name`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.sensitivity_bar_charts` (Line 1765)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.reference_data['3']`
  - `self.reference_data['2']`
  - `self.reference_data['6']`
  - `self.get_certain_objective_function`
  - `self.required_torque`
  - `self.reference_design['3']`
  - `self.reference_data['4']`
  - `self.rotor_volume`
  - `self.sw.fea_config_dict`
  - `self.reference_design`
  - `self.sw.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
  - `self.reference_data['5']`
  - `self.rotor_weight`
  - `self.reference_design['3'].split`
  - `self.sw`
  - `self.reference_data`
- **Writes (State Changes)**:
  - `self.reference_data`

### `SwarmDataAnalyzer.build_basic_info` (Line 2127)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.stack_length`
  - `self.speed_rpm`
  - `self.sw`
  - `self.Qs`
  - `self.template.SI['GP']['mm_r_ro']`
  - `self.template.SI`
  - `self.weights_name`
  - `self.required_torque`
  - `self.mec_power`
  - `self.template`
  - `self.weights_used`
  - `self.rotor_weight`
  - `self.template.SI['GP']`
  - `self.Omega`
  - `self.spec`
  - `self.stack_length_max`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.Qr`
  - `self.ExcitationFreqSimulated`
  - `self.rotor_volume`
- **Writes (State Changes)**:
  - `self.stack_length_max`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.Qr`
  - `self.ExcitationFreqSimulated`
  - `self.required_torque`
  - `self.stack_length`
  - `self.speed_rpm`
  - `self.rotor_volume`
  - `self.mec_power`
  - `self.weights_used`
  - `self.rotor_weight`
  - `self.Omega`
  - `self.Qs`
  - `self.weights_name`

### `build_sensitivity_bar_charts` (Line 2159)
- **Arguments**: spec, sw
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `build_Pareto_plot` (Line 2164)
- **Arguments**: spec, sw
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\winding_layout.py`

### `infer_Y_layer_phases_from_X_layer_and_coil_pitch_y` (Line 9)
- **Arguments**: layer_X_phases, coil_pitch
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `infer_Y_layer_signs_from_X_layer_and_coil_pitch_y` (Line 11)
- **Arguments**: layer_X_signs, coil_pitch
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `infer_Y_layer_grpAC_from_X_layer_and_coil_pitch_y` (Line 14)
- **Arguments**: grouping_AC, coil_pitch
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_layout_v2.__init__` (Line 21)
- **Arguments**: self, DPNV_or_SEPA, Qs, p, ps, coil_pitch_y, pr, m, Wrap_Around
- **Reads (State Dependencies)**:
  - `self.distributed_or_concentrated`
  - `self.ox_distribution_three_phase`
  - `self.list_layer_motor_signs`
  - `self.kd1`
  - `self.number_parallel_branch`
  - `self.bool_3PhaseCurrentSource`
  - `self.layer_X_phases`
  - `self.pr`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.Qs`
  - `self.CommutatingSequenceD`
  - `self.m`
  - `self.layer_X_signs`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.grouping_AC`
  - `self.list_layer_suspension_signs`
  - `self.kp1`
  - `self.layer_Y_phases`
  - `self.ox_distribution_phase_U`
  - `self.layer_Y_signs`
  - `self.ps`
  - `self.number_winding_layer`
  - `self.p`
  - `self.dict_coil_connection`
  - `self.SIPNV_or_SEPA`
  - `self.SPP`
  - `self.list_layer_suspension_phases`
  - `self.list_layer_motor_phases`
- **Writes (State Changes)**:
  - `self.distributed_or_concentrated`
  - `self.ox_distribution_three_phase`
  - `self.list_layer_motor_signs`
  - `self.kd1`
  - `self.number_parallel_branch`
  - `self.bool_3PhaseCurrentSource`
  - `self.layer_X_phases`
  - `self.pr`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.Qs`
  - `self.CommutatingSequenceD`
  - `self.m`
  - `self.layer_X_signs`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.grouping_AC`
  - `self.list_layer_suspension_signs`
  - `self.kp1`
  - `self.layer_Y_phases`
  - `self.ox_distribution_phase_U`
  - `self.layer_Y_signs`
  - `self.ps`
  - `self.number_winding_layer`
  - `self.p`
  - `self.dict_coil_connection`
  - `self.SIPNV_or_SEPA`
  - `self.SPP`
  - `self.list_layer_suspension_phases`
  - `self.list_layer_motor_phases`

### `pole_specific_winding_with_neutral.__init__` (Line 1187)
- **Arguments**: self, Qr, p, ps, coil_pitch_y
- **Reads (State Dependencies)**:
  - `self.pairs`
- **Writes (State Changes)**:
  - `self.pairs`

### `nextpow2` (Line 1277)
- **Arguments**: L
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `periodic2pi` (Line 1283)
- **Arguments**: x
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `segmented_func` (Line 1294)
- **Arguments**: x, lst_x, lst_y
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `PhaseWinding.__init__` (Line 1319)
- **Arguments**: self, Qs, m, turns_per_slot, ox_distribution_phase_U, desc_type
- **Reads (State Dependencies)**:
  - `self.turns_per_slot`
  - `self.avg_val_of_turn_func`
  - `self.sym_turn_func`
  - `self.winding_func`
  - `self.ox_distribution_phase_U`
  - `self.turn_func`
  - `self.slot_per_phase`
  - `self.radian_between_slots`
  - `self.sym_winding_func`
  - `self.setTurnFuncObject`
  - `self.setSymPos`
  - `self.degree_between_slots`
  - `self.sym_begin_pos`
- **Writes (State Changes)**:
  - `self.turns_per_slot`
  - `self.avg_val_of_turn_func`
  - `self.sym_turn_func`
  - `self.winding_func`
  - `self.ox_distribution_phase_U`
  - `self.slot_per_phase`
  - `self.radian_between_slots`
  - `self.sym_winding_func`
  - `self.degree_between_slots`

### `PhaseWinding.setTurnFuncObject` (Line 1346)
- **Arguments**: self, ox_distribution_phase_U
- **Reads (State Dependencies)**:
  - `self.lst_x`
  - `self.turns_per_slot`
  - `self.lst_y`
  - `self.turn_func`
  - `self.radian_between_slots`
- **Writes (State Changes)**:
  - `self.lst_x`
  - `self.lst_y`
  - `self.turn_func`

### `PhaseWinding.setSymPos` (Line 1378)
- **Arguments**: self, index
- **Reads (State Dependencies)**:
  - `self.lst_x`
  - `self.lst_y.index`
  - `self.lst_y`
  - `self.sym_begin_pos`
  - `self.sym_begin_pos_2`
  - `self.sym_begin_pos_1`
- **Writes (State Changes)**:
  - `self.sym_begin_pos`
  - `self.sym_begin_pos_2`
  - `self.sym_begin_pos_1`

### `PhaseWinding.plot2piFft` (Line 1403)
- **Arguments**: self, func, Fs, L
- **Reads (State Dependencies)**:
  - `self.fig_plot2piFft`
- **Writes (State Changes)**:
  - `self.fig_plot2piFft`

### `PhaseWinding.plotFuncObj` (Line 1457)
- **Arguments**: self, func
- **Reads (State Dependencies)**:
  - `self.fig_plotFuncObj`
- **Writes (State Changes)**:
  - `self.fig_plotFuncObj`

### `winding_layout.__init__` (Line 1521)
- **Arguments**: self, DPNV_or_SEPA, Qs, p, ps
- **Reads (State Dependencies)**:
  - `self.distributed_or_concentrated`
  - `self.number_parallel_branch`
  - `self.bool_3PhaseCurrentSource`
  - `self.l42`
  - `self.Qs`
  - `self.CommutatingSequenceD`
  - `self.l21`
  - `self.l41`
  - `self.l_leftlayer1`
  - `self.coil_pitch`
  - `self.CommutatingSequenceB`
  - `self.grouping_AC`
  - `self.layer_A2`
  - `self.l22`
  - `self.l_rightlayer2`
  - `self.l_leftlayer2`
  - `self.layer_B1`
  - `self.no_winding_layer`
  - `self.initial_excitation_bias_compensation_deg`
  - `self.layer_B2`
  - `self.p`
  - `self.l_rightlayer1`
  - `self.layer_A1`
- **Writes (State Changes)**:
  - `self.distributed_or_concentrated`
  - `self.number_parallel_branch`
  - `self.bool_3PhaseCurrentSource`
  - `self.l42`
  - `self.Qs`
  - `self.CommutatingSequenceD`
  - `self.l21`
  - `self.l41`
  - `self.l_leftlayer1`
  - `self.coil_pitch`
  - `self.CommutatingSequenceB`
  - `self.grouping_AC`
  - `self.layer_A2`
  - `self.l22`
  - `self.l_rightlayer2`
  - `self.l_leftlayer2`
  - `self.layer_B1`
  - `self.no_winding_layer`
  - `self.initial_excitation_bias_compensation_deg`
  - `self.layer_B2`
  - `self.p`
  - `self.l_rightlayer1`
  - `self.layer_A1`

## File: `codes4\winding_layout_derivation_ismb2021_asymetry_no_drawing.py`

### `_print` (Line 18)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `set_verbose` (Line 24)
- **Arguments**: verbose
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `limit_to_360_deg` (Line 30)
- **Arguments**: PHI
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `belong_to_which_phase_belt` (Line 38)
- **Arguments**: PHI, phase_belt
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `belong_to_band` (Line 85)
- **Arguments**: LB, UB, PHI
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `phase_angle_of_slot_i_at_frequency_h` (Line 102)
- **Arguments**: slot_number, h, Q
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `compute_star_of_slots` (Line 108)
- **Arguments**: Q, p, m, verbose
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `compute_connection_star_at_another_frequency` (Line 146)
- **Arguments**: connection_star_raw_dict, frequency_ratio, which_phase, verbose
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_distribution_factor` (Line 206)
- **Arguments**: Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict, Aa
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_short_pitch_factor_v2` (Line 243)
- **Arguments**: h, coil_pitch_y, Q
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Winding_Derivation.__init__` (Line 253)
- **Arguments**: self, slot_pole_comb, bool_double_layer_winding, verbose
- **Reads (State Dependencies)**:
  - `self.turn_func_bias`
  - `self.list_slot_number_of_phase`
  - `self.qs`
  - `self.suspen_kd_at_h`
  - `self.dpnv_grouping_dict_b`
  - `self.Q`
  - `self.bool_double_layer_winding`
  - `self.torque_kd_at_h`
  - `self.q`
  - `self.ts`
  - `self.list_phase_u_slot_number`
  - `self.m`
  - `self.coil_pitch_y`
  - `self.t`
  - `self.verbose`
  - `self.dpnv_grouping_dict_c`
  - `self.list_phase_w_slot_number`
  - `self.ps`
  - `self.list_phase_v_slot_number`
  - `self.torque_kp_at_h`
  - `self.p`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.torque_kw_at_h`
  - `self.suspen_kp_at_h`
  - `self.suspen_kw_at_h`
- **Writes (State Changes)**:
  - `self.turn_func_bias`
  - `self.list_slot_number_of_phase`
  - `self.qs`
  - `self.suspen_kd_at_h`
  - `self.dpnv_grouping_dict_b`
  - `self.Q`
  - `self.bool_double_layer_winding`
  - `self.torque_kd_at_h`
  - `self.q`
  - `self.ts`
  - `self.list_phase_u_slot_number`
  - `self.m`
  - `self.coil_pitch_y`
  - `self.t`
  - `self.verbose`
  - `self.dpnv_grouping_dict_c`
  - `self.list_phase_w_slot_number`
  - `self.ps`
  - `self.list_phase_v_slot_number`
  - `self.torque_kp_at_h`
  - `self.p`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.torque_kw_at_h`
  - `self.suspen_kp_at_h`
  - `self.suspen_kw_at_h`

### `Winding_Derivation.get_complex_number_winding_factor_of_coil_i` (Line 444)
- **Arguments**: self, i, coil_pitch_y, Q, v, p
- **Reads (State Dependencies)**:
  - `self.verbose`
- **Writes**: None state variables detected

### `Winding_Derivation.get_complex_number_kw_per_phase` (Line 458)
- **Arguments**: self, v, p, positive_connected_coils, negative_connected_coils
- **Reads (State Dependencies)**:
  - `self.Q`
  - `self.get_complex_number_winding_factor_of_coil_i`
  - `self.verbose`
  - `self.coil_pitch_y`
- **Writes**: None state variables detected

### `Winding_Derivation.get_complex_number_kw` (Line 486)
- **Arguments**: self, p_or_ps, v, bool_study_suspension_subharmonics
- **Reads (State Dependencies)**:
  - `self.dpnv_grouping_dict_a['GAC']`
  - `self.dpnv_grouping_dict_b['GBD']`
  - `self.get_complex_number_kw_per_phase`
  - `self.verbose`
  - `self.dpnv_grouping_dict_c`
  - `self.dpnv_grouping_dict_b`
  - `self.dpnv_grouping_dict_a['GBD']`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.dpnv_grouping_dict_c['GBD']`
  - `self.ps`
  - `self.dpnv_grouping_dict_c['GAC']`
  - `self.dpnv_grouping_dict_b['GAC']`
- **Writes**: None state variables detected

### `Winding_Derivation.format_print_out_string` (Line 531)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.dpnv_grouping_dict_b`
  - `self.Q`
  - `self.layer_X_phases`
  - `self.bool_double_layer_winding`
  - `self.dict_suspension_kw_els['C_angle']`
  - `self.dict_torque_kw_cjh`
  - `self.dict_suspension_kw_els['A_angle']`
  - `self.dict_torque_kw_els`
  - `self.dict_suspension_kw_cjh`
  - `self.m`
  - `self.layer_X_signs`
  - `self.coil_pitch_y`
  - `self.verbose`
  - `self.grouping_AC`
  - `self.dpnv_grouping_dict_c`
  - `self.ps`
  - `self.get_complex_number_kw`
  - `self.print_out_string`
  - `self.dict_suspension_kw_els`
  - `self.p`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.dict_suspension_kw_els['B_angle']`
- **Writes (State Changes)**:
  - `self.grouping_AC`
  - `self.layer_X_phases`
  - `self.layer_X_signs`
  - `self.coil_pitch_y`
  - `self.print_out_string`

### `main_derivation` (Line 655)
- **Arguments**: m, Qs, p, ps, coil_pitch_y, verbose
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `codes4\winding_layout_derivation_ismb2021_asymetry_ori.py`

### `limit_to_360_deg` (Line 21)
- **Arguments**: PHI
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `belong_to_which_phase_belt` (Line 27)
- **Arguments**: PHI, phase_belt
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `belong_to_band` (Line 79)
- **Arguments**: LB, UB, PHI
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `angular_location` (Line 95)
- **Arguments**: PHI, radius_bias
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `phase_angle_of_slot_i_at_frequency_h` (Line 99)
- **Arguments**: slot_number, h, Q
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `draw_star_of_slots` (Line 107)
- **Arguments**: Q, p, m
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `draw_connection_star` (Line 172)
- **Arguments**: m, phase_belt, connection_star_raw_dict
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `draw_connection_star_at_another_frequency` (Line 239)
- **Arguments**: connection_star_raw_dict, frequency_ratio, which_phase
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_distribution_factor_verPyrhonen` (Line 323)
- **Arguments**: h, Q, p, m
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_distribution_factor` (Line 338)
- **Arguments**: Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict, Aa
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `draw_turn_function` (Line 389)
- **Arguments**: connection_star_raw_dict, coil_pitch_y, Q, phase_Aa_dpnv_grouping_list, turn_func_bias, Aa, bool_double_layer_winding
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_short_pitch_factor` (Line 537)
- **Arguments**: h, coil_pitch_y, Q, npp
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_short_pitch_factor_v2` (Line 581)
- **Arguments**: h, coil_pitch_y, Q
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Winding_Derivation.__init__` (Line 591)
- **Arguments**: self, slot_pole_comb, bool_double_layer_winding
- **Reads (State Dependencies)**:
  - `self.turn_func_bias`
  - `self.list_slot_number_of_phase`
  - `self.qs`
  - `self.drawer_T3b`
  - `self.drawer_T4c`
  - `self.suspen_kd_at_h`
  - `self.drawer_T3c`
  - `self.dpnv_grouping_dict_b`
  - `self.drawer_T2`
  - `self.Q`
  - `self.drawer_T4b`
  - `self.torque_kd_at_h`
  - `self.q`
  - `self.ts`
  - `self.list_phase_u_slot_number`
  - `self.drawer_T3a`
  - `self.m`
  - `self.coil_pitch_y`
  - `self.t`
  - `self.dpnv_grouping_dict_c`
  - `self.list_phase_w_slot_number`
  - `self.drawer_Text`
  - `self.Q_prime`
  - `self.ps`
  - `self.list_phase_v_slot_number`
  - `self.drawer_T4`
  - `self.drawer_T4a`
  - `self.torque_kp_at_h`
  - `self.p`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.torque_kw_at_h`
  - `self.suspen_kp_at_h`
  - `self.drawer_T1`
  - `self.suspen_kw_at_h`
- **Writes (State Changes)**:
  - `self.turn_func_bias`
  - `self.list_slot_number_of_phase`
  - `self.qs`
  - `self.drawer_T3b`
  - `self.drawer_T4c`
  - `self.suspen_kd_at_h`
  - `self.drawer_T3c`
  - `self.dpnv_grouping_dict_b`
  - `self.drawer_T2`
  - `self.Q`
  - `self.drawer_T4b`
  - `self.torque_kd_at_h`
  - `self.q`
  - `self.ts`
  - `self.list_phase_u_slot_number`
  - `self.drawer_T3a`
  - `self.m`
  - `self.coil_pitch_y`
  - `self.t`
  - `self.dpnv_grouping_dict_c`
  - `self.list_phase_w_slot_number`
  - `self.drawer_Text`
  - `self.Q_prime`
  - `self.ps`
  - `self.list_phase_v_slot_number`
  - `self.drawer_T4`
  - `self.drawer_T4a`
  - `self.torque_kp_at_h`
  - `self.p`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.torque_kw_at_h`
  - `self.suspen_kp_at_h`
  - `self.drawer_T1`
  - `self.suspen_kw_at_h`

### `Winding_Derivation.get_complex_number_winding_factor_of_coil_i` (Line 828)
- **Arguments**: self, i, coil_pitch_y, Q, v, p
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Winding_Derivation.get_complex_number_kw_per_phase` (Line 846)
- **Arguments**: self, v, p, positive_connected_coils, negative_connected_coils
- **Reads (State Dependencies)**:
  - `self.Q`
  - `self.get_complex_number_winding_factor_of_coil_i`
  - `self.coil_pitch_y`
- **Writes**: None state variables detected

### `Winding_Derivation.get_complex_number_kw` (Line 882)
- **Arguments**: self, p_or_ps, v, bool_study_suspension_subharmonics
- **Reads (State Dependencies)**:
  - `self.dict_kw_cjh`
  - `self.dpnv_grouping_dict_a['GAC']`
  - `self.dpnv_grouping_dict_b['GBD']`
  - `self.get_complex_number_kw_per_phase`
  - `self.dpnv_grouping_dict_c`
  - `self.dpnv_grouping_dict_b`
  - `self.dpnv_grouping_dict_a['GBD']`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_a`
  - `self.dpnv_grouping_dict_c['GBD']`
  - `self.ps`
  - `self.dpnv_grouping_dict_c['GAC']`
  - `self.dict_kw_els`
  - `self.dpnv_grouping_dict_b['GAC']`
- **Writes (State Changes)**:
  - `self.dict_kw_cjh`
  - `self.dict_kw_els`

### `main_derivation` (Line 927)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_diagram.__init__` (Line 1181)
- **Arguments**: self, layer_X_phases, layer_X_signs, coil_pitch_y, grouping_AC
- **Reads (State Dependencies)**:
  - `self.dl_grouping_AC`
  - `self.grouping_AC`
  - `self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y`
  - `self.dl_grouping_BD`
  - `self.dl_rightlayer`
  - `self.layer_Y_phases`
  - `self.layer_X_phases`
  - `self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y`
  - `self.dl_leftlayer`
  - `self.layer_Y_signs`
  - `self.layer_X_signs`
  - `self.coil_pitch_y`
- **Writes (State Changes)**:
  - `self.dl_grouping_AC`
  - `self.grouping_AC`
  - `self.dl_grouping_BD`
  - `self.dl_rightlayer`
  - `self.layer_Y_phases`
  - `self.layer_X_phases`
  - `self.dl_leftlayer`
  - `self.layer_Y_signs`
  - `self.layer_X_signs`
  - `self.coil_pitch_y`

### `winding_diagram.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y` (Line 1245)
- **Arguments**: layer_X_phases, coil_pitch
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_diagram.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y` (Line 1249)
- **Arguments**: layer_X_signs, coil_pitch
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `winding_diagram.draw` (Line 1253)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.dl_grouping_AC`
  - `self.dl_rightlayer['V']`
  - `self.dl_rightlayer['W']`
  - `self.grouping_AC`
  - `self.dl_grouping_BD`
  - `self.dl_rightlayer`
  - `self.dl_leftlayer['U']`
  - `self.dl_leftlayer['W']`
  - `self.dl_leftlayer`
  - `self.dl_rightlayer['U']`
  - `self.coil_pitch_y`
  - `self.dl_leftlayer['V']`
- **Writes**: None state variables detected

