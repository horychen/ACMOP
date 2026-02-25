# Variables & State Dependencies Report

This report extracts parameters read and modified by functions/methods across the codebase, focusing on state variables like `self.*`, `machine.*`, `EX[...]`, dictionary `.get()`, etc.

## Categorized Variables summary

### Variables Read (Inputs)
- `EX['DriveW_zQ']`
- `EX['Js']`
- `EX['WindingFill']`
- `EX['wily']`
- `EX['wily'].number_parallel_branch`
- `acm_variant.rotorMagnet.notched_rotor.p`
- `im_variant.winding.EX['RatedSpeed']`
- `motor.rotor.OD`
- `motor.rotor.airGap`
- `motor.stator.ID`
- `motor.stator.OD`
- `motor.stator.liner`
- `motor.stator.toothDepth`
- `motor.stator.toothWidth`
- `motor.stator.yoke`
- `rotor.ID`
- `rotor.OD`
- `rotor.magnetDepth`
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
- `self.La`
- `self.Length_HeadNeckRotorSlot`
- `self.Omega`
- `self.PowerFactor`
- `self.Q`
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
- `self.air_gap`
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
- `self.awg`
- `self.b1`
- `self.bFillRegion`
- `self.bMirror`
- `self.basic_info`
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
- `self.calc_bemf_constants`
- `self.calc_bounds`
- `self.calc_motor_losses`
- `self.calculate_excitation_current`
- `self.calculate_slot_area`
- `self.checkGeomApp`
- `self.circuit_current`
- `self.coeff`
- `self.coil_fluxLinkage`
- `self.coil_pitch`
- `self.coil_pitch_y`
- `self.color`
- `self.components_make_region`
- `self.connection_star_raw_dict`
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
- `self.ctx.path_extents`
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
- `self.d_tooth`
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
- `self.dict_suspension_kw_cjh`
- `self.dict_suspension_kw_els`
- `self.dict_suspension_kw_els['A_angle']`
- `self.dict_suspension_kw_els['B_angle']`
- `self.dict_suspension_kw_els['C_angle']`
- `self.dict_torque_kw_cjh`
- `self.dict_torque_kw_els`
- `self.distributed_or_concentrated`
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
- `self.draw_machine_using_CairoDrawer`
- `self.edge4Ref`
- `self.edge4ref`
- `self.estimate_max_wires_in_slot`
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
- `self.geometry.add_part`
- `self.geometry.all_points`
- `self.geometry.machineGeometry`
- `self.geometry.machineGeometry.items`
- `self.geometry.machineGeometry.values`
- `self.geometry.parts`
- `self.geometry.show_geometry_svg`
- `self.geometry.sync`
- `self.getSketch`
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
- `self.get_wire_properties`
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
- `self.l_stack`
- `self.l_stator_copper_loss_in_end_turn`
- `self.l_torque_average`
- `self.layer_A1`
- `self.layer_A2`
- `self.layer_B1`
- `self.layer_B2`
- `self.layer_X_phases`
- `self.layer_X_signs`
- `self.layer_Y_phases`
- `self.layer_Y_signs`
- `self.liner`
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
- `self.lst_x`
- `self.lst_y`
- `self.lst_y.index`
- `self.m`
- `self.machine_data`
- `self.machine_data.append`
- `self.mec_power`
- `self.message`
- `self.mm_r_ro`
- `self.mm_r_ro.value`
- `self.mm_r_si`
- `self.mm_r_si.append`
- `self.mm_w_st`
- `self.mm_w_st.append`
- `self.model`
- `self.motor`
- `self.motor.stator`
- `self.motor.stator.tooth_shape`
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
- `self.options`
- `self.overwritten`
- `self.ox_distribution_phase_U`
- `self.ox_distribution_three_phase`
- `self.p`
- `self.pairs`
- `self.parameter_dict`
- `self.parameter_dict_by_name`
- `self.parameter_dict_by_name.get`
- `self.parts`
- `self.parts.append`
- `self.path2SwarmData`
- `self.payload`
- `self.phase`
- `self.pole_count`
- `self.poles`
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
- `self.required_torque`
- `self.rotated_position`
- `self.rotor_od`
- `self.rotor_tooth_width_b_dr`
- `self.rotor_volume`
- `self.rotor_weight`
- `self.run_integer`
- `self.save`
- `self.scale`
- `self.scalingFactor`
- `self.select_FEA_tool`
- `self.select_fea_config_dict`
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
- `self.slot_count`
- `self.slot_per_phase`
- `self.slots`
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
- `self.stator`
- `self.stator.ID`
- `self.stator.OD`
- `self.stator.yoke`
- `self.stator_id`
- `self.stator_od`
- `self.stator_tooth_width_b_ds`
- `self.str_best_design_details`
- `self.study`
- `self.study_name`
- `self.surface`
- `self.surface.finish`
- `self.suspen_kd_at_h`
- `self.suspen_kp_at_h`
- `self.suspen_kw_at_h`
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
- `self.t`
- `self.target`
- `self.target.machine_class`
- `self.target_j`
- `self.template`
- `self.template.SI`
- `self.template.SI['GP']`
- `self.template.SI['GP']['mm_r_ro']`
- `self.template.SI['GP']['mm_r_ro'].value`
- `self.terminal_voltage`
- `self.time_list`
- `self.to_dict`
- `self.to_dict_full`
- `self.tooth_depth`
- `self.tooth_shape`
- `self.tooth_width`
- `self.torque_average`
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
- `self.value`
- `self.verbose`
- `self.verbose_drawing`
- `self.view`
- `self.visualization_points`
- `self.visualization_points.items`
- `self.weights_name`
- `self.weights_used`
- `self.winding`
- `self.winding.EX`
- `self.winding.EX.copy`
- `self.winding.EX['mm_stack_length_specified']`
- `self.winding.wily`
- `self.winding.wily.to_dict`
- `self.winding['l_stack']`
- `self.winding['pole_count']`
- `self.winding['slot_count']`
- `self.winding_func`
- `self.workDir`
- `self.yoke_thickness`
- `stator.ID`
- `stator.OD`
- `stator.toothDepth`
- `stator.toothWidth`
- `stator.yoke`

### Variables Written (Outputs)
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
- `self.La`
- `self.Length_HeadNeckRotorSlot`
- `self.Omega`
- `self.PowerFactor`
- `self.Q`
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
- `self.air_gap`
- `self.air_gap_length_delta`
- `self.all_points`
- `self.ampl`
- `self.app`
- `self.ass`
- `self.avg_val_of_turn_func`
- `self.awg`
- `self.b1`
- `self.bFillRegion`
- `self.bMirror`
- `self.basic_info`
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
- `self.distributed_or_concentrated`
- `self.dm`
- `self.doc`
- `self.dpnv_grouping_dict_a`
- `self.dpnv_grouping_dict_b`
- `self.dpnv_grouping_dict_c`
- `self.draw_function`
- `self.edge4Ref`
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
- `self.gp`
- `self.grouping_AC`
- `self.horizontal_position`
- `self.iRotateCopy`
- `self.id`
- `self.id_rotorCore`
- `self.id_statorCore`
- `self.imag`
- `self.initial_excitation_bias_compensation_deg`
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
- `self.layer_A1`
- `self.layer_A2`
- `self.layer_B1`
- `self.layer_B2`
- `self.layer_X_phases`
- `self.layer_X_signs`
- `self.layer_Y_phases`
- `self.layer_Y_signs`
- `self.liner`
- `self.list_layer_motor_phases`
- `self.list_layer_motor_signs`
- `self.list_layer_suspension_phases`
- `self.list_layer_suspension_signs`
- `self.list_phase_u_slot_number`
- `self.list_phase_v_slot_number`
- `self.list_phase_w_slot_number`
- `self.list_slot_number_of_phase`
- `self.lst_x`
- `self.lst_y`
- `self.m`
- `self.machine_data`
- `self.mec_power`
- `self.message`
- `self.mm_r_si`
- `self.mm_w_st`
- `self.model`
- `self.motor`
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
- `self.pole_count`
- `self.poles`
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
- `self.rated_speed`
- `self.real`
- `self.reference_data`
- `self.reference_design`
- `self.required_torque`
- `self.rotated_position`
- `self.rotor_od`
- `self.rotor_tooth_width_b_dr`
- `self.rotor_volume`
- `self.rotor_weight`
- `self.run_integer`
- `self.scale`
- `self.scalingFactor`
- `self.sine`
- `self.sketch`
- `self.sketchNameList`
- `self.sketch_color`
- `self.slot_count`
- `self.slot_per_phase`
- `self.slots`
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
- `self.stator_id`
- `self.stator_od`
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
- `self.target_j`
- `self.template.SI['GP']['mm_r_ro'].value`
- `self.time_list`
- `self.tooth_depth`
- `self.tooth_width`
- `self.torque_average`
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
- `self.value`
- `self.verbose`
- `self.verbose_drawing`
- `self.view`
- `self.visualization_points`
- `self.weights_name`
- `self.weights_used`
- `self.winding`
- `self.winding_func`
- `self.workDir`
- `self.yoke_thickness`

## File: `angle_error_nick.py`

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

## File: `JMAG.py`

### `JMAG.__init__` (Line 31)
- **Arguments**: self, fea_config_dict
- **Reads (State Dependencies)**:
  - `self.flag_material_already_loaded`
  - `self.app`
  - `self.consts`
  - `self.geomApp`
  - `self.bMirror`
  - `self.ass`
  - `self.iRotateCopy`
  - `self.model`
  - `self.doc`
  - `self.jd`
  - `self.edge4Ref`
  - `self.fea_config_dict`
  - `self.study`
  - `self.sketchNameList`
  - `self.bool_suppressShaft`
  - `self.verbose_drawing`
  - `self.workDir`
  - `self.sketch`
  - `self.view`
  - `self.projName`
  - `self.defaultUnit`
  - `self.JMAG_version_number`
- **Writes (State Changes)**:
  - `self.flag_material_already_loaded`
  - `self.app`
  - `self.consts`
  - `self.geomApp`
  - `self.bMirror`
  - `self.ass`
  - `self.iRotateCopy`
  - `self.model`
  - `self.doc`
  - `self.jd`
  - `self.edge4Ref`
  - `self.fea_config_dict`
  - `self.study`
  - `self.sketchNameList`
  - `self.bool_suppressShaft`
  - `self.verbose_drawing`
  - `self.workDir`
  - `self.sketch`
  - `self.view`
  - `self.projName`
  - `self.defaultUnit`
  - `self.JMAG_version_number`

### `JMAG.open` (Line 59)
- **Arguments**: self, Steel_name, expected_project_file_path, pc_name, dir_parent, bool_jmagDesignerShow
- **Reads (State Dependencies)**:
  - `self.flag_material_already_loaded`
  - `self.app`
  - `self.fea_config_dict['pc_name']`
  - `self.JMAG_version_number`
  - `self.fea_config_dict`
  - `self.JMAG_version_string`
- **Writes (State Changes)**:
  - `self.flag_material_already_loaded`
  - `self.app`
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
  - `self.doc.SaveModel`
  - `self.app`
  - `self.doc`
  - `self.app.GetCurrentModel`
- **Writes**: None state variables detected

### `JMAG.pre_process_PMSM` (Line 280)
- **Arguments**: self, app, model, acm_variant
- **Reads (State Dependencies)**:
  - `self.fea_config_dict['designer.show']`
  - `self.doc`
  - `self.id_statorCore`
  - `self.id_rotorCore`
  - `self.fea_config_dict`
  - `self.doc.GetSelection`
- **Writes (State Changes)**:
  - `self.id_statorCore`
  - `self.id_rotorCore`

### `JMAG.add_magnetic_transient_study` (Line 456)
- **Arguments**: self, app, model, path2FEACsv, study_name, acm_variant
- **Reads (State Dependencies)**:
  - `self.id_statorCore`
  - `self.JMAG_version_number`
  - `self.fea_config_dict['designer.max_nonlinear_iteration']`
  - `self.id_rotorCore`
  - `self.study_name`
  - `self.add_material`
  - `self.fea_config_dict`
  - `self.add_circuit`
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
  - `self.doc`
  - `self.doc.CreateReferenceFromItem`
  - `self.sketch.CreateBiConstraint`
  - `self.sketch`
  - `self.sketch.GetItem`
- **Writes**: None state variables detected

### `JMAG.drawLine` (Line 1060)
- **Arguments**: self, startxy, endxy, returnVertexName
- **Reads (State Dependencies)**:
  - `self.getSketch`
  - `self.sketch.OpenSketch`
  - `self.sketch`
  - `self.sketch.CreateLine`
- **Writes (State Changes)**:
  - `self.sketch`

### `JMAG.drawArc` (Line 1078)
- **Arguments**: self, centerxy, startxy, endxy, returnVertexName
- **Reads (State Dependencies)**:
  - `self.getSketch`
  - `self.sketch.OpenSketch`
  - `self.sketch`
  - `self.sketch.CreateArc`
- **Writes (State Changes)**:
  - `self.sketch`

### `JMAG.drawCircle` (Line 1095)
- **Arguments**: self, centerxy, radius, returnVertexName
- **Reads (State Dependencies)**:
  - `self.getSketch`
  - `self.sketch.OpenSketch`
  - `self.sketch`
  - `self.sketch.CreateCircle`
- **Writes (State Changes)**:
  - `self.sketch`

### `JMAG.checkGeomApp` (Line 1109)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.app`
  - `self.geomApp`
  - `self.doc`
  - `self.app.LaunchGeometryEditor`
  - `self.geomApp.NewDocument`
  - `self.app.CreateGeometryEditor`
- **Writes (State Changes)**:
  - `self.geomApp`
  - `self.doc`

### `JMAG.getSketch` (Line 1117)
- **Arguments**: self, sketchName, color
- **Reads (State Dependencies)**:
  - `self.doc.GetAssembly`
  - `self.geomApp`
  - `self.geomApp.GetDocument`
  - `self.ass`
  - `self.checkGeomApp`
  - `self.doc`
  - `self.ass.GetItem`
  - `self.doc.CreateReferenceFromItem`
  - `self.sketchNameList`
  - `self.sketchNameList.append`
  - `self.sketch.SetProperty`
  - `self.ass.CreateSketch`
  - `self.sketch.OpenSketch`
  - `self.sketch`
- **Writes (State Changes)**:
  - `self.geomApp`
  - `self.sketch`
  - `self.ass`
  - `self.doc`

### `JMAG.prepareSection` (Line 1142)
- **Arguments**: self, token, bMirrorMerge, bRotateMerge
- **Reads (State Dependencies)**:
  - `self.sketch.GetItem`
  - `self.bMirror`
  - `self.doc`
  - `self.sketch.CloseSketch`
  - `self.edge4Ref`
  - `self.regionMirrorCopy`
  - `self.sketch.CreateRegions`
  - `self.iRotateCopy`
  - `self.sketch`
  - `self.edge4ref`
  - `self.doc.GetSelection`
  - `self.regionCircularPattern360Origin`
- **Writes**: None state variables detected

### `JMAG.regionMirrorCopy` (Line 1197)
- **Arguments**: self, region, edge4Ref, symmetryType, bMerge
- **Reads (State Dependencies)**:
  - `self.ass`
  - `self.doc`
  - `self.doc.CreateReferenceFromItem`
  - `self.ass.GetItem`
  - `self.sketch.CreateRegionMirrorCopy`
  - `self.sketch`
  - `self.sketch.GetItem`
- **Writes**: None state variables detected

### `JMAG.regionCircularPattern360Origin` (Line 1219)
- **Arguments**: self, idx, region, Q_float, bMerge
- **Reads (State Dependencies)**:
  - `self.sketch.CreateRegionCircularPattern`
  - `self.sketch`
  - `self.doc`
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
  - `self.save`
  - `self.show`
  - `self.bMirror`
  - `acm_variant.rotorMagnet.notched_rotor.p`
  - `self.calculate_excitation_current`
  - `self.iRotateCopy`
  - `self.prepareSection`
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
  - `self.Current_dict['Time(s)']`
  - `self.FluxLinkage_dict`
  - `self.ForConY_list`
  - `self.basic_info`
  - `self.DisplacementAngle_list`
  - `self.myvoltage`
  - `self.ForConAbs_list`
  - `self.terminal_voltage`
  - `self.mytime`
  - `self.ForConX_list`
  - `self.get_voltage_and_current`
  - `self.Current_dict`
  - `self.coil_fluxLinkage`
  - `self.TorCon_list`
  - `self.mycurrent`
  - `self.ui_info`
  - `self.femm_loss_list`
  - `self.time_list`
  - `self.jmag_loss_list`
  - `self.circuit_current`
- **Writes (State Changes)**:
  - `self.TorCon_list`
  - `self.basic_info`
  - `self.ForConX_list`
  - `self.mycurrent`
  - `self.myvoltage`
  - `self.ui_info`
  - `self.ForConAbs_list`
  - `self.femm_loss_list`
  - `self.time_list`
  - `self.jmag_loss_list`
  - `self.mytime`
  - `self.ForConY_list`

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

## File: `machine.py`

### `Machine.__init__` (Line 6)
- **Arguments**: self, motor
- **Reads (State Dependencies)**:
  - `self.motor`
  - `self.geometry`
- **Writes (State Changes)**:
  - `self.motor`
  - `self.geometry`

### `Machine.sync` (Line 10)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.motor.stator.tooth_shape`
  - `self.motor.stator`
  - `self.geometry.sync`
  - `self.motor`
  - `self.geometry.add_part`
  - `self.geometry`
- **Writes**: None state variables detected

### `Machine.draw_svg` (Line 22)
- **Arguments**: self, filename
- **Reads (State Dependencies)**:
  - `self.geometry`
  - `self.geometry.show_geometry_svg`
- **Writes**: None state variables detected

### `Machine.draw_jmag` (Line 25)
- **Arguments**: self, filename
- **Reads (State Dependencies)**:
  - `self.geometry.all_points`
  - `self.geometry.parts`
  - `self.geometry`
- **Writes**: None state variables detected

## File: `machine_analyzer.py`

### `MotorPerformanceAnalyzer.__init__` (Line 16)
- **Arguments**: self, motor
- **Reads (State Dependencies)**:
  - `self.tooth_width`
  - `self.motor`
  - `self.target_j`
  - `motor.stator.liner`
  - `self.stator_id`
  - `self.stator_od`
  - `motor.stator.yoke`
  - `self.liner`
  - `self.yoke_thickness`
  - `motor.stator.OD`
  - `motor.rotor.OD`
  - `self.rotor_od`
  - `self.slots`
  - `motor.stator.ID`
  - `motor.stator.toothWidth`
  - `self.air_gap`
  - `motor.rotor.airGap`
  - `motor.stator.toothDepth`
  - `self.rated_speed`
  - `self.La`
  - `self.poles`
  - `self.awg`
  - `self.tooth_depth`
- **Writes (State Changes)**:
  - `self.rotor_od`
  - `self.slots`
  - `self.stator_od`
  - `self.tooth_width`
  - `self.rated_speed`
  - `self.motor`
  - `self.target_j`
  - `self.liner`
  - `self.air_gap`
  - `self.stator_id`
  - `self.La`
  - `self.poles`
  - `self.yoke_thickness`
  - `self.awg`
  - `self.tooth_depth`

### `MotorPerformanceAnalyzer.get_wire_properties` (Line 35)
- **Arguments**: self, awg_size
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `MotorPerformanceAnalyzer.calculate_slot_area` (Line 42)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.stator_id`
  - `self.tooth_depth`
  - `self.slots`
  - `self.tooth_width`
- **Writes**: None state variables detected

### `MotorPerformanceAnalyzer.estimate_max_wires_in_slot` (Line 52)
- **Arguments**: self, d_od
- **Reads (State Dependencies)**:
  - `self.slots`
  - `self.tooth_width`
  - `self.liner`
  - `self.stator_id`
  - `self.tooth_depth`
- **Writes**: None state variables detected

### `MotorPerformanceAnalyzer.calc_bemf_constants` (Line 89)
- **Arguments**: self, z_slot, B_sat
- **Reads (State Dependencies)**:
  - `self.La`
  - `self.poles`
  - `self.slots`
  - `self.yoke_thickness`
- **Writes**: None state variables detected

### `MotorPerformanceAnalyzer.calc_motor_losses` (Line 113)
- **Arguments**: self, turns_per_phase, current, rpm
- **Reads (State Dependencies)**:
  - `self.tooth_width`
  - `self.La`
  - `self.poles`
  - `self.awg`
  - `self.get_wire_properties`
- **Writes**: None state variables detected

### `MotorPerformanceAnalyzer.analyze` (Line 141)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.calculate_slot_area`
  - `self.slots`
  - `self.rated_speed`
  - `self.calc_bemf_constants`
  - `self.target_j`
  - `self.stator_id`
  - `self.estimate_max_wires_in_slot`
  - `self.calc_motor_losses`
  - `self.yoke_thickness`
  - `self.awg`
  - `self.get_wire_properties`
- **Writes**: None state variables detected

## File: `machine_core.py`

### `MotorParameters.rOuter` (Line 37)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.stator`
  - `self.stator.OD`
  - `self.stator.yoke`
- **Writes**: None state variables detected

### `MotorParameters.rInner` (Line 41)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.stator`
  - `self.stator.ID`
- **Writes**: None state variables detected

## File: `machine_geometry.py`

### `rotate_point` (Line 11)
- **Arguments**: p, deg
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `AllPoints.__post_init__` (Line 23)
- **Arguments**: self, motor
- **Reads (State Dependencies)**:
  - `self.HP['4']['1']`
  - `rotor.OD`
  - `stator.toothDepth`
  - `self.HP['7']['1']`
  - `self.HP['9']`
  - `self.HP['9']['1']`
  - `stator.toothWidth`
  - `stator.ID`
  - `self.HP['4']['0']`
  - `self.horizontal_position`
  - `self.HP['4']`
  - `self.HP['6']['0']`
  - `rotor.magnetDepth`
  - `self.HP_mirror`
  - `self.HP`
  - `self.slot_count`
  - `self.HP['7']['0']`
  - `self.rotated_position`
  - `stator.OD`
  - `self.RP`
  - `self.HP['9']['0']`
  - `self.HP['6']['1']`
  - `self.HP['7']`
  - `self.pole_count`
  - `self.HP['6']`
  - `stator.yoke`
  - `rotor.ID`
- **Writes (State Changes)**:
  - `self.HP_mirror`
  - `self.HP`
  - `self.rotated_position`
  - `self.slot_count`
  - `self.RP`
  - `self.horizontal_position`
  - `self.pole_count`

### `parse_point_name` (Line 87)
- **Arguments**: name, all_points, rotation_deg
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `draw_instruction_parser` (Line 119)
- **Arguments**: part, all_points, drawer
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `StatorCore.draw_instruction` (Line 211)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points`
  - `self.all_points.slot_count`
- **Writes**: None state variables detected

### `RotorCore.draw_instruction` (Line 230)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points`
  - `self.options`
  - `self.all_points.pole_count`
- **Writes**: None state variables detected

### `Magnet.draw_instruction` (Line 248)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points`
  - `self.all_points.pole_count`
- **Writes**: None state variables detected

### `Coil.draw_instruction` (Line 268)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points`
  - `self.all_points.slot_count`
- **Writes**: None state variables detected

### `MachineGeometry.add_part` (Line 298)
- **Arguments**: self, part
- **Reads (State Dependencies)**:
  - `self.parts.append`
  - `self._next_index`
  - `self.all_points`
  - `self.parts`
- **Writes**: None state variables detected

### `MachineGeometry.sync` (Line 305)
- **Arguments**: self, motor
- **Reads (State Dependencies)**:
  - `self.all_points`
- **Writes (State Changes)**:
  - `self.all_points`

### `MachineGeometry.draw_machine_using_CairoDrawer` (Line 308)
- **Arguments**: self, drawer
- **Reads (State Dependencies)**:
  - `self.all_points`
  - `self.parts`
- **Writes**: None state variables detected

### `MachineGeometry.show_geometry_svg` (Line 324)
- **Arguments**: self, filename, scale
- **Reads (State Dependencies)**:
  - `self.draw_machine_using_CairoDrawer`
- **Writes**: None state variables detected

## File: `machine_geometry_utils.py`

### `rotate_point` (Line 9)
- **Arguments**: p, deg
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `AllPoints.__post_init__` (Line 21)
- **Arguments**: self, user_input
- **Reads (State Dependencies)**:
  - `self.HP['4']['1']`
  - `self.HP['7']['1']`
  - `self.HP['9']`
  - `self.HP['9']['1']`
  - `self.HP['3']['0']`
  - `self.HP['3']['1']`
  - `self.HP['2']['1']`
  - `self.HP['2']`
  - `self.HP['4']['0']`
  - `self.horizontal_position`
  - `self.HP['4']`
  - `self.GP`
  - `self.HP['6']['0']`
  - `self.HP['1']['1']`
  - `self.HP_mirror`
  - `self.HP`
  - `self.slot_count`
  - `self.HP['7']['0']`
  - `self.HP['3']`
  - `self.HP['1']`
  - `self.rotated_position`
  - `self.HP['1']['0']`
  - `self.HP['2']['0']`
  - `self.RP`
  - `self.HP['9']['0']`
  - `self.HP['6']['1']`
  - `self.HP['7']`
  - `self.pole_count`
  - `self.HP['6']`
- **Writes (State Changes)**:
  - `self.HP_mirror`
  - `self.HP`
  - `self.rotated_position`
  - `self.slot_count`
  - `self.RP`
  - `self.horizontal_position`
  - `self.pole_count`
  - `self.GP`

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
  - `self.all_points`
  - `self.options`
  - `self.all_points.pole_count`
- **Writes**: None state variables detected

### `StatorCore.draw_instruction` (Line 269)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.options`
  - `self.all_points`
  - `self.all_points.slot_count`
- **Writes**: None state variables detected

### `Magnet.draw_instruction` (Line 305)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points`
  - `self.all_points.pole_count`
- **Writes**: None state variables detected

### `Coil.draw_instruction` (Line 325)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.all_points`
  - `self.all_points.slot_count`
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
  - `self.winding`
  - `self.winding['l_stack']`
- **Writes**: None state variables detected

### `MachineGeometry.d_air_gap` (Line 372)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['d_air_gap']`
- **Writes**: None state variables detected

### `MachineGeometry.r_stator_outer` (Line 376)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['r_stator_outer']`
- **Writes**: None state variables detected

### `MachineGeometry.r_rotor_outer` (Line 380)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['r_rotor_outer']`
- **Writes**: None state variables detected

### `MachineGeometry.r_shaft` (Line 384)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['r_shaft']`
- **Writes**: None state variables detected

### `MachineGeometry.w_tooth` (Line 388)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['w_tooth']`
- **Writes**: None state variables detected

### `MachineGeometry.d_tooth` (Line 392)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['d_tooth']`
- **Writes**: None state variables detected

### `MachineGeometry.d_magnet` (Line 396)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['d_magnet']`
- **Writes**: None state variables detected

### `MachineGeometry.tooth_shape` (Line 400)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.gp['tooth_shape']`
- **Writes**: None state variables detected

### `MachineGeometry.deg_alpha_rm` (Line 404)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.pole_count`
- **Writes**: None state variables detected

### `MachineGeometry.split_ratio` (Line 409)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.r_stator_outer`
  - `self.r_rotor_outer`
  - `self.d_air_gap`
- **Writes**: None state variables detected

### `MachineGeometry.d_stator_yoke` (Line 413)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.d_tooth`
  - `self.r_stator_outer`
  - `self.r_rotor_outer`
- **Writes**: None state variables detected

### `MachineGeometry.d_tooth_shoe` (Line 417)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.gp`
  - `self.tooth_shape`
  - `self.gp['d_tooth_shoe']`
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
  - `self.all_points`
  - `self.parts`
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
  - `self.gp`
  - `self.winding`
- **Writes (State Changes)**:
  - `self.gp`
  - `self.winding`

## File: `main.py`

### `get_motor_parameters` (Line 26)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `update_motor_parameters` (Line 30)
- **Arguments**: params
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `modern_machine_designer_utility.py`

### `Parameter.__init__` (Line 6)
- **Arguments**: self, name, type, value, bounds, calc, calc_bounds, unit, parameter_dict
- **Reads (State Dependencies)**:
  - `self.calc_bounds`
  - `self.calc`
  - `self.unit`
  - `self.bounds['0']`
  - `self.name`
  - `self.bounds`
  - `self.initialized`
  - `self.parameter_dict`
  - `self.overwritten`
  - `self.value`
  - `self.type`
  - `self.bounds['1']`
- **Writes (State Changes)**:
  - `self.calc_bounds`
  - `self.calc`
  - `self.unit`
  - `self.name`
  - `self.bounds`
  - `self.initialized`
  - `self.parameter_dict`
  - `self.value`
  - `self.type`

### `Parameter.__repr__` (Line 43)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.value`
  - `self.name`
  - `self.unit`
  - `self.type`
- **Writes**: None state variables detected

### `Parameter.sensitivity` (Line 46)
- **Arguments**: self, param_name
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Parameter.to_dict` (Line 49)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.unit`
  - `self.name`
  - `self.bounds`
  - `self.value`
  - `self.type`
- **Writes**: None state variables detected

### `Parameter.from_dict` (Line 68)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Winding.__init__` (Line 92)
- **Arguments**: self, phase_number_m, stator_slot_number_Qs, pole_pair_number_p, suspension_pole_pair_number_ps, coil_pitch_y, bool_DPNVorSEPA, number_of_parallel_branch
- **Reads (State Dependencies)**:
  - `self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.ps`
  - `self.CommutatingSequenceD`
  - `self.layer_Y_phases`
  - `self.number_of_winding_layer`
  - `self.dict_coil_connection`
  - `self.bool_distributed_or_concentrated`
  - `self.bool_CustomizedCircuit`
  - `self.number_of_parallel_branch`
  - `self.layer_Y_signs`
  - `self.SPP`
  - `self.grouping_AC`
  - `self.layer_X_signs`
  - `self.m`
  - `self.bool_3PhaseCurrentSource`
  - `self.bool_DPNVorSEPA`
  - `self.get_winding_factor`
  - `self.layer_X_phases`
  - `self.kw1`
  - `self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.Qs`
  - `self.p`
- **Writes (State Changes)**:
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.ps`
  - `self.CommutatingSequenceD`
  - `self.layer_Y_phases`
  - `self.number_of_winding_layer`
  - `self.number_of_parallel_branch`
  - `self.dict_coil_connection`
  - `self.bool_CustomizedCircuit`
  - `self.layer_Y_signs`
  - `self.SPP`
  - `self.grouping_AC`
  - `self.layer_X_signs`
  - `self.m`
  - `self.bool_3PhaseCurrentSource`
  - `self.bool_DPNVorSEPA`
  - `self.layer_X_phases`
  - `self.kw1`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.Qs`
  - `self.p`

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
  - `self.ps`
  - `self.m`
  - `self.coil_pitch_y`
  - `self.Qs`
  - `self.p`
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
  - `self.kw1`
  - `self.ps`
  - `self.m`
  - `self.derivation`
  - `self.p`
  - `self.derivation.__dict__.items`
  - `self.coil_pitch_y`
  - `self.Qs`
  - `self.number_of_parallel_branch`
  - `self.derivation.__dict__`
- **Writes**: None state variables detected

### `Winding.from_dict` (Line 317)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Geometry.__init__` (Line 362)
- **Arguments**: self, name, GP, draw_function, color
- **Reads (State Dependencies)**:
  - `self.GP.items`
  - `self.color`
  - `self.name`
  - `self.draw_function`
  - `self.GP`
- **Writes (State Changes)**:
  - `self.draw_function`
  - `self.color`
  - `self.GP`
  - `self.name`

### `Geometry.update_from_GP` (Line 371)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.GP.items`
  - `self.GP`
- **Writes**: None state variables detected

### `Geometry.__repr__` (Line 382)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.GP.items`
  - `self.name`
  - `self.color`
  - `self.GP`
- **Writes**: None state variables detected

### `Geometry.print_parameters` (Line 393)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.GP.items`
  - `self.color`
  - `self.name`
  - `self.__dict__.items`
  - `self.visualization_points.items`
  - `self.GP`
  - `self.visualization_points`
  - `self.__dict__`
- **Writes**: None state variables detected

### `Geometry.to_dict` (Line 439)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.__dict__.items`
  - `self.__dict__`
- **Writes**: None state variables detected

### `Geometry.from_dict` (Line 484)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Geometry.draw` (Line 490)
- **Arguments**: self, drawer
- **Reads (State Dependencies)**:
  - `self.components_make_region`
  - `self.draw_function`
  - `self.visualization_points`
- **Writes (State Changes)**:
  - `self.components_make_region`
  - `self.visualization_points`

### `CairoDrawer.__init__` (Line 514)
- **Arguments**: self, width_in_points, height_in_points, filename, verbose_drawing, scale, bFillRegion
- **Reads (State Dependencies)**:
  - `self.ctx.transform`
  - `self.surface`
  - `self.bMirror`
  - `self.ctx.restore`
  - `self.ctx.paint`
  - `self.ctx.save`
  - `self.verbose_drawing`
  - `self.ctx.set_source_rgb`
  - `self.bFillRegion`
  - `self.filename`
  - `self.scale`
  - `self.ctx.scale`
  - `self.iRotateCopy`
  - `self.ctx`
- **Writes (State Changes)**:
  - `self.surface`
  - `self.bMirror`
  - `self.verbose_drawing`
  - `self.filename`
  - `self.bFillRegion`
  - `self.scale`
  - `self.iRotateCopy`
  - `self.ctx`

### `CairoDrawer.apply_stroke` (Line 537)
- **Arguments**: self, lw
- **Reads (State Dependencies)**:
  - `self.ctx.stroke`
  - `self.ctx.set_line_cap`
  - `self.ctx.set_source_rgba`
  - `self.ctx.set_line_width`
  - `self.ctx`
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
  - `self.ctx.stroke`
  - `self.ctx.arc`
  - `self.ctx.move_to`
  - `self.ctx.restore`
  - `self.ctx.arc_negative`
  - `self.ctx.save`
  - `self.ctx.set_line_width`
  - `self.ctx.fill_preserve`
  - `self.ctx.set_source_rgba`
  - `self.ctx.rotate`
  - `self.hex_to_rgb`
  - `self.bFillRegion`
  - `self.scale`
  - `self.ctx.scale`
  - `self.ctx.path_extents`
  - `self.ctx.line_to`
  - `self.ctx`
- **Writes**: None state variables detected

### `CairoDrawer.finalize_part` (Line 675)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.__init__` (Line 680)
- **Arguments**: self, specs
- **Reads (State Dependencies)**:
  - `self.winding`
  - `self.specs`
  - `self.flag_do_not_evaluate_when_init_pop`
  - `self.geometry`
  - `self.target`
- **Writes (State Changes)**:
  - `self.winding`
  - `self.specs`
  - `self.flag_do_not_evaluate_when_init_pop`
  - `self.geometry`
  - `self.target`

### `Modern_Machine_Designer_Utility._get_parameter_logger` (Line 690)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.name`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.apply_parameter_dict` (Line 711)
- **Arguments**: self, prev_params, key_map
- **Reads (State Dependencies)**:
  - `self.get_parameters_by_type`
  - `self.get_parameter_dict_by_name`
  - `self.parameter_dict_by_name.get`
  - `self.geometry.machineGeometry.values`
  - `self._get_parameter_logger`
  - `self.geometry.machineGeometry`
  - `self.parameter_dict_by_name`
  - `self.geometry`
- **Writes (State Changes)**:
  - `self.parameter_dict_by_name`

### `Modern_Machine_Designer_Utility.update_geometric_parameters` (Line 781)
- **Arguments**: self, x_denorm, x_denorm_dict
- **Reads (State Dependencies)**:
  - `self.get_parameters_by_type`
  - `self.get_parameter_dict_by_name`
  - `self.geometry.machineGeometry.items`
  - `self.geometry.machineGeometry`
  - `self.parameter_dict_by_name`
  - `self.geometry`
  - `self.get_free_variables`
- **Writes (State Changes)**:
  - `self.parameter_dict_by_name`

### `Modern_Machine_Designer_Utility.get_pc_name` (Line 809)
- **Arguments**: None
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.learn_about_the_archive` (Line 822)
- **Arguments**: self, prob, swarm_data, popsize, bool_plot_and_show, bool_more_info
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.write_swarm_survivor` (Line 901)
- **Arguments**: self, pop, counter_fitness_return
- **Reads (State Dependencies)**:
  - `self.path2SwarmData`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_bad_fintess_values` (Line 917)
- **Arguments**: self, machine_type, ref
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_rotor_volume` (Line 941)
- **Arguments**: self, stack_length
- **Reads (State Dependencies)**:
  - `self.mm_r_ro.value`
  - `self.mm_r_ro`
  - `self.winding`
  - `self.winding.EX['mm_stack_length_specified']`
  - `self.winding.EX`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_rotor_weight` (Line 949)
- **Arguments**: self, gravity, stack_length
- **Reads (State Dependencies)**:
  - `self.get_rotor_volume`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_free_variables` (Line 966)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_free_variables_as_dict` (Line 969)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_free_variables`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_free_variable_bounds_dict` (Line 980)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_free_variables`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.set_free_variables_from_dict` (Line 992)
- **Arguments**: self, free_variables_dict
- **Reads (State Dependencies)**:
  - `self.get_parameter`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.update_derived_parameters` (Line 999)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter_fields` (Line 1017)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameters_by_type` (Line 1032)
- **Arguments**: self, param_type
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter` (Line 1045)
- **Arguments**: self, name
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.set_parameter_value` (Line 1058)
- **Arguments**: self, name, value
- **Reads (State Dependencies)**:
  - `self.get_parameter`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter_dict` (Line 1075)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameter_dict_by_name` (Line 1084)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.list_parameters` (Line 1093)
- **Arguments**: self, param_type
- **Reads (State Dependencies)**:
  - `self.get_parameters_by_type`
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.get_parameters_summary` (Line 1107)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.validate_parameters` (Line 1136)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.get_parameter_fields`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.__repr__` (Line 1168)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.target.machine_class`
  - `self.get_parameters_summary`
  - `self.target`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.to_dict` (Line 1177)
- **Arguments**: self
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.to_json` (Line 1219)
- **Arguments**: self, indent, ensure_ascii
- **Reads (State Dependencies)**:
  - `self.to_dict`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.save_to_file` (Line 1232)
- **Arguments**: self, filepath, indent
- **Reads (State Dependencies)**:
  - `self.to_dict`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.to_dict_full` (Line 1243)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.__class__.__module__`
  - `self.bool_jmagDeleteResultsAfterCalculation`
  - `self.geometry.machineGeometry.items`
  - `self.winding.EX`
  - `self.winding.EX.copy`
  - `self.geometry`
  - `self.__class__`
  - `self.winding.wily`
  - `self.target`
  - `self.target.machine_class`
  - `self.winding`
  - `self.select_fea_config_dict`
  - `self.select_FEA_tool`
  - `self.get_parameter_fields`
  - `self.__class__.__name__`
  - `self.bool_RotorNotched`
  - `self.counter`
  - `self.bool_PermanentMagnet`
  - `self.name`
  - `self.winding.wily.to_dict`
  - `self.geometry.machineGeometry`
  - `self.bool_StatorSlotClosed`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.save_to_file_full` (Line 1399)
- **Arguments**: self, filepath, indent
- **Reads (State Dependencies)**:
  - `self.to_dict_full`
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.from_dict_full` (Line 1411)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.load_from_file_full` (Line 1552)
- **Arguments**: cls, filepath
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.from_dict` (Line 1567)
- **Arguments**: cls, data
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.from_json` (Line 1639)
- **Arguments**: cls, json_str
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.load_from_file` (Line 1653)
- **Arguments**: cls, filepath
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.remove_jfiles_folders` (Line 1669)
- **Arguments**: root_dir
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.myLogger` (Line 1691)
- **Arguments**: dir_log, prefix
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Modern_Machine_Designer_Utility.draw_individual_from_swarm` (Line 1720)
- **Arguments**: self, index
- **Reads (State Dependencies)**:
  - `self.path2SwarmData`
  - `self.show_geometry`
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.decode_py_reduce_ordered_dict` (Line 1764)
- **Arguments**: x_denorm_dict_raw
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.__init__` (Line 1807)
- **Arguments**: self, fname, desired_x_denorm_dict, bool_filter_pareto_front
- **Reads (State Dependencies)**:
  - `self.swarm_data_as_dict`
  - `self.swarm_data_xf`
  - `self.swarm_data_xf.append`
  - `self.decode_py_reduce_ordered_dict`
  - `self.number_of_free_variables`
  - `self.swarm_data_project_names`
  - `self.swarm_data_xf['0']`
  - `self.number_of_chromosome`
  - `self.filter_data`
  - `self.get_metric_of_the_whole_swarm`
- **Writes (State Changes)**:
  - `self.swarm_data_as_dict`
  - `self.swarm_data_xf`
  - `self.number_of_free_variables`
  - `self.swarm_data_project_names`
  - `self.number_of_chromosome`

### `Swarm_Data_Analyzer.filter_data` (Line 1976)
- **Arguments**: self, data, param_type, filter_key, direction, filter_value
- **Reads (State Dependencies)**:
  - `self.decode_py_reduce_ordered_dict`
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.decode` (Line 1999)
- **Arguments**: d
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.get_metric_of_the_whole_swarm` (Line 2005)
- **Arguments**: self, metric
- **Reads (State Dependencies)**:
  - `self.swarm_data_as_dict`
  - `self.swarm_data_as_dict.items`
- **Writes**: None state variables detected

### `Swarm_Data_Analyzer.prepare_data_for_post_processing` (Line 2025)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.rotor_weight`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.Ea`
  - `self.Cost_Fe`
  - `self.f1`
  - `self.f2`
  - `self.l_rated_stack_length`
  - `self.get_metric_of_the_whole_swarm`
  - `self.torque_average`
  - `self.TorqueRipple`
  - `self.Tripple`
  - `self.ss_avg_force_magnitude`
  - `self.TRV`
  - `self.RatedEfficiency`
  - `self.force_error_angle`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.FRW`
  - `self.normalized_force_error_magnitude`
  - `self.f3`
  - `self.RatedStkLen`
  - `self.PowerFactor`
  - `self.Em`
  - `self.swarm_data_xf`
  - `self.l_rated_iron_loss`
  - `self.Cost_PM`
  - `self.l_rated_magnet_Joule_loss`
  - `self.Cost`
  - `self.Cost_Cu`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_rated_stator_copper_loss_along_stack`
- **Writes (State Changes)**:
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.rotor_weight`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.Ea`
  - `self.Cost_Fe`
  - `self.f1`
  - `self.f2`
  - `self.l_rated_stack_length`
  - `self.torque_average`
  - `self.TorqueRipple`
  - `self.Tripple`
  - `self.ss_avg_force_magnitude`
  - `self.TRV`
  - `self.RatedEfficiency`
  - `self.force_error_angle`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.FRW`
  - `self.normalized_force_error_magnitude`
  - `self.f3`
  - `self.RatedStkLen`
  - `self.PowerFactor`
  - `self.Em`
  - `self.l_rated_iron_loss`
  - `self.Cost_PM`
  - `self.l_rated_magnet_Joule_loss`
  - `self.Cost`
  - `self.Cost_Cu`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_rated_stator_copper_loss_along_stack`

### `swarm_data_container.__init__` (Line 2111)
- **Arguments**: self, swarm_data_raw, fea_config_dict, swarm_data_json, swarm_data_json_file_path
- **Reads (State Dependencies)**:
  - `self.swarm_data_raw`
  - `self._initialize_empty`
  - `self.fea_config_dict`
  - `self._load_from_raw`
  - `self._load_from_json`
- **Writes (State Changes)**:
  - `self.fea_config_dict`
  - `self.swarm_data_raw`

### `swarm_data_container._initialize_empty` (Line 2141)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.RatedWeight`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.Ea`
  - `self.mm_w_st`
  - `self.l_original_stack_length`
  - `self.l_rated_stack_length`
  - `self.l_rated_rotor_weight`
  - `self.deg_alpha_st`
  - `self.l_OB`
  - `self.l_rated_efficiency`
  - `self.l_TRV`
  - `self.l_normalized_torque_ripple`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.FRW`
  - `self.rated_data`
  - `self.l_original_rotor_weight`
  - `self.l_ss_avg_force_magnitude`
  - `self.RatedStkLen`
  - `self.l_force_error_angle`
  - `self.Em`
  - `self.RatedVol`
  - `self.swarm_data_xf`
  - `self.l_rated_iron_loss`
  - `self.machine_data`
  - `self.number_of_free_variables`
  - `self.l_efficiency`
  - `self.Trip`
  - `self.mm_r_si`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_OA`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_OC`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_FRW`
  - `self.project_names`
  - `self.l_design_parameters`
- **Writes (State Changes)**:
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.RatedWeight`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.Ea`
  - `self.mm_w_st`
  - `self.l_original_stack_length`
  - `self.l_rated_stack_length`
  - `self.l_rated_rotor_weight`
  - `self.deg_alpha_st`
  - `self.l_OB`
  - `self.l_rated_efficiency`
  - `self.l_TRV`
  - `self.l_normalized_torque_ripple`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.FRW`
  - `self.rated_data`
  - `self.l_original_rotor_weight`
  - `self.l_ss_avg_force_magnitude`
  - `self.RatedStkLen`
  - `self.l_force_error_angle`
  - `self.Em`
  - `self.RatedVol`
  - `self.swarm_data_xf`
  - `self.l_rated_iron_loss`
  - `self.machine_data`
  - `self.number_of_free_variables`
  - `self.l_efficiency`
  - `self.Trip`
  - `self.mm_r_si`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_OA`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_OC`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_FRW`
  - `self.project_names`
  - `self.l_design_parameters`

### `swarm_data_container._load_from_json` (Line 2188)
- **Arguments**: self, swarm_data_json
- **Reads (State Dependencies)**:
  - `self.Ea.append`
  - `self.RatedWeight`
  - `self.swarm_data_xf.append`
  - `self.Ea`
  - `self._extract_performance_lists`
  - `self.Trip.append`
  - `self.project_names.append`
  - `self.machine_data.append`
  - `self.RatedVol.append`
  - `self.FRW`
  - `self.rated_data`
  - `self.RatedStkLen`
  - `self.FRW.append`
  - `self.Em`
  - `self.RatedStkLen.append`
  - `self.RatedVol`
  - `self.swarm_data_xf`
  - `self.machine_data`
  - `self.Em.append`
  - `self._initialize_empty`
  - `self.number_of_free_variables`
  - `self.Trip`
  - `self.rated_data.append`
  - `self.RatedWeight.append`
  - `self.swarm_data_xf['0']`
  - `self.project_names`
- **Writes (State Changes)**:
  - `self.number_of_free_variables`

### `swarm_data_container._extract_performance_lists` (Line 2302)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_original_stack_length`
  - `self.l_rated_stack_length`
  - `self.l_rated_rotor_weight`
  - `self.l_OB`
  - `self.l_rated_efficiency`
  - `self.l_TRV`
  - `self.l_normalized_torque_ripple`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.rated_data`
  - `self.l_original_rotor_weight`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_force_error_angle`
  - `self.swarm_data_xf`
  - `self.l_rated_iron_loss`
  - `self.machine_data`
  - `self.l_efficiency`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_OA`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_OC`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_FRW`
  - `self.l_design_parameters`
- **Writes (State Changes)**:
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.l_original_stack_length`
  - `self.l_rated_stack_length`
  - `self.l_rated_rotor_weight`
  - `self.l_OB`
  - `self.l_rated_efficiency`
  - `self.l_TRV`
  - `self.l_normalized_torque_ripple`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.l_original_rotor_weight`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_force_error_angle`
  - `self.l_rated_iron_loss`
  - `self.l_efficiency`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_OA`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_OC`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_FRW`
  - `self.l_design_parameters`

### `swarm_data_container._load_from_raw` (Line 2347)
- **Arguments**: self, swarm_data_raw
- **Reads (State Dependencies)**:
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.Ea.append`
  - `self.RatedWeight`
  - `self.swarm_data_xf.append`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.Ea`
  - `self.mm_r_si.append`
  - `self._extract_performance_lists`
  - `self.mm_w_st`
  - `self.Trip.append`
  - `self.l_original_stack_length`
  - `self.l_rated_stack_length`
  - `self.project_names.append`
  - `self.l_rated_rotor_weight`
  - `self.deg_alpha_st`
  - `self.l_OB`
  - `self.machine_data.append`
  - `self.l_rated_efficiency`
  - `self.l_TRV`
  - `self.RatedVol.append`
  - `self.l_normalized_torque_ripple`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.FRW`
  - `self.rated_data`
  - `self.l_original_rotor_weight`
  - `self.l_ss_avg_force_magnitude`
  - `self.RatedStkLen`
  - `self.FRW.append`
  - `self.l_force_error_angle`
  - `self.Em`
  - `self.RatedStkLen.append`
  - `self.deg_alpha_st.append`
  - `self.RatedVol`
  - `self.swarm_data_xf`
  - `self.l_rated_iron_loss`
  - `self.machine_data`
  - `self.Em.append`
  - `self._initialize_empty`
  - `self.number_of_free_variables`
  - `self.l_efficiency`
  - `self.mm_w_st.append`
  - `self.Trip`
  - `self.mm_r_si`
  - `self.rated_data.append`
  - `self.RatedWeight.append`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_OA`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_OC`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.swarm_data_xf['0']`
  - `self.l_FRW`
  - `self.project_names`
  - `self.l_design_parameters`
- **Writes (State Changes)**:
  - `self.l_normalized_force_error_magnitude`
  - `self.l_rated_rotor_copper_loss_along_stack`
  - `self.l_rated_total_loss`
  - `self.l_rated_windage_loss`
  - `self.mm_w_st`
  - `self.l_original_stack_length`
  - `self.l_rated_stack_length`
  - `self.l_rated_rotor_weight`
  - `self.deg_alpha_st`
  - `self.l_OB`
  - `self.l_rated_efficiency`
  - `self.l_TRV`
  - `self.l_normalized_torque_ripple`
  - `self.l_stator_copper_loss_in_end_turn`
  - `self.l_torque_average`
  - `self.l_original_rotor_weight`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_force_error_angle`
  - `self.l_rated_iron_loss`
  - `self.number_of_free_variables`
  - `self.l_efficiency`
  - `self.mm_r_si`
  - `self.l_rotor_copper_loss_in_end_turn`
  - `self.l_OA`
  - `self.l_power_factor`
  - `self.l_rated_shaft_power`
  - `self.l_OC`
  - `self.l_rated_rotor_volume`
  - `self.l_rated_stator_copper_loss_along_stack`
  - `self.l_FRW`
  - `self.l_design_parameters`

### `swarm_data_container.get_list_y_data` (Line 2580)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.l_OB`
  - `self.l_TRV`
  - `self.l_force_error_angle`
- **Writes**: None state variables detected

### `swarm_data_container.sensitivity_bar_charts` (Line 2593)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.reference_design`
  - `self.l_normalized_force_error_magnitude`
  - `self.rotor_weight`
  - `self.reference_data`
  - `self.reference_design['3'].split`
  - `self.l_rated_total_loss`
  - `self.reference_data['3']`
  - `self.get_certain_objective_function`
  - `self.reference_data['6']`
  - `self.l_OB`
  - `self.l_normalized_torque_ripple`
  - `self.required_torque`
  - `self.fea_config_dict`
  - `self.l_ss_avg_force_magnitude`
  - `self.l_original_rotor_weight`
  - `self.l_force_error_angle`
  - `self.reference_design['3']`
  - `self.number_of_free_variables`
  - `self.rotor_volume`
  - `self.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
  - `self.reference_data['5']`
  - `self.reference_data['2']`
  - `self.reference_data['4']`
  - `self.l_OA`
  - `self.l_OC`
- **Writes (State Changes)**:
  - `self.reference_data`

## File: `utility.py`

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
  - `self.message`
  - `self.payload`
- **Writes (State Changes)**:
  - `self.message`
  - `self.payload`

### `ExceptionReTry.__str__` (Line 46)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.message`
- **Writes**: None state variables detected

### `ExceptionBadNumberOfParts.__init__` (Line 51)
- **Arguments**: self, message, payload
- **Reads (State Dependencies)**:
  - `self.message`
  - `self.payload`
- **Writes (State Changes)**:
  - `self.message`
  - `self.payload`

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
  - `self.rotor_tooth_width_b_dr`
  - `self.SIesign_parameters_denorm`
  - `self.Angle_StatorSlotOpen`
  - `self.air_gap_length_delta`
  - `self.show_norm`
  - `self.Width_StatorTeethHeadThickness`
  - `self.stator_tooth_width_b_ds`
  - `self.b1`
  - `self.Length_HeadNeckRotorSlot`
- **Writes (State Changes)**:
  - `self.rotor_tooth_width_b_dr`
  - `self.SIesign_parameters_denorm`
  - `self.Angle_StatorSlotOpen`
  - `self.air_gap_length_delta`
  - `self.Width_StatorTeethHeadThickness`
  - `self.stator_tooth_width_b_ds`
  - `self.b1`
  - `self.Length_HeadNeckRotorSlot`

### `Pyrhonen_design.show_denorm` (Line 375)
- **Arguments**: self, bounds, design_parameters_norm
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `Pyrhonen_design.show_norm` (Line 383)
- **Arguments**: self, bounds, design_parameters_denorm
- **Reads (State Dependencies)**:
  - `self.SIesign_parameters_norm.tolist`
  - `self.SIesign_parameters_norm`
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
  - `self.force_x`
  - `self.force_err_abs`
  - `self.force_err_ang_old_way`
  - `self.ss_max_force_err_abs['1']`
  - `self.ss_max_force_err_abs['0']`
  - `self.force_y`
  - `self.ss_avg_force_magnitude`
  - `self.ss_max_force_err_ang['1']`
  - `self.force_ang.append`
  - `self.force_error_angle`
  - `self.normalized_force_error_magnitude`
  - `self.ss_avg_force_vector`
  - `self.ss_max_force_err_ang`
  - `self.ss_max_force_err_ang['0']`
  - `self.force_abs`
  - `self.force_err_ang`
  - `self.ss_avg_force_angle`
  - `self.ss_avg_force_vector['1']`
  - `self.ss_avg_force_vector['0']`
  - `self.force_err_ang_new_way`
  - `self.force_ang`
  - `self.range_ss`
  - `self.ss_max_force_err_abs`
- **Writes (State Changes)**:
  - `self.force_err_ang_new_way`
  - `self.force_x`
  - `self.ss_avg_force_magnitude`
  - `self.normalized_force_error_magnitude`
  - `self.force_abs`
  - `self.force_err_abs`
  - `self.force_err_ang`
  - `self.force_err_ang_old_way`
  - `self.range_ss`
  - `self.force_error_angle`
  - `self.ss_avg_force_angle`
  - `self.force_ang`
  - `self.force_y`
  - `self.ss_avg_force_vector`
  - `self.ss_max_force_err_ang`
  - `self.ss_max_force_err_abs`

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
  - `self.im.template.SI['GP']['mm_r_ro']`
  - `self.im.template.SI['GP']['mm_r_ro'].value`
  - `self.im.rotor_slot_height_h_sr`
  - `EX['WindingFill']`
  - `EX['wily'].number_parallel_branch`
  - `EX['DriveW_zQ']`
  - `self.im`
  - `self.im.design_parameters['2']`
  - `self.list_rotor_current_amp`
  - `self.im.template.SI`
  - `EX['Js']`
  - `self.im.stack_length`
  - `self.im.template.SI['GP']`
  - `self.im.template`
  - `self.im.DriveW_poles`
  - `EX['wily']`
  - `self.im.design_parameters`
  - `self.im.Qr`
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
  - `self.count`
  - `self.imag`
  - `self.sine`
  - `self.k`
  - `self.q2`
  - `self.phase`
  - `self.q`
  - `self.id`
  - `self.real`
  - `self.coeff`
  - `self.scalingFactor`
  - `self.accumSquaredData`
- **Writes (State Changes)**:
  - `self.ampl`
  - `self.bool_initialized`
  - `self.cosine`
  - `self.count`
  - `self.imag`
  - `self.sine`
  - `self.k`
  - `self.q2`
  - `self.phase`
  - `self.q`
  - `self.id`
  - `self.real`
  - `self.coeff`
  - `self.scalingFactor`
  - `self.accumSquaredData`

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
  - `self.reference_design`
  - `self.sw`
  - `self.build_basic_info`
  - `self.SIir_run`
  - `self.run_integer`
  - `self.spec`
  - `self.buf`
  - `self.buf_length`
  - `self.number_of_designs`
- **Writes (State Changes)**:
  - `self.reference_design`
  - `self.sw`
  - `self.SIir_run`
  - `self.run_integer`
  - `self.spec`
  - `self.buf`
  - `self.buf_length`
  - `self.number_of_designs`

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
  - `self.rotor_weight`
  - `self.sw.im.template.SI`
  - `self.Qr`
  - `self.mec_power`
  - `self.spec.Jr`
  - `self.best_design_display`
  - `self.weights_used`
  - `self.stack_length`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.speed_rpm`
  - `self.spec.stator_phase_current_rms`
  - `self.sw.im.BeariW_poles`
  - `self.sw.im.template`
  - `self.template`
  - `self.template.SI`
  - `self.required_torque`
  - `self.weights_name`
  - `self.template.SI['GP']['mm_r_ro']`
  - `self.SIesign_parameters_generator`
  - `self.Omega`
  - `self.sw.fea_config_dict`
  - `self.SIesign_display_generator`
  - `self.best_design_display.split`
  - `self.str_best_design_details`
  - `self.sw.im.template.SI['GP']`
  - `self.spec.Js`
  - `self.run_integer`
  - `self.spec`
  - `self.sw.im.DriveW_poles`
  - `self.sw.im.template.SI['GP']['mm_r_ro']`
  - `self.rotor_volume`
  - `self.template.SI['GP']`
  - `self.best_design_denorm['0']`
  - `self.sw`
  - `self.sw.fea_config_dict['use_weights']`
  - `self.spec.Steel`
  - `self.ExcitationFreqSimulated`
  - `self.spec.VoltageRating`
  - `self.sw.im`
  - `self.sw.im.Radius_OuterStatorYoke`
  - `self.best_design_denorm`
  - `self.Qs`
  - `self.sw.im.template.SI['GP']['mm_r_ro'].value`
  - `self.stack_length_max`
- **Writes (State Changes)**:
  - `self.str_best_design_details`
  - `self.best_design_denorm`
  - `self.best_design_display`

### `SwarmDataAnalyzer.pareto_plot_eta_vs_stack_legnth` (Line 1560)
- **Arguments**: self, fig, ax, marker, bool_filtered
- **Reads (State Dependencies)**:
  - `self.rotor_weight`
  - `self.mec_power`
  - `self.get_certain_objective_function`
  - `self.required_torque`
  - `self.weights_used`
  - `self.my_scatter_plot`
  - `self.rotor_volume`
  - `self.stack_length_max`
  - `self.stack_length`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.pareto_plot_torque_force` (Line 1648)
- **Arguments**: self, fig2, axeses, marker
- **Reads (State Dependencies)**:
  - `self.rotor_weight`
  - `self.get_certain_objective_function`
  - `self.weights_name`
  - `self.weights_used`
  - `self.my_scatter_plot`
  - `self.rotor_volume`
- **Writes**: None state variables detected

### `SwarmDataAnalyzer.sensitivity_bar_charts` (Line 1765)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.reference_design`
  - `self.sw`
  - `self.rotor_weight`
  - `self.reference_data['5']`
  - `self.reference_data`
  - `self.reference_design['3'].split`
  - `self.reference_data['2']`
  - `self.reference_data['4']`
  - `self.reference_data['3']`
  - `self.get_certain_objective_function`
  - `self.reference_design['3']`
  - `self.sw.fea_config_dict['local_sensitivity_analysis_number_of_variants']`
  - `self.required_torque`
  - `self.reference_data['6']`
  - `self.sw.fea_config_dict`
  - `self.rotor_volume`
- **Writes (State Changes)**:
  - `self.reference_data`

### `SwarmDataAnalyzer.build_basic_info` (Line 2127)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.rotor_weight`
  - `self.Qr`
  - `self.mec_power`
  - `self.weights_used`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.stack_length`
  - `self.speed_rpm`
  - `self.template`
  - `self.template.SI`
  - `self.required_torque`
  - `self.weights_name`
  - `self.template.SI['GP']['mm_r_ro']`
  - `self.Omega`
  - `self.spec`
  - `self.template.SI['GP']`
  - `self.rotor_volume`
  - `self.sw`
  - `self.ExcitationFreqSimulated`
  - `self.Qs`
  - `self.stack_length_max`
- **Writes (State Changes)**:
  - `self.speed_rpm`
  - `self.rotor_weight`
  - `self.template.SI['GP']['mm_r_ro'].value`
  - `self.Qr`
  - `self.mec_power`
  - `self.ExcitationFreqSimulated`
  - `self.required_torque`
  - `self.weights_name`
  - `self.weights_used`
  - `self.Qs`
  - `self.Omega`
  - `self.rotor_volume`
  - `self.stack_length_max`
  - `self.stack_length`

### `build_sensitivity_bar_charts` (Line 2159)
- **Arguments**: spec, sw
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `build_Pareto_plot` (Line 2164)
- **Arguments**: spec, sw
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `winding_layout.py`

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
  - `self.kd1`
  - `self.kp1`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.ps`
  - `self.CommutatingSequenceD`
  - `self.layer_Y_phases`
  - `self.dict_coil_connection`
  - `self.list_layer_motor_signs`
  - `self.number_winding_layer`
  - `self.layer_Y_signs`
  - `self.number_parallel_branch`
  - `self.SIPNV_or_SEPA`
  - `self.SPP`
  - `self.grouping_AC`
  - `self.layer_X_signs`
  - `self.distributed_or_concentrated`
  - `self.m`
  - `self.pr`
  - `self.bool_3PhaseCurrentSource`
  - `self.list_layer_motor_phases`
  - `self.ox_distribution_three_phase`
  - `self.ox_distribution_phase_U`
  - `self.layer_X_phases`
  - `self.list_layer_suspension_phases`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.Qs`
  - `self.p`
  - `self.list_layer_suspension_signs`
- **Writes (State Changes)**:
  - `self.kd1`
  - `self.kp1`
  - `self.deg_winding_U_phase_phase_axis_angle`
  - `self.ps`
  - `self.CommutatingSequenceD`
  - `self.layer_Y_phases`
  - `self.dict_coil_connection`
  - `self.list_layer_motor_signs`
  - `self.number_winding_layer`
  - `self.layer_Y_signs`
  - `self.number_parallel_branch`
  - `self.SIPNV_or_SEPA`
  - `self.SPP`
  - `self.grouping_AC`
  - `self.layer_X_signs`
  - `self.distributed_or_concentrated`
  - `self.m`
  - `self.pr`
  - `self.bool_3PhaseCurrentSource`
  - `self.list_layer_motor_phases`
  - `self.ox_distribution_three_phase`
  - `self.ox_distribution_phase_U`
  - `self.layer_X_phases`
  - `self.list_layer_suspension_phases`
  - `self.CommutatingSequenceB`
  - `self.coil_pitch_y`
  - `self.Qs`
  - `self.p`
  - `self.list_layer_suspension_signs`

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
  - `self.sym_winding_func`
  - `self.turns_per_slot`
  - `self.avg_val_of_turn_func`
  - `self.turn_func`
  - `self.degree_between_slots`
  - `self.ox_distribution_phase_U`
  - `self.sym_turn_func`
  - `self.slot_per_phase`
  - `self.sym_begin_pos`
  - `self.radian_between_slots`
  - `self.setTurnFuncObject`
  - `self.winding_func`
  - `self.setSymPos`
- **Writes (State Changes)**:
  - `self.sym_winding_func`
  - `self.turns_per_slot`
  - `self.avg_val_of_turn_func`
  - `self.degree_between_slots`
  - `self.ox_distribution_phase_U`
  - `self.sym_turn_func`
  - `self.slot_per_phase`
  - `self.radian_between_slots`
  - `self.winding_func`

### `PhaseWinding.setTurnFuncObject` (Line 1346)
- **Arguments**: self, ox_distribution_phase_U
- **Reads (State Dependencies)**:
  - `self.lst_y`
  - `self.turns_per_slot`
  - `self.lst_x`
  - `self.radian_between_slots`
  - `self.turn_func`
- **Writes (State Changes)**:
  - `self.lst_y`
  - `self.lst_x`
  - `self.turn_func`

### `PhaseWinding.setSymPos` (Line 1378)
- **Arguments**: self, index
- **Reads (State Dependencies)**:
  - `self.lst_y`
  - `self.lst_y.index`
  - `self.sym_begin_pos_2`
  - `self.sym_begin_pos_1`
  - `self.sym_begin_pos`
  - `self.lst_x`
- **Writes (State Changes)**:
  - `self.sym_begin_pos_2`
  - `self.sym_begin_pos_1`
  - `self.sym_begin_pos`

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
  - `self.Qs`
  - `self.coil_pitch`
  - `self.CommutatingSequenceD`
  - `self.l22`
  - `self.l_rightlayer2`
  - `self.l_leftlayer2`
  - `self.layer_B1`
  - `self.number_parallel_branch`
  - `self.grouping_AC`
  - `self.initial_excitation_bias_compensation_deg`
  - `self.distributed_or_concentrated`
  - `self.l41`
  - `self.layer_B2`
  - `self.bool_3PhaseCurrentSource`
  - `self.no_winding_layer`
  - `self.layer_A2`
  - `self.l_rightlayer1`
  - `self.layer_A1`
  - `self.l21`
  - `self.l42`
  - `self.CommutatingSequenceB`
  - `self.l_leftlayer1`
  - `self.p`
- **Writes (State Changes)**:
  - `self.Qs`
  - `self.coil_pitch`
  - `self.CommutatingSequenceD`
  - `self.l22`
  - `self.l_rightlayer2`
  - `self.l_leftlayer2`
  - `self.layer_B1`
  - `self.number_parallel_branch`
  - `self.grouping_AC`
  - `self.initial_excitation_bias_compensation_deg`
  - `self.distributed_or_concentrated`
  - `self.l41`
  - `self.layer_B2`
  - `self.bool_3PhaseCurrentSource`
  - `self.no_winding_layer`
  - `self.layer_A2`
  - `self.l_rightlayer1`
  - `self.layer_A1`
  - `self.l21`
  - `self.l42`
  - `self.CommutatingSequenceB`
  - `self.l_leftlayer1`
  - `self.p`

## File: `winding_layout_derivation_ismb2021_asymetry_no_drawing.py`

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
  - `self.connection_star_raw_dict`
  - `self.ps`
  - `self.bool_double_layer_winding`
  - `self.suspen_kp_at_h`
  - `self.t`
  - `self.suspen_kd_at_h`
  - `self.Q`
  - `self.verbose`
  - `self.turn_func_bias`
  - `self.torque_kp_at_h`
  - `self.m`
  - `self.torque_kw_at_h`
  - `self.ts`
  - `self.suspen_kw_at_h`
  - `self.dpnv_grouping_dict_a`
  - `self.list_phase_u_slot_number`
  - `self.list_phase_w_slot_number`
  - `self.q`
  - `self.list_phase_v_slot_number`
  - `self.list_slot_number_of_phase`
  - `self.dpnv_grouping_dict_b`
  - `self.qs`
  - `self.torque_kd_at_h`
  - `self.coil_pitch_y`
  - `self.dpnv_grouping_dict_c`
  - `self.p`
- **Writes (State Changes)**:
  - `self.connection_star_raw_dict`
  - `self.ps`
  - `self.bool_double_layer_winding`
  - `self.suspen_kp_at_h`
  - `self.t`
  - `self.suspen_kd_at_h`
  - `self.Q`
  - `self.verbose`
  - `self.turn_func_bias`
  - `self.torque_kp_at_h`
  - `self.m`
  - `self.torque_kw_at_h`
  - `self.ts`
  - `self.suspen_kw_at_h`
  - `self.dpnv_grouping_dict_a`
  - `self.list_phase_u_slot_number`
  - `self.list_phase_w_slot_number`
  - `self.q`
  - `self.list_phase_v_slot_number`
  - `self.list_slot_number_of_phase`
  - `self.dpnv_grouping_dict_b`
  - `self.qs`
  - `self.torque_kd_at_h`
  - `self.coil_pitch_y`
  - `self.dpnv_grouping_dict_c`
  - `self.p`

### `Winding_Derivation.get_complex_number_winding_factor_of_coil_i` (Line 444)
- **Arguments**: self, i, coil_pitch_y, Q, v, p
- **Reads (State Dependencies)**:
  - `self.verbose`
- **Writes**: None state variables detected

### `Winding_Derivation.get_complex_number_kw_per_phase` (Line 458)
- **Arguments**: self, v, p, positive_connected_coils, negative_connected_coils
- **Reads (State Dependencies)**:
  - `self.coil_pitch_y`
  - `self.get_complex_number_winding_factor_of_coil_i`
  - `self.verbose`
  - `self.Q`
- **Writes**: None state variables detected

### `Winding_Derivation.get_complex_number_kw` (Line 486)
- **Arguments**: self, p_or_ps, v, bool_study_suspension_subharmonics
- **Reads (State Dependencies)**:
  - `self.dpnv_grouping_dict_a`
  - `self.connection_star_raw_dict`
  - `self.dpnv_grouping_dict_c['GBD']`
  - `self.verbose`
  - `self.ps`
  - `self.dpnv_grouping_dict_b`
  - `self.dpnv_grouping_dict_c['GAC']`
  - `self.dpnv_grouping_dict_a['GAC']`
  - `self.dpnv_grouping_dict_a['GBD']`
  - `self.dpnv_grouping_dict_b['GBD']`
  - `self.get_complex_number_kw_per_phase`
  - `self.dpnv_grouping_dict_c`
  - `self.dpnv_grouping_dict_b['GAC']`
- **Writes**: None state variables detected

### `Winding_Derivation.format_print_out_string` (Line 531)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.dict_suspension_kw_els['C_angle']`
  - `self.connection_star_raw_dict`
  - `self.ps`
  - `self.bool_double_layer_winding`
  - `self.Q`
  - `self.dict_suspension_kw_cjh`
  - `self.dict_suspension_kw_els['B_angle']`
  - `self.verbose`
  - `self.dict_torque_kw_cjh`
  - `self.grouping_AC`
  - `self.get_complex_number_kw`
  - `self.layer_X_signs`
  - `self.dict_suspension_kw_els['A_angle']`
  - `self.m`
  - `self.dict_torque_kw_els`
  - `self.dpnv_grouping_dict_a`
  - `self.print_out_string`
  - `self.layer_X_phases`
  - `self.dict_suspension_kw_els`
  - `self.dpnv_grouping_dict_b`
  - `self.coil_pitch_y`
  - `self.dpnv_grouping_dict_c`
  - `self.p`
- **Writes (State Changes)**:
  - `self.print_out_string`
  - `self.grouping_AC`
  - `self.layer_X_signs`
  - `self.coil_pitch_y`
  - `self.layer_X_phases`

### `main_derivation` (Line 655)
- **Arguments**: m, Qs, p, ps, coil_pitch_y, verbose
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `WireSlot_v1.py`

### `calc_stator_geometry` (Line 11)
- **Arguments**: OD, ID, tooth_width, tooth_depth, yoke, liner, slots
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `calc_winding_capacity` (Line 32)
- **Arguments**: net_area, gross_area, awg, orthocyclic_factor, fill_heuristic
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `calc_thermal_load` (Line 64)
- **Arguments**: ID, z_slot, a_bare, J, slots
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `calc_bemf_constants` (Line 84)
- **Arguments**: La, yoke, z_slot, B_sat, poles, slots
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `calc_motor_losses` (Line 120)
- **Arguments**: La, tooth_width, turns_per_phase, current, rpm, awg, poles
- **Reads**: None state variables detected
- **Writes**: None state variables detected

## File: `WireSlot_v2.py`

### `MotorThermalAnalyzer.__init__` (Line 9)
- **Arguments**: self, stator_od, rotor_od, air_gap, tooth_depth, tooth_width, slots, liner
- **Reads (State Dependencies)**:
  - `self.rotor_od`
  - `self.slots`
  - `self.stator_od`
  - `self.tooth_width`
  - `self.liner`
  - `self.air_gap`
  - `self.stator_id`
  - `self.yoke_thickness`
  - `self.tooth_depth`
- **Writes (State Changes)**:
  - `self.rotor_od`
  - `self.slots`
  - `self.stator_od`
  - `self.tooth_width`
  - `self.liner`
  - `self.air_gap`
  - `self.stator_id`
  - `self.yoke_thickness`
  - `self.tooth_depth`

### `MotorThermalAnalyzer.calculate_slot_area` (Line 33)
- **Arguments**: self
- **Reads (State Dependencies)**:
  - `self.stator_id`
  - `self.tooth_depth`
  - `self.slots`
  - `self.tooth_width`
- **Writes**: None state variables detected

### `MotorThermalAnalyzer.get_wire_properties` (Line 44)
- **Arguments**: self, awg_size
- **Reads**: None state variables detected
- **Writes**: None state variables detected

### `MotorThermalAnalyzer.estimate_max_wires_in_slot` (Line 61)
- **Arguments**: self, d_od
- **Reads (State Dependencies)**:
  - `self.slots`
  - `self.tooth_width`
  - `self.liner`
  - `self.stator_id`
  - `self.tooth_depth`
- **Writes**: None state variables detected

### `MotorThermalAnalyzer.analyze_thermal_performance` (Line 104)
- **Arguments**: self, awg_size, target_j
- **Reads (State Dependencies)**:
  - `self.calculate_slot_area`
  - `self.slots`
  - `self.stator_id`
  - `self.estimate_max_wires_in_slot`
  - `self.yoke_thickness`
  - `self.get_wire_properties`
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

