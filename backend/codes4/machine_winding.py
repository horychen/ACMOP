import math
from types import SimpleNamespace
import modern_machine_designer_utility

class MachineWinding:
    def __init__(self, winding_dict: dict = None):
        self.winding_dict = winding_dict if winding_dict is not None else {}

    def get_InitialRotationAngle(self, machineGeometry):
        p = self.winding_dict.get('pole_count', 2) // 2
        deg_pole_span = 180 / p
        self.winding_dict['initial_rotation_angle'] = (deg_pole_span - machineGeometry.get_deg_alpha_rm()) * 0.5 + self.winding_dict['wily'].deg_winding_U_phase_phase_axis_angle
        return self.winding_dict['initial_rotation_angle']

    def update_excitations(self, machineGeometry, materials):
        p = self.winding_dict.get('pole_count', 2) // 2
        Qs = self.winding_dict.get('slot_count', 12)
        m = self.winding_dict.get('phase_count', 3)
        
        RatedSpeed = self.winding_dict.get('rated_speed', 1000)
        self.winding_dict['excitation_frequency_simulated'] = RatedSpeed / 60 * p
        ExcitationFreqSimulated = self.winding_dict['excitation_frequency_simulated']
        is_wye_connection = self.winding_dict.get('connection_type', 'wye').lower() == "wye"
        V_stator_phase_voltage_amp = math.sqrt(2) * self.winding_dict.get('dc_bus_voltage', 12) / (math.sqrt(3) if is_wye_connection else 1.0)
        leakage_inductance_factor = 0.95
        V_desired_emf_Em = leakage_inductance_factor * V_stator_phase_voltage_amp
        alpha_i = 2.0/math.pi

        T_air_gap_flux_density_Bg_guessed = materials.material_dict['magnet_bg']

        mm_stack_length_specified = self.winding_dict.get('stack_length_specified', 16.0)
        mm_d_magnetic_air_gap = machineGeometry.gp['d_air_gap']
        mm_stack_length_effective = mm_stack_length_specified + 2 * mm_d_magnetic_air_gap
        r_si = machineGeometry.gp['r_rotor_outer'] + machineGeometry.gp['d_air_gap']
        mm_pole_pitch_tau_p = math.pi * r_si / p
        Wb_air_gap_flux_Phi_m = alpha_i * T_air_gap_flux_density_Bg_guessed * mm_pole_pitch_tau_p*1e-3 * mm_stack_length_effective*1e-3
        
        no_series_coil_turns_N = V_desired_emf_Em / (2*math.pi* ExcitationFreqSimulated * self.winding_dict['wily'].kw1 * Wb_air_gap_flux_Phi_m) if Wb_air_gap_flux_Phi_m > 0 else 1
        no_series_coil_turns_N = round(no_series_coil_turns_N)
        SPP = Qs / (2*p*m)
        bool_weHavePlentyVoltage = True
        if bool_weHavePlentyVoltage:
            no_series_coil_turns_N = min([p*SPP*i for i in range(1000,0,-1)], key=lambda x:abs(x - no_series_coil_turns_N))
        else:
            no_series_coil_turns_N = min([p*SPP*i for i in range(1000)], key=lambda x:abs(x - no_series_coil_turns_N))
        
        if no_series_coil_turns_N > 990:
            no_series_coil_turns_N = 990

        self.winding_dict['series_turns'] = no_series_coil_turns_N
        self.winding_dict['drive_winding_conductors_per_slot'] = 2 * m * no_series_coil_turns_N / Qs * self.winding_dict['wily'].number_of_parallel_branch
        self.winding_dict['conductors_per_slot'] = self.winding_dict['drive_winding_conductors_per_slot']
        
        torque_current_ratio = self.winding_dict.get('torque_current_ratio', 1.0)
        suspension_current_ratio = self.winding_dict.get('suspension_current_ratio', 0.0)

        self.winding_dict['bearing_winding_conductors_per_slot'] = self.winding_dict['drive_winding_conductors_per_slot'] if self.winding_dict['wily'].bool_DPNVorSEPA == True else self.winding_dict['drive_winding_conductors_per_slot'] / torque_current_ratio * suspension_current_ratio if torque_current_ratio > 0 else 0

        mm_r_sy = machineGeometry.gp['r_stator_outer'] - machineGeometry.get_d_stator_yoke()
        mm_r_si = machineGeometry.gp['r_rotor_outer'] + machineGeometry.gp['d_air_gap']
        mm_r_ss = mm_r_si + machineGeometry.get_d_tooth_shoe()
        self.winding_dict['slot_area'] = (math.pi*(mm_r_sy**2 - mm_r_ss**2) / Qs - machineGeometry.gp['w_tooth'] * (machineGeometry.gp['d_tooth'] - machineGeometry.get_d_tooth_shoe()))
        self.winding_dict['slot_current_amplitude']   = self.winding_dict['slot_area'] * self.winding_dict.get('rated_current_density', 5.0) * self.winding_dict.get('fill_factor', 0.5) * math.sqrt(2)
        self.winding_dict['conductor_current_amplitude'] = self.winding_dict['slot_current_amplitude'] / self.winding_dict['drive_winding_conductors_per_slot'] if self.winding_dict['drive_winding_conductors_per_slot'] > 0 else 0
        self.winding_dict['phase_current_amplitude']     = self.winding_dict['conductor_current_amplitude'] * self.winding_dict['wily'].number_of_parallel_branch
        self.winding_dict['drive_winding_current'] = torque_current_ratio     * self.winding_dict['phase_current_amplitude']
        self.winding_dict['bearing_winding_current'] = suspension_current_ratio * self.winding_dict['phase_current_amplitude']
        self.winding_dict['torque_current_utilization_ratio'] = (self.winding_dict['drive_winding_current'] + self.winding_dict['bearing_winding_current']) / self.winding_dict['phase_current_amplitude'] if self.winding_dict['phase_current_amplitude'] > 0 else 0

        self.get_InitialRotationAngle(machineGeometry)

        Rout = machineGeometry.gp['r_rotor_outer']
        Rin  = Rout - machineGeometry.gp['d_magnet']
        deg_alpha_rp = 360 / (2*p)
        deg_alpha_rm = 360 / (2*p)
        self.winding_dict['magnet_area'] = deg_alpha_rm/deg_alpha_rp * math.pi*(Rout**2 - Rin**2)

    def sync(self, machineGeometry, materials, l_stack):
        self.winding_dict['rotor_core_material'] = materials.material_dict['rotor_steel']
        self.winding_dict['stator_core_material'] = materials.material_dict['stator_steel']
        self.winding_dict['lamination_factor'] = materials.material_dict['steel_stack_factor']
        self.winding_dict['magnet_name'] = materials.material_dict['magnet_grade']
        self.winding_dict['magnet_temperature'] = 20
        p = self.winding_dict.get('pole_count', 2) // 2
        deg_alpha_rm = 360 / (2*p)
        self.winding_dict['magnet_start_angle'] = 0.0 * deg_alpha_rm
        self.winding_dict['stack_length_specified'] = l_stack
        
        self.winding_dict['wily'] = modern_machine_designer_utility.Winding(
            phase_number_m=self.winding_dict.get('phase_count', 3),
            stator_slot_number_Qs=self.winding_dict.get('slot_count', 12),
            pole_pair_number_p=p,
            suspension_pole_pair_number_ps=self.winding_dict.get('ps', 1),
            coil_pitch_y=self.winding_dict.get('coil_pitch', 1),
            bool_DPNVorSEPA=self.winding_dict.get('DPNV_or_SEPA', True),
            number_of_parallel_branch=self.winding_dict.get('wily').number_of_parallel_branch if self.winding_dict.get('wily') else 1
        )
        
        self.update_excitations(machineGeometry, materials)

        r_si = machineGeometry.gp['r_rotor_outer'] + machineGeometry.gp['d_air_gap']
        hm = machineGeometry.gp['d_magnet']
        gap_dist = machineGeometry.gp['d_air_gap']
        
        br = materials.material_dict['magnet_br']
        b_gap = br * hm / (hm + 1.05 * gap_dist) if (hm + 1.05 * gap_dist) > 0 else 0
        print(f'b_gap = {b_gap} T')

        wire_diameter = self.winding_dict.get('wire_diameter', 1.0)
        wire_area = math.pi * (wire_diameter/2)**2
        conductor_current = self.winding_dict.get('rated_current_density', 5.0) * wire_area
        self.winding_dict['slot_current_at'] = conductor_current * self.winding_dict.get('conductors_per_slot', 0)
        
        end_winding = math.pi * (2*r_si) / self.winding_dict.get('slot_count', 12) * self.winding_dict.get('coil_pitch', 1)
        n_series = (self.winding_dict.get('conductors_per_slot', 0) * self.winding_dict.get('slot_count', 12)) / (2 * self.winding_dict.get('phase_count', 3) * self.winding_dict['wily'].number_of_parallel_branch)
        turn_length = 2 * (l_stack + end_winding) * 1e-3
        rho_copper = 1.72e-8
        self.winding_dict['phase_resistance'] = rho_copper * (turn_length * n_series) / (wire_area * 1e-6) if wire_area > 0 else 0
        
        pole_count = self.winding_dict.get('pole_count', 2)
        self.winding_dict['excitation_frequency'] = (self.winding_dict.get('rated_speed', 1000) * (pole_count / 2)) / 60.0
        
        area_pole = (2 * math.pi * r_si * l_stack) / pole_count
        estimated_flux = b_gap * area_pole * 1e-6 
        angular_speed = self.winding_dict.get('rated_speed', 1000) * (2 * math.pi / 60.0)
        winding_factor = 0.93
        self.winding_dict['estimated_back_emf'] = estimated_flux * n_series * winding_factor * angular_speed

    def validate_fill_factor(self):
        import math
        zQ = self.winding_dict.get('drive_winding_conductors_per_slot', 0)
        D_wire_insul = self.winding_dict.get('wire_diameter_with_insulation', 1.0)
        D_wire_bare = self.winding_dict.get('wire_diameter', 1.0)
        
        area_one_insulated_strand = math.pi * (D_wire_insul / 2)**2
        total_insulated_area = zQ * area_one_insulated_strand
        
        A_slot = self.winding_dict.get('slot_area', 1e-6)
        true_fill_factor = total_insulated_area / A_slot
        
        target_fill_factor = self.winding_dict.get('fill_factor', 0.5)
        
        rho_copper_hot = 2.0e-8 
        J_A_m2 = self.winding_dict.get('rated_current_density', 5.0) * 1e6
        loss_density_W_m3 = rho_copper_hot * (J_A_m2**2) * true_fill_factor
        loss_density_W_cm3 = loss_density_W_m3 * 1e-6
        temp_rise_estimate = loss_density_W_cm3 * 60.0
        
        is_valid = true_fill_factor <= target_fill_factor
        status_word = "PASSED" if is_valid else "FAILED"
        
        report = []
        report.append(f"====== WIRE GAUGE & FILL FACTOR VALIDATION: {status_word} ======")
        report.append(f"  Slot Area (available): {A_slot:.2f} mm^2")
        report.append(f"  Req. Conductors/Slot (zQ): {zQ:.1f} strands")
        report.append(f"  Wire Dia (insulated): {D_wire_insul:.3f} mm (Area = {area_one_insulated_strand:.3f} mm^2/strand)")
        report.append(f"  Total Insulated Wire Area: {total_insulated_area:.2f} mm^2")
        report.append(f"")
        report.append(f"  Target Fill Factor (User Limit): {target_fill_factor:.3f} ({(target_fill_factor*100):.1f}%)")
        report.append(f"  True Physical Fill Factor:     {true_fill_factor:.3f} ({(true_fill_factor*100):.1f}%)")
        
        report.append(f"")
        report.append(f"  --- Thermal Performance Estimation ---")
        report.append(f"  RMS Current Density (J_s): {self.winding_dict.get('rated_current_density', 5.0):.1f} A/mm^2")
        report.append(f"  Est. Volumetric Heat Load: {loss_density_W_cm3:.2f} W/cm^3 of slot")
        report.append(f"  Est. Steady-State Temp Rise: ~{temp_rise_estimate:.1f} °C")
        if temp_rise_estimate > 120.0:
            report.append("  [WARNING] High Thermal Risk! Consider reducing current density (J_s) or improving cooling.")
        
        if not is_valid:
            report.append(f"")
            report.append(f"  [!] VALIDATION FAILED: Wires cannot physically fit into the slot!")
            report.append(f"  --- ACTIONABLE OPTIONS TO RESOLVE MISTMATCH ---")
            
            required_wire_area = (A_slot * target_fill_factor) / zQ if zQ > 0 else 0
            req_D_insul = math.sqrt(required_wire_area / math.pi) * 2
            req_D_bare = req_D_insul * (D_wire_bare / D_wire_insul) if D_wire_insul > 0 else 0
            req_bare_area = math.pi * (req_D_bare / 2)**2
            orig_bare_area = math.pi * (D_wire_bare / 2)**2
            orig_cond_current = orig_bare_area * self.winding_dict.get('rated_current_density', 5.0)
            new_Js = orig_cond_current / req_bare_area if req_bare_area > 0 else 0
            
            report.append(f"  Option 1: Thinner Wire & Higher Current Density")
            report.append(f"            - Use thinner wire with insulated diameter <= {req_D_insul:.3f} mm")
            report.append(f"            - This will increase required Current Density (J_s) to ~{new_Js:.1f} A/mm^2")
            
            required_A_slot = total_insulated_area / target_fill_factor
            report.append(f"  Option 2: Increase Mechanical Slot Area")
            report.append(f"            - Increase Slot Area from {A_slot:.2f} mm^2 to at least {required_A_slot:.2f} mm^2")
            report.append(f"            - Adjust geometry: increase stator outer radius or decrease tooth width.")
            
            required_zQ = (A_slot * target_fill_factor) / area_one_insulated_strand
            report.append(f"  Option 3: Relax Electrical/Operating Constraints")
            report.append(f"            - Reduce required conductors per slot from {zQ:.1f} to {required_zQ:.1f} turns/slot")
            report.append(f"            - This reduces your achievable back-EMF proportionally (e.g., lower base speed or lower DC bus voltage).")

        report.append("================================================================")
        return is_valid, "\n".join(report)

if __name__ == "__main__":
    w = MachineWinding()
    print("MachineWinding Debug:")
