from dataclasses import dataclass, field
from typing import Any
import math
from types import SimpleNamespace
import modern_machine_designer_utility

@dataclass
class MachineWinding:
    winding_dict: dict = field(default_factory=dict)
    
    # --- Basic Specifications ---
    @property
    def phase_count(self) -> int: return self.winding_dict['phase_count']
    @property
    def slot_count(self) -> int: return self.winding_dict['slot_count']
    @property
    def pole_count(self) -> int: return self.winding_dict['pole_count']
    @property
    def coil_pitch(self) -> int: return self.winding_dict['coil_pitch']
    
    # Aliases/Parameters for legacy compatibility
    @property
    def Qs(self): return self.slot_count
    @property
    def p(self): return self.pole_count // 2
    
    # TODO: use zQ to compute fill factor
    # @property
    # def conductors_per_slot(self) -> int: return self.winding_dict['conductors_per_slot']
    @property
    def wire_diameter_with_insulation(self) -> float: return self.winding_dict['wire_diameter_with_insulation']
    @property
    def wire_diameter(self) -> float: return self.winding_dict['wire_diameter']
    @property
    def rated_current_density(self) -> float: return self.winding_dict['rated_current_density']
    @property
    def parallel_branch_count(self) -> int: return self.wily.number_of_parallel_branch
    
    # # --- Expanded Excitation & Thermal ---
    @property
    def connection_type(self) -> str: return self.winding_dict['connection_type']
    @property
    def is_wye_connection(self) -> bool: return self.connection_type == "Wye"
    @property
    def rated_speed(self) -> float: return self.winding_dict['rated_speed']
    @property
    def rated_power(self) -> float: return self.winding_dict['rated_power']
    @property
    def dc_bus_voltage(self) -> float: return self.winding_dict['dc_bus_voltage']
    @property
    def fill_factor(self) -> float: return self.winding_dict['fill_factor']
    # @property
    # def excitation_frequency(self) -> float: return self.winding_dict['excitation_frequency']
    
    # @property
    # def is_wye_connection(self) -> bool: return self.winding_dict['is_wye_connection']
    
    # # Bearingless specific
    @property
    def torque_current_ratio(self) -> float: return self.winding_dict['torque_current_ratio']
    @property
    def suspension_current_ratio(self) -> float: return self.winding_dict['suspension_current_ratio']
    @property
    def drive_winding_resistance(self) -> float: return self.winding_dict['drive_winding_resistance']
    @property
    def bearing_winding_resistance(self) -> float: return self.winding_dict['bearing_winding_resistance']
    
    # Performance & Design Targets
    # series_turns: int = 0 # (no_series_coil_turns_N)
    # drive_winding_conductors_per_slot: float = 0.0 # (DriveW_zQ)
    # bearing_winding_conductors_per_slot: float = 0.0 # (BeariW_zQ)
    # slot_area: float = 0.0 # (mm2_slot_area)
    # slot_current_amplitude: float = 0.0 # (CurrentAmp_in_the_slot)
    # conductor_current_amplitude: float = 0.0 # (CurrentAmp_per_conductor)
    # phase_current_amplitude: float = 0.0 # (CurrentAmp_per_phase)
    # drive_winding_current: float = 0.0 # (DriveW_CurrentAmp)
    # bearing_winding_current: float = 0.0 # (BeariW_CurrentAmp)
    # torque_current_utilization_ratio: float = 0.0 # (slot_current_utilizing_ratio_for_torque)
    # magnet_area: float = 0.0 # (mm2_magnet_area)
    # initial_rotation_angle: float = 0.0 # (InitialRotationAngle)
    
    # Simulation & Configuration
    bool_DPNVorSEPA: bool = True
    bool_3PhaseCurrentSource: bool = False 
    # Derived results
    phase_resistance: float = 0.0
    estimated_back_emf: float = 0.0
    slot_current_at: float = 0.0
    # wily: Any = field(default=None, init=False)
    ps: int = field(init=False) # suspension_pole_pair_count

    def __post_init__(self):
        from types import SimpleNamespace
        p = self.pole_count // 2
        # Try p+1 first, then p-1, and avoid multiples of phase_count (typically 3)
        for candidate in [p + 1, p - 1]:
            if candidate > 0 and candidate % self.phase_count != 0:
                self.ps = candidate
                break

        # get winding layout
        # val_wily = modern_machine_designer_utility.Winding(phase_number_m=3, stator_slot_number_Qs=24, pole_pair_number_p=2, suspension_pole_pair_number_ps=1, coil_pitch_y=6, bool_DPNVorSEPA=False, number_of_parallel_branch=1)
        # val_wily = modern_machine_designer_utility.Winding(phase_number_m=3, stator_slot_number_Qs=12, pole_pair_number_p=5, suspension_pole_pair_number_ps=1, coil_pitch_y=1, bool_DPNVorSEPA=True, number_of_parallel_branch=2)
        # val_wily = modern_machine_designer_utility.Winding(phase_number_m=3, stator_slot_number_Qs=12, pole_pair_number_p=2, suspension_pole_pair_number_ps=1, coil_pitch_y=3, bool_DPNVorSEPA=None, number_of_parallel_branch=1)
        self.wily = modern_machine_designer_utility.Winding(
            phase_number_m=self.phase_count, 
            stator_slot_number_Qs=self.slot_count, 
            pole_pair_number_p=p, 
            suspension_pole_pair_number_ps=self.ps, 
            coil_pitch_y=self.coil_pitch, 
            bool_DPNVorSEPA=self.bool_DPNVorSEPA
        )

    def get_InitialRotationAngle(self, machineGeometry):
        deg_pole_span = 180 / self.p
        self.initial_rotation_angle = (deg_pole_span - machineGeometry.deg_alpha_rm) * 0.5 + self.wily.deg_winding_U_phase_phase_axis_angle
        return self.initial_rotation_angle

    def update_excitations(self, machineGeometry, materials):
        # This logic is based on the user-provided snippet
        p = self.p
        Qs = self.slot_count
        m = self.phase_count
        
        RatedSpeed = self.rated_speed
        self.excitation_frequency_simulated = RatedSpeed / 60 * p
        ExcitationFreqSimulated = self.excitation_frequency_simulated
        V_stator_phase_voltage_amp = math.sqrt(2) * self.dc_bus_voltage / (math.sqrt(3) if self.is_wye_connection else 1.0) # SVPWM
        leakage_inductance_factor = 0.95
        V_desired_emf_Em = leakage_inductance_factor * V_stator_phase_voltage_amp
        alpha_i = 2.0/math.pi

        T_air_gap_flux_density_Bg_guessed = materials.magnet_bg  # T

        mm_stack_length_specified = self.stack_length_specified
        mm_d_magnetic_air_gap = machineGeometry.d_air_gap # + machineGeometry.d_sleeve
        mm_stack_length_effective = mm_stack_length_specified + 2 * mm_d_magnetic_air_gap
        r_si = machineGeometry.r_rotor_outer + machineGeometry.d_air_gap
        mm_pole_pitch_tau_p = math.pi * r_si / p
        Wb_air_gap_flux_Phi_m = alpha_i * T_air_gap_flux_density_Bg_guessed * mm_pole_pitch_tau_p*1e-3 * mm_stack_length_effective*1e-3 # Wb
        
        # Calculate turns
        no_series_coil_turns_N = V_desired_emf_Em / (2*math.pi* ExcitationFreqSimulated * self.wily.kw1 * Wb_air_gap_flux_Phi_m) if Wb_air_gap_flux_Phi_m > 0 else 1
        no_series_coil_turns_N = round(no_series_coil_turns_N)
        SPP = Qs / (2*p*m)
        bool_weHavePlentyVoltage = True
        if bool_weHavePlentyVoltage:
            no_series_coil_turns_N = min([p*SPP*i for i in range(1000,0,-1)], key=lambda x:abs(x - no_series_coil_turns_N))
        else:
            no_series_coil_turns_N = min([p*SPP*i for i in range(1000)], key=lambda x:abs(x - no_series_coil_turns_N))
        
        if no_series_coil_turns_N > 990:
            no_series_coil_turns_N = 990 # Cap instead of raise for stability during optimization

        self.series_turns = no_series_coil_turns_N
        self.drive_winding_conductors_per_slot = self.conductors_per_slot = 2* m * no_series_coil_turns_N / Qs * self.wily.number_of_parallel_branch
        self.bearing_winding_conductors_per_slot = self.drive_winding_conductors_per_slot if self.wily.bool_DPNVorSEPA == True else self.drive_winding_conductors_per_slot / self.torque_current_ratio * self.suspension_current_ratio

        ''' Excitations Consiering Thermal Capability Limit (Simple) '''
        mm_r_sy = machineGeometry.r_stator_outer - machineGeometry.d_stator_yoke
        mm_r_si = machineGeometry.r_rotor_outer + machineGeometry.d_air_gap
        mm_r_ss = mm_r_si + machineGeometry.d_tooth_shoe
        self.slot_area = (math.pi*(mm_r_sy**2 - mm_r_ss**2) / Qs - machineGeometry.w_tooth * (machineGeometry.d_tooth - machineGeometry.d_tooth_shoe))
        self.slot_current_amplitude   = self.slot_area * self.rated_current_density * self.fill_factor * math.sqrt(2)
        self.conductor_current_amplitude = self.slot_current_amplitude / self.drive_winding_conductors_per_slot if self.drive_winding_conductors_per_slot > 0 else 0
        self.phase_current_amplitude     = self.conductor_current_amplitude * self.wily.number_of_parallel_branch
        self.drive_winding_current = self.torque_current_ratio     * self.phase_current_amplitude
        self.bearing_winding_current = self.suspension_current_ratio * self.phase_current_amplitude
        self.torque_current_utilization_ratio = (self.drive_winding_current + self.bearing_winding_current) / self.phase_current_amplitude if self.phase_current_amplitude > 0 else 0

        self.initial_rotation_angle = self.get_InitialRotationAngle(machineGeometry)

        Rout = machineGeometry.r_rotor_outer
        Rin  = Rout - machineGeometry.d_magnet
        deg_alpha_rp = 360 / (2*self.p)
        deg_alpha_rm = 360 / (2*self.p)
        self.magnet_area = deg_alpha_rm/deg_alpha_rp * math.pi*(Rout**2 - Rin**2)

    def sync(self, machineGeometry, materials, l_stack):
        """Update derived performance parameters using geometry and materials."""
        # 0.1 Update Materials individually
        self.rotor_core_material = materials.rotor_steel
        self.stator_core_material = materials.stator_steel
        self.lamination_factor = materials.steel_stack_factor
        self.magnet_name = materials.magnet_grade
        self.magnet_temperature = 20 # Default
        deg_alpha_rm = 360 / (2*self.p)
        self.magnet_start_angle = 0.0 * deg_alpha_rm # Default 永磁体励磁默认是N极朝右，S极朝左；如果是1*deg_alpha_rm，则是N极朝左，S极朝右
        self.stack_length_specified = l_stack
        
        # Update wily
        # Update wily
        self.wily = modern_machine_designer_utility.Winding(
            phase_number_m=self.phase_count,
            stator_slot_number_Qs=self.slot_count,
            pole_pair_number_p=self.p,
            suspension_pole_pair_number_ps=self.ps,
            coil_pitch_y=self.coil_pitch,
            bool_DPNVorSEPA=self.bool_DPNVorSEPA,
            number_of_parallel_branch=self.wily.number_of_parallel_branch
        )
        
        # 1. Update Excitations (turns, current, etc.)
        self.update_excitations(machineGeometry, materials)

        m = materials
        
        # 1. Geometry-derived values
        r_si = machineGeometry.r_rotor_outer + machineGeometry.d_air_gap
        r_ro = machineGeometry.r_rotor_outer
        hm = machineGeometry.d_magnet
        gap_dist = machineGeometry.d_air_gap
        
        # 2. Bg Estimation
        br = m.magnet_br
        b_gap = br * hm / (hm + 1.05 * gap_dist) if (hm + 1.05 * gap_dist) > 0 else 0
        print(f'b_gap = {b_gap} T')
        # quit()

        # 3. Electrical Metrics
        wire_area = math.pi * (self.wire_diameter/2)**2
        conductor_current = self.rated_current_density * wire_area # [A]
        self.slot_current_at = conductor_current * self.conductors_per_slot
        
        # 4. Resistance Estimation
        end_winding = math.pi * (2*r_si) / self.slot_count * self.coil_pitch
        n_series = (self.conductors_per_slot * self.slot_count) / (2 * self.phase_count * self.parallel_branch_count)
        turn_length = 2 * (l_stack + end_winding) * 1e-3 # [m]
        rho_copper = 1.72e-8 # [Ohm-m]
        self.phase_resistance = rho_copper * (turn_length * n_series) / (wire_area * 1e-6)
        
        # 5. Back-EMF
        pole_count = self.pole_count
        self.excitation_frequency = (self.rated_speed * (pole_count / 2)) / 60.0
        
        area_pole = (2 * math.pi * r_si * l_stack) / pole_count
        estimated_flux = b_gap * area_pole * 1e-6 
        angular_speed = self.rated_speed * (2 * math.pi / 60.0)
        winding_factor = 0.93
        self.estimated_back_emf = estimated_flux * n_series * winding_factor * angular_speed

    def validate_fill_factor(self):
        """
        Validate if the chosen wire gauge can physically fit into the mechanical slot area
        based on the desired conductors per slot (zQ) and assigned fill factor limit.
        Generates a detailed diagnostic report with actionable options and thermal estimation.
        """
        import math
        
        # 1. Physical Calculations
        zQ = getattr(self, 'drive_winding_conductors_per_slot', 0)
        D_wire_insul = self.wire_diameter_with_insulation
        D_wire_bare = self.wire_diameter
        
        # Area of ONE insulated strand
        area_one_insulated_strand = math.pi * (D_wire_insul / 2)**2
        # Total area required for all strands
        total_insulated_area = zQ * area_one_insulated_strand
        
        # True Fill Factor computation
        A_slot = getattr(self, 'slot_area', 1e-6) # avoid div by zero
        true_fill_factor = total_insulated_area / A_slot
        
        # Target Fill Factor assigned by user
        target_fill_factor = self.fill_factor
        
        # Thermal Estimation (Very simplified steady-state Copper Loss density)
        # Assuming typical slot perimeter and air cooling. 
        # Heat load (W/m^3 of slot) \propto J_s^2 * k_fill
        # rho_Cu at 20C roughly 1.72e-8 Ohm-m, typically 2.2e-8 at 100C. We use ~2.0e-8
        rho_copper_hot = 2.0e-8 
        # Copper loss density in slot [W/m^3] = rho * (J * 1e6)^2 * true_fill_factor
        # where J is in A/mm^2 = A/(1e-6 m^2) = 1e6 A/m^2
        J_A_m2 = self.rated_current_density * 1e6
        loss_density_W_m3 = rho_copper_hot * (J_A_m2**2) * true_fill_factor
        # Rough empirical conversion back to surface Heat Flux or Temperature Rise
        # Very rough rule of thumb for small enclosed motors: 1 W/cm^3 \approx 50-80 degC rise.
        loss_density_W_cm3 = loss_density_W_m3 * 1e-6
        temp_rise_estimate = loss_density_W_cm3 * 60.0 # degC roughly
        
        # 2. Report Generation
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
        
        # Thermal Report
        report.append(f"")
        report.append(f"  --- Thermal Performance Estimation ---")
        report.append(f"  RMS Current Density (J_s): {self.rated_current_density:.1f} A/mm^2")
        report.append(f"  Est. Volumetric Heat Load: {loss_density_W_cm3:.2f} W/cm^3 of slot")
        report.append(f"  Est. Steady-State Temp Rise: ~{temp_rise_estimate:.1f} °C")
        if temp_rise_estimate > 120.0:
            report.append("  [WARNING] High Thermal Risk! Consider reducing current density (J_s) or improving cooling.")
        
        if not is_valid:
            report.append(f"")
            report.append(f"  [!] VALIDATION FAILED: Wires cannot physically fit into the slot!")
            report.append(f"  --- ACTIONABLE OPTIONS TO RESOLVE MISTMATCH ---")
            
            # Option 1: Find Required J_s to drop wire size
            # If we keep A_slot and zQ, we need the insulated area to drop to match target_fill_factor
            required_wire_area = (A_slot * target_fill_factor) / zQ if zQ > 0 else 0
            req_D_insul = math.sqrt(required_wire_area / math.pi) * 2
            # Assuming bare wire ratio follows same proportionality
            req_D_bare = req_D_insul * (D_wire_bare / D_wire_insul) if D_wire_insul > 0 else 0
            req_bare_area = math.pi * (req_D_bare / 2)**2
            # New J_s needed = original conductor current / new bare wire area
            # original conductor current = original bare area * original J_s
            orig_bare_area = math.pi * (D_wire_bare / 2)**2
            orig_cond_current = orig_bare_area * self.rated_current_density
            new_Js = orig_cond_current / req_bare_area if req_bare_area > 0 else 0
            
            report.append(f"  Option 1: Thinner Wire & Higher Current Density")
            report.append(f"            - Use thinner wire with insulated diameter <= {req_D_insul:.3f} mm")
            report.append(f"            - This will increase required Current Density (J_s) to ~{new_Js:.1f} A/mm^2")
            
            # Option 2: Increase Slot Area
            required_A_slot = total_insulated_area / target_fill_factor
            report.append(f"  Option 2: Increase Mechanical Slot Area")
            report.append(f"            - Increase Slot Area from {A_slot:.2f} mm^2 to at least {required_A_slot:.2f} mm^2")
            report.append(f"            - Adjust geometry: increase stator outer radius or decrease tooth width.")
            
            # Option 3: Lower DC Bus Constraint / Speed to drop required turns (zQ)
            required_zQ = (A_slot * target_fill_factor) / area_one_insulated_strand
            report.append(f"  Option 3: Relax Electrical/Operating Constraints")
            report.append(f"            - Reduce required conductors per slot from {zQ:.1f} to {required_zQ:.1f} turns/slot")
            report.append(f"            - This reduces your achievable back-EMF proportionally (e.g., lower base speed or lower DC bus voltage).")

        report.append("================================================================")
        return is_valid, "\n".join(report)

if __name__ == "__main__":
    w = MachineWinding()
    print("MachineWinding Debug:")
    print(f"Num Slots: {w.slot_count}, Num Poles: {w.pole_count}")
    # Partial sync test would require mock geometry/materials
