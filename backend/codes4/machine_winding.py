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
    def phase_count(self) -> int: return self.winding_dict.get('phase_count', 3)
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
    @property
    def conductors_per_slot(self) -> int: return self.winding_dict.get('conductors_per_slot', 42)
    @property
    def wire_diameter_with_insulation(self) -> float: return self.winding_dict.get('wire_diameter_with_insulation', 0.23)
    @property
    def wire_diameter(self) -> float: return self.winding_dict.get('wire_diameter', 0.21)
    @property
    def rated_current_density(self) -> float: return self.winding_dict.get('rated_current_density', 5.0)
    @property
    def parallel_branch_count(self) -> int: return self.wily.number_of_parallel_branch
    
    # # --- Expanded Excitation & Thermal ---
    @property
    def connection_type(self) -> str: return self.winding_dict.get('connection_type', "Wye")
    @property
    def is_wye_connection(self) -> bool: return self.connection_type == "Wye"
    @property
    def rated_speed(self) -> float: return self.winding_dict.get('rated_speed', 10000.0)
    # @property
    # def rated_power(self) -> float: return self.winding_dict.get('rated_power', 50.0)
    @property
    def dc_bus_voltage(self) -> float: return self.winding_dict.get('dc_bus_voltage', 24.0)
    @property
    def fill_factor(self) -> float: return self.winding_dict.get('fill_factor', 0.4)
    # @property
    # def excitation_frequency(self) -> float: return self.winding_dict.get('excitation_frequency', 250.0)
    
    # @property
    # def is_wye_connection(self) -> bool: return self.winding_dict.get('is_wye_connection', True)
    
    # # Bearingless specific
    @property
    def torque_current_ratio(self) -> float: return self.winding_dict.get('torque_current_ratio', 1.0)
    @property
    def suspension_current_ratio(self) -> float: return self.winding_dict.get('suspension_current_ratio', 0.0)
    @property
    def drive_winding_resistance(self) -> float: return self.winding_dict.get('drive_winding_resistance', 1.0)
    @property
    def bearing_winding_resistance(self) -> float: return self.winding_dict.get('bearing_winding_resistance', 1.0)
    
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
    bool_3PhaseCurrentSource: bool = True
    EX: dict = field(default_factory=dict)
    
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

    def update_excitations(self, machineGeometry):
        # This logic is based on the user-provided snippet
        p = self.p
        Qs = self.slot_count
        m = self.phase_count
        EX = self.EX
        
        # Ensure essential keys exist in EX with defaults
        EX.setdefault('DCBusVoltage', self.dc_bus_voltage)
        EX.setdefault('mm_stack_length_specified', 0.0)
        EX.setdefault('Js', self.rated_current_density)
        EX.setdefault('WindingFill', self.fill_factor)
        EX.setdefault('TORQUE_CURRENT_RATIO', self.torque_current_ratio)
        EX.setdefault('SUSPENSION_CURRENT_RATIO', self.suspension_current_ratio)
        EX.setdefault('RatedSpeed', self.rated_speed)

        RatedSpeed = self.rated_speed
        ExcitationFreqSimulated: float = RatedSpeed / 60 * p
        EX['ExcitationFreqSimulated'] = ExcitationFreqSimulated
        V_stator_phase_voltage_amp = math.sqrt(2) * EX['DCBusVoltage'] / (math.sqrt(3) if self.is_wye_connection else 1.0) 
        V_desired_emf_Em = 0.95 * V_stator_phase_voltage_amp
        alpha_i = 2.0/math.pi
        T_air_gap_flux_density_Bg_guessed = 0.9 # T
        mm_stack_length_specified = EX['mm_stack_length_specified']
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

        EX['no_series_coil_turns_N'] = no_series_coil_turns_N
        self.series_turns = no_series_coil_turns_N
        
        EX['DriveW_zQ'] = no_conductors_per_slot_zQ = 2* m * no_series_coil_turns_N / Qs * self.wily.number_of_parallel_branch
        self.drive_winding_conductors_per_slot = no_conductors_per_slot_zQ
        
        EX['BeariW_zQ'] = EX['DriveW_zQ'] if self.wily.bool_DPNVorSEPA == True else EX['DriveW_zQ'] / EX['TORQUE_CURRENT_RATIO'] * EX['SUSPENSION_CURRENT_RATIO']
        self.bearing_winding_conductors_per_slot = EX['BeariW_zQ']

        ''' Excitations Consiering Thermal Capability Limit (Simple) '''
        mm_r_sy = machineGeometry.r_stator_outer - machineGeometry.d_stator_yoke
        r_si = machineGeometry.r_rotor_outer + machineGeometry.d_air_gap
        mm_r_ss = r_si + machineGeometry.d_tooth_shoe
        EX['mm2_slot_area'] = (math.pi*(mm_r_sy**2 - mm_r_ss**2) / Qs - machineGeometry.w_tooth * machineGeometry.d_tooth)
        
        EX['CurrentAmp_in_the_slot']   = EX['mm2_slot_area'] * 1e-6 * EX['Js'] * EX['WindingFill'] * math.sqrt(2)
        
        EX['CurrentAmp_per_conductor'] = EX['CurrentAmp_in_the_slot'] / EX['DriveW_zQ'] if EX['DriveW_zQ'] > 0 else 0
        
        EX['CurrentAmp_per_phase']     = EX['CurrentAmp_per_conductor'] * self.wily.number_of_parallel_branch
        
        EX['DriveW_CurrentAmp'] = EX['TORQUE_CURRENT_RATIO']     * EX['CurrentAmp_per_phase']
        
        EX['BeariW_CurrentAmp'] = EX['SUSPENSION_CURRENT_RATIO'] * EX['CurrentAmp_per_phase']
        
        EX['slot_current_utilizing_ratio_for_torque'] = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase'] if EX['CurrentAmp_per_phase'] > 0 else 0

        self.slot_area = EX['mm2_slot_area']
        self.slot_current_amplitude = EX['CurrentAmp_in_the_slot']
        self.conductor_current_amplitude = EX['CurrentAmp_per_conductor']
        self.phase_current_amplitude = EX['CurrentAmp_per_phase']
        self.drive_winding_current = EX['DriveW_CurrentAmp']
        self.bearing_winding_current = EX['BeariW_CurrentAmp']
        self.torque_current_utilization_ratio = EX['slot_current_utilizing_ratio_for_torque']

        EX['InitialRotationAngle'] = self.get_InitialRotationAngle(machineGeometry)

        Rout = machineGeometry.r_rotor_outer
        Rin  = Rout - machineGeometry.d_magnet
        deg_alpha_rp = 360 / (2*self.p)
        EX['mm2_magnet_area'] = machineGeometry.deg_alpha_rm/deg_alpha_rp * math.pi*(Rout**2 - Rin**2)
        self.magnet_area = EX['mm2_magnet_area']

    def sync(self, machineGeometry, materials, l_stack):
        """Update derived performance parameters using geometry and materials."""
        # 0.1 Update Materials in EX
        self.EX['RotorCore_Material'] = materials.rotor_steel
        self.EX['StatorCore_Material'] = materials.stator_steel
        self.EX['LaminationFactor'] = materials.steel_stack_factor
        self.EX['Magnet_Name'] = materials.magnet_grade
        self.EX['Magnet_Temperature'] = 20 # Default
        self.EX['Magnet_StartAngle'] = 0.5 * machineGeometry.deg_alpha_rm # Default
        self.EX['mm_stack_length_specified'] = l_stack
        
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
        self.update_excitations(machineGeometry)

        m = materials
        
        # 1. Geometry-derived values
        r_si = machineGeometry.r_rotor_outer + machineGeometry.d_air_gap
        r_ro = machineGeometry.r_rotor_outer
        hm = machineGeometry.d_magnet
        gap_dist = machineGeometry.d_air_gap
        
        # 2. Bg Estimation
        br = m.magnet_br
        b_gap = br * hm / (hm + 1.05 * gap_dist) if (hm + 1.05 * gap_dist) > 0 else 0
        
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
        self.EX['DriveW_Rs'] = self.phase_resistance
        self.EX['BeariW_Rs'] = self.phase_resistance
        
        # 5. Back-EMF
        pole_count = self.pole_count
        self.excitation_frequency = (self.rated_speed * (pole_count / 2)) / 60.0
        
        area_pole = (2 * math.pi * r_si * l_stack) / pole_count
        estimated_flux = b_gap * area_pole * 1e-6 
        angular_speed = self.rated_speed * (2 * math.pi / 60.0)
        winding_factor = 0.95 
        self.estimated_back_emf = estimated_flux * n_series * winding_factor * angular_speed

if __name__ == "__main__":
    w = MachineWinding()
    print("MachineWinding Debug:")
    print(f"Num Slots: {w.slot_count}, Num Poles: {w.pole_count}")
    # Partial sync test would require mock geometry/materials
