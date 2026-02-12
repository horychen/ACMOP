from dataclasses import dataclass
import math

@dataclass
class MachineWinding:
    # --- Basic Specifications ---
    phase_count: int = 3 # (num_phases)
    slot_count: int = 12 # (num_slots)
    pole_count: int = 10 # (num_poles)
    coil_pitch: int = 1 # (coil_pitch_y)
    conductors_per_slot: int = 42 # (conductors_per_slot)
    wire_diameter: float = 0.21 # (wire_diameter)
    connection_type: str = "Wye" # (connection)
    rated_speed: float = 3000.0 # (rated_speed_rpm)
    rated_current_density: float = 5.0 # [A/mm^2] (rated_current_density_Js)
    parallel_branch_count: int = 1 # (number_of_parallel_branch)
    
    # --- Expanded Excitation & Thermal ---
    rated_power: float = 0.0 # [W] (RatedPower)
    dc_bus_voltage: float = 0.0 # [V] (DCBusVoltage)
    fill_factor: float = 0.4 # (WindingFill)
    excitation_frequency: float = 0.0 # [Hz] (ExcitationFreqSimulated)
    is_wye_connection: bool = True # (True: Wye) (bool_WyeConnectOrDeltaConnect)
    
    # Bearingless specific
    torque_current_ratio: float = 1.0 # (TORQUE_CURRENT_RATIO)
    suspension_current_ratio: float = 0.0 # (SUSPENSION_CURRENT_RATIO)
    drive_winding_resistance: float = 0.0 # [Ohm] (DriveW_Rs)
    bearing_winding_resistance: float = 0.0 # [Ohm] (BeariW_Rs)
    
    # Performance & Design Targets
    series_turns: int = 0 # (no_series_coil_turns_N)
    drive_winding_conductors_per_slot: int = 0 # (DriveW_zQ)
    bearing_winding_conductors_per_slot: int = 0 # (BeariW_zQ)
    slot_area: float = 0.0 # (mm2_slot_area)
    slot_current_amplitude: float = 0.0 # (CurrentAmp_in_the_slot)
    conductor_current_amplitude: float = 0.0 # (CurrentAmp_per_conductor)
    phase_current_amplitude: float = 0.0 # (CurrentAmp_per_phase)
    drive_winding_current: float = 0.0 # (DriveW_CurrentAmp)
    bearing_winding_current: float = 0.0 # (BeariW_CurrentAmp)
    torque_current_utilization_ratio: float = 0.0 # (slot_current_utilizing_ratio_for_torque)
    magnet_area: float = 0.0 # (mm2_magnet_area)
    initial_rotation_angle: float = 0.0 # (InitialRotationAngle)

    # Derived results
    phase_resistance: float = 0.0
    estimated_back_emf: float = 0.0
    slot_current_at: float = 0.0

    def sync(self, geometry_points, materials, l_stack):
        """Update derived performance parameters using geometry and materials."""
        # geometry_points is expected to be an AllPoints.HP or similar dict
        g = geometry_points
        m = materials
        
        # 1. Geometry-derived values
        r_si = g[4][0]
        r_ro = g[3][0]
        hm = g[3][0] - g[2][0]
        gap_dist = r_si - r_ro
        
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
