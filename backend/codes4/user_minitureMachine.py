from dataclasses import dataclass, field
from typing import List, Optional
import math
import sys
import os

# Ensure the current directory is in the path to import from sibling files if needed
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

@dataclass
class MachineDesignInputOri:
    # Winding
    m: int = 3
    Qs: int = 12
    p: int = 4
    ps: int = 5
    coil_pitch_y: int = 1

    # Nameplate
    RatedPower: float = 50e3 # W
    RatedSpeed: float = 30000 # rpm

    # Excitation
    bool_WyeConnectOrDeltaConnect: bool = True
    bool_weHavePlentyVoltage: bool = True
    Temperature: float = 75
    TORQUE_CURRENT_RATIO: float = 0.95
    SUSPENSION_CURRENT_RATIO: float = 0.05
    DCBusVoltage: float = 600
    Js: float = 4e6
    WindingFill: float = 0.3882
    DriveW_Rs: float = 1.0 # [Ohm]
    BeariW_Rs: float = 1.0 # [Ohm]

    # Materials
    SteelMaterial: str = "M-19 Steel Gauge-29"
    Magnet_Name: str = "Arnold/Reversible/N40H"
    LaminationFactor: float = 95
    StatorCore_Material: Optional[str] = None
    RotorCore_Material: Optional[str] = None

    # Geometry (Fixed/Initial)
    mm_stack_length_specified: float = 50 # mm
    SR: float = (8+2*0.15) / 13 # 1.0 - 0.35
    mm_r_so: float = 123.5
    mm_d_mech_air_gap: float = 0.5

    # Machine Geometry
    bool_PermanentMagnet: bool = True
    bool_StatorSlotClosed: bool = True
    bool_RotorNotched: bool = True

    # Bounds
    split_ratio_r_si_slash_r_so_bounds: List[float] = field(default_factory=lambda: [0.55, 0.75])
    yoke_split_ratio_bounds: List[float] = field(default_factory=lambda: [0.15, 0.35])
    tooth_split_ratio_at_middle_slot_bounds: List[float] = field(default_factory=lambda: [0.4, 0.75])

    def __post_init__(self):
        if self.StatorCore_Material is None:
            self.StatorCore_Material = self.SteelMaterial
        if self.RotorCore_Material is None:
            self.RotorCore_Material = self.SteelMaterial


@dataclass(frozen=True)
class MaterialSpecs:
    """Specifications for magnetic and conductive materials."""
    magnet_grade: str = "N42SH"
    magnet_br: float = 1.3  # Remanence [T]
    magnet_h_cj: float = 1592.0  # Coercivity [kA/m] (SH grade >= 20kOe)
    magnet_temp_max: float = 150.0  # [°C]
    
    stator_steel: str = "20JNEH1200"
    steel_thickness: float = 0.20  # [mm]
    steel_stack_factor: float = 0.95
    steel_max_flux_density: float = 1.9  # [T] (Saturation knee)

@dataclass
class GeometrySpecs:
    """Physical dimensions and tolerances of the 13mm actuator."""
    d_stator_outer: float = 13.0  # [mm]
    d_rotor_outer: float = 8.0   # [mm]
    d_shaft: float = 2.0         # [mm]
    magnet_thickness: float = 3.0 # [mm]
    air_gap: float = 0.15        # [mm]
    stack_length_options: List[float] = field(default_factory=lambda: [6.0, 10.0, 16.0]) # [mm]

    # Stator Tooth Geometry
    tooth_width: float = 1.2  # [mm]
    tooth_depth: float = 2.0  # [mm]
    tooth_shape_options: List[str] = field(default_factory=lambda: ['wide-open', 'semi-open', 'closed'])
    
    # Manufacturing Tolerance
    nominal_eccentricity: float = 0.05  # [mm] (5-絲)

    # Rotor Geometry
    has_rotor_magnet: bool = True
    has_rotor_notch: bool = False

    def __post_init__(self):
        self.stack_length = self.stack_length_options[-1]
        self.tooth_shape = self.tooth_shape_options[-1]

@dataclass
class WindingSpecs:
    """Winding and electrical loading parameters."""
    num_slots: int = 12
    num_poles: int = 10
    coil_pitch_y: int = 1
    winding_type: str = "Concentrated Double-Layer"
    winding_factor: float = 0.933  # k_w for 12S10P
    
    conductors_per_slot: int = 42  # z_Q
    wire_diameter_with_insulation: float = 0.21  # [mm]
    wire_gauge_approx: str = "AWG 32"
    
    # Current Density
    rated_current_density: float = 14.0  # [A/mm^2] (Peak/Pulse)

@dataclass
class PerformanceTargets:
    """Calculated physics limits and performance goals."""
    magnetic_loading_target: float = 1.14  # B_g [T]
    maxwell_stress_radial: float = 520.0  # [kPa]
    torque_constant_kt: float = 38.2      # [mNm/A]
    ump_threshold_max: float = 25.0       # [N]
    temp_limit: float = 120.0             # [°C] (Class F operating limit)

@dataclass
class MotorSpecs:
    """The master configuration for the JIAHAO-DEX-13."""
    geometry: GeometrySpecs = field(default_factory=GeometrySpecs)
    materials: MaterialSpecs = field(default_factory=MaterialSpecs)
    winding: WindingSpecs = field(default_factory=WindingSpecs)
    targets: PerformanceTargets = field(default_factory=PerformanceTargets)
    
    def get_eccentricity_ratio(self) -> float:
        """Calculate epsilon (e/g)."""
        return self.geometry.nominal_eccentricity / self.geometry.air_gap

    def calculate_copper_area(self) -> float:
        """Calculate total copper cross-section area in mm^2."""
        # Assuming 0.01mm insulation thickness
        d_bare = self.winding.wire_diameter_with_insulation - 0.02
        return (math.pi * (d_bare**2) / 4.0) * self.winding.conductors_per_slot

def convert_to_machine_design_input(specs: MotorSpecs) -> MachineDesignInputOri:
    """
    Convert MotorSpecs to MachineDesignInputOri for machine design workflow.
    Fills in known parameters and uses placeholders for others.
    """
    # Geometry Calculations
    r_so = specs.geometry.d_stator_outer / 2.0
    r_ro = specs.geometry.d_rotor_outer / 2.0
    # Assuming air gap is mechanical air gap and no sleeve for now
    r_si = r_ro + specs.geometry.air_gap 
    SR = r_si / r_so

    return MachineDesignInputOri(
        # Winding
        Qs = specs.winding.num_slots,
        p = specs.winding.num_poles // 2,
        coil_pitch_y = specs.winding.coil_pitch_y,
        m = 3, # Placeholder: Phase number (Usually 3)
        ps = 4, # Placeholder: Suspension pole pair number (For bearingless motors)

        # Nameplate
        RatedPower = 10.0, # Placeholder: [W]
        RatedSpeed = 10000.0, # Placeholder: [rpm]

        # Excitation
        Temperature = 60.0, # Placeholder: Operating temperature [C] (Limit is specs.targets.temp_limit)
        DCBusVoltage = 24.0, # Placeholder: [V]
        Js = specs.winding.rated_current_density * 1e6, # Convert [A/mm^2] to [A/m^2]
        WindingFill = 0.4, # Placeholder: Slot fill factor
        DriveW_Rs = 1.0, # Placeholder: Drive winding resistance [Ohm]
        BeariW_Rs = 1.0, # Placeholder: Bearing winding resistance [Ohm]
        bool_WyeConnectOrDeltaConnect = True, # Placeholder: Connection type
        bool_weHavePlentyVoltage = True, # Placeholder: Voltage constraint assumption
        TORQUE_CURRENT_RATIO = 1.0, # Placeholder
        SUSPENSION_CURRENT_RATIO = 0.0, # Placeholder

        # Materials
        SteelMaterial = specs.materials.stator_steel,
        Magnet_Name = specs.materials.magnet_grade,
        LaminationFactor = specs.materials.steel_stack_factor * 100, # Convert to percentage
        StatorCore_Material = specs.materials.stator_steel,
        RotorCore_Material = specs.materials.stator_steel, # Assuming same material for rotor back iron

        # Geometry (Fixed/Initial)
        mm_stack_length_specified = specs.geometry.stack_length_options[-1], # Defaulting to middle option (16.0mm)
        SR = SR,
        mm_r_so = r_so,
        mm_d_mech_air_gap = specs.geometry.air_gap,

        # Machine Geometry (from GeometrySpecs)
        bool_PermanentMagnet = specs.geometry.bool_PermanentMagnet,
        bool_StatorSlotClosed = specs.geometry.bool_StatorSlotClosed,
        bool_RotorNotched = specs.geometry.bool_RotorNotched,

        # Bounds (Using defaults from MachineDesignInputOri)
    )

# Example Usage:
if __name__ == "__main__":
    dex13 = MotorSpecs()
    
    print(f"Project: JIAHAO-DEX-13 Automation")
    print(f"Eccentricity Ratio: {dex13.get_eccentricity_ratio():.2f}")
    print(f"Total Slot Cu Area: {dex13.calculate_copper_area():.3f} mm^2")
    
    if dex13.get_eccentricity_ratio() > 0.3:
        print("WARNING: High UMP Risk Detected. Mechanical stiffness check mandatory.")

    # Convert to MachineDesignInputOri
    machine_input = convert_to_machine_design_input(dex13)
    print("\n--- Converted Machine Design Input ---")
    print(f"Stator Outer Radius (mm_r_so): {machine_input.mm_r_so:.2f} mm")
    print(f"Split Ratio (SR): {machine_input.SR:.4f}")
    print(f"Stack Length: {machine_input.mm_stack_length_specified} mm")
    print(f"Pole Pairs (p): {machine_input.p}")
    print(f"Slots (Qs): {machine_input.Qs}")
    print(f"Current Density (Js): {machine_input.Js:.2e} A/m^2")
    print(f"Steel Material: {machine_input.SteelMaterial}")


