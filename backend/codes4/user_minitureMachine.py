from dataclasses import dataclass, field
from typing import List, Optional
import math
import sys
import os

# Ensure the current directory is in the path to import from sibling files if needed
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

@dataclass
class Parameter:
    name: str
    type: str  # 'free', 'fixed', 'derived'
    value: float
    bounds: Optional[List[float]] = None
    unit: str = 'mm'

    def __post_init__(self):
        # Ensure value is within bounds if free
        if self.type == 'free' and self.bounds and len(self.bounds) == 2:
            if self.value < self.bounds[0]: self.value = self.bounds[0]
            if self.value > self.bounds[1]: self.value = self.bounds[1]


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
    """Physical dimensions and exploration space using structured Parameters."""
    d_stator_outer: Parameter = field(default_factory=lambda: Parameter("Stator OD", "fixed", 13.0))
    d_rotor_outer: Parameter = field(default_factory=lambda: Parameter("Rotor OD", "free", 8.0, [6.0, 10.0]))
    d_shaft: Parameter = field(default_factory=lambda: Parameter("Shaft", "fixed", 2.0))
    magnet_thickness: Parameter = field(default_factory=lambda: Parameter("Magnet thickness", "free", 3.0, [1.0, 5.0]))
    air_gap: Parameter = field(default_factory=lambda: Parameter("Air gap", "fixed", 0.15))
    
    # Stator Tooth Geometry
    tooth_width: Parameter = field(default_factory=lambda: Parameter("Tooth width", "free", 1.2, [0.8, 2.5]))
    tooth_depth: Parameter = field(default_factory=lambda: Parameter("Tooth depth", "free", 2.0, [1.0, 4.0]))
    tooth_shoe_depth: Parameter = field(default_factory=lambda: Parameter("Tooth shoe", "fixed", 0.5))
    
    tooth_shape: str = "closed"  # Options: "open", "semi-closed", "closed"
    tooth_shape_options: List[str] = field(default_factory=lambda: ["open", "semi-closed", "closed"])

    # Options for dropdowns
    stack_length: Parameter = field(default_factory=lambda: Parameter("Stack length", "fixed", 16.0))
    stack_length_options: List[float] = field(default_factory=lambda: [6.0, 10.0, 16.0])

    # Derived parameter example
    split_ratio: Parameter = field(default_factory=lambda: Parameter("Split ratio", "derived", 0.615))
    
    def update_derived(self):
        """Update derived parameters based on free/fixed values."""
        self.split_ratio.value = self.d_rotor_outer.value / self.d_stator_outer.value

@dataclass
class WindingSpecs:
    """Winding, electrical loading, and Back-EMF estimation."""
    num_slots: int = 12
    num_poles: int = 10
    coil_pitch_y: int = 1
    
    conductors_per_slot: int = 42  # z_Q
    wire_diameter_with_insulation: float = 0.21  # [mm]
    
    # Excitation & Nameplate
    rated_speed: float = 10000.0  # [rpm]
    dc_bus_voltage: float = 24.0  # [V]
    rated_current_density: float = 14.0  # [A/mm^2]
    
    winding_factor: float = 0.933  # k_w for 12S10P
    
    def calculate_suggested_turns(self, flux_per_pole: float) -> int:
        """
        Estimate series turns N based on V_dc and speed.
        V_emf = 4.44 * f * k_w * N * Phi
        """
        if flux_per_pole <= 0: return 0
        freq = (self.rated_speed / 60.0) * (self.num_poles / 2.0)
        # Target Back-EMF is ~90% of phase voltage
        v_phase_max = self.dc_bus_voltage / math.sqrt(3) # Star connection
        v_target = 0.9 * v_phase_max
        
        # Simplified N calculation
        n_series = v_target / (2 * math.pi * freq * self.winding_factor * flux_per_pole)
        return int(round(n_series))

    def estimate_back_emf(self, flux_per_pole: float) -> float:
        """Calculate peak phase Back-EMF [V]."""
        freq = (self.rated_speed / 60.0) * (self.num_poles / 2.0)
        n_series = (self.conductors_per_slot * self.num_slots) / (2 * 3) # Simplified for 3-phase
        return 2 * math.pi * freq * self.winding_factor * n_series * flux_per_pole

@dataclass
class PerformanceTargets:
    """Optimization objectives and FEA result containers."""
    torque_target: float = 50.0 # [mNm]
    efficiency_target: float = 0.85
    
    # FEA Post-processing results (Placeholders for time-domain waves)
    last_fea_torque_wave: List[float] = field(default_factory=list)
    last_fea_time_steps: List[float] = field(default_factory=list)
    torque_ripple: float = 0.0
    thd_voltage: float = 0.0

    # Meta Data
    individual_name_format: str = 'gen-0-ind-0'
    product_name: str = 'miniture_machine'

    # FEA Config
    select_FEA_tool: str = 'JMAG Designer' # FEMM
    select_fea_config_dict: str = '#0301 JMAG Non-Bearingless'
    fea_config_dict: dict = field(default_factory=lambda: {
        "circuit.TORQUE_CURRENT_RATIO" : 1.0,
        "circuit.SUSPENSION_CURRENT_RATIO" : 0.0,

        "local_sensitivity_analysis": False,
        "local_sensitivity_analysis_number_of_variants": 20,

        "bool_post_processing": False,
        "delete_results_after_calculation": False, 

        "designer.Show": True,
        "designer.MultipleCPUs": False,
        "designer.OnlyTableResults": False, 
        "designer.Restart": False, 
        "designer.JMAG_Scheduler": False, 

        "designer.max_nonlinear_iteration": 50,
        "designer.AddIronLossCondition": True,
        "designer.number_of_steps_1stTSS": 24, 
        "designer.number_of_steps_2ndTSS": 32, 
        "designer.StepPerCycle_3rdTSS": 64, 
        "designer.number_cycles_in_1stTSS": 3, 
        "designer.number_cycles_in_2ndTSS": 0.5, 
        "designer.number_cycles_in_3rdTSS": 0, 
        "designer.number_cycles_prolonged": 0, 
        "designer.TranRef-StepPerCycle": 64, 

        "designer.CircumferentialDivision": 720,
        "designer.meshSize_Magnet": 8,
        "designer.meshSize_Shaft": 40,
        "designer.meshSizeAir": 4.0,
        "designer.meshSize_General": 10.0,

        "femm.Coarse_Mesh": False,
        "femm.deg_per_step": None,
        "moo.popsize": 78,
        "moo.apply_constraints": False,
        "moo.fitness_OA": "TorqueDensityOverSquireRootCopperLoss",
        "moo.fitness_OB": "Efficiency",
        "moo.fitness_OC": "TorqueRipple",
        "moo.fitness_OD": None
    })
    bool_jmagDeleteResultsAfterCalculation: bool = False

    # Optimization
    generation: int = 0
    counter: int = 0
    counter_fitness_called: int = 0
    counter_fitness_return: int = 0
    # toolJd: JMAG.JMAG = None

@dataclass
class MotorSpecs:
    """The master configuration for the JIAHAO-DEX-13."""
    geometry: GeometrySpecs = field(default_factory=GeometrySpecs)
    materials: MaterialSpecs = field(default_factory=MaterialSpecs)
    winding: WindingSpecs = field(default_factory=WindingSpecs)
    targets: PerformanceTargets = field(default_factory=PerformanceTargets)
    
    def calculate_copper_area(self) -> float:
        """Calculate total copper area in one slot [mm^2]."""
        wire_radius = self.winding.wire_diameter_with_insulation / 2.0
        return self.winding.conductors_per_slot * (math.pi * wire_radius**2)

    def validate_inputs(self, step: str = "geometry") -> dict:
        """Validate inputs and return physics metrics for the given step."""
        # Sync derived parameters first
        self.geometry.update_derived()
        errors = []
        warnings = []
        metrics = {}
        
        if step == "materials":
            m = self.materials
            if m.magnet_br < 0 or m.magnet_br > 2.0:
                errors.append("Magnet Br out of physical range (0-2T)")
            if m.steel_max_flux_density < 1.0:
                warnings.append("Saturation target B_sat is unusually low")
            
            # (BH)max estimate [MGOe]
            mu0 = 4 * math.pi * 1e-7
            mu_rec = 1.05
            bh_max_si = (m.magnet_br**2) / (4 * mu0 * mu_rec) # [J/m^3]
            bh_max_mgoe = bh_max_si / 7957.7 # Convert to MGOe approx
            
            metrics["Magnet Br"] = f"{m.magnet_br} T"
            metrics["Magnet (BH)max"] = f"{bh_max_mgoe:.1f} MGOe"
            metrics["Steel B_sat"] = f"{m.steel_max_flux_density} T"
            metrics["Steel Grade"] = m.stator_steel

        elif step == "geometry":
            g = self.geometry
            if g.d_rotor_outer.value >= g.d_stator_outer.value:
                errors.append("Rotor OD must be smaller than Stator OD")
            
            r_si = (g.d_rotor_outer.value / 2.0) + g.air_gap.value
            if r_si >= (g.d_stator_outer.value / 2.0):
                errors.append("Air gap extends beyond stator outer diameter!")

            # Split Ratio
            split_ratio = g.split_ratio.value # Using pre-updated derived param
            metrics["Split Ratio"] = f"{split_ratio:.3f}"
            
            # Equivalent Air Gap
            mu_rec = 1.05
            metrics["Equivalent Gap"] = f"{g.air_gap.value + g.magnet_thickness.value/mu_rec:.3f} mm"

            # Slot & Area Ratios
            stator_inner_radius = r_si
            slot_bottom_radius = stator_inner_radius + g.tooth_depth.value
            avg_slot_width = (math.pi * (stator_inner_radius + slot_bottom_radius) / self.winding.num_slots) - g.tooth_width.value
            slot_area = avg_slot_width * g.tooth_depth.value
            copper_area = self.calculate_copper_area()
            fill_factor = copper_area / slot_area if slot_area > 0 else 0
            
            metrics["Slot Area"] = f"{slot_area:.2f} mm²"
            metrics["Copper Area"] = f"{copper_area:.2f} mm²"
            metrics["Fill Factor"] = f"{fill_factor*100:.1f}%"
            metrics["Conductors/Slot"] = str(self.winding.conductors_per_slot)
            
            # Volume & Weight Est (Very rough)
            vol_stator = (math.pi * (g.d_stator_outer.value**2 - (2*r_si)**2) / 4.0) * g.stack_length.value * 1e-3 # [cm^3]
            weight_stator = vol_stator * 7.7 * 1e-3 # [kg] (Iron + Winding approx)
            metrics["Est. Core Weight"] = f"{weight_stator*1000:.1f} g"

        elif step == "winding":
            w = self.winding
            g_specs = self.geometry
            if w.num_slots % 3 != 0:
                errors.append("Slot number must be multiple of 3")
            
            # 1. Bg Estimation
            br = self.materials.magnet_br
            hm = g_specs.magnet_thickness.value
            gap_dist = g_specs.air_gap.value
            b_gap = br * hm / (hm + 1.05 * gap_dist) if (hm + 1.05 * gap_dist) > 0 else 0
            metrics["Est. Gap Flux Bg"] = f"{b_gap:.3f} T"

            # 2. Electrical Metrics
            wire_area = math.pi * (w.wire_diameter_with_insulation/2)**2
            conductor_current = w.rated_current_density * wire_area
            slot_current = conductor_current * w.conductors_per_slot
            
            # Resistance Estimation (Phase)
            r_si = (g_specs.d_rotor_outer.value / 2.0) + g_specs.air_gap.value
            end_winding = math.pi * (2*r_si) / w.num_slots * w.coil_pitch_y
            turn_length = 2 * (g_specs.stack_length.value + end_winding) * 1e-3 # [m]
            n_series = (w.conductors_per_slot * w.num_slots) / 6 # Turns per phase
            rho_copper = 1.72e-8 # [Ohm-m]
            resistance = rho_copper * (turn_length * n_series) / (wire_area * 1e-6)
            
            metrics["Magnet Br"] = f"{self.materials.magnet_br} T"
            metrics["Phase Resistance"] = f"{resistance:.3f} Ω"
            metrics["Conductor Current"] = f"{conductor_current:.2f} A"
            metrics["Total Slot Current"] = f"{slot_current:.1f} A-t"
            
            # Copper Loss
            p_cu = 3 * (conductor_current**2) * resistance
            metrics["Copper Loss Pcu"] = f"{p_cu:.2f} W"

            # Fill Factor (repeated for context)
            stator_inner_radius = (g_specs.d_rotor_outer.value / 2.0) + g_specs.air_gap.value
            slot_bottom_radius = stator_inner_radius + g_specs.tooth_depth.value
            avg_slot_width = (math.pi * (stator_inner_radius + slot_bottom_radius) / w.num_slots) - g_specs.tooth_width.value
            slot_area = avg_slot_width * g_specs.tooth_depth.value
            copper_area = self.calculate_copper_area()
            fill_factor = copper_area / slot_area if slot_area > 0 else 0
            
            metrics["Slot Area"] = f"{slot_area:.2f} mm²"
            metrics["Winding Fill Factor"] = f"{fill_factor*100:.1f}%"

            # Back-EMF
            area_pole = (2 * math.pi * r_si * g_specs.stack_length.value) / w.num_poles
            estimated_flux = b_gap * area_pole * 1e-6 
            emf = w.estimate_back_emf(estimated_flux)
            metrics["Estimated Back-EMF"] = f"{emf:.2f} V"

        elif step == "targets":
            metrics["Torque Target"] = f"{self.targets.torque_target} mNm"
            # Power Est
            power = (self.targets.torque_target * 1e-3) * (self.winding.rated_speed * 2 * math.pi / 60.0)
            metrics["Rated Power"] = f"{power:.2f} W"

        status = "pass"
        if errors: status = "fail"
        elif warnings: status = "warn"
            
        return {
            "status": status, 
            "errors": errors, 
            "warnings": warnings, 
            "metrics": metrics,
            "fill_factor": metrics.get("Fill Factor", "0%") # For backward compatibility
        }

    def validate_outputs(self, step: str, results: dict) -> dict:
        """
        Validate outputs/results of a specific step.
        """
        # Placeholder for future steps
        return {"status": "pass", "errors": [], "warnings": []}

    def calculate_copper_area(self) -> float:
        """Calculate total copper cross-section area in mm^2."""
        # Assuming 0.01mm insulation thickness
        d_bare = self.winding.wire_diameter_with_insulation - 0.02
        return (math.pi * (d_bare**2) / 4.0) * self.winding.conductors_per_slot




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
        bool_PermanentMagnet = specs.geometry.has_rotor_magnet,
        bool_StatorSlotClosed = specs.geometry.tooth_shape=='closed',
        bool_RotorNotched = specs.geometry.has_rotor_notch,

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


