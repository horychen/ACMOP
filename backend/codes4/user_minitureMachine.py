import math
from dataclasses import dataclass, field, fields
from typing import List, Dict, Any, Optional
import sys
import os

__version__ = "1.0.2_debug_renames"
print(f"DEBUG: Loading user_minitureMachine.py version {__version__}")

# Ensure the current directory is in the path to import from sibling files if needed
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modern_machine_designer_utility import Parameter, Geometry


@dataclass
class MaterialSpecs:
    """Specifications for magnetic and conductive materials."""
    magnet_grade: str = "N42SH"
    magnet_br: float = 1.3  # Remanence [T]
    magnet_h_cj: float = 1592.0  # Coercivity [kA/m] (SH grade >= 20kOe)
    magnet_temp_max: float = 150.0  # [°C]
    magnet_temperature: float = 20.0 # [°C]
    
    stator_steel: str = "20JNEH1200"
    stator_core_material: str = "20JNEH1200"
    rotor_core_material: str = "20JNEH1200"
    steel_thickness: float = 0.20  # [mm]
    steel_stack_factor: float = 0.95
    lamination_factor: float = 95.0 # [%]
    steel_max_flux_density: float = 1.9  # [T] (Saturation knee)



@dataclass
class ToothGeometry:
    """Base class for tooth geometry parameters."""
    shape: str = field(init=False)

@dataclass
class OpenTooth(ToothGeometry):
    shape: str = "open"
    d_tooth_shoe: Parameter = field(default_factory=lambda: Parameter("stator_tooth_shoe_depth", "fixed", 0.0))

@dataclass
class SemiClosedTooth(ToothGeometry):
    shape: str = "semi-closed"
    alpha_tooth: Parameter = field(default_factory=lambda: Parameter("stator_tooth_span_angle", "fixed", 0.0))
    d_tooth_open: Parameter = field(default_factory=lambda: Parameter("stator_tooth_open_depth", "derived"))
    alpha_tooth_open: Parameter = field(default_factory=lambda: Parameter("stator_tooth_open_angle", "derived"))
    d_tooth_shoe: Parameter = field(default_factory=lambda: Parameter("stator_tooth_shoe_depth", "fixed", 0.5))

@dataclass
class ClosedTooth(ToothGeometry):
    shape: str = "closed"
    d_tooth_shoe: Parameter = field(default_factory=lambda: Parameter("stator_tooth_shoe_depth", "fixed", 0.5))
@dataclass
class GeometrySpecs:
    """Physical dimensions and exploration space using structured Parameters."""
    # Fixed / Input
    r_stator_outer: Parameter = field(default_factory=lambda: Parameter("stator_outer_radius", "fixed", 13.0/2.0))
    r_rotor_outer: Parameter = field(default_factory=lambda: Parameter("rotor_outer_radius", "free", 8.0/2.0, [6.0/2.0, 10.0/2.0]))
    r_shaft: Parameter = field(default_factory=lambda: Parameter("shaft_radius", "fixed", 2.0/2.0))
    d_magnet: Parameter = field(default_factory=lambda: Parameter("magnet_depth", "free", 3.0, [1.0, 5.0]))
    d_air_gap: Parameter = field(default_factory=lambda: Parameter("mechanical_air_gap_depth", "fixed", 0.15))
    l_stack: Parameter = field(default_factory=lambda: Parameter("stack_length", "fixed", 16.0))
    
    # Motor specific flags
    bool_PermanentMagnet: bool = True
    bool_StatorSlotClosed: bool = True
    bool_RotorNotched: bool = True
    bool_use_sleeve: bool = False
    bool_use_shaft: bool = False
    
    tooth_shape: str = "closed" # Options: "open", "semi-closed", "closed"
    tooth_specs: ToothGeometry = field(init=False)

    # Stator Tooth Geometry
    w_tooth: Parameter = field(default_factory=lambda: Parameter("stator_tooth_width", "free", 1.2))
    d_tooth: Parameter = field(default_factory=lambda: Parameter("stator_tooth_depth", "derived"))
    d_stator_yoke: Parameter = field(default_factory=lambda: Parameter("stator_yoke_depth", "free", 0.5))

    # Rotor specific
    r_rotor_inner: Parameter = field(default_factory=lambda: Parameter("rotor_inner_radius", "fixed", 0.0))
    d_sleeve: Parameter = field(default_factory=lambda: Parameter("rotor_sleeve_depth", "fixed", 0.0))
    alpha_magnet_pole: Parameter = field(default_factory=lambda: Parameter("magnet_pole_span_angle", "fixed", 18.0))
    alpha_magnet_segment: Parameter = field(default_factory=lambda: Parameter("magnet_segment_span_angle", "derived"))
    inter_polar_iron_thickness: Parameter = field(default_factory=lambda: Parameter("inter_polar_iron_thickness", "derived"))
    inter_segment_iron_thickness: Parameter = field(default_factory=lambda: Parameter("inter_segment_iron_thickness", "fixed", 0.0))

    # Options for dropdowns
    tooth_shape_options: List[str] = field(default_factory=lambda: ["open", "semi-closed", "closed"])
    stack_length_options: List[float] = field(default_factory=lambda: [6.0, 10.0, 16.0])

    # Derived parameter example
    split_ratio: Parameter = field(default_factory=lambda: Parameter("split_ratio_r_si_slash_r_so", "derived", 0.615))
        
    p: Parameter = field(default_factory=lambda: Parameter("pole_pair_number_p", "fixed", 5.0))
    s: Parameter = field(default_factory=lambda: Parameter("number_of_magnet_segments_per_pole", "fixed", 1.0))

    def __post_init__(self):
        # Initialize tooth_specs based on tooth_shape
        if self.tooth_shape == "open":
            self.tooth_specs = OpenTooth()
        elif self.tooth_shape == "semi-closed":
            self.tooth_specs = SemiClosedTooth()
        elif self.tooth_shape == "closed":
            self.tooth_specs = ClosedTooth()
        else:
            self.tooth_specs = OpenTooth() # Fallback

        self.update_derived()

        # Import needed classes locally to avoid circular dependencies if any
        from modern_machine_designer_utility import CairoDrawer
        import CrossSectInnerNotchedRotor, CrossSectStator

        # Geometry is defined by Key Points and relations among key points.
        # geometric parameters are only used to define the key points
        self.machineGeometry = {
            "rotorCore": Geometry(name='rotorCore',
                GP={
                    'r_stator_outer': self.r_stator_outer,
                    'r_rotor_outer': self.r_rotor_outer,
                    'd_magnet': self.d_magnet,
                    'r_rotor_inner': self.r_rotor_inner,
                    'inter_polar_iron_thickness': self.inter_polar_iron_thickness,
                    'inter_segment_iron_thickness': self.inter_segment_iron_thickness,
                    'p': self.p,
                    's': self.s,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                        name="rotorCore",
                        color="#FE840E",
                        mm_d_pm=self.d_magnet.value,
                        deg_alpha_rm=self.alpha_magnet_pole.value,
                        deg_alpha_rs=self.alpha_magnet_pole.value if self.s.value==1 else self.alpha_magnet_segment.value,
                        mm_d_ri=self.r_rotor_outer.value - self.d_magnet.value,
                        mm_r_ri=self.r_rotor_inner.value,
                        mm_d_rp=self.inter_polar_iron_thickness.value,
                        mm_d_rs=self.inter_segment_iron_thickness.value,
                        p=int(self.p.value),
                        s=int(self.s.value)
                    ).draw(drawer, **kwargs)
                )
            ),
            "rotorMagnet": Geometry(name='rotorMagnet',
                GP={
                    'd_magnet': self.d_magnet,
                    'r_rotor_outer': self.r_rotor_outer,
                    'r_rotor_inner': self.r_rotor_inner,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet(
                        name="rotorMagnet",
                        color="#1C96E0",
                        rotorCore=CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                            mm_d_pm=self.d_magnet.value,
                            deg_alpha_rm=self.alpha_magnet_pole.value,
                            deg_alpha_rs=self.alpha_magnet_pole.value if self.s.value==1 else self.alpha_magnet_segment.value,
                            mm_d_ri=self.r_rotor_outer.value - self.d_magnet.value,
                            mm_r_ri=self.r_rotor_inner.value,
                            mm_d_rp=self.inter_polar_iron_thickness.value,
                            mm_d_rs=self.inter_segment_iron_thickness.value,
                            p=int(self.p.value),
                            s=int(self.s.value)
                        )
                    ).draw(drawer, **kwargs)
                ),
            ),
            "statorCore": Geometry(name='statorCore',
                GP={
                    'r_stator_outer': self.r_stator_outer,
                    'd_tooth_shoe': self.tooth_specs.d_tooth_shoe,
                    'd_tooth': self.d_tooth,
                    'd_stator_yoke': self.d_stator_yoke,
                    'w_tooth': self.w_tooth,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectStator.CrossSectInnerRotorClosedSlotStator(
                        name="statorCore",
                        color="#BAFA01",
                        mm_r_so=self.r_stator_outer.value,
                        mm_d_sts=self.tooth_specs.d_tooth_shoe.value,
                        mm_r_si=self.r_stator_outer.value * self.split_ratio.value,
                        mm_d_st=self.d_tooth.value,
                        mm_d_sy=self.d_stator_yoke.value,
                        mm_w_st=self.w_tooth.value,
                        Q=int(Parameter("stator_slot_number_Qs", "fixed", 12).value), 
                    ).draw(drawer, **kwargs)
                ),
            ),
            "coils": None
        }
        self.machineGeometry['coils'] = Geometry(name='coils',
            GP={
                'r_stator_outer': self.r_stator_outer,
                'd_stator_yoke': self.d_stator_yoke,
                'w_tooth': self.w_tooth,
                'd_tooth': self.d_tooth,
            },
            draw_function=lambda drawer, **kwargs: (
                CrossSectStator.CrossSectInnerRotorClosedSlotStatorWinding(
                    stator_core=self.machineGeometry['statorCore'],
                ).draw(drawer, **kwargs)
            )
        )

    def show_geometry(self, filename=None):
        from modern_machine_designer_utility import CairoDrawer
        def draw_spmsm(lw, width_in_points, height_in_points, filename='machine_geometry.svg', bool_draw_whole_model=True):
            self.drawer = drawer = CairoDrawer(width_in_points, height_in_points, filename=filename, verbose_drawing=getattr(self, 'verbose_drawing', False))

            self.machineGeometry['rotorCore'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)
            # self.machineGeometry['rotorMagnet'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)
            # self.machineGeometry['statorCore'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)
            # self.machineGeometry['coils'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)

            drawer.apply_stroke(lw=lw)
            drawer.convert_to_pdf()

        lw = 0.1 if self.r_rotor_outer.value < 15 else 0.5
        width_in_points  = self.r_stator_outer.value*2.1
        height_in_points = self.r_stator_outer.value*2.1
        bool_draw_whole_model = True

        draw_spmsm(lw, width_in_points, height_in_points, filename=filename or 'machine_geometry.svg', bool_draw_whole_model=bool_draw_whole_model)

    def get_parameter_dict(self):
        d = {}
        for f in fields(self):
            attr = getattr(self, f.name)
            if isinstance(attr, Parameter):
                d[attr.name] = attr
        
        # Add tooth_specs parameters
        if hasattr(self, 'tooth_specs'):
            for f in fields(self.tooth_specs):
                attr = getattr(self.tooth_specs, f.name)
                if isinstance(attr, Parameter):
                    d[attr.name] = attr
        return d

    _derivation_errors: List[str] = field(default_factory=list, init=False, repr=False)

    def _validate_positive(self, param_name, value, dependencies):
        """Check if a derived value is positive. If not, log dependencies for debugging and store error."""
        if value <= 0:
            msg = f"ERROR: Derived parameter '{param_name}' is non-positive (value: {value:.3f})."
            msg += " components: " + ", ".join([f"{k}={v:.3f}" for k, v in dependencies.items()])
            print("\n" + msg)
            if msg not in self._derivation_errors:
                self._derivation_errors.append(msg)
        return value

    def update_derived(self, motor_specs=None):
        """Update derived parameters based on free/fixed values."""
        self._derivation_errors = []
        # 1. Stator Inner Radius (derived from outer and split ratio)
        r_so = self.r_stator_outer.value
        r_si = r_so * self.split_ratio.value
        
        # 2. Tooth Depth
        # Logic: r_so - r_si - yoke - shoe
        val_td = r_so - r_si - self.d_stator_yoke.value - self.tooth_specs.d_tooth_shoe.value
        self.d_tooth.value = self._validate_positive("d_tooth", val_td, {
            "r_stator_outer": r_so,
            "r_si": r_si,
            "d_stator_yoke": self.d_stator_yoke.value,
            "d_tooth_shoe": self.tooth_specs.d_tooth_shoe.value
        })
        
        # 3. Rotor Outer Radius
        r_ro = r_si - self.d_air_gap.value - self.d_sleeve.value
        self.r_rotor_outer.value = r_ro
        self._validate_positive("r_rotor_outer", self.r_rotor_outer.value, {
            "r_si": r_si,
            "d_air_gap": self.d_air_gap.value,
            "d_sleeve": self.d_sleeve.value
        })
        
        # 4. Rotor Iron Depth
        self.inter_polar_iron_thickness.value = self.d_magnet.value
        
        # 5. Stator Tooth Open (OO parameters)
        if hasattr(self.tooth_specs, 'd_tooth_open'):
            self.tooth_specs.d_tooth_open.value = self.tooth_specs.d_tooth_shoe.value * 0.667
            self._validate_positive("d_tooth_open", self.tooth_specs.d_tooth_open.value, {
                "d_tooth_shoe": self.tooth_specs.d_tooth_shoe.value
            })
        
        if hasattr(self.tooth_specs, 'alpha_tooth_open') and hasattr(self.tooth_specs, 'alpha_tooth'):
            self.tooth_specs.alpha_tooth_open.value = self.tooth_specs.alpha_tooth.value * 0.5
        
        # 6. Magnet Segments
        self.alpha_magnet_segment.value = self.alpha_magnet_pole.value
        
        # Update split_ratio for reporting if it's derived from final r_si
        if r_so != 0:
            self.split_ratio.value = r_si / r_so

@dataclass
class WindingSpecs:
    """Winding, electrical loading, and Back-EMF estimation."""
    m: int = 3
    slot_count: int = 12
    pole_count: int = 10
    ps: int = 5 # suspension_pole_pair_number
    coil_pitch_y: int = 1

    @property
    def p(self): return self.pole_count // 2
    @property
    def Qs(self): return self.slot_count
    bool_DPNVorSEPA: bool = True
    number_of_parallel_branch: int = 2
    
    conductors_per_slot: int = 42  # z_Q
    wire_diameter_with_insulation: float = 0.21  # [mm]
    
    # Excitation & Nameplate
    rated_speed: float = 10000.0  # [rpm]
    rated_power: float = 50.0  # [W]
    dc_bus_voltage: float = 24.0  # [V]
    rated_current_density_Js: float = 14.0e6  # [A/m^2]
    winding_fill_factor: float = 0.4
    
    winding_factor: float = 0.933  # k_w for 12S10P
    
    # Connection and Ratios
    bool_WyeConnectOrDeltaConnect: bool = True
    torque_current_ratio: float = 1.0
    suspension_current_ratio: float = 0.0
    drive_winding_resistance: float = 1.0 # [Ohm]
    bearing_winding_resistance: float = 1.0 # [Ohm]
    
    # Dynamic values (to be populated by sync)
    wily: Any = None # Modern Winding object
    no_series_coil_turns_N: int = 0
    DriveW_zQ: float = 0.0
    BeariW_zQ: float = 0.0
    mm2_slot_area: float = 0.0
    CurrentAmp_per_phase: float = 0.0
    DriveW_CurrentAmp: float = 0.0
    BeariW_CurrentAmp: float = 0.0

    def calculate_suggested_turns(self, flux_per_pole: float) -> int:
        if flux_per_pole <= 0: return 0
        freq = (self.rated_speed / 60.0) * (self.pole_count / 2.0)
        v_phase_max = self.dc_bus_voltage / math.sqrt(3) # Star connection
        v_target = 0.95 * v_phase_max * math.sqrt(2) # Ampere value
        
        n_series = v_target / (2 * math.pi * freq * self.winding_factor * flux_per_pole)
        return int(round(n_series))

    def estimate_back_emf(self, flux_per_pole: float) -> float:
        freq = (self.rated_speed / 60.0) * (self.pole_count / 2.0)
        n_series = self.no_series_coil_turns_N
        return 2 * math.pi * freq * self.winding_factor * n_series * flux_per_pole

@dataclass
class PerformanceTargets:
    """Optimization objectives and FEA result containers."""
    torque_target: float = 50.0 # [mNm]
    efficiency_target: float = 0.85
    
    # Extra Objectives
    objectives: Dict[str, float] = field(default_factory=lambda: {
        'TorqueDensity': 0.0,
        'Cost': 0.0,
        'TorqueDensityOverSquireRootCopperLoss': 0.0,
        'Efficiency': 0.0,
        'TorqueRipple': 0.0,
        'ForceErrorMagnitude': 0.0,
        'ForceErrorAngle': 0.0,
        'IronLoss': 0.0,
        'CoggingTorque': 0.0
    })

    # FEA Post-processing results
    last_fea_torque_wave: List[float] = field(default_factory=list)
    last_fea_time_steps: List[float] = field(default_factory=list)
    torque_ripple: float = 0.0
    thd_voltage: float = 0.0
    initial_rotation_angle: float = 0.0
    mm2_magnet_area: float = 0.0
    
    # Excitation Results
    no_series_coil_turns_N: int = 0
    no_conductors_per_slot_zQ: float = 0.0
    mm2_slot_area: float = 0.0
    CurrentAmp_in_the_slot: float = 0.0
    CurrentAmp_per_conductor: float = 0.0
    CurrentAmp_per_phase: float = 0.0
    DriveW_CurrentAmp: float = 0.0
    BeariW_CurrentAmp: float = 0.0
    slot_current_utilizing_ratio_for_torque: float = 0.0

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
    
    def get_eccentricity_ratio(self) -> float:
        """Calculate eccentricity ratio (offset/airgap). Placeholder implementation."""
        # For a standard non-bearingless motor, this is 0
        return 0.0

    def validate_inputs(self, step: str = "geometry") -> dict:
        """Validate inputs and return physics metrics for the given step."""
        try:
            # Sync derived parameters first
            self.geometry.update_derived(self)
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
                # Add derivation errors found during update_derived
                if hasattr(g, "_derivation_errors") and g._derivation_errors:
                    errors.extend(g._derivation_errors)
                
                r_so = g.r_stator_outer.value
                r_ro = g.r_rotor_outer.value
                
                if r_ro >= r_so:
                    errors.append("Rotor radius must be smaller than Stator radius")
                
                r_si = r_so * g.split_ratio.value
                if r_si >= r_so:
                    errors.append("Air gap extends beyond stator outer diameter!")

                # Split Ratio
                metrics["Split Ratio"] = f"{g.split_ratio.value:.3f} (r_si/r_so)"
                
                # Equivalent Air Gap
                mu_rec = 1.05
                metrics["Equivalent Gap"] = f"{g.d_air_gap.value + g.d_magnet.value/mu_rec:.3f} mm"

                # Slot & Area Ratios
                stator_inner_radius = r_si
                slot_bottom_radius = r_si + g.d_tooth.value
                avg_slot_width = (math.pi * (stator_inner_radius + slot_bottom_radius) / self.winding.slot_count) - g.w_tooth.value
                slot_area = avg_slot_width * g.d_tooth.value
                copper_area = self.calculate_copper_area()
                fill_factor = copper_area / slot_area if slot_area > 0 else 0
                
                metrics["Slot Area"] = f"{slot_area:.2f} mm²"
                metrics["Copper Area"] = f"{copper_area:.2f} mm²"
                metrics["Fill Factor"] = f"{fill_factor*100:.1f}%"
                metrics["Conductors/Slot"] = str(self.winding.conductors_per_slot)
                
                # Volume & Weight Est (Very rough)
                vol_stator = (math.pi * (r_so**2 - r_si**2)) * g.l_stack.value * 1e-3 # [cm^3]
                weight_stator = vol_stator * 7.7 * 1e-3 # [kg]
                metrics["Est. Core Weight"] = f"{weight_stator*1000:.1f} g"

            elif step == "winding":
                w = self.winding
                g_specs = self.geometry
                if w.slot_count % 3 != 0:
                    errors.append("Slot number must be multiple of 3")
                
                # 1. Bg Estimation
                br = self.materials.magnet_br
                hm = g_specs.d_magnet.value
                gap_dist = g_specs.d_air_gap.value
                b_gap = br * hm / (hm + 1.05 * gap_dist) if (hm + 1.05 * gap_dist) > 0 else 0
                metrics["Est. Gap Flux Bg"] = f"{b_gap:.3f} T"

                # 2. Electrical Metrics
                wire_area = math.pi * (w.wire_diameter_with_insulation/2)**2
                conductor_current = w.rated_current_density_Js * wire_area * 1e-6 
                slot_current = conductor_current * w.conductors_per_slot
                
                # Resistance Estimation (Phase)
                r_si = g_specs.r_stator_outer.value * g_specs.split_ratio.value
                end_winding = math.pi * (2*r_si) / w.slot_count * w.coil_pitch_y
                n_series = (w.conductors_per_slot * w.slot_count) / (2 * w.m * w.number_of_parallel_branch)
                turn_length = 2 * (g_specs.l_stack.value + end_winding) * 1e-3 # [m]
                rho_copper = 1.72e-8 # [Ohm-m]
                resistance = rho_copper * (turn_length * n_series) / (wire_area * 1e-6)
                
                metrics["Magnet Br"] = f"{self.materials.magnet_br} T"
                metrics["Phase Resistance"] = f"{resistance:.3f} Ω"
                metrics["Conductor Current"] = f"{conductor_current:.2f} A"
                metrics["Total Slot Current"] = f"{slot_current:.1f} A-t"
                
                # Copper Loss
                p_cu = 3 * (conductor_current**2) * resistance
                metrics["Copper Loss Pcu"] = f"{p_cu:.2f} W"

                # Fill Factor
                stator_inner_radius = r_si
                slot_bottom_radius = r_si + g_specs.d_tooth.value
                avg_slot_width = (math.pi * (stator_inner_radius + slot_bottom_radius) / w.slot_count) - g_specs.w_tooth.value
                slot_area = avg_slot_width * g_specs.d_tooth.value
                copper_area = self.calculate_copper_area()
                fill_factor = copper_area / slot_area if slot_area > 0 else 0
                
                metrics["Slot Area"] = f"{slot_area:.2f} mm²"
                metrics["Winding Fill Factor"] = f"{fill_factor*100:.1f}%"

                # Back-EMF
                area_pole = (2 * math.pi * r_si * g_specs.l_stack.value) / w.pole_count
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
                "fill_factor": metrics.get("Fill Factor", "0%")
            }
        except Exception as e:
            return {
                "status": "fail",
                "errors": [f"Validation Crash: {str(e)}"],
                "warnings": [],
                "metrics": {},
                "fill_factor": "0%"
            }

    def validate_outputs(self, step: str, results: dict) -> dict:
        """Validate outputs/results of a specific step."""
        return {"status": "pass", "errors": [], "warnings": []}

    def sync(self):
        """Orchestrate initialization logic matching Modern_Machine_Designer."""
        from modern_machine_designer_utility import Winding
        
        # 1. Winding initialization
        w = self.winding
        w.wily = Winding(w.m, w.slot_count, w.pole_count // 2, w.ps, w.coil_pitch_y, bool_DPNVorSEPA=w.bool_DPNVorSEPA)

        # 2. Material Temperature Adjustment
        available_temperature_list = [-40, 20, 60, 80, 100, 120, 150, 180, 200, 220]
        self.materials.magnet_temperature = min(available_temperature_list, key=lambda x: abs(x - self.materials.magnet_temperature))

        # 3. Parameter and Geometry Setup
        g = self.geometry
        p_dict = g.get_parameter_dict()
        
        # Link parameter_dicts
        for p in p_dict.values():
            p.parameter_dict = p_dict
        
        # Derived Variables (Sync logic)
        g.update_derived(self)
        
        # 4. Turns and Excitation Logic
        V_stator_phase_voltage_amp = math.sqrt(2) * self.winding.dc_bus_voltage / math.sqrt(3) # Assume Wye
        V_desired_emf_Em = 0.95 * V_stator_phase_voltage_amp
        
        p = w.pole_count // 2
        r_si = g.r_stator_outer.value * g.split_ratio.value
        tau_p = math.pi * r_si / p
        Wb_flux = (2.0/math.pi) * 0.9 * tau_p * 1e-3 * (g.l_stack.value + 2 * g.d_air_gap.value) * 1e-3
        
        freq = (w.rated_speed / 60) * p
        no_series_coil_turns_N = V_desired_emf_Em / (2*math.pi* freq * w.wily.kw1 * Wb_flux) if Wb_flux > 0 else 0
        w.no_series_coil_turns_N = round(no_series_coil_turns_N)
        
        # Voltage priority
        SPP = w.slot_count / (2*p*w.m)
        w.no_series_coil_turns_N = min([int(p*SPP*i) for i in range(100,0,-1)], key=lambda x:abs(x - w.no_series_coil_turns_N))
        
        w.DriveW_zQ = 2 * w.m * w.no_series_coil_turns_N / w.slot_count * w.number_of_parallel_branch
        w.BeariW_zQ = w.DriveW_zQ

        # Thermal / Current Density
        r_sy = g.r_stator_outer.value - g.d_stator_yoke.value
        r_ss = r_si + g.tooth_specs.d_tooth_shoe.value
        w.mm2_slot_area = (math.pi*(r_sy**2 - r_ss**2) / w.slot_count - g.w_tooth.value * (r_sy - r_ss))
        i_slot = w.mm2_slot_area * 1e-6 * w.rated_current_density_Js * w.winding_fill_factor * math.sqrt(2)
        i_cond = i_slot / w.DriveW_zQ if w.DriveW_zQ > 0 else 0
        w.CurrentAmp_per_phase = i_cond * w.number_of_parallel_branch
        w.DriveW_CurrentAmp = 1.0 * w.CurrentAmp_per_phase 
        w.BeariW_CurrentAmp = 0.0 * w.CurrentAmp_per_phase 
        
        # Targets
        t = self.targets
        t.no_series_coil_turns_N = w.no_series_coil_turns_N
        t.no_conductors_per_slot_zQ = w.DriveW_zQ
        t.mm2_slot_area = w.mm2_slot_area
        t.CurrentAmp_in_the_slot = i_slot
        t.CurrentAmp_per_conductor = i_cond
        t.CurrentAmp_per_phase = w.CurrentAmp_per_phase
        t.DriveW_CurrentAmp = w.DriveW_CurrentAmp
        t.BeariW_CurrentAmp = w.BeariW_CurrentAmp

    def print_summary(self):
        """Print a summary of MotorSpecs."""
        print("-" * 50)
        print(f"{'Motor Specification Summary':^50}")
        print("-" * 50)
        g = self.geometry
        w = self.winding
        print(f"Geometry: {g.r_stator_outer.value*2}x{g.l_stack.value} mm, Split: {g.split_ratio.value:.3f}")
        print(f"Winding:  {w.slot_count}S/{w.pole_count}P, Turns: {w.no_series_coil_turns_N}, zQ: {w.DriveW_zQ:.1f}")
        print(f"Materials: Stator:{self.materials.stator_core_material}, Magnet:{self.materials.magnet_grade}")
        print(f"Excitation: Js:{w.rated_current_density_Js/1e6:.1f} A/mm2, Current:{w.DriveW_CurrentAmp:.1f} A")
        print("-" * 50)

    def calculate_copper_area(self) -> float:
        """Calculate total copper area in one slot [mm^2]."""
        d_bare = self.winding.wire_diameter_with_insulation - 0.02
        return self.winding.conductors_per_slot * (math.pi * (d_bare / 2.0)**2)

# Global helper for command line
def run_visualization():
    specs = MotorSpecs()
    specs.sync()
    specs.geometry.show_geometry()
    print("Geometry visualization generated as machine_geometry.svg/pdf")
    specs.print_summary()

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--viz":
        run_visualization()
    else:
        dex13 = MotorSpecs()
        print(f"Project: JIAHAO-DEX-13 Automation")
        print(f"Eccentricity Ratio: {dex13.get_eccentricity_ratio():.2f}")
        dex13.sync()
        dex13.print_summary()
