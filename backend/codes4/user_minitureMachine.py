from dataclasses import dataclass, field
from typing import List, Optional
import math

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
    air_gap: float = 0.15        # [mm]
    magnet_thickness: float = 3.0 # [mm]
    stack_length_options: List[float] = field(default_factory=lambda: [6.0, 10.0, 16.0]) # [mm]
    
    # Stator Tooth Geometry
    tooth_width: float = 2.0  # [mm]
    tooth_depth: float = 2.0  # [mm]
    tooth_shape_options: List[str] = field(default_factory=lambda: ['wide-open', 'semi-open', 'closed'])
    
    # Manufacturing Tolerance
    nominal_eccentricity: float = 0.05  # [mm] (5-絲)

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

# Example Usage:
if __name__ == "__main__":
    dex13 = MotorSpecs()
    
    print(f"Project: JIAHAO-DEX-13 Automation")
    print(f"Eccentricity Ratio: {dex13.get_eccentricity_ratio():.2f}")
    print(f"Total Slot Cu Area: {dex13.calculate_copper_area():.3f} mm^2")
    
    if dex13.get_eccentricity_ratio() > 0.3:
        print("WARNING: High UMP Risk Detected. Mechanical stiffness check mandatory.")
