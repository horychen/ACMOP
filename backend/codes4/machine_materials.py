from dataclasses import dataclass

@dataclass
class MachineMaterial:
    magnet_grade: str = "N42SH"
    magnet_br: float = 1.3
    magnet_h_cj: float = 1592.0
    stator_steel: str = "20JNEH1200"
    rotor_steel: str = "20JNEH1200"
    steel_stack_factor: float = 0.95
    steel_max_flux_density: float = 1.9
    copper_fill_factor: float = 0.4

if __name__ == "__main__":
    m = MachineMaterial()
    print("MachineMaterial Debug:")
    print(f"Magnet Grade: {m.magnet_grade}")
    print(f"Stator Steel: {m.stator_steel}")
