from dataclasses import dataclass
import builtins

@dataclass
class MachineMaterial:
    magnet_grade: str = "N42SH"
    magnet_br: float = 1.3
    magnet_h_cj: float = 1592.0
    d_magnet: float = 3.0
    stator_steel: str = "20JNEH1200"
    rotor_steel: str = "20JNEH1200"
    steel_stack_factor: float = 0.95
    steel_max_flux_density: float = 1.9
    l_stack: float = 16.0 # [mm]
    l_lamination: float = 0.2 # [mm]

    # copper_fill_factor: float = 0.4

    def __post_init__(self):
        self.lamination_count = self.l_stack * self.steel_stack_factor / self.l_lamination
        if builtins.verbose:
            print('- 设计目标：我们的目标是二维场电机的电磁设计与分析')
            print(f'- 永磁体{self.magnet_grade}的厚度{self.d_magnet=}，退磁目标')
            print(f'- 硅钢片{self.stator_steel}叠长，影响电机的aspect ratio，铁芯由：{self.lamination_count=}组成')
