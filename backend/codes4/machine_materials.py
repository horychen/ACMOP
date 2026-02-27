import builtins

class MachineMaterial:
    def __init__(self, material_dict: dict = None):
        self.material_dict = material_dict if material_dict is not None else {}
        
        defaults = {
            'magnet_grade': "N42SH",
            'magnet_bg': 1.23,
            'magnet_br': 1.3,
            'magnet_h_cj': 1592.0,
            'd_magnet': 3.0,
            'stator_steel': "20JNEH1200",
            'rotor_steel': "20JNEH1200",
            'steel_stack_factor': 0.95,
            'steel_max_flux_density': 1.9,
            'l_stack': 16.0,
            'l_lamination': 0.2,
        }
        for k, v in defaults.items():
            if k not in self.material_dict:
                self.material_dict[k] = v

        self.material_dict['lamination_count'] = self.material_dict['l_stack'] * self.material_dict['steel_stack_factor'] / self.material_dict['l_lamination']
        if builtins.verbose:
            print('- 设计目标：我们的目标是二维场电机的电磁设计与分析')
            print(f"- 永磁体{self.material_dict['magnet_grade']}的厚度d_magnet={self.material_dict['d_magnet']}，退磁目标")
            print(f"- 硅钢片{self.material_dict['stator_steel']}叠长，影响电机的aspect ratio，铁芯由：lamination_count={self.material_dict['lamination_count']}组成")
