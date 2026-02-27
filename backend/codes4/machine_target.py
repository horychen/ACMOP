from typing import Dict, Any

class MachineTarget:
    def __init__(self, target_dict: dict = None):
        self.target_dict = target_dict if target_dict is not None else {}
        
        if 'fea_config_dict' not in self.target_dict:
            self.target_dict['fea_config_dict'] = {
                'pc_name': 'localhost',
                'designer.show': True,
                'designer.max_nonlinear_iteration': 50,
                'mesh.average_size': 0.002, # 2mm
                'delete_results_after_calculation': False,
                'designer.JMAG_Scheduler': False,
                'designer.MultipleCPUs': True,
                'designer.AddIronLossCondition': True,
                'designer.OnlyTableResults': True,
                'designer.number_cycles_in_1stTSS': 0.5,
                'designer.number_cycles_in_2ndTSS': 0.5,
                'designer.number_cycles_in_3rdTSS': 0.0,
                'designer.number_cycles_prolonged': 0.0,
                'designer.number_of_steps_1stTSS': 20,
                'designer.number_of_steps_2ndTSS': 20,
                'designer.StepPerCycle_3rdTSS': 40,
                'designer.TranRef-StepPerCycle': 40,
                'designer.CircumferentialDivision': 720,
                'designer.meshSize_Magnet': 2.0,
                'designer.meshSize_Shaft': 2.0,
                'designer.meshSizeAir': 2.0,
                'designer.meshSize_General': 2.0,
            }
        if 'select_FEA_tool' not in self.target_dict:
            self.target_dict['select_FEA_tool'] = "JMAG"
        if 'machine_class' not in self.target_dict:
            self.target_dict['machine_class'] = "PMSM"

    def update_free_parameters(self, user_input, x: list):
        """Update geometric parameters in user_input from an optimization vector x."""
        keys = ['r_rotor_outer', 'd_magnet'] # Example search space
        for i, key in enumerate(keys):
            user_input['geometry'][key] = x[i]

    def get_required_GP(self, user_input):
        """Retrieve geometry dict from user_input."""
        return user_input['geometry']

if __name__ == "__main__":
    target = MachineTarget({})
    print("MachineTarget simplified.")
