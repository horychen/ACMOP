from dataclasses import dataclass, field
from typing import Dict, Any, Optional
import os
import math
import numpy as np
from machine_geometry import RotorCore, StatorCore, Magnet

@dataclass
class MachineTarget:
    # FEA Configuration
    select_FEA_tool: str = "JMAG"
    machine_class: str = "PMSM"
    fea_config_dict: Dict[str, Any] = field(default_factory=lambda: {
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
    })
    swarm_data_json_file_path: str = "SwarmData.json"
    path2Data: str = ""
    path2SwarmData: str = ""
    results_for_optimization: tuple = field(default=(), init=False)
    counter: int = field(default=0, init=False)

    def update_free_parameters(self, user_input, x: list):
        """Update geometric parameters in user_input from an optimization vector x."""
        # This now targets user_input directly
        keys = ['r_rotor_outer', 'd_magnet'] # Example search space
        for i, key in enumerate(keys):
            user_input['geometry'][key] = x[i]

    def get_required_GP(self, user_input):
        """Retrieve geometry dict from user_input."""
        return user_input['geometry']

if __name__ == "__main__":
    target = MachineTarget()
    print("MachineTarget simplified.")




if __name__ == "__main__":
    target = MachineTarget()
    print("MachineTarget Debug:")
    print(f"Fixed: {target.fixed_parameters}")
    print(f"Free: {target.free_parameters}")
    print(f"Combined GP: {target.get_required_GP()}")
