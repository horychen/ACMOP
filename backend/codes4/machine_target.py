from dataclasses import dataclass, field
from typing import Dict, Any, Optional
import os
import math
import numpy as np
from machine_geometry import RotorCore, StatorCore, Magnet

@dataclass
class MachineTarget:
    # Fixed Parameters
    fixed_parameters: Dict[str, Any] = field(default_factory=lambda: {
        'r_shaft': 0.0,
        'd_air_gap': 0.15,
        'r_stator_outer': 6.5,
        'd_stator_yoke': 0.3,
        'd_stator_tooth': 2.0,
        'w_stator_width': 1.2,
        'num_slots': 12,
        'num_poles': 10
    })

    # Free Parameters (Search Space)
    free_parameters: Dict[str, Any] = field(default_factory=lambda: {
        'r_rotor_outer': 4.0,
        'd_magnet': 3.0,
    })

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
    # project_name: str = field(default="", init=False)
    results_for_optimization: tuple = field(default=(), init=False)
    counter: int = field(default=0, init=False)

    def update_free_parameters(self, x: list):
        """Update free parameters from an optimization vector x."""
        # Mapping needs to be consistent
        keys = sorted(self.free_parameters.keys())
        for i, key in enumerate(keys):
            self.free_parameters[key] = x[i]

    def get_required_GP(self):
        """Combine fixed and free parameters into required_GP for AllPoints."""
        gp = self.fixed_parameters.copy()
        gp.update(self.free_parameters)
        
        # Calculate derived values if needed for AllPoints initialization
        # The user said "We no longer need derived parameters for drawing. 
        # We use updated free parameters and fixed parameters to update AllPoints object"
        # However, d_stator_tooth_shoe was derived in machine.py:
        # 'd_stator_tooth_shoe': 13/2 - 8/2 - 0.15 - 0.3 - 2
        
        r_so = gp['r_stator_outer']
        r_ro = gp['r_rotor_outer']
        g = gp['d_air_gap']
        dy = gp['d_stator_yoke']
        dt = gp['d_stator_tooth']
        
        # Explicitly calculate d_stator_tooth_shoe to satisfy AllPoints
        gp['d_stator_tooth_shoe'] = r_so - r_ro - g - dy - dt
        
        return gp




if __name__ == "__main__":
    target = MachineTarget()
    print("MachineTarget Debug:")
    print(f"Fixed: {target.fixed_parameters}")
    print(f"Free: {target.free_parameters}")
    print(f"Combined GP: {target.get_required_GP()}")
