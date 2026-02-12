from dataclasses import dataclass, field
from typing import Dict, Any, Optional
import os
import math
import numpy as np
import jsonpickle
import JMAG
from time import time as clock_time
from machine_geometry import draw_instruction_parser, RotorCore, StatorCore, Magnet
from modern_machine_designer_utility import Winding

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
    project_name: str = field(default="", init=False)
    results_for_optimization: tuple = field(default=(), init=False)

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
        # However, d_stator_tooth_shoe was derived in main.py:
        # 'd_stator_tooth_shoe': 13/2 - 8/2 - 0.15 - 0.3 - 2
        
        r_so = gp['r_stator_outer']
        r_ro = gp['r_rotor_outer']
        g = gp['d_air_gap']
        dy = gp['d_stator_yoke']
        dt = gp['d_stator_tooth']
        
        # Explicitly calculate d_stator_tooth_shoe to satisfy AllPoints
        gp['d_stator_tooth_shoe'] = r_so - r_ro - g - dy - dt
        
        return gp

    def FEA_evaluate(self, machine, project_loc, bool_jmagDesignerShow: bool = True, x_denorm=None, counter=0, counter_loop=0):
        """Perform FEA evaluation (JMAG simulation) and compile results."""
        if x_denorm is not None:
            self.update_free_parameters(x_denorm)
        
        machine.sync() # Ensure all components are updated

        # 1. Project and Path Setup
        self.project_name = f"machine-ind{counter}"
        self.project_name += f"-redo{counter_loop}" if counter_loop > 0 else ""

        jmag_temp_dir = os.path.join(project_loc, "jmag_temp")
        jmag_screenshots_dir = os.path.join(project_loc, "jmag_screenshots")
        # Ensure path2FEACsv is something like C:\_Codes\ACMOP\backend\_default\main_test\csv\0
        path2FEACsv = os.path.join(project_loc, "csv", f"{counter}")
        
        for d in [jmag_temp_dir, jmag_screenshots_dir, path2FEACsv]:
            if not os.path.exists(d): os.makedirs(d)
        
        expected_project_file = os.path.join(jmag_temp_dir, f"{self.project_name}.jproj")
        self.swarm_data_json_file_path = os.path.join(project_loc, "SwarmData.json")

        # 2. Compatibility for legacy JMAG.py: Populate machine attributes
        w = machine.winding
        m = machine.materials
        g = machine.geometry
        
        # Add pole pair and slot info for JMAG.py
        from types import SimpleNamespace
        machine.p = SimpleNamespace(value=w.pole_count // 2)
        machine.ps = SimpleNamespace(value=2) 
        machine.Qs = SimpleNamespace(value=w.slot_count)
        machine.s = SimpleNamespace(value=1) # Adjusted for 10-magnet configuration
        machine.bool_PermanentMagnet = True
        machine.machine_class = self.machine_class
        machine.path2SwarmData = project_loc
        machine.path2FEACsv = path2FEACsv
        machine.name = self.project_name
        machine.project_name = self.project_name
        
        # Geometry compatibility for JMAG.py (expects .value)
        gp = machine.geometry.all_points.required_GP
        machine.mm_r_ro = SimpleNamespace(value=gp['r_rotor_outer'])
        machine.mm_r_si = SimpleNamespace(value=gp['r_rotor_outer'] + gp['d_air_gap'])
        machine.mm_d_air_gap = SimpleNamespace(value=gp['d_air_gap'])

        # Additional JMAG compatibility attributes
        fp = machine.target.fixed_parameters
        
        def get_val(key, default=0):
            val = fp.get(key, default)
            if hasattr(val, 'value'): return val
            return SimpleNamespace(value=val)

        machine.mm_d_sleeve = get_val('mm_d_sleeve', 0)
        machine.mm_d_mech_air_gap = get_val('mm_d_mech_air_gap', 0)
        machine.mm_d_pm = get_val('mm_d_pm', 0)
        machine.deg_alpha_rs = get_val('deg_alpha_rs', 0)
        machine.deg_alpha_rm = get_val('deg_alpha_rm', 0)
        
        # ID dummies
        machine.id_shaft = 0
        machine.id_sleeve = 0
        
        # machineGeometry mapping for legacy JMAG.py compatibility
        machine.machineGeometry = {}
        for part in machine.geometry.parts:
            if 'rotor' in part.name.lower() and 'core' in part.name.lower():
                machine.machineGeometry['rotorCore'] = part
            elif 'stator' in part.name.lower() and 'core' in part.name.lower():
                machine.machineGeometry['statorCore'] = part

        # Populate EX dictionary
        machine.EX = {
            'mm_stack_length_specified': g.l_stack,
            'RatedSpeed': w.rated_speed,
            'ExcitationFreqSimulated': w.excitation_frequency,
            'Magnet_Name': m.magnet_grade,
            'StatorCore_Material': m.stator_steel,
            'RotorCore_Material': m.rotor_steel,
            'LaminationFactor': m.steel_stack_factor * 100.0,
            'Magnet_Temperature': 20.0, # Default
            'InitialRotationAngle': w.initial_rotation_angle,
            'DriveW_CurrentAmp': w.drive_winding_current,
            'BeariW_CurrentAmp': w.bearing_winding_current,
            'DriveW_Rs': w.drive_winding_resistance,
            'BeariW_Rs': w.bearing_winding_resistance,
            'DriveW_zQ': w.drive_winding_conductors_per_slot,
            'BeariW_zQ': w.bearing_winding_conductors_per_slot,
            'Magnet_StartAngle': 0.5 * 360 / (2 * machine.p.value),
            'SteelMaterial': m.stator_steel
        }
        print(f"DEBUG: w.excitation_frequency = {w.excitation_frequency}")
        print(f"DEBUG: machine.EX['ExcitationFreqSimulated'] = {machine.EX['ExcitationFreqSimulated']}")

        # Initialize wily (Winding object from utility)
        machine.wily = Winding(
            phase_number_m=w.phase_count,
            stator_slot_number_Qs=w.slot_count,
            pole_pair_number_p=machine.p.value,
            suspension_pole_pair_number_ps=machine.ps.value,
            coil_pitch_y=w.coil_pitch,
            bool_DPNVorSEPA=True, # Default
            number_of_parallel_branch=w.parallel_branch_count
        )
        # Placeholder for winding layout derivation if needed by JMAG
        machine.wily.bool_CustomizedCircuit = False
        machine.wily.bool_3PhaseCurrentSource = True
        machine.wily.bool_DPNVorSEPA = False # Adjusting based on circuit logic in JMAG.py

        if 'JMAG' in self.select_FEA_tool:
            study_name = "Transient"

            # 1. Initialize JMAG Project
            toolJd = JMAG.JMAG(fea_config_dict=self.fea_config_dict)
            machine.fea_config_dict = toolJd.fea_config_dict
            toolJd.open(Steel_name=machine.materials.stator_steel, 
                       expected_project_file_path=expected_project_file, 
                       pc_name=self.fea_config_dict['pc_name'], 
                       dir_parent='../', 
                       bool_jmagDesignerShow=bool_jmagDesignerShow)

            # 2. Draw Machine Parts
            color_rgb_iron = np.array([236,236,236])/255
            color_rgb_magnet = np.array([226,226,226])/255

            for part in machine.geometry.parts:
                region_dict = draw_instruction_parser(part, machine.geometry.all_points, toolJd)
                
                toolJd.bMirror = region_dict.get('bMirror', False)
                toolJd.iRotateCopy = region_dict.get('iRotateCopy', 0)
                
                bRotateMerge = True
                color = None

                if isinstance(part, (RotorCore, StatorCore)):
                    color = color_rgb_iron
                elif isinstance(part, Magnet):
                    bRotateMerge = False
                    color = color_rgb_magnet

                part.inner_coords = region_dict['inner_coords']
                toolJd.prepareSection(region_dict, bRotateMerge=bRotateMerge, color=color)

            # Import Model into Designer
            toolJd.save(self.project_name, "Machine design generated by ACMOP modular architecture")

            # 3. Solver Setup & Run
            app = toolJd.app
            model = app.GetModel(self.project_name)

            if 'PMSM' in self.machine_class:
                toolJd.pre_process_PMSM(app, model, machine)

            study = toolJd.add_magnetic_transient_study(app, model, path2FEACsv, study_name, machine)
            toolJd.mesh_study(machine, app, model, study, output_dir=project_loc)
            toolJd.run_study(machine, app, study, self.fea_config_dict, clock_time())

            # 4. Data Export
            if self.fea_config_dict['delete_results_after_calculation'] == False:
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(os.path.join(path2FEACsv, f"{study_name}_EXPORT_CIRCUIT_VOLTAGE.csv"))

            # 5. Compile & Save Results
            self.compile_results(machine, toolJd, study_name, path2FEACsv, counter)

        else:
            raise Exception(f"FEA tool {self.select_FEA_tool} not implemented or supported.")

    def compile_results(self, machine, toolJd, study_name, path2FEACsv, counter):
        """Pack results into dictionary and save to JSON."""
        results = toolJd.build_str_results(machine, self.project_name, study_name, path2FEACsv, self.fea_config_dict, femm_solver=None)
        
        # Unpack results indices based on build_str_results contract
        # (Simplified mapping for this refactoring, assuming build_str_results returns a list compatible with legacy expectations)
        (cost_function, f1, f2, f3, FRW, \
         normalized_torque_ripple, \
         normalized_force_error_magnitude, \
         force_error_angle, \
         project_name, individual_name, \
         number_current_generation, individual_index,\
         power_factor, \
         rated_ratio, \
         rated_stack_length_mm, \
         rated_total_loss, \
         rated_stator_copper_loss_along_stack, \
         rated_magnet_Joule_loss, \
         rated_rotor_copper_loss_along_stack, \
         stator_copper_loss_in_end_turn, \
         rotor_copper_loss_in_end_turn, \
         rated_iron_loss, \
         rated_windage_loss, \
         str_results, \
         mm2_slot_area, \
         coil_flux_linkage_peak2peak_value, \
         TRV, Cost, Cost_Fe, Cost_Cu, Cost_PM, \
         ss_avg_force_magnitude, rotor_weight, torque_average) = results

        spec_performance_dict = {
            'x_denorm_dict': self.free_parameters.copy(),
            'project_name': project_name,
            'individual_name': individual_name,
            'f1': float(f1),
            'f2': float(f2),
            'f3': float(f3),
            'TRV': float(TRV),
            'FRW': float(FRW),
            'torque_average': torque_average,
            'ss_avg_force_magnitude': ss_avg_force_magnitude,
            'rotor_weight': float(rotor_weight),
            'normalized_torque_ripple': float(normalized_torque_ripple),
            'normalized_force_error_magnitude': float(normalized_force_error_magnitude),
            'force_error_angle': float(force_error_angle),
            'mm2_slot_area': mm2_slot_area,
            'Cost': float(Cost),
            'rated_total_loss': float(rated_total_loss),
            'select_FEA_tool': self.select_FEA_tool,
        }

        # Save to JSON via jsonpickle
        results2file = {f'spec_performance_dict-ind{counter}': spec_performance_dict}
        
        try:
            if os.path.exists(self.swarm_data_json_file_path) and os.path.getsize(self.swarm_data_json_file_path) > 0:
                with open(self.swarm_data_json_file_path, 'r') as rf:
                    existing_data = jsonpickle.decode(rf.read())
            else:
                existing_data = {}
        except:
            existing_data = {}

        existing_data.update(results2file)
        with open(self.swarm_data_json_file_path, 'w') as wf:
            wf.write(jsonpickle.encode(existing_data, indent=4))

        self.results_for_optimization = (cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle)

if __name__ == "__main__":
    target = MachineTarget()
    print("MachineTarget Debug:")
    print(f"Fixed: {target.fixed_parameters}")
    print(f"Free: {target.free_parameters}")
    print(f"Combined GP: {target.get_required_GP()}")
