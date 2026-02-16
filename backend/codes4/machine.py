from machine_geometry import AllPoints, RotorCore, StatorCore, Magnet, Coil, MachineGeometry
from machine_materials import MachineMaterial
from machine_winding import MachineWinding
from machine_target import MachineTarget
from modern_machine_designer_utility import CairoDrawer
import JMAG
import os, logging, jsonpickle, numpy as np, builtins
from time import time as clock_time
from collections import OrderedDict
from rich import print

# control printing out machine design pipeline
builtins.verbose = True
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Data Management
DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '_default', 'main_test'))
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)

class Machine:
    def __init__(self, user_input: dict):
        self.user_input = user_input
        self.materials = MachineMaterial()
        self.geometry = MachineGeometry()
        self.winding = MachineWinding(user_input['winding'])
        self.target = MachineTarget()

        if builtins.verbose:
            print('四个部件之间是相互依赖的，需要运行sync来同步，处理依赖关系。')

    def get_rotor_volume(self, stack_length=None):
        import math
        if stack_length is None:
            stack_length = self.winding.EX.get('mm_stack_length_specified', 0)
        return math.pi * (self.geometry.r_rotor_outer.value * 1e-3)**2 * (stack_length * 1e-3)

    def get_rotor_weight(self, gravity=9.8, stack_length=None):
        material_density_rho = 7860 # kg/m^3
        return gravity * self.get_rotor_volume(stack_length=stack_length) * material_density_rho

    def get_free_variables_as_dict(self):
        return self.target.free_parameters

    def sync(self):
        """Update AllPoints and sub-components based on user_input."""
        # Update derived parameters in user_input
        GP = self.user_input['geometry']
        GP.update({
            # derived geometric parameters
            'd_yoke': GP['r_stator_outer'] - GP['r_rotor_outer'] - GP['d_tooth'],
            'd_tooth_shoe': None if GP['tooth_shape']=='open' else GP['d_tooth_shoe'],
            'alpha_stator_tooth_span': 360 / self.user_input['winding']['slot_count'] * 0.7 if GP['tooth_shape']=='semi-closed' else None,
            'split_ratio': (GP['r_rotor_outer']+GP['d_air_gap']) / GP['r_stator_outer'],
            'aspect_ratio': GP['r_stator_outer'] * 2 / self.user_input['winding']['l_stack']
        })

        # Sync geometry points
        self.geometry.all_points = AllPoints(user_input=self.user_input)
        
        # Update parts' reference to all_points
        for part in self.geometry.parts:
            part.all_points = self.geometry.all_points

        # Sync sub-components
        self.geometry.sync(self.user_input)
        self.winding.sync(
            machineGeometry=self.geometry,
            materials=self.materials,
            l_stack=self.geometry.l_stack
        )

    def draw_machine_using_CairoDrawer(self, drawer):
        from machine_geometry import draw_instruction_parser
        if self.geometry.all_points is None:
            raise ValueError("Machine must have all_points initialized for drawing.")
        
        drawer.regions = []
        for part in self.geometry.parts:
            region_dict = draw_instruction_parser(part, self.geometry.all_points, drawer)
            drawer.regions.append(region_dict)
            drawer.prepareSection(region_dict, color=part.color)
        
        drawer.surface.finish()
        print(f"Machine geometry drawn to {drawer.filename or 'SVG surface'}")
    
    def draw_machine_using_JMAG(self):
        # JMAG project configuration
        expected_project_file = os.path.join(DATA_DIR, 'jmag_machine_test.jproj')
        
        # 1. Initialize JMAG Project
        toolJd = self.open_jmag(expected_project_file, self.materials.stator_steel)

        # 2. Draw Machine Parts
        self.draw_jmag(toolJd)
        print(f"JMAG geometry drawn and saved to {expected_project_file}")

    def open_jmag(self, expected_project_file, Steel_name, bool_jmagDesignerShow=True):
        self.target.project_name = 'default-project-name'

        toolJd = JMAG.JMAG(fea_config_dict=self.target.fea_config_dict)
        self.target.fea_config_dict = toolJd.fea_config_dict
        toolJd.open(Steel_name=Steel_name, 
                    expected_project_file_path=expected_project_file, 
                    pc_name=self.target.fea_config_dict['pc_name'], 
                    dir_parent='../', 
                    bool_jmagDesignerShow=bool_jmagDesignerShow)
        return toolJd

    def draw_jmag(self, toolJd):
        from machine_geometry import draw_instruction_parser, RotorCore, StatorCore, Magnet, Coil
        color_rgb_iron = np.array([236,236,236])/255
        color_rgb_magnet = np.array([226,226,226])/255

        for part in self.geometry.parts:
            region_dict = draw_instruction_parser(part, self.geometry.all_points, toolJd)
            
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
        toolJd.save(self.target.project_name, "Machine design generated by ACMOP modular architecture")

    def FEA_evaluate(self, project_loc, bool_jmagDesignerShow: bool = True, x_denorm=None, counter=0, counter_loop=0):
        """Perform FEA evaluation (JMAG simulation) and compile results."""
        self.target.counter = counter
        if x_denorm is not None:
            self.target.update_free_parameters(x_denorm)
        
        self.sync() # Ensure all components are updated

        # Project and Path Setup
        self.target.project_name = f"machine-ind{counter}"
        self.target.project_name += f"-redo{counter_loop}" if counter_loop > 0 else ""
        jmag_temp_dir = os.path.join(project_loc, "jmag_temp")
        jmag_screenshots_dir = os.path.join(project_loc, "jmag_screenshots")
        # Ensure path2FEACsv is something like C:\_Codes\ACMOP\backend\_default\main_test\csv\0
        path2FEACsv = os.path.join(project_loc, "csv", f"{counter}")
        
        for d in [jmag_temp_dir, jmag_screenshots_dir, path2FEACsv]:
            if not os.path.exists(d): os.makedirs(d)
        
        expected_project_file = os.path.join(jmag_temp_dir, f"{self.target.project_name}.jproj")
        self.target.swarm_data_json_file_path = os.path.join(project_loc, "SwarmData.json")
        self.target.path2Data = project_loc
        self.target.path2SwarmData = project_loc

        if 'JMAG' in self.target.select_FEA_tool:
            study_name = "Transient"

            # 1. Initialize JMAG Project
            toolJd = self.open_jmag(expected_project_file, self.materials.stator_steel)

            # 2. Draw Machine Parts
            self.draw_jmag(toolJd)

            # 3. Solver Setup 
            app = toolJd.app
            model = app.GetModel(self.target.project_name)
            toolJd.pre_process_PMSM(app, model, self)
            study = toolJd.add_magnetic_transient_study(app, model, path2FEACsv, study_name, self)

            # # 4. Mesh and Run
            # toolJd.mesh_study(self, app, model, study, output_dir=project_loc)
            # toolJd.run_study(self, app, study, self.target.fea_config_dict, clock_time())

            # # 5. Data Export
            # if self.target.fea_config_dict['delete_results_after_calculation'] == False:
            #     ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
            #     app.GetDataManager().CreateGraphModel(ref1)
            #     app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(os.path.join(path2FEACsv, f"{study_name}_EXPORT_CIRCUIT_VOLTAGE.csv"))

            # # 6. Compile & Save Results
            # self.compile_results(toolJd, study_name, path2FEACsv, counter)

        else:
            raise Exception(f"FEA tool {self.target.select_FEA_tool} not implemented or supported.")

    def compile_results(self, toolJd, study_name, path2FEACsv, counter):
        """Pack results into dictionary and save to JSON."""
        results = toolJd.build_str_results(self, self.target.project_name, study_name, path2FEACsv, self.target.fea_config_dict, femm_solver=None)
        
        # Unpack results indices based on build_str_results contract
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
            'x_denorm_dict': self.target.free_parameters.copy(),
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
            'select_FEA_tool': self.target.select_FEA_tool,
        }

        # Save to JSON via jsonpickle
        results2file = {f'spec_performance_dict-ind{counter}': spec_performance_dict}
        
        try:
            if os.path.exists(self.target.swarm_data_json_file_path) and os.path.getsize(self.target.swarm_data_json_file_path) > 0:
                with open(self.target.swarm_data_json_file_path, 'r') as rf:
                    existing_data = jsonpickle.decode(rf.read())
            else:
                existing_data = {}
        except:
            existing_data = {}

        existing_data.update(results2file)
        with open(self.target.swarm_data_json_file_path, 'w') as wf:
            wf.write(jsonpickle.encode(existing_data, indent=4))

        self.target.results_for_optimization = (cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle)

def run_step_by_step():

    # 用户需要输入的信息
    user_input = OrderedDict({
        'winding':{
            'slot_count': 12,
            'pole_count': 10,
            'coil_pitch': 1,
            'l_stack': 16.0,
        },
        'geometry':{
            'tooth_shape': 'closed',
            'd_tooth_shoe': 0.3, # this makes effective tooth depth becomes: d_tooth - d_tooth_shoe
            'd_air_gap': 0.15,
            'r_stator_outer': 13/2,
            'r_rotor_outer': 8/2,
            'r_shaft': 0.0,
            'd_tooth': 2.0,
            'w_width': 1.2,
            'd_magnet': 1.0,
        },
        'target': [
            'search w_stator_width within [1.0, 1.5]',
            'search d_stator_tooth within [1.6, 2.0]',
            'search d_magnet within [1.0, 3.0]',
        ]
    })



    print("--- Step 0: Initialize Composite Machine ---")
    dex13 = Machine(user_input)

    print(dex13.materials)
    print(dex13.geometry)
    print(dex13.winding)
    print(dex13.target)

    # Initialize geometry points from user_input
    # dex13.sync_all_points()

    print("--- Step 2: Adding Rotor Core ---")
    dex13.geometry.add_part(RotorCore(name="rotorCore", options="cylinder"))

    print(f"--- Step 3: Adding Magnets ---")
    dex13.geometry.add_part(Magnet(name="magnet", options="arc"))

    print("--- Step 4: Adding Stator Core ---")
    dex13.geometry.add_part(StatorCore(name="statorCore", options="closed-slot"))

    print(f"--- Step 5: Adding Coils ---")
    dex13.geometry.add_part(Coil(name="coil", options="standard"))

    print("\n--- Step 6: Synchronizing Performance Metrics ---")
    dex13.sync()
    print(f"Estimated Back-EMF: {dex13.winding.estimated_back_emf:.2f} V")
    print(f"Phase Resistance: {dex13.winding.phase_resistance:.3f} Ohm")
    print(f"Slot Current density: {dex13.winding.rated_current_density} A/mm2")

    print("\n--- Step 7: Generating SVG and Verifying Regions ---")
    svg_path = os.path.join(DATA_DIR, 'machine_geometry.svg')
    drawer = CairoDrawer(filename=svg_path, scale=30.0, bFillRegion=False)
    dex13.draw_machine_using_CairoDrawer(drawer)

    # Verify the first few regions
    if drawer.regions:
        print(f"Total Regions Collected: {len(drawer.regions)}")
        print(f"Rotor innerCoord: {drawer.regions[0]['innerCoord']}")

    if False:
        print('--- Step 8: Drawing using JMAG Designer ---')
        dex13.draw_machine_using_JMAG() 
    else:
        print('--- Step 9: Evaluate FEA simulaion in JMAG Designer ---')
        dex13.FEA_evaluate(DATA_DIR, bool_jmagDesignerShow=True)

if __name__ == "__main__":

    run_step_by_step()
