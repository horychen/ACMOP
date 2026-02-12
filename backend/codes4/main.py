from machine_geometry import AllPoints, RotorCore, StatorCore, Magnet, Coil, MachineGeometry
from machine_materials import MachineMaterial
from machine_winding import MachineWinding
from machine_target import MachineTarget
from modern_machine_designer_utility import CairoDrawer
import JMAG
import numpy as np
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Data Management
DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '_default', 'main_test'))
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)

class Machine:
    def __init__(self):
        self.materials = MachineMaterial()
        self.target = MachineTarget()
        self.winding = MachineWinding()
        self.geometry = MachineGeometry()

    @property
    def machineGeometry(self):
        return self.geometry.machineGeometry

    @property
    def p(self):
        from types import SimpleNamespace
        return SimpleNamespace(value=self.target.fixed_parameters['num_poles'] // 2)

    @property
    def Qs(self):
        from types import SimpleNamespace
        return SimpleNamespace(value=self.target.fixed_parameters['num_slots'])
        
    def sync_all_points(self):
        """Update AllPoints based on current target parameters."""
        gp = self.target.get_required_GP()
        self.geometry.gp = gp
        self.geometry.all_points = AllPoints(required_GP=gp)
        # Update winding/poles if they changed in target
        self.winding.slot_count = self.target.fixed_parameters['num_slots']
        self.winding.pole_count = self.target.fixed_parameters['num_poles']
        self.geometry.num_slots = self.winding.slot_count
        self.geometry.num_poles = self.winding.pole_count

    def sync(self):
        """Synchronize performance metrics."""
        self.sync_all_points()
        self.winding.sync(
            geometry_points=self.geometry.all_points.HP,
            materials=self.materials,
            l_stack=self.geometry.l_stack
        )

def run_step_by_step():
    print("--- Step 0: Initialize Composite Machine ---")
    my_machine = Machine()
    
    # Initialize geometry points from default target GP
    my_machine.sync_all_points()

    print("--- Step 1: Geometry Setup ---")
    print(f"Target: {my_machine.winding.slot_count} slots, {my_machine.winding.pole_count} poles.\n")

    print("--- Step 2: Adding Rotor Core ---")
    my_machine.geometry.add_part(RotorCore(name="rotorCore", options="cylinder"))

    print(f"--- Step 3: Adding Magnets ---")
    my_machine.geometry.add_part(Magnet(name="magnet", options="arc"))

    print("--- Step 4: Adding Stator Core ---")
    my_machine.geometry.add_part(StatorCore(name="statorCore", options="closed-slot"))

    print(f"--- Step 5: Adding Coils ---")
    my_machine.geometry.add_part(Coil(name="coil", options="standard"))

    print("\n--- Step 6: Synchronizing Performance Metrics ---")
    my_machine.sync()
    print(f"Estimated Back-EMF: {my_machine.winding.estimated_back_emf:.2f} V")
    print(f"Phase Resistance: {my_machine.winding.phase_resistance:.3f} Ohm")
    print(f"Slot Current density: {my_machine.winding.rated_current_density} A/mm2")

    print("\n--- Step 7: Generating SVG and Verifying Regions ---")
    svg_path = os.path.join(DATA_DIR, 'machine_geometry.svg')
    drawer = CairoDrawer(filename=svg_path, scale=30.0, bFillRegion=False)
    
    def draw_machine_using_CairoDrawer(machine):
        from machine_geometry import draw_instruction_parser
        if machine.geometry.all_points is None:
            raise ValueError("Machine must have all_points initialized for drawing.")
        
        drawer.regions = []
        for part in machine.geometry.parts:
            region_dict = draw_instruction_parser(part, machine.geometry.all_points, drawer)
            drawer.regions.append(region_dict)
            drawer.prepareSection(region_dict, color=part.color)
        
        drawer.surface.finish()
        print(f"Machine geometry drawn to {drawer.filename or 'SVG surface'}")
    
    draw_machine_using_CairoDrawer(my_machine)

    # Verify the first few regions
    if drawer.regions:
        print(f"Total Regions Collected: {len(drawer.regions)}")
        print(f"Rotor innerCoord: {drawer.regions[0]['innerCoord']}")

    print('--- Step 8: Drawing using JMAG Designer ---')
    def draw_machine_using_JMAG(machine):
        from machine_geometry import draw_instruction_parser, RotorCore, StatorCore, Magnet, Coil
        
        # JMAG project configuration
        project_file = os.path.join(DATA_DIR, 'jmag_machine_test.jproj')
        fea_config_dict = {
            'pc_name': 'localhost',
            'designer.show': True,
            'designer.max_nonlinear_iteration': 50,
            'mesh.average_size': 0.002, # 2mm
            'delete_results_after_calculation': False,
            'designer.JMAG_Scheduler': False,
            'designer.MultipleCPUs': True,
        }
        
        toolJd = JMAG.JMAG(fea_config_dict=fea_config_dict)
        toolJd.open(Steel_name=machine.materials.stator_steel, 
                    expected_project_file_path=project_file, 
                    pc_name='localhost', 
                    dir_parent='../', 
                    bool_jmagDesignerShow=True)

        print(f"--- Drawing machine to JMAG ---")
        
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

            toolJd.prepareSection(region_dict, bRotateMerge=bRotateMerge, color=color)

        toolJd.save("MachineModel", "Generated using optimized machine.parts architecture")
        print(f"JMAG geometry drawn and saved to {project_file}")

    try:
        draw_machine_using_JMAG(my_machine) 
    except Exception as e:
        print(f"JMAG Drawing failed: {e}")

    print('--- Step 9: Evaluate FEA simulaion in JMAG Designer ---')
    my_machine.target.FEA_evaluate(my_machine, DATA_DIR, bool_jmagDesignerShow=True)

if __name__ == "__main__":
    run_step_by_step()
