from machine_designer_v2 import Machine, RotorCore, StatorCore, Magnet, Coil, AllPoints
from modern_machine_designer_utility import CairoDrawer
import JMAG
import numpy as np

def run_step_by_step():
    num_slots = 12
    num_poles = 10
    
    # Define Geometric Parameters (required_GP)
    GP_as_dict = {
        'r_shaft': 0,
        'r_rotor_outer': 8/2,
        'd_magnet': 3,
        'd_air_gap': 0.15,
        'r_stator_outer': 13/2,
        'd_stator_yoke': 0.3,
        'd_stator_tooth': 2,
        'd_stator_tooth_shoe': 13/2 - 8/2 - 0.15 - 0.3 - 2,
        'w_stator_width': 1.2,
        'num_slots': 12,
        'num_poles': 10
    }

    points = AllPoints(required_GP=GP_as_dict)

    print("--- Step 1: Initialize Machine with AllPoints ---")
    my_machine = Machine(all_points=points)
    print(f"Target: {num_slots} slots, {num_poles} poles.\n")

    print("--- Step 2: Adding Rotor Core (cylinder) ---")
    my_machine.add_part(RotorCore(name="rotorCore", options="cylinder"))

    print(f"--- Step 3: Adding Magnets ---")
    my_machine.add_part(Magnet(name="magnet", options="arc"))

    print("--- Step 4: Adding Stator Core (closed-slot) ---")
    my_machine.add_part(StatorCore(name="statorCore", options="closed-slot"))

    # print(f"--- Step 5: Adding Coils ---")
    # my_machine.add_part(Coil(name="coil", options="standard"))

    print("\n--- Step 6: Synchronizing Performance Metrics ---")
    my_machine.sync()
    print(f"Estimated Back-EMF: {my_machine.winding.estimated_back_emf:.2f} V")
    print(f"Phase Resistance: {my_machine.winding.phase_resistance:.3f} Ohm")
    print(f"Slot Current density: {my_machine.winding.rated_current_density_Js} A/mm2")

    print("\n--- Step 7: Generating SVG and Verifying Regions ---")
    drawer = CairoDrawer(filename='machine_geometry.svg', scale=30.0, bFillRegion=False)
    def draw_machine_using_CairoDrawer(machine):
        from machine_designer_v2 import draw_instruction_parser
        if machine.all_points is None:
            raise ValueError("Machine must have all_points initialized for drawing.")
        
        drawer.regions = []
        for part in machine.parts:

            # if isinstance(part, Magnet):
            #     part.iRotateCopy = 1

            # draw_instruction_parser internally calls drawer.drawArc / drawLine
            # which now return segment metadata instead of drawing immediately.
            region_dict = draw_instruction_parser(part, machine.all_points, drawer)
            drawer.regions.append(region_dict)

            # Now explicitly render and fill the region
            drawer.prepareSection(region_dict, color=part.color)
        
        drawer.surface.finish()
        print(f"Machine geometry drawn to {drawer.filename or 'SVG surface'}")
    draw_machine_using_CairoDrawer(my_machine)

    # Verify the first few regions
    if drawer.regions:
        print(f"Total Regions Collected: {len(drawer.regions)}")
        # Part 0 is usually RotorCore
        print(f"Rotor innerCoord: {drawer.regions[0]['innerCoord']}")




    print('--- Step 8: Drawing using JMAG Designer ---')
    def draw_machine_using_JMAG(machine):
        from machine_designer_v2 import draw_instruction_parser, RotorCore, StatorCore, Magnet, Coil
        
        # JMAG project configuration
        project_file = r'./jmag_machine_test.jproj'
        fea_config_dict = {
            'pc_name': 'localhost',
            'designer.Show': True
        }
        
        toolJd = JMAG.JMAG(fea_config_dict=fea_config_dict)
        # Open JMAG with dummy material and path info
        # Note: In a real environment, Steel_name and paths would come from machine specs.
        toolJd.open(Steel_name="M-19 Steel", 
                    expected_project_file_path=project_file, 
                    pc_name='localhost', 
                    dir_parent='../', 
                    bool_jmagDesignerShow=True)

        print(f"--- Drawing machine to JMAG ---")
        
        # Default gray colors from legacy code
        color_rgb_iron = np.array([236,236,236])/255
        color_rgb_magnet = np.array([226,226,226])/255

        # Get machine parameters for JMAG configuration
        for part in machine.parts:
            # Resolve geometry using unified parser
            region_dict = draw_instruction_parser(part, machine.all_points, toolJd)
            
            # Use metadata from the parser for JMAG optimization
            toolJd.bMirror = region_dict.get('bMirror', False)
            toolJd.iRotateCopy = region_dict.get('iRotateCopy', 0)
            
            # Configuration based on part type
            bRotateMerge = True
            color = None

            if isinstance(part, (RotorCore, StatorCore)):
                color = color_rgb_iron
            elif isinstance(part, Magnet):
                bRotateMerge = False
                color = color_rgb_magnet

            # Create regions in JMAG
            toolJd.prepareSection(region_dict, bRotateMerge=bRotateMerge, color=color)

        toolJd.save("MachineModel", "Generated using optimized machine.parts architecture")
        print(f"JMAG geometry drawn and saved to {project_file}")

    try:
        draw_machine_using_JMAG(my_machine) # Uncomment this if you want to run JMAG
    except Exception as e:
        print(f"JMAG Drawing failed (Expected if JMAG is not installed): {e}")




if __name__ == "__main__":
    run_step_by_step()
