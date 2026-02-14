from machine import Machine
from collections import OrderedDict

def test():
    user_input = OrderedDict({
        'winding':{
            'slot_count': 12,
            'pole_count': 10,
            'coil_pitch': 1,
            'l_stack': 16.0,
        },
        'geometry':{
            'tooth_shape': 'closed',
            'd_tooth_shoe': 0.3,
            'd_air_gap': 0.15,
            'r_stator_outer': 13/2,
            'r_rotor_outer': 8/2,
            'r_shaft': 0.0,
            'd_tooth': 2.0,
            'w_width': 1.2,
            'd_magnet': 1.0,
        }
    })
    
    print("Initializing Machine...")
    m = Machine(user_input)
    
    from machine_geometry import RotorCore, Magnet, StatorCore, Coil
    print("Adding parts...")
    m.geometry.add_part(RotorCore(name="rotorCore", options="cylinder"))
    m.geometry.add_part(Magnet(name="magnet", options="arc"))
    m.geometry.add_part(StatorCore(name="statorCore", options="closed-slot"))
    m.geometry.add_part(Coil(name="coil", options="standard"))

    print("Syncing Machine...")
    m.sync()
    print("Sync successful!")
    print(f"Back-EMF: {m.winding.estimated_back_emf:.2f}")

    from modern_machine_designer_utility import CairoDrawer
    import os
    DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '_default', 'main_test'))
    if not os.path.exists(DATA_DIR): os.makedirs(DATA_DIR)
    svg_path = os.path.join(DATA_DIR, 'verify_geometry.svg')
    print(f"Drawing to {svg_path}...")
    drawer = CairoDrawer(filename=svg_path, scale=30.0, bFillRegion=False)
    m.draw_machine_using_CairoDrawer(drawer)
    print("Drawing successful!")

if __name__ == "__main__":
    test()
