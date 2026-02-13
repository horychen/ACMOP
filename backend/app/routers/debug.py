from fastapi import APIRouter
from machine_geometry import AllPoints, MachineGeometry as Machine, RotorCore, StatorCore, Magnet, Coil
from ReactDrawer import ReactDrawer
import main

router = APIRouter()

@router.get("/points")
async def get_debug_points():
    # Default Geometric Parameters (required_GP) matching main.py example
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
    
    # Calculate RP_mirror dynamically like parse_point_name does
    RP_mirror = {idx: (p[0], -p[1]) for idx, p in points.RP.items()}

    # Return HP, HP_mirror, RP and RP_mirror dictionaries
    return {
        "HP": points.HP,
        "HP_mirror": getattr(points, 'HP_mirror', {}),
        "RP": points.RP,
        "RP_mirror": RP_mirror,
        "parameters": GP_as_dict,
        "num_slots": GP_as_dict['num_slots'],
        "num_poles": GP_as_dict['num_poles'],
        "version": "geometry-debugger-v1"
    }

@router.get("/geometry")
async def get_debug_geometry():
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
    machine = Machine(all_points=points)
    
    # Add parts matching the current main.py flow
    machine.add_part(RotorCore(name="rotorCore", options="notched"))
    machine.add_part(StatorCore(name="statorCore", options="closed-slot"))
    # Add Magnets
    machine.add_part(Magnet(name="magnet", color="#FF0000"))
    # Add Coils
    machine.add_part(Coil(name="coil", color="#0000FF"))
    
    drawer = ReactDrawer()
    regions = drawer.draw_machine(machine)
    
    return {
        "regions": regions,
        "parameters": GP_as_dict
    }

@router.get("/inspection")
async def get_inspection_data():
    my_machine = main.Machine()
    my_machine.sync() # Ensure all points and winding are updated
    
    return {
        "m_spec": {
            "fixed_parameters": my_machine.target.fixed_parameters,
            "materials": {
                "stator_steel": my_machine.materials.stator_steel,
                "rotor_steel": my_machine.materials.rotor_steel,
                "magnet_grade": my_machine.materials.magnet_grade,
                "magnet_br": my_machine.materials.magnet_br,
                "copper_fill_factor": my_machine.materials.copper_fill_factor
            }
        },
        "m_para": {
            "winding_parameters": {
                "slot_count": my_machine.winding.slot_count,
                "pole_count": my_machine.winding.pole_count,
                "estimated_back_emf": my_machine.winding.estimated_back_emf,
                "phase_resistance": my_machine.winding.phase_resistance,
                "rated_current_density": my_machine.winding.rated_current_density,
                "l_stack": my_machine.geometry.l_stack
            },
            "other_derived": {
                "rotor_volume": my_machine.get_rotor_volume(),
                "rotor_weight": my_machine.get_rotor_weight()
            }
        },
        "design_parameters": my_machine.target.free_parameters,
        "jmag_study": my_machine.target.fea_config_dict
    }
