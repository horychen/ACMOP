from fastapi import APIRouter
from machine_designer_v2 import AllPoints, Machine, RotorCore, StatorCore, Magnet, Coil
from ReactDrawer import ReactDrawer

router = APIRouter()
print("--- DEBUG ROUTER LOADED ---")

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
    machine.add_part(RotorCore(name="rotorCore", options="inner_notched"))
    machine.add_part(StatorCore(name="statorCore", options="closed-slot"))
    
    drawer = ReactDrawer()
    regions = drawer.draw_machine(machine)
    
    return {
        "regions": regions,
        "parameters": GP_as_dict
    }
