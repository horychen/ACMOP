from fastapi import APIRouter
from machine_designer_v2 import AllPoints

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
    
    # Return HP, HP_mirror and RP dictionaries
    # Keys in HP/RP are ints, but FastAPI/JSON will convert them to strings
    return {
        "HP": points.HP,
        "HP_mirror": getattr(points, 'HP_mirror', {}),
        "RP": points.RP,
        "parameters": GP_as_dict
    }
