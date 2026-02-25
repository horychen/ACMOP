from fastapi import APIRouter
from machine_geometry import AllPoints, MachineGeometry as Machine, RotorCore, StatorCore, Magnet, Coil
from ReactDrawer import ReactDrawer
import machine
from collections import OrderedDict

router = APIRouter()

def get_default_user_input():
    return OrderedDict({
        'winding':{
            'phase_count': 3,
            'slot_count': 12,
            'pole_count': 10,
            'coil_pitch': 1,
            'l_stack': 16.0,
            'rated_current_density': 14,
            'dc_bus_voltage': 12,
            'fill_factor': 0.58,
            'wire_diameter_with_insulation': 0.226,
            'wire_diameter': 0.179,
            'connection_type': 'wye',
            'rated_speed': 20000,
            'torque_current_ratio': 1.0,
            'suspension_current_ratio': 0.0,
        },
        'geometry':{
            'tooth_shape': 'closed',
            'd_tooth_shoe': 0.3,
            'd_air_gap': 0.15,
            'r_stator_outer': 13/2,
            'r_rotor_outer': 8/2,
            'r_shaft': 0.0,
            'd_tooth': 2.0,
            'w_tooth': 1.2,
            'd_magnet': 3.0,
            'd_yoke': 13/2 - 8/2 - 2.0,
        },
        'target': []
    })

@router.get("/points")
async def get_debug_points():
    user_input = get_default_user_input()
    points = AllPoints(user_input=user_input)
    
    # Calculate RP_mirror dynamically like parse_point_name does
    RP_mirror = {idx: (p[0], -p[1]) for idx, p in points.RP.items()}

    # Return HP, HP_mirror, RP and RP_mirror dictionaries
    return {
        "HP": points.HP,
        "HP_mirror": getattr(points, 'HP_mirror', {}),
        "RP": points.RP,
        "RP_mirror": RP_mirror,
        "parameters": user_input["geometry"],
        "slot_count": user_input["winding"]["slot_count"],
        "pole_count": user_input["winding"]["pole_count"],
        "version": "geometry-debugger-v1"
    }

@router.get("/geometry")
async def get_debug_geometry():
    user_input = get_default_user_input()
    points = AllPoints(user_input=user_input)
    machine_geom = Machine(all_points=points)
    machine_geom.winding = user_input['winding']
    
    # Add parts matching the current main.py flow
    machine_geom.add_part(RotorCore(name="rotorCore", options="notched"))
    machine_geom.add_part(StatorCore(name="statorCore", options="closed-slot"))
    # Add Magnets
    machine_geom.add_part(Magnet(name="magnet", color="#FF0000"))
    # Add Coils
    machine_geom.add_part(Coil(name="coil", color="#0000FF"))
    
    drawer = ReactDrawer()
    regions = drawer.draw_machine(machine_geom)
    
    return {
        "regions": regions,
        "parameters": user_input["geometry"]
    }

@router.get("/inspection")
async def get_inspection_data():
    user_input = get_default_user_input()
    my_machine = machine.Machine(user_input)
    my_machine.sync() # Ensure all points and winding are updated
    
    return {
        "m_spec": {
            "fixed_parameters": {"geometry": user_input["geometry"], "winding": user_input["winding"]},
            "materials": {
                "stator_steel": my_machine.materials.stator_steel,
                "rotor_steel": my_machine.materials.rotor_steel,
                "magnet_grade": my_machine.materials.magnet_grade,
                "magnet_br": my_machine.materials.magnet_br,
                "copper_fill_factor": user_input['winding']['fill_factor']
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
        "design_parameters": {"target": user_input["target"]},
        "jmag_study": my_machine.target.fea_config_dict
    }
