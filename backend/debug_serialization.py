import sys
import json
import traceback
from dataclasses import asdict, is_dataclass

sys.path.append(r'c:\Users\lenovo\Codes\ACMOP\backend\codes4')
import user_minitureMachine

def _specs_to_api_response(specs):
    g = specs.geometry
    w = specs.winding
    
    test_p = g.d_stator_outer
    print(f"DEBUG: type={type(test_p)}, is_dataclass={is_dataclass(test_p)}")
    
    # Mimic machine_specs.py logic
    stack_opts = getattr(g, "stack_length_options", [6.0, 10.0, 16.0])
    stack_length = float(stack_opts[len(stack_opts) // 2] if stack_opts else 10.0)
    
    print("Attempting asdict(g)...")
    geom = asdict(g)
    geom["stack_length"] = stack_length
    
    print("Attempting asdict(w)...")
    wind = asdict(w)
    
    print("Attempting asdict(targets)...")
    targets = asdict(specs.targets)
    
    print("Attempting asdict(materials)...")
    materials = asdict(specs.materials)
    
    response = {
        "geometry": geom,
        "winding": wind,
        "materials": materials,
        "targets": targets,
        "validations": {
            "geometry": specs.validate_inputs("geometry"),
            "materials": specs.validate_inputs("materials"),
            "winding": specs.validate_inputs("winding"),
            "targets": specs.validate_inputs("targets"),
        }
    }
    return response

try:
    specs = user_minitureMachine.MotorSpecs()
    res = _specs_to_api_response(specs)
    print("Transform successful. Attempting JSON dump...")
    json_str = json.dumps(res)
    print("JSON dump successful.")
except Exception as e:
    print("\nFAILURE DETECTED:")
    traceback.print_exc()
