from app.routers.machine_specs import MotorSpecs
import json

try:
    specs = MotorSpecs()
    print("MotorSpecs initialized.")
    
    # Test validate_inputs
    print("\nTesting validate_inputs('geometry')...")
    res = specs.validate_inputs("geometry")
    print("Validation success:", res['status'])
    
    # Test serialization
    print("\nTesting serialization...")
    from app.routers.machine_specs import _robust_dict
    serialized = _robust_dict(specs)
    print("Serialization success.")
    
except Exception as e:
    import traceback
    print("\nFATAL ERROR:")
    traceback.print_exc()
