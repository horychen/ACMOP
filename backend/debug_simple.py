import sys
import os

# Add codes4 to path
sys.path.append(os.path.join(os.getcwd(), 'backend', 'codes4'))

from user_minitureMachine import MotorSpecs
import math

try:
    print("Initializing MotorSpecs...")
    specs = MotorSpecs()
    print("Initialization success.")
    
    print("\nValidating 'geometry' step...")
    res = specs.validate_inputs("geometry")
    print("Geometry validation status:", res['status'])
    
    print("\nValidating 'winding' step...")
    res = specs.validate_inputs("winding")
    print("Winding validation status:", res['status'])
    
    print("\nValidating 'materials' step...")
    res = specs.validate_inputs("materials")
    print("Materials validation status:", res['status'])
    
    print("\nValidation process completed without crash.")
except Exception as e:
    import traceback
    print("\nCRASH DETECTED:")
    traceback.print_exc()
