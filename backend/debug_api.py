import sys
import traceback
sys.path.append(r'c:\Users\lenovo\Codes\ACMOP\backend\codes4')
import user_minitureMachine

try:
    specs = user_minitureMachine.MotorSpecs()
    print("Testing materials...")
    specs.validate_inputs("materials")
    print("Testing geometry...")
    specs.validate_inputs("geometry")
    print("Testing winding...")
    specs.validate_inputs("winding")
    print("Testing targets...")
    specs.validate_inputs("targets")
    print("ALL OK")
except Exception as e:
    print("FAILURE")
    traceback.print_exc()
