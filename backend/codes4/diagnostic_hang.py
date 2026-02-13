
import sys
import os
import time


try:
    print("Simulating router's dynamic load...")
    import types
    path = os.path.abspath("machine_geometry.py")
    codes4_dir = os.path.dirname(path)
    if codes4_dir not in sys.path:
        sys.path.insert(0, codes4_dir)
    
    module = types.ModuleType("machine_geometry")
    module.__file__ = path
    with open(path, "r", encoding="utf-8") as f:
        source = f.read()
    
    print("Compiling...")
    code = compile(source, path, "exec")
    print("Executing...")
    exec(code, module.__dict__)
    print("Execution complete.")
    
    print("Instantiating MotorSpecs from dynamic module...")
    specs = module.MotorSpecs()
    print(f"Created MotorSpecs in {time.time() - start:.2f}s")
    
except Exception as e:
    import traceback
    print(f"Error: {e}")
    traceback.print_exc()

print("Diagnostic complete.")
