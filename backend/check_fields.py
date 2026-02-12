import dataclasses
import sys
import os

# Import user_minitureMachine
sys.path.append(os.path.join(os.getcwd(), "backend", "codes4"))
import user_minitureMachine

print(f"File path: {user_minitureMachine.__file__}")

g = user_minitureMachine.GeometrySpecs()
print(f"Instance attributes: {list(g.__dict__.keys())}")
print(f"Class fields: {[f.name for f in dataclasses.fields(user_minitureMachine.GeometrySpecs)]}")

# Try asdict
try:
    d = dataclasses.asdict(g)
    print("asdict success")
except Exception as e:
    import traceback
    print(f"asdict failed: {e}")
    traceback.print_exc()
