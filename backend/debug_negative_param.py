from codes4.user_minitureMachine import MotorSpecs
import math

def test_negative_tooth_depth():
    print("--- Testing Negative Tooth Depth ---")
    specs = MotorSpecs()
    
    # Intentionally set a very large yoke depth to make tooth depth negative
    # r_so is approx 10.0, r_si is approx 6.15. 
    # tooth_depth = 10 - 6.15 - yoke - shoe
    specs.geometry.stator_yoke_depth.value = 5.0 
    
    print("Calling validate_inputs('geometry')...")
    specs.validate_inputs("geometry")
    print("--- Test Complete ---")

if __name__ == "__main__":
    test_negative_tooth_depth()
