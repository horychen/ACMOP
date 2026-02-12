import sys
from unittest.mock import MagicMock

# Mock problematic dependencies
sys.modules['JMAG'] = MagicMock()
sys.modules['utility'] = MagicMock()
sys.modules['win32com'] = MagicMock()
sys.modules['win32com.client'] = MagicMock()

# Now we can import user_minitureMachine
from codes4.user_minitureMachine import MotorSpecs

def test_negative_tooth_depth_logic():
    print("--- Testing Negative Tooth Depth Logic ---")
    specs = MotorSpecs()
    
    # r_so (d_stator_outer) is 10.0 by default in MotorSpecs? No, let's check
    print(f"Initial d_stator_outer: {specs.geometry.d_stator_outer.value}")
    
    # Intentionally set a very large yoke depth to make tooth depth negative
    # tooth_depth = r_so - r_si - yoke - shoe
    specs.geometry.stator_yoke_depth.value = 10.0 # Very large
    
    print("\nCalling validate_inputs('geometry')...")
    result = specs.validate_inputs("geometry")
    
    print(f"\nAPI Status: {result['status']}")
    print(f"API Errors: {result['errors']}")
    
    found_expected_error = any("tooth_depth" in err and "non-positive" in err for err in result['errors'])
    if found_expected_error:
        print("\nSUCCESS: Found expected derivation error in API response.")
    else:
        print("\nFAILURE: Did not find expected derivation error in API response.")
        
    print("\n--- Test Complete ---")

if __name__ == "__main__":
    test_negative_tooth_depth_logic()
