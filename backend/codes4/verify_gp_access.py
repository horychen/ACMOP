from machine import Machine
import os

def verify():
    print("Initializing Machine...")
    m = Machine()
    print("Syncing Machine...")
    m.sync()
    
    print("\nVerifying machine.geometry.gp...")
    if hasattr(m.geometry, 'gp'):
        print(f"PASS: machine.geometry has 'gp' attribute.")
        print(f"GP values: {m.geometry.gp}")
        assert 'r_rotor_outer' in m.geometry.gp
        assert 'd_magnet' in m.geometry.gp
        print("PASS: Required keys found in gp.")
    else:
        print("FAIL: machine.geometry missing 'gp' attribute.")
        return False

    print("\nVerifying JMAG.py compatibility (simulated)...")
    # Simulate the calculation in JMAG.py line 452
    try:
        R = m.geometry.gp['r_rotor_outer'] - 0.5 * m.geometry.gp['d_magnet']
        print(f"PASS: Successfully calculated R = {R}")
    except Exception as e:
        print(f"FAIL: Calculation failed with error: {e}")
        return False

    return True

if __name__ == "__main__":
    if verify():
        print("\nVerification SUCCESSFUL!")
    else:
        print("\nVerification FAILED!")
