
from machine import Machine
import os

def test_sync():
    print("--- Initializing Machine ---")
    m = Machine()
    
    print("--- Running sync_all_points ---")
    m.sync_all_points()
    
    print("--- Running sync ---")
    # This calls winding.sync which we refactored
    m.sync()
    
    print("--- Verification Results ---")
    print(f"Slot Count: {m.winding.slot_count}")
    print(f"Pole Count: {m.winding.pole_count}")
    print(f"Estimated Back-EMF: {m.winding.estimated_back_emf:.2f} V")
    print(f"Phase Resistance: {m.winding.phase_resistance:.3f} Ohm")
    print(f"Slot Area: {m.winding.slot_area:.2f} mm2")
    
    # Check if some values are calculated (not zero/default)
    assert m.winding.estimated_back_emf > 0
    assert m.winding.phase_resistance > 0
    assert m.winding.slot_area > 0
    print("Sync verification successful!")

if __name__ == "__main__":
    try:
        test_sync()
    except Exception as e:
        import traceback
        traceback.print_exc()
