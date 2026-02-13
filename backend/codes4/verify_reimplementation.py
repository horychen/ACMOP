import sys
import os
import math
import csv
from typing import Dict, Any
from unittest.mock import MagicMock

# Mock win32com and pythoncom if missing
try:
    import win32com
except ImportError:
    sys.modules["win32com"] = MagicMock()
    sys.modules["win32com.client"] = MagicMock()
try:
    import pythoncom
except ImportError:
    sys.modules["pythoncom"] = MagicMock()

# Mock cairo and cairosvg if missing
try:
    import cairo
except ImportError:
    sys.modules["cairo"] = MagicMock()
try:
    import cairosvg
except ImportError:
    sys.modules["cairosvg"] = MagicMock()

# Mock pygmo if missing
try:
    import pygmo
except ImportError:
    sys.modules["pygmo"] = MagicMock()

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from machine_design_guide import Modern_Machine_Designer
from user_minitureMachine import MotorSpecs

def capture_state(designer):
    state = {}
    # Capture Parameter values
    for attr_name, attr_val in designer.__dict__.items():
        if hasattr(attr_val, 'value') and hasattr(attr_val, 'name'):
            state[f"param_{attr_name}"] = attr_val.value
    
    # Capture EX dictionary
    EX = getattr(designer, 'EX', None)
    if EX is None and hasattr(designer, 'winding'):
        EX = getattr(designer.winding, 'EX', None)

    if EX:
        for k, v in EX.items():
            if k == 'mm2_magnet_area': continue # Skip if not in old
            state[f"EX_{k}"] = v
            
    return state

def run_verification():
    print("Running Verification...")
    
    # Initialize MotorSpecs (default)
    specs = MotorSpecs()
    
    # Initialize designer
    # Modern_Machine_Designer(specs=specs)
    designer = Modern_Machine_Designer(specs)
    
    state = capture_state(designer)
    
    csv_file = 'reimplementation_comparison_phase2.csv'
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Variable', 'New Value', 'Status'])
        for k, v in sorted(state.items()):
            writer.writerow([k, v, 'Captured'])
    
    print(f"Comparison CSV generated: {os.path.abspath(csv_file)}")
    
    # Basic assertions to ensure logic is working
    assert designer.mm_r_so.value == specs.geometry.r_stator_outer.value
    assert designer.EX['DCBusVoltage'] == specs.winding.dc_bus_voltage
    
    print("Verification passed (basic checks)!")

if __name__ == "__main__":
    run_verification()
