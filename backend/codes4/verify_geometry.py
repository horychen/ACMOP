import sys
import os
import math
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

def verify_geometry():
    print("Initializing MotorSpecs...")
    specs = MotorSpecs()
    
    # Initialize designer with specs
    print("Initializing Modern_Machine_Designer...")
    designer = Modern_Machine_Designer(specs)
    
    # Generate geometry SVG
    print("Generating geometry SVG...")
    svg_filename = 'machine_geometry_test.svg'
    # We need to make sure show_geometry uses the filename we want
    # The current show_geometry has a hardcoded filename in draw_spmsm inner function:
    # def draw_spmsm(lw, width_in_points, height_in_points, filename='machine_geometry.svg', bool_draw_whole_model=True):
    # But it accepts filename in show_geometry(self, filename=None, x_denorm_dict=None)
    # However, look at the implementation:
    # draw_spmsm(lw, width_in_points, height_in_points, bool_draw_whole_model=bool_draw_whole_model)
    # It doesn't pass filename down!
    
    designer.show_geometry()
    
    if os.path.exists('machine_geometry.svg'):
        print(f"Success! Geometry SVG generated: {os.path.abspath('machine_geometry.svg')}")
    else:
        print("Error: machine_geometry.svg was not generated.")

if __name__ == "__main__":
    verify_geometry()
