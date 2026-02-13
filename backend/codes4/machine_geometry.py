# This script is used to generate the geometry of a machine
# other performance and excitation calculations are done in other scripts
# 
# IMPORTANT: All arcs (marked with "~") must be drawn in Counter-Clockwise (CCW) fashion.
# This means the first point should have a smaller angle than the second point 
# in the context of the arc segment being drawn (e.g., lower_angle ~ higher_angle).

from dataclasses import dataclass, field, InitVar
from typing import List, Dict, Any, Optional
import math
def rotate_point(p, deg):
    rad = math.radians(deg)
    x, y = p
    return (x * math.cos(rad) - y * math.sin(rad),
            x * math.sin(rad) + y * math.cos(rad))

@dataclass
class AllPoints:
    required_GP: InitVar[dict]
    horizontal_position: Dict[int, tuple] = field(default_factory=dict, init=False)
    rotated_position: Dict[int, tuple] = field(default_factory=dict, init=False)

    def __post_init__(self, required_GP: dict):
        # Helper for rotation around origin
        def rotate(r, deg):
            rad = math.radians(deg)
            return (r * math.cos(rad), r * math.sin(rad))

        self.required_GP = required_GP
        r_shaft = required_GP['r_shaft']
        r_rotor_outer = required_GP['r_rotor_outer']
        d_magnet = required_GP['d_magnet']
        d_air_gap = required_GP['d_air_gap']
        r_stator_outer = required_GP['r_stator_outer']
        d_stator_yoke = required_GP['d_stator_yoke']
        d_stator_tooth = required_GP['d_stator_tooth']
        d_stator_tooth_shoe = required_GP['d_stator_tooth_shoe']
        w_stator_width = required_GP['w_stator_width']
        self.num_slots = num_slots = required_GP['num_slots']
        self.num_poles = num_poles = required_GP['num_poles']
        
        # Calculate key angles
        alpha_slot_pitch = 360.0 / num_slots
        alpha_pole_pitch = 360.0 / num_poles
        
        # Stator radii
        r_si = r_rotor_outer + d_air_gap
        r_ss = r_si + d_stator_tooth_shoe
        r_sy = r_stator_outer - d_stator_yoke
        r_so = r_stator_outer

        # Angles for the stalk boundaries (where the straight tooth meets the arcs)
        alpha_si_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_si)) if r_si > 0 else 0.0
        alpha_ss_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_ss)) if r_ss > 0 else 0.0
        alpha_sy_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_sy)) if r_sy > 0 else 0.0
        alpha_so_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_so)) if r_so > 0 else 0.0

        # HP: Horizontal Points (Usually angle 0, but here used for corners on the tooth axis)
        self.HP = self.horizontal_position = {
            0: (0.0, 0.0),
            1: (r_shaft, 0.0),
            2: (r_rotor_outer - d_magnet, 0.0),
            3: (r_rotor_outer, 0.0),
            4: (r_si, 0.0),
            5: (r_ss, 0.0),
            6: rotate(r_ss, alpha_ss_deg), # Shoe corner
            7: rotate(r_sy, alpha_sy_deg), # Yoke corner
            8: (r_sy, 0.0),
            9: (r_so, 0.0),
        }

        # Mirrored points for symmetry and circle drawing
        self.HP_mirror = {
            1: (-self.HP[1][0], self.HP[1][1]), # 180-deg for circles
            2: (-self.HP[2][0], self.HP[2][1]),
            3: (-self.HP[3][0], self.HP[3][1]),
            4: (-self.HP[4][0], self.HP[4][1]),
            6: (self.HP[6][0], -self.HP[6][1]), # Tooth stalk symmetry
            7: (self.HP[7][0], -self.HP[7][1]), # Tooth stalk symmetry
            9: (-self.HP[9][0], self.HP[9][1]),
        }

        # RP: Rotated Points (at slot center axis, angle half_pitch)
        stator_half_pitch = alpha_slot_pitch / 2.0
        rotor_half_pitch = alpha_pole_pitch / 2.0
        self.RP = self.rotated_position = {
            0: rotate(0.0, stator_half_pitch),
            1: rotate(r_shaft, stator_half_pitch),
            2: rotate(r_rotor_outer - d_magnet, rotor_half_pitch),
            3: rotate(r_rotor_outer, rotor_half_pitch),
            4: rotate(r_si, stator_half_pitch),
            5: rotate(r_ss, stator_half_pitch),
            6: rotate(r_ss, stator_half_pitch),
            7: rotate(r_sy, stator_half_pitch),
            8: rotate(r_sy, stator_half_pitch),
            9: rotate(r_so, stator_half_pitch),

            14: rotate(r_si, alpha_si_deg),
            15: rotate(r_ss, alpha_ss_deg),
            17: rotate(r_sy, alpha_sy_deg),
            19: rotate(r_so, alpha_so_deg),
        }


def parse_point_name(name, all_points, rotation_deg=0):
    name = name.split('(')[0].strip()
    is_mirror = False
    if name.endswith("_M"):
        is_mirror = True
        name = name[:-2]
    
    # Support HP_M[...] and RP_M[...] syntax
    if name.startswith("HP_M["):
        is_mirror = True
        name = "HP[" + name[5:]
    elif name.startswith("RP_M["):
        is_mirror = True
        name = "RP[" + name[5:]

    if name.startswith("HP[") and name.endswith("]"):
        idx = int(name[3:-1])
        p = all_points.HP[idx]
        if is_mirror:
            if hasattr(all_points, 'HP_mirror') and idx in all_points.HP_mirror:
                p = all_points.HP_mirror[idx]
            else:
                p = (p[0], -p[1])
    elif name.startswith("RP[") and name.endswith("]"):
        idx = int(name[3:-1])
        p = all_points.RP[idx]
        if is_mirror:
            # We don't typically need a separate RP_mirror dictionary
            # as RP points can be mirrored across the Tooth Axis via Y-negation
            p = (p[0], -p[1])
    else:
        raise ValueError(f"Unknown point name: {name}")
    
    return rotate_point(p, rotation_deg)

def draw_instruction_parser(part, all_points, drawer):
    result = part.draw_instruction()
    
    # Initialize bMirror, mirrorAxis, copyCount from result first
    if isinstance(result, dict):
        bMirror, mirrorAxis = result.get("镜像与否以及镜像轴", (False, None))
        copyCount = result.get("旋转拷贝的个数", 1)
    else:
        bMirror, mirrorAxis = (False, None)
        copyCount = 1

    list_regions = []
    region_names = []
    if isinstance(result, dict) and "区域的几何绘制指导" in result:
        dict_regions = result["区域的几何绘制指导"]
        instructions_list = list(dict_regions.values())
        region_names = list(dict_regions.keys())
    else:
        # Fallback to legacy single list format
        scalar_instrs = result.get("几何绘制字符串", []) if isinstance(result, dict) else result
        instructions_list = [scalar_instrs]
        region_names = ["region1"]

    list_coords = []
    region_centroids = {}
    
    drawer.getSketch(part.name, part.color)
    if getattr(drawer, 'verbose_drawing', False):
        print(f"--- Parsing instructions for {part.name} (Base Rotation: {part.rotation_deg} deg, Copies: {copyCount}) ---")

    # Priority for duplication: part attributes (overrides) > draw_instruction return values
    part_mirror = getattr(part, 'bMirror', None)
    part_rotate = getattr(part, 'iRotateCopy', None)
    
    if part_mirror is not None: bMirror = part_mirror
    if part_rotate is not None: copyCount = part_rotate

    # If it's a visualization drawer (like ReactDrawer), we handle duplication manually
    # If it's a CAD drawer (like JMAG), we let the CAD tool handle it
    is_fea = hasattr(drawer, 'prepareSection')
    # ReactDrawer specifically needs the data for each copy to be provided if it doesn't handle them
    # But ReactDrawer handles copies currently via collecting all regions.
    # Actually, CairoDrawer and JMAG prefer not to have duplicate geometry definitions.
    # We'll treat ReactDrawer as a CAD drawer for this logic to keep it consistent.
    loop_count = 1 if (is_fea or copyCount <= 1) else copyCount
    
    for i in range(loop_count):
        rotation_offset = i * (360.0 / copyCount)
        current_deg = part.rotation_deg + rotation_offset
        
        for idx, instructions in enumerate(instructions_list):
            region_segments = []
            region_coords = []
            for instr in instructions:
                # Strip comments if any
                instr = instr.split('#')[0].strip()
                if not instr: continue
                
                if getattr(drawer, 'verbose_drawing', False) and i == 0:
                    print(f"  Instr: {instr}")
                if " ~ " in instr:
                    p1_name, p2_name = instr.split(" ~ ")
                    p1 = parse_point_name(p1_name, all_points, current_deg)
                    p2 = parse_point_name(p2_name, all_points, current_deg)
                    region_segments += drawer.drawArc((0,0), p1, p2)
                    region_coords.extend([p1, p2])
                elif " - " in instr:
                    p1_name, p2_name = instr.split(" - ")
                    p1 = parse_point_name(p1_name, all_points, current_deg)
                    p2 = parse_point_name(p2_name, all_points, current_deg)
                    region_segments += drawer.drawLine(p1, p2)
                    region_coords.extend([p1, p2])
            
            list_regions.append(region_segments)
            list_coords.extend(region_coords)
            
            # Calculate centroid for this region (in the first copy)
            if i == 0 and region_coords:
                unique_reg_coords = list(set([tuple(p) for p in region_coords]))
                rx = sum(p[0] for p in unique_reg_coords) / len(unique_reg_coords)
                ry = sum(p[1] for p in unique_reg_coords) / len(unique_reg_coords)
                region_centroids[region_names[idx]] = (rx, ry)
    
    if not list_coords:
        return {'innerCoord': (0,0), 'inner_coords': {}, 'list_regions': [[]], 'mirrorAxis': mirrorAxis, 'bMirror': bMirror, 'iRotateCopy': copyCount}
        
    # Calculate overall innerCoord
    unique_coords = list(set([tuple(p) for p in list_coords]))
    innerCoord = (sum(p[0] for p in unique_coords) / len(unique_coords), 
                  sum(p[1] for p in unique_coords) / len(unique_coords))
    
    return {
        'innerCoord': innerCoord, 
        'inner_coords': region_centroids,
        'list_regions': list_regions,
        'mirrorAxis': mirrorAxis,
        'bMirror': bMirror,
        'iRotateCopy': copyCount
    }

@dataclass
class MachinePart:
    name: str
    all_points: Optional['AllPoints'] = field(default=None, repr=False)
    index: int = field(init=False, default=0)
    rotation_deg: float = 0.0
    options: str = ""
    color: str = "#000000"
    inner_coords: Optional[Dict[str, tuple]] = field(default_factory=dict) # Centroids for each region

@dataclass
class RotorCore(MachinePart):
    color: str = "#555555" # Iron gray
    def draw_instruction(self):
        instrs = []
        if self.options == "notched":
            instrs = [
                "RP[1]_M - RP[2]_M", "RP[2]_M ~ RP[2]", "RP[2] - RP[1]", "RP[1]_M ~ RP[1]"
            ]
        elif self.options == "cylinder":
            instrs = ["HP[2] ~ HP[2]_M", "HP[2]_M ~ HP[2]"]
            
        return {
            "区域的几何绘制指导": {"region1": instrs},
            "镜像与否以及镜像轴": (False, None),
            "旋转拷贝的个数": 1 if self.options == "cylinder" else (self.all_points.num_poles if self.all_points else 1)
        }

@dataclass
class StatorCore(MachinePart):
    color: str = "#666666" # Slightly lighter gray
    def draw_instruction(self):
        instrs = []
        num_slots = self.all_points.num_slots if self.all_points else 1
        bFullCircle = False
        
        if self.options == "closed-slot" or self.options == "semi-closed-slot":
            instrs = [
                # Start from tooth center outer back iron and follow a continuous CCW loop for half a sector
                "HP[9] - HP[4]",  # Radial center line: Back iron outer -> Gap inner (IN)
                "HP[4] ~ RP[4]",  # Air gap arc (CCW)
                "RP[4] - RP[6]",  # Slot side radial line (OUT)
                "HP[6] ~ RP[6]",  # Shoe arc (CCW)
                "HP[6] - HP[7]",  # Tooth stalk side radial line (OUT)
                "HP[7] ~ RP[8]",  # Back iron inner arc (CCW)
                "RP[8] - RP[9]",  # Back iron slot center radial line (OUT)
                "HP[9] ~ RP[9]"   # Back iron outer arc (CCW)
            ]
        elif self.options == "open-slot":
            instrs = [
                "HP[9] - HP[4]",  # Radial center line (IN)
                "HP[4] ~ RP[4]",  # Air gap arc (CCW)
                "RP[4] - HP[7]",  # Open slot side (OUT) - Note: slightly different than closed
                "HP[7] ~ RP[8]",  # Back iron inner arc (CCW)
                "RP[8] - RP[9]",  # Radial (OUT)
                "HP[9] ~ RP[9]"   # Outer Arc (CCW)
            ]

        return {
            "区域的几何绘制指导": {"region1": instrs},
            "镜像与否以及镜像轴": (True, None),
            "旋转拷贝的个数": num_slots
        }

@dataclass
class Magnet(MachinePart):
    color: str = "#2222BB" # Blue
    def draw_instruction(self):
        num_poles = self.all_points.num_poles
        return {
            "区域的几何绘制指导":
            {
                "region1":
                [
                    "RP[2]_M - RP[3]_M", 
                    "RP[3]_M ~ RP[3]",
                    "RP[3] - RP[2]", 
                    "RP[2]_M ~ RP[2]"
                ]
            },
            "镜像与否以及镜像轴": (False, None),
            "旋转拷贝的个数": num_poles,
        }

@dataclass
class Coil(MachinePart):
    color: str = "#B87333" # Copper
    def draw_instruction(self):
        num_slots = self.all_points.num_slots
        return {
            "区域的几何绘制指导":
            {
                "region1":
                [
                    "HP[6] - HP[7]", 
                    "HP[7] ~ RP[8]", 
                    "RP[8] - RP[6]", 
                    "HP[6] ~ RP[6]",
                ],
                "region2":
                [
                    "HP_M[6] - HP_M[7]", 
                    "RP_M[8] ~ HP_M[7]", 
                    "RP_M[8] - RP_M[6]", 
                    "RP_M[6] ~ HP_M[6]", 

                ],
            },
            "镜像与否以及镜像轴": (False, None),
            "旋转拷贝的个数": num_slots,
        }

@dataclass
class MachineGeometry:
    parts: List[MachinePart] = field(default_factory=list)
    _next_index: int = field(init=False, default=0)
    all_points: Optional[AllPoints] = None
    gp: dict = field(default_factory=dict) # Parameters from target
    l_stack: float = 16.0 # [mm]
    num_slots: int = 0
    num_poles: int = 0

    # These will be added dynamically in sync() as Parameter objects for compatibility
    r_stator_outer: Any = field(init=False, default=None)
    r_rotor_outer: Any = field(init=False, default=None)
    d_magnet: Any = field(init=False, default=None)
    d_air_gap: Any = field(init=False, default=None)
    split_ratio: Any = field(init=False, default=None)
    w_tooth: Any = field(init=False, default=None)
    d_tooth: Any = field(init=False, default=None)
    d_tooth_shoe: Any = field(init=False, default=None)
    d_stator_yoke: Any = field(init=False, default=None)
    tooth_shape: str = "closed"
    r_shaft: Any = field(init=False, default=None)

    @property
    def machineGeometry(self):
        """Legacy compatibility: returns a dict of parts by name."""
        return {part.name: part for part in self.parts}

    def add_part(self, part: MachinePart):
        part.index = self._next_index
        part.all_points = self.all_points # Link global points
        self.parts.append(part)
        self._next_index += 1
        return part

    def add_radial_array(self, part_type, base_name, count, options="", color=None):
        created_parts = []
        for i in range(count):
            rotation = i * 360.0 / count
            part = part_type(name=f"{base_name}_{i}", rotation_deg=rotation, options=options)
            if color: part.color = color
            self.add_part(part)
            created_parts.append(part)
        return created_parts

    def show_geometry_svg(self, filename='machine_geometry.svg', scale=10.0):
        from modern_machine_designer_utility import CairoDrawer
        drawer = CairoDrawer(filename=filename, scale=scale)
        drawer.draw_machine(self)
        print(f"Geometry drawn to {filename} with scale {scale}")

    def sync(self, gp: dict):
        """Populate compatibility Parameter objects from GP dict."""
        from modern_machine_designer_utility import Parameter
        self.gp = gp
        self.r_stator_outer = Parameter("stator_outer_radius", "fixed", gp['r_stator_outer'])
        self.r_rotor_outer = Parameter("rotor_outer_radius", "free", gp['r_rotor_outer'])
        self.d_magnet = Parameter("magnet_depth", "free", gp['d_magnet'])
        self.d_air_gap = Parameter("mechanical_air_gap_depth", "fixed", gp['d_air_gap'])
        self.split_ratio = Parameter("split_ratio", "derived", gp['r_rotor_outer'] / gp['r_stator_outer'] if gp['r_stator_outer'] else 0.615)
        
        self.w_tooth = Parameter("stator_tooth_width", "free", gp['w_stator_width'])
        self.d_tooth_shoe = Parameter("stator_tooth_shoe_depth", "fixed", gp['d_stator_tooth_shoe'])
        self.d_stator_yoke = Parameter("stator_yoke_depth", "fixed", gp['d_stator_yoke'])
        self.tooth_shape = "closed"
        
        # d_tooth formula matches AllPoints formula: r_so - r_ro - g - dy - dt_shoe
        self.d_tooth = Parameter("stator_tooth_depth", "derived", gp['r_stator_outer'] - (gp['r_rotor_outer'] + gp['d_air_gap']) - gp['d_stator_yoke'] - gp['d_stator_tooth_shoe'])
        self.r_shaft = Parameter("shaft_radius", "fixed", gp['r_shaft'])
        self.deg_alpha_rm = Parameter("magnet_span_angle", "fixed", gp.get('deg_alpha_rm', 18.0))
        self.d_sleeve = Parameter("rotor_sleeve_depth", "fixed", gp.get('d_sleeve', 0.0))

@dataclass
class MotorSpecs:
    """Compatibility class to match legacy MotorSpecs interface."""
    def __init__(self):
        from machine_materials import MachineMaterial
        from machine_target import MachineTarget
        from machine_winding import MachineWinding
        
        self.materials = MachineMaterial()
        self.target = MachineTarget()
        self.winding = MachineWinding()
        self.geometry = MachineGeometry()
        
        # Sync immediately
        self.sync()

    @property
    def targets(self):
        """Legacy alias for target."""
        return self.target

    def sync(self):
        """Synchronize all components."""
        gp = self.target.get_required_GP()
        self.geometry.all_points = AllPoints(required_GP=gp)
        # Update winding/poles if they changed in target
        self.winding.slot_count = self.target.fixed_parameters['num_slots']
        self.winding.pole_count = self.target.fixed_parameters['num_poles']
        self.geometry.num_slots = self.winding.slot_count
        self.geometry.num_poles = self.winding.pole_count
        
        # Use geometry's own sync
        self.geometry.sync(gp)
        
        # Populate parts
        self.geometry.parts = []
        self.geometry.add_part(RotorCore(name="rotorCore"))
        self.geometry.add_part(StatorCore(name="statorCore"))
        self.geometry.add_part(Magnet(name="rotorMagnet"))
        self.geometry.add_part(Coil(name="coils"))
        
        # Update winding aliases
        self.winding.num_slots = self.winding.slot_count
        self.winding.num_poles = self.winding.pole_count
        self.winding.wire_diameter_with_insulation = self.winding.wire_diameter + 0.02
        
        self.winding.sync(
            machineGeometry=self.geometry,
            materials=self.materials,
            l_stack=self.geometry.l_stack
        )

    def validate_inputs(self, step: str = "geometry") -> dict:
        """Compatibility wrapper for validation."""
        # Simple implementation for now, mirroring legacy structure
        self.sync()
        
        # Basic validation logic based on the step
        errors = []
        warnings = []
        metrics = {}
        
        if step == "geometry":
            gp = self.target.get_required_GP()
            metrics["Stator OR"] = f"{gp['r_stator_outer']} mm"
            metrics["Rotor OR"] = f"{gp['r_rotor_outer']} mm"
            metrics["Air Gap"] = f"{gp['d_air_gap']} mm"
            
        return {
            "status": "fail" if errors else ("warn" if warnings else "pass"),
            "errors": errors,
            "warnings": warnings,
            "metrics": metrics,
            "fill_factor": "40%" # Placeholder
        }
