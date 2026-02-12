# This script is used to generate the geometry of a machine
# other performance and excitation calculations are done in other scripts

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
        alpha_slot_pitch = 360.0 / num_slots if num_slots != 0 else 0.0
        alpha_pole_pitch = 360.0 / num_poles if num_poles != 0 else 0.0
        
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
    
    if isinstance(result, dict):
        instructions = result.get("几何绘制字符串", [])
        bMirror, mirrorAxis = result.get("镜像与否以及镜像轴", (False, None))
        copyCount = result.get("旋转拷贝的个数", 1)
    else:
        instructions = result
        bMirror, mirrorAxis = (False, None)
        copyCount = 1

    list_segments = []
    list_coords = []
    
    drawer.getSketch(part.name, part.color)
    if getattr(drawer, 'verbose_drawing', False):
        print(f"--- Parsing instructions for {part.name} (Base Rotation: {part.rotation_deg} deg, Copies: {copyCount}) ---")

    # Priority for duplication: part attributes (overrides) > draw_instruction return values
    part_mirror = getattr(part, 'bMirror', None)
    part_rotate = getattr(part, 'iRotateCopy', None)
    
    if part_mirror is not None: bMirror = part_mirror
    if part_rotate is not None: copyCount = part_rotate

    # If it's not an FEA drawer (like JMAG), we must handle the duplication manually for visualization
    is_fea = hasattr(drawer, 'prepareSection')
    loop_count = 1 if (is_fea or copyCount <= 1) else copyCount
    
    for i in range(loop_count):
        rotation_offset = i * (360.0 / copyCount)
        current_deg = part.rotation_deg + rotation_offset
        
        for instr in instructions:
            if getattr(drawer, 'verbose_drawing', False) and i == 0:
                print(f"  Instr: {instr}")
            if " ~ " in instr:
                p1_name, p2_name = instr.split(" ~ ")
                p1 = parse_point_name(p1_name, all_points, current_deg)
                p2 = parse_point_name(p2_name, all_points, current_deg)
                if getattr(drawer, 'verbose_drawing', False) and i == 0:
                    print(f"    Resolved Arc: {p1_name}->{p1}, {p2_name}->{p2}")
                list_segments += drawer.drawArc((0,0), p1, p2)
                list_coords.extend([p1, p2])
            elif " - " in instr:
                p1_name, p2_name = instr.split(" - ")
                p1 = parse_point_name(p1_name, all_points, current_deg)
                p2 = parse_point_name(p2_name, all_points, current_deg)
                if getattr(drawer, 'verbose_drawing', False) and i == 0:
                    print(f"    Resolved Line: {p1_name}->{p1}, {p2_name}->{p2}")
                list_segments += drawer.drawLine(p1, p2)
                list_coords.extend([p1, p2])
            
    if not list_coords:
        return {'innerCoord': (0,0), 'list_regions': [[]], 'mirrorAxis': mirrorAxis, 'bMirror': bMirror, 'iRotateCopy': copyCount}
        
    # Calculate innerCoord as centroid of unique points
    unique_coords = list(set([tuple(p) for p in list_coords]))
    sum_x = sum(p[0] for p in unique_coords)
    sum_y = sum(p[1] for p in unique_coords)
    innerCoord = (sum_x / len(unique_coords), sum_y / len(unique_coords))
    
    return {
        'innerCoord': innerCoord, 
        'list_regions': [list_segments], 
        'mirrorAxis': mirrorAxis,
        'bMirror': bMirror,
        'iRotateCopy': copyCount
    }

@dataclass
class MaterialSpecs:
    magnet_grade: str = "N42SH"
    magnet_br: float = 1.3
    magnet_h_cj: float = 1592.0
    stator_steel: str = "20JNEH1200"
    rotor_steel: str = "20JNEH1200"
    steel_stack_factor: float = 0.95
    steel_max_flux_density: float = 1.9
    copper_fill_factor: float = 0.4

@dataclass
class WindingSpecs:
    num_phases: int = 3
    num_slots: int = 12
    num_poles: int = 10
    coil_pitch_y: int = 1
    conductors_per_slot: int = 42
    wire_diameter: float = 0.21
    connection: str = "Wye"
    rated_speed_rpm: float = 3000.0
    rated_current_density_Js: float = 5.0 # [A/mm^2]
    number_of_parallel_branch: int = 1
    
    # Derived results
    phase_resistance: float = 0.0
    estimated_back_emf: float = 0.0
    slot_current_at: float = 0.0

@dataclass
class MachinePart:
    name: str
    all_points: Optional['AllPoints'] = field(default=None, repr=False)
    index: int = field(init=False, default=0)
    rotation_deg: float = 0.0
    options: str = ""
    color: str = "#000000"

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
            "几何绘制字符串": instrs,
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
                "HP[4] - HP[9]",
                "HP[9] ~ RP[9]",
                "RP[9] - RP[8]",
                "RP[8] ~ HP[7]",
                "HP[7] - HP[6]",
                "HP[6] ~ RP[6]",
                "RP[6] - RP[4]",
                "RP[4] ~ HP[4]"
            ]
        elif self.options == "open-slot":
            instrs = [
                "HP[4] ~ RP[4]",
                "RP[4] - HP[7]",
                "HP[7] ~ RP[8]",
                "RP[8] - RP[9]",
                "HP[9] ~ RP[9]",
                "HP[9] - HP[4]"
            ]

        return {
            "几何绘制字符串": instrs,
            "镜像与否以及镜像轴": (True, None),
            "旋转拷贝的个数": num_slots
        }

@dataclass
class Magnet(MachinePart):
    color: str = "#2222BB" # Blue
    def draw_instruction(self):
        num_poles = self.all_points.num_poles
        return {
            "几何绘制字符串":
            [
            "RP[2]_M - RP[3]_M", 
            "RP[3]_M ~ RP[3]",
            "RP[3] - RP[2]", 
            "RP[2]_M ~ RP[2]"
        ],
        "镜像与否以及镜像轴": (False, None),
        "旋转拷贝的个数": num_poles,
        }

@dataclass
class Coil(MachinePart):
    color: str = "#B87333" # Copper
    def draw_instruction(self):
        num_slots = self.all_points.num_slots
        return {
            "几何绘制字符串":
            [
            "HP[6] - HP[7]", "HP[7] ~ RP[8]", "RP[8] - RP[6]", "RP[6] ~ HP[6]",
            "HP[4] ~ RP[4]" # Slot reference
        ],
        "镜像与否以及镜像轴": (False, None),
        "旋转拷贝的个数": num_slots,
        }

class Machine:
    def __init__(self, all_points: Optional[AllPoints] = None):
        self.parts = []
        self._next_index = 0
        self.all_points = all_points
        self.materials = MaterialSpecs()
        self.winding = WindingSpecs()
        self.l_stack = 16.0 # [mm]

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

    def sync(self):
        """Update derived performance parameters."""
        if not self.all_points: return
        
        g = self.all_points.HP
        w = self.winding
        m = self.materials
        
        # 1. Geometry-derived values
        r_si = g[4][0]
        r_ro = g[3][0]
        hm = g[3][0] - g[2][0]
        gap_dist = r_si - r_ro
        
        # 2. Bg Estimation
        br = m.magnet_br
        b_gap = br * hm / (hm + 1.05 * gap_dist) if (hm + 1.05 * gap_dist) > 0 else 0
        
        # 3. Electrical Metrics
        wire_area = math.pi * (w.wire_diameter/2)**2
        conductor_current = w.rated_current_density_Js * wire_area # [A]
        w.slot_current_at = conductor_current * w.conductors_per_slot
        
        # 4. Resistance Estimation
        end_winding = math.pi * (2*r_si) / w.num_slots * w.coil_pitch_y
        n_series = (w.conductors_per_slot * w.num_slots) / (2 * w.num_phases * w.number_of_parallel_branch)
        turn_length = 2 * (self.l_stack + end_winding) * 1e-3 # [m]
        rho_copper = 1.72e-8 # [Ohm-m]
        w.phase_resistance = rho_copper * (turn_length * n_series) / (wire_area * 1e-6)
        
        # 5. Back-EMF
        area_pole = (2 * math.pi * r_si * self.l_stack) / w.num_poles
        estimated_flux = b_gap * area_pole * 1e-6 
        angular_speed = w.rated_speed_rpm * (2 * math.pi / 60.0)
        winding_factor = 0.95 
        w.estimated_back_emf = estimated_flux * n_series * winding_factor * angular_speed

    def show_geometry_svg(self, filename='machine_geometry.svg', scale=10.0):
        from modern_machine_designer_utility import CairoDrawer
        drawer = CairoDrawer(filename=filename, scale=scale)
        drawer.draw_machine(self)
        print(f"Geometry drawn to {filename} with scale {scale}")
