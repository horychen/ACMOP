import sys
import os
import math
from collections import OrderedDict
from typing import List, Dict, Any, Optional

# Import CairoDrawer from the local copy
from modern_machine_designer_utility import CairoDrawer

def rotate_point(p, deg):
    """
    二维平面坐标旋转函数。
    将一个点 p=(x,y) 绕原点(0,0) 逆时针旋转 deg 度。
    主要用于将基准槽位(或者磁极位置)的几何坐标旋转复制到整个圆周上。
    """
    rad = math.radians(deg)
    x, y = p
    return (x * math.cos(rad) - y * math.sin(rad),
            x * math.sin(rad) + y * math.cos(rad))

def calculate_all_points(machine_dict: OrderedDict) -> OrderedDict:
    """
    几何特征点提取与计算核心引擎。
    目标：根据输入的电磁与几何设计参数 (GP / WP)，计算出电机截面上关键的骨架点(Nodes)。
    返回: 包含水平参考点(HP), 旋转参考点(RP)及其镜像点的坐标字典，供渲染器或有限元软件在绘制线段与圆弧时引用。
    """
    def rotate(r, deg):
        # 极坐标转换为笛卡尔坐标的快捷函数
        rad = math.radians(deg)
        return (r * math.cos(rad), r * math.sin(rad))

    gp = machine_dict['geometry']
    wp = machine_dict['winding']
    
    # 解析各项径向尺寸
    r_shaft = gp.get('r_shaft', 0.0) 
    r_rotor_outer = gp['r_rotor_outer'] # 转子外径
    d_magnet = gp['d_magnet']           # 磁钢厚度
    
    # 定子尺寸推导
    r_si = r_rotor_outer + gp['d_air_gap'] # 定子内径(定子牙齿末端)
    r_so = gp['r_stator_outer']            # 定子外径
    d_stator_tooth = gp['d_tooth']         # 定子齿牙深度
    
    # 轭部径向厚度推导 = 定子外径 - (定子内径 + 齿深)
    d_stator_yoke = r_so - (r_si + d_stator_tooth)
    tooth_shape = gp.get('tooth_shape', 'open')
    d_stator_tooth_shoe = None if tooth_shape == 'open' else gp.get('d_tooth_shoe')
    r_ss = r_si + d_stator_tooth_shoe if d_stator_tooth_shoe is not None else r_si
    r_sy = r_so - d_stator_yoke            # 轭部内缘(即槽底)所在的半径
    w_stator_width = gp['w_tooth']         # 齿宽
    
    slot_count = wp['slot_count']
    pole_count = wp['pole_count']
    
    # 槽距角与极距角 (机械角度)
    alpha_slot_pitch = 360.0 / slot_count
    alpha_pole_pitch = 360.0 / pole_count

    # 基于定子平直等宽齿设计，逆向推导每一个分段点所在的弦角
    # w_stator_width的半宽所在坐标正好对应到角度：asin( half_width / R )
    alpha_si_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_si)) if r_si > 0 else 0.0
    alpha_ss_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_ss)) if r_ss > 0 else 0.0
    alpha_sy_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_sy)) if r_sy > 0 else 0.0
    alpha_so_deg = math.degrees(math.asin((w_stator_width * 0.5) / r_so)) if r_so > 0 else 0.0

    # HP (Horizontal Points) - 定义在X轴正向或其附近的特征点集
    # 数字代号参考了传统电机二维建模中的节点标号惯例
    HP = {
        0: (0.0, 0.0),                 # 原点
        1: (r_shaft, 0.0),             # 轴孔内圈的起始点
        2: (r_rotor_outer - d_magnet, 0.0), # 磁钢底部的起始点
        3: (r_rotor_outer, 0.0),       # 转子最外侧边缘(气隙下方)的第一点
        4: (r_si, 0.0),                # 定子内孔(齿端中心)气隙上方的起始点
        5: (r_ss, 0.0),                # 极靴弧结束点
        6: rotate(r_ss, alpha_ss_deg), # 齿宽侧边在定子内圆处的交点
        7: rotate(r_sy, alpha_sy_deg), # 齿宽侧边在槽底所在半径的交点
        8: (r_sy, 0.0),                # 槽中心线(X轴正向)与槽底背铁的交点
        9: (r_so, 0.0),                # 电机最外侧背线定子外圆中心点
    }

    # HP_mirror - HP点集关于X轴的翻转镜像点
    # 对于带有对称性的直缝，可以使用镜像迅速构建槽形
    HP_mirror = {
        1: (-r_shaft, 0.0),
        2: (-(r_rotor_outer - d_magnet), 0.0),
        3: (-r_rotor_outer, 0.0),
        4: (-HP[4][0], HP[4][1]),
        6: (HP[6][0], -HP[6][1]),
        7: (HP[7][0], -HP[7][1]),
        9: (-HP[9][0], HP[9][1]),
    }

    # RP (Rotated Points) - 将HP集中的特征点沿圆周旋转一半的槽/极距，到达边界缝隙
    # 这是由于通常抽取一个完整的槽极模型进行分析，其边界为 0 到 pitch 角度或 symmetric pitch
    stator_half_pitch = alpha_slot_pitch / 2.0
    rotor_half_pitch = alpha_pole_pitch / 2.0
    
    RP = {
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

    all_points = OrderedDict()
    all_points['HP'] = HP
    all_points['HP_mirror'] = HP_mirror
    all_points['RP'] = RP
    all_points['slot_count'] = slot_count
    all_points['pole_count'] = pole_count

    return all_points

def parse_point_name(name, all_points, rotation_deg=0):
    """
    点名解析器: 将形如 "HP[4]", "RP_M[6]" 这样的伪代码表达转换成带有坐标偏移的具体数学点位(x,y)
    同时可以叠加上 component 本身的旋转角度(用于整机循环画图)。
    _M 表示取镜像。
    """
    name = name.split('(')[0].strip()
    is_mirror = False
    
    # 处理不同写法的镜像后缀
    if name.endswith("_M"):
        is_mirror = True
        name = name[:-2]
    
    if name.startswith("HP_M["):
        is_mirror = True
        name = "HP[" + name[5:]
    elif name.startswith("RP_M["):
        is_mirror = True
        name = "RP[" + name[5:]

    if name.startswith("HP[") and name.endswith("]"):
        idx = int(name[3:-1])
        p = all_points['HP'][idx]
        if is_mirror:
            if 'HP_mirror' in all_points and idx in all_points['HP_mirror']:
                p = all_points['HP_mirror'][idx]
            else:
                p = (p[0], -p[1])
    elif name.startswith("RP[") and name.endswith("]"):
        idx = int(name[3:-1])
        p = all_points['RP'][idx]
        if is_mirror:
            p = (p[0], -p[1])
    else:
        raise ValueError(f"Unknown point name: {name}")
    
    # 解析出点位后，赋予当前的圆周偏移旋转操作
    return rotate_point(p, rotation_deg)

# ----------------------------------------------------------------------------------------------------------
# 下列函数定义了部件的几何“绘制笔画 (Draw Instructions)”。
# 使用基于我们自定义点位命名的伪代码，配合特定操作符来描述形状:
#   "-" : 代表两点间有一条直线
#   "~" : 代表两点间有一条以原点为圆心的弧线
# 
# 返回结构:
# - "区域的几何绘制指导": 定义内部由哪些region闭合连线组成
# - "镜像与否以及镜像轴": 是否需要镜像复制来构成完整部件
# - "旋转拷贝的个数": 复制数量 (例如定子为槽数，转子磁钢为极数)
# ----------------------------------------------------------------------------------------------------------

def get_stator_core_instruction(part_dict, all_points):
    slot_count = all_points['slot_count']
    options = part_dict.get('options', '')
    
    if options == "closed-slot" or options == "semi-closed-slot":
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
    elif options == "open-slot":
        instrs = [
            "HP[9] - HP[4]",  # Radial center line (IN)
            "HP[4] ~ RP[4]",  # Air gap arc (CCW)
            "RP[4] - HP[7]",  # Open slot side (OUT) - Note: slightly different than closed
            "HP[7] ~ RP[8]",  # Back iron inner arc (CCW)
            "RP[8] - RP[9]",  # Radial (OUT)
            "HP[9] ~ RP[9]"   # Outer Arc (CCW)
        ]
    else:
        # Default fallback to open-slot
        instrs = [
            "HP[9] - HP[4]",
            "HP[4] ~ RP[4]",
            "RP[4] - HP[7]",
            "HP[7] ~ RP[8]",
            "RP[8] - RP[9]",
            "HP[9] ~ RP[9]"
        ]

    return {
        "区域的几何绘制指导": {"region1": instrs},
        "镜像与否以及镜像轴": (True, None), # 定子结构是对称的，直接通过镜像合并
        "旋转拷贝的个数": slot_count
    }

def get_rotor_core_instruction(part_dict, all_points):
    options = part_dict.get('options', '')
    instrs = []
    if options == "notched":
        instrs = [
            "RP_M[1] - RP_M[2]", "RP_M[2] ~ RP[2]", "RP[2] - RP[1]", "RP_M[1] ~ RP[1]"
        ]
    elif options == "cylinder" or options == "":
        instrs = ["HP[2] ~ HP_M[2]", "HP_M[2] ~ HP[2]"] # 画一个完整的内圆柱体
        
    return {
        "区域的几何绘制指导": {"region1": instrs},
        "镜像与否以及镜像轴": (False, None),
        "旋转拷贝的个数": 1 if (options == "cylinder" or options == "") else all_points['pole_count']
    }

def get_magnet_instruction(part_dict, all_points):
    pole_count = all_points['pole_count']
    return {
        "区域的几何绘制指导":
        {
            "region1":
            [
                "RP_M[2] - RP_M[3]",  # 磁钢左侧直线
                "RP_M[3] ~ RP[3]",    # 磁钢贴着气隙外围的圆弧
                "RP[3] - RP[2]",      # 磁钢右侧直线
                "RP_M[2] ~ RP[2]"     # 磁钢底部贴着铁芯位置的封闭圆弧
            ]
        },
        "镜像与否以及镜像轴": (False, None),
        "旋转拷贝的个数": pole_count,
    }

def get_coil_instruction(part_dict, all_points):
    slot_count = all_points['slot_count']
    return {
        "区域的几何绘制指导":
        {
            # 定义分布在前一侧的一半导线捆绑截面
            "region1":
            [
                "HP[6] - HP[7]", 
                "HP[7] ~ RP[8]", 
                "RP[8] - RP[6]", 
                "HP[6] ~ RP[6]",
            ],
            # 定义分布在后侧(镜像侧)的一半导线捆绑截面
            "region2":
            [
                "HP_M[6] - HP_M[7]", 
                "RP_M[8] ~ HP_M[7]", 
                "RP_M[8] - RP_M[6]", 
                "RP_M[6] ~ HP_M[6]", 
            ],
        },
        "镜像与否以及镜像轴": (False, None),
        "旋转拷贝的个数": slot_count,
    }

def draw_instruction_parser(part_dict, all_points, drawer):
    """
    将伪代码指令解析为具体绘图接口可接受的点段集合。
    兼容有限元工具 (例如 JMAG的 prepareSection方法) 和本地图形输出 (CairoDrawer)。
    输出一个规范化的 `region_dict` 包含各种描点。
    """
    type_name = part_dict.get('type')
    
    # 按照 part 的种类获取各自具体的 draw_instruction ()
    if type_name == "statorCore":
        result = get_stator_core_instruction(part_dict, all_points)
    elif type_name == "rotorCore":
        result = get_rotor_core_instruction(part_dict, all_points)
    elif type_name == "magnet":
        result = get_magnet_instruction(part_dict, all_points)
    elif type_name == "coil":
        result = get_coil_instruction(part_dict, all_points)
    else:
        raise ValueError(f"Unknown component type: {type_name}")
    
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
        scalar_instrs = result.get("几何绘制字符串", []) if isinstance(result, dict) else result
        instructions_list = [scalar_instrs]
        region_names = ["region1"]

    list_coords = []
    region_centroids = {}
    
    drawer.getSketch(part_dict['name'], part_dict['color'])

    part_mirror = part_dict.get('bMirror', None)
    part_rotate = part_dict.get('iRotateCopy', None)
    
    if part_mirror is not None: bMirror = part_mirror
    if part_rotate is not None: copyCount = part_rotate

    # 如果Drawer是有限元调用FEA，则往往只需要把这1个part丢给软件配置复制属性即可
    # 如果Drawer是本地SVG/Cairo渲染引擎，需要计算真实的旋转将其画完
    is_fea = hasattr(drawer, 'prepareSection')
    loop_count = 1 if (is_fea or copyCount <= 1) else copyCount
    
    for i in range(loop_count):
        # 计算整机圆形阵列每一份对应的旋转偏移
        rotation_offset = i * (360.0 / copyCount)
        current_deg = part_dict['rotation_deg'] + rotation_offset
        
        # 依次解析直线和弧线请求指令并调用底层API
        for idx, instructions in enumerate(instructions_list):
            region_segments = []
            region_coords = []
            for instr in instructions:
                instr = instr.split('#')[0].strip()
                if not instr: continue
                
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
            
            # 记录第一个Region段落的形心，FEA常需要区域内部质心坐标来识别材料属性域
            if i == 0 and region_coords:
                unique_reg_coords = list(set([tuple(p) for p in region_coords]))
                rx = sum(p[0] for p in unique_reg_coords) / len(unique_reg_coords)
                ry = sum(p[1] for p in unique_reg_coords) / len(unique_reg_coords)
                region_centroids[region_names[idx]] = (rx, ry)
    
    if not list_coords:
        return {'innerCoord': (0,0), 'inner_coords': {}, 'list_regions': [[]], 'mirrorAxis': mirrorAxis, 'bMirror': bMirror, 'iRotateCopy': copyCount}
        
    return {
        'inner_coords': region_centroids,
        'list_regions': list_regions,
        'mirrorAxis': mirrorAxis,
        'bMirror': bMirror,
        'iRotateCopy': copyCount
    }

class MachineGeometry:
    """
    负责统管生成、计算与储存机器截面的所有核心部件集合及关键参数点位。
    """
    def __init__(self):
        self.parts = []
        self._next_index = 0
        self.all_points = None

    def add_part(self, part_dict: OrderedDict):
        part_dict['index'] = self._next_index
        self.parts.append(part_dict)
        self._next_index += 1
        return part_dict
    
    def sync(self, machine_dict: OrderedDict):
        self.all_points = calculate_all_points(machine_dict)

    def draw_machine_using_CairoDrawer(self, drawer):
        """
        调用底层渲染引擎绘制所有的部件 (parts)。本框架利用 CairoDrawer 提供跨域统一绘图逻辑。
        """
        for part_dict in self.parts:
            # print(f"Drawing {part_dict['name']}...")
            res = draw_instruction_parser(part_dict, self.all_points, drawer)
            drawer.prepareSection(res, color=part_dict['color'])
            
            if hasattr(drawer, '_close_sketch_and_add_regions_to_list'):
                drawer._close_sketch_and_add_regions_to_list(
                    res.get('inner_coords', {}),
                    res.get('bMirror', False),
                    res.get('mirrorAxis', None),
                    res.get('iRotateCopy', 1),
                    len(res.get('list_regions', []))
                )

    def show_geometry_svg(self, filename='machine_geometry.svg', scale=30.0):
        # 以SVG格式导出绘图用于Web或其他文档系统的可视化
        drawer = CairoDrawer(filename=filename, scale=scale, bFillRegion=False)
        self.draw_machine_using_CairoDrawer(drawer)
        drawer.surface.finish()
        print(f"Geometry drawn to {filename} with scale {scale}")
