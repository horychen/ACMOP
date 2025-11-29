from dataclasses import dataclass, fields
from typing import Dict, List, Optional, Any
from collections import OrderedDict
import json, math, base64, pickle, cairo, os, jsonpickle

class Parameter(object):
    def __init__(self, name, type, value=None, bounds=None, calc=None, calc_bounds=None, unit='mm', comment=None, args=None) -> None:
        self.name = name
        self.type = type
        self.value = value
        self.bounds = bounds
        self.unit = unit
        self.calc = calc
        self.calc_bounds = calc_bounds
        self.comment = comment
        self.args = args
        # todo: add validation for the type, value, bounds, calc, unit, comment
        if self.calc is not None:
            if self.args is not None:
                self.value = self.calc(*self.args)
            else:
                self.value = self.calc()

        if self.calc_bounds is not None:
            if self.args is not None:
                self.bounds = self.calc_bounds(*self.args)
            else:
                self.bounds = self.calc_bounds()

        if type == 'free' and value is None and self.bounds is not None:
            if isinstance(self.bounds, (list, tuple)) and len(self.bounds) == 2:
                if self.bounds[0] is not None and self.bounds[1] is not None:
                    self.value = self.bounds[0] + (self.bounds[1] - self.bounds[0]) * 0.5

    def __repr__(self):
        return f"Parameter(name='{self.name}', type='{self.type}', value={self.value}, unit='{self.unit}')"

    def sensitivity(self, param_name: str) -> float:
        return 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """
        将 Parameter 对象转换为字典（用于 JSON 序列化）
        
        Returns:
            Dict: 参数字典
        """
        return {
            'name': self.name,
            'type': self.type,
            'value': self.value,
            'bounds': self.bounds,
            'unit': self.unit,
            'comment': self.comment,
            # calc 函数无法序列化，保存为 None，前端可以设置 calc_dependencies
            'calc_dependencies': None,  # 可以扩展为保存依赖的参数名列表
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Parameter':
        """
        从字典创建 Parameter 对象（用于 JSON 反序列化）
        
        Args:
            data: 参数字典
            
        Returns:
            Parameter: 参数对象
        """
        return cls(
            name=data['name'],
            type=data['type'],
            value=data.get('value'),
            bounds=data.get('bounds'),
            unit=data.get('unit', 'mm'),
            comment=data.get('comment'),
            calc=None,  # calc 函数需要从其他地方重建
            calc_bounds=None,  # calc_bounds 函数需要从其他地方重建
            args=data.get('args')  # 保存 args，用于重建 lambda 函数
        )

class Winding(object):
    def __init__(self, phase_number_m: int, stator_slot_number_Qs: int, pole_pair_number_p: int, suspension_pole_pair_number_ps: int, coil_pitch_y :int, bool_DPNVorSEPA: bool=True, number_of_parallel_branch: int=2) -> None:
        self.m = phase_number_m
        self.Qs = stator_slot_number_Qs
        self.p = pole_pair_number_p
        self.ps = suspension_pole_pair_number_ps
        self.SPP = self.Qs / (2*self.p * self.m)

        'Pay attention to the coil setup. The codes are written assuming the first coil is in the 12th slot. In other words PCoil[1] should negative, but {PCoil[1]=}'
        self.deg_winding_U_phase_phase_axis_angle = 0.0

        self.number_of_parallel_branch = number_of_parallel_branch
        self.number_of_winding_layer   = 2
        if self.number_of_winding_layer == 2:
            self.coil_pitch_y = coil_pitch_y
            self.bool_distributed_or_concentrated: bool = False if abs(coil_pitch_y) == 1 else True
        else:
            self.coil_pitch_y = self.Qs/self.p/2.0
            self.bool_distributed_or_concentrated: bool = True

        # TODO: calculate winding factor
        self.kw1 = 0.933
        derivation = self.get_winding_factor()
        # Get all winding factor information and phase/sign/grouping data from self.derivation
        # If self.derivation is a dict or an object with these attributes, extract them.
        # Fallback to defaults if not present.

        if self.derivation is not None:
            # Try both dict and object attribute access
            derivation = self.derivation
            # Winding factor information
            if hasattr(derivation, 'kw1'):
                self.kw1 = derivation.kw1
            elif isinstance(derivation, dict) and 'kw1' in derivation:
                self.kw1 = derivation['kw1']
            
            # Layer X phases
            if hasattr(derivation, 'layer_X_phases'):
                self.layer_X_phases = derivation.layer_X_phases
            elif isinstance(derivation, dict) and 'layer_X_phases' in derivation:
                self.layer_X_phases = derivation['layer_X_phases']
            
            # Layer X signs
            if hasattr(derivation, 'layer_X_signs'):
                self.layer_X_signs = derivation.layer_X_signs
            elif isinstance(derivation, dict) and 'layer_X_signs' in derivation:
                self.layer_X_signs = derivation['layer_X_signs']

            # Layer Y phases
            if hasattr(derivation, 'layer_Y_phases'):
                self.layer_Y_phases = derivation.layer_Y_phases
            elif isinstance(derivation, dict) and 'layer_Y_phases' in derivation:
                self.layer_Y_phases = derivation['layer_Y_phases']

            # Layer Y signs
            if hasattr(derivation, 'layer_Y_signs'):
                self.layer_Y_signs = derivation.layer_Y_signs
            elif isinstance(derivation, dict) and 'layer_Y_signs' in derivation:
                self.layer_Y_signs = derivation['layer_Y_signs']

            # grouping_AC
            if hasattr(derivation, 'grouping_AC'):
                self.grouping_AC = derivation.grouping_AC
            elif isinstance(derivation, dict) and 'grouping_AC' in derivation:
                self.grouping_AC = derivation['grouping_AC']


        # Excitation for DPNV
        self.bool_DPNVorSEPA = bool_DPNVorSEPA
        if self.bool_DPNVorSEPA == True:

            self.bool_3PhaseCurrentSource = False
            self.bool_CustomizedCircuit = False

            # the first coil in layer Y is assigned to phase W then this code is correct.
            self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W']
            self.layer_X_signs  = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
            self.layer_Y_phases = self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]

            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

    def infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self, layer_X_phases, coil_pitch):
        return layer_X_phases[-coil_pitch:] + layer_X_phases[:-coil_pitch]
    def infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self, layer_X_signs, coil_pitch):
        temp = layer_X_signs[-coil_pitch:] + layer_X_signs[:-coil_pitch]
        return [('-' if el == '+' else '+') for el in temp]

    def get_winding_factor(self):
        import winding_layout_derivation_ismb2021_asymetry_no_drawing
        self.derivation = winding_layout_derivation_ismb2021_asymetry_no_drawing.main_derivation(m=self.m, Qs=self.Qs, p=self.p, ps=self.ps, coil_pitch_y=self.coil_pitch_y)

    @staticmethod
    def draw_winding_in_the_slot(u, Qs, list_layer_phases, list_layer_signs, text=''):

        for i in range(Qs):
            radius_slot = 30
            LRIF = layer_radius_incremental_factor = 0.1
            angular_loc = 2*math.pi/Qs*i
            x_slot = radius_slot*math.cos(angular_loc)
            y_slot = radius_slot*math.sin(angular_loc)
            u.pyx_text(   [ x_slot*(1.0+LRIF), 
                            y_slot*(1.0+LRIF)],
                            str(i+1) )
            u.pyx_marker( [ x_slot*(1.0+2*LRIF), 
                            y_slot*(1.0+2*LRIF)], size=0.05)

            radius_tooth = radius_slot + 5 
            x_tooth = radius_tooth*math.cos(angular_loc + math.pi/Qs)
            y_tooth = radius_tooth*math.sin(angular_loc + math.pi/Qs)
            radius_airgap = radius_slot - 5
            x_toothtip = radius_airgap*math.cos(angular_loc+math.pi/Qs)
            y_toothtip = radius_airgap*math.sin(angular_loc+math.pi/Qs)
            u.pyx_line([x_toothtip, y_toothtip], [x_tooth, y_tooth])

            for ind, phases in enumerate(list_layer_phases):
                signs = list_layer_signs[ind]
                u.pyx_text(   [ x_slot*(1.0-ind*LRIF), 
                                y_slot*(1.0-ind*LRIF)],
                                '$' + phases[i].lower() + '^' + signs[i]
                                + '$' )

        u.pyx_text([0,0], (r'DPNV Winding' if wily.DPNV_or_SEPA else r'Separate Winding') + text)

        u = PyX_Utility.PyX_Utility()
        draw_winding_in_the_slot(u, wily.Qs, wily.list_layer_motor_phases, wily.list_layer_motor_signs, text=' Motor Mode' )
        u.cvs.writePDFfile(self.path2SwarmData + 'part_winding_pyx_output_M')

        u = PyX_Utility.PyX_Utility()
        draw_winding_in_the_slot(u, wily.Qs, wily.list_layer_suspension_phases, wily.list_layer_suspension_signs, text=' Suspension Mode' )
        u.cvs.writePDFfile(self.path2SwarmData + 'part_winding_pyx_output_S')
        # u.cvs.writeSVGfile(r'C:\Users\horyc\Desktop\pyx_output')
        # u.cvs.writeEPSfile(r'C:\Users\horyc\Desktop\pyx_output')
        # quit()

    @staticmethod
    def plot_winding_function(wily):
        '[1.2] Winding function / Current Linkage waveform'
        from pylab import plt, np
        zQ = 100 # number of conductors/turns per slot (Assume to be 100 for now)
        turns_per_layer = zQ / wily.number_winding_layer
        U_phase = winding_layout.PhaseWinding(wily.Qs, wily.m, turns_per_layer, wily.ox_distribution_phase_U)
        U_phase.plotFuncObj(U_phase.winding_func)
        U_phase.fig_plotFuncObj.savefig(self.path2SwarmData + 'part_winding_winding_function.png')
        U_phase.plot2piFft(U_phase.winding_func, Fs=1/(2*np.pi/3600), L=32000*2**4) # 在2pi的周期内取360个点
        U_phase.fig_plot2piFft.savefig(self.path2SwarmData + 'part_winding_winding_function_·.png')
        plt.show()

    def to_dict(self) -> Dict[str, Any]:
        """将 Winding 对象转换为字典"""
        return {
            'phase_number_m': self.m,
            'stator_slot_number_Qs': self.Qs,
            'pole_pair_number_p': self.p,
            'suspension_pole_pair_number_ps': self.ps,
            'coil_pitch_y': self.coil_pitch_y,
            'number_of_parallel_branch': self.number_of_parallel_branch,
            'kw1': self.kw1,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Winding':
        """从字典创建 Winding 对象"""
        # 创建对象实例
        instance = cls(
            phase_number_m=data.get('phase_number_m', data.get('m', 3)),
            stator_slot_number_Qs=data.get('stator_slot_number_Qs', data.get('Qs', 12)),
            pole_pair_number_p=data.get('pole_pair_number_p', data.get('p', 4)),
            suspension_pole_pair_number_ps=data.get('suspension_pole_pair_number_ps', data.get('ps', 5)),
            coil_pitch_y=data.get('coil_pitch_y', 1),
            number_of_parallel_branch=data.get('number_of_parallel_branch', 2)
        )
        
        # 设置其他属性（如果 JSON 中有的话）
        if 'kw1' in data:
            instance.kw1 = data['kw1']
        if 'bool_DPNVorSEPA' in data:
            instance.bool_DPNVorSEPA = data['bool_DPNVorSEPA']
        if 'layer_X_phases' in data:
            instance.layer_X_phases = data['layer_X_phases']
        if 'layer_X_signs' in data:
            instance.layer_X_signs = data['layer_X_signs']
        if 'grouping_AC' in data:
            instance.grouping_AC = data['grouping_AC']
        if 'number_parallel_branch' in data:
            instance.number_parallel_branch = data['number_parallel_branch']
        if 'number_winding_layer' in data:
            instance.number_winding_layer = data['number_winding_layer']
        if 'bool_3PhaseCurrentSource' in data:
            instance.bool_3PhaseCurrentSource = data['bool_3PhaseCurrentSource']
        if 'CommutatingSequenceD' in data:
            instance.CommutatingSequenceD = data['CommutatingSequenceD']
        if 'CommutatingSequenceB' in data:
            instance.CommutatingSequenceB = data['CommutatingSequenceB']
        
        # 重新计算 layer_Y_phases 和 layer_Y_signs（因为它们是通过方法计算的）
        instance.layer_Y_phases = instance.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(
            instance.layer_X_phases, instance.coil_pitch_y
        )
        instance.layer_Y_signs = instance.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(
            instance.layer_X_signs, instance.coil_pitch_y
        )
        
        return instance

class Geometry(object):
    def __init__(self, name, GP: dict = None, draw_function: callable = None, color: str = None):
        self.name = name
        self.color = color
        self.GP = GP or []
        self.draw_function = draw_function
        for i, (name, gp) in enumerate(self.GP.items()):
            if isinstance(gp, Parameter):
                exec(f"self.{name} = {gp.value}")
    def to_dict(self) -> Dict[str, Any]:
        # INSERT_YOUR_CODE
        """
        扩展序列化方法：将所有当前成员变量都存入字典（不局限于定义时的变量）
        排除不能序列化的 draw_function 属性
        """
        result = {}
        # 收集所有实例属性
        for k, v in self.__dict__.items():
            if k in ['draw_function']:
                # 无法序列化，略去
                continue
            elif k == 'GP':
                # 单独处理 GP
                if isinstance(v, list):
                    result['GP'] = [
                        gp.name if isinstance(gp, Parameter) else str(gp) if gp is not None else None
                        for gp in v
                    ]
                elif isinstance(v, dict):
                    result['GP'] = {kk: (vv.name if isinstance(vv, Parameter) else str(vv) if vv is not None else None)
                                    for kk, vv in v.items()}
                else:
                    result['GP'] = v
            elif k == 'visualization_points':
                # 特殊处理 visualization_points：确保它是可序列化的字典
                if isinstance(v, dict):
                    result['visualization_points'] = v
                else:
                    result['visualization_points'] = {}
            else:
                # 其它属性，直接存储基础类型，否则转为字符串
                if isinstance(v, (int, float, str, bool, type(None))):
                    result[k] = v
                else:
                    try:
                        # 尝试用 .to_dict()
                        result[k] = v.to_dict()
                    except Exception:
                        result[k] = str(v)
        # 必须保证 _needs_rebuild 存在
        result['_needs_rebuild'] = True
        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Geometry':
        return cls(
            color=data['color'],
            GP=data['GP'],
            draw_function=data['draw_function'],
        )
    def draw(self, drawer, *args, **kwargs):
        # 记录调用 draw 之前的 visualization_points 键，以便识别新添加的点
        old_keys = set(getattr(drawer, 'visualization_points', {}).keys())
        
        # 调用 draw_function
        self.components_make_region = self.draw_function(drawer, *args, **kwargs)
        # print(self.components_make_region)
        # print(drawer)
        # print() 

        # 从 drawer.visualization_points 中提取新添加的点坐标
        if hasattr(drawer, 'visualization_points'):
            new_keys = set(drawer.visualization_points.keys()) - old_keys
            # 如果有新添加的点，将它们存储到 self.visualization_points
            if new_keys:
                # 通常只有一个新键，但为了安全起见，我们存储所有新添加的点
                for key in new_keys:
                    if key in drawer.visualization_points:
                        self.visualization_points = drawer.visualization_points[key]
                        break  # 通常只需要第一个匹配的键
        
        return self.components_make_region

class CairoDrawer(object):
    def __init__(self, width_in_points=500, height_in_points=500, filename=None):
        self.filename = filename
        self.surface = cairo.SVGSurface(self.filename, width_in_points, height_in_points)
        self.ctx = cairo.Context(self.surface)
        # self.ctx.scale(width_in_points, height_in_points)
        # m = cairo.Matrix(yy=-1, y0=height_in_points) # Cartetian Coordinate
        m = cairo.Matrix(yy=-1, y0=0.5*height_in_points, x0=+0.5*width_in_points) # Offset to center
        self.ctx.transform(m)
        # Set a background color
        self.ctx.save()
        self.ctx.set_source_rgb(0.95, 0.95, 0.95)
        self.ctx.paint()
        self.ctx.restore()
    def apply_stroke(self, lw=0.5):
        self.ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        self.ctx.set_line_width(lw)
        # setting color of the context
        self.ctx.set_source_rgba(0.0, 0.0, 0.0, 1)
        # stroke out the color and width property
        self.ctx.stroke()
    def convert_to_pdf(self, bool_open_pdf=False, filename=None): # 这个代码只是把SVG转换为PDF而已
        # self.surface.write_to_svg()
        self.surface.finish()
        import cairosvg
        cairosvg.svg2pdf(url=f'machine_geometry.svg', write_to=f'machine_geometry.pdf')
        if bool_open_pdf:
            import os
            os.system('sumatraPDF2.exe ' + 'machine_geometry.pdf')
        print('[machine_design_guide.py] Find the file machine_geometry.pdf in the current folder.')

    def getSketch(self, name, color):
        self.sketch_name = name
        self.sketch_color = color

    def drawLine(self, p1, p2):
        self.ctx.move_to(p1[0], p1[1])
        self.ctx.line_to(p2[0], p2[1])
        return [{'move_to': (p1[0], p1[1]), 'line_to': (p2[0], p2[1])}]

    def drawArc(self, centerxy, startxy, endxy):
        EPS = 1e-3
        v1 = [startxy[0] - centerxy[0], startxy[1] - centerxy[1]]
        v2 = [endxy[0]   - centerxy[0], endxy[1]   - centerxy[1]]
        cos夹角 = (v1[0]*v2[0] + v1[1]*v2[1]) / (math.sqrt(v1[0]*v1[0] + v1[1]*v1[1])*math.sqrt(v2[0]*v2[0] + v2[1]*v2[1]))
        if 1.0 < cos夹角 < 1.0+EPS:
            cos夹角 = 1.0
        elif -1.0-EPS < cos夹角 < -1.0:
            cos夹角 = -1.0
        angle_between = math.acos(cos夹角)

        radius = math.sqrt(v1[0]*v1[0] + v1[1]*v1[1])
        angle_start = math.atan2(v1[1], v1[0])
        angle_end = angle_start + angle_between

        self.ctx.move_to(startxy[0], startxy[1])
        self.ctx.arc(centerxy[0], centerxy[1], radius, angle_start, angle_end)
        # self.ctx.arc_negative(centerxy[0], centerxy[1], radius, angle_end, angle_start)
        return [{'move_to': (centerxy[0], centerxy[1]), 'arc': (radius, angle_start, angle_end)}]

@dataclass
class Modern_Machine_Designer(object):

    # Meta Data
    name: str = 'SPMSM'
    machine_class: str = 'bearingless_spmsm_heart.bearingless_spmsm_design_variant'

    # Machine Geometry
    bool_PermanentMagnet: bool = True
    bool_StatorSlotClosed: bool = False
    bool_RotorNotched: bool = True

    select_FEA_tool: str = 'JMAG Designer' # FEMM
    select_fea_config_dict: str = '#0213 JMAG Bearingless Sub-hamonics'
    fea_config_dict: dict = None
    bool_jmagDeleteResultsAfterCalculation: bool = False

    # Optimization
    counter: int = 0
    counter_fitness_called: int = 0
    counter_fitness_return: int = 0

    def __post_init__(self):

        # 绕组
        m : int = 3
        Qs: int = 12
        p : int = 4
        ps: int = 5
        coil_pitch_y: int = 1


        self.wily = Winding(m, Qs, p, ps, coil_pitch_y, bool_DPNVorSEPA=True)

        # 激励（含热负荷）
        bool_WyeConnectOrDeltaConnect: bool = True
        bool_weHavePlentyVoltage: bool = True

        Temperature : float = 75
        available_temperature_list = [-40, 20, 60, 80, 100, 120, 150, 180, 200, 220] # according to JMAG
        Magnet_Temperature = min(available_temperature_list, key=lambda x:abs(x - Temperature))

        RatedPower: float = 50e3 # W
        RatedSpeed: float = 30000 # rpm
        ExcitationFreqSimulated: float = RatedSpeed / 60 * p

        TORQUE_CURRENT_RATIO: float = 0.95
        SUSPENSION_CURRENT_RATIO: float = 0.05

        SteelMaterial = "M-19 Steel Gauge-29"

        self.EX = EX = {
            # 3D
            'mm_stack_length_specified': 50, # mm
            # Materials
            'Magnet_Name': u"Arnold/Reversible/N40H",
            'Magnet_StartAngle': 0.5* 360/(2*p),
            'Magnet_Temperature': Magnet_Temperature,
            'SteelMaterial': SteelMaterial,
            'StatorCore_Material': SteelMaterial, # "M-15 Steel", "Arnon5-final", u"35CS250", "DCMagnetic Type/50A1000",
            'RotorCore_Material': SteelMaterial, # "M-15 Steel", "Arnon5-final", u"35CS250", "DCMagnetic Type/50A1000",
            'LaminationFactor': 95,
            # Thermal
            'RatedPower': RatedPower,
            'RatedSpeed': RatedSpeed,
            'ExcitationFreqSimulated': ExcitationFreqSimulated,
            'bool_WyeConnectOrDeltaConnect' : bool_WyeConnectOrDeltaConnect,
            'DCBusVoltage': 600,
            'Js' : 4e6,
            'Temperature': 75,
            'WindingFill': 0.3882,
            'TORQUE_CURRENT_RATIO': TORQUE_CURRENT_RATIO,
            'SUSPENSION_CURRENT_RATIO': SUSPENSION_CURRENT_RATIO,
            'DriveW_Rs': 1.0, # [Ohm]
            'BeariW_Rs': 1.0, # [Ohm]
        }

        # 定子裂比和外径
        SR: float = 0.35
        mm_r_so: float = 123.5

        # 利用不同的裂比去估算合理的边界值
        yoke_split_ratio_bounds: list[float] = [0.2, 0.45]
        tooth_split_ratio_at_middle_slot: list[float] = [0.25, 0.50] # TODO: 下界需要考虑到w_st的宽度和半径的值

        '''Fixed variables'''
        self.m: Parameter            = Parameter('phase_number_m', 'fixed', m)
        self.Qs: Parameter           = Parameter('stator_slot_number_Qs', 'fixed', Qs)
        self.p: Parameter            = Parameter('pole_pair_number_p', 'fixed', p)
        self.ps: Parameter           = Parameter('suspension_pole_pair_number_ps', 'fixed', ps)
        self.coil_pitch_y: Parameter = Parameter('coil_pitch_y', 'fixed', coil_pitch_y)
        self.mm_r_so: Parameter      = Parameter('stator_outer_radius', 'fixed', mm_r_so)
        self.mm_d_mech_air_gap: Parameter = Parameter('mechanical_air_gap_depth', 'fixed', 0.5)

        if self.bool_PermanentMagnet:
            self.mm_d_pm: Parameter      = Parameter('magnet_depth', 'free', 3, bounds=[2, 6])
            self.mm_d_sleeve: Parameter  = Parameter('rotor_sleeve_depth', 'fixed', 1.0)
            self.s: Parameter            = Parameter('number_of_magnet_segments_per_pole', 'fixed', 1)

        '''Free variables'''
        self.split_ratio: Parameter  = Parameter('split_ratio_r_si_slash_r_so', 'free', SR, calc_bounds=lambda p: [0.2, 0.5] if p < 10 else [0.15, 0.35], args=[p])
            # "split_ratio":  [0.4, 0.6], # Binder-2020-MLMS-0953@Fig.7
            # "split_ratio":  [0.35, 0.5], # Q12p4优化的时候，轭部经常不够用，所以就把split_ratio减小——Exception: ('Error: Negative derived parameter', "acmop_parameter(type='derived', name='stator_yoke_depth', value=-1.362043443071423, bounds=[None, None], calc=<function template_machine_as_numbers.__init__.<locals>.<lambda> at 0x00000237CC403D30>)")

        # 计算 mm_r_si 的初始值用于 calc_bounds
        mm_r_si_initial = mm_r_so * SR
        self.mm_w_st: Parameter = Parameter('stator_tooth_width', 'free', calc_bounds=lambda mm_r_so, mm_r_si, tooth_split_ratio_at_middle_slot, Q: [el / Q * math.pi * (mm_r_so + mm_r_si) for el in tooth_split_ratio_at_middle_slot], args=[mm_r_so, mm_r_si_initial, tooth_split_ratio_at_middle_slot, Qs])
        self.mm_d_sy: Parameter = Parameter('stator_yoke_depth', 'free', calc_bounds=lambda mm_r_so, mm_r_si, yoke_split_ratio_bounds: [el * (mm_r_so - mm_r_si) for el in yoke_split_ratio_bounds], args=[mm_r_so, mm_r_si_initial, yoke_split_ratio_bounds])
        self.mm_d_sts: Parameter = Parameter('stator_tooth_shoe_depth', 'free', bounds=[1, 5])

        if not self.bool_StatorSlotClosed:
            self.deg_alpha_st: Parameter = Parameter('stator_tooth_span_angle', 'free', calc_bounds=lambda Qs: [360/Qs*0.1, 360/Qs], unit='deg', args=[self.Qs.value])

        '''derived variables have dependency on the other geometric parameters'''
        self.mm_r_si: Parameter      = Parameter('stator_inner_radius', 'derived', calc=lambda mm_r_so, split_ratio: mm_r_so * split_ratio, args=[self.mm_r_so.value, self.split_ratio.value])
        self.mm_d_st: Parameter      = Parameter('stator_tooth_depth', 'derived', calc=lambda mm_r_so, mm_r_si, mm_d_sy, mm_d_sts: mm_r_so - mm_r_si - mm_d_sy - mm_d_sts, args=[self.mm_r_so.value, self.mm_r_si.value, self.mm_d_sy.value, self.mm_d_sts.value])
        if self.bool_PermanentMagnet:
            self.mm_d_ri: Parameter      = Parameter('rotor_iron (back iron) depth', 'derived', calc=lambda mm_d_pm: 4 if mm_d_pm < 4 else mm_d_pm, args=[self.mm_d_pm.value])
            self.mm_r_ro: Parameter      = Parameter('rotor_outer_radius', 'derived', calc=lambda mm_r_si, mm_d_mech_air_gap, mm_d_sleeve: mm_r_si - mm_d_mech_air_gap - mm_d_sleeve, args=[self.mm_r_si.value, self.mm_d_mech_air_gap.value, self.mm_d_sleeve.value])
            self.mm_r_ri: Parameter      = Parameter('rotor_inner_radius', 'derived', calc=lambda r_ro, mm_d_pm, mm_d_ri: r_ro-mm_d_pm-mm_d_ri, args=[self.mm_r_ro.value, self.mm_d_pm.value, self.mm_d_ri.value])

        if not self.bool_StatorSlotClosed:
            self.mm_d_sto: Parameter = Parameter('stator_tooth_open_depth', 'derived', calc=lambda mm_d_sts: mm_d_sts*0.667, args=[self.mm_d_sts.value])
            # deg_alpha_sto 依赖于 deg_alpha_st，使用 deg_alpha_st 的当前值（如果已计算）或使用 bounds 的中间值
            deg_alpha_st_value = self.deg_alpha_st.value if self.deg_alpha_st.value is not None else (self.deg_alpha_st.bounds[0] + self.deg_alpha_st.bounds[1]) / 2 if self.deg_alpha_st.bounds else 360/12*0.1*0.5
            self.deg_alpha_sto: Parameter = Parameter('stator_tooth_open_angle', 'derived', calc=lambda deg_alpha_st: deg_alpha_st*0.5, args=[deg_alpha_st_value])

        if self.bool_RotorNotched:
            self.deg_alpha_rm: Parameter = Parameter('magnet_pole_span_angle', 'free', bounds=[180/self.p.value*0.7, 180/self.p.value])

            # deg_alpha_rs 依赖于 deg_alpha_rm，使用 deg_alpha_rm 的当前值（如果已计算）或使用 bounds 的中间值
            deg_alpha_rm_value = self.deg_alpha_rm.value if self.deg_alpha_rm.value is not None else (self.deg_alpha_rm.bounds[0] + self.deg_alpha_rm.bounds[1]) / 2 if self.deg_alpha_rm.bounds else 360/12*0.1
            self.deg_alpha_rs: Parameter = Parameter('magnet_segment_span_angle', 'derived', calc=lambda deg_alpha_rm: deg_alpha_rm, args=[deg_alpha_rm_value])

            self.mm_d_rp: Parameter      = Parameter('inter_polar_iron_thickness', 'derived', calc=lambda mm_d_pm: mm_d_pm, args=[self.mm_d_pm.value])
            self.mm_d_rs: Parameter      = Parameter('inter_segment_iron_thickness', 'fixed', 0.0)


        V_stator_phase_voltage_amp = math.sqrt(2) *EX['DCBusVoltage'] / (math.sqrt(3) if bool_WyeConnectOrDeltaConnect else 1.0) 
        V_desired_emf_Em = 0.95 * V_stator_phase_voltage_amp
        alpha_i = 2.0/math.pi # ideal sinusoidal flux density distribusion, when the saturation happens in teeth, alpha_i becomes higher.
        T_air_gap_flux_density_Bg_guessed = 0.9 # T
        mm_stack_length_specified = EX['mm_stack_length_specified']
        mm_d_magnetic_air_gap = self.mm_d_mech_air_gap.value + self.mm_d_sleeve.value
        mm_stack_length_effective = mm_stack_length_specified + 2 * mm_d_magnetic_air_gap
        mm_pole_pitch_tau_p = math.pi *self.mm_r_si.value / p
        Wb_air_gap_flux_Phi_m = alpha_i * T_air_gap_flux_density_Bg_guessed * mm_pole_pitch_tau_p*1e-3 * mm_stack_length_effective*1e-3 # Wb
        no_series_coil_turns_N = V_desired_emf_Em / (2*math.pi* ExcitationFreqSimulated * self.wily.kw1 * Wb_air_gap_flux_Phi_m)
        no_series_coil_turns_N = round(no_series_coil_turns_N)
        SPP = Qs / (2*p*m) # slot per pole per phase
        print(f"[DEBUG] m={m}")
        print(f"[DEBUG] Qs={Qs}")
        print(f"[DEBUG] p={p}")
        print(f"[DEBUG] ps={ps}")
        print(f"[DEBUG] coil_pitch_y={coil_pitch_y}")
        print(f"[DEBUG] V_stator_phase_voltage_amp={V_stator_phase_voltage_amp}")
        print(f"[DEBUG] V_desired_emf_Em={V_desired_emf_Em}")
        print(f"[DEBUG] alpha_i={alpha_i}")
        print(f"[DEBUG] T_air_gap_flux_density_Bg_guessed={T_air_gap_flux_density_Bg_guessed}")
        print(f"[DEBUG] mm_stack_length_specified={mm_stack_length_specified}")
        print(f"[DEBUG] mm_d_magnetic_air_gap={mm_d_magnetic_air_gap}")
        print(f"[DEBUG] mm_stack_length_effective={mm_stack_length_effective}")
        print(f"[DEBUG] mm_pole_pitch_tau_p={mm_pole_pitch_tau_p}")
        print(f"[DEBUG] Wb_air_gap_flux_Phi_m={Wb_air_gap_flux_Phi_m}")
        print(f"[DEBUG] no_series_coil_turns_N={no_series_coil_turns_N}")
        print(f"[DEBUG] SPP={SPP}")
        if bool_weHavePlentyVoltage:
            no_series_coil_turns_N = min([p*SPP*i for i in range(1000,0,-1)], key=lambda x:abs(x - no_series_coil_turns_N)) # using larger turns value has priority
        else:
            no_series_coil_turns_N = min([p*SPP*i for i in range(1000)], key=lambda x:abs(x - no_series_coil_turns_N))  # using lower turns value has priority # https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
        if no_series_coil_turns_N > 990:
            raise Exception(f'What? no_series_coil_turns_N is too large: {no_series_coil_turns_N=}')
        # print(f'[zQ] We need {no_series_coil_turns_N=} to reach the desired voltage: {V_desired_emf_Em=} V when {EX["DCBusVoltage"]=} V')
        # print(f'[zQ] {no_series_coil_turns_N=} should be multiple of pq: q * p = {SPP} * {p}')
        EX['no_series_coil_turns_N'] = no_series_coil_turns_N
        EX['DriveW_zQ'] = no_conductors_per_slot_zQ = 2* m * no_series_coil_turns_N / Qs * self.wily.number_of_parallel_branch
        EX['BeariW_zQ'] = EX['DriveW_zQ'] if self.wily.bool_DPNVorSEPA == True else EX['DriveW_zQ'] / EX['TORQUE_CURRENT_RATIO'] * EX['SUSPENSION_CURRENT_RATIO']

        ''' Excitations Consiering Thermal Capability Limit (Simple) '''
        mm_r_sy = self.mm_r_so.value - self.mm_d_sy.value  # radius stator yoke
        mm_r_ss = self.mm_r_si.value + self.mm_d_sts.value # radius stator slot
        EX['mm2_slot_area']            = (math.pi*(mm_r_sy**2 - mm_r_ss**2) / Qs - self.mm_w_st.value * self.mm_d_st.value) # 计算槽面积
        EX['CurrentAmp_in_the_slot']   = EX['mm2_slot_area'] * 1e-6 * EX['Js'] * EX['WindingFill'] * math.sqrt(2)
        EX['CurrentAmp_per_conductor'] = EX['CurrentAmp_in_the_slot'] / EX['DriveW_zQ']
        EX['CurrentAmp_per_phase']     = EX['CurrentAmp_per_conductor'] * self.wily.number_of_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
        EX['DriveW_CurrentAmp'] = EX['TORQUE_CURRENT_RATIO']     * EX['CurrentAmp_per_phase']
        EX['BeariW_CurrentAmp'] = EX['SUSPENSION_CURRENT_RATIO'] * EX['CurrentAmp_per_phase']
        EX['slot_current_utilizing_ratio_for_torque'] = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']

        EX['InitialRotationAngle'] = self.get_InitialRotationAngle()

        Rout = self.mm_r_ri.value+ self.mm_d_ri.value+ self.mm_d_pm.value
        Rin  = self.mm_r_ri.value+ self.mm_d_ri.value
        deg_alpha_rp = 360 / (2*self.p.value)
        EX['mm2_magnet_area'] = self.deg_alpha_rm.value/deg_alpha_rp * math.pi*(Rout**2 - Rin**2)

        # INSERT_YOUR_CODE
        print(f"[DEBUG] mm_r_sy={mm_r_sy}")
        print(f"[DEBUG] mm_r_ss={mm_r_ss}")
        print(f"[DEBUG] EX['mm2_slot_area']={EX['mm2_slot_area']}")
        print(f"[DEBUG] EX['CurrentAmp_in_the_slot']={EX['CurrentAmp_in_the_slot']}")
        print(f"[DEBUG] EX['CurrentAmp_per_conductor']={EX['CurrentAmp_per_conductor']}")
        print(f"[DEBUG] EX['CurrentAmp_per_phase']={EX['CurrentAmp_per_phase']}")
        print(f"[DEBUG] EX['DriveW_CurrentAmp']={EX['DriveW_CurrentAmp']}")
        print(f"[DEBUG] EX['BeariW_CurrentAmp']={EX['BeariW_CurrentAmp']}")
        print(f"[DEBUG] EX['slot_current_utilizing_ratio_for_torque']={EX['slot_current_utilizing_ratio_for_torque']}")
        print(f"[DEBUG] EX['InitialRotationAngle']={EX['InitialRotationAngle']}")
        print(f"[DEBUG] Rout={Rout}")
        print(f"[DEBUG] Rin={Rin}")
        print(f"[DEBUG] deg_alpha_rp={deg_alpha_rp}")
        print(f"[DEBUG] EX['mm2_magnet_area']={EX['mm2_magnet_area']}")

        # raise KeyboardInterrupt

        import CrossSectInnerNotchedRotor, CrossSectStator
        self.machineGeometry = {
            "rotorCore": Geometry(name='rotorCore',
                GP={
                    'mm_r_ro': self.mm_r_ro,
                    'mm_d_ri': self.mm_d_ri,
                    'mm_d_pm': self.mm_d_pm,
                    'mm_d_rp': self.mm_d_rp,
                    'mm_d_rs': self.mm_d_rs,
                    'p': self.p,
                    's': self.s,

                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                        name="rotorCore",
                        color="#FE840E",
                        mm_d_pm=self.mm_d_pm.value if hasattr(self, "mm_d_pm") and self.mm_d_pm.value is not None else 6,
                        deg_alpha_rm=self.deg_alpha_rm.value if hasattr(self, "deg_alpha_rm") and self.deg_alpha_rm.value is not None else 60,
                        deg_alpha_rs=self.deg_alpha_rs.value if hasattr(self, "deg_alpha_rs") and self.deg_alpha_rs.value is not None else 10,
                        mm_d_ri=self.mm_d_ri.value if hasattr(self, "mm_d_ri") and self.mm_d_ri.value is not None else 8,
                        mm_r_ri=self.mm_r_ri.value if hasattr(self, "mm_r_ri") and self.mm_r_ri.value is not None else 40,
                        mm_d_rp=self.mm_d_rp.value if hasattr(self, "mm_d_rp") and self.mm_d_rp.value is not None else 5,
                        mm_d_rs=self.mm_d_rs.value if hasattr(self, "mm_d_rs") and self.mm_d_rs.value is not None else 3,
                        p=self.p.value if hasattr(self, "p") and self.p.value is not None else 2,
                        s=self.s.value if hasattr(self, "s") and self.s.value is not None else 4
                    ).draw(drawer, **kwargs)
                )
            ),
            "shaft": Geometry(name='shaft',
                GP={'mm_r_ri': self.mm_r_ri},
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectShaft(
                        name="shaft",
                        color="#0EE0E2",
                        rotorCore=CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                            mm_d_pm=self.mm_d_pm.value if hasattr(self, "mm_d_pm") and self.mm_d_pm.value is not None else 6,
                            deg_alpha_rm=self.deg_alpha_rm.value if hasattr(self, "deg_alpha_rm") and self.deg_alpha_rm.value is not None else 60,
                            deg_alpha_rs=self.deg_alpha_rs.value if hasattr(self, "deg_alpha_rs") and self.deg_alpha_rs.value is not None else 10,
                            mm_d_ri=self.mm_d_ri.value if hasattr(self, "mm_d_ri") and self.mm_d_ri.value is not None else 8,
                            mm_r_ri=self.mm_r_ri.value if hasattr(self, "mm_r_ri") and self.mm_r_ri.value is not None else 40,
                            mm_d_rp=self.mm_d_rp.value if hasattr(self, "mm_d_rp") and self.mm_d_rp.value is not None else 5,
                            mm_d_rs=self.mm_d_rs.value if hasattr(self, "mm_d_rs") and self.mm_d_rs.value is not None else 3,
                            p=self.p.value if hasattr(self, "p") and self.p.value is not None else 2,
                            s=self.s.value if hasattr(self, "s") and self.s.value is not None else 4
                        )
                    ).draw(drawer, **kwargs)
                )
            ),
            "rotorMagnet": Geometry(name='rotorMagnet',
                GP={
                    'mm_d_pm': self.mm_d_pm,
                    'mm_d_ri': self.mm_d_ri,
                    'mm_r_ri': self.mm_r_ri,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet(
                        name="rotorMagnet",
                        color="#1C96E0",
                        rotorCore=CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                            mm_d_pm=self.mm_d_pm.value if hasattr(self, "mm_d_pm") and self.mm_d_pm.value is not None else 6,
                            deg_alpha_rm=self.deg_alpha_rm.value if hasattr(self, "deg_alpha_rm") and self.deg_alpha_rm.value is not None else 60,
                            deg_alpha_rs=self.deg_alpha_rs.value if hasattr(self, "deg_alpha_rs") and self.deg_alpha_rs.value is not None else 10,
                            mm_d_ri=self.mm_d_ri.value if hasattr(self, "mm_d_ri") and self.mm_d_ri.value is not None else 8,
                            mm_r_ri=self.mm_r_ri.value if hasattr(self, "mm_r_ri") and self.mm_r_ri.value is not None else 40,
                            mm_d_rp=self.mm_d_rp.value if hasattr(self, "mm_d_rp") and self.mm_d_rp.value is not None else 5,
                            mm_d_rs=self.mm_d_rs.value if hasattr(self, "mm_d_rs") and self.mm_d_rs.value is not None else 3,
                            p=self.p.value if hasattr(self, "p") and self.p.value is not None else 2,
                            s=self.s.value if hasattr(self, "s") and self.s.value is not None else 4
                        )
                    ).draw(drawer, **kwargs)
                ),
            ),
            "sleeve": Geometry(name='sleeve',
                GP={
                    'mm_r_ri': self.mm_r_ri,
                    'mm_d_ri': self.mm_d_ri,
                    'mm_d_pm': self.mm_d_pm,
                    'p': self.p,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectSleeve(
                        mm_r_ri=self.mm_r_ri.value if hasattr(self, "mm_r_ri") and self.mm_r_ri.value is not None else 5,
                        mm_d_ri=self.mm_d_ri.value if hasattr(self, "mm_d_ri") and self.mm_d_ri.value is not None else 5,
                        mm_d_pm=self.mm_d_pm.value if hasattr(self, "mm_d_pm") and self.mm_d_pm.value is not None else 3,
                        p=self.p.value if hasattr(self, "p") and self.p.value is not None else 4,
                        d_sleeve=self.d_sleeve.value if hasattr(self, "d_sleeve") and self.d_sleeve.value is not None else 1
                    ).draw(drawer, **kwargs)
                )
            ),
            "statorCore": Geometry(name='statorCore',
                GP={
                    'mm_r_si': self.mm_r_si,
                    'mm_d_sto': self.mm_d_sto,
                    'mm_d_sts': self.mm_d_sts,
                    'mm_d_st': self.mm_d_st,
                    'mm_d_sy': self.mm_d_sy,
                    'mm_w_st': self.mm_w_st,
                    'deg_alpha_st': self.deg_alpha_st,
                    'deg_alpha_sto': self.deg_alpha_sto,
                    'mm_r_si': self.mm_r_si,
                    'mm_d_sto': self.mm_d_sto,
                    'Q': self.Qs,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectStator.CrossSectInnerRotorStator(
                        name="statorCore",
                        color="#BAFA01",
                        deg_alpha_st=self.deg_alpha_st.value if hasattr(self, "deg_alpha_st") and self.deg_alpha_st.value is not None else 40,
                        deg_alpha_sto=self.deg_alpha_sto.value if hasattr(self, "deg_alpha_sto") and self.deg_alpha_sto.value is not None else 20,
                        mm_r_si=self.mm_r_si.value if hasattr(self, "mm_r_si") and self.mm_r_si.value is not None else 40,
                        mm_d_sto=self.mm_d_sto.value if hasattr(self, "mm_d_sto") and self.mm_d_sto.value is not None else 5,
                        mm_d_sts=self.mm_d_sts.value if hasattr(self, "mm_d_sts") and self.mm_d_sts.value is not None else 10,
                        mm_d_st=self.mm_d_st.value if hasattr(self, "mm_d_st") and self.mm_d_st.value is not None else 15,
                        mm_d_sy=self.mm_d_sy.value if hasattr(self, "mm_d_sy") and self.mm_d_sy.value is not None else 15,
                        mm_w_st=self.mm_w_st.value if hasattr(self, "mm_w_st") and self.mm_w_st.value is not None else 13,
                        Q=self.Q.value if hasattr(self, "Q") and self.Q.value is not None else 6,
                    ).draw(drawer, **kwargs)
                ),
            ),
            "coils": None
        }
        self.machineGeometry['coils'] = Geometry(name='coils',
                GP={
                    'mm_r_so': self.mm_r_so,
                    'mm_d_sy': self.mm_d_sy,
                    'mm_w_st': self.mm_w_st,
                    'mm_d_st': self.mm_d_st,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectStator.CrossSectInnerRotorStatorWinding(
                        stator_core=self.machineGeometry['statorCore'],
                    ).draw(drawer, **kwargs)
                )
            )


    def get_InitialRotationAngle(self):
        # 设置转子角度的初始位置条件，以使得在t=0时刻，转子的q轴与U相绕组的相轴重合，并且此时U相电流应该为交流最大（需要同步调整circuit中的激励正弦信号的相位）。
        # 不仅如此，初始转子角度还影响着永磁体的励磁角度是否对齐，最好手动确认一下： study.GetMaterial(u"Magnet").SetValue(u"StartAngle", 0.5* 360/(2*acm_variant.['p']) ) # 半个极距
        # Implementation of id=0 control:
        #   After rotate the rotor by half the inter-pole notch span, The d-axis initial position is at pole pitch angle divided by 2.
        #   The U-phase current is sin(omega_syn*t) = 0 at t=0 and requires the d-axis to be at the winding phase axis (to obtain id=0 control)
        deg_pole_span = 180/self.p.value
        #                           inter-pole notch is rotated to x-axis (0.5 for half)  winding placing bias (one slot angle)        align with q-axis
        self.InitialRotationAngle = (deg_pole_span-self.deg_alpha_rm.value)*0.5 + self.wily.deg_winding_U_phase_phase_axis_angle  # is made by set phase U current maximum at t=0, that is the current is a cosine function.
        print(f"[bearingless_spmsm_design.py] [PMSM JMAG] {self.InitialRotationAngle} deg = ", (deg_pole_span-self.deg_alpha_rm.value)*0.5,  self.wily.deg_winding_U_phase_phase_axis_angle,  deg_pole_span*0.5)
        print(f"[bearingless_spmsm_design.py] [PMSM JMAG] {self.InitialRotationAngle} deg")

        return self.InitialRotationAngle


    def show_geometry(self, filename=None) -> None:
        bool_draw_whole_model = True
        
        # 确保 machineGeometry 已初始化
        if not hasattr(self, 'machineGeometry') or self.machineGeometry is None:
            # 如果 machineGeometry 不存在，调用 __post_init__ 来创建
            self.__post_init__()

        def draw_spmsm(lw, width_in_points, height_in_points, filename='machine_geometry.svg', bool_draw_whole_model=True):
            self.drawer = drawer = CairoDrawer(width_in_points, height_in_points, filename=filename)

            # 检查 machineGeometry 是否存在且包含必要的键
            if not hasattr(self, 'machineGeometry') or self.machineGeometry is None:
                raise ValueError("machineGeometry is not initialized. Please ensure __post_init__ was called.")

            # 安全地调用 draw 方法
            if  'rotorCore' in self.machineGeometry and self.machineGeometry['rotorCore'] is not None\
            and 'shaft' in self.machineGeometry and self.machineGeometry['shaft'] is not None\
            and 'rotorMagnet' in self.machineGeometry and self.machineGeometry['rotorMagnet'] is not None\
            and 'statorCore' in self.machineGeometry and self.machineGeometry['statorCore'] is not None\
            and 'coils' in self.machineGeometry and self.machineGeometry['coils'] is not None:
                list_regions = self.machineGeometry['rotorCore'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)
                list_regions = self.machineGeometry['shaft'].draw(drawer)
                list_regions = self.machineGeometry['rotorMagnet'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)
                list_regions = self.machineGeometry['statorCore'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)
                list_regions = self.machineGeometry['coils'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)

            drawer.apply_stroke(lw=lw)
            drawer.convert_to_pdf()

            # import builtins
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['rotorCore'] = acm_variant.rotorCore
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['shaft'] = acm_variant.shaft
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['rotorMagnet'] = acm_variant.rotorMagnet
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['sleeve'] = acm_variant.sleeve
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['statorCore'] = acm_variant.statorCore
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['coils'] = acm_variant.coils


        # mm_r_ro 只在 bool_PermanentMagnet 为 True 时存在
        if hasattr(self, 'mm_r_ro') and self.mm_r_ro.value is not None:
            lw = 0.1 if self.mm_r_ro.value < 15 else 0.5
        else:
            # 使用默认值
            lw = 0.5
        # 检查必要的参数是否存在
        if not hasattr(self, 'mm_r_so') or self.mm_r_so.value is None:
            raise ValueError("mm_r_so parameter is required but not found or has no value")
        width_in_points  = self.mm_r_so.value*2.1
        height_in_points = self.mm_r_so.value*2.1
        draw_spmsm(lw, width_in_points, height_in_points)

    def FEA_evaluate(self, project_loc=fr'../_default/', bool_jmagDesignerShow: bool = True):

        if self.fea_config_dict is None:
            with open((os.path.dirname(__file__))+'/machine_simulation.json', 'r') as f:
                raw_fea_config_dicts = json.load(f)
                self.fea_config_dict = OrderedDict(raw_fea_config_dicts[self.select_fea_config_dict])

        ''' 工程和文件路径 '''
        def get_pc_name():
            import platform, socket
            n1 = platform.node()
            n2 = socket.gethostname()
            n3 = os.environ["COMPUTERNAME"]
            if n1 == n2 == n3:
                return n1
            else:
                raise Exception(f"Computer names are not equal to each other. {n1,n2,n3}")
        dir_parent = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/'
        dir_codes  = os.path.abspath(os.path.dirname(__file__)) + '/'
        pc_name = get_pc_name()
        # os.chdir(dir_codes)
        self.path2SwarmData = project_loc + self.name.replace(' ', '_')+'/'
        if not os.path.isdir(self.path2SwarmData): os.makedirs(self.path2SwarmData)

        self.project_name = self.name + '-' + str(self.counter)
        self.expected_project_file = self.path2SwarmData + "temp/%s.jproj"%(self.project_name)

        self.path2FEACsv = os.path.abspath(self.path2SwarmData + 'csv/') + '/'

        if not os.path.isdir(self.path2FEACsv): os.makedirs(self.path2FEACsv)

        if 'JMAG' in self.select_FEA_tool:
            study_name = "Transient" # Change here and there 

            # Leave the solving task to JMAG
            def build_jmag_project(study_name):
                import JMAG
                toolJd = JMAG.JMAG()
                toolJd.open(Steel_name=self.EX['SteelMaterial'], expected_project_file_path=self.expected_project_file, pc_name=pc_name, dir_parent=dir_parent, bool_jmagDesignerShow=bool_jmagDesignerShow)
                return toolJd

            def draw_spmsm(toolJd):
                import numpy as np
                # gray
                color_rgb_A = np.array([236,236,236])/255
                color_rgb_B = np.array([226,226,226])/255

                # Rotor Core
                list_regions_1 = self.machineGeometry['rotorCore'].draw(toolJd)
                if hasattr(toolJd, 'visualization_points') and 'rotorCore' in toolJd.visualization_points:
                    self.machineGeometry['rotorCore'].visualization_points = toolJd.visualization_points['rotorCore']
                toolJd.bMirror = False
                toolJd.iRotateCopy = self.machineGeometry['rotorCore'].p*2
                region1 = toolJd.prepareSection(list_regions_1, color=color_rgb_A)

                # print(list_regions_1)
                # raise KeyboardInterrupt

                # Shaft
                list_regions = self.machineGeometry['shaft'].draw(toolJd)
                if hasattr(toolJd, 'visualization_points') and 'shaft' in toolJd.visualization_points:
                    self.machineGeometry['shaft'].visualization_points = toolJd.visualization_points['shaft']
                toolJd.bMirror = False
                toolJd.iRotateCopy = 1
                region0 = toolJd.prepareSection(list_regions)

                # Rotor Magnet
                list_regions = self.machineGeometry['rotorMagnet'].draw(toolJd)
                if hasattr(toolJd, 'visualization_points') and 'rotorMagnet' in toolJd.visualization_points:
                    self.machineGeometry['rotorMagnet'].visualization_points = toolJd.visualization_points['rotorMagnet']
                toolJd.bMirror = False
                toolJd.iRotateCopy = self.machineGeometry['rotorCore'].p*2
                region2 = toolJd.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

                # Sleeve
                list_regions = self.machineGeometry['sleeve'].draw(toolJd)
                if hasattr(toolJd, 'visualization_points') and 'Sleeve' in toolJd.visualization_points:
                    self.machineGeometry['sleeve'].visualization_points = toolJd.visualization_points['Sleeve']
                toolJd.bMirror = False
                toolJd.iRotateCopy = self.machineGeometry['rotorCore'].p*2
                regionS = toolJd.prepareSection(list_regions)

                # Stator Core
                list_regions = self.machineGeometry['statorCore'].draw(toolJd)
                if hasattr(toolJd, 'visualization_points') and 'statorCore' in toolJd.visualization_points:
                    self.machineGeometry['statorCore'].visualization_points = toolJd.visualization_points['statorCore']
                toolJd.bMirror = True
                toolJd.iRotateCopy = self.machineGeometry['statorCore'].Q
                region3 = toolJd.prepareSection(list_regions, color=color_rgb_A)

                # Stator Winding
                list_regions = self.machineGeometry['coils'].draw(toolJd)
                if hasattr(toolJd, 'visualization_points') and 'Coils' in toolJd.visualization_points:
                    self.machineGeometry['coils'].visualization_points = toolJd.visualization_points['Coils']
                toolJd.bMirror = False
                toolJd.iRotateCopy = self.machineGeometry['statorCore'].Q
                region4 = toolJd.prepareSection(list_regions)

                # self.calculate_excitation_current(acm_variant)

                # Import Model into Designer
                toolJd.save(self.name, self.to_json())

            self.toolJd = toolJd = build_jmag_project(study_name)
            if 'PMSM' in self.name:
                draw_spmsm(self.toolJd)

            # JMAG
            app = toolJd.app
            model = app.GetModel(self.name)

            if 'PMSM' in self.name:
                toolJd.pre_process_PMSM(app, model, self)

            study = toolJd.add_magnetic_transient_study(app, model, self.path2FEACsv, study_name, self)
            toolJd.mesh_study(self, app, model, study, output_dir=self.path2SwarmData)
            # raise KeyboardInterrupt
            from time import time as clock_time
            toolJd.run_study(self, app, study, self.fea_config_dict, clock_time())

            # export Voltage if field data exists.
            if self.fea_config_dict['delete_results_after_calculation'] == False:
                # Export Circuit Voltage
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(self.path2FEACsv + study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")



            def compile_results(study_name, toolJd):
                self.results_to_be_unpacked = self.toolJd.build_str_results(self, self.project_name, study_name, self.path2FEACsv, self.fea_config_dict, femm_solver=None)
                cost_function, f1, f2, f3, FRW, \
                normalized_torque_ripple, \
                normalized_force_error_magnitude, \
                force_error_angle, \
                project_name, individual_name, \
                number_current_generation, individual_index,\
                power_factor, \
                rated_ratio, \
                rated_stack_length_mm, \
                rated_total_loss, \
                rated_stator_copper_loss_along_stack, \
                rated_magnet_Joule_loss, \
                rated_rotor_copper_loss_along_stack, \
                stator_copper_loss_in_end_turn, \
                rotor_copper_loss_in_end_turn, \
                rated_iron_loss, \
                rated_windage_loss, \
                str_results, \
                mm2_slot_area, \
                coil_flux_linkage_peak2peak_value, \
                TRV, Cost, Cost_Fe, Cost_Cu, Cost_PM, \
                ss_avg_force_magnitude, rotor_weight, torque_average = self.results_to_be_unpacked

                # acm_variant.spec_geometry_dict['x_denorm'] = list(x_denorm)

                self.spec_performance_dict = spec_performance_dict = dict()
                spec_performance_dict['x_denorm_dict'] = self.get_free_variables_as_dict() # ['x_denorm_dict']
                spec_performance_dict['project_name'] = project_name
                spec_performance_dict['individual_name'] = individual_name
                spec_performance_dict['number_current_generation'] = number_current_generation
                spec_performance_dict['individual_index'] = individual_index
                # spec_performance_dict['cost_function'] = cost_function
                spec_performance_dict['f1'] = f1
                spec_performance_dict['f2'] = f2
                spec_performance_dict['f3'] = float(f3)
                spec_performance_dict['TRV'] = TRV
                spec_performance_dict['FRW'] = FRW
                spec_performance_dict['torque_average'] = torque_average
                spec_performance_dict['ss_avg_force_magnitude'] = ss_avg_force_magnitude
                spec_performance_dict['rotor_weight'] = rotor_weight
                spec_performance_dict['normalized_torque_ripple'] = float(normalized_torque_ripple)
                spec_performance_dict['normalized_force_error_magnitude'] = float(normalized_force_error_magnitude)
                spec_performance_dict['force_error_angle'] = float(force_error_angle)
                spec_performance_dict['coil_flux_linkage_peak2peak_value'] = float(coil_flux_linkage_peak2peak_value)
                spec_performance_dict['mm2_slot_area'] = mm2_slot_area
                spec_performance_dict['Cost'] = Cost
                spec_performance_dict['Cost_Fe'] = Cost_Fe
                spec_performance_dict['Cost_Cu'] = Cost_Cu
                spec_performance_dict['Cost_PM'] = Cost_PM
                spec_performance_dict['power_factor'] = power_factor
                spec_performance_dict['rated_ratio'] = rated_ratio
                spec_performance_dict['rated_stack_length_mm'] = rated_stack_length_mm
                spec_performance_dict['rated_total_loss'] = rated_total_loss
                spec_performance_dict['rated_stator_copper_loss_along_stack'] = rated_stator_copper_loss_along_stack
                spec_performance_dict['rated_rotor_copper_loss_along_stack'] = rated_rotor_copper_loss_along_stack
                spec_performance_dict['rated_magnet_Joule_loss'] = rated_magnet_Joule_loss
                spec_performance_dict['stator_copper_loss_in_end_turn'] = stator_copper_loss_in_end_turn
                spec_performance_dict['rotor_copper_loss_in_end_turn'] = rotor_copper_loss_in_end_turn
                spec_performance_dict['rated_iron_loss'] = rated_iron_loss
                spec_performance_dict['rated_windage_loss'] = rated_windage_loss
                # spec_performance_dict['str_results'] = str_results
                spec_performance_dict['select_FEA_tool'] = self.select_FEA_tool
                spec_performance_dict['moo.fitness_OA'] = self.fea_config_dict['moo.fitness_OA']
                spec_performance_dict['moo.fitness_OB'] = self.fea_config_dict['moo.fitness_OB']
                spec_performance_dict['moo.fitness_OC'] = self.fea_config_dict['moo.fitness_OC']

                # GP = acm_variant.template.SI['GP']
                # EX = acm_variant.template.SI['EX']

                # Save to disk
                # self.save_to_disk(acm_variant, spec_performance_dict, GP, EX)

                number_current_generation = spec_performance_dict['number_current_generation'] #= int(acm_variant.counter//popsize), 
                individual_index = spec_performance_dict['individual_index'] #= acm_variant.counter
                # self.visualize_dict[f'FEA_Evaluated_Performance-{number_current_generation}-{individual_index}'] = spec_performance_dict
                json_file_path = self.path2SwarmData + self.name + f'-gen{number_current_generation}-ind{individual_index}.json'

                # Read the possibly-existing current json data
                try:
                    if os.path.getsize(json_file_path) > 0:
                        with open(json_file_path, 'r') as rf:
                            loaded_json = json.load(rf)
                    else:
                        loaded_json = {}
                except Exception:
                    loaded_json = {}

                # Compose new key
                # key = f'gen{number_current_generation}-ind{individual_index}'
                # loaded_json[key] = self.visualize_dict

                json_string = jsonpickle.encode(loaded_json, indent=4)
                with open(json_file_path, 'w+') as f:
                    f.write(json_string)

                number_current_generation = spec_performance_dict['number_current_generation'] #= int(acm_variant.counter//popsize), 
                individual_index = spec_performance_dict['individual_index'] #= acm_variant.counter

                # this is for optimization
                self.results_for_optimization = (cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle)

            compile_results(study_name, self.toolJd)

        elif 'FEMM' in self.select_FEA_tool:
            pass 
            # self.toolFEMM = self.build_femm_project(acm_variant)
            # acm_variant.results_to_be_unpacked = results_to_be_unpacked = toolFEMM.build_str_results(self.axeses, acm_variant, self.project_name, study_name, self.dir_csv_output_folder, self.fea_config_dict, femm_solver=None)
            # return acm_variant
        else:
            raise Exception('[acm_designer.py] Wrong string of select_FEA_tool:', self.select_FEA_tool)

    def start_optimization(self):

        import logging, datetime, os
        import builtins, utility_moo
        import pygmo as pg

        builtins.ad = self # share global variable between modules # https://stackoverflow.com/questions/142545/how-to-make-a-cross-module-variable
        ad = self
        import Problem_BearinglessSynchronousDesign # must import this after __builtins__.ad = ad


        def myLogger(dir_log, prefix='default_prefix_'): # This works even when the module is reloaded (which is not the case of the other answers) https://stackoverflow.com/questions/7173033/duplicate-log-output-when-using-python-logging-module

            logger = logging.getLogger()
            if not len(logger.handlers):
                logger.setLevel(logging.DEBUG)
                now = datetime.datetime.now()

                if not os.path.isdir(dir_log):
                    os.makedirs(dir_log)

                # create a file handler
                handler=logging.FileHandler(dir_log + prefix + '-' + now.strftime("%Y-%m-%d") +'.log')
                handler.setLevel(logging.DEBUG)

                # create a logging format
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)

                # add the handlers to the logger
                logger.addHandler(handler)
            return logger

        self.logger = myLogger(self.path2SwarmData, prefix='acmdm')
        logger = logging.getLogger(__name__)

        ################################################################
        # MOO Step 1:
        #   Create UserDefinedProblem and create population
        #   The magic method __init__ cannot be fined for UDP class
        ################################################################
        # [4.3.1] Basic setup
        _, prob = Problem_BearinglessSynchronousDesign.get_prob()
        popsize = self.fea_config_dict["moo.popsize"]
        logger.info(f'Pop size is {popsize}')
        # print('[acmop.py]', '-'*40 + '\n[acmop.py] Pop size is', popsize)

        # [4.3.2] Generate the pop
        if False:
            pop = pg.population(prob, size=popsize) 
        # Add Restarting Feature when generating pop
        else:

            # 检查swarm_data.txt，如果有至少一个数据，返回就不是None。
            # logger.info(f'Check for swarm data from: {self.select_spec}.json ...')
            self.ad.acm_template.build_x_denorm()
            # quit()
            swarm_data_file = ad.   read_swarm_data_json(self.select_spec, self.ad.acm_template.x_denorm_dict)
            
            number_of_chromosome = ad.analyzer.number_of_chromosome
            # print(number_of_chromosome)
            # quit()
            for index, (k, v) in enumerate(self.ad.acm_template.x_denorm_dict.items()):
                logger.info(f'x_denorm_dict variable no. {index} is {k} = {v} with bounds: {ad.acm_template.bounds_denorm[index]}')
            # quit()
            # case 1: swarm_data.txt exists # Restarting feature related codes
            if number_of_chromosome != 0:

                number_of_finished_iterations                       = number_of_chromosome // popsize
                number_of_finished_chromosome_in_current_generation = number_of_chromosome % popsize

                # 如果刚好整除，把余数0改为popsize
                if number_of_finished_chromosome_in_current_generation == 0:
                    number_of_finished_chromosome_in_current_generation = popsize
                    logger.info(f'\tThere are {number_of_chromosome} chromosomes found in {ad.swarm_data_file}.')
                    logger.info('\tWhat is the odds! The script just stopped when the evaluation of the whole pop is finished.')
                    logger.info('\tSet number_of_finished_chromosome_in_current_generation to popsize %d'%(number_of_finished_chromosome_in_current_generation))

                logger.info('This is a restart of '+ self.path2SwarmData)
                logger.info('\tNumber of finished iterations is %d'%(number_of_finished_iterations))
                # print('This means the initialization of the population class is interrupted. So the pop in swarm_data.txt is used as the survivor.')

                # 这些计数器的值永远都是评估过的chromosome的个数。
                ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome
                logger.info('ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome = %d', number_of_chromosome)

                # 禁止在初始化pop时运行有限元
                ad.flag_do_not_evaluate_when_init_pop = True

                # 初始化population，如果ad.flag_do_not_evaluate_when_init_pop是False，那么就说明是 new run，否则，整代个体的fitness都是[0,0,0]。
                pop = pg.population(prob, size=popsize)
                # quit()
                # 如果整代个体的fitness都是[0,0,0]，那就需要调用set_xf，把txt文件中的数据写入pop。如果发现数据的个数不够，那就调用set_x()来产生数据，形成初代个体。
                if ad.flag_do_not_evaluate_when_init_pop == True:
                    pop_array = pop.get_x()
                    # print(pop_array)
                    # quit()
                    if number_of_chromosome <= popsize: # 个体数不够一代的情况
                        for i in range(popsize):
                            if i < number_of_chromosome: #number_of_finished_chromosome_in_current_generation:
                                pop.set_xf(i, ad.   swarm_data[i][:-3], ad.   swarm_data[i][-3:])
                                # print(pop.set_xf(i, ad.   swarm_data[i][:-3], ad.   swarm_data[i][-3:]))
                                # quit()
                            else:
                                logger.info('Set "ad.flag_do_not_evaluate_when_init_pop" to False...')
                                ad.flag_do_not_evaluate_when_init_pop = False
                                logger.info('Calling pop.set_x()---this is a restart for individual#%d during pop initialization.', i)
                                logger.info('i=%d: call get_fevals: %s', i, prob.get_fevals()) # https://esa.github.io/pygmo2/problem.html?highlight=get_fevals#pygmo.problem.get_fevals
                                pop.set_x(i, pop_array[i]) # evaluate this guy
                    else:
                        # 新办法，直接从swarm_data.txt（相当于archive）中判断出当前最棒的群体
                        swarm_data_on_pareto_front = utility_moo.learn_about_the_archive(prob, ad.   swarm_data, popsize, self.fea_config_dict)
                        # print(swarm_data_on_pareto_front)
                        # quit()
                        for i in range(popsize):
                            pop.set_xf(i, swarm_data_on_pareto_front[i][:-3], swarm_data_on_pareto_front[i][-3:])
                            # quit()
                    # 必须放到这个if的最后，因为在 learn_about_the_archive 中是有初始化一个 pop_archive 的，会调用fitness方法。
                    ad.flag_do_not_evaluate_when_init_pop = False

            # case 2: swarm_data.txt does not exist
            else:
                number_of_finished_chromosome_in_current_generation = None
                number_of_finished_iterations = 0 # 实际上跑起来它不是零，而是一，因为我们认为初始化的一代也是一代。或者，我们定义number_of_finished_iterations = number_of_chromosome // popsize

                # case 2-A: swarm_data.txt does not exist and this is a whole new run.
                logger.info('Nothing exists in the archival json file. This is a whole new run.')
                ad.flag_do_not_evaluate_when_init_pop = False
                pop = pg.population(prob, size=popsize)

            # this flag must be false before moving on
            ad.flag_do_not_evaluate_when_init_pop = False

        logger.info(f'Pop is initialized:\n {pop}')
        # hv = pg.hypervolume(pop)
        # quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
        # logger.info('[acmop.py] quality_measure: %g'%(quality_measure))
        # raise KeyboardInterrupt

        # 初始化以后，pop.problem.get_fevals()就是popsize，但是如果大于popsize，说明“pop.set_x(i, pop_array[i]) # evaluate this guy”被调用了，说明还没输出过 survivors 数据，那么就写一下。
        if pop.problem.get_fevals() > popsize:
            logger.info('Write survivors.')
            ad.   write_swarm_survivor(pop, ad.counter_fitness_return)


        ################################################################
        # MOO Step 2:
        #   Select algorithm (another option is pg.nsga2())
        ################################################################
        # [4.3.3] Selecting algorithm
        # Don't forget to change neighbours to be below popsize (default is 20) decomposition="bi"
        algo = pg.algorithm(pg.moead(gen=1, weight_generation="grid", decomposition="tchebycheff", 
                                     neighbours=int(popsize/4), 
                                     CR=1, F=0.5, eta_m=20, 
                                     realb=0.9, 
                                     limit=2, preserve_diversity=True)) # https://esa.github.io/pagmo2/docs/python/algorithms/py_algorithms.html#pygmo.moead
        logger.info(f'{algo}')
        logger.info(f'\t MOEA/D neighbourhood size is set to 1/4 of the popsize as {int(popsize/4)}')
        # quit()

        ################################################################
        # MOO Step 3:
        #   Begin optimization
        ################################################################
        # [4.3.4] Begin optimization
        # number_of_chromosome = ad.   read_swarm_data(self.select_spec)
        # swarm_data_file = ad.   read_swarm_data_json(self.select_spec, self.ad.acm_template.x_denorm_dict)
        number_of_chromosome = ad.analyzer.number_of_chromosome
        number_of_finished_iterations = number_of_chromosome // popsize
        number_of_iterations = 500

        for _ in range(number_of_finished_iterations, number_of_iterations):
            msg = '[acmop.py] This is iteration #%d. '%(_)
            # print(msg)
            logger.info(msg)
            pop = algo.evolve(pop)

            msg += 'Write survivors to file. '
            ad.   write_swarm_survivor(pop, ad.counter_fitness_return)

            hv = pg.hypervolume(pop)
            quality_measure = hv.compute(ref_point=get_bad_fintess_values(machine_type='PMSM', ref=True)) # ref_point must be dominated by the pop's pareto front
            msg += 'Quality measure by hyper-volume: %g'% (quality_measure)
            # print('[acmop.py]', msg)
            logger.info(msg)

            utility_moo.my_print(ad, pop, _)
            # my_plot(fits, vectors, ndf)
        pass





    ''' 实用
    '''
    def get_bad_fintess_values(self, machine_type='IM', ref=False):
        # define bad values for different MOO objectives

        if ref == False:
            if 'IM' in machine_type:
                return 0, 0, 99
            elif 'PM' in machine_type:
                return 9999, 0, 999
        else:
            if 'IM' in machine_type:
                return 1,     10, 100
            elif 'PM' in machine_type:
                return 10000, 10, 1000

    def get_rotor_volume(self, stack_length=None):
        if stack_length is None:
            return math.pi*(self.mm_r_ro.value*1e-3)**2 * (self.EX['mm_stack_length_specified']*1e-3)
        else:
            return math.pi*(self.mm_r_ro.value*1e-3)**2 * (stack_length*1e-3)
    def get_rotor_weight(self, gravity=9.8, stack_length=None):
        
        def get_material_data():
            material_density_rho = 7860 # kg/m^3
            Poisson_ratio_nu = 0.29 # (i.e. the ratio of lateral contraction to longitudinal extension in the direction of the stretching force)
            Youngs_modulus_of_elasticity = 190 * 1e9 # Young's modulus for steel is in [190, 210] GPa
            return material_density_rho, Poisson_ratio_nu, Youngs_modulus_of_elasticity

        material_density_rho = get_material_data()[0]
        if stack_length is None:
            return gravity * self.get_rotor_volume() * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg
        else:
            return gravity * self.get_rotor_volume(stack_length=stack_length) * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg


    # =========== 变量管理 ===========

    def get_free_variables(self) -> List[Parameter]:
        return [param for param in self.get_parameter_fields().values() if param.type == 'free']

    def get_free_variables_as_dict(self) -> OrderedDict:
        """
        获取所有 free 类型参数的值字典（有序）
        
        Returns:
            OrderedDict: 参数字典，键为参数名，值为参数值
            顺序与 get_free_variable_bounds_dict() 保持一致
        """
        free_vars = self.get_free_variables()
        return OrderedDict((param.name, param.value) for param in free_vars)

    def get_free_variable_bounds_dict(self) -> OrderedDict:
        """
        获取所有 free 类型参数的边界值字典（有序）
        
        Returns:
            OrderedDict: 参数字典，键为参数名，值为边界值（bounds）
            如果参数没有边界值，则值为 None
            顺序与 get_free_variables_as_dict() 保持一致
        """
        free_vars = self.get_free_variables()
        return OrderedDict((param.name, param.bounds) for param in free_vars)

    def set_free_variables_from_dict(self, free_variables_dict: Dict[str, Any]) -> None:
        for name, value in free_variables_dict.items():
            param = self.get_parameter(name)
            if param is None:
                raise Exception(f"Parameter '{name}' not found")
            param.value = value

    def update_derived_parameters(self) -> None:
        """
        更新所有 derived 类型的参数值
        注意：这需要 calc 函数存在。如果从 JSON 反序列化，calc 函数可能为 None。
        在这种情况下，需要在反序列化后手动设置 calc 函数。
        """
        for name, param in self.get_parameter_fields().items():
            if param.type == 'derived' and param.calc is not None:
                try:
                    # calc 函数将在未来提供正确的定义
                    # 这里暂时不执行计算，等待用户提供 calc 函数的正确定义
                    pass
                except Exception as e:
                    # 如果计算失败，保持原值或设置为 None
                    print(f"Warning: Failed to calculate derived parameter '{name}': {e}")
    
    # ========== 自省方法：参数管理 ==========
    
    def get_parameter_fields(self) -> Dict[str, 'Modern_Machine_Designer.Parameter']:
        """
        获取所有 Parameter 类型的字段
        
        Returns:
            Dict[str, Parameter]: 字段名到 Parameter 对象的映射
        """
        param_fields = {}
        # 使用 vars(self) 或 self.__dict__ 来获取所有实例属性
        # 因为参数是在 __post_init__ 中动态添加的，不是 dataclass 字段
        for attr_name, attr_value in vars(self).items():
            if isinstance(attr_value, Parameter):
                param_fields[attr_name] = attr_value
        return param_fields
    
    def get_parameters_by_type(self, param_type: str) -> Dict[str, 'Modern_Machine_Designer.Parameter']:
        """
        根据参数类型筛选参数
        
        Args:
            param_type: 参数类型 ('fixed', 'derived', 'free' 等)
            
        Returns:
            Dict[str, Parameter]: 符合条件的参数字典
        """
        return {name: param for name, param in self.get_parameter_fields().items() 
                if param.type == param_type}
    
    def get_parameter(self, name: str) -> Optional['Modern_Machine_Designer.Parameter']:
        """
        根据字段名获取参数
        
        Args:
            name: 字段名
            
        Returns:
            Parameter 对象，如果不存在则返回 None
        """
        params = self.get_parameter_fields()
        return params.get(name)
    
    def set_parameter_value(self, name: str, value: Any) -> bool:
        """
        设置参数值
        
        Args:
            name: 字段名
            value: 要设置的值
            
        Returns:
            bool: 是否设置成功
        """
        param = self.get_parameter(name)
        if param is None:
            return False
        param.value = value
        return True
    
    def get_parameter_dict(self) -> OrderedDict:
        """
        获取所有参数的 OrderedDict（兼容旧代码风格）
        
        Returns:
            OrderedDict: 参数的有序字典
        """
        return OrderedDict(self.get_parameter_fields())
    
    def list_parameters(self, param_type: Optional[str] = None) -> List[str]:
        """
        列出所有参数名
        
        Args:
            param_type: 可选，如果指定则只返回该类型的参数名
            
        Returns:
            List[str]: 参数名列表
        """
        if param_type:
            return list(self.get_parameters_by_type(param_type).keys())
        return list(self.get_parameter_fields().keys())
    
    def get_parameters_summary(self) -> Dict[str, Any]:
        """
        获取参数摘要信息
        
        Returns:
            Dict: 包含参数统计信息的字典
        """
        params = self.get_parameter_fields()
        summary = {
            'total_count': len(params),
            'by_type': {},
            'by_unit': {},
            'with_values': 0,
            'without_values': 0
        }
        
        for param in params.values():
            # 按类型统计
            summary['by_type'][param.type] = summary['by_type'].get(param.type, 0) + 1
            # 按单位统计
            summary['by_unit'][param.unit] = summary['by_unit'].get(param.unit, 0) + 1
            # 统计有值/无值
            if param.value is not None:
                summary['with_values'] += 1
            else:
                summary['without_values'] += 1
        
        return summary
    
    def validate_parameters(self) -> Dict[str, List[str]]:
        """
        验证参数的有效性
        
        Returns:
            Dict: 包含 'errors' 和 'warnings' 的字典
        """
        errors = []
        warnings = []
        
        for name, param in self.get_parameter_fields().items():
            # 检查 derived 类型是否有 calc 方法
            if param.type == 'derived' and param.calc is None:
                warnings.append(f"Derived parameter '{name}' has no calculation method")
            
            # 检查 bounds 是否有效
            if param.bounds is not None:
                if isinstance(param.bounds, (list, tuple)) and len(param.bounds) == 2:
                    if param.bounds[0] is not None and param.bounds[1] is not None:
                        if param.bounds[0] > param.bounds[1]:
                            errors.append(f"Parameter '{name}' has invalid bounds: {param.bounds}")
            
            # 检查值是否在范围内
            if param.value is not None and param.bounds is not None:
                if isinstance(param.bounds, (list, tuple)) and len(param.bounds) == 2:
                    if param.bounds[0] is not None and param.value < param.bounds[0]:
                        errors.append(f"Parameter '{name}' value {param.value} is below lower bound {param.bounds[0]}")
                    if param.bounds[1] is not None and param.value > param.bounds[1]:
                        errors.append(f"Parameter '{name}' value {param.value} is above upper bound {param.bounds[1]}")
        
        return {'errors': errors, 'warnings': warnings}
    
    def __repr__(self):
        """提供类的字符串表示，包含参数摘要"""
        summary = self.get_parameters_summary()
        return (f"Modern_Machine_Designer(machine_class='{self.machine_class}', "
                f"parameters={summary['total_count']}, "
                f"with_values={summary['with_values']})")
    
    # ========== JSON 序列化和反序列化方法 ==========
    
    def to_dict(self) -> Dict[str, Any]:
        """
        将 Modern_Machine_Designer 对象的所有成员变量（包括自定义和全部属性）转换为字典。
        对自定义对象（如 Parameter/Winding/Geometry）自动调用其 to_dict 方法。
        """
        def serialize_value(val):
            # 基础类型直接返回
            if isinstance(val, (int, float, str, bool, type(None))):
                return val
            # Parameter类
            if isinstance(val, Parameter):
                return val.to_dict()
            # 假如有Winding类
            if hasattr(val, "to_dict") and callable(val.to_dict):
                return val.to_dict()
            # dict递归
            if isinstance(val, dict):
                return {k: serialize_value(v) for k, v in val.items()}
            # list/tuple递归
            if isinstance(val, (list, tuple)):
                return [serialize_value(x) for x in val]
            # Geometry特判（常见于 machineGeometry）
            if 'Geometry' in type(val).__name__ or hasattr(val, 'GP'):
                # 只序列化字段和参数名称
                geo_dict = {
                    'color': getattr(val, 'color', None),
                    'GP': [
                        gp.name if isinstance(gp, Parameter) else str(getattr(gp, 'name', gp))
                        for gp in getattr(val, 'GP', [])
                    ]
                }
                geo_dict['_needs_rebuild'] = True
                return geo_dict
            # 其它类型尝试转为字符串
            return str(val)

        result = {}
        for attr_name in vars(self):
            attr_val = getattr(self, attr_name)
            result[attr_name] = serialize_value(attr_val)
        return result

    def to_json(self, indent: Optional[int] = 2, ensure_ascii: bool = False) -> str:
        """
        将 Modern_Machine_Designer 对象转换为 JSON 字符串
        
        Args:
            indent: JSON 缩进空格数，None 表示不缩进
            ensure_ascii: 是否确保 ASCII 编码
            
        Returns:
            str: JSON 字符串
        """
        return json.dumps(self.to_dict(), indent=indent, ensure_ascii=ensure_ascii)
    
    def save_to_file(self, filepath: str, indent: Optional[int] = 2) -> None:
        """
        将对象保存到 JSON 文件
        
        Args:
            filepath: 文件路径
            indent: JSON 缩进空格数
        """
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(self.to_dict(), f, indent=indent, ensure_ascii=False)
    
    def to_dict_full(self) -> Dict[str, Any]:
        """
        将 Modern_Machine_Designer 对象转换为完整字典（包含所有信息）
        使用纯 JSON 序列化，不依赖 pickle，确保所有数据都能正确保存
        
        Returns:
            Dict: 包含完整对象信息的字典
        """
        import inspect
        
        def serialize_parameter_full(param: Parameter) -> Dict[str, Any]:
            """完整序列化 Parameter 对象，包括 lambda 函数的信息"""
            param_dict = param.to_dict()
            
            # 保存 args（用于重建 lambda 函数）
            if param.args is not None:
                # 序列化 args，将 Parameter 对象转换为名称引用
                serialized_args = []
                for arg in param.args:
                    if isinstance(arg, Parameter):
                        serialized_args.append({'_type': 'Parameter', 'name': arg.name})
                    elif isinstance(arg, (int, float, str, bool, type(None))):
                        serialized_args.append(arg)
                    elif isinstance(arg, (list, tuple)):
                        serialized_args.append([serialize_parameter_full(a) if isinstance(a, Parameter) else a for a in arg])
                    else:
                        serialized_args.append(str(arg))
                param_dict['args'] = serialized_args
            
            # 尝试保存 lambda 函数的源代码（如果可能）
            calc_info = {}
            if param.calc is not None and callable(param.calc):
                try:
                    # 尝试获取源代码
                    source_lines = inspect.getsourcelines(param.calc)
                    if source_lines and len(source_lines) > 0:
                        source = ''.join(source_lines[0]).strip()
                        calc_info['source'] = source
                        calc_info['has_calc'] = True
                except Exception:
                    calc_info['has_calc'] = True
                    calc_info['source'] = None
            else:
                calc_info['has_calc'] = False
            
            if param.calc_bounds is not None and callable(param.calc_bounds):
                try:
                    source_lines = inspect.getsourcelines(param.calc_bounds)
                    if source_lines and len(source_lines) > 0:
                        source = ''.join(source_lines[0]).strip()
                        calc_info['calc_bounds_source'] = source
                        calc_info['has_calc_bounds'] = True
                except Exception:
                    calc_info['has_calc_bounds'] = True
                    calc_info['calc_bounds_source'] = None
            else:
                calc_info['has_calc_bounds'] = False
            
            if calc_info:
                param_dict['_calc_info'] = calc_info
            
            return param_dict
        
        def serialize_geometry_full(geo: Geometry) -> Dict[str, Any]:
            """完整序列化 Geometry 对象"""
            geo_dict = geo.to_dict()
            
            # 确保 visualization_points 被保存
            if hasattr(geo, 'visualization_points') and geo.visualization_points:
                geo_dict['visualization_points'] = geo.visualization_points
            
            # 保存 GP 中的 Parameter 对象引用
            if hasattr(geo, 'GP') and geo.GP:
                if isinstance(geo.GP, dict):
                    gp_dict = {}
                    for k, v in geo.GP.items():
                        if isinstance(v, Parameter):
                            gp_dict[k] = {'_type': 'Parameter', 'name': v.name}
                        else:
                            gp_dict[k] = v
                    geo_dict['GP'] = gp_dict
                elif isinstance(geo.GP, list):
                    geo_dict['GP'] = [
                        {'_type': 'Parameter', 'name': gp.name} if isinstance(gp, Parameter) else gp
                        for gp in geo.GP
                    ]
            
            return geo_dict
        
        # 开始构建完整字典
        full_dict = {}
        
        # 1. 保存基本属性
        full_dict['name'] = self.name
        full_dict['bool_PermanentMagnet'] = self.bool_PermanentMagnet
        full_dict['bool_StatorSlotClosed'] = self.bool_StatorSlotClosed
        full_dict['bool_RotorNotched'] = self.bool_RotorNotched
        full_dict['select_FEA_tool'] = self.select_FEA_tool
        full_dict['select_fea_config_dict'] = self.select_fea_config_dict
        full_dict['bool_jmagDeleteResultsAfterCalculation'] = self.bool_jmagDeleteResultsAfterCalculation
        full_dict['counter'] = self.counter
        
        # 2. 保存所有 Parameter 对象（完整版本）
        parameters_dict = {}
        for field_name, param in self.get_parameter_fields().items():
            parameters_dict[field_name] = serialize_parameter_full(param)
        full_dict['parameters'] = parameters_dict
        
        # 3. 保存 Winding 对象
        if hasattr(self, 'wily') and self.wily is not None:
            full_dict['wily'] = self.wily.to_dict()
        
        # 4. 保存 EX 字典（激励参数）
        if hasattr(self, 'EX') and self.EX is not None:
            full_dict['EX'] = self.EX.copy()
        
        # 5. 保存 machineGeometry（完整版本）
        if hasattr(self, 'machineGeometry') and self.machineGeometry is not None:
            machine_geometry_dict = {}
            for key, geo in self.machineGeometry.items():
                if geo is None:
                    machine_geometry_dict[key] = None
                elif isinstance(geo, Geometry):
                    machine_geometry_dict[key] = serialize_geometry_full(geo)
                else:
                    machine_geometry_dict[key] = str(geo)
            full_dict['machineGeometry'] = machine_geometry_dict
        
        # 6. 保存其他可能存在的属性
        other_attrs = ['path2SwarmData', 'project_name', 'expected_project_file', 
                      'path2FEACsv', 'toolJd', 'results_to_be_unpacked', 
                      'spec_performance_dict', 'results_for_optimization',
                      'InitialRotationAngle', 'drawer']
        for attr in other_attrs:
            if hasattr(self, attr):
                attr_val = getattr(self, attr)
                # 只保存可序列化的属性
                if isinstance(attr_val, (int, float, str, bool, type(None), dict, list)):
                    full_dict[attr] = attr_val
                elif isinstance(attr_val, Parameter):
                    full_dict[attr] = serialize_parameter_full(attr_val)
                elif hasattr(attr_val, 'to_dict'):
                    try:
                        full_dict[attr] = attr_val.to_dict()
                    except Exception:
                        pass  # 跳过无法序列化的对象
        
        # 7. 保存元数据
        full_dict['_metadata'] = {
            'class_name': self.__class__.__name__,
            'module': self.__class__.__module__,
            'has_machineGeometry': hasattr(self, 'machineGeometry') and self.machineGeometry is not None,
            'has_wily': hasattr(self, 'wily') and self.wily is not None,
            'has_EX': hasattr(self, 'EX') and self.EX is not None,
            'serialization_version': '2.0',  # 版本号，用于未来兼容性
            # 添加路径和项目信息到 metadata
            'path2SwarmData': getattr(self, 'path2SwarmData', None),
            'project_name': getattr(self, 'project_name', None),
            'expected_project_file': getattr(self, 'expected_project_file', None),
            'path2FEACsv': getattr(self, 'path2FEACsv', None),
        }
        
        return full_dict
    
    def save_to_file_full(self, filepath: str = 'machine_designer_full.json', indent: Optional[int] = 2) -> None:
        """
        将对象完整信息保存到 JSON 文件（类似 pickle 保存）
        
        Args:
            filepath: 文件路径，默认为 'machine_designer_full.json'
            indent: JSON 缩进空格数
        """
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(self.to_dict_full(), f, indent=indent, ensure_ascii=False)
    
    @classmethod
    def from_dict_full(cls, data: Dict[str, Any]) -> 'Modern_Machine_Designer':
        """
        从完整字典创建 Modern_Machine_Designer 对象（从完整 JSON 数据恢复）
        
        Args:
            data: 包含完整对象数据的字典
            
        Returns:
            Modern_Machine_Designer: 重建的对象
        """
        # 检查是否是旧格式（包含 pickle 数据）
        if '_pickle_data' in data and data['_pickle_data'] is not None:
            try:
                # 尝试从旧格式恢复
                pickled_base64 = data['_pickle_data']
                pickled_data = base64.b64decode(pickled_base64.encode('utf-8'))
                instance = pickle.loads(pickled_data)
                instance.__post_init__()
                return instance
            except Exception as e:
                print(f"Warning: Failed to unpickle object, using standard deserialization: {e}")
                # 回退到新格式
                pass
        
        # 新格式：从纯 JSON 数据恢复
        # 创建对象实例
        instance = cls.__new__(cls)
        
        # 恢复基本属性
        instance.name = data.get('name', 'SPMSM')
        instance.bool_PermanentMagnet = data.get('bool_PermanentMagnet', True)
        instance.bool_StatorSlotClosed = data.get('bool_StatorSlotClosed', False)
        instance.bool_RotorNotched = data.get('bool_RotorNotched', True)
        instance.select_FEA_tool = data.get('select_FEA_tool', 'JMAG Designer')
        instance.select_fea_config_dict = data.get('select_fea_config_dict', '#0213 JMAG Bearingless Sub-hamonics')
        instance.bool_jmagDeleteResultsAfterCalculation = data.get('bool_jmagDeleteResultsAfterCalculation', False)
        instance.counter = data.get('counter', 0)
        instance.fea_config_dict = data.get('fea_config_dict', None)
        
        # 恢复 Parameter 对象
        parameters_data = data.get('parameters', {})
        for field_name, param_data in parameters_data.items():
            # 从字典恢复 Parameter
            param = Parameter.from_dict(param_data)
            # 注意：calc 和 calc_bounds 函数需要在 __post_init__ 中重建
            setattr(instance, field_name, param)
        
        # 恢复 Winding 对象
        if 'wily' in data and data['wily'] is not None:
            instance.wily = Winding.from_dict(data['wily'])
        
        # 恢复 EX 字典
        if 'EX' in data:
            instance.EX = data['EX'].copy()
        
        # 恢复其他属性
        for attr in ['path2SwarmData', 'project_name', 'expected_project_file', 
                    'path2FEACsv', 'results_to_be_unpacked', 
                    'spec_performance_dict', 'results_for_optimization',
                    'InitialRotationAngle']:
            if attr in data:
                setattr(instance, attr, data[attr])
        
        # 调用 __post_init__ 来重建 lambda 函数、machineGeometry 和其他依赖项
        # 注意：__post_init__ 会使用已设置的参数值
        # 但是，如果参数值已经存在，我们需要确保它们不会被覆盖
        # 所以我们需要在调用 __post_init__ 之前保存当前值，然后在之后恢复
        
        # 保存当前参数值
        saved_param_values = {}
        for field_name in parameters_data.keys():
            param = getattr(instance, field_name, None)
            if param is not None:
                saved_param_values[field_name] = param.value
        
        # 调用 __post_init__ 重建所有内容
        instance.__post_init__()
        
        # 恢复保存的参数值（如果 __post_init__ 覆盖了它们）
        for field_name, saved_value in saved_param_values.items():
            param = getattr(instance, field_name, None)
            if param is not None and param.value != saved_value:
                param.value = saved_value
        
        # 恢复 machineGeometry（如果 JSON 中有保存）
        if 'machineGeometry' in data and data['machineGeometry'] is not None:
            # 注意：machineGeometry 中的 visualization_points 会被保留
            # 但 draw_function 需要在 __post_init__ 中重建
            for key, geo_data in data['machineGeometry'].items():
                if geo_data is None:
                    instance.machineGeometry[key] = None
                elif isinstance(geo_data, dict):
                    # 尝试恢复 Geometry 对象
                    # 由于 draw_function 无法序列化，我们只恢复可序列化的部分
                    if hasattr(instance, 'machineGeometry') and key in instance.machineGeometry:
                        geo = instance.machineGeometry[key]
                        if geo is not None and isinstance(geo, Geometry):
                            # 恢复 visualization_points
                            if 'visualization_points' in geo_data:
                                geo.visualization_points = geo_data['visualization_points']
        
        return instance
    
    @classmethod
    def load_from_file_full(cls, filepath: str = 'machine_designer_full.json') -> 'Modern_Machine_Designer':
        """
        从完整 JSON 文件加载对象（从 pickle 数据恢复）
        
        Args:
            filepath: 文件路径，默认为 'machine_designer_full.json'
            
        Returns:
            Modern_Machine_Designer: 重建的对象
        """
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return cls.from_dict_full(data)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Modern_Machine_Designer':
        """
        从字典创建 Modern_Machine_Designer 对象（用于 JSON 反序列化）
        
        Args:
            data: 包含对象数据的字典
            
        Returns:
            Modern_Machine_Designer: 重建的对象
        """
        # 首先获取所有字段的默认值
        field_defaults = {}
        for field in fields(cls):
            if hasattr(cls, field.name):
                field_defaults[field.name] = getattr(cls, field.name)
        
        # 创建参数字典，从 JSON 数据中恢复，缺失的使用默认值
        param_dict = {}
        json_params = data.get('parameters', {})
        
        # 先设置所有默认字段
        for field_name, default_value in field_defaults.items():
            if isinstance(default_value, Parameter):
                # 如果 JSON 中有这个参数，使用 JSON 数据；否则使用默认值
                if field_name in json_params:
                    param_dict[field_name] = Parameter.from_dict(json_params[field_name])
                else:
                    # 创建默认参数的副本
                    param_dict[field_name] = Parameter(
                        name=default_value.name,
                        type=default_value.type,
                        value=default_value.value,
                        bounds=default_value.bounds,
                        unit=default_value.unit,
                        comment=default_value.comment,
                        calc=default_value.calc,
                        calc_bounds=default_value.calc_bounds
                    )
        
        # 创建对象实例
        instance = cls.__new__(cls)
        
        # 设置基本字段
        instance.machine_class = data.get('machine_class', field_defaults.get('machine_class', cls.machine_class))
        instance.bool_PermanentMagnet = data.get('bool_PermanentMagnet', field_defaults.get('bool_PermanentMagnet', True))
        instance.bool_StatorSlotClosed = data.get('bool_StatorSlotClosed', field_defaults.get('bool_StatorSlotClosed', False))
        instance.bool_RotorNotched = data.get('bool_RotorNotched', field_defaults.get('bool_RotorNotched', True))
        
        # 设置所有参数字段
        for field_name, param in param_dict.items():
            setattr(instance, field_name, param)
        
        # 重建 Winding 对象
        if 'winding' in data:
            instance.wily = Winding.from_dict(data['winding'])
        
        # 调用 __post_init__ 来创建 machineGeometry 和其他依赖项
        # 注意：__post_init__ 会使用已设置的参数值来创建 machineGeometry
        instance.__post_init__()
        
        return instance
    
    @classmethod
    def from_json(cls, json_str: str) -> 'Modern_Machine_Designer':
        """
        从 JSON 字符串创建 Modern_Machine_Designer 对象
        
        Args:
            json_str: JSON 字符串
            
        Returns:
            Modern_Machine_Designer: 重建的对象
        """
        data = json.loads(json_str)
        return cls.from_dict(data)
    
    @classmethod
    def load_from_file(cls, filepath: str) -> 'Modern_Machine_Designer':
        """
        从 JSON 文件加载对象
        
        Args:
            filepath: 文件路径
            
        Returns:
            Modern_Machine_Designer: 重建的对象
        """
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return cls.from_dict(data)

if __name__ == "__main__":
    # 创建对象并导出为 JSON
    mmd = Modern_Machine_Designer()

    mmd.show_geometry()
    # print(dir(mmd.machineGeometry['statorCore']))
    mmd.drawer.visualization_points['Coils']['PCoil']

    mmd.FEA_evaluate()
    mmd.save_to_file('machine_designer.json')
    mmd.save_to_file_full('machine_designer_full.json') # 保存完整信息到文件（类似 pickle）

    mmd.start_optimization()
    quit()


    # 从 JSON 文件恢复对象
    mmd2 = Modern_Machine_Designer.load_from_file('machine_designer.json')

    # 从完整文件恢复对象（包含所有信息，包括 lambda 函数）
    mmd3 = Modern_Machine_Designer.load_from_file_full('machine_designer_full.json')
    print("=== 已从完整文件恢复对象 ===")


