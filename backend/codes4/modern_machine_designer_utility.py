from typing import Dict, List, Optional, Any
from collections import OrderedDict
import json, math, base64, pickle, cairo, os, jsonpickle, logging, utility, JMAG

class Parameter(object):
    def __init__(self, name, type, value=None, bounds=None, calc=None, calc_bounds=None, unit='mm', parameter_dict=None) -> None:
        self.name = name
        self.type = type
        self.value = value
        self.bounds = bounds
        self.unit = unit
        self.calc = calc
        self.calc_bounds = calc_bounds
        self.parameter_dict = parameter_dict
        self.overwritten: bool = False  # 标记是否被外部强制覆盖
        # todo: add validation for the type, value, bounds, calc, unit, comment
        if self.calc is not None and self.parameter_dict is not None:
            try: 
                self.value = self.calc(self.parameter_dict)
                self.initialized = True
            except KeyError as e:
                print(f"The derivation of {self.name} failed due to KeyError: {e}. Need to run calc method again later.")
                self.initialized = False

        if self.calc_bounds is not None:
            try:
                # 如果 parameter_dict 是 None，lambda 函数可能使用闭包中的 self
                # 在这种情况下，传递 None 作为参数，lambda 函数会使用闭包中的 self
                if self.parameter_dict is None:
                    # 尝试调用 calc_bounds，lambda 函数会使用闭包中的 self
                    self.bounds = self.calc_bounds(None)
                else:
                    self.bounds = self.calc_bounds(self.parameter_dict)
            except Exception as e:
                # 如果 calc_bounds 失败，保持原有 bounds
                print(f"Warning: calc_bounds failed for parameter '{self.name}': {e}")

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
            # calc 函数无法序列化，保存为 None，前端可以设置 calc_dependencies
            'calc_dependencies': None,  # 可以扩展为保存依赖的参数名列表
            'overwritten': getattr(self, 'overwritten', False),
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
        param = cls(
            name=data['name'],
            type=data['type'],
            value=data.get('value'),
            bounds=data.get('bounds'),
            unit=data.get('unit', 'mm'),
            calc=None,  # calc 函数需要从其他地方重建
            calc_bounds=None,  # calc_bounds 函数需要从其他地方重建
            parameter_dict=data.get('parameter_dict')  # 保存 parameter_dict，用于重建 lambda 函数
        )
        param.overwritten = data.get('overwritten', False)
        return param

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
        # Default layer attributes
        self.layer_X_phases = None
        self.layer_X_signs = None
        self.layer_Y_phases = None
        self.layer_Y_signs = None

        derivation = self.get_winding_factor()
        # Get all winding factor information and phase/sign/grouping data from self.derivation
        # If self.derivation is a dict or an object with these attributes, extract them.
        # Fallback to defaults if not present.

        if derivation is not None:
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

            # Populate dict_coil_connection for JMAG.py compatibility
            self.dict_coil_connection = {
                'layer X phases': self.layer_X_phases,
                'layer X signs': self.layer_X_signs,
                'layer Y phases': self.layer_Y_phases,
                'layer Y signs': self.layer_Y_signs,
            }


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
            
            # Update dict_coil_connection for DPNV defaults
            self.dict_coil_connection = {
                'layer X phases': self.layer_X_phases,
                'layer X signs': self.layer_X_signs,
                'layer Y phases': self.layer_Y_phases,
                'layer Y signs': self.layer_Y_signs,
            }

    def infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self, layer_X_phases, coil_pitch):
        return layer_X_phases[-coil_pitch:] + layer_X_phases[:-coil_pitch]
    def infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self, layer_X_signs, coil_pitch):
        temp = layer_X_signs[-coil_pitch:] + layer_X_signs[:-coil_pitch]
        return [('-' if el == '+' else '+') for el in temp]

    def get_winding_factor(self):
        import winding_layout_derivation_ismb2021_asymetry_no_drawing
        derivation = winding_layout_derivation_ismb2021_asymetry_no_drawing.main_derivation(m=self.m, Qs=self.Qs, p=self.p, ps=self.ps, coil_pitch_y=self.coil_pitch_y)
        return derivation

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
        u.cvs.writePDFfile(self.path2SwarmData + '/part_winding_pyx_output_M')

        u = PyX_Utility.PyX_Utility()
        draw_winding_in_the_slot(u, wily.Qs, wily.list_layer_suspension_phases, wily.list_layer_suspension_signs, text=' Suspension Mode' )
        u.cvs.writePDFfile(self.path2SwarmData + '/part_winding_pyx_output_S')
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
        """将 Winding 对象转换为字典，导出所有成员变量"""
        result = {
            'phase_number_m': self.m,
            'stator_slot_number_Qs': self.Qs,
            'pole_pair_number_p': self.p,
            'suspension_pole_pair_number_ps': self.ps,
            'coil_pitch_y': self.coil_pitch_y,
            'number_of_parallel_branch': self.number_of_parallel_branch,
            'kw1': self.kw1,
        }
        
        # 导出所有其他成员变量
        additional_attrs = [
            'SPP', 'deg_winding_U_phase_phase_axis_angle', 'number_of_winding_layer',
            'bool_distributed_or_concentrated', 'bool_DPNVorSEPA', 
            'bool_3PhaseCurrentSource', 'bool_CustomizedCircuit',
            'CommutatingSequenceD', 'CommutatingSequenceB',
            'layer_X_phases', 'layer_X_signs', 'layer_Y_phases', 'layer_Y_signs',
            'grouping_AC'
        ]
        
        for attr in additional_attrs:
            if hasattr(self, attr):
                value = getattr(self, attr)
                # 只序列化可序列化的类型
                if isinstance(value, (int, float, str, bool, list, tuple, type(None))):
                    result[attr] = value
                elif isinstance(value, dict):
                    result[attr] = value
                else:
                    # 对于其他类型，尝试转换为字符串或跳过
                    try:
                        result[attr] = str(value)
                    except Exception:
                        pass  # 跳过无法序列化的属性
        
        # 导出 derivation（如果存在且可序列化）
        if hasattr(self, 'derivation') and self.derivation is not None:
            try:
                if isinstance(self.derivation, dict):
                    result['derivation'] = self.derivation
                elif hasattr(self.derivation, '__dict__'):
                    # 尝试将对象转换为字典
                    result['derivation'] = {k: v for k, v in self.derivation.__dict__.items() 
                                          if isinstance(v, (int, float, str, bool, list, tuple, dict, type(None)))}
                else:
                    result['derivation'] = str(self.derivation)
            except Exception:
                result['derivation'] = None
        
        return result
    
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
        self.GP = GP or {}
        self.draw_function = draw_function
        for name, gp in self.GP.items():
            if isinstance(gp, Parameter):
                setattr(self, name, gp.value)
    
    def update_from_GP(self):
        """
        从 GP 字典中的 Parameter 对象更新所有实例属性值
        当 Parameter 的值更新后，调用此方法来同步 Geometry 对象的属性值
        """
        if isinstance(self.GP, dict):
            for name, gp in self.GP.items():
                if isinstance(gp, Parameter):
                    # 更新实例属性值为 Parameter 的当前值
                    setattr(self, name, gp.value)
    
    def __repr__(self):
        """返回 Geometry 对象的基本信息"""
        gp_info = {}
        if isinstance(self.GP, dict):
            for k, v in self.GP.items():
                if isinstance(v, Parameter):
                    gp_info[k] = f"Parameter(name='{v.name}', value={v.value}, unit='{v.unit}')"
                else:
                    gp_info[k] = str(v)
        return f"Geometry(name='{self.name}', color='{self.color}', GP={gp_info})"
    
    def print_parameters(self):
        """打印所有参数值，包括 GP 中的参数和实例属性"""
        print(f"\n=== Geometry: {self.name} ===")
        print(f"Color: {self.color}")
        
        # 打印 GP 中的参数
        if self.GP:
            print("\nGP Parameters:")
            if isinstance(self.GP, dict):
                for param_name, param_value in self.GP.items():
                    if isinstance(param_value, Parameter):
                        print(f"  {param_name}: Parameter(name='{param_value.name}', type='{param_value.type}', value={param_value.value}, bounds={param_value.bounds}, unit='{param_value.unit}')")
                    else:
                        print(f"  {param_name}: {param_value}")
            elif isinstance(self.GP, list):
                for i, gp in enumerate(self.GP):
                    if isinstance(gp, Parameter):
                        print(f"  [{i}]: Parameter(name='{gp.name}', type='{gp.type}', value={gp.value}, bounds={gp.bounds}, unit='{gp.unit}')")
                    else:
                        print(f"  [{i}]: {gp}")
        
        # 打印所有实例属性（排除 GP 和 draw_function）
        print("\nInstance Attributes:")
        for attr_name, attr_value in self.__dict__.items():
            if attr_name not in ['GP', 'draw_function']:
                if isinstance(attr_value, Parameter):
                    print(f"  {attr_name}: Parameter(name='{attr_value.name}', type='{attr_value.type}', value={attr_value.value}, bounds={attr_value.bounds}, unit='{attr_value.unit}')")
                elif isinstance(attr_value, (int, float, str, bool, type(None))):
                    print(f"  {attr_name}: {attr_value}")
                elif isinstance(attr_value, (list, tuple)):
                    print(f"  {attr_name}: {attr_value}")
                elif isinstance(attr_value, dict):
                    print(f"  {attr_name}: {attr_value}")
                else:
                    print(f"  {attr_name}: {type(attr_value).__name__} object")
        
        # 打印 visualization_points（如果存在）
        if hasattr(self, 'visualization_points') and self.visualization_points:
            print("\nVisualization Points:")
            if isinstance(self.visualization_points, dict):
                for key, value in self.visualization_points.items():
                    print(f"  {key}: {value}")
            else:
                print(f"  {self.visualization_points}")
        
        print("=" * 50)
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
    def __init__(self, width_in_points=500, height_in_points=500, filename=None, verbose_drawing=False, scale=1.0, bFillRegion=True):
        self.filename = filename
        self.verbose_drawing = verbose_drawing
        self.scale = scale
        self.bFillRegion = bFillRegion
        self.iRotateCopy = 0
        self.bMirror = False
        self.ctx = None
        self.surface = None
        if self.filename:
            self.surface = cairo.SVGSurface(self.filename, width_in_points, height_in_points)
            self.ctx = cairo.Context(self.surface)
        # self.ctx.scale(width_in_points, height_in_points)
        # m = cairo.Matrix(yy=-1, y0=height_in_points) # Cartetian Coordinate
        m = cairo.Matrix(yy=-1, y0=0.5*height_in_points, x0=+0.5*width_in_points) # Offset to center
        self.ctx.transform(m)
        self.ctx.scale(self.scale, self.scale)
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
        try:
            import cairosvg
            cairosvg.svg2pdf(url=filename or f'machine_geometry.svg', write_to=f'machine_geometry.pdf')
            if bool_open_pdf:
                import os
                os.system('sumatraPDF2.exe ' + 'machine_geometry.pdf')
            print('[machine_design_guide.py] Find the file machine_geometry.pdf in the current folder.')
        except ImportError:
            print('[machine_design_guide.py] cairosvg not found, skipping PDF conversion.')
        except Exception as e:
            print(f'[machine_design_guide.py] PDF conversion failed: {e}')

        self.sketch_color = color

    def hex_to_rgb(self, hex_color):
        if not hex_color or not isinstance(hex_color, str):
            return (0.0, 0.0, 0.0)
        hex_color = hex_color.lstrip('#')
        lv = len(hex_color)
        try:
            if lv == 3:
                rgb = tuple(int(hex_color[i:i+1]*2, 16) for i in range(0, 3))
            elif lv == 6:
                rgb = tuple(int(hex_color[i:i+2], 16) for i in range(0, 6, 2))
            else:
                return (0.0, 0.0, 0.0)
            return tuple(c/255.0 for c in rgb)
        except ValueError:
            return (0.0, 0.0, 0.0)

    def drawLine(self, p1, p2):
        if self.verbose_drawing:
            print(f'[CairoDrawer.py] drawLine({p1=}, {p2=})')
        return [{'type': 'line', 'p1': p1, 'p2': p2}]

    def drawArc(self, centerxy, startxy, endxy):
        if self.verbose_drawing:
            print(f'[CairoDrawer.py] drawArc({centerxy=}, {startxy=}, {endxy=})')
        return [{'type': 'arc', 'center': centerxy, 'p1': startxy, 'p2': endxy}]

    def getSketch(self, name, color):
        # We handle sketches (paths) in prepareSection now
        pass

    def prepareSection(self, region_dict, color=None, **kwargs):
        list_regions = region_dict.get('list_regions', [])
        
        # Priority: kwargs > region_dict (from part) > Drawer Attributes
        bMirror = kwargs.get('bMirror', region_dict.get('bMirror', getattr(self, 'bMirror', False)))
        copyCount = kwargs.get('iRotateCopy', region_dict.get('iRotateCopy', getattr(self, 'iRotateCopy', 0)))
        
        rgb = self.hex_to_rgb(color or "#888888")
        
        # If copyCount is 0 or 1, we just do one rotation (base)
        loop_count = copyCount if copyCount >= 2 else 1
        
        if self.ctx is None:
            return region_dict

        for i in range(loop_count):
            rotation_offset = i * (360.0 / loop_count) if loop_count > 1 else 0.0
            self.ctx.save()
            self.ctx.rotate(math.radians(rotation_offset))
            
            # Helper to draw the segments
            def draw_region_content(segments_list):
                for segments in segments_list:
                    if not segments: continue
                    
                    # Track current position to avoid unnecessary move_to or unwanted lines
                    current_pos = None
                    
                    for seg in segments:
                        p1 = seg['p1']
                        p2 = seg['p2']
                        
                        # If not at the start of the segment, move_to p1
                        # Use a small epsilon for float comparison if needed, but simple compare is usually fine for exact coordinates
                        if current_pos is None or current_pos != p1:
                            self.ctx.move_to(p1[0], p1[1])
                        
                        if seg['type'] == 'line':
                            self.ctx.line_to(p2[0], p2[1])
                            current_pos = p2
                        elif seg['type'] == 'arc':
                            center = seg['center']
                            v1 = [p1[0] - center[0], p1[1] - center[1]]
                            v2 = [p2[0] - center[0], p2[1] - center[1]]
                            radius = math.sqrt(v1[0]**2 + v1[1]**2)
                            angle_start = math.atan2(v1[1], v1[0])
                            
                            dot = v1[0]*v2[0] + v1[1]*v2[1]
                            mag = radius * math.sqrt(v2[0]**2 + v2[1]**2)
                            cos_val = max(-1.0, min(1.0, dot / mag))
                            angle_diff = math.acos(cos_val)
                            cross_prod = v1[0]*v2[1] - v1[1]*v2[0]
                            if cross_prod < 0:
                                self.ctx.arc_negative(center[0], center[1], radius, angle_start, angle_start - angle_diff)
                            else:
                                self.ctx.arc(center[0], center[1], radius, angle_start, angle_start + angle_diff)
                            current_pos = p2

            # Draw original
            self.ctx.new_path()
            draw_region_content(list_regions)
            
            # Fill and Stroke
            if self.bFillRegion:
                self.ctx.set_source_rgba(*rgb, 0.6)
                self.ctx.fill_preserve()
            
            self.ctx.set_source_rgba(*rgb, 1.0)
            self.ctx.set_line_width(0.01)
            self.ctx.stroke()
            
            if bMirror:
                self.ctx.scale(1, -1)
                self.ctx.new_path()
                draw_region_content(list_regions)
                if self.bFillRegion:
                    self.ctx.set_source_rgba(*rgb, 0.6)
                    self.ctx.fill_preserve()
                self.ctx.set_source_rgba(*rgb, 1.0)
                self.ctx.stroke()

            self.ctx.restore()

    def finalize_part(self):
        # Decommissioned in favor of prepareSection
        pass

class Modern_Machine_Designer_Utility(object):
    def __init__(self, specs: 'MotorSpecs'):
        self.specs = specs
        self.geometry = specs.geometry
        self.winding = specs.winding
        self.target = specs.targets
        # Keep machine_class as an attribute for now if it's used directly, 
        # but let's try to transition to self.target.machine_class
        # self.machine_class = specs.targets.machine_class 
        self.flag_do_not_evaluate_when_init_pop = False

    def _get_parameter_logger(self):
        """
        获取用于参数覆盖的日志记录器，确保日志写入文件。
        """
        logger_name = f"{__name__}.parameter_override"
        logger = logging.getLogger(logger_name)
        if not getattr(logger, "_acm_param_handler", False):
            log_dir = getattr(self, "path2SwarmData", os.getcwd())
            try:
                os.makedirs(log_dir, exist_ok=True)
            except Exception:
                pass
            log_file = os.path.join(log_dir, f"{self.name}-parameter_override.log")
            handler = logging.FileHandler(log_file, encoding='utf-8')
            handler.setLevel(logging.INFO)
            handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
            logger._acm_param_handler = True
        return logger

    def apply_parameter_dict(self, prev_params: Dict[str, Any], key_map: Optional[Dict[str, str]] = None):
        """
        根据提供的参数字典更新 mmd 参数。
        - free 参数直接更新。
        - 非 free 参数被覆盖时会记录日志并标记 overwritten。
        - 最后刷新 derived 参数（跳过已被标记 overwritten 的）。
        """
        if not prev_params:
            return

        logger = self._get_parameter_logger()
        default_key_map = {
            "rotor_sleeve_depth": "mm_d_sleeve",
            "magnet_depth": "mm_d_pm",
            "magnet_pole_span_angle": "deg_alpha_rm",
            "stator_tooth_width": "mm_w_st",
            "stator_yoke_depth": "mm_d_sy",
            "stator_tooth_shoe_depth": "mm_d_sts",
            "stator_tooth_span_angle": "deg_alpha_st",
            "split_ratio_r_si_slash_r_so": "split_ratio",
        }
        key_map = key_map or default_key_map

        # 确保 parameter_dict_by_name 最新
        self.parameter_dict_by_name = self.get_parameter_dict_by_name()

        for raw_key, new_value in prev_params.items():
            attr_name = key_map.get(raw_key, raw_key)
            param = self.parameter_dict_by_name.get(raw_key)
            if param is None:
                candidate = getattr(self, attr_name, None)
                param = candidate if isinstance(candidate, Parameter) else None

            if param is None:
                logger.warning("prev_params key '%s' 未匹配到已知参数，已忽略。", raw_key)
                continue

            old_value = param.value
            if param.type != 'free':
                param.overwritten = True
                param.value = new_value
                logger.warning(
                    "覆盖非free参数 '%s' (type=%s) 从 %s 改为 %s，来源 key='%s'",
                    param.name, param.type, old_value, new_value, raw_key
                )
            else:
                param.value = new_value
                param.overwritten = False
                logger.info("设置 free 参数 '%s' 为 %s（来源 key='%s'）", param.name, new_value, raw_key)

        # 刷新参数映射，确保 calc 中引用的是最新对象
        self.parameter_dict_by_name = self.get_parameter_dict_by_name()

        # 更新 derived 参数，跳过被覆盖的
        for name, param in self.get_parameters_by_type('derived').items():
            if getattr(param, 'overwritten', False):
                logger.info("跳过已标记覆盖的导出参数 '%s'", param.name)
                continue
            if param.calc is not None and param.parameter_dict is not None:
                try:
                    param.parameter_dict = self.get_parameter_dict_by_name()
                    param.value = param.calc(param.parameter_dict)
                except Exception as e:
                    logger.warning("计算导出参数 '%s' 失败: %s", param.name, e)

        # 更新几何对象中的值
        if hasattr(self, "geometry") and hasattr(self.geometry, "machineGeometry"):
            for geo in self.geometry.machineGeometry.values():
                geo.update_from_GP()

    def update_geometric_parameters(self, x_denorm=None, x_denorm_dict=None):

        # 更新决策变量
        if x_denorm is not None:
            for i, param in enumerate(self.get_free_variables()):
                param.value = x_denorm[i]
        elif x_denorm_dict is not None:
            # 如果 x_denorm_dict 是 OrderedDict 或普通字典，直接使用
            for param in self.get_free_variables():
                if param.name in x_denorm_dict:
                    param.value = x_denorm_dict[param.name]

        # 更新 parameter_dict_by_name 以确保所有引用都是最新的
        self.parameter_dict_by_name = self.get_parameter_dict_by_name()

        # 同时刷新依赖于决策变量的导出参数。
        for i, param in enumerate(self.get_parameters_by_type('derived').values()):
            if param.calc is not None and param.parameter_dict is not None:
                param.value = param.calc(param.parameter_dict)

        # 当几何参数更新后，调用此方法来同步 machineGeometry 中的值
        # 更新 machineGeometry 中所有 Geometry 对象的参数值
        if hasattr(self, "geometry") and hasattr(self.geometry, "machineGeometry"):
            for geo_name, geo in self.geometry.machineGeometry.items():
                if geo is not None:
                    geo.update_from_GP()

    @staticmethod
    def get_pc_name():
        import platform, socket
        n1 = platform.node()
        # n2 = socket.gethostname()
        # n3 = os.environ["COMPUTERNAME"]
        # if n1 == n2 == n3:
        #     return n1
        # else:
        #     raise Exception(f"Computer names are not equal to each other. {n1,n2,n3}")
        return n1

    ''' 实用
    '''
    def learn_about_the_archive(self, prob, swarm_data, popsize, bool_plot_and_show=False, bool_more_info=False):
        """
        从群体数据中学习帕累托前沿，并基于支配排序和拥挤距离选择最优个体
        
        Args:
            prob: pygmo problem 对象
            swarm_data: 群体数据列表，每个元素是 [x_denorm..., f1, f2, f3]
            popsize: 种群大小
            bool_plot_and_show: 是否绘制并显示帕累托前沿
            bool_more_info: 是否返回额外信息
            
        Returns:
            如果 bool_more_info=False: 返回排序后的群体数据（前 popsize 个个体）
            如果 bool_more_info=True: 返回 (排序后的群体数据, 额外信息)
        """
        import pygmo as pg
        logger = logging.getLogger(__name__)
        
        number_of_chromosome = len(swarm_data)
        logger.info('Archive size: %d', number_of_chromosome)
        
        # 创建 archive 种群
        pop_archive = pg.population(prob, size=number_of_chromosome)
        for i in range(number_of_chromosome):
            pop_archive.set_xf(i, swarm_data[i][:-3], swarm_data[i][-3:])
        
        # 使用 sort_population_mo 对种群进行排序（基于支配排序和拥挤距离）
        sorted_index = pg.sort_population_mo(points=pop_archive.get_f())
        logger.info('Sorted by domination rank and crowding distance: %d', len(sorted_index))
        logger.debug('\t %s', sorted_index)
        
        # 获取非支配排序信息
        fits, vectors = pop_archive.get_f(), pop_archive.get_x()
        ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(fits)
        
        more_info = []
        ind1, ind2 = 0, 0
        for rank_minus_1, front in enumerate(ndf):
            ind2 += len(front)
            sorted_index_at_this_front = sorted_index[ind1:ind2]
            fits_at_this_front = [fits[point] for point in sorted_index_at_this_front]
            
            # Rank 1 Pareto Front
            if ind1 == 0:
                rank1_ParetoPoints = fits_at_this_front
                if len(front) < popsize:
                    logger.warning('There are not enough chromosomes (%d) belonging to domination rank 1 (the best Pareto front). Will use rank 2 or lower to reach popsize of %d.', len(front), popsize)
            
            # 计算拥挤距离
            if len(fits_at_this_front) >= 2:
                crwdst = pg.crowding_distance(fits_at_this_front)
            else:
                logger.warning('A non dominated front must contain at least two points: 1 detected.')
                crwdst = [999999]
            
            more_info.append((rank_minus_1+1, len(front), len(sorted_index_at_this_front)))
            ind1 = ind2
        
        # 获取排序后的向量和适应度值
        sorted_vectors = [vectors[index].tolist() for index in sorted_index]
        sorted_fits = [fits[index].tolist() for index in sorted_index]
        
        # 组合成完整的群体数据格式 [x_denorm..., f1, f2, f3]
        swarm_data_on_pareto_front = [design_parameters_denorm + fits 
                                      for design_parameters_denorm, fits in zip(sorted_vectors, sorted_fits)]
        
        # 只返回前 popsize 个个体（这些是具有高拥挤距离值的个体）
        swarm_data_selected = swarm_data_on_pareto_front[:popsize]
        
        if bool_plot_and_show:
            from pylab import plt
            # 可以在这里添加绘图代码
            pass
        
        if bool_more_info:
            return swarm_data_selected, more_info
        else:
            return swarm_data_selected
    
    def write_swarm_survivor(self, pop, counter_fitness_return):
        """
        将种群中的幸存者写入文件
        
        Args:
            pop: pygmo population 对象
            counter_fitness_return: 计数器值
        """

        survivor_file = self.path2SwarmData + '/swarm_survivor.txt'
        with open(survivor_file, 'a', encoding='utf-8') as f:
            f.write('\n---------%d\n' % counter_fitness_return)
            for el in zip(pop.get_x(), pop.get_f()):
                line = ','.join('%.4f' % x for x in el[0].tolist() + el[1].tolist())
                f.write(line + '\n')
    
    def get_bad_fintess_values(self, machine_type='IM', ref=False):
        # define bad values for different MOO objectives

        # return bad fitness values according to the moo.objective_functions defined in fea_config_dict
        # if 'moo.objective_functions' in self.fea_config_dict:
        #     objective_functions = self.fea_config_dict['moo.objective_functions']
        #     if 'f1' in objective_functions:
        #         return 0, 0, 99
        #     elif 'f2' in objective_functions:
        #         return 9999, 0, 999
        #     elif 'f3' in objective_functions:
        #         return 10000, 10, 1000

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
            if hasattr(self, 'winding') and hasattr(self.winding, 'EX'):
                return math.pi*(self.mm_r_ro.value*1e-3)**2 * (self.winding.EX['mm_stack_length_specified']*1e-3)
            else:
                return 0.0 # Or raise an error, depending on desired behavior
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
            return gravity * self.get_rotor_volume(stack_length=stack_length) * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 N/kg


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
        # 因为参数是在 __post_init__ 中动态添加的，不是 dataclass 字段，所以用fields是不全的。
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
            OrderedDict: 参数的有序字典，键为字段名
        """
        return OrderedDict(self.get_parameter_fields())
    
    def get_parameter_dict_by_name(self) -> OrderedDict:
        """
        获取所有参数的 OrderedDict，以参数名（Parameter.name）为键
        
        Returns:
            OrderedDict: 参数的有序字典，键为参数名（Parameter.name）
        """
        return OrderedDict((param.name, param) for param in self.get_parameter_fields().values())
    
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
        return (f"Modern_Machine_Designer(machine_class='{self.target.machine_class}', "
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
            
            # 保存 parameter_dict 的引用标记（用于重建 lambda 函数）
            # 注意：parameter_dict 本身不序列化，因为它包含循环引用
            # 在反序列化时，parameter_dict 会从对象的 get_parameter_dict_by_name() 方法获取
            if hasattr(param, 'parameter_dict') and param.parameter_dict is not None:
                param_dict['_has_parameter_dict'] = True
            
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
        full_dict['machine_class'] = self.target.machine_class
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
        if hasattr(self, 'winding') and hasattr(self.winding, 'wily') and self.winding.wily is not None:
            full_dict['wily'] = self.winding.wily.to_dict()
        
        # 4. 保存 EX 字典（激励参数）
        if hasattr(self, 'winding') and hasattr(self.winding, 'EX') and self.winding.EX is not None:
            full_dict['EX'] = self.winding.EX.copy()
        
        # 5. 保存 machineGeometry（完整版本）
        if hasattr(self, 'geometry') and hasattr(self.geometry, 'machineGeometry') and self.geometry.machineGeometry is not None:
            machine_geometry_dict = {}
            for key, geo in self.geometry.machineGeometry.items():
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
            'has_machineGeometry': hasattr(self, 'geometry') and hasattr(self.geometry, 'machineGeometry') and self.geometry.machineGeometry is not None,
            'has_wily': hasattr(self, 'winding') and hasattr(self.winding, 'wily') and self.winding.wily is not None,
            'has_EX': hasattr(self, 'winding') and hasattr(self.winding, 'EX') and self.winding.EX is not None,
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
        
        # New: Import necessary classes for deserialization
        from .parameter import Parameter
        from .winding import Winding
        from .geometry import Geometry
        from .winding import MachineWinding
        from .geometry import MachineGeometry

        # 新格式：从纯 JSON 数据恢复
        # 创建对象实例
        instance = cls.__new__(cls)
        
        # 恢复基本属性
        instance.name = data.get('name', 'SPMSM')
        instance.target.machine_class = data.get('machine_class', 'SPMSM')
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
            if not hasattr(instance, 'winding'):
                from machine_winding import MachineWinding
                instance.winding = MachineWinding()
            instance.winding.wily = Winding.from_dict(data['wily'])
        
        if 'EX' in data and data['EX'] is not None:
            if not hasattr(instance, 'winding'):
                from machine_winding import MachineWinding
                instance.winding = MachineWinding()
            instance.winding.EX = data['EX'].copy()
        
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
        
        # 更新所有参数的 parameter_dict 引用
        parameter_dict = instance.get_parameter_dict_by_name()
        parameter_dict_by_name = instance.get_parameter_dict_by_name()
        for field_name, param in instance.get_parameter_fields().items():
            # 对于需要参数名作为键的参数（derived 参数和 calc_bounds），使用 parameter_dict_by_name
            if param.calc is not None or param.calc_bounds is not None:
                param.parameter_dict = parameter_dict_by_name
            elif hasattr(param, 'parameter_dict'):
                # 其他情况也使用参数名作为键的字典
                param.parameter_dict = parameter_dict
        
        # 恢复保存的参数值（如果 __post_init__ 覆盖了它们）
        for field_name, saved_value in saved_param_values.items():
            param = getattr(instance, field_name, None)
            if param is not None and param.value != saved_value:
                param.value = saved_value
        
        # 更新所有 derived 参数的值（基于恢复的 free 参数值）
        parameter_dict_by_name = instance.get_parameter_dict_by_name()
        for field_name, param in instance.get_parameter_fields().items():
            if param.type == 'derived' and param.calc is not None:
                param.parameter_dict = parameter_dict_by_name
                try:
                    param.value = param.calc(parameter_dict_by_name)
                except Exception as e:
                    print(f"Warning: Failed to recalculate derived parameter '{field_name}': {e}")
        
        # 恢复 machineGeometry（如果 JSON 中有保存）
        if 'machineGeometry' in data and data['machineGeometry'] is not None:
            # 注意：machineGeometry 中的 visualization_points 会被保留
            if not hasattr(instance, 'geometry'):
                from machine_geometry import MachineGeometry
                instance.geometry = MachineGeometry()
            
            for key, geo_data in data['machineGeometry'].items():
                if geo_data is None:
                    if hasattr(instance.geometry, 'machineGeometry'):
                        instance.geometry.machineGeometry[key] = None
                elif isinstance(geo_data, dict):
                    # 只有当 machineGeometry 中已经存在该几何体时才恢复
                    if hasattr(instance.geometry, 'machineGeometry') and key in instance.geometry.machineGeometry:
                        geo = instance.geometry.machineGeometry[key]
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
        # New: Import necessary classes for deserialization
        from .parameter import Parameter
        from .winding import Winding
        from .geometry import Geometry
        from .winding import MachineWinding
        from .geometry import MachineGeometry

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
        instance.target.machine_class = data.get('machine_class', field_defaults.get('machine_class', cls.target.machine_class))
        instance.bool_PermanentMagnet = data.get('bool_PermanentMagnet', field_defaults.get('bool_PermanentMagnet', True))
        instance.bool_StatorSlotClosed = data.get('bool_StatorSlotClosed', field_defaults.get('bool_StatorSlotClosed', False))
        instance.bool_RotorNotched = data.get('bool_RotorNotched', field_defaults.get('bool_RotorNotched', True))
        
        # 设置所有参数字段
        for field_name, param in param_dict.items():
            setattr(instance, field_name, param)
        
        if 'winding' in data and data['winding'] is not None:
            if not hasattr(instance, 'winding'):
                from machine_winding import MachineWinding
                instance.winding = MachineWinding()
            instance.winding.wily = Winding.from_dict(data['winding'])
        
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


    @staticmethod
    def remove_jfiles_folders(root_dir):
        """
        Recursively remove all folders whose names end with 'jfiles' in the given directory.

        Args:
            root_dir (str): The root directory to start searching from.
        """
        import os
        import shutil

        for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
            for dirname in dirnames:
                if dirname.endswith('jfiles'):
                    folder_path = os.path.join(dirpath, dirname)
                    try:
                        shutil.rmtree(folder_path)
                        print(f"Deleted folder: {folder_path}")
                    except Exception as e:
                        print(f"Failed to delete {folder_path}: {e}")


    @staticmethod
    def myLogger(dir_log, prefix='/default_prefix_'):
        import logging, datetime, os
        """创建日志记录器"""
        # Disable matplotlib DEBUG logging to reduce log noise
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
        
        logger = logging.getLogger()
        if not len(logger.handlers):
            logger.setLevel(logging.DEBUG)
            now = datetime.datetime.now()

            if not os.path.isdir(dir_log):
                os.makedirs(dir_log)

            # create a file handler
            handler = logging.FileHandler(dir_log + prefix + '-' + now.strftime("%Y-%m-%d") + '.log')
            handler.setLevel(logging.DEBUG)

            # create a logging format
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)

            # add the handlers to the logger
            logger.addHandler(handler)
        return logger



    def draw_individual_from_swarm(self, index):
        """
        Draws the geometry for an individual from SwarmData.json, specified by index.

        Args:
            index (int): The index of the individual (e.g., 3383 for 'ind3383').
        """
        import json
        import os

        # Path to SwarmData.json
        swarm_data_path = os.path.join(self.path2SwarmData, "SwarmData.json")

        # Read the SwarmData.json file
        with open(swarm_data_path, "r", encoding="utf-8") as f:
            swarm_data = json.load(f)

        # Build the individual key
        suffix = f"ind{index}"
        # Find the key that ends with e.g. 'ind3383'
        target_key = None
        for key in swarm_data.keys():
            if key.endswith(suffix):
                target_key = key
                break

        if target_key is None:
            raise ValueError(f"Individual {suffix} not found in SwarmData.json.")

        # Extract and draw
        x_denorm_dict_raw = swarm_data[target_key]["x_denorm_dict"]
        # Decode x_denorm_dict if it's in py/reduce format
        x_denorm_dict = Swarm_Data_Analyzer.decode_py_reduce_ordered_dict(x_denorm_dict_raw)
        self.show_geometry(x_denorm_dict=x_denorm_dict)




class Swarm_Data_Analyzer(object):
    """
    分析群体优化数据的类，用于从 SwarmData.json 文件中读取和分析所有个体的设计参数和性能指标。
    """
    
    @staticmethod
    def decode_py_reduce_ordered_dict(x_denorm_dict_raw):
        """
        解码 jsonpickle 序列化的 OrderedDict (py/reduce 格式) 或普通字典
        
        Args:
            x_denorm_dict_raw: 包含 py/reduce 格式的字典或普通字典
            
        Returns:
            OrderedDict: 解码后的有序字典
        """
        # 如果是空字典，直接返回
        if not isinstance(x_denorm_dict_raw, dict):
            return OrderedDict()
        
        if len(x_denorm_dict_raw) == 0:
            return OrderedDict()
        
        # 如果是普通字典（不是 py/reduce 格式），直接转换为 OrderedDict
        if 'py/reduce' not in x_denorm_dict_raw:
            return OrderedDict(x_denorm_dict_raw)
        
        # 处理 py/reduce 格式
        try:
            # py/reduce 格式: [type_info, tuple_info, None, None, data_tuple]
            reduce_data = x_denorm_dict_raw['py/reduce']
            if len(reduce_data) >= 5 and isinstance(reduce_data[4], dict) and 'py/tuple' in reduce_data[4]:
                tuples = reduce_data[4]['py/tuple']
                result = OrderedDict()
                for item in tuples:
                    if isinstance(item, dict) and 'py/tuple' in item:
                        key_value_pair = item['py/tuple']
                        if len(key_value_pair) >= 2:
                            key = key_value_pair[0]
                            value = key_value_pair[1]
                            result[key] = value
                return result
        except (KeyError, IndexError, TypeError) as e:
            logger = logging.getLogger(__name__)
            logger.warning(f'Failed to decode py/reduce OrderedDict: {e}')
            return OrderedDict()
        
        return OrderedDict()
    
    def __init__(self, fname, desired_x_denorm_dict=None, bool_filter_pareto_front=False):
        """
        初始化群体数据分析器
        
        Args:
            fname: SwarmData.json 文件路径
            desired_x_denorm_dict: 期望的设计参数字典（用于排序），如果为 None 则使用所有参数
            bool_filter_pareto_front: 是否过滤帕累托前沿
        """
        # logger = logging.getLogger(__name__)
        # logger.info('Swarm_Data_Analyzer: %s', fname)
        
        if not os.path.exists(fname):
            self.number_of_chromosome = 0
            self.swarm_data_xf = None
            self.swarm_data_as_dict = {}
            self.swarm_data_project_names = []
            return
        
        ''' 1. Load json file
        '''
        # print(f'[Swarm_Data_Analyzer] read in {fname=}')
        with open(fname, 'r', encoding='utf-8') as f:
            try:
                # 尝试使用 jsonpickle 解码（如果文件是用 jsonpickle 保存的）
                swarm_data_as_dict = jsonpickle.decode(f.read())
            except Exception:
                # 如果 jsonpickle 解码失败，尝试使用标准 json 加载
                f.seek(0)  # 重置文件指针
                swarm_data_as_dict = json.load(f)
        
        if bool_filter_pareto_front:
            swarm_data_as_dict = self.filter_data(swarm_data_as_dict, 'Geometric parameters', 'split_ratio', 'bigger', 0.45)
        
        # for el in swarm_data_as_dict.keys():
        #     print(el)
        # print(len(swarm_data_as_dict.keys()))
        # quit()

        self.number_of_chromosome = len(swarm_data_as_dict)
        self.swarm_data_as_dict = swarm_data_as_dict

        ''' 2. Extract x_denorm_dict and performance metrics for each individual
        '''
        def sort_as_desired(x_denorm_dict, desired_x_denorm_dict=None):
            """根据 desired_x_denorm_dict 的顺序提取参数值"""
            if desired_x_denorm_dict is None:
                # 如果没有指定顺序，按字典的原始顺序返回所有值
                return list(x_denorm_dict.values())
            else:
                try:
                    # 按照 desired_x_denorm_dict 的键顺序提取值
                    return [x_denorm_dict[key] for key in desired_x_denorm_dict.keys()]
                except KeyError as e:
                    logger = logging.getLogger(__name__)
                    logger.warning(f'Error: some geometric parameters are renamed. Missing key: {e}')
                    raise KeyError
                    # 返回所有可用的值
                    return list(x_denorm_dict.values())

        # 处理每个个体
        self.swarm_data_xf = []
        for key, individual_data in swarm_data_as_dict.items():
            # 解码 x_denorm_dict
            x_denorm_dict_raw = individual_data.get('x_denorm_dict', {})
            x_denorm_dict = self.decode_py_reduce_ordered_dict(x_denorm_dict_raw)

            # 打印当前与期望的key情况
            logger = logging.getLogger(__name__)
            # logger.debug(f'[BEFORE] x_denorm_dict keys: {list(x_denorm_dict.keys())}')
            if desired_x_denorm_dict is not None:
                # logger.debug(f'[BEFORE] desired_x_denorm_dict keys: {list(desired_x_denorm_dict.keys())}')

                # 删除x_denorm_dict中多余的key
                keys_in_x = set(x_denorm_dict.keys())
                keys_in_desired = set(desired_x_denorm_dict.keys())

                # 1. 删去x_denorm_dict中多余的键
                redundant_keys = keys_in_x - keys_in_desired
                if redundant_keys:
                    logger.warning(f"x_denorm_dict has redundant keys {redundant_keys}, removing them.")
                for k in redundant_keys:
                    x_denorm_dict.pop(k, None)

                # 2. 补充x_denorm_dict中缺少的键
                missing_keys = keys_in_desired - keys_in_x
                if missing_keys:
                    logger.warning(f"x_denorm_dict missing keys {missing_keys}. Try to fill from other individuals.")
                for k in missing_keys:
                    # 改为在builtins.ad中的成员变量中找缺失的键
                    import builtins
                    value_to_add = None
                    ad = getattr(builtins, 'ad', None)
                    if ad is not None and hasattr(ad, 'x_denorm_dict'):
                        builtins_x_denorm_dict = getattr(ad, 'x_denorm_dict', {})
                        if k in builtins_x_denorm_dict:
                            value_to_add = builtins_x_denorm_dict[k]
                    if value_to_add is None and ad is not None and hasattr(ad, 'get_free_variables_as_dict'):
                        # 尝试通过接口函数获取
                        try:
                            from_free_vars = ad.get_free_variables_as_dict()
                            if k in from_free_vars:
                                value_to_add = from_free_vars[k]
                        except Exception:
                            pass
                    if value_to_add is None:
                        logger.warning(f"Could not find value for {k} in builtins.ad, will use current desired_x_denorm_dict value instead")
                        value_to_add = desired_x_denorm_dict[k]
                    x_denorm_dict[k] = value_to_add

                # 最后校正顺序
                x_denorm_dict = {k: x_denorm_dict[k] for k in desired_x_denorm_dict.keys()}

            # 如果 x_denorm_dict 为空，尝试从 desired_x_denorm_dict 获取默认值
            if not x_denorm_dict and desired_x_denorm_dict is not None:
                # 如果字典为空，使用 desired_x_denorm_dict 的当前值作为占位符
                # 这通常发生在数据保存时 x_denorm_dict 没有被正确序列化
                logger.warning(f'x_denorm_dict is empty for {key}, using current parameter values as placeholder')
                x_denorm_dict = desired_x_denorm_dict.copy()
            
            # 提取设计参数值
            x_denorm = sort_as_desired(x_denorm_dict, desired_x_denorm_dict)
            
            # 提取性能指标 f1, f2, f3
            # 处理 f1, f2, f3 可能是空字典 {} 的情况
            def get_value_or_default(value, default=0.0):
                """如果值是空字典，返回默认值；否则返回实际值"""
                if isinstance(value, dict) and len(value) == 0:
                    return default
                elif value is None:
                    return default
                else:
                    return float(value) if value != {} else default
            
            f1 = get_value_or_default(individual_data.get('f1'), 0.0)
            f2 = get_value_or_default(individual_data.get('f2'), 0.0)
            f3 = get_value_or_default(individual_data.get('f3'), 0.0)
            
            # 组合成 [x_denorm..., f1, f2, f3]
            # 确保 x_denorm 的长度与 desired_x_denorm_dict 一致
            if desired_x_denorm_dict is not None and len(x_denorm) != len(desired_x_denorm_dict):
                logger = logging.getLogger(__name__)
                logger.warning(f'x_denorm length mismatch for {key}: expected {len(desired_x_denorm_dict)}, got {len(x_denorm)}')
                # 如果长度不匹配，使用 desired_x_denorm_dict 的当前值填充
                if len(x_denorm) == 0:
                    x_denorm = list(desired_x_denorm_dict.values())
                elif len(x_denorm) < len(desired_x_denorm_dict):
                    # 如果长度不足，用当前值填充
                    x_denorm = x_denorm + list(desired_x_denorm_dict.values())[len(x_denorm):]
                else:
                    # 如果长度超出，截断
                    x_denorm = x_denorm[:len(desired_x_denorm_dict)]
            
            self.swarm_data_xf.append(x_denorm + [f1, f2, f3])
        
        if len(self.swarm_data_xf) > 0:
            self.number_of_free_variables = len(self.swarm_data_xf[0]) - 3
        else:
            self.number_of_free_variables = 0

        # DEBUG
        # print('[Swarm_Data_Analyzer]')
        # for ind, xf in enumerate(self.swarm_data_xf):
        #     print(f'{ind:04d}', ',\t'.join([f'{el:.2f}' for el in xf]))

        ''' 3. Get the list of other attribute by individuals (not needed for optimization)
        '''
        self.swarm_data_project_names = self.get_metric_of_the_whole_swarm('project_name')
    
    def filter_data(self, data, param_type, filter_key, direction, filter_value):
        """
        过滤数据（保留旧接口以兼容）
        """
        filtered_data = {}
        for key in data.keys():
            individual_data = data[key]
            # 尝试从 x_denorm_dict 中获取参数值
            x_denorm_dict_raw = individual_data.get('x_denorm_dict', {})
            x_denorm_dict = self.decode_py_reduce_ordered_dict(x_denorm_dict_raw)
            
            if filter_key in x_denorm_dict:
                value = x_denorm_dict[filter_key]
                flag = False
                if direction == "bigger" and value > filter_value:
                    filtered_data[key] = individual_data
                    flag = True
                elif direction == "smaller" and value <= filter_value:
                    filtered_data[key] = individual_data
                    flag = True
        return filtered_data

    @staticmethod
    def decode(d):
        """解码函数（保留旧接口以兼容）"""
        if isinstance(d, dict):
            return list(d.values())[0] if len(d) > 0 else d
        return d

    def get_metric_of_the_whole_swarm(self, metric):
        """
        获取整个群体中所有个体的某个指标
        
        Args:
            metric: 指标名称（如 'f1', 'f2', 'Cost', 'TRV' 等）
            
        Returns:
            list: 所有个体的该指标值列表
        """
        result = []
        for key, individual_data in self.swarm_data_as_dict.items():
            value = individual_data.get(metric, None)
            if value is not None:
                result.append(value)
            else:
                logger = logging.getLogger(__name__)
                logger.warning(f'Metric "{metric}" not found in individual {key}')
                result.append(0.0)  # 默认值
        return result
    def prepare_data_for_post_processing(self):
        """
        准备后处理数据，收集所有个体的各种性能指标
        """
        # 基本性能指标
        self.FRW = self.get_metric_of_the_whole_swarm('FRW')
        self.Em = self.get_metric_of_the_whole_swarm('normalized_force_error_magnitude')
        self.Ea = self.get_metric_of_the_whole_swarm('force_error_angle')
        self.Tripple = self.get_metric_of_the_whole_swarm('normalized_torque_ripple')
        self.RatedStkLen = self.get_metric_of_the_whole_swarm('rated_stack_length_mm')

        # 从 swarm_data_xf 提取 f1, f2, f3
        if self.swarm_data_xf is not None and len(self.swarm_data_xf) > 0:
            self.f1 = [raw[-3] for raw in self.swarm_data_xf]
            self.f2 = [raw[-2] for raw in self.swarm_data_xf]
            self.f3 = [raw[-1] for raw in self.swarm_data_xf]
        else:
            self.f1 = []
            self.f2 = []
            self.f3 = []

        # 功率因数和效率相关
        self.PowerFactor = self.get_metric_of_the_whole_swarm('power_factor')
        
        # Cost 和 TRV
        try:
            self.Cost = self.get_metric_of_the_whole_swarm('Cost')
        except:
            import numpy as np
            self.Cost = np.array(self.get_metric_of_the_whole_swarm('f1'))
        
        try:
            self.TRV = self.get_metric_of_the_whole_swarm('TRV')
        except:
            import numpy as np
            self.TRV = -np.array(self.get_metric_of_the_whole_swarm('f1'))
        
        # 效率
        try:
            self.RatedEfficiency = self.get_metric_of_the_whole_swarm('RatedEfficiency')
        except:
            import numpy as np
            self.RatedEfficiency = -np.array(self.get_metric_of_the_whole_swarm('f2'))
        
        self.TorqueRipple = self.get_metric_of_the_whole_swarm('normalized_torque_ripple')
        
        # 其他性能指标
        self.torque_average = self.get_metric_of_the_whole_swarm('torque_average')
        self.ss_avg_force_magnitude = self.get_metric_of_the_whole_swarm('ss_avg_force_magnitude')
        self.normalized_force_error_magnitude = self.get_metric_of_the_whole_swarm('normalized_force_error_magnitude')
        self.force_error_angle = self.get_metric_of_the_whole_swarm('force_error_angle')

        # 损耗相关
        self.l_rated_total_loss = self.get_metric_of_the_whole_swarm('rated_total_loss')
        self.l_rated_stator_copper_loss_along_stack = self.get_metric_of_the_whole_swarm('rated_stator_copper_loss_along_stack')
        self.l_rated_rotor_copper_loss_along_stack = self.get_metric_of_the_whole_swarm('rated_rotor_copper_loss_along_stack')
        self.l_stator_copper_loss_in_end_turn = self.get_metric_of_the_whole_swarm('stator_copper_loss_in_end_turn')
        self.l_rotor_copper_loss_in_end_turn = self.get_metric_of_the_whole_swarm('rotor_copper_loss_in_end_turn')
        self.l_rated_iron_loss = self.get_metric_of_the_whole_swarm('rated_iron_loss')
        self.l_rated_windage_loss = self.get_metric_of_the_whole_swarm('rated_windage_loss')
        self.l_rated_magnet_Joule_loss = self.get_metric_of_the_whole_swarm('rated_magnet_Joule_loss')
        self.l_rated_stack_length = self.get_metric_of_the_whole_swarm('rated_stack_length_mm')
        
        # 成本和重量相关
        try:
            self.Cost_Fe = self.get_metric_of_the_whole_swarm('Cost_Fe')
            self.Cost_Cu = self.get_metric_of_the_whole_swarm('Cost_Cu')
            self.Cost_PM = self.get_metric_of_the_whole_swarm('Cost_PM')
        except:
            self.Cost_Fe = []
            self.Cost_Cu = []
            self.Cost_PM = []
        
        try:
            self.rotor_weight = self.get_metric_of_the_whole_swarm('rotor_weight')
        except:
            self.rotor_weight = []

class swarm_data_container(object):
    """
    群体数据容器类，用于从原始文本数据或 JSON 文件读取和分析群体优化数据。
    支持两种数据源：
    1. 原始文本格式（swarm_data_raw）
    2. JSON 格式（swarm_data_json 或从 JSON 文件读取）
    """
    
    def __init__(self, swarm_data_raw=None, fea_config_dict=None, swarm_data_json=None, swarm_data_json_file_path=None):
        """
        初始化群体数据容器
        
        Args:
            swarm_data_raw: 原始文本数据（列表格式，向后兼容）
            fea_config_dict: FEA 配置字典
            swarm_data_json: JSON 格式的数据字典（如果提供，将优先使用）
            swarm_data_json_file_path: JSON 文件路径（如果提供，将从文件读取）
        """
        self.swarm_data_raw = swarm_data_raw
        self.fea_config_dict = fea_config_dict or {}
        
        # 如果提供了 JSON 文件路径，从文件读取
        if swarm_data_json_file_path is not None and os.path.exists(swarm_data_json_file_path):
            logger = logging.getLogger(__name__)
            logger.info(f'Loading swarm data from JSON file: {swarm_data_json_file_path}')
            with open(swarm_data_json_file_path, 'r', encoding='utf-8') as f:
                swarm_data_json = json.load(f)
        
        # 如果提供了 JSON 数据，使用 JSON 数据源
        if swarm_data_json is not None:
            self._load_from_json(swarm_data_json)
        elif swarm_data_raw is not None:
            self._load_from_raw(swarm_data_raw)
        else:
            logger = logging.getLogger(__name__)
            logger.warning('No data source provided (neither swarm_data_raw nor swarm_data_json/swarm_data_json_file_path)')
            self._initialize_empty()
    
    def _initialize_empty(self):
        """初始化空的数据结构"""
        self.swarm_data_xf = []
        self.project_names = []
        self.machine_data = []
        self.rated_data = []
        self.Trip = []
        self.FRW = []
        self.Em = []
        self.Ea = []
        self.RatedVol = []
        self.RatedWeight = []
        self.RatedStkLen = []
        self.deg_alpha_st = []
        self.mm_w_st = []
        self.mm_r_si = []
        self.number_of_free_variables = 0
        
        # 性能指标列表
        self.l_OA = []
        self.l_OB = []
        self.l_OC = []
        self.l_design_parameters = []
        self.l_power_factor = []
        self.l_efficiency = []
        self.l_torque_average = []
        self.l_normalized_torque_ripple = []
        self.l_ss_avg_force_magnitude = []
        self.l_normalized_force_error_magnitude = []
        self.l_force_error_angle = []
        self.l_rated_shaft_power = []
        self.l_rated_efficiency = []
        self.l_rated_total_loss = []
        self.l_rated_stator_copper_loss_along_stack = []
        self.l_rated_rotor_copper_loss_along_stack = []
        self.l_stator_copper_loss_in_end_turn = []
        self.l_rotor_copper_loss_in_end_turn = []
        self.l_rated_iron_loss = []
        self.l_rated_windage_loss = []
        self.l_rated_rotor_volume = []
        self.l_rated_rotor_weight = []
        self.l_rated_stack_length = []
        self.l_original_stack_length = []
        self.l_original_rotor_weight = []
        self.l_TRV = []
        self.l_FRW = []
    
    def _load_from_json(self, swarm_data_json):
        """
        从 JSON 格式的数据加载
        
        Args:
            swarm_data_json: JSON 格式的群体数据字典
        """
        logger = logging.getLogger(__name__)
        logger.info(f'Loading {len(swarm_data_json)} individuals from JSON data')
        
        # 初始化数据结构
        self._initialize_empty()
        
        # 使用 Swarm_Data_Analyzer 的解码方法
        decoder = Swarm_Data_Analyzer.__new__(Swarm_Data_Analyzer)  # 创建实例但不调用 __init__
        
        for key, individual_data in swarm_data_json.items():
            # 解码 x_denorm_dict
            x_denorm_dict_raw = individual_data.get('x_denorm_dict', {})
            x_denorm_dict = Swarm_Data_Analyzer.decode_py_reduce_ordered_dict(x_denorm_dict_raw)
            
            # 提取设计参数值（按顺序）
            x_denorm = list(x_denorm_dict.values())
            
            # 提取性能指标
            f1 = individual_data.get('f1', 0.0)
            f2 = individual_data.get('f2', 0.0)
            f3 = individual_data.get('f3', 0.0)
            
            # 组合成 [x_denorm..., f1, f2, f3]
            self.swarm_data_xf.append(x_denorm + [f1, f2, f3])
            
            # 提取项目名称
            project_name = individual_data.get('project_name', key)
            self.project_names.append(project_name)
            
            # 构建 machine_data（模拟原始格式）
            # [power_factor, efficiency, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
            machine_data_item = [
                individual_data.get('power_factor', 0.0),
                individual_data.get('RatedEfficiency', -f2 if f2 < 0 else 0.0),  # 效率可能是负的 f2
                individual_data.get('torque_average', 0.0),
                individual_data.get('normalized_torque_ripple', 0.0),
                individual_data.get('ss_avg_force_magnitude', 0.0),
                individual_data.get('normalized_force_error_magnitude', 0.0),
                individual_data.get('force_error_angle', 0.0),
            ]
            self.machine_data.append(machine_data_item)
            
            # 构建 rated_data（模拟原始格式）
            rated_data_item = [
                individual_data.get('rated_shaft_power', 0.0),
                individual_data.get('RatedEfficiency', -f2 if f2 < 0 else 0.0),
                individual_data.get('rated_total_loss', 0.0),
                individual_data.get('rated_stator_copper_loss_along_stack', 0.0),
                individual_data.get('rated_rotor_copper_loss_along_stack', 0.0),
                individual_data.get('stator_copper_loss_in_end_turn', 0.0),
                individual_data.get('rotor_copper_loss_in_end_turn', 0),
                individual_data.get('rated_iron_loss', 0.0),
                individual_data.get('rated_windage_loss', 0.0),
                individual_data.get('rated_rotor_volume', 0.0),
                individual_data.get('rated_stack_length_mm', 0.0),
                individual_data.get('original_stack_length', 0.0),
            ]
            self.rated_data.append(rated_data_item)
            
            # 提取其他指标
            self.Trip.append(individual_data.get('normalized_torque_ripple', 0.0))
            
            # 计算 FRW（如果需要）
            ss_avg_force = individual_data.get('ss_avg_force_magnitude', 0.0)
            rotor_weight = individual_data.get('rotor_weight', 0.0)
            if rotor_weight > 0:
                individual_FRW = ss_avg_force / rotor_weight
            else:
                individual_FRW = 0.0
            self.FRW.append(individual_FRW)
            
            self.Em.append(individual_data.get('normalized_force_error_magnitude', 0.0))
            self.Ea.append(individual_data.get('force_error_angle', 0.0))
            
            # 体积和重量
            rated_rotor_volume = individual_data.get('rated_rotor_volume', 0.0)
            if rated_rotor_volume == 0.0:
                # 如果没有直接提供，尝试从其他数据计算
                rated_rotor_volume = 0.0  # 需要更多信息才能计算
            self.RatedVol.append(rated_rotor_volume)
            
            # 转子重量（密度 8050 kg/m^3，转换为 N）
            if rated_rotor_volume > 0:
                individual_rated_rotor_weight = rated_rotor_volume * 8050 * 9.8  # N
            else:
                individual_rated_rotor_weight = individual_data.get('rotor_weight', 0.0)
            self.RatedWeight.append(individual_rated_rotor_weight)
            
            rated_stack_length = individual_data.get('rated_stack_length_mm', 0.0)
            self.RatedStkLen.append(rated_stack_length)
            
            # 提取特定设计参数（如果存在）
            if 'magnet_depth' in x_denorm_dict:
                # 尝试从 x_denorm_dict 中提取特定参数
                pass  # 这些参数可能不在 x_denorm_dict 中
        
        # 设置自由变量数量
        if len(self.swarm_data_xf) > 0:
            self.number_of_free_variables = len(self.swarm_data_xf[0]) - 3
        else:
            self.number_of_free_variables = 0
        
        # 提取所有性能指标列表
        self._extract_performance_lists()
        
        logger.info(f'Loaded {len(self.swarm_data_xf)} individuals, {self.number_of_free_variables} free variables')
    
    def _extract_performance_lists(self):
        """从 machine_data 和 rated_data 提取所有性能指标列表"""
        # 从 swarm_data_xf 提取目标函数值
        if len(self.swarm_data_xf) > 0:
            self.l_OA = [raw[-3] for raw in self.swarm_data_xf]  # f1
            self.l_OB = [raw[-2] for raw in self.swarm_data_xf]  # f2
            self.l_OC = [raw[-1] for raw in self.swarm_data_xf]  # f3
            self.l_design_parameters = [raw[:-3] for raw in self.swarm_data_xf]
        
        # 从 machine_data 提取
        if len(self.machine_data) > 0:
            self.l_power_factor = [raw[0] for raw in self.machine_data]
            self.l_efficiency = [raw[1] for raw in self.machine_data]
            self.l_torque_average = [raw[2] for raw in self.machine_data]
            self.l_normalized_torque_ripple = [raw[3] for raw in self.machine_data]
            self.l_ss_avg_force_magnitude = [raw[4] for raw in self.machine_data]
            self.l_normalized_force_error_magnitude = [raw[5] for raw in self.machine_data]
            self.l_force_error_angle = [raw[6] for raw in self.machine_data]
        
        # 从 rated_data 提取
        if len(self.rated_data) > 0:
            self.l_rated_shaft_power = [raw[0] for raw in self.rated_data]
            self.l_rated_efficiency = [raw[1] for raw in self.rated_data]
            self.l_rated_total_loss = [raw[2] for raw in self.rated_data]
            self.l_rated_stator_copper_loss_along_stack = [raw[3] for raw in self.rated_data]
            self.l_rated_rotor_copper_loss_along_stack = [raw[4] for raw in self.rated_data]
            self.l_stator_copper_loss_in_end_turn = [raw[5] for raw in self.rated_data]
            self.l_rotor_copper_loss_in_end_turn = [raw[6] for raw in self.rated_data]
            self.l_rated_iron_loss = [raw[7] for raw in self.rated_data]
            self.l_rated_windage_loss = [raw[8] for raw in self.rated_data]
            self.l_rated_rotor_volume = [raw[9] for raw in self.rated_data]
            self.l_rated_rotor_weight = [(V*8050*9.8) for V in self.l_rated_rotor_volume]  # N
            self.l_rated_stack_length = [raw[10] for raw in self.rated_data]
            self.l_original_stack_length = [raw[11] for raw in self.rated_data]
            self.l_original_rotor_weight = [weight/rated*ori if rated > 0 else 0 
                                            for weight, rated, ori in zip(self.l_rated_rotor_weight, 
                                                                          self.l_rated_stack_length, 
                                                                          self.l_original_stack_length)]
        
        # 计算 TRV 和 FRW
        import numpy as np
        required_torque = 50e3 / (30000/60*2*math.pi)  # TODO: 应该使用额定堆叠长度和平均转矩计算
        self.l_TRV = [required_torque/raw if raw > 0 else 0 for raw in self.l_rated_rotor_volume]
        self.l_FRW = [F/W if W > 0 else 0 for W, F in zip(self.l_original_rotor_weight, self.l_ss_avg_force_magnitude)]
    
    def _load_from_raw(self, swarm_data_raw):
        """
        从原始文本数据加载（保持向后兼容）
        
        Args:
            swarm_data_raw: 原始文本数据列表，每个元素是一个列表，包含：
                raw[0]: 索引或其他信息
                raw[1]: 项目名称
                raw[2]: 包含 f1, f2, f3 的字符串
                raw[3]: 机器数据（逗号分隔的字符串）
                raw[4]: 额定数据（逗号分隔的字符串）
                raw[5]: 设计参数（逗号分隔的字符串）
        """
        logger = logging.getLogger(__name__)
        logger.info(f'Loading {len(swarm_data_raw)} individuals from raw text data')
        
        # 初始化数据结构
        self._initialize_empty()
        
        self.deg_alpha_st = []
        self.mm_w_st = []
        self.mm_r_si = []
        
        # 处理每个个体
        for raw in swarm_data_raw:

                    # spmsm_template.design_parameters = [
                    #                                   0 spmsm_template.deg_alpha_st 
                    #                                   1 spmsm_template.deg_alpha_sto 
                    #                                   2 spmsm_template.mm_r_si      
                    #                                   3 spmsm_template.mm_d_sto      
                    #                                   4 spmsm_template.mm_d_sts      
                    #                                   5 spmsm_template.mm_d_st      
                    #                                   6 spmsm_template.mm_d_sy      
                    #                                   7 spmsm_template.mm_w_st      
                    #                                   8 spmsm_template.mm_r_st      
                    #                                   9 spmsm_template.mm_r_sf      
                    #                                  10 spmsm_template.mm_r_sb      
                    #                                  11 spmsm_template.Q            
                    #                                  12 spmsm_template.sleeve_length
                    #                                  13 spmsm_template.fixed_air_gap_length
                    #                                  14 spmsm_template.mm_d_pm      
                    #                                  15 spmsm_template.deg_alpha_rm 
                    #                                  16 spmsm_template.deg_alpha_rs 
                    #                                  17 spmsm_template.mm_d_ri      
                    #                                  18 spmsm_template.mm_r_ri      
                    #                                  19 spmsm_template.mm_d_rp      
                    #                                  20 spmsm_template.mm_d_rs      
                    #                                  21 spmsm_template.p
                    #                                  22 spmsm_template.s
                    #                                 ]
                    design_parameters_denorm = [float(x) for x in raw[5].split(',')]
                    self.deg_alpha_st.append(design_parameters_denorm[0] )
                    self.mm_w_st.append(     design_parameters_denorm[7] )
                    self.mm_r_si.append(   design_parameters_denorm[2])

                    loc1 = raw[2].find('f1')
                    loc2 = raw[2].find('f2')
                    loc3 = raw[2].find('f3')
                    f1 = float(raw[2][loc1+3:loc2-1])
                    f2 = float(raw[2][loc2+3:loc3-1])
                    f3 = float(raw[2][loc3+3:])

                    if len(design_parameters_denorm) > 20:
                        ''' 永磁电机 复古 '''

                        # 在 acmop 里，我们已经放弃了使用 bound_filter 的概念。
                        # x_denorm = self.get_x_denorm_from_design_parameters(design_parameters_denorm, bound_filter)

                        """ This is consistent with bopt-python """
                        # x_denorm = [None]*11
                        # x_denorm[0]  = design_parameters_denorm[0] # spmsm_template.deg_alpha_st 
                        # x_denorm[1]  = design_parameters_denorm[3] # spmsm_template.mm_d_sto         
                        # x_denorm[2]  = design_parameters_denorm[5] # spmsm_template.mm_d_st
                        # x_denorm[3]  = design_parameters_denorm[7] # spmsm_template.mm_w_st         
                        # x_denorm[4]  = design_parameters_denorm[12] # spmsm_template.sleeve_length   
                        # x_denorm[5]  = design_parameters_denorm[14] # spmsm_template.mm_d_pm         
                        # x_denorm[6]  = design_parameters_denorm[15] # spmsm_template.deg_alpha_rm    
                        # x_denorm[7]  = design_parameters_denorm[16] # spmsm_template.deg_alpha_rs    
                        # x_denorm[8]  = design_parameters_denorm[17] # spmsm_template.mm_d_ri         
                        # x_denorm[9]  = design_parameters_denorm[19] # spmsm_template.mm_d_rp         
                        # x_denorm[10] = design_parameters_denorm[20] # spmsm_template.mm_d_rs         

                        if False:
                            """ This is consistent with ACMOP """
                            x_denorm = [None]*11
                            x_denorm[0]  = design_parameters_denorm[0] # spmsm_template.deg_alpha_st 
                            x_denorm[1]  = design_parameters_denorm[3] # spmsm_template.mm_d_sto         
                            x_denorm[2]  = design_parameters_denorm[5] # spmsm_template.mm_d_st
                            x_denorm[3]  = sum([design_parameters_denorm[i] for i in (2,4,5,6)]) # outer_stator_radius mm_r_so
                            x_denorm[4]  = design_parameters_denorm[7] # spmsm_template.mm_w_st   
                            x_denorm[5]  = design_parameters_denorm[12] #            mm_d_sleeve
                            r_si = design_parameters_denorm[2] # 2 spmsm_template.mm_r_si      
                            try:
                                x_denorm[6]  = r_si /  x_denorm[3] # split_ratio     r_is_slash_r_os 
                            except ZeroDivisionError as e:
                                print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                                print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                                print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                                raise e
                            x_denorm[7]  = design_parameters_denorm[14] # spmsm_template.mm_d_pm      
                            x_denorm[8]  = design_parameters_denorm[17] # spmsm_template.mm_d_ri         
                            # childGP
                            x_denorm[9]  = design_parameters_denorm[15] # spmsm_template.deg_alpha_rm    
                            x_denorm[10]  = design_parameters_denorm[19] # spmsm_template.mm_d_rp         
                            # x_denorm[11]  = design_parameters_denorm[16] # spmsm_template.deg_alpha_rs    
                            # x_denorm[12] = design_parameters_denorm[20] # spmsm_template.mm_d_rs         

                            # DEBUG
                            # odict_keys(['deg_alpha_st', 'mm_d_sto', 'mm_d_st', 'mm_r_so', 'mm_w_st', 'mm_d_sleeve', 'split_ratio', 'mm_d_pm', 'mm_d_ri', 'deg_alpha_rm', 'mm_d_rp'])
                            # deg_alpha_st 11.1183
                            # mm_d_sto 1.50079
                            # mm_d_st 42.9701
                            # mm_r_so 16.099
                            # mm_w_st 5.89091
                            # mm_d_sleeve 5.19948
                            # split_ratio 44.9638
                            # mm_d_pm 44.9638
                            # mm_d_ri 3.67901
                            # deg_alpha_rm 5.19948
                            # mm_d_rp 0.0
                        else:
                            """ This is consistent with ACMOP """
                            x_denorm = [None]*11
                            x_denorm[0]  = design_parameters_denorm[0] # spmsm_template.deg_alpha_st 
                            x_denorm[1]  = design_parameters_denorm[3] # spmsm_template.mm_d_sto         
                            x_denorm[2]  = design_parameters_denorm[5] # spmsm_template.mm_d_st
                            x_denorm[3]  = sum([design_parameters_denorm[i] for i in (2,4,5,6)]) # outer_stator_radius mm_r_so
                            print(f'{x_denorm[1]=}, {x_denorm[0]/2=}')
                            x_denorm[4]  = design_parameters_denorm[7] # spmsm_template.mm_w_st   
                            x_denorm[5]  = design_parameters_denorm[12] #            mm_d_sleeve
                            r_si = design_parameters_denorm[2] # 2 spmsm_template.mm_r_si      
                            try:
                                x_denorm[6]  = r_si /  x_denorm[3] # split_ratio     r_is_slash_r_os 
                                # print(f'{x_denorm[6]=}, {r_si=}')
                            except ZeroDivisionError as e:
                                print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                                print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                                print('Error: You need to clean up the swarm_data.txt file. There is a design with zero element in design_parameters (which is intended with ACMOP).')
                                raise e
                            x_denorm[7]  = design_parameters_denorm[14] # spmsm_template.mm_d_pm      
                            x_denorm[8]  = design_parameters_denorm[17] # spmsm_template.mm_d_ri         
                            # childGP
                            x_denorm[9]  = design_parameters_denorm[15] # spmsm_template.deg_alpha_rm    
                            x_denorm[10]  = design_parameters_denorm[19] # spmsm_template.mm_d_rp         

                    else:
                        '''感应电机 复古 '''
                        raise Exception('not implemented')

                    # print(design_parameters_denorm, f1, f2, f3)
                    # THERE IS A BUT HERE: slot_tip_open_ratio is less than 0.2---Not possible
                        # free_variables[0]  = design_parameters[0] # spmsm_template.deg_alpha_st 
                        # free_variables[4]  = design_parameters[7] # spmsm_template.mm_w_st         
                        # free_variables[10] = sum([design_parameters[i] for i in (18,17,19)]) # spmsm_template.mm_r_ri + spmsm_template.mm_d_ri + spmsm_template.mm_d_rp
                        # self.deg_alpha_st.append(x_denorm[0] ) 
                        # self.mm_w_st.append(x_denorm[4] ) 
                        # self.mm_radius.append(x_denorm[10]) 

                    self.project_names.append(raw[1][:-1])
                    self.machine_data.append([float(x) for x in raw[3].split(',')])
                    self.rated_data.append(  [float(x) for x in raw[4].split(',')])

                    individual_Trip = [float(x) for x in raw[3].split(',')][3]
                    self.Trip.append(individual_Trip)

                    # Get FRW
                    individual_ss_avg_force_magnitude = [float(x) for x in raw[3].split(',')][4]
                    individual_Em                     = [float(x) for x in raw[3].split(',')][5]
                    individual_Ea                     = [float(x) for x in raw[3].split(',')][6]
                    individual_rated_rotor_volume     = [float(x) for x in raw[4].split(',')][9]
                    individual_rated_rotor_weight     = (individual_rated_rotor_volume*8050*9.8)
                    individual_rated_stack_length     = [float(x) for x in raw[4].split(',')][10]
                    individual_original_stack_length  = [float(x) for x in raw[4].split(',')][11]
                    individual_original_rotor_weight  = individual_rated_rotor_weight/individual_rated_stack_length*individual_original_stack_length
                    individual_FRW = individual_ss_avg_force_magnitude/individual_original_rotor_weight
                    self.FRW.append(individual_FRW)
                    self.Em.append(individual_Em)
                    self.Ea.append(individual_Ea)
                    self.RatedVol.append(individual_rated_rotor_volume)
                    self.RatedWeight.append(individual_rated_rotor_weight)
                    self.RatedStkLen.append(individual_rated_stack_length)

                    # Add FRW to xf (This will cause re-starting error)
                    # self.swarm_data_xf.append(x_denorm + [individual_FRW, f1, f2, f3])

                    self.swarm_data_xf.append(x_denorm + [f1, f2, f3])

        # 设置自由变量数量
        if len(self.swarm_data_xf) > 0:
            self.number_of_free_variables = len(self.swarm_data_xf[0]) - 3
        else:
            self.number_of_free_variables = 0
        
        logger.info(f'Loaded {len(self.swarm_data_xf)} individuals from raw data, {self.number_of_free_variables} free variables')
        
        # 提取所有性能指标列表
        self._extract_performance_lists()

        self.l_OA = [raw[-3] for raw in self.swarm_data_xf]
        self.l_OB = [raw[-2] for raw in self.swarm_data_xf]
        self.l_OC = [raw[-1] for raw in self.swarm_data_xf]
        self.l_design_parameters = [raw[:-3] for raw in self.swarm_data_xf]

        # [power_factor, efficiency, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
        self.l_power_factor                     = [raw[0] for raw in self.machine_data]
        self.l_efficiency                       = [raw[1] for raw in self.machine_data]
        self.l_torque_average                   = [raw[2] for raw in self.machine_data]
        self.l_normalized_torque_ripple         = [raw[3] for raw in self.machine_data]
        self.l_ss_avg_force_magnitude           = [raw[4] for raw in self.machine_data]
        self.l_normalized_force_error_magnitude = [raw[5] for raw in self.machine_data]
        self.l_force_error_angle                = [raw[6] for raw in self.machine_data]

        self.l_rated_shaft_power                    = [raw[0] for raw in self.rated_data]
        self.l_rated_efficiency                     = [raw[1] for raw in self.rated_data]
        self.l_rated_total_loss                     = [raw[2] for raw in self.rated_data]
        self.l_rated_stator_copper_loss_along_stack = [raw[3] for raw in self.rated_data]
        self.l_rated_rotor_copper_loss_along_stack  = [raw[4] for raw in self.rated_data]
        self.l_stator_copper_loss_in_end_turn       = [raw[5] for raw in self.rated_data]
        self.l_rotor_copper_loss_in_end_turn        = [raw[6] for raw in self.rated_data]
        self.l_rated_iron_loss                      = [raw[7] for raw in self.rated_data]
        self.l_rated_windage_loss                   = [raw[8] for raw in self.rated_data]
        self.l_rated_rotor_volume                   = [raw[9] for raw in self.rated_data]
        self.l_rated_rotor_weight                   = [(V*8050*9.8) for V in self.l_rated_rotor_volume] # density of rotor is estimated to be that of steraw of 8050 g/cm^3
        self.l_rated_stack_length                   = [raw[10] for raw in self.rated_data] # new!
        self.l_original_stack_length                = [raw[11] for raw in self.rated_data] # new!
        self.l_original_rotor_weight                = [weight/rated*ori for weight, rated, ori in zip(self.l_rated_rotor_weight, self.l_rated_stack_length, self.l_original_stack_length)]

        # TODO: change to EX['mec_power'] and EX['the_speed']
        required_torque = 50e3 / (30000/60*2*math.pi)         # TODO: should use rated stack length and torque average to compute this
        self.l_TRV = [required_torque/raw for raw in self.l_rated_rotor_volume]
        self.l_FRW = [F/W for W, F in zip(self.l_original_rotor_weight, self.l_ss_avg_force_magnitude)] # FRW

    def get_list_y_data(self):

        list_y_data = [ self.l_TRV, ##self.l_rated_stack_length,
                        [100*raw for raw in self.l_OB], 
                        self.l_force_error_angle,
                        ]
        return list_y_data

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Utility
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~


    def sensitivity_bar_charts(self):
        number_of_variant = self.fea_config_dict['local_sensitivity_analysis_number_of_variants'] + 1
        number_of_free_variables = self.number_of_free_variables

        from pylab import subplots, mpl, plt
        mpl.style.use('classic')
        mpl.rcParams['legend.fontsize'] = 12
        # mpl.rcParams['legend.family'] = 'Times New Roman'
        mpl.rcParams['font.family'] = ['Times New Roman']
        # mpl.rcParams['font.size'] = 15.0
        font = {'family' : 'Times New Roman', #'serif',
                'color' : 'darkblue',
                'weight' : 'normal',
                'size' : 14,}
        textfont = {'family' : 'Times New Roman', #'serif',
                    'color' : 'darkblue',
                    'weight' : 'normal',
                    'size' : 11.5,}

        fig, axeses = subplots(4, 2, sharex=True, dpi=150, figsize=(16*0.75, 8*0.75), facecolor='w', edgecolor='k', constrained_layout=True)
        ax_list = []
        for i in range(4):
            ax_list.extend(axeses[i].tolist())
        # O2_prototype_ax.plot(O2_prototype_data[1], 'o-', lw=0.75, alpha=0.5, label=r'$\delta$'         )
        # O2_prototype_ax.plot(O2_prototype_data[0], 'v-', lw=0.75, alpha=0.5, label=r'$b_{\rm tooth,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[3], 's-', lw=0.75, alpha=0.5, label=r'$b_{\rm tooth,r}$')
        # O2_prototype_ax.plot(O2_prototype_data[5], '^-', lw=0.75, alpha=0.5, label=r'$w_{\rm open,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[2], 'd-', lw=0.75, alpha=0.5, label=r'$w_{\rm open,r}$')
        # O2_prototype_ax.plot(O2_prototype_data[6], '*-', lw=0.75, alpha=0.5, label=r'$h_{\rm head,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[4], 'X-', lw=0.75, alpha=0.5, label=r'$h_{\rm head,r}$')

        # Extract data
        free_param_list = [
        r'$L_g$',        
        r'$w_{st}$',     
        r'$w_{rt}$',     
        r'$\theta_{so}$',
        r'$w_{ro}$',
        r'$d_{so}$',
        r'$d_{ro}$']
        y_label_list = ['PF', r'$\eta$ [100%]', r'$T_{em} [N]$', r'$T_{rip}$ [100%]', r'$|F|$ [N]', r'$E_m$ [100%]', r'$E_a$ [deg]', 
                        r'$P_{Cu,s,JMAG}$', r'$P_{Cu,r,JMAG}$', r'$P_{Fe}$ [W]', r'$P_{eddy}$', r'$P_{hyst}$', r'$P_{Cu,s,FEMM}$', r'$P_{Cu,r,FEMM}$', 
                        r'Windage loss', r'Total loss']

        list_y_label = [r'$O_A$ [$\rm kNm/m^3$]', 
                         '$O_C$ [1]', 
                         'FRW [p.u.]',
                         '$E_a$ [deg]', 
                         '$O_B$ [%]', 
                         '$E_m$ [%]', 
                         r'$P_{\rm loss}$ [W]',
                         r'$T_{\rm rip}$ [%]',
                         # 'Rotor Weight [N]', #'Power Factor [1]',
                         ]
        list_y_data_max = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        list_y_data_min = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

        list_y_data = [ [el/1e3 for el in self.l_OA], 
                        self.l_OC,
                        [F/W for W, F in zip(self.l_original_rotor_weight, self.l_ss_avg_force_magnitude)], # FRW
                        self.l_force_error_angle,
                        [100*el for el in self.l_OB], 
                        [100*el for el in self.l_normalized_force_error_magnitude],
                        self.l_rated_total_loss,
                        [100*el for el in self.l_normalized_torque_ripple],
                        # self.l_rated_rotor_weight,
                        ]
        for i in range(len(list_y_label)):
            ax = ax_list[i]
            y_data = list_y_data[i]
            y_value_reference = y_data[0]
            ax.plot(y_value_reference*np.ones(len(y_data)), '-k', alpha=1, zorder=10)
            y_data = y_data[1:]
            number_of_points_per_geometry_variable = len(y_data)/number_of_free_variables
            for idx, part_of_y_data in enumerate([y_data[int(number_of_points_per_geometry_variable*_)\
                                                        :int(number_of_points_per_geometry_variable*(_+1))]\
                                                        for _ in range(number_of_free_variables)]):
                if idx%2 == 0:
                    line_style = '--bo'
                else:
                    line_style = '--ro'
                ax.plot(list(range(len(y_data)))[int(number_of_points_per_geometry_variable*idx)\
                                                :int(number_of_points_per_geometry_variable*(idx+1))], 
                                                part_of_y_data, line_style, alpha=0.33)

            low, high = ax.get_ylim()
            # ax.legend()
            ax.grid()
            ax.set_ylabel(list_y_label[i], **font)
            ax.set_xlim([0,140])
            for j in range(number_of_free_variables):
                if j%2==0:
                    alpha = 0.05
                else:
                    alpha = 0.15
                ax.axvspan(j*number_of_variant-0.5, (j+1)*number_of_variant-0.5, facecolor='k', alpha=alpha)
                ax.text(0.33*number_of_free_variables+j*number_of_variant, high-(high-low)*0.125, free_param_list[j])

            if i == 0:
                ax.set_yticks([-24, -22, -20, -18, -16])
            list_y_data_max[i].append(max(y_data))
            list_y_data_min[i].append(min(y_data))

        ax_list[-2].set_xlabel('Count of design variant', **font)
        ax_list[-1].set_xlabel('Count of design variant', **font)
        fig.savefig(r'C:\Users\horyc\Desktop/'+ 'LSA_curves.png', dpi=300)
        # plt.show()
        return













        INDEX_TOTAL_LOSS = 15 + 4 # index of total loss in the machine_data list

        # ------------------------------------ Sensitivity Analysis Bar Chart Scripts

        # print next(self.get_list_objective_function())
        data_max = []
        data_min = []
        eta_at_50kW_max = []
        eta_at_50kW_min = []
        O1_max   = []
        O1_min   = []
        for ind, i in enumerate(list(range(7))+[INDEX_TOTAL_LOSS]):
            print('\n-----------', y_label_list[i])
            l = list(self.get_certain_objective_function(i))
            y = l
            print('ind=', ind, 'i=', i, 'len(y)=', len(y))

            data_max.append([])
            data_min.append([])

            for j in range(int(len(y)/number_of_variant)): # iterate design parameters
                y_vs_design_parameter = y[j*number_of_variant:(j+1)*number_of_variant]

                try:
                    # if j == 6:
                    ax_list[ind].plot(y_vs_design_parameter, 'o-', lw=0.75, label=str(j)+' '+param_list[j], alpha=0.5)
                except IndexError as e:
                    print('Check the length of y should be 7*(%d+1)=%d, or else you should remove the redundant results in swarm_data.txt (they are produced because of the interrupted/resumed script run.)'%(number_of_variant, 7*number_of_variant))
                    raise e
                print('\tj=', j, param_list[j], '\t\t Max-Min:', max(y_vs_design_parameter) - min(y_vs_design_parameter))

                data_max[ind].append(max(y_vs_design_parameter))
                data_min[ind].append(min(y_vs_design_parameter))            

            if i==1:
                ax_list[ind].legend(prop={'family':'Times New Roman'})
            ax_list[ind].grid()
            ax_list[ind].set_ylabel(y_label_list[i], **font)

        print('\nObjectives vs. geometry variables:')
        for ind, el in enumerate(data_max):
            print(ind, 'Max', el)
        print('\nObjectives vs. geometry variables:')
        for ind, el in enumerate(data_min):
            print(ind, 'Min', el)

        if self.reference_design is not None:
            print('\n-------------------- Here goes the reference design:')
            for el in self.reference_design[1:]:
                print(el, end=' ')
            self.reference_data = [float(el) for el in self.reference_design[3].split(',')]
            O2_ref = fobj_scalar(self.reference_data[2],
                                 self.reference_data[4],
                                 self.reference_data[3],
                                 self.reference_data[5],
                                 self.reference_data[6],
                                 self.reference_data[INDEX_TOTAL_LOSS],
                                 weights=use_weights('O2'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
            O1_ref = fobj_scalar(self.reference_data[2],
                                 self.reference_data[4],
                                 self.reference_data[3],
                                 self.reference_data[5],
                                 self.reference_data[6],
                                 self.reference_data[INDEX_TOTAL_LOSS],
                                 weights=use_weights('O1'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
        else:
            raise Exception('self.reference_design is None.')

        print('Objective function 1')
        O1 = fobj_list( list(self.get_certain_objective_function(2)), 
                        list(self.get_certain_objective_function(4)), 
                        list(self.get_certain_objective_function(3)), 
                        list(self.get_certain_objective_function(5)), 
                        list(self.get_certain_objective_function(6)), 
                        np.array(list(self.get_certain_objective_function(9))) + np.array(list(self.get_certain_objective_function(12))) + np.array(list(self.get_certain_objective_function(13))),
                        weights=use_weights('O1'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
        O1_max = []
        O1_min = []
        from pylab import figure
        O1_ax  = figure().gca()
        O2_prototype_data = []
        results_for_refining_bounds = {}
        results_for_refining_bounds['O1'] = []
        for j in range(int(len(O1)/number_of_variant)): # iterate design parameters
            O1_vs_design_parameter = O1[j*number_of_variant:(j+1)*number_of_variant]
            O2_prototype_data.append(O1_vs_design_parameter)

            O1_ax.plot(O1_vs_design_parameter, label=str(j)+' '+param_list[j], alpha=0.5)
            print('\t', j, param_list[j], '\t\t max O1 - min O1:', max(O1_vs_design_parameter) - min(O1_vs_design_parameter), '\t\t', end=' ')

            # narrow bounds (refine bounds)
            results_for_refining_bounds['O1'].append( [ind for ind, el in enumerate(O1_vs_design_parameter) if el < O1_ref*1.0] )
            print(results_for_refining_bounds['O1'][j]) #'<- to derive new original_bounds.'

            O1_max.append(max(O1_vs_design_parameter))
            O1_min.append(min(O1_vs_design_parameter))            
        O1_ax.legend()
        O1_ax.grid()
        O1_ax.set_ylabel('O1 [1]', **font)
        O1_ax.set_xlabel('Count of design variants', **font)

        # fig_prototype = figure(500, figsize=(10, 5), facecolor='w', edgecolor='k')
        # O2_prototype_ax = fig_prototype.gca()
        # O2_prototype_ax.plot(list(range(-1, 22)), O1_ref*np.ones(23), 'k--', label='Reference design')
        # O2_prototype_ax.plot(O2_prototype_data[1], 'o-', lw=0.75, alpha=0.5, label=r'$L_g$')
        # O2_prototype_ax.plot(O2_prototype_data[0], 'v-', lw=0.75, alpha=0.5, label=r'$w_{st}$')
        # O2_prototype_ax.plot(O2_prototype_data[3], 's-', lw=0.75, alpha=0.5, label=r'$w_{rt}$')
        # O2_prototype_ax.plot(O2_prototype_data[5], '^-', lw=0.75, alpha=0.5, label=r'$\theta_{so}$')
        # O2_prototype_ax.plot(O2_prototype_data[2], 'd-', lw=0.75, alpha=0.5, label=r'$w_{ro}$')
        # O2_prototype_ax.plot(O2_prototype_data[6], '*-', lw=0.75, alpha=0.5, label=r'$d_{so}$')
        # O2_prototype_ax.plot(O2_prototype_data[4], 'X-', lw=0.75, alpha=0.5, label=r'$d_{ro}$')
        # O2_prototype_ax.legend()
        # O2_prototype_ax.set_ylabel('$O_2(x)$ [1]', **font)

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # O2
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        print('Objective function 2')
        O2 = fobj_list( list(self.get_certain_objective_function(2)), 
                        list(self.get_certain_objective_function(4)), 
                        list(self.get_certain_objective_function(3)), 
                        list(self.get_certain_objective_function(5)), 
                        list(self.get_certain_objective_function(6)), 
                        np.array(list(self.get_certain_objective_function(9))) + np.array(list(self.get_certain_objective_function(12))) + np.array(list(self.get_certain_objective_function(13))),
                        weights=use_weights('O2'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight )
        O2_max = []
        O2_min = []
        O2_ax  = figure().gca()
        O2_ecce_data = []
        results_for_refining_bounds['O2'] = []
        for j in range(int(len(O2)/number_of_variant)): # iterate design parameters: range(7)
            O2_vs_design_parameter = O2[j*number_of_variant:(j+1)*number_of_variant]
            O2_ecce_data.append(O2_vs_design_parameter)

            # narrow bounds (refine bounds)
            O2_ax.plot(O2_vs_design_parameter, 'o-', label=str(j)+' '+param_list[j], alpha=0.5)
            print('\t', j, param_list[j], '\t\t max O2 - min O2:', max(O2_vs_design_parameter) - min(O2_vs_design_parameter), '\t\t', end=' ')
            results_for_refining_bounds['O2'].append( [ind for ind, el in enumerate(O2_vs_design_parameter) if el < O2_ref*1.0] )
            print(results_for_refining_bounds['O2'][j]) #'<- to derive new original_bounds.'

            O2_max.append(max(O2_vs_design_parameter))
            O2_min.append(min(O2_vs_design_parameter))
        O2_ax.legend()
        O2_ax.grid()
        O2_ax.set_ylabel('O2 [1]', **font)
        O2_ax.set_xlabel('Count of design variants', **font)

        # for ecce digest
        fig_ecce = figure(figsize=(10, 5), facecolor='w', edgecolor='k')
        O2_ecce_ax = fig_ecce.gca()
        O2_ecce_ax.plot(list(range(-1, 22)), O2_ref*np.ones(23), 'k--', label='Reference design')
        O2_ecce_ax.plot(O2_ecce_data[1], 'o-', lw=0.75, alpha=0.5,      label=r'$L_g$')
        O2_ecce_ax.plot(O2_ecce_data[0], 'v-', lw=0.75, alpha=0.5,      label=r'$w_{st}$')
        O2_ecce_ax.plot(O2_ecce_data[3], 's-', lw=0.75, alpha=0.5,      label=r'$w_{rt}$')
        O2_ecce_ax.plot(O2_ecce_data[5], '^-', lw=0.75, alpha=0.5,      label=r'$\theta_{so}$')
        O2_ecce_ax.plot(O2_ecce_data[2], 'd-', lw=0.75, alpha=0.5,      label=r'$w_{ro}$')
        O2_ecce_ax.plot(O2_ecce_data[6], '*-', lw=0.75, alpha=0.5,      label=r'$d_{so}$')
        O2_ecce_ax.plot(O2_ecce_data[4], 'X-', lw=0.75, alpha=0.5,      label=r'$d_{ro}$')

        myfontsize = 12.5
        from pylab import plt
        plt.rcParams.update({'font.size': myfontsize})


        # Reference candidate design
        ref = np.zeros(8)
            # ref[0] = 0.635489                                   # PF
            # ref[1] = 0.963698                                   # eta
            # ref[1] = efficiency_at_50kW(1817.22+216.216+224.706)# eta@50kW

        if self.reference_design is not None:
            list_plotting_weights = [8, 3, self.required_torque, 0.1, self.rotor_weight, 0.2, 10, 2500]
            ref[0] = O2_ref                  / list_plotting_weights[0] 
            ref[1] = O1_ref                  / list_plotting_weights[1] 
            ref[2] = self.reference_data[2]  / list_plotting_weights[2]  # 100%
            ref[3] = self.reference_data[3]  / list_plotting_weights[3]  # 100%
            ref[4] = self.reference_data[4]  / list_plotting_weights[4]  # 100% = FRW
            ref[5] = self.reference_data[5]  / list_plotting_weights[5]  # 100%
            ref[6] = self.reference_data[6]  / list_plotting_weights[6]  # deg
            ref[7] = self.reference_data[INDEX_TOTAL_LOSS] / list_plotting_weights[7]  # W

        O1_ax.plot(list(range(-1, 22)), O1_ref*np.ones(23), 'k--')
        O2_ax.plot(list(range(-1, 22)), O2_ref*np.ones(23), 'k--')
        O2_ecce_ax.legend()
        O2_ecce_ax.grid()
        O2_ecce_ax.set_xticks(list(range(21)))
        O2_ecce_ax.annotate('Lower bound', xytext=(0.5, 5.5), xy=(0, 4), xycoords='data', arrowprops=dict(arrowstyle="->"))
        O2_ecce_ax.annotate('Upper bound', xytext=(18.0, 5.5),  xy=(20, 4), xycoords='data', arrowprops=dict(arrowstyle="->"))
        O2_ecce_ax.set_xlim((-0.5,20.5))
        O2_ecce_ax.set_ylim((0,14)) # 4,14
        O2_ecce_ax.set_xlabel(r'Number of design variant', **font)
        O2_ecce_ax.set_ylabel(r'$O_2(x)$ [1]', **font)
        fig_ecce.tight_layout()
        # fig_ecce.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction_full_paper\images\O2_vs_params.png', dpi=150)
        # plt.show()
        # quit() ###################################


        # Maximum
        data_max = np.array(data_max)
        O1_max   = np.array(O1_max)
        O2_max   = np.array(O2_max)
            # data_max[0] = (data_max[0])                   # PF
            # data_max[1] = (data_max[1])                   # eta
            # data_max[1] = efficiency_at_50kW(data_max[7]) # eta@50kW # should use data_min[7] because less loss, higher efficiency
        data_max[0] = O2_max       / list_plotting_weights[0]  
        data_max[1] = O1_max       / list_plotting_weights[1]  
        data_max[2] = (data_max[2])/ list_plotting_weights[2]  # 100%
        data_max[3] = (data_max[3])/ list_plotting_weights[3]  # 100%
        data_max[4] = (data_max[4])/ list_plotting_weights[4]  # 100% = FRW
        data_max[5] = (data_max[5])/ list_plotting_weights[5]  # 100%
        data_max[6] = (data_max[6])/ list_plotting_weights[6]  # deg
        data_max[7] = (data_max[7])/ list_plotting_weights[7]  # W
        y_max_vs_design_parameter_0 = [el[0] for el in data_max]
        y_max_vs_design_parameter_1 = [el[1] for el in data_max]
        y_max_vs_design_parameter_2 = [el[2] for el in data_max]
        y_max_vs_design_parameter_3 = [el[3] for el in data_max]
        y_max_vs_design_parameter_4 = [el[4] for el in data_max]
        y_max_vs_design_parameter_5 = [el[5] for el in data_max]
        y_max_vs_design_parameter_6 = [el[6] for el in data_max]

        # Minimum
        data_min = np.array(data_min)
        O1_min   = np.array(O1_min)
        O2_min   = np.array(O2_min)
            # data_min[0] = (data_min[0])                    # PF
            # data_min[1] = (data_min[1])                    # eta
            # data_min[1] = efficiency_at_50kW(data_min[7])  # eta@50kW
        data_min[0] = O2_min        / list_plotting_weights[0] 
        data_min[1] = O1_min        / list_plotting_weights[1] 
        data_min[2] = (data_min[2]) / list_plotting_weights[2] # 100%
        data_min[3] = (data_min[3]) / list_plotting_weights[3] # 100%
        data_min[4] = (data_min[4]) / list_plotting_weights[4] # 100% = FRW
        data_min[5] = (data_min[5]) / list_plotting_weights[5] # 100%
        data_min[6] = (data_min[6]) / list_plotting_weights[6] # deg
        data_min[7] = (data_min[7]) / list_plotting_weights[7] # W
        y_min_vs_design_parameter_0 = [el[0] for el in data_min]
        y_min_vs_design_parameter_1 = [el[1] for el in data_min]
        y_min_vs_design_parameter_2 = [el[2] for el in data_min]
        y_min_vs_design_parameter_3 = [el[3] for el in data_min]
        y_min_vs_design_parameter_4 = [el[4] for el in data_min]
        y_min_vs_design_parameter_5 = [el[5] for el in data_min]
        y_min_vs_design_parameter_6 = [el[6] for el in data_min]

        count = np.arange(len(y_max_vs_design_parameter_0))  # the x locations for the groups
        width = 1.0  # the width of the bar

        fig = figure(dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = fig.gca()
        # fig, ax = plt.subplots(dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')                                      #  #1034A
        rects1 = ax.bar(count - 3*width/8, y_min_vs_design_parameter_0, width/8, alpha=0.5, label=r'$L_g$, Air gap length', color='#6593F5')
        rects2 = ax.bar(count - 2*width/8, y_min_vs_design_parameter_1, width/8, alpha=0.5, label=r'$b_{st}$, Stator tooth width', color='#1D2951') # https://digitalsynopsis.com/design/beautiful-color-palettes-combinations-schemes/
        rects3 = ax.bar(count - 1*width/8, y_min_vs_design_parameter_2, width/8, alpha=0.5, label=r'$b_{rt}$, Rotor tooth width', color='#03396c')
        rects4 = ax.bar(count - 0*width/8, y_min_vs_design_parameter_3, width/8, alpha=0.5, label=r'$\theta_{so}$, Stator open width', color='#6497b1')
        rects5 = ax.bar(count + 1*width/8, y_min_vs_design_parameter_4, width/8, alpha=0.5, label=r'$w_{ro}$, Rotor open width',  color='#0E4D92')
        rects6 = ax.bar(count + 2*width/8, y_min_vs_design_parameter_5, width/8, alpha=0.5, label=r'$d_{so}$, Stator open depth', color='#005b96')
        rects7 = ax.bar(count + 3*width/8, y_min_vs_design_parameter_6, width/8, alpha=0.5, label=r'$d_{ro}$, Rotor open depth', color='#b3cde0') 
        print('ylim=', ax.get_ylim())
        autolabel(ax, rects1, bias=-0.10, textfont=textfont)
        autolabel(ax, rects2, bias=-0.10, textfont=textfont)
        autolabel(ax, rects3, bias=-0.10, textfont=textfont)
        autolabel(ax, rects4, bias=-0.10, textfont=textfont)
        autolabel(ax, rects5, bias=-0.10, textfont=textfont)
        autolabel(ax, rects6, bias=-0.10, textfont=textfont)
        autolabel(ax, rects7, bias=-0.10, textfont=textfont)
        one_one = np.array([1, 1])
        minus_one_one = np.array([-1, 1])
        ax.plot(rects4[0].get_x() + 0.5*width*minus_one_one, ref[0]*one_one, 'k--', lw=1.0, alpha=0.6, label='Reference design' )
        ax.plot(rects4[1].get_x() + 0.5*width*minus_one_one, ref[1]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[2].get_x() + 0.5*width*minus_one_one, ref[2]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[3].get_x() + 0.5*width*minus_one_one, ref[3]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[4].get_x() + 0.5*width*minus_one_one, ref[4]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[5].get_x() + 0.5*width*minus_one_one, ref[5]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[6].get_x() + 0.5*width*minus_one_one, ref[6]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[7].get_x() + 0.5*width*minus_one_one, ref[7]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.legend(loc='upper right', prop={'family':'Times New Roman'})
        # text for indicating reference values
        ax.text(rects4[0].get_x() - 3.5/8*width, ref[0]*1.01, '%.2f'%(ref[0]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[1].get_x() - 3.5/8*width, ref[1]*1.01, '%.2f'%(ref[1]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[2].get_x() - 3.5/8*width, ref[2]*1.01, '%.2f'%(ref[2]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[3].get_x() - 3.5/8*width, ref[3]*1.01, '%.2f'%(ref[3]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[4].get_x() - 3.5/8*width, ref[4]*1.01, '%.2f'%(ref[4]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[5].get_x() - 3.5/8*width, ref[5]*1.01, '%.2f'%(ref[5]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[6].get_x() - 3.5/8*width, ref[6]*1.01, '%.2f'%(ref[6]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[7].get_x() - 3.5/8*width, ref[7]*1.01, '%.2f'%(ref[7]), ha='center', va='bottom', rotation=90, **textfont)

        rects1 = ax.bar(count - 3*width/8, y_max_vs_design_parameter_0, width/8, alpha=0.5, label=r'$L_g$,         Air gap length', color='#6593F5')    # bottom=y_min_vs_design_parameter_0, 
        rects2 = ax.bar(count - 2*width/8, y_max_vs_design_parameter_1, width/8, alpha=0.5, label=r'$b_{st}$, Stator tooth width', color='#1D2951')     # bottom=y_min_vs_design_parameter_1, 
        rects3 = ax.bar(count - 1*width/8, y_max_vs_design_parameter_2, width/8, alpha=0.5, label=r'$b_{rt}$, Rotor tooth width', color='#03396c')      # bottom=y_min_vs_design_parameter_2, 
        rects4 = ax.bar(count - 0*width/8, y_max_vs_design_parameter_3, width/8, alpha=0.5, label=r'$\theta_{so}$, Stator open width', color='#6497b1') # bottom=y_min_vs_design_parameter_3, 
        rects5 = ax.bar(count + 1*width/8, y_max_vs_design_parameter_4, width/8, alpha=0.5, label=r'$w_{ro}$, Rotor open width',  color='#0E4D92')      # bottom=y_min_vs_design_parameter_4, 
        rects6 = ax.bar(count + 2*width/8, y_max_vs_design_parameter_5, width/8, alpha=0.5, label=r'$d_{so}$, Stator open depth', color='#005b96')      # bottom=y_min_vs_design_parameter_5, 
        rects7 = ax.bar(count + 3*width/8, y_max_vs_design_parameter_6, width/8, alpha=0.5, label=r'$d_{ro}$, Rotor open depth', color='#b3cde0')       # bottom=y_min_vs_design_parameter_6, 
        autolabel(ax, rects1, textfont=textfont)
        autolabel(ax, rects2, textfont=textfont)
        autolabel(ax, rects3, textfont=textfont)
        autolabel(ax, rects4, textfont=textfont)
        autolabel(ax, rects5, textfont=textfont)
        autolabel(ax, rects6, textfont=textfont)
        autolabel(ax, rects7, textfont=textfont)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Normalized Objective Functions', **font)
        ax.set_xticks(count)
        # ax.set_xticklabels(('Power Factor [100%]', r'$\eta$@$T_{em}$ [100%]', r'$T_{em}$ [15.9 N]', r'$T_{rip}$ [10%]', r'$|F|$ [51.2 N]', r'    $E_m$ [20%]', r'      $E_a$ [10 deg]', r'$P_{\rm Cu,Fe}$ [2.5 kW]')))
        # ax.set_xticklabels(('Power Factor [100%]', r'$O_1$ [3]', r'$T_{em}$ [15.9 N]', r'$T_{rip}$ [10%]', r'$|F|$ [51.2 N]', r'    $E_m$ [20%]', r'      $E_a$ [10 deg]', r'$P_{\rm Cu,Fe}$ [2.5 kW]'))
        ax.set_xticklabels(('$O_2$ [%g]'               %(list_plotting_weights[0]), 
                            '$O_1$ [%g]'               %(list_plotting_weights[1]), 
                            '$T_{em}$ [%g Nm]'         %(list_plotting_weights[2]), 
                            '$T_{rip}$ [%g%%]'         %(list_plotting_weights[3]*100), 
                            '$|F|$ [%g N]'             %(list_plotting_weights[4]), 
                            '    $E_m$ [%g%%]'         %(list_plotting_weights[5]*100), 
                            '      $E_a$ [%g deg]'     %(list_plotting_weights[6]), 
                            '$P_{\\rm Cu,Fe}$ [%g kW]' %(list_plotting_weights[7]*1e-3) ), **font)
        ax.grid()
        ax.set_ylim([0,4])
        # fig.tight_layout()
        # fig.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction\images\sensitivity_results.png', dpi=150)

        # plt.show()
        return results_for_refining_bounds


