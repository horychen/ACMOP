from dataclasses import dataclass, fields
from typing import Dict, List, Optional, Any
from collections import OrderedDict
import json, math, base64, pickle

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
            calc_bounds=None  # calc_bounds 函数需要从其他地方重建
        )

class Winding(object):
    def __init__(self, phase_number_m: int, stator_slot_number_Qs: int, pole_pair_number_p: int, suspension_pole_pair_number_ps: int) -> None:
        self.phase_number_m = phase_number_m
        self.stator_slot_number_Qs = stator_slot_number_Qs
        self.pole_pair_number_p = pole_pair_number_p
        self.suspension_pole_pair_number_ps = suspension_pole_pair_number_ps
    
    def to_dict(self) -> Dict[str, Any]:
        """将 Winding 对象转换为字典"""
        return {
            'phase_number_m': self.phase_number_m,
            'stator_slot_number_Qs': self.stator_slot_number_Qs,
            'pole_pair_number_p': self.pole_pair_number_p,
            'suspension_pole_pair_number_ps': self.suspension_pole_pair_number_ps
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Winding':
        """从字典创建 Winding 对象"""
        return cls(
            phase_number_m=data['phase_number_m'],
            stator_slot_number_Qs=data['stator_slot_number_Qs'],
            pole_pair_number_p=data['pole_pair_number_p'],
            suspension_pole_pair_number_ps=data['suspension_pole_pair_number_ps']
        )

class Geometry(object):
    def __init__(self, GP: dict = None, draw_function: callable = None, _calculate_points: callable = None, color: str = None):
        self.color = color
        self.GP = GP or []
        self.draw_function = draw_function
        self._calculate_points = _calculate_points
        for i, (name, gp) in enumerate(self.GP.items()):
            if isinstance(gp, Parameter):
                exec(f"self.{name} = {gp.value}")
    def to_dict(self) -> Dict[str, Any]:
        # INSERT_YOUR_CODE
        """
        扩展序列化方法：将所有当前成员变量都存入字典（不局限于定义时的变量）
        排除不能序列化的 draw_function/_calculate_points 属性
        """
        result = {}
        # 收集所有实例属性
        for k, v in self.__dict__.items():
            if k in ['draw_function', '_calculate_points']:
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
            _calculate_points=data['_calculate_points']
        )
    def draw(self, drawer, *args, **kwargs):
        self.components_make_region = self.draw_function(drawer, *args, **kwargs)

@dataclass
class Modern_Machine_Designer(object):

    machine_class: str = 'bearingless_spmsm_heart.bearingless_spmsm_design_variant'
    bool_PermanentMagnet: bool = True
    bool_StatorSlotClosed: bool = False
    bool_RotorNotched: bool = True

    select_FEA_tool: str = 'JMAG'
    name: str = 'SuperCoolName'

    def __post_init__(self):

        Qs = 12
        p = 4
        ps = 5
        coil_pitch_y = 1

        # 裂比
        SR = 0.35

        # 定子
        mm_r_so = 123.5

        # 利用不同的裂比去估算合理的边界值
        yoke_split_ratio_bounds = [0.2, 0.45]
        tooth_split_ratio_at_middle_slot = [0.25, 0.50]

        '''Fixed variables'''
        self.m: Parameter            = Parameter('phase_number_m', 'fixed', 3)
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

        # TODO: 下界需要考虑到w_st的宽度和半径的值
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

        import CrossSectInnerNotchedRotor, CrossSectStator
        self.machineGeometry = {
            "rotorCore": Geometry(
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
                ),
                _calculate_points=lambda: (
                    CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor.calculate_points 
                )
            ),
            "shaft": Geometry(
                GP={'mm_r_ri': self.mm_r_ri},
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectShaft(
                        name="shaft",
                        color="#0EE0E2",
                        notched_rotor=CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
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
                _calculate_points=lambda: None,
            ),
            "rotorMagnet": Geometry(
                GP={
                    'mm_d_pm': self.mm_d_pm,
                    'mm_d_ri': self.mm_d_ri,
                    'mm_r_ri': self.mm_r_ri,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet(
                        name="rotorMagnet",
                        color="#1C96E0",
                        notched_rotor=CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
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
                _calculate_points=lambda: None,
            ),
            "statorCore": Geometry(
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
                _calculate_points=lambda: None,
            ),
            "coils": None
        }
        self.machineGeometry['coils'] = Geometry(
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
                ),
                _calculate_points=lambda: None
            )

        self.wily = Winding(phase_number_m=self.m.value if hasattr(self, "m") and self.m.value is not None else 3, stator_slot_number_Qs=self.Qs.value if hasattr(self, "Qs") and self.Qs.value is not None else 12, pole_pair_number_p=self.p.value if hasattr(self, "p") and self.p.value is not None else 4, suspension_pole_pair_number_ps=self.ps.value if hasattr(self, "ps") and self.ps.value is not None else 5)

    def getSketch(self, name, color):
        self.name = name
        self.color = color

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
            # logger = logging.getLogger(__name__)
            # logger.debug('cos夹角=%s', cos夹角)
        elif -1.0-EPS < cos夹角 < -1.0:
            cos夹角 = -1.0
            # logger = logging.getLogger(__name__)
            # logger.debug('cos夹角=%s', cos夹角)
        angle_between = math.acos(cos夹角)

        radius = math.sqrt(v1[0]*v1[0] + v1[1]*v1[1])
        angle_start = math.atan2(v1[1], v1[0])
        angle_end = angle_start + angle_between

        self.ctx.move_to(startxy[0], startxy[1])
        self.ctx.arc(centerxy[0], centerxy[1], radius, angle_start, angle_end)
        # self.ctx.arc_negative(centerxy[0], centerxy[1], radius, angle_end, angle_start)
        return [{'move_to': (centerxy[0], centerxy[1]), 'arc': (radius, angle_start, angle_end)}]

    def show_geometry(self, filename=None) -> None:
        # 检查必要的参数是否存在
        if not hasattr(self, 'mm_r_so') or self.mm_r_so.value is None:
            raise ValueError("mm_r_so parameter is required but not found or has no value")
        
        width_in_points  = self.mm_r_so.value*2.1
        height_in_points = self.mm_r_so.value*2.1
        
        # mm_r_ro 只在 bool_PermanentMagnet 为 True 时存在
        if hasattr(self, 'mm_r_ro') and self.mm_r_ro.value is not None:
            lw = 0.1 if self.mm_r_ro.value < 15 else 0.5
        else:
            # 使用默认值
            lw = 0.5
        
        bool_draw_whole_model = True
        
        # 确保 machineGeometry 已初始化
        if not hasattr(self, 'machineGeometry') or self.machineGeometry is None:
            # 如果 machineGeometry 不存在，调用 __post_init__ 来创建
            self.__post_init__()

        def draw_spmsm(lw, width_in_points, height_in_points, bool_draw_whole_model):
            import cairo

            def init_canvas(width_in_points, height_in_points):
                self.surface = cairo.SVGSurface('machine_geometry.svg', width_in_points, height_in_points)
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
            def apply_stroke(lw=0.5):
                self.ctx.set_line_cap(cairo.LINE_CAP_ROUND)
                self.ctx.set_line_width(lw)
                # setting color of the context
                self.ctx.set_source_rgba(0.0, 0.0, 0.0, 1)
                # stroke out the color and width property
                self.ctx.stroke()
            def convert_to_pdf(bool_open_pdf=False, filename=None): # 这个代码只是把SVG转换为PDF而已
                # self.surface.write_to_svg()
                self.surface.finish()
                import cairosvg
                cairosvg.svg2pdf(url=f'machine_geometry.svg', write_to=f'machine_geometry.pdf')
                if bool_open_pdf:
                    import os
                    os.system('sumatraPDF2.exe ' + 'machine_geometry.pdf')
                print('[machine_design_guide.py] Find the file machine_geometry.pdf in the current folder.')

            init_canvas(width_in_points, height_in_points)

            # 检查 machineGeometry 是否存在且包含必要的键
            if not hasattr(self, 'machineGeometry') or self.machineGeometry is None:
                raise ValueError("machineGeometry is not initialized. Please ensure __post_init__ was called.")
            
            # 安全地调用 draw 方法
            if 'rotorCore' in self.machineGeometry and self.machineGeometry['rotorCore'] is not None:
                list_regions = self.machineGeometry['rotorCore'].draw(self, bool_draw_whole_model=bool_draw_whole_model)
            if 'shaft' in self.machineGeometry and self.machineGeometry['shaft'] is not None:
                list_regions = self.machineGeometry['shaft'].draw(self)
            if 'rotorMagnet' in self.machineGeometry and self.machineGeometry['rotorMagnet'] is not None:
                list_regions = self.machineGeometry['rotorMagnet'].draw(self, bool_draw_whole_model=bool_draw_whole_model)
            if 'statorCore' in self.machineGeometry and self.machineGeometry['statorCore'] is not None:
                list_regions = self.machineGeometry['statorCore'].draw(self, bool_draw_whole_model=bool_draw_whole_model)
            if 'coils' in self.machineGeometry and self.machineGeometry['coils'] is not None:
                list_regions = self.machineGeometry['coils'].draw(self, bool_draw_whole_model=bool_draw_whole_model)

            apply_stroke(lw=lw)
            convert_to_pdf()

            # import builtins
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['rotorCore'] = acm_variant.rotorCore
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['shaft'] = acm_variant.shaft
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['rotorMagnet'] = acm_variant.rotorMagnet
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['sleeve'] = acm_variant.sleeve
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['statorCore'] = acm_variant.statorCore
            # builtins.ad.visualize_dict['GeometricComponentsObjects']['coils'] = acm_variant.coils

        draw_spmsm(lw, width_in_points, height_in_points, bool_draw_whole_model=bool_draw_whole_model)

    def FEA_evaluate(self, project_loc=fr'../_default/'):

        import os
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
        if not os.path.isdir(path2SwarmData): os.makedirs(path2SwarmData)

        self.project_name = self.name+'proj'
        self.expected_project_file = self.path2SwarmData + "temp/%s.jproj"%(self.project_name)

        self.path2FEACsv    = self.path2SwarmData + 'csv/'
        if not os.path.isdir(self.dir_csv_output_folder): os.makedirs(self.dir_csv_output_folder)

        if 'JMAG' in self.select_FEA_tool:

            # study name
            study_name = self.project_name + "-Transient" # Change here and there 

            # Leave the solving task to JMAG
            self.toolJd = self.build_jmag_project(acm_variant, self.project_meta_data, bool_re_evaluate=bool_re_evaluate)
            import rich
            rich.print(self.project_meta_data)

            ################################################################
            # Load data for cost function evaluation
            ################################################################
            acm_variant.results_to_be_unpacked = results_to_be_unpacked = self.toolJd.build_str_results(acm_variant, self.project_name, study_name, self.dir_csv_output_folder, self.fea_config_dict, femm_solver=None)
            if results_to_be_unpacked is not None:
                if self.toolJd.fig_main is not None:
                    try:
                        if False:
                            self.toolJd.fig_main.savefig(self.fea_config_dict['output_dir'] + acm_variant.name + 'results.png', dpi=150)
                    except Exception as e:
                        logger = logging.getLogger(__name__)
                        logger.error('Exception in saving figure: %s', e)
                        logger.info('Ignore error and continue.')
                    finally:
                        utility.pyplot_clear(self.toolJd.axeses)
                # show()
                return acm_variant 
            else:
                raise Exception('[acm_designer] results_to_be_unpacked is None.')

        elif 'FEMM' in self.select_FEA_tool:
            self.toolFEMM = self.build_femm_project(acm_variant)
            # acm_variant.results_to_be_unpacked = results_to_be_unpacked = toolFEMM.build_str_results(self.axeses, acm_variant, self.project_name, study_name, self.dir_csv_output_folder, self.fea_config_dict, femm_solver=None)
            return acm_variant
        else:
            raise Exception('[acm_designer.py] Wrong string of select_FEA_tool:', self.select_FEA_tool)


        if 'JMAG' in self.select_FEA_tool:

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
            ss_avg_force_magnitude, rotor_weight, torque_average = acm_variant.results_to_be_unpacked

            # acm_variant.spec_geometry_dict['x_denorm'] = list(x_denorm)

            spec_performance_dict = dict()
            spec_performance_dict['x_denorm_dict'] = self.acm_template.SI['x_denorm_dict']
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

            GP = acm_variant.template.SI['GP']
            EX = acm_variant.template.SI['EX']

            # Save to disk
            # self.save_to_disk(acm_variant, spec_performance_dict, GP, EX)

            number_current_generation = spec_performance_dict['number_current_generation'] #= int(acm_variant.counter//popsize), 
            individual_index = spec_performance_dict['individual_index'] #= acm_variant.counter
            builtins.ad.visualize_dict[f'FEA_Evaluated_Performance-{number_current_generation}-{individual_index}'] = spec_performance_dict
            json_file_path = self.fea_config_dict['output_dir'] + self.select_spec + '.json'

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
            key = f'gen{number_current_generation}-ind{individual_index}'
            loaded_json[key] = builtins.ad.visualize_dict

            json_string = jsonpickle.encode(loaded_json, indent=4)
            with open(json_file_path, 'w+') as f:
                f.write(json_string)

            number_current_generation = spec_performance_dict['number_current_generation'] #= int(acm_variant.counter//popsize), 
            individual_index = spec_performance_dict['individual_index'] #= acm_variant.counter

            # this is for optimization
            acm_variant.results_for_optimization = (cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle)

    def get_free_variables(self) -> List[Parameter]:
        return [param for param in self.get_parameter_fields().values() if param.type == 'free']

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
        for field in fields(self):
            field_value = getattr(self, field.name)
            if isinstance(field_value, Parameter):
                param_fields[field.name] = field_value
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
        将 Modern_Machine_Designer 对象转换为完整字典（包含所有信息，类似 pickle）
        使用 dill/pickle + base64 编码来保存无法直接序列化的部分
        
        Returns:
            Dict: 包含完整对象信息的字典
        """
        # 首先获取标准的字典表示
        standard_dict = self.to_dict()
        
        # 由于 pickle 无法序列化 lambda 函数，我们需要创建一个可序列化的版本
        # 方法：创建一个新实例，复制所有属性值，但不复制 lambda 函数
        # 然后在反序列化后通过 __post_init__ 重建 lambda 函数
        
        # 保存 lambda 函数的状态信息（用于标记需要重建）
        lambda_info = {}
        
        # 创建新实例并复制所有基本属性
        obj_copy = self.__class__.__new__(self.__class__)
        obj_copy.machine_class = self.machine_class
        obj_copy.bool_PermanentMagnet = self.bool_PermanentMagnet
        obj_copy.bool_StatorSlotClosed = self.bool_StatorSlotClosed
        obj_copy.bool_RotorNotched = self.bool_RotorNotched
        
        # 复制所有参数，但移除 lambda 函数
        for field_name, param in self.get_parameter_fields().items():
            # 创建新的 Parameter 对象，复制所有值，但移除 lambda 函数
            new_param = Parameter(
                name=param.name,
                type=param.type,
                value=param.value,
                bounds=param.bounds,
                unit=param.unit,
                comment=param.comment,
                calc=None,  # lambda 函数会被移除
                calc_bounds=None,  # lambda 函数会被移除
                args=param.args
            )
            setattr(obj_copy, field_name, new_param)
            
            # 记录哪些参数有 lambda 函数
            if param.calc is not None and callable(param.calc):
                func_name = getattr(param.calc, '__name__', '')
                if func_name == '<lambda>' or func_name == '':
                    lambda_info[field_name] = {'has_calc_lambda': True}
            
            if param.calc_bounds is not None and callable(param.calc_bounds):
                func_name = getattr(param.calc_bounds, '__name__', '')
                if func_name == '<lambda>' or func_name == '':
                    if field_name not in lambda_info:
                        lambda_info[field_name] = {}
                    lambda_info[field_name]['has_calc_bounds_lambda'] = True
        
        # 复制 Winding 对象（如果存在）
        if hasattr(self, 'wily') and self.wily is not None:
            obj_copy.wily = Winding.from_dict(self.wily.to_dict())
        
        # machineGeometry 会在 __post_init__ 中重建，所以不需要复制
        
        # 现在尝试序列化
        try:
            pickled_data = pickle.dumps(obj_copy)
            pickled_base64 = base64.b64encode(pickled_data).decode('utf-8')
            pickle_success = True
        except Exception as e:
            # 如果还是失败，只保存标准数据
            pickled_base64 = None
            pickle_success = False
            print(f"Warning: Failed to pickle object: {e}")
        
        # 创建完整字典
        full_dict = {
            '_standard_data': standard_dict,  # 保留标准数据以便前端读取
            '_lambda_info': lambda_info,  # 保存 lambda 函数信息，用于重建
            '_metadata': {
                'class_name': self.__class__.__name__,
                'module': self.__class__.__module__,
                'has_machineGeometry': hasattr(self, 'machineGeometry') and self.machineGeometry is not None,
                'has_wily': hasattr(self, 'wily') and self.wily is not None,
                'pickle_success': pickle_success,
            }
        }
        
        if pickle_success:
            full_dict['_pickle_data'] = pickled_base64
            full_dict['_pickle_version'] = pickle.HIGHEST_PROTOCOL if hasattr(pickle, 'HIGHEST_PROTOCOL') else 4
        else:
            full_dict['_pickle_data'] = None
            full_dict['_error'] = 'Failed to pickle object. Use standard serialization instead.'
        
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
        从完整字典创建 Modern_Machine_Designer 对象（从 pickle 数据恢复）
        
        Args:
            data: 包含完整对象数据的字典（包含 _pickle_data）
            
        Returns:
            Modern_Machine_Designer: 重建的对象
        """
        if '_pickle_data' in data and data['_pickle_data'] is not None:
            try:
                # 从 pickle 数据恢复
                pickled_base64 = data['_pickle_data']
                pickled_data = base64.b64decode(pickled_base64.encode('utf-8'))
                instance = pickle.loads(pickled_data)
                
                # 由于 lambda 函数在序列化时被移除了，需要重新调用 __post_init__ 来重建它们
                # __post_init__ 会重新创建所有 lambda 函数和 machineGeometry
                instance.__post_init__()
                
                return instance
            except Exception as e:
                # 如果 pickle 恢复失败，使用标准方法
                print(f"Warning: Failed to unpickle object, using standard deserialization: {e}")
                return cls.from_dict(data.get('_standard_data', data))
        else:
            # 如果没有 pickle 数据，使用标准方法恢复
            return cls.from_dict(data.get('_standard_data', data))
    
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
    mmd.save_to_file('machine_designer.json')

    mmd.show_geometry()
    print(dir(mmd.machineGeometry['statorCore']))
    quit()

    # 保存完整信息到文件（类似 pickle）
    mmd.save_to_file_full('machine_designer_full.json')
    print("=== 已保存完整信息到 machine_designer_full.json ===")

    # 从 JSON 文件恢复对象
    mmd2 = Modern_Machine_Designer.load_from_file('machine_designer.json')

    # 从完整文件恢复对象（包含所有信息，包括 lambda 函数）
    mmd3 = Modern_Machine_Designer.load_from_file_full('machine_designer_full.json')
    print("=== 已从完整文件恢复对象 ===")

quit()

if __name__ == "__main__":
    # 创建实例
    mmd = Modern_Machine_Designer()
    
    # 导出为 JSON
    json_str = mmd.to_json()
    print("=== JSON 输出 ===")
    print(json_str)
    
    # 保存到文件
    mmd.save_to_file('machine_designer.json')
    print("\n=== 已保存到 machine_designer.json ===")
    
    # 从 JSON 字符串重建
    mmd2 = Modern_Machine_Designer.from_json(json_str)
    print(f"\n=== 从 JSON 重建的对象 ===")
    print(mmd2)
    
    # 从文件加载
    mmd3 = Modern_Machine_Designer.load_from_file('machine_designer.json')
    print(f"\n=== 从文件加载的对象 ===")
    print(mmd3)
    
    # 验证参数摘要
    summary = mmd3.get_parameters_summary()
    print(f"\n=== 参数摘要 ===")
    print(f"总参数数: {summary['total_count']}")
    print(f"按类型: {summary['by_type']}")
    print(f"按单位: {summary['by_unit']}")
    print(f"有值: {summary['with_values']}, 无值: {summary['without_values']}")

    mmd.show_geometry()

