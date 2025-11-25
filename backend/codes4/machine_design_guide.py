from dataclasses import dataclass, fields
from typing import Dict, List, Optional, Any
from collections import OrderedDict
import json, math

class Parameter(object):
    def __init__(self, name, type, value=None, bounds=None, calc=None, calc_bounds=None, unit='mm', comment=None, args=None) -> None:
        # todo: add validation for the type, value, bounds, calc, unit, comment
        self.name = name
        self.type = type
        self.value = value
        self.bounds = bounds
        self.unit = unit
        self.calc = calc
        self.calc_bounds = calc_bounds
        self.comment = comment
        self.args = args
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

        if type == 'free' and value is None:
            self.value = self.bounds[0]+ (self.bounds[1]-self.bounds[0])*0.5

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
    def __init__(self, color: str, GP: list[Parameter], draw_function: callable, _calculate_points: callable):
        self.color = color
        self.GP = GP
        self.draw_function = draw_function
        self._calculate_points = _calculate_points
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'color': self.color,
            'GP': self.GP,
            'draw_function': self.draw_function,
            '_calculate_points': self._calculate_points
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Geometry':
        return cls(
            color=data['color'],
            GP=data['GP'],
            draw_function=data['draw_function'],
            _calculate_points=data['_calculate_points']
        )
    def draw(self):
        return self.draw_function(self.GP)
    def _calculate_points(self):
        return self._calculate_points(self.GP)

@dataclass
class Modern_Machine_Designer(object):

    machine_class: str = 'bearingless_spmsm_heart.bearingless_spmsm_design_variant'
    bool_PermanentMagnet: bool = True
    bool_StatorSlotClosed: bool = False
    bool_RotorNotched: bool = True

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

        # Fixed variables
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

        # Free variables
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

        # derived variables have dependency on the other geometric parameters
        self.mm_r_si: Parameter      = Parameter('stator_inner_radius', 'derived', calc=lambda mm_r_so, split_ratio: mm_r_so * split_ratio, args=[self.mm_r_so.value, self.split_ratio.value])
        if self.bool_PermanentMagnet:
            self.mm_d_ri: Parameter      = Parameter('rotor_iron (back iron) depth', 'derived', calc=lambda mm_d_pm: 4 if mm_d_pm < 4 else mm_d_pm, args=[self.mm_d_pm.value])
            self.mm_r_ro: Parameter      = Parameter('rotor_outer_radius', 'derived', calc=lambda mm_r_si, mm_d_mech_air_gap, mm_d_sleeve: mm_r_si - mm_d_mech_air_gap - mm_d_sleeve, args=[self.mm_r_si.value, self.mm_d_mech_air_gap.value, self.mm_d_sleeve.value])
            self.mm_r_ri: Parameter      = Parameter('rotor_inner_radius', 'derived', calc=lambda r_ro, mm_d_pm, mm_d_ri: r_ro-mm_d_pm-mm_d_ri, args=[self.mm_r_ro.value, self.mm_d_pm.value, self.mm_d_ri.value])
        self.mm_d_st: Parameter      = Parameter('stator_tooth_depth', 'derived', calc=lambda mm_r_so, mm_r_si, mm_d_sy, mm_d_sts: mm_r_so - mm_r_si - mm_d_sy - mm_d_sts, args=[self.mm_r_so.value, self.mm_r_si.value, self.mm_d_sy.value, self.mm_d_sts.value])

        if not self.bool_StatorSlotClosed:
            # deg_alpha_sto 依赖于 deg_alpha_st，使用 deg_alpha_st 的当前值（如果已计算）或使用 bounds 的中间值
            deg_alpha_st_value = self.deg_alpha_st.value if self.deg_alpha_st.value is not None else (self.deg_alpha_st.bounds[0] + self.deg_alpha_st.bounds[1]) / 2 if self.deg_alpha_st.bounds else 360/12*0.1*0.5
            self.deg_alpha_sto: Parameter = Parameter('stator_tooth_open_angle', 'derived', calc=lambda deg_alpha_st: deg_alpha_st*0.5, args=[deg_alpha_st_value])

        if self.bool_RotorNotched:
            self.deg_alpha_rm: Parameter = Parameter('magnet_pole_span_angle', 'free', bounds=[360/self.p.value*0.7, 360/self.p.value], args=[self.p.value])

            # deg_alpha_rs 依赖于 deg_alpha_rm，使用 deg_alpha_rm 的当前值（如果已计算）或使用 bounds 的中间值
            deg_alpha_rm_value = self.deg_alpha_rm.value if self.deg_alpha_rm.value is not None else (self.deg_alpha_rm.bounds[0] + self.deg_alpha_rm.bounds[1]) / 2 if self.deg_alpha_rm.bounds else 360/12*0.1
            self.deg_alpha_rs: Parameter = Parameter('magnet_segment_span_angle', 'derived', calc=lambda deg_alpha_rm: deg_alpha_rm, args=[deg_alpha_rm_value])

            self.mm_d_rp: Parameter      = Parameter('inter_polar_iron_thickness', 'derived', calc=lambda mm_d_pm: mm_d_pm, args=[self.mm_d_pm.value])
            self.mm_d_rs: Parameter      = Parameter('inter_segment_iron_thickness', 'derived', calc=lambda mm_d_pm: mm_d_pm, args=[self.mm_d_pm.value])

        import CrossSectInnerNotchedRotor, CrossSectStator
        self.machineGeometry = {
            "rotorCore": Geometry(
                GP=[
                    self.mm_r_ro, 
                    self.mm_d_ri if hasattr(self, 'mm_d_ri') else None, 
                    self.mm_d_pm if hasattr(self, 'mm_d_pm') else None, 
                    self.mm_d_rp if hasattr(self, 'mm_d_rp') else None, 
                    self.mm_d_rs if hasattr(self, 'mm_d_rs') else None, 
                    self.p if hasattr(self, 'p') else None, 
                    self.s if hasattr(self, 's') else None
                ],
                draw_function=lambda drawer, **kwargs: CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
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
                ).draw(drawer, **kwargs),
                _calculate_points=lambda: CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor.calculate_points
            ),
            "shaft": Geometry(
                GP=[self.mm_r_ri if hasattr(self, 'mm_r_ri') else None],
                draw_function=lambda drawer, **kwargs: CrossSectInnerNotchedRotor.CrossSectShaft(
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
                ).draw(drawer, **kwargs),
                _calculate_points=lambda: None,
            ),
            "rotorMagnet": Geometry(
                GP=[
                    self.mm_d_pm if hasattr(self, "mm_d_pm") else None,
                    self.mm_d_ri if hasattr(self, "mm_d_ri") else None,
                    self.mm_r_ri if hasattr(self, "mm_r_ri") else None,
                ],
                draw_function=lambda drawer, **kwargs: CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet(
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
                ).draw(drawer, **kwargs),
                _calculate_points=lambda: None,
            ),
            "statorCore": Geometry(
                GP=[
                    self.mm_r_si if hasattr(self, "mm_r_si") else None,
                    self.mm_d_sto if hasattr(self, "mm_d_sto") else None,
                    self.mm_d_stt if hasattr(self, "mm_d_stt") else None,
                    self.mm_d_st if hasattr(self, "mm_d_st") else None,
                    self.mm_d_sy if hasattr(self, "mm_d_sy") else None,
                    self.mm_w_st if hasattr(self, "mm_w_st") else None,
                    self.mm_r_st if hasattr(self, "mm_r_st") else None,
                    self.mm_r_sf if hasattr(self, "mm_r_sf") else None,
                    self.mm_r_sb if hasattr(self, "mm_r_sb") else None,
                    self.Q if hasattr(self, "Q") else None,
                ],
                draw_function=lambda drawer, **kwargs: CrossSectStator.CrossSectInnerRotorStator(
                    name="statorCore",
                    color="#BAFA01",
                    deg_alpha_st=self.deg_alpha_st.value if hasattr(self, "deg_alpha_st") and self.deg_alpha_st.value is not None else 40,
                    deg_alpha_sto=self.deg_alpha_sto.value if hasattr(self, "deg_alpha_sto") and self.deg_alpha_sto.value is not None else 20,
                    mm_r_si=self.mm_r_si.value if hasattr(self, "mm_r_si") and self.mm_r_si.value is not None else 40,
                    mm_d_sto=self.mm_d_sto.value if hasattr(self, "mm_d_sto") and self.mm_d_sto.value is not None else 5,
                    mm_d_stt=self.mm_d_stt.value if hasattr(self, "mm_d_stt") and self.mm_d_stt.value is not None else 10,
                    mm_d_st=self.mm_d_st.value if hasattr(self, "mm_d_st") and self.mm_d_st.value is not None else 15,
                    mm_d_sy=self.mm_d_sy.value if hasattr(self, "mm_d_sy") and self.mm_d_sy.value is not None else 15,
                    mm_w_st=self.mm_w_st.value if hasattr(self, "mm_w_st") and self.mm_w_st.value is not None else 13,
                    mm_r_st=self.mm_r_st.value if hasattr(self, "mm_r_st") and self.mm_r_st.value is not None else 0,
                    mm_r_sf=self.mm_r_sf.value if hasattr(self, "mm_r_sf") and self.mm_r_sf.value is not None else 0,
                    mm_r_sb=self.mm_r_sb.value if hasattr(self, "mm_r_sb") and self.mm_r_sb.value is not None else 0,
                    Q=self.Q.value if hasattr(self, "Q") and self.Q.value is not None else 6,
                ).draw(drawer, **kwargs),
                _calculate_points=lambda: None,
            ),
            "coils": Geometry(
                GP=[self.mm_r_so, self.mm_d_sy, self.mm_w_st, self.mm_d_st],
                draw_function=lambda drawer, **kwargs: CrossSectStator.CrossSectInnerRotorStatorWinding(
                    mm_r_si=self.mm_r_si.value if hasattr(self, "mm_r_si") and self.mm_r_si.value is not None else 35,
                    mm_d_st=self.mm_d_st.value if hasattr(self, "mm_d_st") and self.mm_d_st.value is not None else 5,
                    mm_d_sy=self.mm_d_sy.value if hasattr(self, "mm_d_sy") and self.mm_d_sy.value is not None else 6,
                    mm_w_st=self.mm_w_st.value if hasattr(self, "mm_w_st") and self.mm_w_st.value is not None else 4,
                    mm_d_stt=self.mm_d_sts.value if hasattr(self, "mm_d_sts") and self.mm_d_sts.value is not None else 2,
                    Q=self.Qs.value if hasattr(self, "Qs") and self.Qs.value is not None else 12,
                ).draw(drawer, **kwargs),
                _calculate_points=lambda: None
            )
        }


        # 从参数中获取winding参数值，如果没有则使用默认值
        m_val = self.m.value if self.m.value is not None else 3
        Qs_val = self.Qs.value if self.Qs.value is not None else 12
        p_val = self.p.value if self.p.value is not None else 4
        ps_val = self.ps.value if self.ps.value is not None else 5
        self.wily = Winding(phase_number_m=m_val, stator_slot_number_Qs=Qs_val, pole_pair_number_p=p_val, suspension_pole_pair_number_ps=ps_val)

        # 以下代码用于实际创建电机设计变体，但在仅需要参数信息时可以跳过
        # 如果 fea_config_dict 和 spec_input_dict 不存在，说明这是用于参数配置的场景，跳过实际设计创建
        if hasattr(self, 'fea_config_dict') and hasattr(self, 'spec_input_dict'):
            try:
                import bearingless_spmsm_heart
                acm_template = bearingless_spmsm_heart.bearingless_spmsm_template(self.fea_config_dict, self.spec_input_dict)
                # 注意：x_denorm, counter, counter_loop 等变量需要在使用前定义
                # acm_variant = bearingless_spmsm_heart.bearingless_spmsm_design_variant(template=acm_template, x_denorm=x_denorm, counter=counter, counter_loop=counter_loop)
                # acm_variant = self.ad.build_acm_variant(self.ad.acm_template, x_denorm, counter=counter) # counter has the same function as filename
            except Exception as e:
                # 如果创建失败，不影响参数配置功能
                pass

    def show_geometry(self) -> None:
        import VanGogh_Cairo
        toolCairo = VanGogh_Cairo.VanGogh_Cairo(acm_variant, width_in_points=acm_variant.template.SI['GP']['mm_r_so'].value*2.1, 
                                                            height_in_points=acm_variant.template.SI['GP']['mm_r_so'].value*2.1,
                                                            filename=filename)
        if 'PMSM' in acm_variant.template.name:
            lw = 0.1 if acm_variant.template.SI['GP']['mm_r_ro'].value < 15 else 0.5
            saved_filename = toolCairo.draw_spmsm(acm_variant, bool_draw_whole_model=True, lw=lw)
            return saved_filename

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
        将 Modern_Machine_Designer 对象转换为字典（用于 JSON 序列化）
        
        Returns:
            Dict: 包含所有字段的字典
        """
        result = {
            'machine_class': self.machine_class,
            'parameters': {}
        }
        
        # 序列化所有 Parameter 字段
        for field_name, param in self.get_parameter_fields().items():
            result['parameters'][field_name] = param.to_dict()
        
        # 序列化 Winding 对象（如果存在）
        if hasattr(self, 'wily') and self.wily is not None:
            result['winding'] = self.wily.to_dict()
        
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
        instance.machine_class = data.get('machine_class', field_defaults.get('machine_class', cls.machine_class))
        
        # 设置所有参数字段
        for field_name, param in param_dict.items():
            setattr(instance, field_name, param)
        
        # 重建 Winding 对象
        if 'winding' in data:
            instance.wily = Winding.from_dict(data['winding'])
        else:
            # 如果没有 winding 数据，调用 __post_init__ 来创建
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

