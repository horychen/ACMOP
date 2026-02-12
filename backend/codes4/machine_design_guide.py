import json, math, base64, pickle, cairo, os, jsonpickle, logging, utility, JMAG, builtins
try:
    import pygmo as pg
except ImportError:
    from unittest.mock import MagicMock
    pg = MagicMock()
from dataclasses import dataclass, fields, field
from typing import Dict, List, Optional, Any
from collections import OrderedDict
from time import time as clock_time
from modern_machine_designer_utility import Modern_Machine_Designer_Utility, Swarm_Data_Analyzer, swarm_data_container, Parameter, Geometry, Winding, CairoDrawer
from user_minitureMachine import MotorSpecs

# Global verbose control for drawing operations
# Set this to True to enable all print statements in CrossSect classes
builtins.VERBOSE_DRAWING = False  # Default to False, can be changed in __post_init__ or elsewhere


def _default_machine_input():
    """Default machine input from user_minitureMachine: MotorSpecs()."""
    return MotorSpecs()


@dataclass
class Modern_Machine_Designer(Modern_Machine_Designer_Utility):

    def __init__(self, specs: MotorSpecs, machine_class: str = 'dex13', select_fea_config_dict: str = '#0301 JMAG Non-Bearingless', verbose_drawing: bool = False):
        """
        Initialize the Modern_Machine_Designer with MotorSpecs.
        """
        super().__init__(specs)
        self.machine_class = machine_class
        self.select_fea_config_dict = select_fea_config_dict
        self.verbose_drawing = verbose_drawing
        
        # Internal control
        self.counter = 0
        self.bool_jmagDeleteResultsAfterCalculation = False
        self.name = 'SPMSM'
        
        # Run initialization logic
        self._init_components()

    def get_path2SwarmData(self, folder_name):
        self.dir_parent = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/'
        self.dir_codes  = os.path.abspath(os.path.dirname(__file__)) + '/'

        self.path2SwarmData = os.path.abspath(fr'../_default/' + folder_name.replace(' ', '_'))
        self.swarm_data_json_file_path = self.path2SwarmData + f'/SwarmData.json'
        if not os.path.isdir(self.path2SwarmData): os.makedirs(self.path2SwarmData)
        os.chdir(self.dir_codes)

    def _init_components(self):
        """
        Setup components, parameters, and FEA configuration.
        """
        # Set global verbose control for drawing operations
        builtins.VERBOSE_DRAWING = self.verbose_drawing

        # Ensure specs are synced
        self.specs.sync()
        specs = self.specs
        
        # Unpack essential winding/geometry info from specs for compatibility
        m = specs.winding.m
        Qs = specs.winding.num_slots
        p = specs.winding.num_poles // 2
        ps = 4 # Default placeholder for ps if not in specs
        coil_pitch_y = specs.winding.coil_pitch_y

        RatedPower = specs.winding.rated_power
        RatedSpeed = specs.winding.rated_speed

        # Setup paths and FEA config
        suffix = f'minitureMachine'
        if not self.machine_class.endswith(suffix):
            self.machine_class = self.machine_class + suffix

        self.get_path2SwarmData(self.machine_class)

        machine_sim_json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'machine_simulation.json')
        with open(machine_sim_json_path, 'r') as f:
            raw_fea_config_dicts = json.load(f)
            self.fea_config_dict = OrderedDict(raw_fea_config_dicts[self.select_fea_config_dict])
            self.fea_config_dict['pc_name'] = self.get_pc_name()

        self.wily = specs.winding.wily

        # Legacy flags and tool selection from specs
        self.bool_PermanentMagnet = specs.geometry.bool_PermanentMagnet
        self.bool_StatorSlotClosed = specs.geometry.bool_StatorSlotClosed
        self.bool_RotorNotched = specs.geometry.bool_RotorNotched
        self.select_FEA_tool = specs.targets.select_FEA_tool
        self.bool_jmagDeleteResultsAfterCalculation = specs.targets.bool_jmagDeleteResultsAfterCalculation

        # Re-attach Parameters to self for existing code to work
        g = specs.geometry
        self.m = Parameter('phase_number_m', 'fixed', m)
        self.Qs = Parameter('stator_slot_number_Qs', 'fixed', Qs)
        self.p = specs.geometry.p
        self.ps = Parameter('suspension_pole_pair_number_ps', 'fixed', ps)
        self.coil_pitch_y = Parameter('coil_pitch_y', 'fixed', coil_pitch_y)
        self.mm_r_so = g.r_stator_outer
        self.mm_d_mech_air_gap = g.d_air_gap
        self.mm_d_sleeve = g.d_sleeve
        self.mm_r_ri = g.r_rotor_inner
        self.mm_w_st = g.w_tooth
        self.mm_d_sy = g.d_stator_yoke
        self.split_ratio = g.split_ratio
        self.mm_d_pm = g.d_magnet
        self.deg_alpha_rm = g.alpha_magnet_pole
        # Attach these others which might be needed
        self.mm_d_sts = g.d_tooth_shoe
        self.deg_alpha_st = g.tooth_specs.alpha_tooth if hasattr(g.tooth_specs, 'alpha_tooth') else Parameter("stator_tooth_span_angle", "fixed", 0.0)
        self.s = specs.geometry.s

        # Now update dicts
        self.parameter_dict = self.get_parameter_dict_by_name()
        self.parameter_dict_by_name = self.get_parameter_dict_by_name()

        self.machineGeometry = OrderedDict()
        
        # Derived variables - Re-attach these using legacy naming and lambdas
        p_dict = self.parameter_dict_by_name
        self.mm_r_si = Parameter('stator_inner_radius', 'derived', calc=lambda d: d['stator_outer_radius'].value * d['split_ratio_r_si_slash_r_so'].value, parameter_dict=p_dict)
        p_dict['stator_inner_radius'] = self.mm_r_si
        
        self.mm_d_st = Parameter('stator_tooth_depth', 'derived', calc=lambda d: d['stator_outer_radius'].value - d['stator_inner_radius'].value - d['stator_yoke_depth'].value - d['stator_tooth_shoe_depth'].value, parameter_dict=p_dict)
        p_dict['stator_tooth_depth'] = self.mm_d_st
        
        self.mm_r_ro = Parameter('rotor_outer_radius', 'derived', calc=lambda d: d['stator_inner_radius'].value - d['mechanical_air_gap_depth'].value - d['rotor_sleeve_depth'].value, parameter_dict=p_dict)
        p_dict['rotor_outer_radius'] = self.mm_r_ro
        
        self.mm_d_ri = Parameter('rotor_iron (back iron) depth', 'derived', calc=lambda d: d['rotor_outer_radius'].value - d['magnet_depth'].value, parameter_dict=p_dict)
        p_dict['rotor_iron (back iron) depth'] = self.mm_d_ri
        
        self.mm_d_sto = Parameter('stator_tooth_open_depth', 'derived', calc=lambda d: d['stator_tooth_shoe_depth'].value * 0.667, parameter_dict=p_dict)
        p_dict['stator_tooth_open_depth'] = self.mm_d_sto
        
        self.deg_alpha_sto = Parameter('stator_tooth_open_angle', 'derived', calc=lambda d: (d['stator_tooth_span_angle'].value if d['stator_tooth_span_angle'].value is not None else 0) * 0.5, parameter_dict=p_dict)
        p_dict['stator_tooth_open_angle'] = self.deg_alpha_sto
        
        self.deg_alpha_rs = Parameter('magnet_segment_span_angle', 'derived', calc=lambda d: d['magnet_pole_span_angle'].value, parameter_dict=p_dict)
        p_dict['magnet_segment_span_angle'] = self.deg_alpha_rs
        
        self.mm_d_rp = Parameter('inter_polar_iron_thickness', 'derived', calc=lambda d: d['magnet_depth'].value, parameter_dict=p_dict)
        p_dict['inter_polar_iron_thickness'] = self.mm_d_rp
        
        self.mm_d_rs = Parameter('inter_segment_iron_thickness', 'fixed', 0.0)
        p_dict['inter_segment_iron_thickness'] = self.mm_d_rs

        # Refresh dicts again after adding derived
        self.parameter_dict = self.get_parameter_dict_by_name()
        self.parameter_dict_by_name = self.get_parameter_dict_by_name()

        # Update all derived values exactly once
        self.update_geometric_parameters()

        self.InitialRotationAngle = specs.targets.initial_rotation_angle = self.get_InitialRotationAngle()

        # [Important] Print summary
        specs.print_summary()

        EX = {
            # 3D
            'mm_stack_length_specified': specs.geometry.l_stack.value, # mm
            # Materials
            'Magnet_Name': specs.materials.magnet_grade,
            'Magnet_StartAngle': 0.5* 360/(2*p),
            'Magnet_Temperature': specs.materials.magnet_temperature,
            'SteelMaterial': specs.materials.stator_core_material, 
            'StatorCore_Material': specs.materials.stator_core_material, 
            'RotorCore_Material': specs.materials.rotor_core_material,
            'LaminationFactor': specs.materials.steel_stack_factor * 100.0,
            # Thermal
            'RatedPower': specs.winding.rated_power,
            'RatedSpeed': specs.winding.rated_speed,
            'ExcitationFreqSimulated': specs.winding.rated_speed / 60.0 * p,
            'bool_WyeConnectOrDeltaConnect' : specs.winding.bool_WyeConnectOrDeltaConnect,
            'DCBusVoltage': specs.winding.dc_bus_voltage,
            'Js' : specs.winding.rated_current_density_Js,
            'Temperature': specs.materials.magnet_temperature, # Using magnet temperature as general temperature
            'WindingFill': specs.winding.winding_fill_factor,
            'TORQUE_CURRENT_RATIO': specs.winding.torque_current_ratio,
            'SUSPENSION_CURRENT_RATIO': specs.winding.suspension_current_ratio,
            'DriveW_Rs': specs.winding.drive_winding_resistance, # [Ohm]
            'BeariW_Rs': specs.winding.bearing_winding_resistance, # [Ohm]
        }
        self.EX = EX

        # Update EX with values from specs.targets
        EX['no_series_coil_turns_N'] = specs.targets.no_series_coil_turns_N
        EX['DriveW_zQ'] = specs.targets.no_conductors_per_slot_zQ
        EX['BeariW_zQ'] = specs.targets.no_conductors_per_slot_zQ if specs.winding.wily.bool_DPNVorSEPA == True else specs.targets.no_conductors_per_slot_zQ / specs.winding.torque_current_ratio * specs.winding.suspension_current_ratio
        EX['mm2_slot_area'] = specs.targets.mm2_slot_area
        EX['CurrentAmp_in_the_slot'] = specs.targets.CurrentAmp_in_the_slot
        EX['CurrentAmp_per_conductor'] = specs.targets.CurrentAmp_per_conductor
        EX['CurrentAmp_per_phase'] = specs.targets.CurrentAmp_per_phase
        EX['DriveW_CurrentAmp'] = specs.targets.DriveW_CurrentAmp
        EX['BeariW_CurrentAmp'] = specs.targets.BeariW_CurrentAmp
        EX['slot_current_utilizing_ratio_for_torque'] = specs.targets.slot_current_utilizing_ratio_for_torque
        EX['mm2_magnet_area'] = specs.targets.mm2_magnet_area

        EX['InitialRotationAngle'] = self.get_InitialRotationAngle()

        Rout = self.mm_r_ri.value+ self.mm_d_ri.value+ self.mm_d_pm.value
        Rin  = self.mm_r_ri.value+ self.mm_d_ri.value
        deg_alpha_rp = 360 / (2*self.p.value)
        EX['mm2_magnet_area'] = self.deg_alpha_rm.value/deg_alpha_rp * math.pi*(Rout**2 - Rin**2)


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
                        mm_d_pm=self.mm_d_pm.value,
                        deg_alpha_rm=self.deg_alpha_rm.value,
                        deg_alpha_rs=self.deg_alpha_rm.value if self.s.value==1 else self.deg_alpha_rs.value,
                        mm_d_ri=self.mm_d_ri.value,
                        mm_r_ri=self.mm_r_ri.value,
                        mm_d_rp=self.mm_d_rp.value,
                        mm_d_rs=self.mm_d_rs.value,
                        p=self.p.value,
                        s=self.s.value
                    ).draw(drawer, **kwargs)
                )
            ),
            # "shaft": Geometry(name='shaft',
            #     GP={'mm_r_ri': self.mm_r_ri},
            #     draw_function=lambda drawer, **kwargs: (
            #         CrossSectInnerNotchedRotor.CrossSectShaft(
            #             name="shaft",
            #             color="#0EE0E2",
            #             rotorCore=CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
            #                 mm_d_pm=self.mm_d_pm.value,
            #                 deg_alpha_rm=self.deg_alpha_rm.value,
            #                 deg_alpha_rs=self.deg_alpha_rm.value if self.s.value==1 else self.deg_alpha_rs.value,
            #                 mm_d_ri=self.mm_d_ri.value,
            #                 mm_r_ri=self.mm_r_ri.value,
            #                 mm_d_rp=self.mm_d_rp.value,
            #                 mm_d_rs=self.mm_d_rs.value,
            #                 p=self.p.value,
            #                 s=self.s.value
            #             )
            #         ).draw(drawer, **kwargs)
            #     )
            # ),
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
                            mm_d_pm=self.mm_d_pm.value,
                            deg_alpha_rm=self.deg_alpha_rm.value,
                            deg_alpha_rs=self.deg_alpha_rm.value if self.s.value==1 else self.deg_alpha_rs.value,
                            mm_d_ri=self.mm_d_ri.value,
                            mm_r_ri=self.mm_r_ri.value,
                            mm_d_rp=self.mm_d_rp.value,
                            mm_d_rs=self.mm_d_rs.value,
                            p=self.p.value,
                            s=self.s.value
                        )
                    ).draw(drawer, **kwargs)
                ),
            ),
            # "sleeve": Geometry(name='sleeve',
            #     GP={
            #         'mm_r_ri': self.mm_r_ri,
            #         'mm_d_ri': self.mm_d_ri,
            #         'mm_d_pm': self.mm_d_pm,
            #         'p': self.p,
            #     },
            #     draw_function=lambda drawer, **kwargs: (
            #         CrossSectInnerNotchedRotor.CrossSectSleeve(
            #             mm_r_ri=self.mm_r_ri.value,
            #             mm_d_ri=self.mm_d_ri.value,
            #             mm_d_pm=self.mm_d_pm.value,
            #             p=self.p.value,
            #             d_sleeve=self.mm_d_sleeve.value
            #         ).draw(drawer, **kwargs)
            #     )
            # ),
            "statorCore": Geometry(name='statorCore',
                GP={
                    'mm_r_so': self.mm_r_so,
                    'mm_r_si': self.mm_r_si,
                    'mm_d_sts': self.mm_d_sts,
                    'mm_d_st': self.mm_d_st,
                    'mm_d_sy': self.mm_d_sy,
                    'mm_w_st': self.mm_w_st,
                    'Q': self.Qs,
                },
                draw_function=lambda drawer, **kwargs: (
                    CrossSectStator.CrossSectInnerRotorClosedSlotStator(
                        name="statorCore",
                        color="#BAFA01",
                        mm_r_so=self.mm_r_so.value,
                        mm_d_sts=self.mm_d_sts.value,
                        mm_r_si=self.mm_r_si.value,
                        mm_d_st=self.mm_d_st.value,
                        mm_d_sy=self.mm_d_sy.value,
                        mm_w_st=self.mm_w_st.value,
                        Q=self.Qs.value,
                    ).draw(drawer, **kwargs)
                ),
            ),
            # "statorCore": Geometry(name='statorCore',
            #     GP={
            #         'mm_r_si': self.mm_r_si,
            #         'mm_d_sto': self.mm_d_sto,
            #         'mm_d_sts': self.mm_d_sts,
            #         'mm_d_st': self.mm_d_st,
            #         'mm_d_sy': self.mm_d_sy,
            #         'mm_w_st': self.mm_w_st,
            #         'deg_alpha_st': self.deg_alpha_st,
            #         'deg_alpha_sto': self.deg_alpha_sto,
            #         'mm_d_sto': self.mm_d_sto,
            #         'Q': self.Qs,
            #     },
            #     draw_function=lambda drawer, **kwargs: (
            #         CrossSectStator.CrossSectInnerRotorStator(
            #             name="statorCore",
            #             color="#BAFA01",
            #             deg_alpha_st=self.deg_alpha_st.value,
            #             deg_alpha_sto=self.deg_alpha_sto.value,
            #             mm_r_si=self.mm_r_si.value,
            #             mm_d_sto=self.mm_d_sto.value,
            #             mm_d_sts=self.mm_d_sts.value,
            #             mm_d_st=self.mm_d_st.value,
            #             mm_d_sy=self.mm_d_sy.value,
            #             mm_w_st=self.mm_w_st.value,
            #             Q=self.Qs.value
            #         ).draw(drawer, **kwargs)
            #     ),
            # ),
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
                # CrossSectStator.CrossSectInnerRotorStatorWinding(
                CrossSectStator.CrossSectInnerRotorClosedSlotStatorWinding(
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
        # print(f"[bearingless_spmsm_design.py] [PMSM JMAG] {self.InitialRotationAngle} deg = ", (deg_pole_span-self.deg_alpha_rm.value)*0.5,  self.wily.deg_winding_U_phase_phase_axis_angle,  deg_pole_span*0.5)
        # print(f"[bearingless_spmsm_design.py] [PMSM JMAG] {self.InitialRotationAngle} deg")

        return self.InitialRotationAngle

    def show_geometry(self, filename=None, x_denorm_dict=None) -> None:
        if x_denorm_dict is not None:
            self.update_geometric_parameters(x_denorm_dict=x_denorm_dict)

        bool_draw_whole_model = True
        
        def draw_spmsm(lw, width_in_points, height_in_points, filename='machine_geometry.svg', bool_draw_whole_model=True):
            self.drawer = drawer = CairoDrawer(width_in_points, height_in_points, filename=filename, verbose_drawing=getattr(self, 'verbose_drawing', False))

            # 直接调用 draw 方法，如果 machineGeometry 不存在或缺少必要的键，会直接报错
            list_regions = self.machineGeometry['rotorCore'].draw(drawer, bool_draw_whole_model=bool_draw_whole_model)
            # list_regions = self.machineGeometry['shaft'].draw(drawer)
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
        lw = 0.1 if self.mm_r_ro.value < 15 else 0.5
        width_in_points  = self.mm_r_so.value*2.1
        height_in_points = self.mm_r_so.value*2.1
        draw_spmsm(lw, width_in_points, height_in_points, filename=filename or 'machine_geometry.svg', bool_draw_whole_model=bool_draw_whole_model)

    def FEA_evaluate(self, project_loc=fr'../_default/', bool_jmagDesignerShow: bool = True, x_denorm=None, counter=None, counter_loop=0):

        self.update_geometric_parameters(x_denorm=x_denorm)

        if counter is None:
            counter = self.counter
        else:
            self.counter = counter

        # define project_name using counter and counter_loop
        self.project_name = self.name + f'-ind{counter}'
        self.project_name += f'-redo{counter_loop}' if counter_loop > 1 else ''

        jmag_temp_dir = self.path2SwarmData + "/jmag_temp"
        if not os.path.exists(jmag_temp_dir):
            os.makedirs(jmag_temp_dir)
        jmag_screenshots_dir = self.path2SwarmData + "/jmag_screenshots"
        if not os.path.exists(jmag_screenshots_dir):
            os.makedirs(jmag_screenshots_dir)
        self.expected_project_file = jmag_temp_dir + "/%s.jproj" % (self.project_name)


        self.path2FEACsv = os.path.abspath(self.path2SwarmData + '/csv/') + f'/{self.counter}/'
        if not os.path.isdir(self.path2FEACsv): os.makedirs(self.path2FEACsv)

        if 'JMAG' in self.select_FEA_tool:
            study_name = "Transient"

            # Leave the solving task to JMAG
            def build_jmag_project(study_name):
                self.toolJd = JMAG.JMAG(fea_config_dict=self.fea_config_dict)
                self.toolJd.open(Steel_name=self.EX['SteelMaterial'], expected_project_file_path=self.expected_project_file, pc_name=self.fea_config_dict['pc_name'], dir_parent=self.dir_parent, bool_jmagDesignerShow=bool_jmagDesignerShow)
                return self.toolJd

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
                if not self.bool_StatorSlotClosed:
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
                if not self.bool_StatorSlotClosed:
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

            if 'PMSM' in self.machine_class:
                draw_spmsm(self.toolJd)

            # JMAG
            app = toolJd.app
            model = app.GetModel(self.name)

            if 'PMSM' in self.machine_class:
                toolJd.pre_process_PMSM(app, model, self)

            study = toolJd.add_magnetic_transient_study(app, model, self.path2FEACsv, study_name, self)
            toolJd.mesh_study(self, app, model, study, output_dir=self.path2SwarmData)
            # raise KeyboardInterrupt
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

                # import rich
                # rich.print(self.results_to_be_unpacked)

                # acm_variant.spec_geometry_dict['x_denorm'] = list(x_denorm)

                self.spec_performance_dict = spec_performance_dict = dict()
                spec_performance_dict['x_denorm_dict'] = dict(self.get_free_variables_as_dict()) # Convert OrderedDict to dict to be serelized to json file
                spec_performance_dict['project_name'] = project_name
                spec_performance_dict['individual_name'] = individual_name
                spec_performance_dict['number_current_generation'] = number_current_generation
                spec_performance_dict['individual_index'] = individual_index
                # spec_performance_dict['cost_function'] = cost_function
                spec_performance_dict['f1'] = float(f1)
                spec_performance_dict['f2'] = float(f2)
                spec_performance_dict['f3'] = float(f3)
                spec_performance_dict['TRV'] = float(TRV)
                spec_performance_dict['FRW'] = float(FRW)
                spec_performance_dict['torque_average'] = torque_average
                spec_performance_dict['ss_avg_force_magnitude'] = ss_avg_force_magnitude
                spec_performance_dict['rotor_weight'] = float(rotor_weight)
                spec_performance_dict['normalized_torque_ripple'] = float(normalized_torque_ripple)
                spec_performance_dict['normalized_force_error_magnitude'] = float(normalized_force_error_magnitude)
                spec_performance_dict['force_error_angle'] = float(force_error_angle)
                spec_performance_dict['coil_flux_linkage_peak2peak_value'] = float(coil_flux_linkage_peak2peak_value)
                spec_performance_dict['mm2_slot_area'] = mm2_slot_area
                spec_performance_dict['Cost'] = float(Cost)
                spec_performance_dict['Cost_Fe'] = float(Cost_Fe)
                spec_performance_dict['Cost_Cu'] = float(Cost_Cu)
                spec_performance_dict['Cost_PM'] = float(Cost_PM)
                spec_performance_dict['power_factor'] = power_factor
                spec_performance_dict['rated_ratio'] = rated_ratio
                spec_performance_dict['rated_stack_length_mm'] = rated_stack_length_mm
                spec_performance_dict['rated_total_loss'] = float(rated_total_loss)
                spec_performance_dict['rated_stator_copper_loss_along_stack'] = rated_stator_copper_loss_along_stack
                spec_performance_dict['rated_rotor_copper_loss_along_stack'] = rated_rotor_copper_loss_along_stack
                spec_performance_dict['rated_magnet_Joule_loss'] = rated_magnet_Joule_loss
                spec_performance_dict['stator_copper_loss_in_end_turn'] = float(stator_copper_loss_in_end_turn)
                spec_performance_dict['rotor_copper_loss_in_end_turn']  = float(rotor_copper_loss_in_end_turn)
                spec_performance_dict['rated_iron_loss'] = rated_iron_loss
                spec_performance_dict['rated_windage_loss'] = float(rated_windage_loss)
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
                results2file = {f'spec_performance_dict-gen{number_current_generation}-ind{individual_index}' : spec_performance_dict}

                # INSERT_YOUR_CODE
                # Load existing JSON file (if any), append the new results and write back using jsonpickle.
                try:
                    # Try loading the existing JSON file
                    if os.path.exists(self.swarm_data_json_file_path) and os.path.getsize(self.swarm_data_json_file_path) > 0:
                        with open(self.swarm_data_json_file_path, 'r') as rf:
                            existing_data = jsonpickle.decode(rf.read())
                            if not isinstance(existing_data, dict):
                                existing_data = {}
                    else:
                        existing_data = {}
                except Exception:
                    existing_data = {}

                # Add/overwrite the results for the current generation/individual
                existing_data.update(results2file)

                # Write the updated data back to the file
                json_string = jsonpickle.encode(existing_data, indent=4)
                with open(self.swarm_data_json_file_path, 'w') as wf:
                    wf.write(json_string)

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

    def evaluate_design_json_wrapper(self, x_denorm, counter=1, counter_loop=1):
        # This is a wrapper for FEA, in order to build up a json profile for the design variant
        builtins.ad = self  # share global variable between modules
        ad = self

        # 这里应该返回新获得的设计，然后可以获得geometry_dict，然后包括x_denorm的信息方便重构设计。
        self.FEA_evaluate(bool_jmagDesignerShow=self.fea_config_dict["designer.Show"], x_denorm=x_denorm, counter=counter, counter_loop=counter_loop)
        # print('[DEBUG] evaluate_design_json_wrapper: counter =', counter, 'counter_loop =', counter_loop)

        if 'FEMM' in self.select_fea_config_dict:
            self.results_for_optimization = self.analyzer.build_results_for_optimization()

            # Save spec_performance_dict and others to disk
            GP = self.template.d['GP']
            EX = self.template.d['EX']
            self.save_to_disk(self, self.analyzer.spec_performance_dict, GP, EX)

            # Save also the object (self) to disk, but this takes a lot of disk space!
            if self.fea_config_dict['moo.save_self_object_as_jsonpickle'] == True:
                utility_json.to_json_recursively(self, self.name, save_here=self.fea_config_dict['output_dir']+'jsonpickle/')

            # Save time domain data to disk
            self.analyzer.save_time_domain_data(self.fea_config_dict['output_dir']+self.select_spec+f'-ind{counter}.pkl') # counter could be string

            return self

        elif 'JMAG' in self.select_fea_config_dict:

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

            # self.spec_geometry_dict['x_denorm'] = list(x_denorm)

            # spec_performance_dict = dict()
            # spec_performance_dict['x_denorm_dict'] = self.get_free_variables_as_dict()
            # spec_performance_dict['project_name'] = project_name
            # spec_performance_dict['individual_name'] = individual_name
            # spec_performance_dict['number_current_generation'] = number_current_generation
            # spec_performance_dict['individual_index'] = individual_index
            # # spec_performance_dict['cost_function'] = cost_function
            # spec_performance_dict['f1'] = f1
            # spec_performance_dict['f2'] = f2
            # spec_performance_dict['f3'] = float(f3)
            # spec_performance_dict['TRV'] = TRV
            # spec_performance_dict['FRW'] = FRW
            # spec_performance_dict['torque_average'] = torque_average
            # spec_performance_dict['ss_avg_force_magnitude'] = ss_avg_force_magnitude
            # spec_performance_dict['rotor_weight'] = rotor_weight
            # spec_performance_dict['normalized_torque_ripple'] = float(normalized_torque_ripple)
            # spec_performance_dict['normalized_force_error_magnitude'] = float(normalized_force_error_magnitude)
            # spec_performance_dict['force_error_angle'] = float(force_error_angle)
            # spec_performance_dict['coil_flux_linkage_peak2peak_value'] = float(coil_flux_linkage_peak2peak_value)
            # spec_performance_dict['mm2_slot_area'] = mm2_slot_area
            # spec_performance_dict['Cost'] = Cost
            # spec_performance_dict['Cost_Fe'] = Cost_Fe
            # spec_performance_dict['Cost_Cu'] = Cost_Cu
            # spec_performance_dict['Cost_PM'] = Cost_PM
            # spec_performance_dict['power_factor'] = power_factor
            # spec_performance_dict['rated_ratio'] = rated_ratio
            # spec_performance_dict['rated_stack_length_mm'] = rated_stack_length_mm
            # spec_performance_dict['rated_total_loss'] = rated_total_loss
            # spec_performance_dict['rated_stator_copper_loss_along_stack'] = rated_stator_copper_loss_along_stack
            # spec_performance_dict['rated_rotor_copper_loss_along_stack'] = rated_rotor_copper_loss_along_stack
            # spec_performance_dict['rated_magnet_Joule_loss'] = rated_magnet_Joule_loss
            # spec_performance_dict['stator_copper_loss_in_end_turn'] = stator_copper_loss_in_end_turn
            # spec_performance_dict['rotor_copper_loss_in_end_turn'] = rotor_copper_loss_in_end_turn
            # spec_performance_dict['rated_iron_loss'] = rated_iron_loss
            # spec_performance_dict['rated_windage_loss'] = rated_windage_loss
            # # spec_performance_dict['str_results'] = str_results
            # spec_performance_dict['select_fea_config_dict'] = self.select_fea_config_dict
            # spec_performance_dict['moo.fitness_OA'] = self.fea_config_dict['moo.fitness_OA']
            # spec_performance_dict['moo.fitness_OB'] = self.fea_config_dict['moo.fitness_OB']
            # spec_performance_dict['moo.fitness_OC'] = self.fea_config_dict['moo.fitness_OC']

            # # Save to disk
            # # self.save_to_disk(self, spec_performance_dict, GP, EX)

            # number_current_generation = spec_performance_dict['number_current_generation'] #= int(self.counter//popsize), 
            # individual_index = spec_performance_dict['individual_index'] #= self.counter
            # # builtins.ad.visualize_dict[f'FEA_Evaluated_Performance-{number_current_generation}-{individual_index}'] = spec_performance_dict


            # Read the possibly-existing current json data
            try:
                if os.path.getsize(self.swarm_data_json_file_path) > 0:
                    with open(self.swarm_data_json_file_path, 'r') as rf:
                        loaded_json = json.load(rf)
                else:
                    loaded_json = {}
            except Exception:
                loaded_json = {}

            # Compose new key
            # key = f'gen{number_current_generation}-ind{individual_index}'
            # loaded_json[key] = builtins.ad.visualize_dict

            json_string = jsonpickle.encode(loaded_json, indent=4)
            with open(self.swarm_data_json_file_path, 'w+') as f:
                f.write(json_string)

            # save object (self) to disk
            # utility_json.to_json_recursively(self, self.name, save_here=self.fea_config_dict['output_dir']+'jsonpickle/')

            # this is for optimization
            self.results_for_optimization = (cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle)

            # this is for comparison to FEMM
            def compare_with_FEMM(self):
                EX = self.template.d['EX']
                self.analyzer = FEMM_SlidingMesh.Individual_Analyzer_FEMM_Edition(p=EX['wily'].p)
                basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list, \
                    DisplacementAngle_list, \
                    circuit_current_GroupACU, \
                    circuit_current_GroupACV, \
                    circuit_current_GroupACW, \
                    circuit_current_GroupBDU, \
                    circuit_current_GroupBDV, \
                    circuit_current_GroupBDW, \
                    terminal_voltage_GroupACU, \
                    terminal_voltage_GroupACV, \
                    terminal_voltage_GroupACW, \
                    terminal_voltage_GroupBDU, \
                    terminal_voltage_GroupBDV, \
                    terminal_voltage_GroupBDW, \
                    coil_fluxLinkage_GroupACU, \
                    coil_fluxLinkage_GroupACV, \
                    coil_fluxLinkage_GroupACW, \
                    coil_fluxLinkage_GroupBDU, \
                    coil_fluxLinkage_GroupBDV, \
                    coil_fluxLinkage_GroupBDW = self.toolJd.dm.unpack(bool_more_info=True)
                electrical_period = self.template.fea_config_dict['designer.number_cycles_in_2ndTSS']/EX['DriveW_Freq']
                number_of_steps   = self.template.fea_config_dict['designer.number_of_steps_2ndTSS']
                step_size_sec = electrical_period / number_of_steps
                step_size_mech_deg = EX['Omega'] * step_size_sec / math.pi * 180

                for index in range(-self.toolJd.dm.number_of_steps_at_steady_state, 0):
                    time                         = float(time_list[index])
                    RotorAngle_MechanicalDegrees = float(DisplacementAngle_list[index])
                    torque = float(TorCon_list[index])
                    forces = ( float(ForConX_list[index]), float(ForConY_list[index]) )
                    energy = 0.0
                    circuitProperties = ( [ circuit_current_GroupACU[index], terminal_voltage_GroupACU[index], coil_fluxLinkage_GroupACU[index] ],
                                        [ circuit_current_GroupACV[index], terminal_voltage_GroupACV[index], coil_fluxLinkage_GroupACV[index] ],
                                        [ circuit_current_GroupACW[index], terminal_voltage_GroupACW[index], coil_fluxLinkage_GroupACW[index] ],
                                        [ circuit_current_GroupBDU[index], terminal_voltage_GroupBDU[index], coil_fluxLinkage_GroupBDU[index] ],
                                        [ circuit_current_GroupBDV[index], terminal_voltage_GroupBDV[index], coil_fluxLinkage_GroupBDV[index] ],
                                        [ circuit_current_GroupBDW[index], terminal_voltage_GroupBDW[index], coil_fluxLinkage_GroupBDW[index] ] )
                    self.analyzer.add(time, RotorAngle_MechanicalDegrees, torque, forces, energy, circuitProperties)
                self.analyzer.get_ss_data()

            # compare_with_FEMM(self)
            # self.analyzer.save_time_domain_data(counter) # TODO

            return self.results_for_optimization

    def start_optimization(self, bool_local_exploration_around_selected_individual=False, std_fraction=0.05):
        """
        启动多目标优化过程
        
        主要功能：
        1. 读取现有的群体数据（从 SwarmData.json）
        2. 评估帕累托前沿并基于拥挤距离选择最优个体
        3. 初始化种群并开始优化迭代
        """
        builtins.ad = self  # share global variable between modules
        ad = self

        # 设置路径和日志
        self.logger = self.myLogger(self.path2SwarmData + '/', prefix=self.machine_class)
        logger = logging.getLogger(__name__)

        self.nobj = sum(1 for k, v in self.fea_config_dict.items() if k.startswith("moo.fitness") and v is not None)
        if self.nobj!=3: raise Exception(f'Number of objectives is {self.nobj} instead of 3. NOT SUPPORTED!!!')
        self.obj_names = [v for k, v in self.fea_config_dict.items() if k.startswith("moo.fitness") and v is not None]
        logger.info(f'Number of objectives is {self.nobj}')
        logger.info(f'Objectives names are {self.obj_names}')

        ################################################################
        # MOO Step 1: 创建问题并初始化种群
        ################################################################
        # [4.3.1] 基本设置
        import Problem_BearinglessSynchronousDesign
        _, prob = Problem_BearinglessSynchronousDesign.get_prob()
        popsize = self.fea_config_dict["moo.popsize"]
        logger.info(f'Pop size is {popsize}')

        # [4.3.2] 准备设计参数信息
        self.x_denorm = list(self.get_free_variables_as_dict().values())
        self.x_denorm_dict = self.get_free_variables_as_dict()
        self.bounds_denorm = list(self.get_free_variable_bounds_dict().values())

        logger.info('bool_local_exploration_around_selected_individual: %s', bool_local_exploration_around_selected_individual)
        for index, (k, v) in enumerate(self.x_denorm_dict.items()):
            logger.info(f'x_denorm_dict variable no. {index} is {k} = {v} with bounds: {self.bounds_denorm[index]}')

        # [4.3.3] 读取现有的群体数据
        logger.info(f'Reading swarm data from: {self.swarm_data_json_file_path}')
        self.analyzer = Swarm_Data_Analyzer(self.swarm_data_json_file_path, self.x_denorm_dict)
        self.swarm_data = self.analyzer.swarm_data_xf

        number_of_chromosome = self.analyzer.number_of_chromosome


        # normal optimization
        # [4.3.4] 根据是否有现有数据来决定初始化策略
        if number_of_chromosome != 0:

            # 初始化种群（此时所有个体的fitness都是[0,0,0]）
            # 禁止在初始化pop时运行有限元
            ad.flag_do_not_evaluate_when_init_pop = True
            pop = pg.population(prob, size=popsize)
            ad.flag_do_not_evaluate_when_init_pop = False

            # Case 1: 存在现有数据 - 从档案中选择最优个体
            logger.info(f'Found {number_of_chromosome} chromosomes in archive. This is a restart.')
            # 设置计数器
            ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome
            logger.info('ad.counter_fitness_called = ad.counter_fitness_return = number_of_chromosome = %d', number_of_chromosome)

            if number_of_chromosome <= popsize:
                pop_array = pop.get_x()
                # 个体数不够一代的情况
                for i in range(popsize):
                    if i < number_of_chromosome:
                        pop.set_xf(i, ad.swarm_data[i][:-3], ad.swarm_data[i][-3:])
                    else:
                        logger.info('Set "ad.flag_do_not_evaluate_when_init_pop" to False...')
                        ad.flag_do_not_evaluate_when_init_pop = False
                        logger.info('Calling pop.set_x()---this is a restart for individual#%d during pop initialization.', i)
                        logger.info('i=%d: call get_fevals: %s', i, prob.get_fevals())
                        pop.set_x(i, pop_array[i])  # evaluate this guy
            else:
                # 使用 learn_about_the_archive 从档案中选择具有高拥挤距离的个体
                logger.info('Using learn_about_the_archive to select individuals with high crowding distance from archive.')
                swarm_data_selected = self.learn_about_the_archive(prob, ad.swarm_data, popsize)
                
                # 将选中的个体设置到种群中
                for i in range(popsize):
                    pop.set_xf(i, swarm_data_selected[i][:-3], swarm_data_selected[i][-3:])
                
                logger.info(f'Selected {popsize} individuals from archive based on Pareto front and crowding distance.')


            number_of_finished_iterations = number_of_chromosome // popsize + 1
        else:
            # Case 2: 没有现有数据 - 全新运行
            number_of_finished_chromosome_in_current_generation = None
            number_of_finished_iterations = 0

            def get_pop_array_around_selected_individual(std_fraction):
                import numpy as np
                param_names = list(self.x_denorm_dict.keys())
                param_values = np.array(list(self.x_denorm_dict.values()), dtype=np.float64)
                # std_fraction = 0.2  # 20% of each parameter, can be adjusted by user
                min_std = 1e-6       # Minimum std to avoid 0 std

                # Get parameter bounds if available
                if hasattr(ad, "get_free_variable_bounds_dict"):
                    bounds_dict = ad.get_free_variable_bounds_dict()
                    lower_bounds = np.array([bounds_dict[name][0] for name in param_names])
                    upper_bounds = np.array([bounds_dict[name][1] for name in param_names])
                else:
                    lower_bounds = None
                    upper_bounds = None

                stds = np.maximum(np.abs(param_values) * std_fraction, min_std)


                for name, value, std, lower_bound, upper_bound in zip(param_names, param_values, stds, lower_bounds, upper_bounds):
                    logger.info(f"{name}: {value}, std: {std}, lower_bound: {lower_bound}, upper_bound: {upper_bound}")
                # INSERT_YOUR_CODE
                if lower_bounds is not None and upper_bounds is not None:
                    for name, value, lower_bound, upper_bound in zip(param_names, param_values, lower_bounds, upper_bounds):
                        if value < lower_bound or value > upper_bound:
                            logger.warning(f"Parameter '{name}' with value {value} is not within bounds [{lower_bound}, {upper_bound}]")

                # Draw popsize = popsize random samples via local exploration
                pop_array = np.zeros((popsize, len(param_values)), dtype=np.float64)
                pop_array[0, :] = param_values # the first individual is the selected individual

                logging.info('generate pop around the selected individual')
                for i in range(1, popsize):
                    candidate = np.random.normal(loc=param_values, scale=stds)
                    if lower_bounds is not None and upper_bounds is not None:
                        candidate = np.clip(candidate, lower_bounds, upper_bounds)
                    pop_array[i, :] = candidate
                logger.info(f"Index {i}: {param_names} = {candidate}")

                return pop_array
            pop_array = get_pop_array_around_selected_individual(std_fraction) if bool_local_exploration_around_selected_individual else None
            logger.info(f"Initialized local exploration population with {popsize} samples, std_fraction={std_fraction}.")
            logger.info(f"pop_array: {pop_array}")
            # raise

            if pop_array is not None:
                # 初始化种群（此时所有个体的fitness都是[0,0,0]）
                # 禁止在初始化pop时运行有限元
                ad.flag_do_not_evaluate_when_init_pop = True
                pop = pg.population(prob, size=popsize)

                ad.flag_do_not_evaluate_when_init_pop = False
                for i in range(popsize):
                    pop.set_x(i, pop_array[i])  # evaluate this guy
            else:

                logger.info('Nothing exists in the archival json file. This is a whole new run.')
                ad.flag_do_not_evaluate_when_init_pop = False
                pop = pg.population(prob, size=popsize)
                ad.counter_fitness_called = ad.counter_fitness_return = 0

        # 确保这个标志在继续之前是 False
        ad.flag_do_not_evaluate_when_init_pop = False

        logger.info(f'Pop is initialized:\n {pop}')

        # 如果初始化后评估次数大于popsize，说明有新的评估，需要写入survivors
        if pop.problem.get_fevals() > popsize:
            logger.info('Write survivors.')
            ad.write_swarm_survivor(pop, ad.counter_fitness_return)

        ################################################################
        # MOO Step 2: 选择算法
        ################################################################
        # [4.3.5] 选择算法
        algo = pg.algorithm(pg.moead(gen=1, weight_generation="grid", decomposition="tchebycheff",
                                    neighbours=int(popsize / 4),
                                    CR=1, F=0.5, eta_m=20,
                                    realb=0.9,
                                    limit=2, preserve_diversity=True))
        logger.info(f'{algo}')
        logger.info(f'\t MOEA/D neighbourhood size is set to 1/4 of the popsize as {int(popsize / 4)}')

        ################################################################
        # MOO Step 3: 开始优化迭代
        ################################################################
        # [4.3.6] 开始优化
        for iteration in range(number_of_finished_iterations, 500):
            msg = '[acmop.py] This is iteration #%d. ' % iteration
            logger.info(msg)
            ad.generation = iteration
            pop = algo.evolve(pop)

            msg += 'Write survivors to file. '
            ad.write_swarm_survivor(pop, ad.counter_fitness_return)

            # 计算超体积指标
            hv = pg.hypervolume(pop)
            quality_measure = hv.compute(ref_point=self.get_bad_fintess_values(machine_type='PMSM', ref=True))
            msg += 'Quality measure by hyper-volume: %g' % quality_measure
            logger.info(msg)

            # 打印当前种群信息（如果 utility_moo 模块可用）
            try:
                import utility_moo
                utility_moo.my_print(ad, pop, iteration)
            except Exception as e:
                raise e
        logger.info('Optimization completed.')

    def sensitivity_analysis(self, study_name, parameter_dict, parameter_percentage_value_list):

        import logging
        import builtins

        builtins.ad = self  # share global variable between modules
        ad = self

        # 设置路径和日志
        self.logger = self.myLogger(self.path2SwarmData + '/', prefix=study_name)
        logger = logging.getLogger(__name__)

        # 为敏感性分析设置单独的临时文件夹
        self.get_path2SwarmData(study_name)

        # 获取当前设计的参数值
        x_denorm_dict = self.get_free_variables_as_dict()
        import copy
        print(f'当前设计的参数值 x_denorm_dict: {x_denorm_dict}')

        # 验证 parameter_dict 中的参数名是否存在于当前设计中
        for parameter_name in parameter_dict.keys():
            if parameter_name not in x_denorm_dict:
                raise ValueError(f"参数 '{parameter_name}' 不在当前设计的自由变量中。可用的参数: {list(x_denorm_dict.keys())}")

        sensitivity_dict = {}
        for parameter_name in parameter_dict.keys():
            sensitivity_dict[parameter_name] = []
            # 获取当前设计中该参数的原始值
            current_parameter_value = x_denorm_dict[parameter_name]
            for parameter_percentage_value in parameter_percentage_value_list:
                x_denorm_dict_copy = copy.deepcopy(x_denorm_dict)
                # 基于当前设计的参数值进行百分比变化
                x_denorm_dict_copy[parameter_name] = current_parameter_value * (1 + parameter_percentage_value)
                sensitivity_dict[parameter_name].append(x_denorm_dict_copy)

        import pprint
        pp = pprint.PrettyPrinter(indent=2, width=120, compact=False, sort_dicts=False)
        print('sensitivity_dict:')
        pp.pprint(sensitivity_dict)

        self.fea_config_dict["designer.Show"] = True



        # 对每个参数的每个百分比变化分别进行评估
        counter = 0
        for parameter_name in sensitivity_dict.keys():
            for percentage_idx, x_denorm_dict_variant in enumerate(sensitivity_dict[parameter_name]):
                x_denorm = list(x_denorm_dict_variant.values())

                # 为每个敏感性分析案例设置唯一的名称
                self.name = f'param-{parameter_name}-pct-{parameter_percentage_value_list[percentage_idx]:.2f}'
                
                try:
                    cost_function, f1, f2, f3, FRW, \
                    normalized_torque_ripple, \
                    normalized_force_error_magnitude, \
                    force_error_angle = self.evaluate_design_json_wrapper(x_denorm, counter)

                    print(f'敏感性分析结果 [{counter}]: 参数={parameter_name}, 百分比变化={parameter_percentage_value_list[percentage_idx]:.2f}, '
                        f'参数值={x_denorm_dict_variant[parameter_name]:.2f}, '
                        f'f1={f1:.2f}, f2={f2:.2f}, f3={f3:.2f}, FRW={FRW:.2f}, '
                        f'normalized_torque_ripple={normalized_torque_ripple:.2f}, '
                        f'normalized_force_error_magnitude={normalized_force_error_magnitude:.2f}, '
                        f'force_error_angle={force_error_angle:.2f}')
                except Exception as e:
                    print(e)
                    print('大概率是参数值超出范围了导致画不出形状报错，继续敏感性分析')

                counter += 1

if __name__ == "__main__":
    # 创建对象并导出为 JSON
    specs = MotorSpecs()
    mmd = Modern_Machine_Designer(specs)
    full_json_path = os.path.join(mmd.path2SwarmData, 'machine_designer_full.json')     # 保存完整信息到文件（类似 pickle）
    mmd.save_to_file_full(full_json_path)

    if False:  # Add previous design parameters for evaluation
        # ====== Inserted: Apply previous design parameters for evaluation ======
        prev_params = { # p4ps5 prototype from PEMD 2020 paper
            # "stator_outer_radius": 123.49969,
            # "mechanical_air_gap_length": 0.75,
            "rotor_sleeve_depth": 5.89091,               # free variable
            "magnet_depth": 5.19948,                     # free variable
            "magnet_pole_span_angle": 44.9638,           # free variable
            "stator_tooth_width": 16.099,                # free variable
            "stator_yoke_depth": 32.7594,                # free variable
            "stator_tooth_shoe_depth": 1.50079,          # free variable
            "stator_tooth_span_angle": 11.1183,          # free variable
            # "stator_tooth_depth": 42.9701,
            # "stator_tooth_open_depth": 1.50079,
            # "stator_tooth_open_angle": 5.55915,
            # "stator_tooth_tip_depth": 2.251185,
            # "stator_yoke_depth": 32.75940500000001,

            "split_ratio_r_si_slash_r_so": 0.36857582395550953,  # free variable
            # "stator_inner_radius": 45.519,
            # "outer_rotor_radius": 38.87809,
            # "inner_rotor_radius": 29.9996,
            # "rotor_iron_back_iron_depth": 3.67901,
        }

        # 使用统一接口设置参数并处理覆盖与导出参数刷新
        mmd.apply_parameter_dict(prev_params)

        mmd.show_geometry()

        # mmd.start_optimization(bool_local_exploration_around_selected_individual=True, std_fraction=0.20)

    else:
        mmd.show_geometry()

        # mmd.FEA_evaluate(counter=2026)
        # mmd.start_optimization()

        pass 
    # mmd.remove_jfiles_folders(mmd.path2SwarmData)

