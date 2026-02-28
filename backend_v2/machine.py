import os
import platform
from collections import OrderedDict
from machine_geometry import MachineGeometry

class Machine:
    """
    电机的顶层组装与调度类(Machine Wrapper)。
    将解析出的 machine_dict 参数派发至各种下游子系统 (例如几何建模生成器 MachineGeometry，以及有限元交互引擎等)。
    """
    def __init__(self, user_input: OrderedDict):
        self.user_input = user_input
        self.geometry = MachineGeometry()

        
    def sync(self):
        """
        同步机器参数字典到集合建模引擎中，提取点位配置并构建所有独立的模块实例 (定子、线圈、转子、磁钢等)。
        以 OrderedDict 作为基础数据层进行管理隔离。
        """
        self.geometry.parts.clear()
        self.geometry._next_index = 0
        
        # 利用用户参数初始化所有骨架点位
        self.geometry.sync(self.user_input)

        # 核心环节：动态注入 WindingLayout 实例，获取绕组排布与系数
        from winding_layout import winding_layout_v2
        # Notice DPNV_or_SEPA must be None for normal standard motors because winding_layout_v2 throws errors otherwise
        wily = winding_layout_v2(
            DPNV_or_SEPA=self.user_input['winding'].get('DPNV_or_SEPA', None),
            Qs=self.user_input['winding']['slot_count'],
            p=self.user_input['winding']['pole_count'] // 2,
            ps=self.user_input['winding']['suspension_pole_count'] // 2,
            coil_pitch_y=self.user_input['winding'].get('coil_pitch', 1)
        )
        # Store wily's attributes directly into the dictionary rather than the object itself
        wily_dict = {k: v for k, v in dict(wily).items() if not k.startswith('_')}
        self.user_input['winding'].update(wily_dict)
        
        # Analyze winding performance to populate necessary features (currents/turns) automatically during sync
        from machine_analyzer import MotorPerformanceAnalyzer
        self.analyzer = MotorPerformanceAnalyzer(self.user_input)
        self.analysis = self.analyzer.analyze()
        
        # FIX: Populate circuit excitaiton info for JMAG right after WindingLayout
        self.user_input['winding']['drive_winding_current'] = self.analysis['current_per_wire_a'] * 1.41421356
        self.user_input['winding']['bearing_winding_current'] = self.analysis['current_per_wire_a'] * 1.41421356
        self.user_input['winding']['wires_per_slot'] = self.analysis['wires_per_slot'] # zQ
        
        # 定义定子的槽部构形种类 (例如 "closed-slot" 表示闭槽结构)
        stator_options = f"{self.user_input['geometry']['tooth_shape']}-slot"
        
        # =========================================================
        # 通过显式定义各个子部件的属性描述字典(OrderedDict)，来规避
        # 使用厚重的类(Class/Dataclass)而实现极高灵活性的结构
        # =========================================================
        
        # 定子铁心配置
        stator_core_dict = OrderedDict([
            ('type', 'statorCore'),
            ('name', 'statorCore'),
            ('options', stator_options),
            ('color', '#666666'), # 深灰色
            ('rotation_deg', 0.0)
        ])
        
        # 绕组线圈配置
        coil_dict = OrderedDict([
            ('type', 'coil'),
            ('name', 'coil'),
            ('options', 'standard'),
            ('color', '#B87333'), # 铜色
            ('rotation_deg', 0.0)
        ])

        # 转子铁心配置
        rotor_core_dict = OrderedDict([
            ('type', 'rotorCore'),
            ('name', 'rotorCore'),
            ('options', 'cylinder'),
            ('color', '#555555'), # 灰色铁芯
            ('rotation_deg', 0.0)
        ])

        # 永磁体配置
        magnet_dict = OrderedDict([
            ('type', 'magnet'),
            ('name', 'magnet'),
            ('options', ''),
            ('color', '#2222BB'), # 代表N极或通用的蓝色
            ('rotation_deg', 0.0)
        ])

        # 依次加载入子部件管理器列表，顺序由于硬编码绑定需要严格符合 JMAG 获取 ID 列表的情况: Rotor, Magnet, Stator, Coil
        self.geometry.add_part(rotor_core_dict)
        self.geometry.add_part(magnet_dict)
        self.geometry.add_part(stator_core_dict)
        self.geometry.add_part(coil_dict)

    def draw_svg(self, filename: str):
        """输出并构建整机二维截面的矢量图(SVG)，此步无需求助外部软件，全靠内置Cairo/本地绘制算法。"""
        self.geometry.show_geometry_svg(filename=filename, scale=30.0)

    def draw_machine_using_JMAG(self):
        """
        全自动化流程: 通过 COM 接口挂起并直接驱动 JMAG Designer 建立整机全真工作工程模型，省去人工介入。
        包含开启文件定义，建立项目及材质设定。
        """
        import os
        # JMAG项目配置
        expected_project_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stator_v2.jproj')
        
        # 1. 启动 JMAG 底层工具柄，并建立 35CS250 硅钢片工程
        material_dict = self.user_input['material']
        steel_name = material_dict['stator_core_steel_name']
        toolJd = self.open_jmag(expected_project_file, steel_name)

        # 2. 依次利用底层脚本将所有的形状绘制到 JMAG 中
        self.draw_jmag(toolJd)
        print(f"JMAG geometry drawn and saved to {expected_project_file}")

    def open_jmag(self, expected_project_file, Steel_name):
        """提供统一的 JMAG SDK 连接与挂载逻辑。"""
        import os
        import platform
        import JMAG
        
        self.project_name = os.path.splitext(os.path.basename(expected_project_file))[0]
        
        # JMAG 启动时的配置项设定 (从 global 传入)
        fea_config_dict = self.user_input['fea_config_dict']

        toolJd = JMAG.JMAG(fea_config_dict=fea_config_dict)
        
        # Open Project
        toolJd.open(Steel_name=Steel_name, 
                    expected_project_file_path=expected_project_file, 
                    pc_name=fea_config_dict['pc_name'], 
                    dir_parent=os.path.abspath(os.path.dirname(__file__)) + '/')
        return toolJd

    def draw_jmag(self, toolJd):
        """使用底层解析器，依照 JMAG 可接受的特征点操作规范逐一调用绘制API，配置区域与镜像状态。"""
        import numpy as np
        from machine_geometry import draw_instruction_parser
        
        color_rgb_iron = np.array([236,236,236])/255
        color_rgb_copper = np.array([184,115,51])/255

        # 深度循环所有的预定义零部件
        for part_dict in self.geometry.parts:
            print(f"Drawing {part_dict['name']} in JMAG...")
            # 利用draw_instruction_parser获得标准的解析数据(包含是否镜像，有多少个需要环切等)
            region_dict = draw_instruction_parser(part_dict, self.geometry.all_points, toolJd)
            
            toolJd.bMirror = region_dict.get('bMirror', False)
            toolJd.iRotateCopy = region_dict.get('iRotateCopy', 0)
            
            bRotateMerge = True
            color = color_rgb_copper
            if part_dict['type'] == 'statorCore' or part_dict['name'].startswith("stator"):
                color = color_rgb_iron
                
            if part_dict['type'] == 'magnet':
                bRotateMerge = False # Prevent magnets from welding together during rotate copy

            part_dict['inner_coords'] = region_dict.get('inner_coords', {})
            # 在绘图板里执行线段和网格划分绘制
            toolJd.prepareSection(region_dict, bRotateMerge=bRotateMerge, color=color)

        # 载入并生成真实模型文件
        # 因为open_jmag内部把project推到了self.project_name
        toolJd.save(getattr(self, 'project_name', 'stator_v2'), "Machine design generated by V2")

    def FEA_evaluate(self):
        """
        Perform FEA evaluation (JMAG simulation) and compile results.
        eval_config replaces dispersed kwargs/dataclass.
        """
        import os
        from time import time as clock_time
        
        user_input = self.user_input
        
        # 将评测设定合并到机器字典以便随处调用统一管理
        eval_config = user_input['eval_config']
        self.user_input['evaluation'] = eval_config

        counter = eval_config['counter']
        counter_loop = eval_config['counter_loop']
        x_denorm = eval_config['x_denorm']
        project_loc = eval_config['project_loc']

        if x_denorm is not None:
            # 兼容：如果有denorm数据更新
            pass
            
        self.sync()

        # Project and Path Setup
        prefix = eval_config['project_name_prefix']
        project_name = f"{prefix}{counter}"
        if counter_loop > 0:
            project_name += f"-redo{counter_loop}"
            
        self.user_input['evaluation']['project_name'] = project_name
            
        jmag_temp_dir = os.path.join(project_loc, eval_config['jmag_temp_dir_name'])
        jmag_screenshots_dir = os.path.join(project_loc, eval_config['jmag_screenshots_dir_name'])
        path2FEACsv = os.path.join(project_loc, eval_config['jmag_csv_dir_name'], f"{counter}")
        
        for d in [jmag_temp_dir, jmag_screenshots_dir, path2FEACsv]:
            if not os.path.exists(d): 
                os.makedirs(d)
                
        expected_project_file = os.path.join(jmag_temp_dir, f"{project_name}.jproj")
        swarm_data_json_file_path = os.path.join(project_loc, eval_config['swarm_data_file_name'])
        self.user_input['evaluation']['swarm_data_json_file_path'] = swarm_data_json_file_path

        target_dict = self.user_input['target']
        select_FEA_tool = target_dict['select_FEA_tool']

        if 'JMAG' in select_FEA_tool:
            study_name = eval_config['study_name']

            # 1. Initialize JMAG Project
            material_dict = self.user_input['material']
            steel_name = material_dict['stator_core_steel_name']
            toolJd = self.open_jmag(expected_project_file, steel_name)

            # 2. Draw Machine Parts
            self.draw_jmag(toolJd)

            # 3. Solver Setup - Placeholder for porting JMAG API calls
            app = toolJd.app
            model = app.GetModel(project_name)
            
            # NOTE: backend_v2 structure passing self to toolJd requires toolJd to be compatible.
            toolJd.pre_process_PMSM(app, model, self)
            study = toolJd.add_magnetic_transient_study(app, model, path2FEACsv, study_name, self)
            toolJd.mesh_study(self, app, model, study, output_dir=project_loc)
            toolJd.run_study(self, app, study, self.user_input['fea_config_dict'], clock_time())
            
            # 4. Monitor the output results are generated in the specified csv folder
            # import time
            # import os
            # expected_csv = os.path.join(path2FEACsv, f"{study_name}_torque.csv")
            # print(f"Monitoring output results generation in: {path2FEACsv}")
            # wait_time = 0
            # timeout_limit = eval_config['timeout_csv_results_s']
            # while not os.path.exists(expected_csv):
            #     time.sleep(2)
            #     wait_time += 2
            #     if wait_time > timeout_limit: # timeout checks
            #         raise Exception(f"Timeout waiting for JMAG CSV results at {expected_csv}.")
            # print("CSV results detected. Waiting a few seconds for file flush...")
            # time.sleep(eval_config['csv_file_flush_sleep_s']) # Give JMAG some time to finish writing all CSV files
            
            # 5. Compile & Save Results
            self.compile_results(toolJd, study_name, path2FEACsv, swarm_data_json_file_path, select_FEA_tool)
        else:
            raise Exception(f"FEA tool {select_FEA_tool} not implemented or supported.")

    def compile_results(self, toolJd, study_name, path2FEACsv, json_path, select_FEA_tool):
        """Pack results into dictionary and save to JSON."""
        import os, jsonpickle
        project_name = self.user_input['evaluation']['project_name']
        counter = self.user_input['evaluation']['counter']
        target_dict = self.user_input['target']
        fea_config_dict = self.user_input['fea_config_dict']
        
        results = toolJd.build_str_results(self, project_name, study_name, path2FEACsv, fea_config_dict, femm_solver=None)
        
        # Unpack results indices based on build_str_results contract
        (cost_function, f1, f2, f3, FRW, \
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
         ss_avg_force_magnitude, rotor_weight, torque_average) = results

        spec_performance_dict = OrderedDict([
            ('x_denorm_dict', target_dict.get('free_parameters', {}).copy() if isinstance(target_dict.get('free_parameters'), dict) else target_dict.get('free_parameters')),
            ('project_name', project_name),
            ('individual_name', individual_name),
            ('f1', float(f1)),
            ('f2', float(f2)),
            ('f3', float(f3)),
            ('TRV', float(TRV)),
            ('FRW', float(FRW)),
            ('torque_average', torque_average),
            ('ss_avg_force_magnitude', ss_avg_force_magnitude),
            ('rotor_weight', float(rotor_weight)),
            ('normalized_torque_ripple', float(normalized_torque_ripple)),
            ('normalized_force_error_magnitude', float(normalized_force_error_magnitude)),
            ('force_error_angle', float(force_error_angle)),
            ('mm2_slot_area', mm2_slot_area),
            ('Cost', float(Cost)),
            ('rated_total_loss', float(rated_total_loss)),
            ('select_FEA_tool', select_FEA_tool),
        ])

        # Save to JSON via jsonpickle
        results2file = {f'spec_performance_dict-ind{counter}': spec_performance_dict}
        
        try:
            if os.path.exists(json_path) and os.path.getsize(json_path) > 0:
                with open(json_path, 'r') as rf:
                    existing_data = jsonpickle.decode(rf.read())
            else:
                existing_data = {}
        except:
            existing_data = {}

        existing_data.update(results2file)
        with open(json_path, 'w') as wf:
            wf.write(jsonpickle.encode(existing_data, indent=4))

        target_dict['results_for_optimization'] = (cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle)
