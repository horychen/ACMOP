# main_dex13.py

from collections import OrderedDict
from machine_analyzer import WIRE_DATA
import os
import platform
from machine import Machine

# 1. 这个字典即等同于全链路的 Single Source of Truth，承载从用户前端接收的所有特征属性
user_input = OrderedDict([
        ('winding', OrderedDict([
            ('phase_count', 3),
            ('slot_count', 12),
            ('pole_count', 10),
            ('coil_pitch', 1),
            ('l_stack', 16.0), # mm
            ('rated_current_density', 14), # 设定热负荷标准与基准线 A/mm^2
            ('dc_bus_voltage', 12), # 直流母线电压,用于判断最高可输出转速限值 V
            ('fill_factor', 0.58),
            ('rated_speed', 20000), # 额定评估转速点
            # ('DPNV_or_SEPA', True), # Winding Layout option # default is DPNV, SEPA is removed in this work
            # ('separate_winding_utilization_ratio_for_torque', 1.0), # 如果是 SEPA绕组会小于1
            # ('separate_winding_utilization_ratio_for_suspension', 0.0), # 如果是 SEPA绕组小于1，这两个数字加起来是1
            ('torque_current_ratio', 1.0), # 如果不是无轴承电机就是1
            ('suspension_current_ratio', 0.0),
            ('suspension_pole_count', 8), # 悬浮绕组极数
            ('wires_per_slot', None), # zQ, e.g., 34 
            ('drive_winding_current', None), # 驱动目的电流
            ('bearing_winding_current', None), # 悬浮目的电流
            ('awg', 31), # 如果未指定使用的线规AWG，默认采用AWG 31
            ('wire_diameter_with_insulation', WIRE_DATA[31]['coated']), # AWG31 带漆外径
            ('wire_diameter', WIRE_DATA[31]['bare']), # AWG31 裸铜直径
            ('liner', 0.1), # 添加默认的槽绝缘衬套厚度（liner），如果未指定，默认给 0.1mm
            ('connection_type', 'wye'), # or 'delta'
            ('kw1', None), 
            ('phase_resistance', None),
        ])),
        ('geometry', OrderedDict([
            ('tooth_shape', 'closed'), # 牙槽结构样式设定
            ('d_tooth_shoe', 0.15),
            ('d_tooth', 2.0),
            ('w_tooth', 1.2),
            ('d_magnet', 1.0),
            ('d_air_gap', 0.15),
            ('r_stator_outer', 13/2),
            ('r_rotor_outer', 8/2),
            ('r_shaft', 0.0), # 转子轴向孔径 (0等于实心)
        ])),
        ('material', OrderedDict([
            ('stator_core_steel_name', '20JNEH1200'), # other options include China steel 35CS250
            ('rotor_core_steel_name', '20JNEH1200'),
            ('lamination_factor', 98), # in [%]
            ('magnet_material_name', 'N42SH'),
            ('magnet_start_angle', 0), # offset the magnetization to the magnet by an angle 径向磁化永磁体的时候转过去一定的角度
            ('magnet_temperature', 80),
        ])),
        ('fea_config_dict', OrderedDict([
            ('pc_name', platform.node()),
            ('JMAG_Designer_Version', '20.0'),
            ('designer.show', True),
            ('designer.max_nonlinear_iteration', 50),
            ('delete_results_after_calculation', False),
            ('designer.JMAG_Scheduler', False),
            ('designer.MultipleCPUs', False),
            ('designer.AddIronLossCondition', True),
            ('designer.OnlyTableResults', True),
            ('designer.number_cycles_in_1stTSS', 3),
            ('designer.number_cycles_in_2ndTSS', 0.5),
            ('designer.number_cycles_in_3rdTSS', 0.0),
            ('designer.number_cycles_prolonged', 0.0),
            ('designer.number_of_steps_1stTSS', 24),
            ('designer.number_of_steps_2ndTSS', 32),
            ('designer.StepPerCycle_3rdTSS', 64),
            ('designer.TranRef-StepPerCycle', 64),
            ('designer.CircumferentialDivision', 720),
            ('designer.meshSize_Stator', 0.15),
            ('designer.meshSize_Rotor', 0.25),
            ('designer.meshSize_Magnet', 0.2),
            ('designer.meshSizeAir', 0.05),
            ('designer.meshSize_General', 2), # Mesh for Coil
        ])),
        ('target', OrderedDict([
            ('initial_rotation_angle', 0.0), # 转子初始位置，需要和电流矢量垂直
            ('initial_rotation_angle_increment', 3.0), # mechanical degrees
            ('bool_multipleCases', None),
            ('machine_class', 'SPMSM'),
            ('select_FEA_tool', 'JMAG'),
            ('rated_power', 10),
            ('free_parameters', [
                'search w_stator_width within [1.0, 1.5]',
                'search d_stator_tooth within [1.6, 2.0]',
                'search d_magnet within [1.0, 3.0]',
            ])
        ])),
        ('eval_config', OrderedDict([
            ('project_loc', os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")),
            ('project_name_prefix', 'machine-ind'),
            ('counter', 0),
            ('counter_loop', 0),
            ('x_denorm', None),
            ('study_name', 'Transient'),
            ('jmag_temp_dir_name', 'jmag_temp'),
            ('jmag_csv_dir_name', 'csv'),
            ('jmag_screenshots_dir_name', 'jmag_screenshots'),
            ('swarm_data_file_name', 'SwarmData.json'),
            ('timeout_csv_results_s', 180),
            ('csv_file_flush_sleep_s', 5),
            ('moo.popsize', 78),
            ('moo.fitness_OA', 'TorqueDensity'),
            ('moo.fitness_OB', 'Efficiency'),
            ('moo.fitness_OC', 'Cost'),
        ]))
    ])

if __name__ == "__main__":
    mac = Machine(user_input)
    mac.sync()
    

    if True:
        import rich
        from rich.tree import Tree
        from rich.panel import Panel

        # Analysis is already computed during mac.sync(), get reference
        analysis = mac.analysis

        input_tree = Tree("[bold blue]Motor Parameters (user_input)[/bold blue]")
        for category, params in mac.user_input.items():
            branch = input_tree.add(f"[bold magenta]{category}[/bold magenta]")
            for key, value in params.items():
                branch.add(f"[cyan]{key}[/cyan]: [green]{value}[/green]")
        
        rich.print(Panel(input_tree, title="Configuration", expand=False))
        rich.print(analysis)

        r_stator_outer = user_input['geometry']['r_stator_outer']
        print(f"--- 12S10P {r_stator_outer*2}mm 微电机后端第一步分析结果 ---")
        print(f"AWG: {analysis['awg']}, 满槽率: {analysis['slot_fill_factor_cu']}%, 槽内总导线: {analysis['wires_per_slot']} 根")
        print(f"安匝(NI): {analysis['ampere_turns_per_slot']}, 热负荷 (AJ): {analysis['thermal_load_aj']} A^2/(cm*mm^2)")
        print(f"KV 值估算: {analysis['KV_rpm_V']} RPM/V")
        print(f"相电阻 (100°C): {analysis['phase_resistance_100C']} Ohm")
        print(f"20k RPM 总损耗: {analysis['total_loss_W']} W (铜损 {analysis['copper_loss_W']}W, 铁损 {analysis['iron_loss_W']}W)")
        print("-" * 60)
                
    # 步骤二： 输出矢量截面图，用于在不启动 FEA 商业软件情况下的独立 Web UI 展示与几何校验
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output_svg_path = os.path.join(output_dir, "stator_v2.svg")
    mac.draw_svg(output_svg_path)
    print(f"Exported {output_svg_path} successfully with user_input parameters.")
    
    # =================================================================
    # 可选环节：执行有限元建模接口 (此指令将拉起第三方仿真软件建立工程)
    # 屏蔽此段即可进行低成本的 CI/CD 或快速参数迭代
    # =================================================================
    
    # 启用FEA_evaluate流水线
    mac.FEA_evaluate()
