import os
from machine_design_guide import Modern_Machine_Designer

if __name__ == "__main__":
    mmd = Modern_Machine_Designer()
    full_json_path = os.path.join(mmd.path2SwarmData, 'machine_designer_full.json')     # 保存完整信息到文件（类似 pickle）
    mmd.save_to_file_full(full_json_path)

    prev_params = { # p4ps5 prototype from PEMD 2020 paper
        "rotor_sleeve_depth": 5.89091,               # free variable
        "magnet_depth": 5.19948,                     # free variable
        "magnet_pole_span_angle": 44.9638,           # free variable
        "stator_tooth_width": 16.099,                # free variable
        "stator_yoke_depth": 32.7594,                # free variable
        "stator_tooth_shoe_depth": 1.50079,          # free variable
        "stator_tooth_span_angle": 11.1183,          # free variable
        "split_ratio_r_si_slash_r_so": 0.36857582395550953,  # free variable
    }
    mmd.apply_parameter_dict(prev_params)

    x_denorm_list = mmd.sensitivity_analysis(
        study_name='TIA_prototype_SA2',
        parameter_dict=prev_params,
        parameter_percentage_value_list=[-0.2, -0.1, 0.1, 0.2],
    )

    # mmd.remove_jfiles_folders(mmd.path2SwarmData)

