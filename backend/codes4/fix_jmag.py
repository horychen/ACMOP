import re

mapping = {
    'RatedSpeed': 'rated_speed',
    'mm_stack_length_specified': 'stack_length_specified',
    'ExcitationFreqSimulated': 'excitation_frequency_simulated',
    'DCBusVoltage': 'dc_bus_voltage',
    'mmJs': 'rated_current_density',
    'WindingFill': 'fill_factor',
    'TORQUE_CURRENT_RATIO': 'torque_current_ratio',
    'SUSPENSION_CURRENT_RATIO': 'suspension_current_ratio',
    'no_series_coil_turns_N': 'series_turns',
    'DriveW_zQ': 'drive_winding_conductors_per_slot',
    'BeariW_zQ': 'bearing_winding_conductors_per_slot',
    'mm2_slot_area': 'slot_area',
    'CurrentAmp_in_the_slot': 'slot_current_amplitude',
    'CurrentAmp_per_conductor': 'conductor_current_amplitude',
    'CurrentAmp_per_phase': 'phase_current_amplitude',
    'DriveW_CurrentAmp': 'drive_winding_current',
    'BeariW_CurrentAmp': 'bearing_winding_current',
    'slot_current_utilizing_ratio_for_torque': 'torque_current_utilization_ratio',
    'InitialRotationAngle': 'initial_rotation_angle',
    'mm2_magnet_area': 'magnet_area',
    'RotorCore_Material': 'rotor_core_material',
    'StatorCore_Material': 'stator_core_material',
    'LaminationFactor': 'lamination_factor',
    'Magnet_Name': 'magnet_name',
    'Magnet_Temperature': 'magnet_temperature',
    'Magnet_StartAngle': 'magnet_start_angle',
    'DriveW_Rs': 'phase_resistance',
    'BeariW_Rs': 'phase_resistance',
    'RatedPower': 'rated_power',
    'Temperature': 'magnet_temperature',
}

file_path = 'c:/Users/lenovo/Codes/ACMOP/backend/codes4/JMAG.py'
with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

def repl(m):
    key = m.group(1)
    if key in mapping:
        return f'acm_variant.winding.{mapping[key]}'
    return m.group(0)

def repl_EX(m):
    key = m.group(1)
    if key in mapping:
        return f'acm_variant.winding.{mapping[key]}'
    return m.group(0)

content = re.sub(r"acm_variant\.winding\.EX\['([^']+)'\]", repl, content)
content = re.sub(r'acm_variant\.winding\.EX\["([^"]+)"\]', repl, content)

content = re.sub(r"EX\['([^']+)'\]", repl_EX, content)
content = re.sub(r'EX\["([^"]+)"\]', repl_EX, content)

content = re.sub(r'^\s*EX = acm_variant\.winding\.EX\s*$', '', content, flags=re.MULTILINE)

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(content)
print("JMAG fixed")
