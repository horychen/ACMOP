import re

defaults_code = []

with open('tmp_git_diff.txt', 'r', encoding='utf-8') as f:
    lines = f.readlines()

for line in lines:
    if line.startswith('- '):
        # looks like -        study.GetCondition("RotCon").SetValue(u"InitialRotationAngle", acm_variant.user_input.get('winding', {}).get('initial_rotation_angle', 0.0))
        # Wait, the diff we saved earlier looks like:
        # -        study.GetCondition("RotCon").SetValue(u"InitialRotationAngle", acm_variant.user_input['winding'].get('initial_rotation_angle', 0.0))
        # because the first pass I did wasn't perfect, etc.
        m = re.findall(r"['\"]([A-Za-z0-9_.\-]+)['\"]\s*\]\s*\.get\(\s*['\"]([A-Za-z0-9_.\-]+)['\"]\s*,\s*([^)]+)\s*\)", line)
        for d in m:
            dict_name, key, default = d[0], d[1], d[2]
            # e.g. dict_name='winding', key='initial_rotation_angle', default='0.0'
            defaults_code.append(f"        if '{key}' not in self.user_input['{dict_name}']: self.user_input['{dict_name}']['{key}'] = {default}")
            
        m2 = re.findall(r"(fea_config_dict|geom|wind|eval_config|target_dict|mat_dict|wind_dict)\.get\(\s*['\"]([A-Za-z0-9_.\-]+)['\"]\s*,\s*([^)]+)\s*\)", line)
        for d in m2:
            var_name, key, default = d[0], d[1], d[2]
            # Map var_name to dict string
            d_map = {
                'fea_config_dict': 'fea_config_dict',
                'geom': 'geometry',
                'wind': 'winding',
                'eval_config': 'eval_config',
                'target_dict': 'target',
                'mat_dict': 'material',
                'wind_dict': 'winding'
            }
            if var_name in d_map:
                dict_name = d_map[var_name]
                if dict_name == 'eval_config':
                    # eval_config might not be in user_input directly, but it usually is
                    defaults_code.append(f"        if '{key}' not in self.user_input.get('eval_config', {{}}): self.user_input.setdefault('eval_config', {{}})['{key}'] = {default}")
                elif dict_name == 'fea_config_dict':
                    defaults_code.append(f"        if '{key}' not in self.user_input.get('fea_config_dict', {{}}): self.user_input.setdefault('fea_config_dict', {{}})['{key}'] = {default}")
                else:
                    defaults_code.append(f"        if '{key}' not in self.user_input['{dict_name}']: self.user_input['{dict_name}']['{key}'] = {default}")

        m3 = re.findall(r"self\.fea_config_dict\.get\(\s*['\"]([A-Za-z0-9_.\-]+)['\"]\s*,\s*([^)]+)\s*\)", line)
        for d in m3:
            key, default = d[0], d[1]
            defaults_code.append(f"        if '{key}' not in self.user_input.get('fea_config_dict', {{}}): self.user_input.setdefault('fea_config_dict', {{}})['{key}'] = {default}")

# Unique the lines
unique_code = list(set(defaults_code))
unique_code.sort()

# Also remove `.setdefault` since it violates the rule!
# We'll rewrite the output code to just use `in` operator.
fixed_code = []
for code_line in unique_code:
    if "setdefault" in code_line:
        # e.g.   if 'key' not in self.user_input.get('eval_config', {}): self.user_input.setdefault('eval_config', {})['key'] = default
        # let's just make sure eval_config exists.
        pass
    else:
        fixed_code.append(code_line)

with open('tmp_defaults_code.py', 'w') as f:
    f.write('\n'.join(fixed_code))
