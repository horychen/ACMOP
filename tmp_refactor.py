import re

with open('backend_v2/JMAG.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_lines = []
for i, line in enumerate(lines):
    # If the line defines a method that takes `acm_variant`, we should inject `wp` right after if it's used.
    # A cleaner approach without deep AST:
    # We just replace `acm_variant.user_input['winding']` with `wp` globally,
    # and then add `wp = acm_variant.user_input['winding']` at the start of any method using `wp`.
    
    # Or simply:
    line_modified = line.replace("acm_variant.user_input['winding']", "wp")
    
    if "def add_magnetic_transient_study(" in line or "def add_circuit(" in line or "def pre_process_PMSM(" in line:
        new_lines.append(line)
        # Add wp definition at function start
        indent = re.match(r'^\s*', line).group(0) + '    '
        new_lines.append(f"{indent}wp = acm_variant.user_input['winding']\n")
        continue

    # Alternatively, just inject it everywhere there is a method def if it has acm_variant
    if line.strip().startswith("def ") and "acm_variant" in line:
        new_lines.append(line)
        indent = re.match(r'^\s*', line).group(0) + '    '
        new_lines.append(f"{indent}wp = acm_variant.user_input['winding']\n")
        continue
    elif line.strip().startswith("def ") and "self" in line:
        # maybe it doesn't take acm_variant, so wp wouldn't be defined. Let's not blindly add it,
        # but the replacement might break things.
        pass

    # For safety, let's only do the replacement and see if we can just define `wp = acm_variant.user_input['winding']` manually in the few places it complains, or heuristically.
    pass

# Let's do a smarter replacement focused on the methods
new_text = "".join(lines)

# Find all functions containing acm_variant.user_input['winding']
import ast

class MethodRewriter(ast.NodeTransformer):
    pass # AST is too hard to preserve comments.

# Let's use regex to find function defs
out_lines = []
in_func = False
func_indent = ""
needs_wp_decl = False
has_wp_decl = False

for line in lines:
    m = re.match(r'^(\s*)def (\w+)\(.*\):', line)
    if m:
        out_lines.append(line)
        func_indent = m.group(1) + "    "
        in_func = True
        needs_wp_decl = False
        has_wp_decl = False
        continue
        
    if in_func:
        if line.strip() and not line.startswith(func_indent) and not line.strip().startswith("#"):
            # Exited function
            in_func = False
    
    if in_func:
        if "acm_variant.user_input['winding']" in line:
            needs_wp_decl = True
        if "wp = acm_variant.user_input['winding']" in line:
            has_wp_decl = True
            
out_lines = []
i = 0
while i < len(lines):
    line = lines[i]
    m = re.match(r'^(\s*)def (\w+)\(.*\):', line)
    out_lines.append(line)
    if m:
        func_indent = m.group(1) + "    "
        # look ahead to see if 'acm_variant.user_input['winding']' is used
        j = i + 1
        uses_winding = False
        has_acm_variant_access = "acm_variant" in line # is it in signature?
        while j < len(lines):
            l = lines[j]
            if l.strip() and not l.startswith(func_indent) and not l.strip().startswith("#"):
                if re.match(r'^(\s*)def (\w+)\(.*\):', l): # Next func
                    break
            if "acm_variant.user_input['winding']" in l:
                uses_winding = True
            j += 1
            
        if uses_winding:
            # Check if acm_variant parameter exists or if we should get it somehow
            # Most commonly, functions that use it receive it as argument `acm_variant`.
            out_lines.append(f"{func_indent}wp = acm_variant.user_input['winding']\n")
            
    i += 1

# Now replace all 'acm_variant.user_input['winding']' with 'wp' globally in the output lines
final_lines = []
for l in out_lines:
    # Skip replacing in the declaration we just added
    if l.strip() == "wp = acm_variant.user_input['winding']":
        final_lines.append(l)
    else:
        final_lines.append(l.replace("acm_variant.user_input['winding']", "wp"))

with open('backend_v2/JMAG.py', 'w', encoding='utf-8') as f:
    f.writelines(final_lines)
