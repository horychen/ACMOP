import re

def replace_get(match):
    # match.group(0) is the full match, e.g. .get('abc', 10)
    # match.group(1) is the key, e.g. 'abc'
    # match.group(2) is the default, e.g. 10 or None
    key = match.group(1)
    return f"[{key}]"

with open('backend_v2/JMAG.py', 'r', encoding='utf-8') as f:
    text = f.read()

# Pattern to match .get('key') or .get("key") or .get('key', default)
# It handles nested parentheses in default by NOT matching it if it's too complex, but for simple ones it works.
# Actually, since default could be a function call, a simple regex might fail on nested parens.
# Let's use a slightly more robust regex or an ast-based approach.

import ast

class GetReplacer(ast.NodeTransformer):
    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'get':
            if len(node.args) >= 1 and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
                # Replace dict.get('key', ...) with dict['key']
                new_node = ast.Subscript(
                    value=node.func.value,
                    slice=ast.Index(value=node.args[0]) if hasattr(ast, 'Index') else node.args[0],
                    ctx=ast.Load()
                )
                return ast.copy_location(new_node, node)
        return node

# AST unparsing removes comments, so we should try regex.

def regex_replace_get():
    new_text = text
    # Let's repeatedly replace .get('key', default) until no more matches
    # This regex matches .get( "string", ... )
    # We only care about string keys.
    pattern = re.compile(r'\.get\s*\(\s*([\'"][a-zA-Z0-9_ .\-]+[\'"])\s*(?:,[^)]*)?\)')
    
    # We will iterate lines to be safe against multi-line
    new_lines = []
    for line in new_text.split('\n'):
        # simple replacement for `\.get\('key', default\)`
        # this won't handle nested parens in the default argument if there is a closing paren inside it.
        # But looking at JMAG.py, defaults are usually simple: {}, 10, None, 0.0, True, False, etc.
        # Let's do a smart regex: non-greedy up to the next matching parenthesis.
        
        # We can implement a simple bracket matcher
        line_out = line
        while True:
            m = re.search(r'\.get\s*\(\s*([\'"][a-zA-Z0-9_ .\-]+[\'"])', line_out)
            if not m:
                break
            start = m.start()
            # find matching closing parenthesis
            parens = 0
            idx = start + 4
            while idx < len(line_out):
                if line_out[idx] == '(': parens += 1
                elif line_out[idx] == ')':
                    parens -= 1
                    if parens == 0:
                        break
                idx += 1
            if parens == 0: # found matching paren
                end = idx
                key = m.group(1)
                line_out = line_out[:start] + f"[{key}]" + line_out[end+1:]
            else:
                break # couldn't find, skip
        new_lines.append(line_out)
    
    return '\n'.join(new_lines)

final_text = regex_replace_get()
with open('backend_v2/JMAG.py', 'w', encoding='utf-8') as f:
    f.write(final_text)

print("Replacement complete.")
