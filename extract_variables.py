import ast
import os
import sys
from collections import defaultdict

def parse_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        try:
            tree = ast.parse(f.read(), filename=filepath)
        except Exception as e:
            return None, str(e)

    class Analyzer(ast.NodeVisitor):
        def __init__(self):
            self.functions = []
            self.current_class = None

        def visit_ClassDef(self, node):
            prev_class = self.current_class
            self.current_class = node.name
            self.generic_visit(node)
            self.current_class = prev_class

        def visit_FunctionDef(self, node):
            func_name = node.name
            if self.current_class:
                func_name = f"{self.current_class}.{func_name}"
            args = [arg.arg for arg in node.args.args]
            
            reads = set()
            writes = set()
            
            class BodyAnalyzer(ast.NodeVisitor):
                def get_full_name(self, node):
                    if isinstance(node, ast.Name):
                        return node.id
                    elif isinstance(node, ast.Attribute):
                        val = self.get_full_name(node.value)
                        if val:
                            return f"{val}.{node.attr}"
                    elif isinstance(node, ast.Subscript):
                        val = self.get_full_name(node.value)
                        if isinstance(node.slice, ast.Constant):
                            slice_val = node.slice.value
                            if val:
                                return f"{val}['{slice_val}']"
                        elif isinstance(node.slice, ast.Index) and hasattr(node.slice, 'value'): # Python < 3.9
                            if isinstance(node.slice.value, ast.Str):
                                if val:
                                    return f"{val}['{node.slice.value.s}']"
                    elif isinstance(node, ast.Call):
                        if isinstance(node.func, ast.Attribute) and node.func.attr == 'get':
                            val = self.get_full_name(node.func.value)
                            if val and len(node.args) >= 1 and isinstance(node.args[0], ast.Constant):
                                return f"{val}['{node.args[0].value}']"
                    return None

                def visit_Call(self, node):
                    name = self.get_full_name(node)
                    if name:
                        reads.add(name)
                    self.generic_visit(node)

                def visit_Attribute(self, node):
                    name = self.get_full_name(node)
                    if name:
                        reads.add(name)
                    self.generic_visit(node)
                    
                def visit_Subscript(self, node):
                    name = self.get_full_name(node)
                    if name:
                        reads.add(name)
                    self.generic_visit(node)

                def visit_Name(self, node):
                    if isinstance(node.ctx, ast.Store):
                        writes.add(node.id)
                    elif isinstance(node.ctx, ast.Load):
                        reads.add(node.id)
                    self.generic_visit(node)
                
                def visit_Assign(self, node):
                    for target in node.targets:
                        name = self.get_full_name(target)
                        if name:
                            writes.add(name)
                    self.generic_visit(node)

            analyzer = BodyAnalyzer()
            for stmt in node.body:
                analyzer.visit(stmt)
                
            self.functions.append({
                'name': func_name,
                'args': args,
                'reads': sorted(list(reads)),
                'writes': sorted(list(writes)),
                'lineno': node.lineno
            })
            
    visitor = Analyzer()
    visitor.visit(tree)
    return visitor.functions, None

def analyze_directory(directory):
    # Dictionaries to track which variable is read/written by which functions
    var_read_by = defaultdict(list)
    var_written_by = defaultdict(list)
    
    file_markdown = ""
    
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                filepath = os.path.join(root, file)
                rel_path = os.path.relpath(filepath, directory).replace('\\', '/')
                funcs, err = parse_file(filepath)
                if err or not funcs:
                    continue
                
                file_markdown += f"\n## File: `{rel_path}`\n\n"
                for func in funcs:
                    func_id = f"{rel_path}::{func['name']}"
                    file_markdown += f"### `{func['name']}` (Line {func['lineno']})\n"
                    file_markdown += f"- **Arguments**: {', '.join(func['args']) if func['args'] else 'None'}\n"
                    
                    # Heuristic filters for global/state variables
                    def filter_var(x):
                        if not x: return False
                        if 'self.' in x or 'EX[' in x or 'machine.' in x or 'stator.' in x or 'rotor.' in x or 'app' in x or 'wily' in x:
                            return True
                        return False
                    
                    global_reads = [r for r in func['reads'] if filter_var(r)]
                    global_writes = [w for w in func['writes'] if filter_var(w)]
                    
                    for r in global_reads: var_read_by[r].append(func_id)
                    for w in global_writes: var_written_by[w].append(func_id)
                    
                    if global_reads:
                        file_markdown += f"- **State Inputs (Reads)**:\n"
                        for r in sorted(list(set(global_reads))):
                            file_markdown += f"  - `{r}`\n"
                    else:
                        file_markdown += f"- **State Inputs**: None detected\n"
                        
                    if global_writes:
                        file_markdown += f"- **State Outputs (Writes)**:\n"
                        for w in sorted(list(set(global_writes))):
                            file_markdown += f"  - `{w}`\n"
                    else:
                        file_markdown += f"- **State Outputs**: None detected\n"
                    
                    file_markdown += "\n"

    # Now generate the category section
    markdown = "# Backend Dependency Report\n\n"
    markdown += "This report serves to help with refactoring the codebase away from the global object pattern. It lists the state variables and the functions that depend on them, followed by an exhaustive list of functions and their implicit inputs/outputs.\n\n"
    
    markdown += "## Variable Directory (Categorized by where they are needed)\n\n"
    markdown += "This section lists global/state variables and the functions that read or write to them, making it easy to see the dependency graph of each variable.\n\n"
    
    all_vars = sorted(list(set(var_read_by.keys()).union(set(var_written_by.keys()))))
    
    for v in all_vars:
        markdown += f"### Variables: `{v}`\n"
        if v in var_read_by:
            markdown += "- **Read by (Depends on this)**:\n"
            for f in sorted(list(set(var_read_by[v]))):
                markdown += f"  - `{f}`\n"
        if v in var_written_by:
            markdown += "- **Written by (Modifies this)**:\n"
            for f in sorted(list(set(var_written_by[v]))):
                markdown += f"  - `{f}`\n"
        markdown += "\n"
        
    markdown += "## Function Input/Output List\n\n"
    markdown += file_markdown
    
    return markdown

if __name__ == '__main__':
    target_dir = sys.argv[1]
    report = analyze_directory(target_dir)
    with open('backend_dependency_report.md', 'w', encoding='utf-8') as f:
        f.write(report)
    print("Report generated at backend_dependency_report.md")
