import re

with open('backend_v2/JMAG.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

count = 0
for i, line in enumerate(lines):
    if '.get(' in line:
        print(f"L{i+1}: {line.strip()}")
        count += 1
print(f"Total: {count}")
