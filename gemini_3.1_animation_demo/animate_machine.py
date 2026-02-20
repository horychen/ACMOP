import xml.etree.ElementTree as ET
import os
import sys

def animate_machine_svg(input_svg_path, output_html_path):
    # Parse the SVG
    try:
        tree = ET.parse(input_svg_path)
        root = tree.getroot()
    except Exception as e:
        print(f"Error reading SVG: {e}")
        return

    # Handle namespaces
    ns = {'svg': 'http://www.w3.org/2000/svg'}
    ET.register_namespace('', ns['svg'])
    
    # Heuristics for stator vs rotor based on stroke color
    # Typical colors in ACMOP:
    # Stator Core: rgb(40%, 40%, 40%)
    # Coils: rgb(72.156863%, 45.098039%, 20%)
    # Rotor Shaft/Core: rgb(33.333333%, 33.333333%, 33.333333%)
    # Magnets: rgb(13.333333%, 13.333333%, 73.333333%)

    stator_colors = ['rgb(40%, 40%, 40%)', 'rgb(72.156863%, 45.098039%, 20%)']
    rotor_colors = ['rgb(33.333333%, 33.333333%, 33.333333%)', 'rgb(13.333333%, 13.333333%, 73.333333%)']

    stator_paths = []
    rotor_paths = []

    for elem in root.findall('.//svg:path', ns):
        stroke = elem.get('stroke', '')
        # Remove any extra spaces from color for robust matching
        stroke_clean = stroke.replace(" ", "")
        
        # Determine part type
        if stroke_clean in [c.replace(" ", "") for c in stator_colors]:
            stator_paths.append(elem)
        elif stroke_clean in [c.replace(" ", "") for c in rotor_colors]:
            rotor_paths.append(elem)
        else:
            # Fallback based on colors or just default to stator
            # You could also use bbox coordinates to guess, but color is safer here
            stator_paths.append(elem)
            
    # Now create the new HTML content
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Animated Electric Machine</title>
    <style>
        :root {{
            --bg-color: #0f172a;
            --text-color: #f8fafc;
            --accent-1: #818cf8;
            --accent-2: #c084fc;
            --accent-3: #38bdf8;
            --accent-4: #2dd4bf;
        }}

        body {{
            margin: 0;
            padding: 0;
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
            background-color: var(--bg-color);
            color: var(--text-color);
            font-family: system-ui, -apple-system, sans-serif;
            overflow: hidden;
        }}

        .header {{
            text-align: center;
            margin-bottom: 2rem;
            z-index: 10;
        }}

        h1 {{
            font-size: 2.5rem;
            margin: 0 0 0.5rem 0;
            background: linear-gradient(90deg, var(--accent-1), var(--accent-3));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            filter: drop-shadow(0 0 10px rgba(129, 140, 248, 0.3));
        }}

        p {{
            color: #94a3b8;
            font-size: 1.1rem;
            max-width: 600px;
            margin: 0 auto;
            line-height: 1.5;
        }}

        .animation-container {{
            width: 100%;
            max-width: 600px;
            aspect-ratio: 1;
            position: relative;
        }}

        svg {{
            width: 100%;
            height: 100%;
            filter: drop-shadow(0 0 30px rgba(56, 189, 248, 0.15));
        }}

        /* Apply infinite spinning to rotor */
        @keyframes spin-cw {{
            from {{ transform: rotate(0deg); }}
            to {{ transform: rotate(360deg); }}
        }}

        /* Based on transform scale, origin is exactly at 250,250 */
        .rotor-group {{
            transform-origin: 250px 250px;
            animation: spin-cw 4s linear infinite;
        }}
        
        /* Make styling pop for dark mode */
        path {{
            stroke-width: 0.05 !important; /* Thinner strokes */
            filter: url(#glow);
        }}
        
        .pulse {{
            animation: pulse-op 2s ease-in-out infinite alternate;
        }}
        
        @keyframes pulse-op {{
            0% {{ stroke-width: 0.05; stroke-opacity: 0.8; }}
            100% {{ stroke-width: 0.15; stroke-opacity: 1; }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Animated Electric Machine</h1>
        <p>Dynamic Python-generated SVG architecture mapping static paths to an animated, hardware-accelerated rotor.</p>
    </div>

    <div class="animation-container">
        <svg viewBox="0 0 500 500" xmlns="http://www.w3.org/2000/svg">
            <defs>
                <filter id="glow" x="-20%" y="-20%" width="140%" height="140%">
                    <feGaussianBlur stdDeviation="2" result="blur"/>
                    <feMerge>
                        <feMergeNode in="blur"/>
                        <feMergeNode in="SourceGraphic"/>
                    </feMerge>
                </filter>
            </defs>
            
            <!-- Dark background instead of white -->
            <rect x="0" y="0" width="500" height="500" fill="var(--bg-color)" />
            
            <!-- Stator (Stationary) -->
            <g class="stator-group">
                {stator_svg}
            </g>
            
            <!-- Rotor (Spinning) -->
            <g class="rotor-group">
                {rotor_svg}
            </g>
        </svg>
    </div>
</body>
</html>
"""

    def process_elem(elem, is_rotor):
        # We process the strokes to modern glowing dark mode colors
        stroke = elem.get('stroke', '').replace(" ", "")
        
        # Reset original fill since we will only stroke
        elem.set('fill', 'none')
        
        if stroke == 'rgb(40%,40%,40%)':
            elem.set('stroke', '#64748b') # Stator core (slate-500)
        elif stroke == 'rgb(72.156863%,45.098039%,20%)':
            elem.set('stroke', '#f59e0b') # Copper coils (amber-500)
            elem.set('class', 'pulse')    # Make coils pulse slightly
        elif stroke == 'rgb(33.333333%,33.333333%,33.333333%)':
            elem.set('stroke', '#94a3b8') # Rotor core (slate-400)
        elif stroke == 'rgb(13.333333%,13.333333%,73.333333%)':
            elem.set('stroke', '#06b6d4') # Magnets (cyan-500)
        else:
            elem.set('stroke', '#818cf8') 
            
        return ET.tostring(elem, encoding='unicode')

    stator_str = '\n                '.join(process_elem(e, False) for e in stator_paths)
    rotor_str = '\n                '.join(process_elem(e, True) for e in rotor_paths)
    
    final_html = html_template.format(stator_svg=stator_str, rotor_svg=rotor_str)
    
    with open(output_html_path, 'w', encoding='utf-8') as f:
        f.write(final_html)
        
    print(f"Animation successfully written to {output_html_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python animate_machine.py <input_svg> <output_html>")
        sys.exit(1)
        
    animate_machine_svg(sys.argv[1], sys.argv[2])
