import math

class ReactDrawer:
    def __init__(self, verbose_drawing=False):
        self.verbose_drawing = verbose_drawing
        self.regions = []

    def hex_to_rgb(self, hex_color):
        if not hex_color or not isinstance(hex_color, str):
            return (0.0, 0.0, 0.0)
        hex_color = hex_color.lstrip('#')
        lv = len(hex_color)
        try:
            if lv == 3:
                rgb = tuple(int(hex_color[i:i+1]*2, 16) for i in range(0, 3))
            elif lv == 6:
                rgb = tuple(int(hex_color[i:i+2], 16) for i in range(0, 6, 2))
            else:
                return (0.0, 0.0, 0.0)
            return tuple(c/255.0 for c in rgb)
        except ValueError:
            return (0.0, 0.0, 0.0)

    def drawLine(self, p1, p2):
        if self.verbose_drawing:
            print(f'[ReactDrawer] drawLine({p1=}, {p2=})')
        return [{'type': 'line', 'p1': p1, 'p2': p2}]

    def drawArc(self, centerxy, startxy, endxy):
        if self.verbose_drawing:
            print(f'[ReactDrawer] drawArc({centerxy=}, {startxy=}, {endxy=})')
        
        # Calculate flags matching Cairo and JMAG logic
        v1 = [startxy[0] - centerxy[0], startxy[1] - centerxy[1]]
        v2 = [endxy[0] - centerxy[0], endxy[1] - centerxy[1]]
        
        # Cross product to determine direction
        cross_prod = v1[0]*v2[1] - v1[1]*v2[0]
        sweep_flag = 1 if cross_prod >= 0 else 0
        
        # For motor design segments, we almost always use small arcs (large_arc_flag=0)
        # Cairo's implementation uses acos, so it's always <= PI.
        large_arc_flag = 0
            
        return [{
            'type': 'arc', 
            'center': centerxy, 
            'p1': startxy, 
            'p2': endxy,
            'sweep_flag': sweep_flag,
            'large_arc_flag': large_arc_flag
        }]

    def getSketch(self, name, color):
        pass

    def prepareSection(self, region_dict, color=None, **kwargs):
        # We just store the region data to be sent to frontend
        region_dict['color'] = color
        self.regions.append(region_dict)

    def draw_machine(self, machine):
        from machine_designer_v2 import draw_instruction_parser
        self.regions = []
        for part in machine.parts:
            # draw_instruction_parser internally calls drawer.drawArc / drawLine
            region_dict = draw_instruction_parser(part, machine.all_points, self)
            self.prepareSection(region_dict, color=part.color)
        return self.regions
