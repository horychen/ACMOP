import math
import cairo
import os

# Global settings to match PyX_Utility
global_settings = {
    'wscale': 1,
    'xscale': 1.5,
    'linewidth': 0.035, # Default thin line
    'linecolor': (0, 0, 0), # Black
}

# Mapping dict for colors and styles (simplified for Cairo)
mapping_dict = {
    'dashed': [4.0, 4.0], # Simplified dash pattern
    'dense-dashed': [1.0, 1.0],
    'my-thick-line': 0.2,
    'red': (1, 0, 0),
    'blue': (0, 0, 1),
    'RawSienna': (0.78, 0.57, 0.22), # Approx
    'darkgreen': (0, 0.39, 0),
    'tint-red': (1, 0.74, 0.68),
    'tint-blue': (0.7, 0.96, 1),
    'tint-yellow': (1, 0.94, 0.7),
    'warm-grey': (0.98, 0.98, 0.97),
    'night-blue': (0.14, 0.22, 0.36),
    'purple-blue': (0.4, 0.4, 0.6),
}

# Mocking pyx.trafo
class Trafo:
    @staticmethod
    def translate(x, y):
        return ('translate', x, y)
    
    @staticmethod
    def scale(s):
        return ('scale', s)

trafo = Trafo()

class PyX_Utility:
    def __init__(self):
        self.commands = []
        self.cvs = self.CanvasWrapper(self)
        self.scale_factor = 28.35 # 1 cm = 28.35 pts approx
        
        # Bounding box tracking
        self.min_x = float('inf')
        self.min_y = float('inf')
        self.max_x = float('-inf')
        self.max_y = float('-inf')

    def _update_bounds(self, x, y):
        if x < self.min_x: self.min_x = x
        if x > self.max_x: self.max_x = x
        if y < self.min_y: self.min_y = y
        if y > self.max_y: self.max_y = y

    class CanvasWrapper:
        def __init__(self, parent):
            self.parent = parent

        def writePDFfile(self, filename):
            self.parent.render(filename, 'pdf')

        def writeSVGfile(self, filename):
            self.parent.render(filename, 'svg')
            
        def insert(self, other_cvs_wrapper, attrs=[]):
            # Support for collage (inserting another canvas)
            dx, dy = 0, 0
            scale = 1.0
            
            for attr in attrs:
                if isinstance(attr, tuple):
                    if attr[0] == 'translate':
                        dx += attr[1]
                        dy += attr[2]
                    elif attr[0] == 'scale':
                        scale *= attr[1]
            
            # Update bounds based on the inserted canvas's bounds and transformation
            other = other_cvs_wrapper.parent
            if other.min_x != float('inf'):
                # Transform corners of the other bounding box
                corners = [
                    (other.min_x, other.min_y),
                    (other.max_x, other.min_y),
                    (other.min_x, other.max_y),
                    (other.max_x, other.max_y)
                ]
                for cx, cy in corners:
                    tx = cx * scale + dx
                    ty = cy * scale + dy
                    self.parent._update_bounds(tx, ty)

            # Capture commands
            other_commands = other.commands
            
            def draw_inserted(ctx):
                ctx.save()
                ctx.translate(dx, dy)
                ctx.scale(scale, scale)
                
                for cmd in other_commands:
                    ctx.save()
                    cmd(ctx)
                    ctx.restore()
                ctx.restore()
                
            self.parent.commands.append(draw_inserted)

    def _parse_settings(self, settings):
        style = {
            'linewidth': global_settings['linewidth'],
            'linecolor': global_settings['linecolor'],
            'dash': [],
            'rgb': None
        }
        
        for s in settings:
            if s in mapping_dict:
                val = mapping_dict[s]
                if isinstance(val, tuple): # Color
                    style['linecolor'] = val
                    style['rgb'] = val
                elif isinstance(val, list): # Dash
                    style['dash'] = val
                elif isinstance(val, float) or isinstance(val, int): # Linewidth
                    style['linewidth'] = val
            elif isinstance(s, list): 
                 pass 
        return style

    def pyx_line(self, p1, p2, bool_track=False, settings=[], **kwarg):
        style = self._parse_settings(settings)
        self._update_bounds(p1[0], p1[1])
        self._update_bounds(p2[0], p2[1])
        
        def draw(ctx):
            ctx.set_source_rgb(*style['linecolor'])
            ctx.set_line_width(style['linewidth'])
            ctx.set_dash(style['dash'])
            ctx.move_to(p1[0], p1[1])
            ctx.line_to(p2[0], p2[1])
            ctx.stroke()
            
        self.commands.append(draw)
        if bool_track:
            return [*p1, *p2]

    def pyx_arc(self, startxy, endxy, centerxy=(0,0), bool_track=False, settings=[], **kwarg):
        style = self._parse_settings(settings)
        self._update_bounds(startxy[0], startxy[1])
        self._update_bounds(endxy[0], endxy[1])
        # Curve control points might go outside, but start/end is a good approximation for now
        
        def draw(ctx):
            ctx.set_source_rgb(*style['linecolor'])
            ctx.set_line_width(style['linewidth'])
            ctx.set_dash(style['dash'])
            
            ctx.move_to(startxy[0], startxy[1])
            mx, my = (startxy[0] + endxy[0])/2, (startxy[1] + endxy[1])/2
            dx, dy = endxy[0] - startxy[0], endxy[1] - startxy[1]
            length = math.sqrt(dx*dx + dy*dy)
            if length > 0:
                offset = length * 0.2 
                cx, cy = mx - dy/length * offset, my + dx/length * offset
                ctx.curve_to(cx, cy, cx, cy, endxy[0], endxy[1])
            else:
                ctx.line_to(endxy[0], endxy[1])
                
            ctx.stroke()

        self.commands.append(draw)

    def pyx_text(self, loc, text, size=5, scale=1.0, BoxColor=None, settings=[]):
        clean_text = text.replace(r'\textbf{', '').replace('}', '').replace('$', '')
        
        # Estimate bounds (rough)
        # Assuming size is roughly height in user units
        font_size = size * 0.35 * scale
        est_width = len(clean_text) * font_size * 0.6
        est_height = font_size
        
        # Centered
        self._update_bounds(loc[0] - est_width/2, loc[1] - est_height/2)
        self._update_bounds(loc[0] + est_width/2, loc[1] + est_height/2)

        def draw(ctx):
            ctx.set_source_rgb(0, 0, 0) 
            # ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.select_font_face("Times New Roman", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            # Re-calculate font size inside draw context to be sure
            ctx.set_font_size(font_size)
            
            x_bearing, y_bearing, width, height, x_advance, y_advance = ctx.text_extents(clean_text)
            
            x = loc[0] - width / 2 - x_bearing
            y = loc[1] - height / 2 - y_bearing 
            
            if BoxColor and BoxColor in mapping_dict:
                box_rgb = mapping_dict[BoxColor]
                ctx.set_source_rgb(*box_rgb)
                ctx.rectangle(x + x_bearing - 2, y + y_bearing - 2, width + 4, height + 4)
                ctx.fill()
                ctx.set_source_rgb(0, 0, 0)
                
            ctx.move_to(x, y)
            
            matrix = ctx.get_matrix()
            ctx.save()
            ctx.translate(x, y)
            ctx.scale(1, -1) 
            ctx.move_to(0, 0)
            ctx.show_text(clean_text)
            ctx.restore()
            
        self.commands.append(draw)
        return None

    def pyx_arrow(self, PA, PB=None, settings=[]):
        if PB is None:
            PB = PA
            PA = [0,0]
            
        style = self._parse_settings(settings)
        self._update_bounds(PA[0], PA[1])
        self._update_bounds(PB[0], PB[1])
        
        def draw(ctx):
            ctx.set_source_rgb(*style['linecolor'])
            ctx.set_line_width(style['linewidth'])
            ctx.set_dash(style['dash'])
            
            ctx.move_to(PA[0], PA[1])
            ctx.line_to(PB[0], PB[1])
            ctx.stroke()
            
            angle = math.atan2(PB[1] - PA[1], PB[0] - PA[0])
            arrow_len = 0.5 
            arrow_angle = 0.5 
            
            ctx.move_to(PB[0], PB[1])
            ctx.line_to(PB[0] - arrow_len * math.cos(angle - arrow_angle), PB[1] - arrow_len * math.sin(angle - arrow_angle))
            ctx.line_to(PB[0] - arrow_len * math.cos(angle + arrow_angle), PB[1] - arrow_len * math.sin(angle + arrow_angle))
            ctx.close_path()
            ctx.fill()
            
        self.commands.append(draw)

    def pyx_arrow_both_ends(self, PA, PB=None, settings=[]):
        if PB is None:
            PB = PA
            PA = [0,0]
        self.pyx_arrow(PA, PB, settings)
        self.pyx_arrow(PB, PA, settings)

    def pyx_horizental_arrow_both_ends_with_text_inside(self, text, PA, PB=None, BoxColor=None, size=5, scale=1.0, settings=[]):
        if PB is None:
            PB = PA
            PA = [0,0]
        
        mid = [(PA[0]+PB[0])/2, (PA[1]+PB[1])/2]
        self.pyx_text(mid, text, size=size, scale=scale, BoxColor=BoxColor, settings=settings)
        self.pyx_arrow_both_ends(PA, PB, settings)

    def pyx_draw_sector(self, origin, angle_begin, angle_end, radius, bool_stroke=False):
        # Update bounds: origin and points on arc
        self._update_bounds(origin[0], origin[1])
        # Check 4 quadrants + start/end
        angles = [angle_begin, angle_end]
        # Add cardinal directions if they are within the range
        for a in range(0, 361, 90):
            # Normalize angles to handle wrap around if needed, but for now simple check
            if angle_begin <= a <= angle_end:
                angles.append(a)
        
        for a in angles:
            rad = a * math.pi / 180
            self._update_bounds(origin[0] + radius * math.cos(rad), origin[1] + radius * math.sin(rad))

        def draw(ctx):
            ctx.set_source_rgb(0.9, 0.9, 0.9) 
            ctx.move_to(origin[0], origin[1])
            
            a1 = angle_begin * math.pi / 180
            a2 = angle_end * math.pi / 180
            
            ctx.arc(origin[0], origin[1], radius, a1, a2)
            ctx.close_path()
            ctx.fill()
            
            if bool_stroke:
                ctx.set_source_rgb(0, 0, 0)
                ctx.set_line_width(global_settings['linewidth'])
                ctx.move_to(origin[0], origin[1])
                ctx.arc(origin[0], origin[1], radius, a1, a2)
                ctx.close_path()
                ctx.stroke()
                
        self.commands.append(draw)

    def pyx_circle(self, radius, center=[0,0], bool_dashed=False, dash_list=[0,4], linewidth=0.035, settings=[]):
        style = self._parse_settings(settings)
        if bool_dashed:
            style['dash'] = dash_list
            style['linewidth'] = linewidth
            
        self._update_bounds(center[0] - radius, center[1] - radius)
        self._update_bounds(center[0] + radius, center[1] + radius)

        def draw(ctx):
            ctx.set_source_rgb(*style['linecolor'])
            ctx.set_line_width(style['linewidth'])
            ctx.set_dash(style['dash'])
            
            ctx.arc(center[0], center[1], radius, 0, 2 * math.pi)
            ctx.stroke()
            
        self.commands.append(draw)

    def pyx_marker(self, loc, size=0.15, rgb=[0,0,0], settings=[]):
        self._update_bounds(loc[0] - size, loc[1] - size)
        self._update_bounds(loc[0] + size, loc[1] + size)
        def draw(ctx):
            ctx.set_source_rgb(*rgb)
            ctx.arc(loc[0], loc[1], size, 0, 2 * math.pi)
            ctx.fill()
        self.commands.append(draw)

    def pyx_marker_plus(self, loc, size=1, rgb=[0,0,0]):
        self._update_bounds(loc[0] - size, loc[1] - size)
        self._update_bounds(loc[0] + size, loc[1] + size)
        def draw(ctx):
            ctx.set_source_rgb(*rgb)
            ctx.set_line_width(global_settings['linewidth']) 
            
            ctx.arc(loc[0], loc[1], size, 0, 2 * math.pi)
            ctx.stroke()
            
            ctx.move_to(loc[0] - size, loc[1])
            ctx.line_to(loc[0] + size, loc[1])
            ctx.move_to(loc[0], loc[1] - size)
            ctx.line_to(loc[0], loc[1] + size)
            ctx.stroke()
            
        self.commands.append(draw)

    def pyx_marker_minus(self, loc, size=1, rgb=[0,0,0]):
        self._update_bounds(loc[0] - size, loc[1] - size)
        self._update_bounds(loc[0] + size, loc[1] + size)
        def draw(ctx):
            ctx.set_source_rgb(*rgb)
            ctx.set_line_width(global_settings['linewidth'])
            
            ctx.arc(loc[0], loc[1], size, 0, 2 * math.pi)
            ctx.stroke()
            
            ctx.arc(loc[0], loc[1], size * 0.33, 0, 2 * math.pi)
            ctx.fill()
            
        self.commands.append(draw)

    def render(self, filename, fmt):
        if self.min_x == float('inf'):
            # Empty canvas
            self.min_x, self.min_y, self.max_x, self.max_y = -10, -10, 10, 10

        # Add some padding
        padding = 5
        min_x = self.min_x - padding
        max_x = self.max_x + padding
        min_y = self.min_y - padding
        max_y = self.max_y + padding
        
        width = max_x - min_x
        height = max_y - min_y
        
        # Convert to points
        surface_width = width * self.scale_factor
        surface_height = height * self.scale_factor
        
        if fmt == 'pdf':
            if not filename.endswith('.pdf'):
                filename += '.pdf'
            surface = cairo.PDFSurface(filename, surface_width, surface_height)
        elif fmt == 'svg':
            if not filename.endswith('.svg'):
                filename += '.svg'
            surface = cairo.SVGSurface(filename, surface_width, surface_height)
        else:
            return

        ctx = cairo.Context(surface)
        
        # Transform to fit bounds
        # Origin (0,0) in user space should map to correct pixel
        # We want min_x, min_y to be at bottom-left (or top-left depending on Y flip)
        # Cairo origin is top-left.
        # We flip Y: scale(s, -s).
        # So +y goes up.
        # We want min_y to be at the bottom of the surface (surface_height).
        # max_y to be at top (0).
        # min_x at 0.
        
        # Translate to move min_x, min_y to origin
        # ctx.translate(-min_x * self.scale_factor, -min_y * self.scale_factor) # This would put min_x, min_y at 0,0 (top-left)
        
        # Let's do it carefully.
        # 1. Translate origin to bottom-left of surface? No, top-left is 0,0.
        # We want logical point (min_x, min_y) to appear at (0, surface_height) if we didn't flip?
        # With flip:
        # ctx.scale(s, -s) -> +y is up.
        # (0,0) becomes (0,0). (0, 10) becomes (0, -10).
        # We need to translate the whole thing down so that max_y is at the top?
        # Or rather, we want the visible window [min_x, max_x] x [min_y, max_y] to map to [0, W] x [0, H].
        
        # Standard math coords to Cairo:
        # x_cairo = (x_user - min_x) * scale
        # y_cairo = height_cairo - (y_user - min_y) * scale
        
        # Let's use matrix transform
        # ctx.translate(0, surface_height)
        # ctx.scale(self.scale_factor, -self.scale_factor)
        # ctx.translate(-min_x, -min_y)
        
        ctx.translate(0, surface_height)
        ctx.scale(self.scale_factor, -self.scale_factor)
        ctx.translate(-min_x, -min_y)
        
        for cmd in self.commands:
            ctx.save()
            cmd(ctx)
            ctx.restore()
            
        surface.finish()
        print(f"File saved to {filename}")

if __name__ == '__main__':
    u = PyX_Utility()
    u.pyx_text([0,0], 'Hello Cairo!')
    u.pyx_circle(5, [0,0])
    u.cvs.writePDFfile('test_cairo_output')
