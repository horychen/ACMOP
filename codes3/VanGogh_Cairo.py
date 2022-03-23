# importing pycairo
import cairo
from pylab import np

class VanGogh_Cairo:
    def __init__(self, acm_variant, width_in_points=500, height_in_points=500):
        self.acm_variant = acm_variant
        self.output_fname_no_suffix = acm_variant.template.fea_config_dict['output_dir'] + acm_variant.name
        self.surface = cairo.SVGSurface(self.output_fname_no_suffix+'.svg', width_in_points, height_in_points)
        self.ctx = cairo.Context(self.surface)
        # self.ctx.scale(width_in_points, height_in_points)

        # m = cairo.Matrix(yy=-1, y0=height_in_points) # Cartetian Coordinate
        m = cairo.Matrix(yy=-1, y0=0.5*height_in_points, x0=+0.5*width_in_points) # Offset to center
        self.ctx.transform(m)

        # Set a background color
        if True:
            self.ctx.save()
            self.ctx.set_source_rgb(0.95, 0.95, 0.95)
            self.ctx.paint()
            self.ctx.restore()
        else:
            pat = cairo.LinearGradient(0.0, 0.0, 0.0, 1.0)
            pat.add_color_stop_rgba(1, 0, 0, 0, 1)
            pat.add_color_stop_rgba(0, 1, 1, 1, 1)
            self.ctx.rectangle(0, 0, width_in_points, height_in_points)
            self.ctx.set_source(pat)
            self.ctx.fill()

            pat = cairo.RadialGradient(0.45, 0.4, 0.1,
                                    0.4, 0.4, 0.5)
            pat.add_color_stop_rgba(0, 1, 1, 1, 1)
            pat.add_color_stop_rgba(1, 0, 0, 0, 1)
            self.ctx.set_source(pat)

    def drawLine(self, p1, p2):
        self.ctx.move_to(p1[0], p1[1])
        self.ctx.line_to(p2[0], p2[1])
        return []

    def drawArc(self, centerxy, startxy, endxy):

        v1 = np.array([startxy[0] - centerxy[0], startxy[1] - centerxy[1]])
        v2 = np.array([endxy[0]   - centerxy[0], endxy[1]   - centerxy[1]])
        cos夹角 = (v1[0]*v2[0] + v1[1]*v2[1]) / (np.sqrt(v1.dot(v1))*np.sqrt(v2.dot(v2)))
        angle_between = np.arccos(cos夹角)

        radius = np.sqrt(v1.dot(v1))
        angle_start = np.arctan2(v1[1], v1[0])
        angle_end = angle_start + angle_between

        self.ctx.move_to(startxy[0], startxy[1])
        self.ctx.arc(centerxy[0], centerxy[1], radius, angle_start, angle_end)
        # self.ctx.arc_negative(centerxy[0], centerxy[1], radius, angle_end, angle_start)
        return []

    def draw_doubly_salient(self, acm_variant):
        # Rotor Core
        list_regions_1 = acm_variant.rotorCore.draw(self, bool_draw_whole_model=True)
        # self.bMirror = False
        # self.iRotateCopy = acm_variant.rotorCore.p*2
        # region1 = self.prepareSection(list_regions_1)

        # Shaft
        # list_regions = acm_variant.shaft.draw(self)
        # self.bMirror = False
        # self.iRotateCopy = 1
        # region0 = self.prepareSection(list_regions)

        # Stator Core
        list_regions = acm_variant.stator_core.draw(self, bool_draw_whole_model=True)
        # self.bMirror = True
        # self.iRotateCopy = acm_variant.stator_core.Q
        # region3 = self.prepareSection(list_regions)

        # Stator Winding
        # list_regions = acm_variant.coils.draw(self, bool_draw_whole_model=True)
        # self.bMirror = False
        # self.iRotateCopy = acm_variant.coils.stator_core.Q
        # region4 = self.prepareSection(list_regions)

        self.apply_stroke()
        self.save()

        if False:
            # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
            EX = acm_variant.template.d['EX']
            CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * EX['WindingFill'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
            CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。

            # Maybe there is a bug here... regarding the excitation for suspension winding...
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
            EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
            EX['DriveW_CurrentAmp'] = acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            EX['BeariW_CurrentAmp'] = acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
            # print('[FEMM_SlidingMesh.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

        return True

    def draw_spmsm(self, acm_variant):
        # Rotor Core
        list_regions_1 = acm_variant.rotorCore.draw(self, bool_draw_whole_model=True)
        # self.bMirror = False
        # self.iRotateCopy = acm_variant.rotorCore.p*2
        # region1 = self.prepareSection(list_regions_1)

        # Shaft
        # list_regions = acm_variant.shaft.draw(self)
        # self.bMirror = False
        # self.iRotateCopy = 1
        # region0 = self.prepareSection(list_regions)

        # Rotor Magnet
        list_regions = acm_variant.rotorMagnet.draw(self, bool_draw_whole_model=True)
        # self.bMirror = False
        # self.iRotateCopy = acm_variant.rotorMagnet.notched_rotor.p*2
        # region2 = self.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

        # Sleeve
        # list_regions = acm_variant.sleeve.draw(self)
        # self.bMirror = False
        # self.iRotateCopy = acm_variant.rotorMagnet.notched_rotor.p*2
        # regionS = self.prepareSection(list_regions)

        # Stator Core
        list_regions = acm_variant.stator_core.draw(self, bool_draw_whole_model=True)
        # self.bMirror = True
        # self.iRotateCopy = acm_variant.stator_core.Q
        # region3 = self.prepareSection(list_regions)

        # Stator Winding
        # list_regions = acm_variant.coils.draw(self, bool_draw_whole_model=True)
        # self.bMirror = False
        # self.iRotateCopy = acm_variant.coils.stator_core.Q
        # region4 = self.prepareSection(list_regions)

        self.apply_stroke()
        self.save()

        if False:
            # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
            EX = acm_variant.template.d['EX']
            CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * EX['WindingFill'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
            CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。

            # Maybe there is a bug here... regarding the excitation for suspension winding...
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
            EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
            EX['DriveW_CurrentAmp'] = acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            EX['BeariW_CurrentAmp'] = acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
            # print('[FEMM_SlidingMesh.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

        return True

    def apply_stroke(self):

        self.ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        self.ctx.set_line_width(0.5)

        # setting color of the context
        self.ctx.set_source_rgba(0.0, 0.0, 0.0, 1)

        # stroke out the color and width property
        self.ctx.stroke()

    def save(self):
        # self.surface.write_to_svg()
        self.surface.finish()
        import cairosvg
        cairosvg.svg2pdf(url=self.output_fname_no_suffix+'.svg', write_to=self.output_fname_no_suffix+'.pdf')
        print(f"[Vangogh_Cairo.py] Cairo plot saved to {self.output_fname_no_suffix+'.svg (and .pdf)'}")

    def getSketch(self, name, color):
        self.name = name
        self.color = color
        # x_bearing, y_bearing, width, height, x_advance, y_advance = self.ctx.text_extents(name)
        # cr.move_to(innerCoord[0], innerCoord[1])
        # cr.show_text(self.name)

        # we can apply text and color
        pass

if __name__ == '__main__':
    # creating a SVG surface
    with cairo.SVGSurface("VanGoghCairo-Demo.svg", 700, 700) as surface:

        # creating a cairo context object
        context = cairo.Context(surface)

        # creating a rectangle(square) for left eye
        context.rectangle(100, 100, 100, 100)

        # creating a rectangle(square) for right eye
        context.rectangle(500, 100, 100, 100)

        # creating position for the curves
        x, y, x1, y1 = 0.1, 0.5, 0.4, 0.9
        x2, y2, x3, y3 = 0.4, 0.1, 0.9, 0.6

        # setting scale of the context
        # context.scale(700, 700)

        # setting line width of the context
        context.set_line_width(2)

        # move the context to x,y position
        context.move_to(x, y)

        # draw the curve for smile
        context.curve_to(x1, y1, x2, y2, x3, y3)

        context.move_to(100, 100)
        context.arc(100, 100, 50, 0, 0.667*np.pi)

        context.move_to(500, 100)
        context.arc_negative(500, 100, 50, 0, 0.667*np.pi)

        # setting color of the context
        context.set_source_rgba(0.4, 1, 0.4, 1)

        # stroke out the color and width property
        context.stroke()

    import cairosvg
    cairosvg.svg2pdf(url="VanGoghCairo-Demo.svg", write_to="VanGoghCairo-Demo.pdf")
