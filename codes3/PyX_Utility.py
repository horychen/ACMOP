import pyx
pyx.unit.set(wscale=1) # 同时修改VanGogh.py，搜索本行代码，减小画的点的size（其实可以用线cap属性处理而不需要额外画点的）
# pyx.unit.set(wscale=7) # 同时修改VanGogh.py，搜索本行代码，减小画的点的size（其实可以用线cap属性处理而不需要额外画点的）
pyx.text.set(pyx.text.LatexEngine)  # https://pyx-project.org/examples/text/textengine.html
                                    # https://tex.stackexchange.com/questions/406790/mathbb-plot-label-in-pyx
                                    # obsolete
                                    # pyx.text.set(mode="latex")
                                    # pyx.text.preamble(r"\usepackage{amssymb}") # for mathbb # https://tex.stackexchange.com/questions/406790/mathbb-plot-label-in-pyx
pyx.unit.set(xscale=1.5)


mapping_dict = { 
                'dashed': pyx.style.dash([0,4]),
                'dense-dashed': pyx.style.dash([0,1]),
                'my-thick-line': pyx.style.linewidth(0.2),
               } # introduce the mapping dict for pyx

import math
class PyX_Utility:
    def __init__(self):
        self.cvs = pyx.canvas.canvas()
        self.trace_l = [] # all the line path you have been drawn
        self.trace_a = [] # all the arc path you have been drawn

    def pyx_line(self, p1, p2, bool_track=False, arg=[], **kwarg):

        settings = [mapping_dict[setting] for setting in arg]

        # PyX （必须有pyx.style.linecap.round，才能在Overleaf上显现虚线，否则你在Samantrapdf里看到的，overleaf编译以后没了。)
        self.cvs.stroke(pyx.path.line(p1[0], p1[1], p2[0], p2[1]), [pyx.color.rgb.black, pyx.style.linewidth.THICK, pyx.style.linecap.round]+settings)
        self.cvs.fill(pyx.path.circle(p1[0], p1[1], 0.075)) # use this if THICK is used.
        self.cvs.fill(pyx.path.circle(p2[0], p2[1], 0.075)) # use this if THICK is used.

        if bool_track:
            self.trace_l.append( [*p1, *p2] )
            return self.trace_l[-1]

    def pyx_arc(self, startxy, endxy, centerxy=(0,0), bool_track=False, arg=[], **kwarg):

        settings = [mapping_dict[setting] for setting in arg]

        angle_between_points = math.atan2(endxy[1], endxy[0]) - math.atan2(startxy[1], startxy[0])
        relative_angle = -0.5 * angle_between_points # simple math: https://pyx-project.org/manual/connector.html?highlight=relangle

        # PyX
        textattrs = [pyx.text.halign.center, pyx.text.vshift.middlezero]
        X = pyx.text.text(startxy[0], startxy[1], r"", textattrs) # must have textattrs or else you will have normsubpath cannot close error: AssertionError: normsubpathitems do not match
        Y = pyx.text.text(endxy[0], endxy[1], r"", textattrs)
        self.cvs.stroke(pyx.connector.arc(X, Y, boxdists=[0.0, 0.0], relangle=relative_angle/pi*180), [pyx.color.rgb.black, pyx.style.linewidth.THICK]+settings) # https://pyx-project.org/manual/connector.html?highlight=relangle

        if bool_track:
            self.trace_a.append([*startxy, *endxy, *centerxy, relative_angle])
            return self.trace_a[-1]

    def pyx_text(self, loc, text, size=5, scale=1.0):
        # method 1: insert
        textattrs = [pyx.text.halign.center, # left, center, right
                     pyx.text.vshift.middlezero, 
                     pyx.text.size(size), 
                     pyx.trafo.scale(scale)]
        handle = pyx.text.text(*loc, text, textattrs)
        self.cvs.insert(handle)

        return handle

        # method 2: direct
        # self.cvs.text(*loc, text)

    def pyx_arrow(self, PA, PB=None):
        if PB is None:
            PB = PA
            PA = [0,0]
        self.cvs.stroke(pyx.path.line(PA[0], PA[1], PB[0], PB[1]),
                      [ pyx.style.linewidth.Thick, pyx.style.linestyle.solid, pyx.color.rgb.black,
                        pyx.deco.earrow([ pyx.deco.stroked([pyx.color.rgb.black, pyx.style.linejoin.round]), 
                                          pyx.deco.filled([pyx.color.rgb.black]) 
                                        ] , size=0.5)
                      ])

    def pyx_marker(self, loc, size=0.25, rgb=[0,0,0]):
        r,g,b = rgb[0], rgb[1], rgb[2]
        self.cvs.fill(pyx.path.circle(*loc, size), [pyx.color.rgb(r,g,b)])

    def pyx_marker_minus(self, loc, size=1, rgb=[0,0,0]):
        r,g,b = rgb[0], rgb[1], rgb[2]
        # t = self.pyx_text(loc, r"\PyXMarker{id}")
        # center = t.marker("id")
        # center = [el+0.1 for el in center] # 哎，这个center是歪的
        center = loc
        self.cvs.stroke(pyx.path.circle(center[0], center[1], size), [pyx.color.rgb(r,g,b), pyx.style.linewidth.THICK])
        self.cvs.fill(pyx.path.circle(*center, size*0.33), [pyx.color.rgb(r,g,b)])

    def pyx_marker_plus(self, loc, size=1, rgb=[0,0,0]):
        r,g,b = rgb[0], rgb[1], rgb[2]
        # t = self.pyx_text(loc, r"\PyXMarker{id}")
        # center = t.marker("id")
        center = loc
        self.cvs.stroke(pyx.path.circle(center[0], center[1], size), [pyx.color.rgb(r,g,b), pyx.style.linewidth.THICK])
        self.cvs.stroke(pyx.path.line(center[0]-size, center[1], center[0]+size, center[1]), [pyx.color.rgb(r,g,b), pyx.style.linewidth.THICK])
        self.cvs.stroke(pyx.path.line(center[0], center[1]-size, center[0], center[1]+size), [pyx.color.rgb(r,g,b), pyx.style.linewidth.THICK])

    def pyx_draw_sector(self, origin, angle_begin, angle_end, radius, bool_stroke=False):

        upArc = pyx.path.path(pyx.path.arc(*origin, radius, angle_begin, angle_end))
        right = pyx.path.line(*origin, radius*math.cos(angle_begin/180*math.pi), radius*math.sin(angle_begin/180*math.pi))
        left  = pyx.path.line(radius*math.cos(angle_end/180*math.pi), radius*math.sin(angle_end/180*math.pi), *origin)
        p = right<<upArc<<left
        self.cvs.fill(p, [pyx.color.gray(0.9)])
        if bool_stroke:
            self.cvs.stroke(p, [pyx.style.linewidth.thin]) # comment this for no outline

    def pyx_circle(self, radius, center=[0,0], bool_dashed=False, dash_list=[0,4], linewidth=0.035):
        p = pyx.path.circle(*center, radius)
        if bool_dashed:
            # self.cvs.stroke(p, [pyx.style.linewidth.thin, pyx.style.linestyle.dashed])
            self.cvs.stroke(p, [pyx.style.linecap.round, pyx.style.dash(dash_list), pyx.style.linewidth(linewidth)]) # [0,2] is default # https://stackoverflow.com/questions/50621257/python-pyx-plot-changing-spacing-between-dots-in-dotted-line
            # self.cvs.stroke(p, [pyx.style.linecap.round, pyx.style.dash([0,5]), pyx.style.linestyle.dashdotted, pyx.style.linewidth(0.02)])  # https://pyx-project.org/manual/pathstyles.html
        else:
            self.cvs.stroke(p, [pyx.style.linewidth.thin]) 



if __name__ == '__main__':
    u = PyX_Utility()
    u.pyx_text([0,0], 'Hello world!')
    u.cvs.writePDFfile(r'C:\Users\horyc\Desktop\pyx_output')
    u.cvs.writeSVGfile(r'C:\Users\horyc\Desktop\pyx_output')




















# Used in main_post_processing_pm.py
from pylab import np
import os
from copy import deepcopy
# def pyx_draw_path(tool_tikz, path, sign=1, bool_exclude_path=False):
#     # 这个函数完全可以写成 VanGogh.VanGogh_TikZPlotter() 的成员函数的。。。
#     if bool_exclude_path == False:
#         if len(path) == 4: # line
#             p = tool_tikz.draw_line(path[:2], path[2:4], untrack=True)
#         else: # == 6 for arc
#             p = tool_tikz.draw_arc(path[:2], path[2:4], path[4:6], relangle=sign*path[6], untrack=True, ccw=sign)
#         return p
def rotate(_, x, y):
    return np.cos(_)*x + np.sin(_)*y, -np.sin(_)*x + np.cos(_)*y
def is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):
    return np.sqrt(path[0]**2 + path[1]**2) > mm_rotor_outer_radius + 0.5*mm_air_gap_length

def pyx_script_pmsm(ad, best_chromosome, best_idx, Q, p, proj_name):
    if proj_name is None:
        name = 'Q%dp%didx%d'%(Q,p,best_idx)
    else:
        name = 'Q%dp%didx%d%s'%(Q,p,best_idx, proj_name)
    
    # output location 1
    filename_cross_sectional_view = "../_pemd2020/%s.svg"%("Figure_selected_optimal_design_%s"%(name))
    # output location 2
    filename_cross_sectional_view = ad.solver.output_dir + "%s.svg"%("Figure_selected_optimal_design_%s"%(name))
    # output location 3
    filename_cross_sectional_view = '../release/' + "%s.svg"%("Figure_selected_optimal_design_%s"%(name))

    # Plot cross section view
    import bearingless_spmsm_design
    spmsm_best = bearingless_spmsm_design.bearingless_spmsm_design(
                                        spmsm_template=ad.spec.acm_template,
                                        x_denorm=best_chromosome[:-3],
                                        counter=999,
                                        counter_loop=1
                                        )
    spmsm_best.ID = name
    import VanGogh
    tool_tikz = VanGogh.VanGogh_TikZPlotter()
    spmsm_best.draw_spmsm(tool_tikz, bool_pyx=True) # collecting track_path list for tool_tikz

    # 就算可以直接画出来，也只有最小部分而已，没有做镜像和旋转！
    # tool_tikz.c.writePDFfile("../test"+proj_name)
    # tool_tikz.c.writeEPSfile("Figure_selected_optimal_design_%s"%(name))
    tool_tikz.c.writePDFfile(filename_cross_sectional_view[:-4])
    tool_tikz.c.writeSVGfile("../a")
    # print('Write to pdf file:', "../test"+proj_name)
    # raise Exception('DEBUG HERE')
    tool_tikz.c = pyx.canvas.canvas()

    # Use tool_tikz.track_path to redraw the model
    def redraw_cross_section_outline_with_pyx(tool_tikz, no_repeat_stator, no_repeat_rotor, mm_rotor_outer_radius, mm_air_gap_length, mm_rotor_outer_steel_radius, mm_rotor_inner_radius):
        # PyX
        tool_tikz.c = pyx.canvas.canvas() # clear the canvas because we want to redraw 90 deg with the data tool_tikz.track_path

        print('Index   | Path data')
        p_stator = None #pyx.path.path()
        p_rotor  = None #pyx.path.path()
        for index, path in enumerate(tool_tikz.track_path): # track_path is passed by reference and is changed by mirror
            # Failed to fill the closed path, because there is no arc-like path available.
            # p = pyx.path.line(4, 0, 5, 0) << pyx.path.line(5, 0, 5, 1) << pyx.path.line(5, 1, 4, 1)
            # p.append(path.closepath())
            # tool_tikz.c.stroke(p)
            # tool_tikz.c.stroke(path.rect(0, 0, 1, 1), [pyx.style.linewidth.Thick,
            #                      pyx.color.rgb.red,
            #                      pyx.deco.filled([pyx.color.rgb.green])])

            path_mirror = deepcopy(path)
            # for mirror copy (along x-axis)
            path_mirror[1] = path[1]*-1
            path_mirror[3] = path[3]*-1
            # for mirror copy (along y-axis)
            # path_mirror[0] = path[0]*-1
            # path_mirror[2] = path[2]*-1

            bool_exclude_path = False

            # rotate path and plot
            if is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):
                Q = no_repeat_stator
            else:
                Q = no_repeat_rotor*2



            EPS = 1e-6
            if is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):
                # 按照Eric的要求，把不必要的线给删了。
                if abs(path[1] + path[3]) < EPS: # 镜像对称线
                    bool_exclude_path = True
                if abs(path[0] - path[2]) + np.cos(2*np.pi/Q/2) < EPS: # 旋转对称线（特别情况，tan(90°) = ∞
                    bool_exclude_path = True                
                else:
                    if abs( abs((path[1] - path[3])/(path[0] - path[2])) - abs(np.tan(2*np.pi/Q/2)) ) < EPS: # 旋转对称线
                        bool_exclude_path = True

            if not is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):
                # 按照Eric的要求，把不必要的线给删了。
                if  (abs(np.sqrt(path[0]**2+path[1]**2) - mm_rotor_inner_radius)<EPS or abs(np.sqrt(path[2]**2+path[3]**2) - mm_rotor_inner_radius)<EPS) \
                    and (len(path)==4): # 转子铁芯内径到外径的直线(len(path)==4)
                    bool_exclude_path = True

            #     # 特别的是，画永磁体的时候，边界要闭合哦。
            #     if abs(np.sqrt(path[0]**2+path[1]**2) - mm_rotor_outer_steel_radius) < EPS or abs(np.sqrt(path[2]**2+path[3]**2) - mm_rotor_outer_steel_radius) < EPS:
            #         bool_exclude_path = False

            # A trick that makes sure models with different outer diameters have the same scale.
            # tool_tikz.draw_arc([125,0], [-125,0], relangle=sign*180, untrack=True)
            tool_tikz.c.fill(pyx.path.circle(0, 0, 125), [pyx.color.transparency(1)]) # use this if THICK is used. <- Warn: Transparency not available in PostScript, proprietary ghostscript extension code inserted. (save as eps format)
            # tool_tikz.c.fill(pyx.path.circle(0, 0, 125), [pyx.color.rgb.white]) # use this if THICK is used. <- this will over-write everthing... how to change zorder?


            _ = 2*np.pi/Q
 
            if True: # full model
                for counter in range(Q):

                    # 转子：旋转复制
                    if not is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):
                        path[0], path[1] = rotate(_, path[0], path[1])
                        path[2], path[3] = rotate(_, path[2], path[3])
                        路径 = tool_tikz.pyx_draw_path(path, sign=1, bool_exclude_path=bool_exclude_path, bool_stroke=True)
                        if 路径 is not None:
                            if p_rotor is None:
                                p_rotor = 路径
                            else:
                                p_rotor = p_rotor << 路径

                    # 定子：镜像+旋转复制
                    if is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):

                        path[0], path[1] = rotate(_, path[0], path[1])
                        path[2], path[3] = rotate(_, path[2], path[3])
                        路径 = tool_tikz.pyx_draw_path(path, sign=1, bool_exclude_path=bool_exclude_path, bool_stroke=True)
                        # print(index, '\t|', ',\t'.join(['%g'%(el) for el in path]))
                        if 路径 is not None:
                            if p_stator is None:
                                p_stator = 路径
                            else:
                                p_stator = p_stator << 路径

                        path_mirror[0], path_mirror[1] = rotate(_, path_mirror[0], path_mirror[1])
                        path_mirror[2], path_mirror[3] = rotate(_, path_mirror[2], path_mirror[3])
                        路径 = tool_tikz.pyx_draw_path(path_mirror, sign=-1, bool_exclude_path=bool_exclude_path, bool_stroke=True)
                        if 路径 is not None:
                            if p_stator is None:
                                p_stator = 路径
                            else:
                                p_stator = p_stator << 路径

                    # break
            else: # backup

                # 转子：旋转复制
                if not is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):
                    path[0], path[1] = rotate(0.5*np.pi - 0.5*0.5*_, path[0], path[1])
                    path[2], path[3] = rotate(0.5*np.pi - 0.5*0.5*_, path[2], path[3])
                    pyx_draw_path(tool_tikz, path, sign=1, bool_exclude_path=bool_exclude_path)
                    # print(index, '\t|', ',\t'.join(['%g'%(el) for el in path]))

                    # path[0], path[1] = rotate(0.5*np.pi - 0*0.5*_, path[0], path[1])
                    # path[2], path[3] = rotate(0.5*np.pi - 0*0.5*_, path[2], path[3])
                    # pyx_draw_path(tool_tikz, path, sign=1, bool_exclude_path=bool_exclude_path)

                # 定子：镜像+旋转复制
                if is_at_stator(path, mm_rotor_outer_radius, mm_air_gap_length):

                    path[0], path[1] = rotate(0.5*np.pi - 0.5*_, path[0], path[1])
                    path[2], path[3] = rotate(0.5*np.pi - 0.5*_, path[2], path[3])
                    pyx_draw_path(tool_tikz, path, sign=1, bool_exclude_path=bool_exclude_path)
                    # print(index, '\t|', ',\t'.join(['%g'%(el) for el in path]))

                    path_mirror[0], path_mirror[1] = rotate(0.5*np.pi - 0.5*_, path_mirror[0], path_mirror[1])
                    path_mirror[2], path_mirror[3] = rotate(0.5*np.pi - 0.5*_, path_mirror[2], path_mirror[3])
                    pyx_draw_path(tool_tikz, path_mirror, sign=-1, bool_exclude_path=bool_exclude_path)

                    # 注意，所有 tack_path 中的 path 都已经转动了90度了！
                    # for mirror copy (along y-axis)
                    path[0] *= -1
                    path[2] *= -1
                    pyx_draw_path(tool_tikz, path, sign=-1, bool_exclude_path=bool_exclude_path)

                    path_mirror[0] *= -1
                    path_mirror[2] *= -1
                    pyx_draw_path(tool_tikz, path_mirror, sign=1, bool_exclude_path=bool_exclude_path)

        # You can have a cool logo if you un-comment this...
        # tool_tikz.c.fill(p_stator, [pyx.color.gray(0.8)])
        # tool_tikz.c.fill(p_rotor, [pyx.color.gray(0.4)])
        # tool_tikz.c.stroke(p_stator)
        # tool_tikz.c.stroke(p_rotor)

    if True: # Draw the outline?
        redraw_cross_section_outline_with_pyx(tool_tikz, spmsm_best.Q, spmsm_best.p, spmsm_best.Radius_OuterRotor, spmsm_best.Length_AirGap, spmsm_best.Radius_OuterRotorSteel, spmsm_best.Radius_InnerRotor)
        # tool_tikz.c.writePDFfile("../Test_Fill_Plot" + proj_name)
        tool_tikz.c.writeSVGfile("../b")
        tool_tikz.c.writePDFfile(filename_cross_sectional_view[:-4]+'outline')
        print('Final cross sectional outline files are printed to', filename_cross_sectional_view[:-4]+'outline')
            # tool_tikz.c.writePDFfile("Figure_selected_optimal_design_%s"%(name))
        # tool_tikz.c.writeEPSfile("Figure_selected_optimal_design_%s"%(name))
        # tool_tikz.c.writeSVGfile("Figure_selected_optimal_design_%s"%(name))
        # print('Write to pdf file: Figure_selected_optimal_design_%s.pdf.'%(name))
        # os.system('start %s'%("selected_optimal_design%s.pdf"%(spmsm_best.name)))
        # quit()



        # Option 1
        if True:
            from svgutils.compose import Figure, SVG # https://svgutils.readthedocs.io/en/latest/tutorials/composing_multipanel_figures.html

            Figure("5cm", "5cm",
                    SVG("../a.svg"),
                    SVG("../b.svg"), # the order of a and b matters.
                    ).save(filename_cross_sectional_view)
            os.system('inkscape --file=%s --export-area-drawing --without-gui --export-pdf=%s.pdf'%(filename_cross_sectional_view, filename_cross_sectional_view[:-4]))

        else:

            # Option 2
            import svgutils.transform as sg

            #create new SVG figure
            fig = sg.SVGFigure("16cm", "6.5cm")

            # load matpotlib-generated figures
            fig1 = sg.fromfile('../a.svg')
            fig2 = sg.fromfile('../b.svg')

            # get the plot objects
            plot1 = fig1.getroot()
            plot2 = fig2.getroot()
            # plot2.moveto(280, 0, scale=0.5)

            # append plots and labels to figure
            fig.append([plot1, plot2])
            fig.save("../d.svg")

def pyx_script_im(ad, best_chromosome, best_idx, Q, p, proj_name):
    if proj_name is None:
        name = 'Q%dp%didx%d'%(Q,p,best_idx)
    else:
        name = 'Q%dp%didx%d%s'%(Q,p,best_idx, proj_name)
    
    # output location 1
    filename_cross_sectional_view = "../_pemd2020/%s.svg"%("Figure_selected_optimal_design_IM_%s"%(name))
    # output location 2
    filename_cross_sectional_view = ad.solver.output_dir + "%s.svg"%("Figure_selected_optimal_design_IM_%s"%(name))
    # output location 3
    filename_cross_sectional_view = '../release/' + "%s.svg"%("Figure_selected_optimal_design_IM_%s"%(name))

    # Below: Plot cross section view
    # Below: Plot cross section view
    # Below: Plot cross section view

    import population
    im_best = population.bearingless_induction_motor_design.local_design_variant(ad.spec.acm_template, 0, 999, x_denorm=best_chromosome[:-3])
    im_best.name = name
    im_best.ID = str(best_idx)

    from utility_moo import pyx_draw_model
    pyx_draw_model(im_best)
