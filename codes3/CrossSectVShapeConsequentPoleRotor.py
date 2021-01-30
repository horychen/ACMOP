from pylab import np, cos, sin
EPS = 1e-3 # [mm]

# pole_slot_combinations = {
#     0: ("fixed", "p",   "torque_pole_pair_number", None),
#     1: ("fixed", "ps",  "suspension_pole_pair_number", None),
#     2: ("fixed", "pr",  "rotor_pole_pair_number", None),
# }
# rotor_geometric_parameters = {
#     "pr": ("fixed", "rotor_pole_pair_number",              "pr",           10),
#     0: ("free",  "air_bridge_depth",                    "mm_d_bg_air",  1.5),
#     1: ("free",  "magnet_bridge_depth",                 "mm_d_bg_magnet", 5),
#     2: ("free",  "magnet_depth",                        "mm_d_pm",      9.4),
#     3: ("free",  "magnet_rotation",                     "deg_alpha_pm", 20.3),
#     4: ("free",  "rotor_outer_radius",                  "mm_r_or",      129.8),
#     5: ("free",  "rotor_inter_radius",                  "mm_r_ir",      84.6),    
#     6: ("free",  "rotor_inner_depth (back-iron depth)", "mm_d_ri",      15.86),
#     7: ("derived", "magnet_width"                     , "mm_w_pm",      None)
# }
# "D:\DrH\[00]GetWorking\118 VernierMotor\导出永磁体宽度w_pm.afx"
# magnet_width = rotor_geometric_parameters[8] = ( r_os - d_bg_air - (r_ir + d_ir) ) / cos(alpha_pm) - d_pm * tan(alpha_pm)

class CrossSectVShapeConsequentPoleRotor(object):
    # CrossSectInnerNotchedRotor Describes the inner notched rotor.
    #    Properties are set upon class creation and cannot be modified.
    #    The anchor point for this is the center of the rotor,
    #    with the x-axis directed along the center of one of the rotor poles
    def __init__(self, 
                    name = 'V-Shape Consequent Pole Rotor',
                    color = '#677781', #FE840E', # https://www.color-hex.com/ or https://htmlcolorcodes.com/
                    # mm_d_pm = 6,
                    # deg_alpha_rm = 60,
                    # deg_alpha_rs = 10,
                    # mm_d_ri = 8,
                    # mm_r_ri = 40,
                    # mm_d_rp = 5,
                    # mm_d_rs = 3,
                    p = 2, # Set pole-pairs to 2
                    # s = 4, # Set magnet segments/pole to 4
                    location = None
                    ):
        self.name = name
        self.color = color
        self.mm_d_pm = mm_d_pm # depth of the permanent magnet
        self.deg_alpha_rm = deg_alpha_rm # angular span of the pole: class type DimAngular
        self.deg_alpha_rs = deg_alpha_rs # segment span: class type DimAngular
        self.mm_d_ri = mm_d_ri           # inner radius of rotor: class type DimLinear
        self.mm_r_ri = mm_r_ri           # rotor iron thickness: class type DimLinear
        self.mm_d_rp = mm_d_rp           # interpolar iron thickness: class type DimLinear
        self.mm_d_rs = mm_d_rs           # inter segment iron thickness: class type DimLinear
        self.p = p                       # number of pole pairs
        self.s = s                       # number of segments  
        self.location = location         # move this part to another location other than origin (not supported yet)

        # Validate parameters
        # TODO

    def draw(self, drawer):

        drawer.getSketch(self.name, self.color)

        mm_d_pm  = self.mm_d_pm
        alpha_rm = self.deg_alpha_rm * np.pi/180
        alpha_rs = self.deg_alpha_rs * np.pi/180
        r_ri     = self.mm_r_ri
        d_ri     = self.mm_d_ri
        d_rp     = self.mm_d_rp
        d_rs     = self.mm_d_rs
        p        = self.p
        s        = self.s
        alpha_rp = 2*np.pi/(2*p) # pole span



        P1 = [r_ri, 0]

        r_P2 = r_ri + d_ri + d_rp
        P2 = [r_P2, 0]

        alpha_P3 = alpha_rp - alpha_rm
        P3 = [r_P2*cos(alpha_P3), r_P2*-sin(alpha_P3)]

        r_P4 = r_ri + d_ri # = (r_P2 - d_rp) 
        P4 = [r_P4*cos(alpha_P3), r_P4*-sin(alpha_P3)]

        alpha_P5 = alpha_P3 + alpha_rs # alpha_rs means rotor segment (of PM)
        if abs(alpha_rs*s - alpha_rm)<EPS: # This means the inter-segment notch should span 0 deg, which mans the segmented design is reduced to a non-segmented design such that alpha_rm == alpha_rs*s
            alpha_P5 = alpha_P3 + alpha_rs*s
            s = 1 # this is a bad practice but it will help to re-use the drawing code of case s==1 below
        P5 = [r_P4*cos(alpha_P5), r_P4*-sin(alpha_P5)]

        list_segments = []
        if s == 1:
            # No magnet sement!
            # Then P6 is an extra point for a full rotor
            P6 = [r_ri*cos(alpha_P5), r_ri*-sin(alpha_P5)]

            if alpha_rm >= np.pi/p*0.9800:
                print('Non-NOTCHED ROTOR IS USED.\n')
                print('alpha_P5 is', alpha_P5, alpha_P5/np.pi*180)
                P1p5 = [P2[0] - d_rp, P2[1]]
                list_segments += drawer.drawLine(P1, P1p5)
                list_segments += drawer.drawArc([0,0], P5, P1p5)
                list_segments += drawer.drawLine(P5, P6)
                list_segments += drawer.drawArc([0,0], P6, P1)
            else:
                list_segments += drawer.drawLine(P1, P2)
                list_segments += drawer.drawArc([0,0], P3, P2)
                list_segments += drawer.drawLine(P3, P4)
                list_segments += drawer.drawArc([0,0], P5, P4)
                list_segments += drawer.drawLine(P5, P6)
                list_segments += drawer.drawArc([0,0], P6, P1)
                # debug
                # def print_point(P):
                #     print( '(%g, %g)' % (P[0], P[1]) )
                # print_point(P1)
                # print_point(P2)
                # print_point(P3)
                # print_point(P4)
                # print_point(P5)
                # print_point(P6)
        else:

            alpha_notch  = (alpha_rm - s*alpha_rs) / (s-1) # inter-segment notch占的弧度
            P5

            r_P6 = r_ri + d_ri + d_rs
            P6 = [r_P6*cos(alpha_P5), r_P6*-sin(alpha_P5)]

            alpha_P7 = alpha_P5 + alpha_notch
            P7 = [r_P6*cos(alpha_P7), r_P6*-sin(alpha_P7)]

            P8 = [r_P4*cos(alpha_P7), r_P4*-sin(alpha_P7)]

            if alpha_rm >= np.pi/p*0.9800: # no inter-pole notch
                P1p5 = [r_ri + d_ri, 0]
                list_segments += drawer.drawLine(P1, P1p5)
                list_segments += drawer.drawArc([0,0], P5, P1p5)
                list_segments += drawer.drawLine(P5, P6)
                list_segments += drawer.drawArc([0,0], P7, P6)
                list_segments += drawer.drawLine(P7, P8)
            else:
                list_segments += drawer.drawLine(P1, P2)
                list_segments += drawer.drawArc([0,0], P3, P2)
                list_segments += drawer.drawLine(P3, P4)
                list_segments += drawer.drawArc([0,0], P5, P4)
                list_segments += drawer.drawLine(P5, P6)
                list_segments += drawer.drawArc([0,0], P7, P6)
                list_segments += drawer.drawLine(P7, P8)
            # raise
            for _ in range(1, s-1):
                alpha_temp = (alpha_rs+alpha_notch) # alpha_rs means rotor segment (of PM)

                P4 = P8 #[  cos(alpha_temp)*P4[0] + -sin(alpha_temp)*P4[0],
                        #sin(alpha_temp)*P4[1] + cos(alpha_temp)*P4[1] ]

                P5 = [  cos(alpha_temp)*P5[0] + sin(alpha_temp)*P5[1],
                        -sin(alpha_temp)*P5[0] + cos(alpha_temp)*P5[1] ]

                P6 = [  cos(alpha_temp)*P6[0] + sin(alpha_temp)*P6[1],
                        -sin(alpha_temp)*P6[0] + cos(alpha_temp)*P6[1] ]

                P7 = [  cos(alpha_temp)*P7[0] + sin(alpha_temp)*P7[1],
                        -sin(alpha_temp)*P7[0] + cos(alpha_temp)*P7[1] ]

                P8 = [  cos(alpha_temp)*P8[0] + sin(alpha_temp)*P8[1],
                        -sin(alpha_temp)*P8[0] + cos(alpha_temp)*P8[1] ]

                list_segments += drawer.drawArc([0,0], P5, P4)
                list_segments += drawer.drawLine(P5, P6)
                list_segments += drawer.drawArc([0,0], P7, P6)
                list_segments += drawer.drawLine(P7, P8)

            P9 = [  cos(alpha_rs)*P8[0] + sin(alpha_rs)*P8[1],
                    -sin(alpha_rs)*P8[0] + cos(alpha_rs)*P8[1] ] # alpha_rs means rotor segment (of PM)
            P10 = [ r_ri*cos(alpha_rp), r_ri*-sin(alpha_rp)] # alpha_rp is pole span

            list_segments += drawer.drawArc([0,0], P9, P8)
            list_segments += drawer.drawLine(P9, P10)

            list_segments += drawer.drawArc([0,0], P10, P1)

        innerCoord = ( 0.5*(P1[0]+P3[0]), 0.5*(P1[1]+P3[1]))

        # return [list_segments] # csToken # cross section token
        return {'innerCoord': innerCoord, 'list_regions':[list_segments], 'mirrorAxis': None}

class CrossSectVShapeConsequentMagnet(object):
    def __init__(self, 
                    name = 'Notched Rotor',
                    color = '#0BA0E2',
                    notched_rotor = None,
                    ):
        self.name = name
        self.color = color
        self.notched_rotor = notched_rotor

    def draw(self, drawer, bool_re_evaluate=False):

        if False == bool_re_evaluate:
            drawer.getSketch(self.name, self.color)

        d_pm     = self.notched_rotor.mm_d_pm
        alpha_rm = self.notched_rotor.deg_alpha_rm * np.pi/180
        alpha_rs = self.notched_rotor.deg_alpha_rs * np.pi/180
        r_ri     = self.notched_rotor.mm_r_ri
        d_ri     = self.notched_rotor.mm_d_ri
        d_rp     = self.notched_rotor.mm_d_rp
        d_rs     = self.notched_rotor.mm_d_rs
        p        = self.notched_rotor.p
        s        = self.notched_rotor.s
        alpha_rp = 2*np.pi/(2*p) # pole span

        if abs(alpha_rp - alpha_rm) <= 2 * np.pi/180: # if alpha_rm and alpha_rp has a difference smaller than 2 deg, then let alpha_rm equal to alpha_rp.
            alpha_rm = alpha_rp
            if s == 1:
                alpha_rs = alpha_rm # alpha_rs is the variable actually being used in the following...
            else:
                print('[Warn] s=%d: This is not tested. For now it simply assumes the iron notch between poles becomes the iron notch between the segments of one pole.' % (s))
            print('[Warn] [class CrossSectInnerNotchedMagnet] Magnet is fully spanned.')

        if d_pm + 2*EPS < d_rp:
            print('[Warn]: [class CrossSectInnerNotchedMagnet] Detect d_rp is too close to d_pm. To avoid small line entity error in JMAG, set d_pm equal to d_rp because rotor core is plotted already.')
            raise ExceptionBadDesign('[Error] Magnet depth d_pm is too close to inter-pole notch depth d_rp.')

        P1 = [r_ri, 0]

        r_P2 = r_ri + d_ri + d_rp
        P2 = [r_P2, 0]

        alpha_P3 = alpha_rp - alpha_rm
        # P3 = [r_P2*cos(alpha_P3), r_P2*-sin(alpha_P3)]

        r_P4 = r_ri + d_ri # = (r_P2 - d_rp) 
        P4 = [r_P4*cos(alpha_P3), r_P4*-sin(alpha_P3)]

        P3_extra = [(r_P4+d_pm)*cos(alpha_P3), (r_P4+d_pm)*-sin(alpha_P3)]

        alpha_P5 = alpha_P3 + alpha_rs # alpha_rs means rotor segment (of PM)
        P5 = [r_P4*cos(alpha_P5), r_P4*-sin(alpha_P5)]

        list_regions = []
        list_segments = []
        if True:

            if s>1:
                alpha_notch  = (alpha_rm - s*alpha_rs) / (s-1) # 永磁体占的弧度
            P6_extra = [(r_P4+d_pm)*cos(alpha_P5), (r_P4+d_pm)*-sin(alpha_P5)]

            Rout = r_P4+d_pm
            Rin  = r_P4
            self.mm2_magnet_area = alpha_rm/alpha_rp  *  np.pi*(Rout**2 - Rin**2) # magnet area for all the poles
            print('Magnet area in total is %g mm^2'%(self.mm2_magnet_area))
            if bool_re_evaluate:
                return self.mm2_magnet_area

            # rotor inter-pole notch
            if self.notched_rotor.deg_alpha_rm >= 180/p*0.9800:
                print('FULL POLE PITCH MAGNET IS USED.')
                list_segments += drawer.drawLine(P3_extra, P4)
                list_segments += drawer.drawArc([0,0], P5, P4)
                list_segments += drawer.drawLine(P5, P6_extra)
                list_segments += drawer.drawArc([0,0], P6_extra, P3_extra)
            else:
                list_segments += drawer.drawLine(P3_extra, P4)
                list_segments += drawer.drawArc([0,0], P5, P4)
                list_segments += drawer.drawLine(P5, P6_extra)
                list_segments += drawer.drawArc([0,0], P6_extra, P3_extra)

            list_regions.append(list_segments)
            list_segments = []

            # raise
            for _ in range(s-1):
                alpha_temp = (alpha_rs+alpha_notch) # alpha_rs means rotor segment (of PM)

                P3_extra = [  cos(alpha_temp)*P3_extra[0] + sin(alpha_temp)*P3_extra[1],
                             -sin(alpha_temp)*P3_extra[0] + cos(alpha_temp)*P3_extra[1] ]

                P4 = [  cos(alpha_temp)*P4[0] + sin(alpha_temp)*P4[1],
                        -sin(alpha_temp)*P4[0] + cos(alpha_temp)*P4[1] ]

                P5 = [  cos(alpha_temp)*P5[0] + sin(alpha_temp)*P5[1],
                        -sin(alpha_temp)*P5[0] + cos(alpha_temp)*P5[1] ]

                P6_extra = [  cos(alpha_temp)*P6_extra[0] + sin(alpha_temp)*P6_extra[1],
                             -sin(alpha_temp)*P6_extra[0] + cos(alpha_temp)*P6_extra[1] ]

                list_segments += drawer.drawLine(P3_extra, P4)
                list_segments += drawer.drawArc([0,0], P5, P4)
                list_segments += drawer.drawLine(P5, P6_extra)
                list_segments += drawer.drawArc([0,0], P6_extra, P3_extra)

                list_regions.append(list_segments)
                list_segments = []

        innerCoord = ( 0.5*(P4[0]+P6_extra[0]), 0.5*(P4[1]+P6_extra[1]))
        return {'innerCoord': innerCoord, 'list_regions':list_regions, 'mirrorAxis': None}
        # return list_regions # csToken # cross section token

class CrossSectShaft(object):
    def __init__(self, 
                    name = 'Shaft',
                    color = '#0EE0E2',
                    notched_rotor = None,
                    ):
        self.name = name
        self.color = color
        self.notched_rotor = notched_rotor

    def draw(self, drawer):

        drawer.getSketch(self.name, self.color)

        r_ri     = self.notched_rotor.mm_r_ri

        P1 = [r_ri, 0]
        NP1 = [-r_ri, 0]

        list_regions = []
        list_segments = []
        if True:
            list_segments += drawer.drawArc([0,0], NP1, P1)
            list_segments += drawer.drawArc([0,0], P1, NP1)

            list_regions.append(list_segments)
            list_segments = []

        innerCoord = ( 0, 0 )

        # return list_regions # csToken # cross section token
        return {'innerCoord': innerCoord, 'list_regions':list_regions, 'mirrorAxis': None}

if __name__ == '__main__':
    import JMAG
    import Location2D
    import sys; sys.path.insert(0, '../')
    import acmop

    mop = acmop.AC_Machine_Optiomization_Wrapper(
            select_fea_config_dict = "#03 JMAG Non-Nearingless Motor Evaluation Setting",
            select_spec            = "PMVM p2pr10-Q12y3 Wenbo",
            project_loc            = r'D:/DrH/acmop/_WenboVShapeVernier/'
        )

    toolJd = JMAG.JMAG(mop.fea_config_dict, spec_input_dict=mop.spec_input_dict)

    project_name               = '_test_vernier'
    expected_project_file_path = f"./{project_name}.jproj"
    toolJd.open(expected_project_file_path)

    # %% Define cross sections
    notched_rotor = CrossSectInnerNotchedRotor( name = 'NotchedRotor',
                                                color = '#FE840E',
                                                deg_alpha_rm = 60,
                                                deg_alpha_rs = 10,
                                                mm_d_ri = 8,
                                                mm_r_ri = 40,
                                                mm_d_rp = 5,
                                                mm_d_rs = 3,
                                                p = 2, # Set pole-pairs to 2
                                                s = 4, # Set magnet segments/pole to 4
                                                location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                                )

    list_regions = notched_rotor.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = notched_rotor.pr*2
    region1 = toolJd.prepareSection(list_regions)
    
    # Magnet
    # notched_magnet = CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
    #                                               color = '#0E001E',
    #                                               notched_rotor = notched_rotor
    #                                             )

    # list_regions = notched_magnet.draw(toolJd)
    # toolJd.bMirror = False
    # toolJd.iRotateCopy = notched_rotor.p*2
    # region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('PMVM Modeling')
    model.SetDescription('PMVM Test')

    # Pre-process
    # toolJd.preProcess(makeToken)
    # model.CloseCadLink()



