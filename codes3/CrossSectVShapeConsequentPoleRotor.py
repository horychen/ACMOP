from pylab import np, cos, sin
from utility import EPS

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
    def __init__(self, 
                    name = 'V-Shape Consequent Pole Rotor',
                    color = '#677781', #FE840E', # https://www.color-hex.com/ or https://htmlcolorcodes.com/
                    mm_d_bg_air=1.5,
                    mm_d_bg_magnet=2,
                    mm_d_pm=9.4,
                    deg_alpha_vspm=20.3,
                    mm_r_or=129.8,
                    mm_r_ri=15.86,
                    mm_d_ri=84.6,
                    mm_w_pm=20,
                    p = 2, 
                    pr = 10,
                    location = None
                    ):
        self.name = name
        self.color = color
        # BUG: why passed arguments are tuples?
        self.mm_d_bg_air=mm_d_bg_air
        self.mm_d_bg_magnet=mm_d_bg_magnet
        self.mm_d_pm=mm_d_pm
        self.deg_alpha_vspm=deg_alpha_vspm
        self.mm_r_or=mm_r_or
        self.mm_r_ri=mm_r_ri
        self.mm_d_ri=mm_d_ri
        self.mm_w_pm=mm_w_pm
        self.p = p 
        self.pr = pr
        self.location = location         # move this part to another location other than origin (not supported yet)

        # Validate parameters
        # TODO
        # 永磁体倾斜角度不能使得镂空区域超过18°（360/pr/2）

    def draw(self, drawer):

        drawer.getSketch(self.name, self.color)

        p = self.p
        pr = self.pr

        # print(type(self.mm_d_bg_air))
        # print(self.mm_d_bg_air)
        # quit()

        d_bg_air    = self.mm_d_bg_air
        d_bg_magnet = self.mm_d_bg_magnet
        d_pm        = self.mm_d_pm
        alpha_vspm  = self.deg_alpha_vspm / 180 * np.pi
        r_or        = self.mm_r_or
        r_ri        = self.mm_r_ri
        d_ri        = self.mm_d_ri
        w_pm        = self.mm_w_pm

        P1 = [r_ri, 0]

        P2 = [r_ri + d_ri, 0]

        P3 = [P2[0], d_bg_magnet/2]

        P4 = [r_or - d_bg_air, d_bg_magnet/2]

        l34 = r_or - d_bg_air - d_ri - r_ri
        P5 = [ P4[0], l34/cos(alpha_vspm)]

        l37 = d_pm / cos(alpha_vspm)
        P6 = [ P5[0], P5[1]+l37 ] 

        P7 = [ P3[0], P3[1]+l37 ]

        P8 = [  (r_or - d_bg_air) * cos(360/pr/2 / 180*np.pi),
                (r_or - d_bg_air) *-sin(360/pr/2 / 180*np.pi) ]

        P9 = [  P1[0] * cos(360/pr/2 / 180*np.pi),
                P1[0] * sin(360/pr/2 / 180*np.pi) ]

        P11 = [r_or, 0]
        P12 = [ r_or* cos(360/pr/2 / 180*np.pi),
                r_or* sin(360/pr/2 / 180*np.pi) ]


        P5[0] -= 6
        P6[0] -= 6

        list_segments_Core = []
        list_segments_Core += drawer.drawLine(P11, P1)
        list_segments_Core += drawer.drawArc([0,0], P11, P12)
        list_segments_Core += drawer.drawLine(P12, P9)
        list_segments_Core += drawer.drawArc([0,0], P1, P9)
  
        list_segments_Hole = []
        list_segments_Hole += drawer.drawLine(P3, P5)
        list_segments_Hole += drawer.drawLine(P5, P6)
        list_segments_Hole += drawer.drawLine(P6, P7)
        list_segments_Hole += drawer.drawLine(P7, P3)

        innerCoord = ( 0.5*(P4[0]+P3[0]), 0.5*(P4[1]+P3[1]) )

        # return [list_segments] # csToken # cross section token
        return {    'innerCoord': innerCoord, 
                    'list_regions':[list_segments_Core, list_segments_Hole], 
                    'list_regions_to_remove': [0, 1],
                    'mirrorAxis': None,
                }

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
            print('[CrossSectVShapeConsequentPoleRotor.py] Magnet area in total is %g mm^2'%(self.mm2_magnet_area))
            if bool_re_evaluate:
                return self.mm2_magnet_area

            # rotor inter-pole notch
            if self.notched_rotor.deg_alpha_rm >= 180/p*0.9800:
                print('[CrossSectVShapeConsequentPoleRotor.py] FULL POLE PITCH MAGNET IS USED.')
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

    mop.fea_config_dict['designer.Show'] = True
    toolJd = JMAG.JMAG(mop.fea_config_dict, spec_input_dict=mop.spec_input_dict)

    project_name               = '_test_vernier'
    expected_project_file_path = f"./{project_name}.jproj"
    toolJd.open(expected_project_file_path)

    # Rotor core
    vshape_vernier_rotor = CrossSectVShapeConsequentPoleRotor(
                                                name = 'V-Shape Consequent Pole Rotor',
                                                color = '#677781', #FE840E', # https://www.color-hex.com/ or https://htmlcolorcodes.com/
                                                mm_d_bg_air=1.5,
                                                mm_d_bg_magnet=2,
                                                mm_d_pm=9.4,
                                                deg_alpha_vspm=15.3,
                                                mm_r_or=129.8,
                                                mm_r_ri=15.86,
                                                mm_d_ri=84.6,
                                                mm_w_pm=20,
                                                p = 2, 
                                                pr = 10,
                                                location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                                )
    token = vshape_vernier_rotor.draw(toolJd) # drawer

    toolJd.bMirror = True
    toolJd.iRotateCopy = vshape_vernier_rotor.pr
    region1 = toolJd.prepareSection(token)

    # Magnet
    # notched_magnet = CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
    #                                               color = '#0E001E',
    #                                               vshape_vernier_rotor = vshape_vernier_rotor
    #                                             )
    # token = notched_magnet.draw(toolJd)
    # toolJd.bMirror = False
    # toolJd.iRotateCopy = vshape_vernier_rotor.p*2
    # region2 = toolJd.prepareSection(token)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('PMVM Modeling')
    model.SetDescription('PMVM Test')

    # Pre-process
    # toolJd.preProcess(makeToken)
    # model.CloseCadLink()

