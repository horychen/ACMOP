from pylab import np, cos, sin
EPS = 1e-3 # [mm]

class ExceptionBadDesign(Exception):
    """Exception for notifying bad design."""
    def __init__(self, message, payload=None):
        self.message = message
        self.payload = 'you could add more args here'
    def __str__(self):
        return str(self.message)

class CrossSectInnerSalientPoleRotorV1(object):
    def __init__(self, 
                    name = 'SalientPoleRotor',
                    color = '#FE840E',
                    mm_r_ro = 40,
                    mm_d_sleeve = 1.5,
                    split_ratio_rotor_salient = 0.2,
                    deg_alpha_rsp = 10,
                    pm = 2,
                    mm_r_ri = 5,
                    location = None
                    ):
        self.name = name
        self.color = color
        self.mm_r_ro = mm_r_ro
        self.mm_d_sleeve = mm_d_sleeve
        self.split_ratio_rotor_salient = split_ratio_rotor_salient
        self.deg_alpha_rsp = deg_alpha_rsp
        self.pm = pm                     # number of rotor modulator pole pairs
        self.mm_r_ri = mm_r_ri
        self.location = location         # move this part to another location other than origin (not supported yet)

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        mm_r_ro = self.mm_r_ro
        mm_d_sleeve = self.mm_d_sleeve
        split_ratio_rotor_salient = self.split_ratio_rotor_salient
        mm_d_rsp = split_ratio_rotor_salient * mm_r_ro
        # print('DEBUG: mm_d_rsp=', mm_d_rsp)
        deg_alpha_rsp = self.deg_alpha_rsp
        pm = self.pm
        alpha_rp = 2*np.pi/(2*pm) # pole span

        P1 = [mm_r_ro * cos(0.5*deg_alpha_rsp/180*np.pi), - mm_r_ro * sin(0.5*deg_alpha_rsp/180*np.pi)]
        P2 = [mm_r_ro * cos(0.5*deg_alpha_rsp/180*np.pi), + mm_r_ro * sin(0.5*deg_alpha_rsp/180*np.pi)]
        P3 = [P2[0] - mm_d_rsp, P2[1]]
        P0 = [P1[0] - mm_d_rsp, P1[1]]
        def iPark(P, theta):
            return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
        P4 = iPark(P0, alpha_rp)

        list_segments = []
        if bool_draw_whole_model:
            def draw_fraction(list_segments, P1, P2, P3, P4):
                # P1 = iPark(P1, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                # P2 = iPark(P2, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                # P3 = iPark(P3, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                # P4 = iPark(P4, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                list_segments += drawer.drawArc([0,0], P1, P2)
                list_segments += drawer.drawLine(P2, P3)
                list_segments += drawer.drawArc([0,0], P3, P4)
                P1_CCW = iPark(P1, alpha_rp)
                list_segments += drawer.drawLine(P1_CCW, P4)
            for i in range(2*pm):
                draw_fraction(list_segments, iPark(P1, i*alpha_rp), 
                                             iPark(P2, i*alpha_rp), 
                                             iPark(P3, i*alpha_rp), 
                                             iPark(P4, i*alpha_rp))
            # draw a circle (this is officially suggested by FEMM)
            PRI = [self.mm_r_ri, 0 ]
            list_segments += drawer.drawArc([0,0], PRI, [-PRI[0], PRI[1]])
            list_segments += drawer.drawArc([0,0],      [-PRI[0], PRI[1]], PRI)

            innerCoord = ( 0.5*(P1[0]+P4[0]), 0.5*(P1[1]+P4[1]))

            # return [list_segments] # csToken # cross section token
            return {'innerCoord': innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': None,
                    'inner_or_outer_region_to_remove': [True, False]
                    }
        else:
            raise Exception('Not supported. Please set bool_draw_whole_model to True.')

class CrossSectInnerSalientPoleRotorV2(object):
    def __init__(self, 
                    name = 'SalientPoleRotor',
                    color = '#FE840E',
                    mm_r_ro = 40,
                    mm_d_sleeve = 1.5,
                    split_ratio_rotor_salient = 0.2,
                    deg_alpha_rsp = 10,
                    pm = 2,
                    mm_r_ri = 5,
                    location = None
                    ):
        self.name = name
        self.color = color
        self.mm_r_ro = mm_r_ro
        self.mm_d_sleeve = mm_d_sleeve
        self.split_ratio_rotor_salient = split_ratio_rotor_salient
        self.deg_alpha_rsp = deg_alpha_rsp
        self.pm = pm                     # number of rotor modulator pole pairs
        self.mm_r_ri = mm_r_ri
        self.location = location         # move this part to another location other than origin (not supported yet)

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        mm_r_ro = self.mm_r_ro
        mm_d_sleeve = self.mm_d_sleeve
        split_ratio_rotor_salient = self.split_ratio_rotor_salient
        mm_d_rsp = split_ratio_rotor_salient * mm_r_ro
        # print('DEBUG: mm_d_rsp=', mm_d_rsp)
        deg_alpha_rsp = self.deg_alpha_rsp
        pm = self.pm

        mm_r_ri = self.mm_r_ri

        alpha_rp = 2*np.pi/pm # pole span
        deg_span_angle_rotor_pole_pair = 360/pm

        P0 = [mm_r_ro, 0]
        P1 = [mm_r_ro*cos(deg_alpha_rsp/180*np.pi/2), 
              mm_r_ro*-sin(deg_alpha_rsp/180*np.pi/2)]
        mm_r_groove = mm_r_ro * (1-split_ratio_rotor_salient)
        P3 = [  mm_r_groove *  cos(deg_span_angle_rotor_pole_pair/180*np.pi),
                mm_r_groove * -sin(deg_span_angle_rotor_pole_pair/180*np.pi) ]
        rad_temp = (mm_r_ro - mm_r_groove) * np.tan(8/180*np.pi) / mm_r_groove
        rad_P3_rotate_angle = deg_span_angle_rotor_pole_pair/180*np.pi - (deg_alpha_rsp / 2 / 180*np.pi + rad_temp)
        def iPark(P, theta):
            return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
        P2 = iPark(P3, rad_P3_rotate_angle)

        # P4 = iPark(P0, alpha_rp)

        P1_Mirror_CCW = P1_Mirror = P1[0], -P1[1]
        # print(P1_Mirror, P1)
        P2_Mirror_CCW = P2_Mirror = P2[0], -P2[1]
        P2_Rotate = iPark(P2, alpha_rp)

        angle_at_P2 = np.arctan2(P2[1], P2[0])

        P2_InnerEdge = [self.mm_r_ri* cos(angle_at_P2),
                        self.mm_r_ri* sin(angle_at_P2) ]
        P2_InnerEdge_Rotate =  iPark(P2_InnerEdge, alpha_rp)

        list_segments = []
        if bool_draw_whole_model:
            def draw_fraction(list_segments, P1, P2, P2_Rotate, P1_Mirror_CCW, P2_Mirror_CCW):
                # P1 = iPark(P1, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                # P2 = iPark(P2, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                # P2_Rotate = iPark(P2_Rotate, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                # P4 = iPark(P4, 45/2/180*np.pi) # rotate the rotor by an angle of 22.5 deg
                list_segments += drawer.drawArc([0,0], P2_Mirror_CCW, P2_Rotate)
                list_segments += drawer.drawLine(P2_Mirror_CCW, P1_Mirror_CCW)
                list_segments += drawer.drawArc([0,0], P1, P1_Mirror_CCW)
                list_segments += drawer.drawLine(P2, P1)
                # list_segments += drawer.drawLine(P1_CCW, P4)
            for i in range(pm):
                draw_fraction(list_segments, iPark(P1, i*alpha_rp), 
                                             iPark(P2, i*alpha_rp), 
                                             iPark(P2_Rotate, i*alpha_rp),
                                             iPark(P1_Mirror_CCW, i*alpha_rp),
                                             iPark(P2_Mirror_CCW, i*alpha_rp),
                                             )
                # if i==1:
                #     break
            # draw a circle (this is officially suggested by FEMM)
            PRI = [self.mm_r_ri, 0 ]
            list_segments += drawer.drawArc([0,0], PRI, [-PRI[0], PRI[1]])
            list_segments += drawer.drawArc([0,0],      [-PRI[0], PRI[1]], PRI)

            # innerCoord = ( 0.5*(P1[0]+P4[0]), 0.5*(P1[1]+P4[1]))

            # return [list_segments] # csToken # cross section token
            return {#'innerCoord': innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': None,
                    'inner_or_outer_region_to_remove': [True, False]
                    }
        else:
            list_segments += drawer.drawArc([0,0], P2_Mirror, P2_Rotate)
            list_segments += drawer.drawLine(P2_Mirror, P1_Mirror)
            list_segments += drawer.drawArc([0,0], P1, P1_Mirror)
            list_segments += drawer.drawLine(P2, P1)
            list_segments += drawer.drawLine(P2, P2_InnerEdge)
            list_segments += drawer.drawArc([0,0], P2_InnerEdge, P2_InnerEdge_Rotate)
            list_segments += drawer.drawLine(P2_InnerEdge_Rotate, P2_Rotate)

            # draw a circle (this is officially suggested by FEMM)
            PRI = [self.mm_r_ri, 0 ]
            # list_segments += drawer.drawArc([0,0], PRI, [-PRI[0], PRI[1]])
            # list_segments += drawer.drawArc([0,0],      [-PRI[0], PRI[1]], PRI)

            # innerCoord = ( 0.5*(P1[0]+P4[0]), 0.5*(P1[1]+P4[1]))

            # return [list_segments] # csToken # cross section token
            return {#'innerCoord': innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': None,
                    # 'inner_or_outer_region_to_remove': [True, False]
                    }
            raise Exception('Not supported. Please set bool_draw_whole_model to True.')


    @staticmethod
    def get_node_at_intersection(c,l):
        # this works for c and l having one or two intersections
        if c[0][0] != 0 or c[0][1] != 0:
            raise Exception('Not implemented for non-origin centered circle.')
        r = c[1]
        c = None
        x1, y1 = l[0][0], l[0][1]
        x2, y2 = l[1][0], l[1][1]
        if x1 == x2:
            raise Exception('Not implemented.')
        a = -(y2-y1)/(x2-x1)
        b = 1
        c = y1-(y2-y1)/(x2-x1)*x1
        Delta = sqrt(r**2*(a**2+b**2)-c**2)
        if Delta < 0:
            raise Exception('No intersection for given line and circle')
        x_solutions = (a*c + b*Delta)/(a**2+b**2), (a*c - b*Delta)/(a**2+b**2)
        y_solutions = (b*c - a*Delta)/(a**2+b**2), (b*c + a*Delta)/(a**2+b**2)
        return (x_solutions[0], y_solutions[0]), (x_solutions[1], y_solutions[1])
