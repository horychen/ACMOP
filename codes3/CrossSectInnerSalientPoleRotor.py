from pylab import np, cos, sin
EPS = 1e-3 # [mm]

class ExceptionBadDesign(Exception):
    """Exception for notifying bad design."""
    def __init__(self, message, payload=None):
        self.message = message
        self.payload = 'you could add more args here'
    def __str__(self):
        return str(self.message)

class CrossSectInnerSalientPoleRotor(object):
    def __init__(self, 
                    name = 'SalientPoleRotor',
                    color = '#FE840E',
                    mm_r_ro = 40,
                    mm_d_sleeve = 1.5,
                    split_ratio_rotor_salient = 0.2,
                    deg_alpha_rsp = 10,
                    pm = 2,
                    mm_r_ri = 10,
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
