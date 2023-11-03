from pylab import np, cos, sin
EPS = 1e-3 # [mm]

class ExceptionBadDesign(Exception):
    """Exception for notifying bad design."""
    def __init__(self, message, payload=None):
        self.message = message
        self.payload = 'you could add more args here'
    def __str__(self):
        return str(self.message)

class CrossSectInnerNotchedRotor(object):
    # CrossSectInnerNotchedRotor Describes the inner notched rotor.
    #    Properties are set upon class creation and cannot be modified.
    #    The anchor point for this is the center of the rotor,
    #    with the x-axis directed along the center of one of the rotor poles
    def __init__(self, 
                    name = 'Notched Rotor',
                    color = '#FE840E',
                    mm_d_pm = 6,
                    deg_alpha_rm = 60,
                    deg_alpha_rs = 10,
                    mm_d_ri = 8,
                    mm_r_ri = 40,
                    mm_d_rp = 5,
                    mm_d_rs = 3,
                    p = 2, # Set pole-pairs to 2
                    s = 4, # Set magnet segments/pole to 4
                    location = None
                    ):
        self.name = name
        self.color = color
        self.mm_d_pm = mm_d_pm # depth of the permanent magnet
        self.deg_alpha_rm = deg_alpha_rm # angular span of the pole: class type DimAngular
        self.deg_alpha_rs = deg_alpha_rs # segment span: class type DimAngular
        self.mm_d_ri = mm_d_ri           # rotor iron thickness: class type DimLinear
        self.mm_r_ri = mm_r_ri           # inner radius of rotor: class type DimLinear
        self.mm_d_rp = mm_d_rp           # interpolar iron thickness: class type DimLinear
        self.mm_d_rs = mm_d_rs           # inter segment iron thickness: class type DimLinear
        self.p = p                       # number of pole pairs
        self.s = s                       # number of segments  
        self.location = location         # move this part to another location other than origin (not supported yet)

        # Validate that magnet spans only one pole pitch  
        if self.deg_alpha_rm>(180/self.p):
            raise Exception('Invalid alpha_rm. Check that it is less than 180/p')

        if self.s>1:
            # Validate that d_rs is non zero if there are segments  
            if self.mm_d_rs==0:
                raise Exception('Invalid d_rs. Check that it is positive for s>1')

            # Validate that segment span is legitimate
            if not (self.deg_alpha_rs<=self.deg_alpha_rm/self.s):
                raise Exception('Invalid deg_alpha_rs=%g. Check that it is less than alpha_rm/s=%g'%(self.deg_alpha_rs, self.deg_alpha_rm/self.s))
        elif self.s==1:
            # Validate that alpha_rs and alpha_rm are set equal for s =1 
            if not (self.deg_alpha_rs==(self.deg_alpha_rm/self.s)):
                raise Exception('Invalid alpha_rs. Check that it is equal to alpha_rm for s=1', self.deg_alpha_rs, self.deg_alpha_rm, self.s)

            # Validate that d_rs is set zero for s=1
            if not (self.mm_d_rs==0):
                raise Exception('Invalid d_rs. Check that it is equal to 0 for s =1')

    def draw(self, drawer, bool_draw_whole_model=False):

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

        # if abs(d_pm - d_rp) < 2*EPS: # d_pm is not defined
        #     print('Warn: [class CrossSectInnerNotchedRotor] Detect d_rp is too close to d_pm. To avoid small line entity error in JMAG, will set d_pm equal to d_rp in CrossSectInnerNotchedMagnet.') # d_pm is not defined here so we cannot set d_rp to d_pm.
        if abs(alpha_rp - alpha_rm) <= 2 * np.pi/180: # if alpha_rm and alpha_rp has a difference smaller than 2 deg, then let alpha_rm equal to alpha_rp.
            alpha_rm = alpha_rp
            if s == 1:
                alpha_rs = alpha_rm # alpha_rs is the variable actually being used in the following...
            else:
                print('DEBUG s>1 notched rotor')
            #     print('s=%d: This is not tested. For now it simply assumes the iron notch between poles becomes the iron notch between the segments of one pole.' % (s))
            # print('[class CrossSectInnerNotchedRotor] Rotor has no notch, i.e., there is no P2 or P3.')
            # print('alpha_rp is', alpha_rp)
            # print('alpha_rm is', alpha_rm)
            # print('alpha_rs is', alpha_rs)

        P1 = [r_ri, 0]

        r_P2 = r_ri + d_ri + d_rp
        # print('[CrossSectInnerNotchedRotor.py] DEBUG: ', r_P2, mm_r_ro)
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
        # debug

        # list_segments += drawer.drawLine([5, 10], [0, 0])
        # list_segments += drawer.drawLine([5, 10], [-10, 20])
        # list_segments += drawer.drawLine([5, 10], [-10, 200])
        # print(P1, P2, P3, P4, P5)
        # print(P1, P2, P3, P4, P5)
        # print(P1, P2, P3, P4, P5)



        def print_point(P):
            print( '(%g, %g)' % (P[0], P[1]) )
            print_point(P1)
            print_point(P2)
            print_point(P3)
            print_point(P4)
            print_point(P5)
            print_point(P6)
            quit()

        list_segments = []
        if s == 1:
            # No magnet sement!
            # Then P6 is an extra point for a full rotor
            P6 = [r_ri*cos(alpha_P5), r_ri*-sin(alpha_P5)]
            # list_segments += drawer.drawLine([5, 10], [0, 0])
            # list_segments += drawer.drawLine([5, 10], [-10, 20])
            # list_segments += drawer.drawLine([5, 10], [-10, 200])
            # print(P1, P2, P3, P4, P5, P6)
            # print(P1, P2, P3, P4, P5, P6)
            # print(P1, P2, P3, P4, P5, P6)
            if alpha_rm >= alpha_rp*0.9800:
                print('[CrossSectInnerNotchedRotor.py] Non-NOTCHED ROTOR IS USED.\n')
                print('[CrossSectInnerNotchedRotor.py] alpha_P5 is', alpha_P5, alpha_P5/np.pi*180)
                P1p5 = [P2[0] - d_rp, P2[1]]
                if bool_draw_whole_model:
                    list_segments += drawer.drawArc([0,0], P1, [-P1[0], P1[1]])
                    list_segments += drawer.drawArc([0,0], [-P1[0], P1[1]], P1)
                    list_segments += drawer.drawArc([0,0], P1p5, [-P1p5[0], P1p5[1]])
                    list_segments += drawer.drawArc([0,0], [-P1p5[0], P1p5[1]], P1p5)
                    # raise
                else:
                    list_segments += drawer.drawLine(P1, P1p5)
                    list_segments += drawer.drawArc([0,0], P5, P1p5)
                    list_segments += drawer.drawLine(P5, P6)
                    list_segments += drawer.drawArc([0,0], P6, P1)
            else:
                if bool_draw_whole_model:
                    def iPark(P, theta):
                        return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
                    def draw_fraction(list_segments, P2, P3, P4, P5):
                        list_segments += drawer.drawArc([0,0], P3, P2)
                        list_segments += drawer.drawLine(P3, P4)
                        list_segments += drawer.drawArc([0,0], P5, P4)
                        P5_CCW = iPark(P5, alpha_rp)
                        list_segments += drawer.drawLine(P5_CCW, P2)
                    for i in range(2*p):
                        draw_fraction(list_segments, iPark(P2, i*alpha_rp), iPark(P3, i*alpha_rp), iPark(P4, i*alpha_rp), iPark(P5, i*alpha_rp))
                    # draw a circle (this is officially suggested by FEMM)
                    list_segments += drawer.drawArc([0,0], P1, [-P1[0], P1[1]])
                    list_segments += drawer.drawArc([0,0],     [-P1[0], P1[1]], P1)
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
            if bool_draw_whole_model == True:
                raise Exception('NOT IMPLEMENTED for s>1.')

            alpha_notch  = (alpha_rm - s*alpha_rs) / (s-1) # inter-segment notch占的弧度
            P5

            r_P6 = r_ri + d_ri + d_rs
            P6 = [r_P6*cos(alpha_P5), r_P6*-sin(alpha_P5)]

            alpha_P7 = alpha_P5 + alpha_notch
            P7 = [r_P6*cos(alpha_P7), r_P6*-sin(alpha_P7)]

            P8 = [r_P4*cos(alpha_P7), r_P4*-sin(alpha_P7)]

            if alpha_rm >= alpha_rp*0.9800: # no inter-pole notch
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

        innerCoord = ( 0.5*(P1[0]+P4[0]), 0.5*(P1[1]+P4[1]))

        # return [list_segments] # csToken # cross section token
        return {'innerCoord': innerCoord, 'list_regions':[list_segments], 'mirrorAxis': None}

class CrossSectInnerNotchedMagnet(object):
    def __init__(self, 
                    name = 'Notched Rotor',
                    color = '#0BA0E2',
                    notched_rotor = None,
                    ):
        self.name = name
        self.color = color
        self.notched_rotor = notched_rotor

    def draw(self, drawer, bool_re_evaluate=False, bool_draw_whole_model=False):

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

        # rotor inter-pole notch being too small
        if alpha_rm >= alpha_rp*0.9800:
            print('[CrossSectInnerNotchedRotor.py] FULL POLE PITCH MAGNET IS USED.')
            alpha_rm = alpha_rp

        # print(alpha_rm/np.pi*180)
        # print(alpha_rm/np.pi*180)
        # print(alpha_rm/np.pi*180)

        if abs(alpha_rp - alpha_rm) <= 2 * np.pi/180: # if alpha_rm and alpha_rp has a difference smaller than 2 deg, then let alpha_rm equal to alpha_rp.
            alpha_rm = alpha_rp
            if s == 1:
                alpha_rs = alpha_rm # alpha_rs is the variable actually being used in the following...
            else:
                print('[Warn] s=%d: This is not tested. For now it simply assumes the iron notch between poles becomes the iron notch between the segments of one pole.' % (s))
            print('[Warn] [class CrossSectInnerNotchedMagnet] Magnet is fully spanned.')


        # print(alpha_rm/np.pi*180)
        # print(alpha_rm/np.pi*180)
        # print(alpha_rm/np.pi*180)

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
            print('[CrossSectInnerNotchedRotor.py] Magnet area in total is %g mm^2'%(self.mm2_magnet_area))
            if bool_re_evaluate:
                return self.mm2_magnet_area

            if bool_draw_whole_model:
                def iPark(P, theta):
                    return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
                def draw_fraction(list_segments, P3_extra, P4, P5, P6_extra):
                    list_segments += drawer.drawLine(P3_extra, P4)
                    list_segments += drawer.drawArc([0,0], P5, P4)
                    list_segments += drawer.drawLine(P5, P6_extra)
                    list_segments += drawer.drawArc([0,0], P6_extra, P3_extra)
                for i in range(2*p):
                    draw_fraction(list_segments, iPark(P3_extra, i*alpha_rp), 
                                                 iPark(P4, i*alpha_rp), 
                                                 iPark(P5, i*alpha_rp), 
                                                 iPark(P6_extra, i*alpha_rp))
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

class CrossSectSleeve(object):
    def __init__(self, 
                    name = 'Sleeve',
                    color = '#11E322',
                    notched_magnet = None,
                    d_sleeve = None
                    ):
        self.name = name
        self.color = color
        self.notched_magnet = notched_magnet

        self.d_sleeve = d_sleeve

    def draw(self, drawer):

        drawer.getSketch(self.name, self.color)

        r_ri  = self.notched_magnet.notched_rotor.mm_r_ri
        d_ri  = self.notched_magnet.notched_rotor.mm_d_ri
        d_pm  = self.notched_magnet.notched_rotor.mm_d_pm
        p     = self.notched_magnet.notched_rotor.p

        r_or = r_ri + d_ri + d_pm 
        d_sleeve = self.d_sleeve

        P1 = [r_or, 0]
        P2 = [r_or+d_sleeve, 0]

        P3 = [  cos(np.pi/p)*P1[0] + sin(np.pi/p)*P1[1],
               -sin(np.pi/p)*P1[0] + cos(np.pi/p)*P1[1] ]        
        P4 = [  cos(np.pi/p)*P2[0] + sin(np.pi/p)*P2[1],
               -sin(np.pi/p)*P2[0] + cos(np.pi/p)*P2[1] ]        

        list_regions = []
        list_segments = []
        list_segments += drawer.drawLine(P1, P2)
        list_segments += drawer.drawArc([0,0], P4, P2)
        list_segments += drawer.drawLine(P4, P3)
        list_segments += drawer.drawArc([0,0], P3, P1)

        list_regions.append(list_segments)
        list_segments = []
        # raise KeyboardInterrupt

        # this is on the edge
        innerCoord = ( 0.5*(P1[0]+P2[0]), 0.5*(P1[1]+P2[1]))

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

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        # d_pm     = self.notched_rotor.mm_d_pm
        # alpha_rm = self.notched_rotor.deg_alpha_rm * np.pi/180
        # alpha_rs = self.notched_rotor.deg_alpha_rs * np.pi/180
        r_ri     = self.notched_rotor.mm_r_ri
        # d_ri     = self.notched_rotor.mm_d_ri
        # d_rp     = self.notched_rotor.mm_d_rp
        # d_rs     = self.notched_rotor.mm_d_rs
        # p        = self.notched_rotor.p
        # s        = self.notched_rotor.s
        # alpha_rp = 2*np.pi/(2*p) # pole span

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
    if True:
        from utility import my_execfile
        my_execfile('./default_setting.py', g=globals(), l=locals())
        fea_config_dict

        toolJd = JMAG.JMAG(fea_config_dict)

        project_name          = 'proj%d'%(0)
        expected_project_file_path = './' + "%s.jproj"%(project_name)
        toolJd.open(expected_project_file_path)

    if True:
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
    toolJd.iRotateCopy = notched_rotor.p*2
    region1 = toolJd.prepareSection(list_regions)
    
    if True:
        notched_magnet = CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
                                                      color = '#0E001E',
                                                      notched_rotor = notched_rotor
                                                    )

    list_regions = notched_magnet.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = notched_rotor.p*2
    region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')

    # Pre-process
    # toolJd.preProcess(makeToken)
    # model.CloseCadLink()



