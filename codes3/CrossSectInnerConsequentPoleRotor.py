from pylab import np, cos, sin
EPS = 1e-3 # [mm]

# class ExceptionBadDesign(Exception):
#     """Exception for notifying bad design."""
#     def __init__(self, message, payload=None):
#         self.message = message
#         self.payload = 'you could add more args here'
#     def __str__(self):
#         return str(self.message)

class CrossSectConsequentPoleRotor(object):
    # CrossSectInnerConsequentPoleRotor Describes the inner ConsequentPole rotor.
    #    Properties are set upon class creation and cannot be modified.
    #    The anchor point for this is the center of the rotor,
    #    with the x-axis directed along the center of one of the rotor poles
    def __init__(self, 
                    name = 'ConsequentPole Rotor',
                    color = '#FE840E',
                    mm_d_pm = 6,
                    deg_alpha_rm = 60,
                    # deg_alpha_rs = 10,
                    mm_d_ri = 8,
                    mm_r_ri = 40,
                    mm_d_rp = 5,
                    # mm_d_rs = 3,
                    p = 2, # Set pole-pairs to 2
                    s = 4, # Set magnet segments/pole to 4
                    location = None
                    ):
        self.name = name
        self.color = color

        self.mm_d_pm = mm_d_pm # depth of the permanent magnet

        self.deg_alpha_rm = deg_alpha_rm # angular span of the pole: class type DimAngular
        # self.deg_alpha_rs = deg_alpha_rs # segment span: class type DimAngular
        self.mm_d_ri = mm_d_ri           # rotor iron thickness: class type DimLinear
        self.mm_r_ri = mm_r_ri           # inner radius of rotor: class type DimLinear
        self.mm_d_rp = mm_d_rp           # interpolar iron thickness: class type DimLinear   
        # self.mm_d_rs = mm_d_rs           # inter segment iron thickness: class type DimLinear
        self.p = p                       # number of pole pairs
        self.s = s                       # number of segments  
        self.location = location         # move this part to another location other than origin (not supported yet)

        # Validate that magnet spans only one pole pitch  
        if self.deg_alpha_rm>(180/self.p):
            raise Exception('Invalid alpha_rm. Check that it is less than 180/p')

        # if self.s>1:
        #     # Validate that d_rs is non zero if there are segments  
        #     if self.mm_d_rs==0:
        #         raise Exception('Invalid d_rs. Check that it is positive for s>1')

        #     # Validate that segment span is legitimate
        #     if not (self.deg_alpha_rs<=self.deg_alpha_rm/self.s):
        #         raise Exception('Invalid deg_alpha_rs=%g. Check that it is less than alpha_rm/s=%g'%(self.deg_alpha_rs, self.deg_alpha_rm/self.s))
        # elif self.s==1:
        #     # Validate that alpha_rs and alpha_rm are set equal for s =1 
        #     if not (self.deg_alpha_rs==(self.deg_alpha_rm/self.s)):
        #         raise Exception('Invalid alpha_rs. Check that it is equal to alpha_rm for s=1', self.deg_alpha_rs, self.deg_alpha_rm, self.s)

        #     # Validate that d_rs is set zero for s=1
        #     if not (self.mm_d_rs==0):
        #         raise Exception('Invalid d_rs. Check that it is equal to 0 for s =1')

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        d_pm  = self.mm_d_pm
        alpha_rm = self.deg_alpha_rm * np.pi/180
        # alpha_rs = self.deg_alpha_rs * np.pi/180
        r_ri     = self.mm_r_ri
        d_ri     = self.mm_d_ri
        d_rp     = self.mm_d_rp
        # d_rs     = self.mm_d_rs
        p        = self.p
        s        = self.s
        alpha_rp = 2*np.pi/(2*p) # pole span    
        alpha_rs = 0
        # if abs(d_pm - d_rp) < 2*EPS: # d_pm is not defined
        #     print('Warn: [class CrossSectInnerNotchedRotor] Detect d_rp is too close to d_pm. To avoid small line entity error in JMAG, will set d_pm equal to d_rp in CrossSectInnerNotchedMagnet.') # d_pm is not defined here so we cannot set d_rp to d_pm.
        if abs(alpha_rp - alpha_rm) <= 2 * np.pi/180: # if alpha_rm and alpha_rp has a difference smaller than 2 deg, then let alpha_rm equal to alpha_rp.
            alpha_rm = alpha_rp
            if s == 1:
                alpha_rs = alpha_rm # alpha_rs is the variable actually being used in the following...
            else:
                raise 


        P1 = [r_ri, 0]

        r_P2 = r_ri + d_ri + d_pm
        # print('[CrossSectInnerNotchedRotor.py] DEBUG: ', r_P2, mm_r_ro)
        P2 = [r_P2, 0]


        # For consequent-pole iron span is equal to pole span

        alpha_P3 = alpha_rp   
        P3 = [r_P2*cos(alpha_P3), r_P2*sin(alpha_P3)]
        
        # print(alpha_rm)
        # print(alpha_rm)        
        # print(alpha_rm)
        # print(alpha_rm)
        # quit()        

        alpha_P4 = alpha_rp  
        r_P4 = r_ri + d_ri # = (r_P2 - d_rp) 
        P4 = [r_P4*cos(alpha_P4), r_P4*sin(alpha_P4)]
        print(alpha_P4)
        print(P4)
        
        r_P5 = r_P4
        P5 = [0, r_P5]
        print(P5)

        # alpha_P5 = alpha_P3 + alpha_rs # alpha_rs means rotor segment (of PM)
        # if abs(alpha_rs*s - alpha_rm)<EPS: # This means the inter-segment notch should span 0 deg, which mans the segmented design is reduced to a non-segmented design such that alpha_rm == alpha_rs*s
        #     alpha_P5 = alpha_P3 + alpha_rs*s
        #     s = 1 # this is a bad practice but it will help to re-use the drawing code of case s==1 below
        # P5 = [r_P4*cos(alpha_P5), r_P4*-sin(alpha_P5)]

        list_segments = []
        if s == 1:
            # No magnet segment!
            
            P6 = [0, r_ri]
            # P6 = [r_ri*cos(alpha_P5), r_ri*-sin(alpha_P5)]

            if alpha_rm >= alpha_rp*0.9800:
                print('[CrossSectInnerNotchedRotor.py] Non-NOTCHED ROTOR IS USED.\n')
                print('[CrossSectInnerNotchedRotor.py] alpha_P5 is', alpha_P5, alpha_P5/np.pi*180)
                P1p5 = [P2[0] - d_rp, P2[1]]
                if bool_draw_whole_model:
                    list_segments += drawer.drawArc([0,0], P1, [-P1[0], P1[1]])
                    list_segments += drawer.drawArc([0,0], [-P1[0], P1[1]], P1)
                    list_segments += drawer.drawArc([0,0], P1p5, [-P1p5[0], P1p5[1]])
                    list_segments += drawer.drawArc([0,0], [-P1p5[0], P1p5[1]], P1p5)
                    
                    # list_segments += drawer.drawLine([5, 10], [0, 0])
                    # list_segments += drawer.drawLine([5, 10], [-10, 20])
                    # list_segments += drawer.drawLine([5, 10], [-10, 200])
                    # print(P1, P2, P3, P4, P5, P6)
                    # print(P1, P2, P3, P4, P5, P6)
                    # print(P1, P2, P3, P4, P5, P6)

                    # useless

                    # raise
                else:
                    list_segments += drawer.drawLine(P1, P1p5)
                    list_segments += drawer.drawArc([0,0], P5, P1p5)
                    list_segments += drawer.drawLine(P5, P6)
                    list_segments += drawer.drawArc([0,0], P6, P1)
                    # list_segments += drawer.drawLine([5, 10], [0, 0])
                    # list_segments += drawer.drawLine([5, 10], [-10, 20])
                    # list_segments += drawer.drawLine([5, 10], [-10, 200])
                    # print(P1, P2, P3, P4, P5, P6)
                    # print(P1, P2, P3, P4, P5, P6)
                    # print(P1, P2, P3, P4, P5, P6)

                    # useless

            else:
                print(F'{alpha_rm=}, {alpha_rp=}, {s=}')
                print(F'{alpha_rm=}, {alpha_rp=}, {s=}')
                print(F'{alpha_rm=}, {alpha_rp=}, {s=}')
                if bool_draw_whole_model:
                    def iPark(P, theta):
                        return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
                    

                    def draw_fraction(list_segments, P1, P2, P3, P4, P5, P6):
                        if (i % 2) == 0:
                            list_segments += drawer.drawLine(P1, P2)
                            list_segments += drawer.drawArc([0,0], P2, P3)
                        else :
                            list_segments += drawer.drawLine(P3, P4)
                            list_segments += drawer.drawArc([0,0], P4, P5)
                            list_segments += drawer.drawLine(P5, P6)
                            list_segments += drawer.drawArc([0,0], P1, P6)
######################### 非常有用！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！############################
                    for i in range(2*p):
                        draw_fraction(list_segments, P1, P2, P3, P4, P5, P6)
                        # break





                    # # draw a circle (this is officially suggested by FEMM)
                    # list_segments += drawer.drawArc([0,0], P1, [-P1[0], P1[1]])
                    # list_segments += drawer.drawArc([0,0],     [-P1[0], P1[1]], P1)
######################### 非常有用！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！############################

        innerCoord = ( 0.5*(P1[0]+P4[0]), 0.5*(P1[1]+P4[1]))

        # return [list_segments] # csToken # cross section token
        return {'innerCoord': innerCoord, 'list_regions':[list_segments], 'mirrorAxis': None,}
                # 'list_regions_to_remove': }

class CrossSectConsequentPoleMagnet(object):
    def __init__(self, 
                    name = 'ConsequentPole',
                    color = '#0BA0E2',
                    ConsequentPole_rotor = None,
                    ):
        self.name = name
        self.color = color
        self.ConsequentPole_rotor = ConsequentPole_rotor

    def draw(self, drawer, bool_re_evaluate= False, bool_draw_whole_model= False):

        if False == bool_re_evaluate:
            drawer.getSketch(self.name, self.color)

        d_pm     = self.ConsequentPole_rotor.mm_d_pm
        alpha_rm = self.ConsequentPole_rotor.deg_alpha_rm * np.pi/180
        # alpha_rs = self.ConsequentPole_rotor.deg_alpha_rs * np.pi/180
        r_ri     = self.ConsequentPole_rotor.mm_r_ri
        d_ri     = self.ConsequentPole_rotor.mm_d_ri
        d_rp     = self.ConsequentPole_rotor.mm_d_rp
        # d_rs     = self.ConsequentPole_rotor.mm_d_rs
        p        = self.ConsequentPole_rotor.p
        s        = self.ConsequentPole_rotor.s
        alpha_rp = 2*np.pi/(2*p) # pole span
        alpha_rs = 0


        # rotor inter-pole notch being too small
        if alpha_rm >= alpha_rp*0.9800:
            print('[CrossSectInnerConsequentPoleRotor.py] FULL POLE PITCH MAGNET IS USED.')
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
            print('[Warn] [class CrossSectInnerConsequentPoleMagnet] Magnet is fully spanned.')


        # print(alpha_rm/np.pi*180)
        # print(alpha_rm/np.pi*180)
        # print(alpha_rm/np.pi*180)

        if d_pm + 2*EPS < d_rp:
            print('[Warn]: [class CrossSectInnerNotchedMagnet] Detect d_rp is too close to d_pm. To avoid small line entity error in JMAG, set d_pm equal to d_rp because rotor core is plotted already.')
            raise ExceptionBadDesign('[Error] Magnet depth d_pm is too close to inter-pole notch depth d_rp.')
        


        
        P1 = [r_ri, 0]

        r_P2 = r_ri + d_ri + d_pm
        P2 = [r_P2, 0]
        P2_extra = [0, r_P2]

        alpha_P3 = alpha_rp           # For consequent-pole iron span is equal to pole span
        P3 = [r_P2*cos(alpha_P3), r_P2*sin(alpha_P3)]

        alpha_P4 = alpha_rp  
        r_P4 = r_ri + d_ri # = (r_P2 - d_rp) 
        P4 = [r_P4*cos(alpha_P4), r_P4*sin(alpha_P4)]
        # print(alpha_P4)
        # print(P4)
        P4_extra = [r_P4, 0]

        P3_extra = [(r_P4+d_pm)*cos(alpha_P3), (r_P4+d_pm)*-sin(alpha_P3)]

        alpha_P5 = alpha_P3 + alpha_rs # alpha_rs means rotor segment (of PM)
        P5 = [r_P4*cos(alpha_P5), r_P4*-sin(alpha_P5)]
        r_P5 = r_P4
        P5 = [0, r_P5]
        # print(P5)

        list_regions = []
        list_segments = []
        if True:

            if s>1:
                alpha_notch  = (alpha_rm - s*alpha_rs) / (s-1) # 永磁体占的弧度
            P6_extra = [(r_P4+d_pm)*cos(alpha_P5), (r_P4+d_pm)*-sin(alpha_P5)]

            Rout = r_P4+d_pm
            Rin  = r_P4
            self.mm2_magnet_area = alpha_rm/alpha_rp  *  np.pi*(Rout**2 - Rin**2) # magnet area for all the poles
            print('[CrossSectInnerConsequentPoleRotor.py] Magnet area in total is %g mm^2'%(self.mm2_magnet_area))
            if bool_re_evaluate:
                return self.mm2_magnet_area


                    ################## 只有我有用 ###############################
                    ################## 只有我有用 ###############################

            if bool_draw_whole_model:
                def iPark(P, theta):
                    return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
                def draw_fraction(list_segments, P3, P4, P5, P2_extra):
                   
                    if (i % 2) == 1:
                        list_segments += drawer.drawLine(P3, P4)
                        list_segments += drawer.drawArc([0,0], P4, P5)
                        list_segments += drawer.drawArc([0,0], P3, P2_extra)
                        list_segments += drawer.drawLine(P5, P2_extra)

                for i in range(2*p):
                    draw_fraction(list_segments, P3, 
                                                 P4, 
                                                 P5, 
                                                 P2_extra)

                #     list_segments += drawer.drawLine(P2, P4_extra)
                #     if (i % 2) == 1:
                #         P2_CCW = iPark(P2, alpha_rp)
                #         list_segments += drawer.drawArc([0,0], P2, P2_CCW)
                #         P4_CCW = iPark(P4_extra, alpha_rp)
                #         list_segments += drawer.drawArc([0,0], P4_extra, P4_CCW)

                # for i in range(2*p):
                #     draw_fraction(list_segments, iPark(P2, i*alpha_rp), 
                #                                  iPark(P3_extra, i*alpha_rp), 
                #                                  iPark(P4_extra, i*alpha_rp), 
                #                                  iPark(P6_extra, i*alpha_rp))


                    ################## 只有我有用 ###############################
                    ################## 只有我有用 ###############################
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


####################### Shaft ########################
class CrossSectConsequentPoleShaft(object):

    def __init__(self, 
                    name = 'Shaft',
                    color = '#0EE0E2',
                    ConsequentPole_rotor = None,
                    ):
        self.name = name
        self.color = color
        self.ConsequentPole_rotor = ConsequentPole_rotor

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        d_pm     = self.ConsequentPole_rotor.mm_d_pm
        alpha_rm = self.ConsequentPole_rotor.deg_alpha_rm * np.pi/180
        # alpha_rs = self.notched_rotor.deg_alpha_rs * np.pi/180
        r_ri     = self.ConsequentPole_rotor.mm_r_ri
        d_ri     = self.ConsequentPole_rotor.mm_d_ri
        d_rp     = self.ConsequentPole_rotor.mm_d_rp
        # d_rs     = self.notched_rotor.mm_d_rs
        p        = self.ConsequentPole_rotor.p
        s        = self.ConsequentPole_rotor.s
        alpha_rp = 2*np.pi/(2*p) # pole span

        P1 = [r_ri, 0]
        NP1 = [-r_ri, 0]

        list_regions = []
        list_segments = []
        if True:
            list_segments += drawer.drawArc([0,0], NP1, P1)
            list_segments += drawer.drawArc([0,0], P1, NP1)

            # list_segments += drawer.drawLine([5, 10], [0, 0])
            # list_segments += drawer.drawLine([5, 10], [-10, 20])
            # list_segments += drawer.drawLine([5, 10], [-10, 200])
            # print(P1)
            # print(P1)
            # print(P1)


            list_regions.append(list_segments)
            list_segments = []

        innerCoord = ( 0, 0 )

        # return list_regions # csToken # cross section token
        return {'innerCoord': innerCoord, 'list_regions':list_regions, 'mirrorAxis': None}


### testrun

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



