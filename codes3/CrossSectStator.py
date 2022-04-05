from pylab import np, cos, sin, arctan
class CrossSectInnerRotorStator:
    # CrossSectInnerRotorStator Describes the inner rotor motor stator.
    #    Properties are set upon class creation and cannot be modified.
    #    The anchor point for this is the center of the stator,
    #    with the x-axis directed down the center of one of the stator teeth.
    def __init__(self,
                    name  = 'StatorCore',
                    color = '#BAFD01',
                    deg_alpha_st = 40, # span angle of tooth: class type DimAngular
                    deg_alpha_sto = 20, # angle of tooth edge: class type DimAngular
                    mm_r_si      = 40, # inner radius of stator teeth: class type DimLinear
                    mm_d_sto      = 5,  # tooth edge length: class type DimLinear
                    mm_d_stt      = 10, # tooth tip length: class type DimLinear
                    mm_d_st      = 15, # tooth base length: class type DimLinear
                    mm_d_sy      = 15, # back iron thickness: class type DimLinear
                    mm_w_st      = 13, # tooth base width: class type DimLinear
                    mm_r_st      = 0,  # fillet on outter tooth: class type DimLinear
                    mm_r_sf      = 0,  # fillet between tooth tip and base: class type DimLinear
                    mm_r_sb      = 0,  # fillet at tooth base: class type DimLinear
                    Q            = 6,  # number of stator slots (integer)
                    location = None
                ):

        self.name = name
        self.color = color
        self.deg_alpha_st = deg_alpha_st
        self.deg_alpha_sto = deg_alpha_sto
        self.mm_r_si      = mm_r_si     
        self.mm_d_sto      = mm_d_sto     
        self.mm_d_stt      = mm_d_stt     
        self.mm_d_st      = mm_d_st     
        self.mm_d_sy      = mm_d_sy     
        self.mm_w_st      = mm_w_st     
        self.mm_r_st      = mm_r_st     
        self.mm_r_sf      = mm_r_sf     
        self.mm_r_sb      = mm_r_sb     
        self.Q = Q               
        self.location = location 

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        alpha_st = self.deg_alpha_st * np.pi/180
        alpha_so = -self.deg_alpha_sto * np.pi/180
        r_si = self.mm_r_si
        d_so = self.mm_d_sto
        d_sp = self.mm_d_stt
        d_st = self.mm_d_st
        d_sy = self.mm_d_sy
        w_st = self.mm_w_st
        r_st = self.mm_r_st
        r_sf = self.mm_r_sf
        r_sb = self.mm_r_sb
        Q    = self.Q

        alpha_slot_span = 360/Q * np.pi/180

        P1 = [r_si, 0]
        P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]
        P3_temp = [ d_so*cos(alpha_st*0.5), 
                    d_so*-sin(alpha_st*0.5)]
        P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
                             -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        P3 = [  P3_local_rotate[0] + P2[0],
                P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_sp
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度)]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        # Option 1
        # 6 [94.01649113009418, -25.191642873516543] 97.33303383272235
        # Radius_InnerStatorYoke = r_si+d_sp+d_st
        # print('Radius_InnerStatorYoke #1:', Radius_InnerStatorYoke)
        # Option 2
        Radius_InnerStatorYoke = np.sqrt(P5[0]**2 + P5[1]**2)
        # print('Radius_InnerStatorYoke #2:', Radius_InnerStatorYoke)
        P6 = [ Radius_InnerStatorYoke *  cos(alpha_slot_span*0.5),
               Radius_InnerStatorYoke * -sin(alpha_slot_span*0.5) ]

        P7 = [ (r_si+d_sp+d_st+d_sy)*cos(alpha_slot_span*0.5),
               (r_si+d_sp+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
        P8 = [  r_si+d_sp+d_st+d_sy, 0]

        list_segments = []
        if bool_draw_whole_model:
            P2_Mirror = [P2[0], -P2[1]] # = iPark(P2, alpha_st)
            P3_Mirror = [P3[0], -P3[1]]
            P4_Mirror = [P4[0], -P4[1]]
            P5_Mirror = [P5[0], -P5[1]]
            def iPark(P, theta):
                return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
            def draw_fraction(list_segments, P2, P3, P4, P5,
                                            P2_Mirror, P3_Mirror, P4_Mirror, P5_Mirror):
                P5_Rotate = iPark(P5, alpha_slot_span)
                list_segments += drawer.drawArc([0,0], P2, P2_Mirror)
                list_segments += drawer.drawLine(P2, P3)
                list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
                list_segments += drawer.drawLine(P3, P4)
                list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
                list_segments += drawer.drawLine(P4, P5)
                list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)
                list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate)
            for i in range(Q):
                draw_fraction(list_segments, iPark(P2, i*alpha_slot_span), 
                                             iPark(P3, i*alpha_slot_span), 
                                             iPark(P4, i*alpha_slot_span), 
                                             iPark(P5, i*alpha_slot_span),
                                             iPark(P2_Mirror, i*alpha_slot_span), 
                                             iPark(P3_Mirror, i*alpha_slot_span), 
                                             iPark(P4_Mirror, i*alpha_slot_span), 
                                             iPark(P5_Mirror, i*alpha_slot_span), )
                # raise
            # draw a circle (this is officially suggested)
            list_segments += drawer.drawArc([0,0], P8, [-P8[0], P8[1]])
            list_segments += drawer.drawArc([0,0],     [-P8[0], P8[1]], P8)
        else:
            list_segments += drawer.drawArc([0,0], P2, P1)
            list_segments += drawer.drawLine(P2, P3)
            list_segments += drawer.drawLine(P3, P4)
            list_segments += drawer.drawLine(P4, P5)
            list_segments += drawer.drawArc([0,0], P6, P5)
            list_segments += drawer.drawLine(P6, P7)
            list_segments += drawer.drawArc([0,0], P7, P8)
            list_segments += drawer.drawLine(P8, P1)

        # DEBUG
        # for ind, point in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
        #     print(ind+1, point, np.sqrt(point[0]**2+point[1]**2))

        self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

        # return [list_segments]
        return {'innerCoord': self.innerCoord, 'list_regions':[list_segments], 'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])]}

def get_area_polygon(a,b,c,d):
    x1, x2, x3, x4 = a[0], b[0], c[0], d[0]
    y1, y2, y3, y4 = a[1], b[1], c[1], d[1]

    return 0.5*abs( x1*y2-y1*x2 + x2*y3-y2*x3 + x3*y4-y3*x4 + x4*y1-y4*x1 )

class CrossSectInnerRotorStatorWinding(object):
    def __init__(self, 
                    name = 'Coils',
                    color = '#3D9970',
                    stator_core = None
                    ):
        self.type = 'Cynlinder'
        self.name = name
        self.color = color
        self.stator_core = stator_core

    def draw(self, drawer, bool_re_evaluate=False, bool_draw_whole_model=False):

        if False == bool_re_evaluate:
            drawer.getSketch(self.name, self.color)

        alpha_st = self.stator_core.deg_alpha_st * np.pi/180
        alpha_so = self.stator_core.deg_alpha_sto * np.pi/180
        r_si     = self.stator_core.mm_r_si
        d_so     = self.stator_core.mm_d_sto
        d_sp     = self.stator_core.mm_d_stt
        d_st     = self.stator_core.mm_d_st
        d_sy     = self.stator_core.mm_d_sy
        w_st     = self.stator_core.mm_w_st
        r_st     = self.stator_core.mm_r_st
        r_sf     = self.stator_core.mm_r_sf
        r_sb     = self.stator_core.mm_r_sb
        Q        = self.stator_core.Q

        alpha_slot_span = 360/Q * np.pi/180

        P1 = [r_si, 0]

            # 乘以0.99或0.95避免上层导体和下层导体重合导致导入Designer时产生多余的Parts。
        POpen = [(r_si+d_sp)*cos(alpha_slot_span*0.5*1.00), (r_si+d_sp)*-sin(alpha_slot_span*0.5*1.00)]

        # P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]

        # P3_temp = [ d_so*cos(alpha_st*0.5), 
        #             d_so*-sin(alpha_st*0.5)]
        # P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
        #                      -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        # P3 = [  P3_local_rotate[0] + P2[0],
        #         P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_sp
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度) ]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        PMiddle45 = [0.5*(P4[0] + P5[0]), P4[1]]
        TheRadius = (P5[0] - P4[0])*0.45

            # 为了使得槽和导体之间不要接触，试着添加5%的clearance？
        P6 = [ (r_si+d_sp+d_st)*cos(alpha_slot_span*0.5) *1.00,
               (r_si+d_sp+d_st)*-sin(alpha_slot_span*0.5) *1.00 ]

        self.mm2_slot_area = 2 * get_area_polygon(P4, P5, P6, POpen)
        print('[CrossSectStator.py] Stator slot area is %g mm^2'%(self.mm2_slot_area))
        if bool_re_evaluate:
            return self.mm2_slot_area

        PMiddle6Open = [ 0.5*(P6[0]+POpen[0]), 0.5*(P6[1]+POpen[1])]
        self.PCoil = PCoil = [ 0.5*(PMiddle45[0]+PMiddle6Open[0]), 0.5*(PMiddle45[1]+PMiddle6Open[1])]

        # P7 = [ (r_si+d_sp+d_st+d_sy)*cos(alpha_slot_span*0.5),
        #        (r_si+d_sp+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
        # P8 = [  r_si+d_sp+d_st+d_sy, 0]

        # Compute the vector starting from PCoil to one of the corner of the polygon.
        def shrink(PC, P):
            vector = [ P[0] - PC[0], P[1] - PC[1]]
            return [ PC[0]+0.900*vector[0], PC[1]+0.900*vector[1] ]
        P6_Shrink = shrink(PCoil, P6)
        P5_Shrink = shrink(PCoil, P5)
        P4_Shrink = shrink(PCoil, P4)
        POpen_Shrink = shrink(PCoil, POpen)

        list_regions = []

        list_segments = []
        if bool_draw_whole_model:
            P6_Shrink_Mirror = [P6_Shrink[0], -P6_Shrink[1]]
            P5_Shrink_Mirror = [P5_Shrink[0], -P5_Shrink[1]]
            P4_Shrink_Mirror = [P4_Shrink[0], -P4_Shrink[1]]
            POpen_Shrink_Mirror = [POpen_Shrink[0], -POpen_Shrink[1]]
            def iPark(P, theta):
                return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
            def draw_fraction(list_segments, P6_Shrink, P5_Shrink, P4_Shrink, POpen_Shrink,
                                            P6_Shrink_Mirror, P5_Shrink_Mirror, P4_Shrink_Mirror,  POpen_Shrink_Mirror):
                list_segments += drawer.drawArc([0,0], P6_Shrink, P5_Shrink)
                list_segments += drawer.drawLine(P5_Shrink, P4_Shrink)
                list_segments += drawer.drawLine(P4_Shrink, POpen_Shrink)
                list_segments += drawer.drawLine(POpen_Shrink, P6_Shrink)
                list_regions.append(list_segments)
                list_segments = []

                list_segments += drawer.drawArc([0,0], P5_Shrink_Mirror, P6_Shrink_Mirror)
                list_segments += drawer.drawLine(P5_Shrink_Mirror, P4_Shrink_Mirror)
                list_segments += drawer.drawLine(P4_Shrink_Mirror, POpen_Shrink_Mirror)
                list_segments += drawer.drawLine(POpen_Shrink_Mirror, P6_Shrink_Mirror)
                list_regions.append(list_segments)
                list_segments = []

            for i in range(Q):
                draw_fraction(list_segments, iPark(P6_Shrink, i*alpha_slot_span), 
                                             iPark(P5_Shrink, i*alpha_slot_span), 
                                             iPark(P4_Shrink, i*alpha_slot_span), 
                                             iPark(POpen_Shrink, i*alpha_slot_span),
                                             iPark(P6_Shrink_Mirror, i*alpha_slot_span), 
                                             iPark(P5_Shrink_Mirror, i*alpha_slot_span), 
                                             iPark(P4_Shrink_Mirror, i*alpha_slot_span), 
                                             iPark(POpen_Shrink_Mirror, i*alpha_slot_span), )
                # raise
        else:
            # list_segments += drawer.drawCircle(PCoil, TheRadius)
            list_segments += drawer.drawArc([0,0], P6_Shrink, P5_Shrink)
            list_segments += drawer.drawLine(P5_Shrink, P4_Shrink)
            list_segments += drawer.drawLine(P4_Shrink, POpen_Shrink)
            list_segments += drawer.drawLine(POpen_Shrink, P6_Shrink)
            list_regions.append(list_segments)

            PCoil[1] *= -1
            P6_Shrink[1] *= -1
            P5_Shrink[1] *= -1
            P4_Shrink[1] *= -1
            POpen_Shrink[1] *= -1
            list_segments = []
            # list_segments += drawer.drawCircle(PCoil, TheRadius)
            list_segments += drawer.drawArc([0,0], P5_Shrink, P6_Shrink)
            list_segments += drawer.drawLine(P5_Shrink, P4_Shrink)
            list_segments += drawer.drawLine(P4_Shrink, POpen_Shrink)
            list_segments += drawer.drawLine(POpen_Shrink, P6_Shrink)
            list_regions.append(list_segments)
            list_segments = []

        # 我乱给的
        self.innerCoord = ( 0.5*(POpen[0]+P6[0]), 0.5*(POpen[1]+P6[1]))

        # return [list_segments]
        return {'innerCoord': self.innerCoord, 'list_regions':list_regions, 'mirrorAxis': None}



class CrossSectInnerRotorStator_PMAtYoke:
    # CrossSectInnerRotorStator_PMAtYoke Describes the inner-rotor motor stator 
    # that is dedicated for doubly salient PM motor.
    def __init__(self,
                    name  = 'StatorCore_PMAtYoke',
                    color = '#BAFD01',
                    deg_alpha_st = 40, # span angle of tooth: class type DimAngular
                    deg_alpha_sto = 20, # angle of tooth edge: class type DimAngular
                    mm_r_si      = 40, # inner radius of stator teeth: class type DimLinear
                    mm_d_sto      = 5,  # tooth edge length: class type DimLinear
                    mm_d_stt      = 10, # tooth tip length: class type DimLinear
                    mm_d_st      = 15, # tooth base length: class type DimLinear
                    mm_d_sy      = 15, # back iron thickness: class type DimLinear
                    mm_w_st      = 13, # tooth base width: class type DimLinear
                    mm_r_st      = 0,  # fillet on outter tooth: class type DimLinear
                    mm_r_sf      = 0,  # fillet between tooth tip and base: class type DimLinear
                    mm_r_sb      = 0,  # fillet at tooth base: class type DimLinear
                    Q            = 6,  # number of stator slots (integer)
                    mm_d_pm      = 5,
                    mm_difference_pm_yoke = 0.1, # mm
                    location = None
                ):


        self.name = name
        self.color = color
        # self.deg_alpha_st = deg_alpha_st
        self.Q = Q * 2 # 先搞起来

        self.alpha_slot_span = 2*np.pi / self.Q

        self.deg_alpha_st = (self.alpha_slot_span/2) /np.pi*180 *0.9 # 360 / self.Q / 2
        self.deg_alpha_sto = self.deg_alpha_st/2 #deg_alpha_sto
        self.mm_r_si      = mm_r_si     
        self.mm_d_sto     = mm_d_sto # * 0
        self.mm_d_stt     = mm_d_stt # * 0 
        self.mm_d_st      = mm_d_st     
        self.mm_d_sy      = mm_d_sy     
        # self.mm_w_st      = mm_w_st     
        # print(self.deg_alpha_st, self.alpha_slot_span/np.pi*180)
        self.mm_w_st = 2*(self.mm_r_si+self.mm_d_stt)*sin(0.5*(self.alpha_slot_span/2))
        self.mm_r_st      = mm_r_st     
        self.mm_r_sf      = mm_r_sf     
        self.mm_r_sb      = mm_r_sb  
        self.mm_d_pm =  mm_d_pm
        self.mm_difference_pm_yoke = mm_difference_pm_yoke
        self.location = location 

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        alpha_st = self.deg_alpha_st * np.pi/180
        alpha_so = -self.deg_alpha_sto * np.pi/180
        r_si = self.mm_r_si
        d_sto = self.mm_d_sto
        d_stt = self.mm_d_stt
        d_st = self.mm_d_st
        d_sy = self.mm_d_sy
        w_st = self.mm_w_st
        r_st = self.mm_r_st
        r_sf = self.mm_r_sf
        r_sb = self.mm_r_sb
        Q    = self.Q
        mm_d_pm = self.mm_d_pm

        alpha_slot_span = self.alpha_slot_span
        P1 = [r_si, 0]
        P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]
        P3_temp = [ d_sto*cos(alpha_st*0.5), 
                    d_sto*-sin(alpha_st*0.5)]
        P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
                             -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        P3 = [  P3_local_rotate[0] + P2[0],
                P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_stt
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度)]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        # Option 1
        # 6 [94.01649113009418, -25.191642873516543] 97.33303383272235
        # Radius_InnerStatorYoke = r_si+d_sp+d_st
        # print('Radius_InnerStatorYoke #1:', Radius_InnerStatorYoke)
        # Option 2
        r_sy = Radius_InnerStatorYoke = np.sqrt(P5[0]**2 + P5[1]**2)
        # print('Radius_InnerStatorYoke #2:', Radius_InnerStatorYoke)
        P6 = [ r_sy *  cos(alpha_slot_span*0.5),
               r_sy * -sin(alpha_slot_span*0.5) ]

        if False:
            P7 = [  (r_si+d_stt+d_st+d_sy)*cos(alpha_slot_span*0.5),
                    (r_si+d_stt+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
            P8 = [  r_si+d_stt+d_st+d_sy, 0]
        else:
            P7 = [  (r_sy+d_sy)*cos(alpha_slot_span*0.5),
                    (r_sy+d_sy)*-sin(alpha_slot_span*0.5) ]
            P8 = [   r_sy+d_sy, 0]

        # Draw slot for inserting PM
        P0 = [r_sy, 0]
        # print('||DEBUG', P0, np.sqrt(P5[0]**2 + P5[1]**2))
        print(mm_d_pm, r_sy, P5)
        alpha_pm_depth = 2*np.arcsin(0.5*mm_d_pm / r_sy)
        # print('magnet angle', alpha_pm_depth/np.pi*180)
        mm_w_pm = d_sy - self.mm_difference_pm_yoke

        list_segments = []
        if bool_draw_whole_model:
            P2_Mirror = [P2[0], -P2[1]] # = iPark(P2, alpha_st)
            P3_Mirror = [P3[0], -P3[1]]
            P4_Mirror = [P4[0], -P4[1]]
            P5_Mirror = [P5[0], -P5[1]]
            def iPark(P, theta):
                return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
            def draw_fraction(list_segments, i, P0, P2, P3, P4, P5,
                                            P2_Mirror, P3_Mirror, P4_Mirror, P5_Mirror):
                P5_Rotate = iPark(P5, alpha_slot_span)
                list_segments += drawer.drawArc([0,0], P2, P2_Mirror)
                list_segments += drawer.drawLine(P2, P3)
                list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
                list_segments += drawer.drawLine(P3, P4)
                list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
                list_segments += drawer.drawLine(P4, P5)
                list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)
                if i % 2 == 1:
                    list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate)
                else:
                    # print(alpha_slot_span/2/np.pi*180,  alpha_pm_depth/2/np.pi*180)
                    P5_Rotate2PM = iPark(P0, alpha_slot_span/2 - alpha_pm_depth/2)
                    list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate2PM)
                    P5_RotatePassPM = iPark(P5_Rotate2PM, alpha_pm_depth)
                    P51 = [ P5_Rotate2PM[0]    + mm_w_pm*cos(alpha_slot_span/2+alpha_slot_span*i), P5_Rotate2PM[1]    + mm_w_pm*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P52 = [ P5_RotatePassPM[0] + mm_w_pm*cos(alpha_slot_span/2+alpha_slot_span*i), P5_RotatePassPM[1] + mm_w_pm*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    list_segments += drawer.drawLine(P5_Rotate2PM, P51)
                    list_segments += drawer.drawLine(P51, P52)
                    list_segments += drawer.drawLine(P52, P5_RotatePassPM)
                    list_segments += drawer.drawArc([0,0], P5_RotatePassPM, P5_Rotate)
            for i in range(Q):
                draw_fraction(list_segments, i,
                                             iPark(P0, i*alpha_slot_span), 
                                             iPark(P2, i*alpha_slot_span), 
                                             iPark(P3, i*alpha_slot_span), 
                                             iPark(P4, i*alpha_slot_span), 
                                             iPark(P5, i*alpha_slot_span),
                                             iPark(P2_Mirror, i*alpha_slot_span), 
                                             iPark(P3_Mirror, i*alpha_slot_span), 
                                             iPark(P4_Mirror, i*alpha_slot_span), 
                                             iPark(P5_Mirror, i*alpha_slot_span), )
                # if i >=1:
                #     break
                # raise
            # draw a circle (this is officially suggested)
            list_segments += drawer.drawArc([0,0], P8, [-P8[0], P8[1]])
            list_segments += drawer.drawArc([0,0],     [-P8[0], P8[1]], P8)


            # DEBUG
            # for ind, point in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
            #     print(ind+1, point, np.sqrt(point[0]**2+point[1]**2))

            self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

            # return [list_segments]
            return {'innerCoord': self.innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])],
                    'inner_or_outer_region_to_remove': [False, True]
                    }
        else:
            raise Exception('Not supported. Please set bool_draw_whole_model to True.')

class CrossSectStatorMagnetAtYoke:
    def __init__(self,
                    name  = 'StatorMagnetAtYoke',
                    color = '#0BA0E2',
                    stator_core = None,
                ):


        self.name = name
        self.color = color
        self.stator_core = stator_core
    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        alpha_st = self.stator_core.deg_alpha_st * np.pi/180
        alpha_so = -self.stator_core.deg_alpha_sto * np.pi/180
        r_si = self.stator_core.mm_r_si
        d_so = self.stator_core.mm_d_sto
        d_sp = self.stator_core.mm_d_stt
        d_st = self.stator_core.mm_d_st
        d_sy = self.stator_core.mm_d_sy
        w_st = self.stator_core.mm_w_st
        r_st = self.stator_core.mm_r_st
        r_sf = self.stator_core.mm_r_sf
        r_sb = self.stator_core.mm_r_sb
        Q    = self.stator_core.Q
        mm_d_pm = self.stator_core.mm_d_pm

        alpha_slot_span = self.stator_core.alpha_slot_span
        P1 = [r_si, 0]
        P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]
        P3_temp = [ d_so*cos(alpha_st*0.5), 
                    d_so*-sin(alpha_st*0.5)]
        P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
                             -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        P3 = [  P3_local_rotate[0] + P2[0],
                P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_sp
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度)]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        # Option 1
        # 6 [94.01649113009418, -25.191642873516543] 97.33303383272235
        # Radius_InnerStatorYoke = r_si+d_sp+d_st
        # print('Radius_InnerStatorYoke #1:', Radius_InnerStatorYoke)
        # Option 2
        Radius_InnerStatorYoke = np.sqrt(P5[0]**2 + P5[1]**2)
        # print('Radius_InnerStatorYoke #2:', Radius_InnerStatorYoke)
        P6 = [ Radius_InnerStatorYoke *  cos(alpha_slot_span*0.5),
               Radius_InnerStatorYoke * -sin(alpha_slot_span*0.5) ]

        P7 = [ (r_si+d_sp+d_st+d_sy)*cos(alpha_slot_span*0.5),
               (r_si+d_sp+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
        P8 = [  r_si+d_sp+d_st+d_sy, 0]

        # Draw slot for inserting PM
        r_sy = Radius_InnerStatorYoke
        P0 = [r_sy, 0]
        # print('||DEBUG', P0, np.sqrt(P5[0]**2 + P5[1]**2))
        alpha_pm_depth = 2*np.arcsin(0.5*mm_d_pm / r_sy)
        # print('magnet angle', alpha_pm_depth/np.pi*180)
        mm_w_pm = d_sy - self.stator_core.mm_difference_pm_yoke

        self.maly=[]
        list_segments = []
        if bool_draw_whole_model:
            P2_Mirror = [P2[0], -P2[1]] # = iPark(P2, alpha_st)
            P3_Mirror = [P3[0], -P3[1]]
            P4_Mirror = [P4[0], -P4[1]]
            P5_Mirror = [P5[0], -P5[1]]
            def iPark(P, theta):
                return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
            def draw_fraction(list_segments, i, P0, P2, P3, P4, P5,
                                            P2_Mirror, P3_Mirror, P4_Mirror, P5_Mirror):
                P5_Rotate = iPark(P5, alpha_slot_span)
                # list_segments += drawer.drawArc([0,0], P2, P2_Mirror)
                # list_segments += drawer.drawLine(P2, P3)
                # list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
                # list_segments += drawer.drawLine(P3, P4)
                # list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
                # list_segments += drawer.drawLine(P4, P5)
                # list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)
                if i % 2 == 1:
                    # list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate)
                    pass
                else:
                    magnet = dict()
                    # print(alpha_slot_span/2/np.pi*180,  alpha_pm_depth/2/np.pi*180)
                    P5_Rotate2PM = iPark(P0, alpha_slot_span/2 - alpha_pm_depth/2)
                    # list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate2PM)
                    P5_RotatePassPM = iPark(P5_Rotate2PM, alpha_pm_depth)
                    P51 = [ P5_Rotate2PM[0]    + mm_w_pm*cos(alpha_slot_span/2+alpha_slot_span*i), P5_Rotate2PM[1]    + mm_w_pm*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P52 = [ P5_RotatePassPM[0] + mm_w_pm*cos(alpha_slot_span/2+alpha_slot_span*i), P5_RotatePassPM[1] + mm_w_pm*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    list_segments += drawer.drawLine(P5_Rotate2PM, P51)
                    list_segments += drawer.drawLine(P51, P52)
                    list_segments += drawer.drawLine(P52, P5_RotatePassPM)
                    list_segments += drawer.drawLine(P5_RotatePassPM, P5_Rotate2PM)
                    # list_segments += drawer.drawArc([0,0], P5_RotatePassPM, P5_Rotate)
                    magnet['X'] = 0.5*(P5_Rotate2PM[0]+P52[0])
                    magnet['Y'] = 0.5*(P5_Rotate2PM[1]+P52[1])
                    magnet['Direction-Theta'] = np.arctan2(P52[1]-P51[1], P52[0]-P51[0])
                    self.maly.append(magnet)
            for i in range(Q):
                draw_fraction(list_segments, i,
                                             iPark(P0, i*alpha_slot_span), 
                                             iPark(P2, i*alpha_slot_span), 
                                             iPark(P3, i*alpha_slot_span), 
                                             iPark(P4, i*alpha_slot_span), 
                                             iPark(P5, i*alpha_slot_span),
                                             iPark(P2_Mirror, i*alpha_slot_span), 
                                             iPark(P3_Mirror, i*alpha_slot_span), 
                                             iPark(P4_Mirror, i*alpha_slot_span), 
                                             iPark(P5_Mirror, i*alpha_slot_span), )
                # if i >=1:
                #     break
                # raise
            # draw a circle (this is officially suggested)
            # list_segments += drawer.drawArc([0,0], P8, [-P8[0], P8[1]])
            # list_segments += drawer.drawArc([0,0],     [-P8[0], P8[1]], P8)


            # DEBUG
            # for ind, point in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
            #     print(ind+1, point, np.sqrt(point[0]**2+point[1]**2))

            self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

            # return [list_segments]
            return {'innerCoord': self.innerCoord, 
                    'list_regions':[list_segments], 
                    # 'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])],
                    # 'inner_or_outer_region_to_remove': [False, True]
                    }
        else:
            raise Exception('Not supported. Please set bool_draw_whole_model to True.')

class CrossSectToroidalWiniding:
    def __init__(self, 
                    name = 'CoilToroidal',
                    color = '#3D9970',
                    stator_core = None
                    ):
        self.type = 'Toroidal'
        self.name = name
        self.color = color
        self.stator_core = stator_core

    def draw(self, drawer, bool_draw_whole_model=True):

        drawer.getSketch(self.name, self.color)

        alpha_st = self.stator_core.deg_alpha_st * np.pi/180
        alpha_so = -self.stator_core.deg_alpha_sto * np.pi/180
        r_si = self.stator_core.mm_r_si
        d_so = self.stator_core.mm_d_sto
        d_sp = self.stator_core.mm_d_stt
        d_st = self.stator_core.mm_d_st
        d_sy = self.stator_core.mm_d_sy
        w_st = self.stator_core.mm_w_st
        r_st = self.stator_core.mm_r_st
        r_sf = self.stator_core.mm_r_sf
        r_sb = self.stator_core.mm_r_sb
        Q    = self.stator_core.Q
        mm_d_pm = self.stator_core.mm_d_pm

        alpha_slot_span = self.stator_core.alpha_slot_span
        P1 = [r_si, 0]
        P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]
        P3_temp = [ d_so*cos(alpha_st*0.5), 
                    d_so*-sin(alpha_st*0.5)]
        P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
                             -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        P3 = [  P3_local_rotate[0] + P2[0],
                P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_sp
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度)]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        # Option 1
        # 6 [94.01649113009418, -25.191642873516543] 97.33303383272235
        # Radius_InnerStatorYoke = r_si+d_sp+d_st
        # print('Radius_InnerStatorYoke #1:', Radius_InnerStatorYoke)
        # Option 2
        Radius_InnerStatorYoke = np.sqrt(P5[0]**2 + P5[1]**2)
        # print('Radius_InnerStatorYoke #2:', Radius_InnerStatorYoke)
        P6 = [ Radius_InnerStatorYoke *  cos(alpha_slot_span*0.5),
               Radius_InnerStatorYoke * -sin(alpha_slot_span*0.5) ]

        P7 = [ (r_si+d_sp+d_st+d_sy)*cos(alpha_slot_span*0.5),
               (r_si+d_sp+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
        P8 = [  r_si+d_sp+d_st+d_sy, 0]

        # Draw slot for inserting PM
        r_sy = Radius_InnerStatorYoke
        P0 = [r_sy, 0]
        # print('||DEBUG', P0, np.sqrt(P5[0]**2 + P5[1]**2))
        alpha_pm_depth = 2*np.arcsin(0.5*mm_d_pm / r_sy)
        # print('magnet angle', alpha_pm_depth/np.pi*180)
        mm_w_pm = d_sy - self.stator_core.mm_difference_pm_yoke

        self.wily = list()
        list_segments = []
        if bool_draw_whole_model:
            P2_Mirror = [P2[0], -P2[1]] # = iPark(P2, alpha_st)
            P3_Mirror = [P3[0], -P3[1]]
            P4_Mirror = [P4[0], -P4[1]]
            P5_Mirror = [P5[0], -P5[1]]
            def iPark(P, theta):
                return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]
            def draw_fraction(list_segments, i, P0, P2, P3, P4, P5,
                                            P2_Mirror, P3_Mirror, P4_Mirror, P5_Mirror):
                P5_Rotate = iPark(P5, alpha_slot_span)
                # list_segments += drawer.drawArc([0,0], P2, P2_Mirror)
                # list_segments += drawer.drawLine(P2, P3)
                # list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
                # list_segments += drawer.drawLine(P3, P4)
                # list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
                # list_segments += drawer.drawLine(P4, P5)
                # list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)
                # if i % 2 == 1:
                #     # list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate)
                #     pass
                # else:
                if True:
                    # print(alpha_slot_span/2/np.pi*180,  alpha_pm_depth/2/np.pi*180)
                    P5_Rotate2PM = iPark(P0, alpha_slot_span/2 - alpha_pm_depth/2)
                    # list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate2PM)
                    P5_RotatePassPM = iPark(P5_Rotate2PM, alpha_pm_depth)
                    P53  = [ P5_Rotate2PM[0]    + d_sy* 1.2*cos(alpha_slot_span/2+alpha_slot_span*i), P5_Rotate2PM[1]    + d_sy* 1.2*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P54  = [ P5_RotatePassPM[0] + d_sy* 1.2*cos(alpha_slot_span/2+alpha_slot_span*i), P5_RotatePassPM[1] + d_sy* 1.2*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P55  = [ P5_Rotate2PM[0]    + d_sy* 1.4*cos(alpha_slot_span/2+alpha_slot_span*i), P5_Rotate2PM[1]    + d_sy* 1.4*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P56  = [ P5_RotatePassPM[0] + d_sy* 1.4*cos(alpha_slot_span/2+alpha_slot_span*i), P5_RotatePassPM[1] + d_sy* 1.4*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P57  = [ P5_Rotate2PM[0]    + d_sy*-0.2*cos(alpha_slot_span/2+alpha_slot_span*i), P5_Rotate2PM[1]    + d_sy*-0.2*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P58  = [ P5_RotatePassPM[0] + d_sy*-0.2*cos(alpha_slot_span/2+alpha_slot_span*i), P5_RotatePassPM[1] + d_sy*-0.2*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P59  = [ P5_Rotate2PM[0]    + d_sy*-0.4*cos(alpha_slot_span/2+alpha_slot_span*i), P5_Rotate2PM[1]    + d_sy*-0.4*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    P510 = [ P5_RotatePassPM[0] + d_sy*-0.4*cos(alpha_slot_span/2+alpha_slot_span*i), P5_RotatePassPM[1] + d_sy*-0.4*sin(alpha_slot_span/2+alpha_slot_span*i) ]
                    list_segments += drawer.drawLine(P53, P54)
                    list_segments += drawer.drawLine(P54, P56)
                    list_segments += drawer.drawLine(P56, P55)
                    list_segments += drawer.drawLine(P55, P53)
                    coil=dict()
                    coil['LayerX-Phase'] = 'B' if i%2==0 else 'M'
                    coil['LayerX-X'] = 0.5*(P53[0] + P56[0])
                    coil['LayerX-Y'] = 0.5*(P53[1] + P56[1])
                    coil['LayerX-R'] = np.sqrt(coil['LayerX-X']**2 + coil['LayerX-Y']**2)
                    coil['LayerX-Theta'] = np.arctan2(coil['LayerX-Y'], coil['LayerX-X'])
                    coil['LayerX-Direction'] = '+'
                    coil['LayerX-Purpose'] = 'SuspensionOnly' if i%2==0 else 'DualPurpose'
                    list_segments += drawer.drawLine(P57, P58)
                    list_segments += drawer.drawLine(P58, P510)
                    list_segments += drawer.drawLine(P510, P59)
                    list_segments += drawer.drawLine(P59, P57)
                    coil['LayerY-Phase'] = 'B' if i%2==0 else 'M'
                    coil['LayerY-X'] = 0.5*(P57[0] + P510[0])
                    coil['LayerY-Y'] = 0.5*(P57[1] + P510[1])
                    coil['LayerY-R'] = np.sqrt(coil['LayerY-X']**2 + coil['LayerY-Y']**2)
                    coil['LayerY-Theta'] = np.arctan2(coil['LayerY-Y'], coil['LayerY-X'])
                    coil['LayerY-Direction'] = '-'
                    coil['LayerY-Purpose'] = 'SuspensionOnly' if i%2==0 else 'DualPurpose'
                    self.wily.append(coil)
            for i in range(Q):
                draw_fraction(list_segments, i,
                                             iPark(P0, i*alpha_slot_span), 
                                             iPark(P2, i*alpha_slot_span), 
                                             iPark(P3, i*alpha_slot_span), 
                                             iPark(P4, i*alpha_slot_span), 
                                             iPark(P5, i*alpha_slot_span),
                                             iPark(P2_Mirror, i*alpha_slot_span), 
                                             iPark(P3_Mirror, i*alpha_slot_span), 
                                             iPark(P4_Mirror, i*alpha_slot_span), 
                                             iPark(P5_Mirror, i*alpha_slot_span), )
                # if i >=1:
                #     break
                # raise
            # draw a circle (this is officially suggested)
            # list_segments += drawer.drawArc([0,0], P8, [-P8[0], P8[1]])
            # list_segments += drawer.drawArc([0,0],     [-P8[0], P8[1]], P8)

            # DEBUG
            # for ind, point in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
            #     print(ind+1, point, np.sqrt(point[0]**2+point[1]**2))

            # self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

            # return [list_segments]
            return {
                    # 'innerCoord': self.innerCoord, 
                    'list_regions':[list_segments], 
                    # 'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])],
                    # 'inner_or_outer_region_to_remove': [False, True]
                    }
        else:
            raise Exception('Not supported. Please set bool_draw_whole_model to True.')




class CrossSectInnerRotorStator_PMAtToothBody:
    def __init__(self,
                    name  = 'StatorCore_PMAtToothBody',
                    color = '#BAFD01',
                    deg_alpha_st = 40, # span angle of tooth: class type DimAngular
                    deg_alpha_sto = 20, # angle of tooth edge: class type DimAngular
                    mm_r_si      = 40, # inner radius of stator teeth: class type DimLinear
                    mm_d_sto      = 5,  # tooth edge length: class type DimLinear
                    mm_d_stt      = 10, # tooth tip length: class type DimLinear
                    mm_d_st      = 15, # tooth base length: class type DimLinear
                    mm_d_sy      = 15, # back iron thickness: class type DimLinear
                    mm_w_st      = 13, # tooth base width: class type DimLinear
                    mm_r_st      = 0,  # fillet on outter tooth: class type DimLinear
                    mm_r_sf      = 0,  # fillet between tooth tip and base: class type DimLinear
                    mm_r_sb      = 0,  # fillet at tooth base: class type DimLinear
                    Q            = 6,  # number of stator slots (integer)
                    mm_d_pm      = 5,
                    deg_alpha_pm_at_airgap = None,
                    mm_difference_pm_yoke = 1, # mm
                    location = None
                ):


        self.name = name
        self.color = color
        # self.deg_alpha_st = deg_alpha_st
        self.Q = Q

        self.alpha_slot_span = 2*np.pi / self.Q

        self.deg_alpha_st = deg_alpha_st # (self.alpha_slot_span/2) /np.pi*180 *0.9 # 360 / self.Q / 2
        self.deg_alpha_sto = self.deg_alpha_st/2 #deg_alpha_sto
        self.mm_r_si      = mm_r_si     
        self.mm_d_sto     = mm_d_sto # * 0
        self.mm_d_stt     = mm_d_stt # * 0 
        self.mm_d_st      = mm_d_st     
        self.mm_d_sy      = mm_d_sy     
        # self.mm_w_st      = mm_w_st     
        # print(self.deg_alpha_st, self.alpha_slot_span/np.pi*180)
        self.mm_w_st = 2*(self.mm_r_si+self.mm_d_stt)*sin(0.5*(self.alpha_slot_span/2))
        self.mm_r_st      = mm_r_st     
        self.mm_r_sf      = mm_r_sf     
        self.mm_r_sb      = mm_r_sb  
        self.mm_d_pm =  mm_d_pm
        self.deg_alpha_pm_at_airgap = deg_alpha_pm_at_airgap
        self.mm_difference_pm_yoke = mm_difference_pm_yoke
        self.location = location 

    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        alpha_st = self.deg_alpha_st   * np.pi/180
        # print('DEBUBG', self.deg_alpha_st)

        alpha_so = -self.deg_alpha_sto * np.pi/180
        r_si = self.mm_r_si
        d_sto = self.mm_d_sto
        d_stt = self.mm_d_stt
        d_st = self.mm_d_st
        d_sy = self.mm_d_sy
        w_st = self.mm_w_st
        r_st = self.mm_r_st
        r_sf = self.mm_r_sf
        r_sb = self.mm_r_sb
        Q    = self.Q
        mm_d_pm = self.mm_d_pm
        alpha_pm_at_airgap = self.deg_alpha_pm_at_airgap/180*np.pi

        alpha_slot_span = self.alpha_slot_span
        print('DEUBG', alpha_slot_span/np.pi*180)

        P0 = [r_si, 0]
        P1 = [r_si* cos(alpha_pm_at_airgap/2), 
              r_si*-sin(alpha_pm_at_airgap/2) ]
        P2 = [r_si* cos(alpha_st/2), 
              r_si*-sin(alpha_st/2) ]

        # P1 = [r_si, 0]
        # P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]
        P3_temp = [ d_sto*cos(alpha_st*0.5), 
                    d_sto*-sin(alpha_st*0.5)]
        P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
                             -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        P3 = [  P3_local_rotate[0] + P2[0],
                P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_stt
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度)]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        # Option 1
        # 6 [94.01649113009418, -25.191642873516543] 97.33303383272235
        # Radius_InnerStatorYoke = r_si+d_sp+d_st
        # print('Radius_InnerStatorYoke #1:', Radius_InnerStatorYoke)
        # Option 2
        r_sy = Radius_InnerStatorYoke = np.sqrt(P5[0]**2 + P5[1]**2)
        # print('Radius_InnerStatorYoke #2:', Radius_InnerStatorYoke)
        P6 = [ r_sy *  cos(alpha_slot_span*0.5),
               r_sy * -sin(alpha_slot_span*0.5) ]

        if False:
            P7 = [  (r_si+d_stt+d_st+d_sy)*cos(alpha_slot_span*0.5),
                    (r_si+d_stt+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
            P8 = [  r_si+d_stt+d_st+d_sy, 0]
        else:
            P7 = [  (r_sy+d_sy)*cos(alpha_slot_span*0.5),
                    (r_sy+d_sy)*-sin(alpha_slot_span*0.5) ]
            P8 = [   r_sy+d_sy, 0]

        r_so = r_sy + d_sy
        P_PM = [np.sqrt((r_so)**2 - P1[1]**2) - self.mm_difference_pm_yoke, P1[1]]

        斜边 = r_so
        高 = np.abs(P5[1])
        底边 = np.sqrt(斜边**2 - 高**2)
        P5_OuterEdge = [底边, -高]

        P1_Mirror = [P1[0], -P1[1]]
        P2_Mirror = [P2[0], -P2[1]] # = iPark(P2, alpha_st)
        P3_Mirror = [P3[0], -P3[1]]
        P4_Mirror = [P4[0], -P4[1]]
        P5_Mirror = [P5[0], -P5[1]]
        P_PM_Mirror = [P_PM[0], -P_PM[1]]
        def iPark(P, theta):
            return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]

        list_segments = []
        if bool_draw_whole_model:
            def draw_fraction(list_segments, i, P1, P2, P3, P4, P5,
                                            P1_Mirror,
                                            P2_Mirror, P3_Mirror, P4_Mirror, P5_Mirror,
                                            P_PM, P_PM_Mirror):
                P5_Rotate = iPark(P5, alpha_slot_span)
                list_segments += drawer.drawArc([0,0], P2, P1)
                list_segments += drawer.drawLine(P2, P3)
                list_segments += drawer.drawLine(P3, P4)
                list_segments += drawer.drawLine(P4, P5)

                list_segments += drawer.drawArc([0,0], P1_Mirror, P2_Mirror)
                list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
                list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
                list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)

                list_segments += drawer.drawLine(P1, P_PM)
                list_segments += drawer.drawLine(P1_Mirror, P_PM_Mirror)
                list_segments += drawer.drawLine(P_PM, P_PM_Mirror)

                list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate)

            for i in range(Q):
                draw_fraction(list_segments, i,
                                             iPark(P1, i*alpha_slot_span), 
                                             iPark(P2, i*alpha_slot_span), 
                                             iPark(P3, i*alpha_slot_span), 
                                             iPark(P4, i*alpha_slot_span), 
                                             iPark(P5, i*alpha_slot_span),
                                             iPark(P1_Mirror, i*alpha_slot_span), 
                                             iPark(P2_Mirror, i*alpha_slot_span), 
                                             iPark(P3_Mirror, i*alpha_slot_span), 
                                             iPark(P4_Mirror, i*alpha_slot_span), 
                                             iPark(P5_Mirror, i*alpha_slot_span),
                                             iPark(P_PM, i*alpha_slot_span),
                                             iPark(P_PM_Mirror, i*alpha_slot_span), )
                # if i >=1:
                #     break
                # raise
                # break
            # draw a circle (this is officially suggested)
            list_segments += drawer.drawArc([0,0], P8, [-P8[0], P8[1]])
            list_segments += drawer.drawArc([0,0],     [-P8[0], P8[1]], P8)


            # DEBUG
            # for ind, point in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
            #     print(ind+1, point, np.sqrt(point[0]**2+point[1]**2))

            self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

            # return [list_segments]
            return {'innerCoord': self.innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])],
                    'inner_or_outer_region_to_remove': [False, True]
                    }
        else:
            P5_Rotate = iPark(P5, alpha_slot_span)
            P5_OuterEdge_Rotate = iPark(P5_OuterEdge, alpha_slot_span)
            list_segments += drawer.drawArc([0,0], P2, P1)
            list_segments += drawer.drawLine(P2, P3)
            list_segments += drawer.drawLine(P3, P4)
            list_segments += drawer.drawLine(P4, P5_OuterEdge)
            list_segments += drawer.drawLine(P5_Rotate, P5_OuterEdge_Rotate)
            list_segments += drawer.drawArc([0,0], P5_OuterEdge, P5_OuterEdge_Rotate)

            list_segments += drawer.drawArc([0,0], P1_Mirror, P2_Mirror)
            list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
            list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
            list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)

            list_segments += drawer.drawLine(P1, P_PM)
            list_segments += drawer.drawLine(P1_Mirror, P_PM_Mirror)
            list_segments += drawer.drawLine(P_PM, P_PM_Mirror)

            list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate)

            self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

            # return [list_segments]
            return {'innerCoord': self.innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])],
                    # 'inner_or_outer_region_to_remove': [False, True]
                    }
            # raise Exception('Not supported. Please set bool_draw_whole_model to True.')

class CrossSectStatorMagnetAtToothBody:
    def __init__(self,
                    name  = 'StatorMagnetAtToothBody',
                    color = '#0BA0E2',
                    stator_core = None,
                    mm_d_air_pm = 5,
                ):
        self.name = name
        self.color = color
        self.stator_core = stator_core
        self.mm_d_air_pm = mm_d_air_pm
    def draw(self, drawer, bool_draw_whole_model=False):

        drawer.getSketch(self.name, self.color)

        alpha_st = self.stator_core.deg_alpha_st * np.pi/180
        alpha_so = -self.stator_core.deg_alpha_sto * np.pi/180
        r_si = self.stator_core.mm_r_si
        d_sto = self.stator_core.mm_d_sto
        d_stt = self.stator_core.mm_d_stt
        d_st = self.stator_core.mm_d_st
        d_sy = self.stator_core.mm_d_sy
        w_st = self.stator_core.mm_w_st
        r_st = self.stator_core.mm_r_st
        r_sf = self.stator_core.mm_r_sf
        r_sb = self.stator_core.mm_r_sb
        Q    = self.stator_core.Q
        mm_d_pm = self.stator_core.mm_d_pm

        mm_d_air_pm = self.mm_d_air_pm

        alpha_slot_span = self.stator_core.alpha_slot_span

        alpha_pm_at_airgap = self.stator_core.deg_alpha_pm_at_airgap/180*np.pi

        P0 = [r_si + mm_d_air_pm, 0]
        P1 = [r_si* cos(alpha_pm_at_airgap/2) + mm_d_air_pm, 
              r_si*-sin(alpha_pm_at_airgap/2) ]
        P2 = [r_si* cos(alpha_st/2), 
              r_si*-sin(alpha_st/2) ]

        # P1 = [r_si, 0]
        # P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]
        P3_temp = [ d_sto*cos(alpha_st*0.5), 
                    d_sto*-sin(alpha_st*0.5)]
        P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
                             -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        P3 = [  P3_local_rotate[0] + P2[0],
                P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_stt
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度)]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        # Option 1
        # 6 [94.01649113009418, -25.191642873516543] 97.33303383272235
        # Radius_InnerStatorYoke = r_si+d_sp+d_st
        # print('Radius_InnerStatorYoke #1:', Radius_InnerStatorYoke)
        # Option 2
        r_sy = Radius_InnerStatorYoke = np.sqrt(P5[0]**2 + P5[1]**2)
        # print('Radius_InnerStatorYoke #2:', Radius_InnerStatorYoke)
        P6 = [ r_sy *  cos(alpha_slot_span*0.5),
               r_sy * -sin(alpha_slot_span*0.5) ]

        if False:
            P7 = [  (r_si+d_stt+d_st+d_sy)*cos(alpha_slot_span*0.5),
                    (r_si+d_stt+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
            P8 = [  r_si+d_stt+d_st+d_sy, 0]
        else:
            P7 = [  (r_sy+d_sy)*cos(alpha_slot_span*0.5),
                    (r_sy+d_sy)*-sin(alpha_slot_span*0.5) ]
            P8 = [   r_sy+d_sy, 0]

        r_so = r_sy + d_sy
        P_PM = [np.sqrt((r_so)**2 - P1[1]**2) - self.stator_core.mm_difference_pm_yoke, P1[1]]

        P1_Mirror = [P1[0], -P1[1]]
        P2_Mirror = [P2[0], -P2[1]] # = iPark(P2, alpha_st)
        P3_Mirror = [P3[0], -P3[1]]
        P4_Mirror = [P4[0], -P4[1]]
        P5_Mirror = [P5[0], -P5[1]]
        P_PM_Mirror = [P_PM[0], -P_PM[1]]
        def iPark(P, theta):
            return [P[0]*np.cos(theta)+P[1]*-np.sin(theta), P[0]*np.sin(theta)+P[1]*np.cos(theta)]

        self.mm2_magnet_area = np.sqrt( (         P1[0]-P_PM[0])**2 + (         P1[1]-P_PM[1])**2 ) \
                             * np.sqrt( (P_PM_Mirror[0]-P_PM[0])**2 + (P_PM_Mirror[1]-P_PM[1])**2 ) * Q # Note Q = 2*pe

        list_segments = []
        if bool_draw_whole_model:
            def draw_fraction(list_segments, i, P1, P2, P3, P4, P5,
                                            P1_Mirror,
                                            P2_Mirror, P3_Mirror, P4_Mirror, P5_Mirror,
                                            P_PM, P_PM_Mirror):
                P5_Rotate = iPark(P5, alpha_slot_span)
                # list_segments += drawer.drawArc([0,0], P2, P1)
                # list_segments += drawer.drawLine(P2, P3)
                # list_segments += drawer.drawLine(P3, P4)
                # list_segments += drawer.drawLine(P4, P5)

                # list_segments += drawer.drawArc([0,0], P1_Mirror, P2_Mirror)
                # list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
                # list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
                # list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)

                list_segments += drawer.drawLine(P1, P1_Mirror)
                list_segments += drawer.drawLine(P1, P_PM)
                list_segments += drawer.drawLine(P_PM, P_PM_Mirror)
                list_segments += drawer.drawLine(P1_Mirror, P_PM_Mirror)

                # list_segments += drawer.drawArc([0,0], P5_Mirror, P5_Rotate)

            for i in range(Q):
                draw_fraction(list_segments, i,
                                             iPark(P1, i*alpha_slot_span), 
                                             iPark(P2, i*alpha_slot_span), 
                                             iPark(P3, i*alpha_slot_span), 
                                             iPark(P4, i*alpha_slot_span), 
                                             iPark(P5, i*alpha_slot_span),
                                             iPark(P1_Mirror, i*alpha_slot_span), 
                                             iPark(P2_Mirror, i*alpha_slot_span), 
                                             iPark(P3_Mirror, i*alpha_slot_span), 
                                             iPark(P4_Mirror, i*alpha_slot_span), 
                                             iPark(P5_Mirror, i*alpha_slot_span),
                                             iPark(P_PM, i*alpha_slot_span),
                                             iPark(P_PM_Mirror, i*alpha_slot_span), )
                # if i >=1:
                #     break
                # raise
                # break
            # draw a circle (this is officially suggested)
            # list_segments += drawer.drawArc([0,0], P8, [-P8[0], P8[1]])
            # list_segments += drawer.drawArc([0,0],     [-P8[0], P8[1]], P8)


            # DEBUG
            # for ind, point in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
            #     print(ind+1, point, np.sqrt(point[0]**2+point[1]**2))

            self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

            # return [list_segments]
            return {'innerCoord': self.innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])],
                    'inner_or_outer_region_to_remove': [False, True]
                    }
        else:
            P5_Rotate = iPark(P5, alpha_slot_span)
            # list_segments += drawer.drawArc([0,0], P2, P1)
            # list_segments += drawer.drawLine(P2, P3)
            # list_segments += drawer.drawLine(P3, P4)
            # list_segments += drawer.drawLine(P4, P5)

            # list_segments += drawer.drawArc([0,0], P1_Mirror, P2_Mirror)
            # list_segments += drawer.drawLine(P2_Mirror, P3_Mirror)
            # list_segments += drawer.drawLine(P3_Mirror, P4_Mirror)
            # list_segments += drawer.drawLine(P4_Mirror, P5_Mirror)

            list_segments += drawer.drawLine(P1, P1_Mirror)
            list_segments += drawer.drawLine(P1, P_PM)
            list_segments += drawer.drawLine(P_PM, P_PM_Mirror)
            list_segments += drawer.drawLine(P1_Mirror, P_PM_Mirror)
            self.innerCoord = ( 0.5*(P1[0]+P5[0]), 0.5*(P1[1]+P5[1]))

            # return [list_segments]
            return {'innerCoord': self.innerCoord, 
                    'list_regions':[list_segments], 
                    'mirrorAxis': [(P8[0]+5, P8[1]), (P8[0]+15, P8[1])],
                    # 'inner_or_outer_region_to_remove': [False, True]
                    }
            raise Exception('Not supported. Please set bool_draw_whole_model to True.')

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

    stator_core = CrossSectInnerRotorStator( name = 'StatorCore',
                                        deg_alpha_st = 40,
                                        deg_alpha_sto = 20,
                                        mm_r_si = 40,
                                        mm_d_sto = 5,
                                        mm_d_stt = 10,
                                        mm_d_st = 15,
                                        mm_d_sy = 15,
                                        mm_w_st = 13,
                                        mm_r_st = 0,
                                        mm_r_sf = 0,
                                        mm_r_sb = 0,
                                        Q = 6,
                                        location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                        )

    list_regions = stator_core.draw(toolJd)
    toolJd.bMirror = True
    toolJd.iRotateCopy = stator_core.Q
    region1 = toolJd.prepareSection(list_regions)

    coils = CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                stator_core = stator_core)

    list_regions = coils.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = coils.stator_core.Q
    region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')

