from math import cos, sin
import math
import logging
EPS = 1e-3 # [mm]
import numpy as np

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
        """Calculate all point coordinates based on the geometric parameters."""
        mm_d_pm  = self.mm_d_pm
        alpha_rm = self.deg_alpha_rm * math.pi/180
        alpha_rs = self.deg_alpha_rs * math.pi/180
        r_ri     = self.mm_r_ri
        d_ri     = self.mm_d_ri
        d_rp     = self.mm_d_rp
        d_rs     = self.mm_d_rs
        p        = self.p
        s        = self.s
        alpha_rp = 2*math.pi/(2*p) # pole span

        # Adjust alpha_rm if it's too close to alpha_rp
        if abs(alpha_rp - alpha_rm) <= 2 * math.pi/180:
            alpha_rm = alpha_rp
            if s == 1:
                alpha_rs = alpha_rm

        # Basic points
        P1 = [r_ri, 0]

        r_P2 = r_ri + d_ri + d_rp
        P2 = [r_P2, 0]

        alpha_P3 = alpha_rp - alpha_rm
        P3 = [r_P2*cos(alpha_P3), r_P2*-sin(alpha_P3)]

        r_P4 = r_ri + d_ri
        P4 = [r_P4*cos(alpha_P3), r_P4*-sin(alpha_P3)]

        alpha_P5 = alpha_P3 + alpha_rs
        if abs(alpha_rs*s - alpha_rm) < EPS:
            alpha_P5 = alpha_P3 + alpha_rs*s
        P5 = [r_P4*cos(alpha_P5), r_P4*-sin(alpha_P5)]

        # Points for s == 1 case
        P6 = [r_ri*cos(alpha_P5), r_ri*-sin(alpha_P5)]
        P1p5 = [P2[0] - d_rp, P2[1]]

        # Points for s > 1 case
        if s > 1:
            alpha_notch = (alpha_rm - s*alpha_rs) / (s-1)
            r_P6 = r_ri + d_ri + d_rs
            P6 = [r_P6*cos(alpha_P5), r_P6*-sin(alpha_P5)]
            
            alpha_P7 = alpha_P5 + alpha_notch
            P7 = [r_P6*cos(alpha_P7), r_P6*-sin(alpha_P7)]
            P8 = [r_P4*cos(alpha_P7), r_P4*-sin(alpha_P7)]
            
            # P9 and P10 for the last segment
            P9 = [cos(alpha_rs)*P8[0] + sin(alpha_rs)*P8[1],
                       -sin(alpha_rs)*P8[0] + cos(alpha_rs)*P8[1]]
            P10 = [r_ri*cos(alpha_rp), r_ri*-sin(alpha_rp)]
        else:
            # For s == 1, set these to None or default values
            P7 = None
            P8 = None
            P9 = None
            P10 = None
            alpha_notch = None

        drawer.getSketch(self.name, self.color)

        list_segments = []
        if s == 1:
            if alpha_rm >= alpha_rp*0.9800:
                logger = logging.getLogger(__name__)
                logger.info('Non-NOTCHED ROTOR IS USED.')
                logger.info('alpha_P5 is %s, %s deg', alpha_P5, alpha_P5/math.pi*180)
                if bool_draw_whole_model:
                    list_segments += drawer.drawArc([0,0], P1, [-P1[0], P1[1]])
                    list_segments += drawer.drawArc([0,0], [-P1[0], P1[1]], P1)
                    list_segments += drawer.drawArc([0,0], P1p5, [-P1p5[0], P1p5[1]])
                    list_segments += drawer.drawArc([0,0], [-P1p5[0], P1p5[1]], P1p5)
                else:
                    list_segments += drawer.drawLine(P1, P1p5)
                    list_segments += drawer.drawArc([0,0], P5, P1p5)
                    list_segments += drawer.drawLine(P5, P6)
                    list_segments += drawer.drawArc([0,0], P6, P1)
            else:
                if bool_draw_whole_model:
                    def iPark(P, theta):
                        return [P[0]*math.cos(theta)+P[1]*-math.sin(theta), P[0]*math.sin(theta)+P[1]*math.cos(theta)]
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
        else:
            if bool_draw_whole_model == True:
                raise Exception('NOT IMPLEMENTED for s>1.')

            P7 = self.P7
            P8 = self.P8
            P9 = self.P9
            P10 = self.P10
            alpha_notch = self.alpha_notch
            alpha_rs = self.alpha_rs

            if alpha_rm >= alpha_rp*0.9800: # no inter-pole notch
                P1p5 = [self.mm_r_ri + self.mm_d_ri, 0]
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
            
            # Rotate points for additional segments
            current_P4 = P4
            current_P5 = P5
            current_P6 = P6
            current_P7 = P7
            current_P8 = P8
            
            for _ in range(1, s-1):
                alpha_temp = (alpha_rs + alpha_notch)

                current_P4 = current_P8
                current_P5 = [cos(alpha_temp)*current_P5[0] + sin(alpha_temp)*current_P5[1],
                             -sin(alpha_temp)*current_P5[0] + cos(alpha_temp)*current_P5[1]]
                current_P6 = [cos(alpha_temp)*current_P6[0] + sin(alpha_temp)*current_P6[1],
                             -sin(alpha_temp)*current_P6[0] + cos(alpha_temp)*current_P6[1]]
                current_P7 = [cos(alpha_temp)*current_P7[0] + sin(alpha_temp)*current_P7[1],
                             -sin(alpha_temp)*current_P7[0] + cos(alpha_temp)*current_P7[1]]
                current_P8 = [cos(alpha_temp)*current_P8[0] + sin(alpha_temp)*current_P8[1],
                             -sin(alpha_temp)*current_P8[0] + cos(alpha_temp)*current_P8[1]]

                list_segments += drawer.drawArc([0,0], current_P5, current_P4)
                list_segments += drawer.drawLine(current_P5, current_P6)
                list_segments += drawer.drawArc([0,0], current_P7, current_P6)
                list_segments += drawer.drawLine(current_P7, current_P8)

            list_segments += drawer.drawArc([0,0], P9, current_P8)
            list_segments += drawer.drawLine(P9, P10)
            list_segments += drawer.drawArc([0,0], P10, P1)

        innerCoord = ( 0.5*(P1[0]+P4[0]), 0.5*(P1[1]+P4[1]))

        self.list_region=[list_segments]
        
        # Pass point coordinates to drawer for frontend visualization
        if not hasattr(drawer, 'visualization_points'):
            drawer.visualization_points = {}
        points_dict = {'P1': P1, 'P2': P2, 'P3': P3, 'P4': P4, 'P5': P5, 'P6': P6, 'P1p5': P1p5}
        if s > 1:
            points_dict.update({'P7': P7, 'P8': P8, 'P9': P9, 'P10': P10})
        drawer.visualization_points[self.name] = points_dict
        
        return {'innerCoord': innerCoord, 'list_regions':[list_segments], 'mirrorAxis': None}

class CrossSectInnerNotchedMagnet(object):
    def __init__(self, 
                    name = 'Notched Rotor',
                    color = '#0BA0E2',
                    rotorCore = None,
                    ):
        self.name = name
        self.color = color
        self.rotorCore = rotorCore

    def draw(self, drawer, bool_re_evaluate=False, bool_draw_whole_model=False):
        """Calculate all point coordinates based on the rotorCore."""
        d_pm     = self.rotorCore.mm_d_pm
        alpha_rm = self.rotorCore.deg_alpha_rm * math.pi/180
        alpha_rs = self.rotorCore.deg_alpha_rs * math.pi/180
        r_ri     = self.rotorCore.mm_r_ri
        d_ri     = self.rotorCore.mm_d_ri
        d_rp     = self.rotorCore.mm_d_rp
        d_rs     = self.rotorCore.mm_d_rs
        p        = self.rotorCore.p
        s        = self.rotorCore.s
        alpha_rp = 2*math.pi/(2*p) # pole span

        # rotor inter-pole notch being too small
        if alpha_rm >= alpha_rp*0.9800:
            logger = logging.getLogger(__name__)
            logger.info('FULL POLE PITCH MAGNET IS USED.')
            alpha_rm = alpha_rp

        if abs(alpha_rp - alpha_rm) <= 2 * math.pi/180:
            alpha_rm = alpha_rp
            if s == 1:
                alpha_rs = alpha_rm
            else:
                logger = logging.getLogger(__name__)
                logger.warning('s=%d: This is not tested. For now it simply assumes the iron notch between poles becomes the iron notch between the segments of one pole.', s)
            logger = logging.getLogger(__name__)
            logger.warning('[class CrossSectInnerNotchedMagnet] Magnet is fully spanned.')

        if d_pm + 2*EPS < d_rp:
            logger = logging.getLogger(__name__)
            logger.warning('[class CrossSectInnerNotchedMagnet] Detect d_rp is too close to d_pm. To avoid small line entity error in JMAG, set d_pm equal to d_rp because rotor core is plotted already.')
            raise ExceptionBadDesign('[Error] Magnet depth d_pm is too close to inter-pole notch depth d_rp.')

        P1 = [r_ri, 0]

        r_P2 = r_ri + d_ri + d_rp
        P2 = [r_P2, 0]

        alpha_P3 = alpha_rp - alpha_rm
        r_P4 = r_ri + d_ri
        P4 = [r_P4*cos(alpha_P3), r_P4*-sin(alpha_P3)]

        P3_extra = [(r_P4+d_pm)*cos(alpha_P3), (r_P4+d_pm)*-sin(alpha_P3)]

        alpha_P5 = alpha_P3 + alpha_rs
        P5 = [r_P4*cos(alpha_P5), r_P4*-sin(alpha_P5)]

        if s > 1:
            alpha_notch = (alpha_rm - s*alpha_rs) / (s-1)
        else:
            alpha_notch = None
            
        P6_extra = [(r_P4+d_pm)*cos(alpha_P5), (r_P4+d_pm)*-sin(alpha_P5)]

        Rout = r_P4+d_pm
        Rin  = r_P4
        mm2_magnet_area = alpha_rm/alpha_rp  *  math.pi*(Rout**2 - Rin**2)
        logger = logging.getLogger(__name__)
        logger.info('Magnet area in total is %g mm^2', mm2_magnet_area)

        if False == bool_re_evaluate:
            drawer.getSketch(self.name, self.color)

        if bool_re_evaluate:
            return mm2_magnet_area

        list_regions = []
        list_segments = []

        if bool_draw_whole_model:
            def iPark(P, theta):
                return [P[0]*math.cos(theta)+P[1]*-math.sin(theta), P[0]*math.sin(theta)+P[1]*math.cos(theta)]
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

        # Rotate points for additional segments
        current_P3_extra = P3_extra
        current_P4 = P4
        current_P5 = P5
        current_P6_extra = P6_extra
        
        for _ in range(s-1):
            alpha_temp = (alpha_rs + alpha_notch) if alpha_notch is not None else alpha_rs

            current_P3_extra = [cos(alpha_temp)*current_P3_extra[0] + sin(alpha_temp)*current_P3_extra[1],
                               -sin(alpha_temp)*current_P3_extra[0] + cos(alpha_temp)*current_P3_extra[1]]
            current_P4 = [cos(alpha_temp)*current_P4[0] + sin(alpha_temp)*current_P4[1],
                         -sin(alpha_temp)*current_P4[0] + cos(alpha_temp)*current_P4[1]]
            current_P5 = [cos(alpha_temp)*current_P5[0] + sin(alpha_temp)*current_P5[1],
                         -sin(alpha_temp)*current_P5[0] + cos(alpha_temp)*current_P5[1]]
            current_P6_extra = [cos(alpha_temp)*current_P6_extra[0] + sin(alpha_temp)*current_P6_extra[1],
                               -sin(alpha_temp)*current_P6_extra[0] + cos(alpha_temp)*current_P6_extra[1]]

            list_segments += drawer.drawLine(current_P3_extra, current_P4)
            list_segments += drawer.drawArc([0,0], current_P5, current_P4)
            list_segments += drawer.drawLine(current_P5, current_P6_extra)
            list_segments += drawer.drawArc([0,0], current_P6_extra, current_P3_extra)

            list_regions.append(list_segments)
            list_segments = []

        innerCoord = ( 0.5*(P4[0]+P6_extra[0]), 0.5*(P4[1]+P6_extra[1]))
        self.list_region = list_regions
        
        # Pass point coordinates to drawer for frontend visualization
        if not hasattr(drawer, 'visualization_points'):
            drawer.visualization_points = {}
        drawer.visualization_points[self.name] = {
            'P1': P1, 'P2': P2, 'P3_extra': P3_extra, 'P4': P4, 'P5': P5, 'P6_extra': P6_extra
        }
        
        return {'innerCoord': innerCoord, 'list_regions':list_regions, 'mirrorAxis': None}

class CrossSectSleeve(object):
    def __init__(self, 
                    name = 'Sleeve',
                    color = '#11E322',
                    mm_r_ri=5,
                    mm_d_ri=5,
                    mm_d_pm=3,
                    p=4,
                    d_sleeve=1
                    ):
        self.name = name
        self.color = color

        self.mm_r_ri = mm_r_ri
        self.mm_d_ri = mm_d_ri
        self.mm_d_pm = mm_d_pm
        self.p = p
        self.d_sleeve = d_sleeve

    def draw(self, drawer):
        """Calculate all point coordinates."""
        r_ri  = self.mm_r_ri
        d_ri  = self.mm_d_ri
        d_pm  = self.mm_d_pm
        p     = self.p

        r_or = r_ri + d_ri + d_pm 
        d_sleeve = self.d_sleeve

        P1 = [r_or, 0]
        P2 = [r_or+d_sleeve, 0]

        P3 = [cos(math.pi/p)*P1[0] + sin(math.pi/p)*P1[1],
                  -sin(math.pi/p)*P1[0] + cos(math.pi/p)*P1[1]]
        P4 = [cos(math.pi/p)*P2[0] + sin(math.pi/p)*P2[1],
                  -sin(math.pi/p)*P2[0] + cos(math.pi/p)*P2[1]]

        drawer.getSketch(self.name, self.color)

        list_regions = []
        list_segments = []
        list_segments += drawer.drawLine(P1, P2)
        list_segments += drawer.drawArc([0,0], P4, P2)
        list_segments += drawer.drawLine(P4, P3)
        list_segments += drawer.drawArc([0,0], P3, P1)

        list_regions.append(list_segments)
        list_segments = []

        innerCoord = ( 0.5*(P1[0]+P2[0]), 0.5*(P1[1]+P2[1]))
        self.list_region = list_regions
        
        # Pass point coordinates to drawer for frontend visualization
        if not hasattr(drawer, 'visualization_points'):
            drawer.visualization_points = {}
        drawer.visualization_points[self.name] = {
            'P1': P1, 'P2': P2, 'P3': P3, 'P4': P4
        }
        
        return {'innerCoord': innerCoord, 'list_regions':list_regions, 'mirrorAxis': None}

class CrossSectShaft(object):
    def __init__(self, 
                    name = 'Shaft',
                    color = '#0EE0E2',
                    rotorCore = None,
                    ):
        self.name = name
        self.color = color
        self.rotorCore = rotorCore

    def draw(self, drawer, bool_draw_whole_model=False):
        """Calculate all point coordinates."""
        r_ri = self.rotorCore.mm_r_ri

        P1 = [r_ri, 0]
        NP1 = [-r_ri, 0]

        drawer.getSketch(self.name, self.color)

        list_regions = []
        list_segments = []
        list_segments += drawer.drawArc([0,0], NP1, P1)
        list_segments += drawer.drawArc([0,0], P1, NP1)

        list_regions.append(list_segments)
        list_segments = []

        innerCoord = ( 0, 0 )

        self.list_region = list_regions
        
        # Pass point coordinates to drawer for frontend visualization
        if not hasattr(drawer, 'visualization_points'):
            drawer.visualization_points = {}
        drawer.visualization_points[self.name] = {
            'P1': P1, 'NP1': NP1
        }
        
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
        rotorCore = CrossSectInnerNotchedRotor( name = 'NotchedRotor',
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

    list_regions = rotorCore.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = rotorCore.p*2
    region1 = toolJd.prepareSection(list_regions)
    
    if True:
        notched_magnet = CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
                                                      color = '#0E001E',
                                                      rotorCore = rotorCore
                                                    )

    list_regions = notched_magnet.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = rotorCore.p*2
    region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')
