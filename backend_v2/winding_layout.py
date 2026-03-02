import logging

def infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(layer_X_phases, coil_pitch):
    return layer_X_phases[-coil_pitch:] + layer_X_phases[:-coil_pitch]
def infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(layer_X_signs, coil_pitch):
    temp = layer_X_signs[-coil_pitch:] + layer_X_signs[:-coil_pitch]
    return [('-' if el == '+' else '+') for el in temp]
def infer_Y_layer_grpAC_from_X_layer_and_coil_pitch_y(grouping_AC, coil_pitch):
    return grouping_AC[-coil_pitch:] + grouping_AC[:-coil_pitch]
    # >>> a = [1,2,3,4,5,6,7,8,0]
    # >>> a[-7:] + a[:-7]
    # [3, 4, 5, 6, 7, 8, 0, 1, 2]
from collections import OrderedDict

class winding_layout_v3(OrderedDict):
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        if name.startswith('_'):
            super().__setattr__(name, value)
        else:
            self[name] = value

    def __init__(self, Qs, p, ps=None, coil_pitch_y=None, pr=None, m=3, Wrap_Around=None):
        ''' Naming convention:
        # right layer = 1st layer = X layer = torque layer for separate winding
        # left layer  = 2nd layer = Y layer = suspension layer for separate winding
        '''
        super().__init__()
        self.bool_distributed_or_concentrated = (coil_pitch_y != 1)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Combined Winding for Bearingless Motor
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

        if Qs == 12 \
        and p == 4 \
        and ps == 5:

            self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W']
            self.layer_X_signs  = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if Qs == 12 \
            and p == 5 \
            and ps == 4:

                self.layer_X_phases = ['U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U']
                self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
                self.coil_pitch_y   = coil_pitch_y
                self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
                self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

                self.grouping_AC            = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0]
                self.number_parallel_branch = 2
                self.number_winding_layer   = 2

                self.bool_3PhaseCurrentSource = False
                self.CommutatingSequenceD = 1
                self.CommutatingSequenceB = 0

        if Qs == 12 \
            and p == 5 \
            and ps == 1:

                self.layer_X_phases = ['U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U']
                self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
                self.coil_pitch_y   = coil_pitch_y
                self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
                self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

                self.grouping_AC            = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
                self.number_parallel_branch = 2
                self.number_winding_layer   = 2

                self.bool_3PhaseCurrentSource = False
                self.CommutatingSequenceD = 1
                self.CommutatingSequenceB = 0
        
        if Qs == 12 \
            and p == 4 \
            and ps == 1 \
            and coil_pitch_y == 1:

                self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W']
                self.layer_X_signs  = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
                self.coil_pitch_y   = coil_pitch_y
                self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
                self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

                self.grouping_AC            = [0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1]
                self.number_parallel_branch = 2
                self.number_winding_layer   = 2

                self.bool_3PhaseCurrentSource = False
                self.CommutatingSequenceD = 1
                self.CommutatingSequenceB = 0

        if Qs == 12 \
            and p == 2 \
            and ps == 1:

                self.layer_X_phases = ['U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V']
                self.layer_X_signs  = ['+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-']
                self.coil_pitch_y   = coil_pitch_y
                self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
                self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

                self.grouping_AC            = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0]
                self.number_parallel_branch = 2
                self.number_winding_layer   = 2

                self.bool_3PhaseCurrentSource = False
                self.CommutatingSequenceD = 1
                self.CommutatingSequenceB = 0

        if Qs == 12 \
        and p == 1 \
        and ps == 2:

            self.layer_X_phases = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if Qs == 12 \
        and p == 4 \
        and ps == 5 \
        and coil_pitch_y == 1:

            self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W']
            self.layer_X_signs  = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

            self.kd1 = 1.0
            self.kp1 = 0.866

        if Qs == 18 \
        and p == 2 \
        and ps == 1:

            self.layer_X_phases = ['U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V', 'U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if Qs == 24 \
        and p == 4 \
        and ps == 1 \
        and coil_pitch_y == 9:
            self.layer_X_phases = ['U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V']
            self.layer_X_signs  = ['+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if Qs == 24 \
        and p == 2 \
        and ps == 1:
            self.layer_X_phases = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V'] # ExampleQ24p2m3ps1: torque winding outer layer
            self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.coil_pitch_y = int(Qs/(2*p)) if coil_pitch_y is None else coil_pitch_y # Y layer can be inferred from coil pitch and X layer diagram
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            # grouping AC is valid for the X layer signs, the Y layer can always be inferred from the X layer signs and the coil pitch.
            self.grouping_AC    = [  0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0] # 只取决于1s tlayer/outer layer/right layer的反相情况
            self.number_parallel_branch = 2.
            self.number_winding_layer = 2 # for torque winding and this means there could be a short pitch

            # Excitation options
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.CommutatingSequenceD = 1 # D stands for Drive winding (i.e., torque winding)
            self.CommutatingSequenceB = 0 # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field

        if Qs == 6 \
            and p == 4 \
            and ps == 1:

                self.layer_X_phases = ['U', 'W', 'V', 'U', 'W', 'V']
                self.layer_X_signs  = ['+', '+', '+', '+', '+', '+']
                self.coil_pitch_y   = coil_pitch_y
                self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
                self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

                self.grouping_AC            = [0, 0, 1, 1, 1, 0]
                self.number_parallel_branch = 2
                self.number_winding_layer   = 2

                self.bool_3PhaseCurrentSource = False
                self.CommutatingSequenceD = 1
                self.CommutatingSequenceB = 0
        
        if Qs == 6 \
            and p ==  2 \
            and ps == 1:
                self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W']
                self.layer_X_signs  = ['+', '+', '+', '+', '+', '+']
                self.coil_pitch_y   = coil_pitch_y
                self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
                self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

                self.grouping_AC            = [0, 1, 0, 1, 0, 1]
                self.number_parallel_branch = 2
                self.number_winding_layer   = 2

                self.bool_3PhaseCurrentSource = False
                self.CommutatingSequenceD = 1
                self.CommutatingSequenceB = 0

        if Qs == 6 \
        and p == 5 \
        and ps == 4:

            self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W']
            self.layer_X_signs  = ['+', '-', '+', '-', '+', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 1, 0, 1, 0, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if Qs == 6 \
        and p == 1 \
        and ps == 2:

            self.layer_X_phases = ['U', 'W', 'V', 'U', 'W', 'V']
            self.layer_X_signs  = ['+', '-', '+', '-', '+', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 1, 0, 1, 0, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # End of winding definition
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        try: 
            self.coil_pitch_y
            self.distributed_or_concentrated = False if abs(self.coil_pitch_y) == 1 else True

            # below is valid for PMSM only

            q = SPP = Qs / (2*p * m)

            # 寻找初始d轴激励角度
            if q%1 == 0:
                # integral slot

                phase_U_belt = self.layer_X_phases[:int(q)]
                number_of_U = sum([1 for el in phase_U_belt if el =='U'])
                if number_of_U<q:
                    # print('目前只有%d个字母U，需要寻找一共q(=%d)个字母U。'%(number_of_U, q))
                    for phase_U_starting_slot_number in range(-1, -int(q)-1, -1):
                        new_phase_U_belt = self.layer_X_phases[phase_U_starting_slot_number:] + phase_U_belt
                        number_of_U = sum([1 for el in new_phase_U_belt if el =='U'])
                        # print(phase_U_starting_slot_number, number_of_U)
                        if number_of_U < q:
                            continue
                        else:
                            break
                else:
                    phase_U_starting_slot_number = 1 # 没有等于0的哈，等于0你得以槽的中线作图绘制定子铁芯；现在的代码是以齿的中线作图绘制定子铁芯的哦。

                self.deg_winding_U_phase_phase_axis_angle = 360/Qs*0.5 * (phase_U_starting_slot_number + self.coil_pitch_y+(SPP-1) )
                logger = logging.getLogger(__name__)
                logger.info('[wily] self.deg_winding_U_phase_phase_axis_angle=%s', self.deg_winding_U_phase_phase_axis_angle)
                logger.info('[wily] q = SPP = %s', SPP)
            else:
                # This clause includes fractional slot winding with an SPP value below 1.

                for index, UVW in enumerate(self.layer_X_phases):
                    if UVW == 'V':
                        phase_V_starting_slot_number = index+1
                        break
                for index, UVW in enumerate(self.layer_X_phases):
                    if UVW == 'W':
                        phase_W_starting_slot_number = index+1
                        break
                # print(self.layer_X_phases)
                # print(phase_V_starting_slot_number, phase_W_starting_slot_number, phase_W_starting_slot_number-phase_V_starting_slot_number)

                尾 = phase_W_starting_slot_number-1 + self.coil_pitch_y
                头 = phase_V_starting_slot_number
                deg_winding_V_phase_phase_axis_angle = 360/Qs * (尾 - 头) + 360/Qs*0.5
                相邻的属于同一相的线圈的个数 = phase_W_starting_slot_number - phase_V_starting_slot_number

                if True:
                    # The winding axis is defined in JMAG.preProcess where the upper/lower coils are assigned to phase connections.
                    self.deg_winding_U_phase_phase_axis_angle = None
                else:
                    if q<1:
                        # fractional slot and q<1
                        self.deg_winding_U_phase_phase_axis_angle = deg_winding_V_phase_phase_axis_angle - 360/Qs*相邻的属于同一相的线圈的个数
                    else:
                        # fractional slot and q>1
                        self.deg_winding_U_phase_phase_axis_angle = deg_winding_V_phase_phase_axis_angle - 360/Qs*相邻的属于同一相的线圈的个数
                        msg = '[winding_layout.py] [Warning] This case (q=%g) is not thought thorough, so you must inspect the initial excitation angle and initial rotor position manually to make sure it is id=0 control.'%(q)
                        logger = logging.getLogger(__name__)
                        logger.warning(msg)
                        if q>2:
                            raise Exception(msg)
                    logger = logging.getLogger(__name__)
                    logger.info('[wily] self.deg_winding_U_phase_phase_axis_angle=%s, deg_winding_V_phase_phase_axis_angle=%s', self.deg_winding_U_phase_phase_axis_angle, deg_winding_V_phase_phase_axis_angle)
                    logger.info('[wily] q = SPP = %s', SPP)

                # print(self.deg_winding_U_phase_phase_axis_angle, deg_winding_V_phase_phase_axis_angle, 尾, 头, 相邻的属于同一相的线圈的个数)
                # quit()

                    # 2019/12/25: Q18, p2, ps1 case is verified to be okay with q=1.5.
            # quit()
        except:
            raise Exception(f'Error: This winding is not implemented to get deg_winding_U_phase_phase_axis_angle: Qs,p,ps,coil_pitch_y = {Qs,p,ps,coil_pitch_y}', Wrap_Around, DPNV_or_SEPA)


        # 这是实际在pre_procee中调用的字典
        self.dict_coil_connection = {'layer X phases': self.layer_X_phases, 'layer X signs':self.layer_X_signs,   # 这里的命名规则是按照seprate winding的情况来的。
                                     'layer Y phases': self.layer_Y_phases, 'layer Y signs':self.layer_Y_signs}   # 这里的命名规则是按照seprate winding的情况来的。


from pylab import np, plt, fft, linspace
import math
import scipy.integrate as integrate

def nextpow2(L):
    n = 0
    while 2**n < L:
        n += 1
    return n

def periodic2pi(x):
    # make x periodic positive
    while(x<0):
        x += 2*math.pi

    # make it within 2*pi
    if int(x / (2*math.pi)) != 0:
        x = x % (2*math.pi)

    return x

def segmented_func(x, lst_x, lst_y):
    ''' The value in lst_x must increase monotonously.
        根据判断，x在lst_x中的位置，返回对应的阶梯函数值。
    '''
    x = periodic2pi(x)
    
    if len(lst_x) != len(lst_y):
        raise Exception('Error: lst_x & y must has same length.')

    for i in range(0,len(lst_x)):
        if x > lst_x[i]: # x 超过了lst_x[i]的位置，
            if i == len(lst_x)-1: # last i。比最后一个位置lst_x[-1]大，比2pi小。
                return lst_y[-1]

        else: # x 比lst_x[i]小，但是比lst_x[i-1]大。
            if i-1 < 0: # first i。比第一个位置lst_x[0]小，比0大。
                return lst_y[-1] # 注意lst_y[0]的值不是0，但是lst_y[-1]是0.
            else: # 普通情况：
                return lst_y[i-1]

class PhaseWinding:
    '''[PhaseWinding]
    
    [The winding function and turn function are consistent with Lipo's 2012 book.]
    '''
    def __init__(self, Qs, m, turns_per_slot, ox_distribution_phase_U, desc_type = 'STEP_DISTRIBUTED'):
        if desc_type == 'STEP_DISTRIBUTED':
            pass
        else:
            raise Exception('Not implemented for other type of descripsion.')

        self.slot_per_phase = Qs / m # Only upper layer counts for one slot.
        self.turns_per_slot = turns_per_slot
        self.degree_between_slots = 360.0 / Qs # mechanical deg.
        self.radian_between_slots = self.degree_between_slots / 180.0 * math.pi # mechanical rad.
        self.ox_distribution_phase_U = ox_distribution_phase_U
        print(self.ox_distribution_phase_U)

        # turn function 
        self.setTurnFuncObject(self.ox_distribution_phase_U) # define self.turn_func
        # <turn function> defined by Lipo 2012
        self.avg_val_of_turn_func = integrate.quad(self.turn_func, 0, 2*math.pi)[0] / (2*math.pi) # [1] is error of integration
        # winding function
        self.winding_func = lambda x: self.turn_func(x)-self.avg_val_of_turn_func

        # get self.sym_begin_pos 
        self.setSymPos()
        # symmetric turn function
        self.sym_turn_func = lambda x: self.turn_func(x + self.sym_begin_pos)
        # symmetric winding function
        self.sym_winding_func = lambda x: self.winding_func(x + self.sym_begin_pos)

    def setTurnFuncObject(self, ox_distribution_phase_U):
        '''  ['x','x','x','x',
              'n','n','n','n',
              'o','o','o','o',
              'n','n','n','n']
              非n地方才记录，x增，o减。
        '''
        x_val, y_val, lst_x, lst_y = 0, 0, [], []
        for el in ox_distribution_phase_U:
            x_val += self.radian_between_slots
            if el == 'x':
                y_val += self.turns_per_slot
                lst_y.append(y_val)
                lst_x.append(x_val)
            elif el == 'xx':
                y_val += 2* self.turns_per_slot                
                lst_y.append(y_val)
                lst_x.append(x_val)
            elif el == 'o':
                y_val -= self.turns_per_slot                
                lst_y.append(y_val)
                lst_x.append(x_val)
            elif el == 'oo':
                y_val -= 2 * self.turns_per_slot                
                lst_y.append(y_val)
                lst_x.append(x_val)
            else:
                pass
        self.lst_x = lst_x
        self.lst_y = lst_y
        self.turn_func = lambda x: segmented_func(x, self.lst_x, self.lst_y)

    def setSymPos(self, index=2):
        '''
            可选择两种不同的对称位置。用作原始turn_func和wind_func的偏移。举例：
            和lst_x一一对应的lst_y:
                50
                100
                150
                200      <- pos2 @ max y 
                150x
                100x
                50x
                0        <- pos1 @ min y
        '''
        # symmetrical pos_2 = pos_200 + (pos_150x - pos_200) / 2
        index = self.lst_y.index(max(self.lst_y))
        self.sym_begin_pos_2 = self.lst_x[index] + (self.lst_x[index+1] - self.lst_x[index]) / 2.

        # symmetrical pos_1
        self.sym_begin_pos_1 = self.sym_begin_pos_2 + math.pi

        if index == 1:
            self.sym_begin_pos = self.sym_begin_pos_1
        else:
            self.sym_begin_pos = self.sym_begin_pos_2

    def plot2piFft(self, func, Fs, L):
        ''' Fs is the sampling freq. 
            L is length of signal list.
            This plot is for a func that has period of 2pi.

            If you found the time domain wave is not very accurate,
            that is because you set too small Fs, which leads to
            to big step Ts.
        '''
        base_freq = 1.0/(2*math.pi) #频域横坐标除以基频，即以基频为单位，此处的基频为 2*pi rad/s
        Ts = 1.0/Fs
        t = [el*Ts for el in range(0,L)]
        x = [func(el) for el in t]

        # https://www.ritchievink.com/blog/2017/04/23/understanding-the-fourier-transform-by-example/

        # 小明给的代码：
        # sampleF = Fs
        # print('小明：')
        # for f, Y in zip(
        #                 np.arange(0, len(x)*sampleF,1) * 1/len(x) * sampleF, 
        #                 np.log10(np.abs(np.fft.fft(x) / len(x))) 
        #              ):
            # print('\t', f, Y)


        L_4pi = int(4*math.pi / Ts) +1 # 画前两个周期的
        
        self.fig_plot2piFft = plt.figure(7)
        plt.subplot(211)
        plt.plot(t[:L_4pi], x[:L_4pi])
        #title('Signal in Time Domain')
        #xlabel('Time / s')
        #ylabel('x(t)')
        plt.title('Winding Function')
        plt.xlabel('Angular location along air gap [mech. rad.]')
        plt.ylabel('Current Linkage by unit current [Ampere]')

        NFFT = 2**nextpow2(L)
        print('NFFT =', NFFT, '= 2^%g' % (nextpow2(L)), '>= L =', L)
        y = fft(x,NFFT) # y is a COMPLEX defined in numpy
        Y = [2 * el.__abs__() / L for el in y] # /L for spectrum aplitude consistent with actual signal. 2* for single-sided. abs for amplitude.
        f = Fs/2.0/base_freq*linspace(0,1,int(NFFT/2+1)) # unit is base_freq Hz
        #f = Fs/2.0*linspace(0,1,NFFT/2+1) # unit is Hz

        plt.subplot(212)
        plt.plot(f, Y[0:int(NFFT/2+1)])
        plt.title('Single-Sided Amplitude Spectrum of x(t)')
        plt.xlabel('Frequency divided by base_freq [base freq * Hz]')
        #plt.ylabel('|Y(f)|')
        plt.ylabel('Amplitude [1]')
        plt.xlim([0,50])
        # plt.show()

    def plotFuncObj(self, func):
        x = np.arange(0, 2*math.pi, 0.5/180*math.pi)
        y = [func(el) for el in x]
        x = [el/math.pi for el in x] # x for plot. unit is pi
        
        self.fig_plotFuncObj = plt.figure()
        ax = plt.subplot(111) #注意:一般都在ax中设置,不在plot中设置
        ax.plot(x, y)

        xmajorLocator   = plt.MultipleLocator(0.25) #将x主刻度标签设置为20的倍数
        xmajorFormatter = plt.FormatStrFormatter('%.2fπ') #设置x轴标签文本的格式
        ax.xaxis.set_major_locator(xmajorLocator)  
        ax.xaxis.set_major_formatter(xmajorFormatter)

        ##matplotlib.pyplot.minorticks_on()
        ##xminorLocator   = MultipleLocator(0.25)
        ##xminorFormatter = FormatStrFormatter(u'%.2fπ')
        ##ax.xaxis.set_minor_locator(xminorLocator)  
        ##ax.xaxis.set_minor_formatter(xminorFormatter)

        plt.xlabel('Angular location along the gap [mech. rad.]')
        plt.ylabel('Turns of winding [1]')
        # plt.title('Turn Function or Winding Function')
        plt.grid(True) # or ax.grid(True)
        # plt.gcf().savefig("turn_function.png")
        # plt.show()

if __name__ == '__main__':

    # wily = winding_layout_v2(DPNV_or_SEPA=False, Qs=24, p=2, ps=1)
    # wily = winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=2, ps=1, coil_pitch_y=6)
    wily = winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=1, ps=2, coil_pitch_y=9)

    zQ = 100 # number of conductors/turns per slot
    turns_per_layer = zQ / wily.number_winding_layer
    U_phase = PhaseWinding(wily.Qs, wily.m, turns_per_layer, wily.ox_distribution_phase_U)
    U_phase.plotFuncObj(U_phase.winding_func)
    U_phase.plot2piFft(U_phase.winding_func, Fs=1/(2*math.pi/3600), L=65536*2**4) # 采样频率：在2pi的周期内取720个点

    plt.show()
