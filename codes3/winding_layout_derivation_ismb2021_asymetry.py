ABCDEFGHIJKLMNOPQRSTUVWXYZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
import os
output_dir = os.path.dirname(__file__) + r'/../_wily/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
print('Output directory is:', output_dir)

# Globals for drawing 
RADIUS = 8
BELT_BIAS = 5 # deg. elec.
distance_between_label_layers = 1.33 # 1.0, 1.25
angle_between_arrow_and_label_star_of_slots = 6 # deg
angle_between_arrow_and_label = 4 # deg
PLOT_SPACING = 35 # 35


from pylab import np
def limit_to_360_deg(PHI):
    while PHI < 0:
        PHI += 360
    while PHI >= 360:
        PHI -= 360
    return PHI
def belong_to_which_phase_belt(PHI, phase_belt):

    PHI = limit_to_360_deg(PHI)

    if phase_belt == 60:
        if PHI   <= phase_belt*0.5 + BELT_BIAS or PHI > 360 - phase_belt*0.5 + BELT_BIAS:
            return 'A' # 'u+'
        elif 180 - phase_belt*0.5 + BELT_BIAS < PHI <= 180 + phase_belt*0.5 + BELT_BIAS: 
            return 'a' # 'u-'

        elif 120 - phase_belt*0.5 + BELT_BIAS < PHI <= 120 + phase_belt*0.5 + BELT_BIAS: 
            return 'B' # 'v+'
        elif 300 - phase_belt*0.5 + BELT_BIAS < PHI <= 300 + phase_belt*0.5 + BELT_BIAS: 
            return 'b' # 'v-'

        elif 240 - phase_belt*0.5 + BELT_BIAS < PHI <= 240 + phase_belt*0.5 + BELT_BIAS: 
            return 'C' # 'w+'
        elif 60  - phase_belt*0.5 + BELT_BIAS < PHI <= 60  + phase_belt*0.5 + BELT_BIAS: 
            return 'c' # 'w-'
        else:
            raise Exception('Unexpected PHI=%g'%(PHI))
    elif phase_belt == 120:
        if PHI   <= phase_belt*0.5 + BELT_BIAS or PHI > 360 - phase_belt*0.5 + BELT_BIAS:
            return 'A' # 'u+'
        elif 120 - phase_belt*0.5 + BELT_BIAS < PHI <= 120 + phase_belt*0.5 + BELT_BIAS: 
            return 'B' # 'v+'
        elif 240 - phase_belt*0.5 + BELT_BIAS < PHI <= 240 + phase_belt*0.5 + BELT_BIAS: 
            return 'C' # 'w+'
        else:
            raise Exception('Unexpected PHI=%g'%(PHI))
    else:
        phase_name_positive_list = ABCDEFGHIJKLMNOPQRSTUVWXYZ
        phase_name_negative_list = ABCDEFGHIJKLMNOPQRSTUVWXYZ.lower()

        i = 0
        # print('PHASE BELT:', i, ':', PHI)
        if PHI   <= phase_belt*0.5 + BELT_BIAS or PHI > 360 - phase_belt*0.5 + BELT_BIAS:
            return phase_name_positive_list[0] # 'u+'
        elif 180 - phase_belt*0.5 + BELT_BIAS < PHI <= 180 + phase_belt*0.5 + BELT_BIAS: 
            return phase_name_negative_list[0] # 'u-'
        else:
            while True:
                i += 1
                # print('PHASE BELT:', i, ':', PHI, end='in')
                # print('[', i*2*phase_belt       - phase_belt*0.5 + BELT_BIAS, ',', i*2*phase_belt     + phase_belt*0.5 + BELT_BIAS, ']')
                if   i*2*phase_belt     - phase_belt*0.5 + BELT_BIAS < PHI <= i*2*phase_belt     + phase_belt*0.5 + BELT_BIAS: 
                    return phase_name_positive_list[i]
                elif i*2*phase_belt+180 - phase_belt*0.5 + BELT_BIAS < PHI <= i*2*phase_belt+180 + phase_belt*0.5 + BELT_BIAS: 
                    return phase_name_negative_list[i]
                if i>100:
                    raise Exception('Dead loop.')
        # raise Exception('Unexpected phase_belt: %g'%(phase_belt))
def belong_to_band(LB, UB, PHI):

    LB = limit_to_360_deg(LB)
    UB = limit_to_360_deg(UB)
    PHI = limit_to_360_deg(PHI)

    if LB > UB:
        LB -= 360

    if LB <= PHI <= UB:
        return True
    else:
        if LB <= PHI-360 <= UB:
            return True
        else:
            return False
def angular_location(PHI, radius_bias=0):
    X = (RADIUS+radius_bias) * np.cos(PHI/180*np.pi)
    Y = (RADIUS+radius_bias) * np.sin(PHI/180*np.pi)
    return X, Y
def phase_angle_of_slot_i_at_frequency_h(slot_number, h, Q):
    # - \alpha_{i,h}^e
    # - fundamental frequency <=> h=p
    槽距角 = deg_elec_angle_between_adjacent_slots = 2*np.pi * h / Q 
    return deg_elec_angle_between_adjacent_slots * (slot_number-1) # 水平向右定义为1号槽

import PyX_Utility
from utility import gcd
def draw_star_of_slots(Q, p, m):

    电角度 = electrical_degree = 360 * p
    槽距角 = deg_elec_angle_between_adjacent_slots = 360 * p / Q 
    # 极距   = pole_pitch = Q / (2*p)
    # 节距   = coil_pitch = pole_pitch - 0 # -0, -1, -2, -3
    相带   = phase_belt = 360/(2*m) # 60 or 120 deg. elec. | /(2*m) or /m
    print('\t Number of Zones = 2*p*m = %d'%(2*p*m))
    print('\t 槽距角：\t', deg_elec_angle_between_adjacent_slots, 'deg. elec.')
    print('\t 相带：\t', phase_belt, 'deg. elec.')
    # 槽电势星形图圈数 = t = gcd(Q,p)
    # 槽电势星形图每圈箭头数 = Q / t

    u = PyX_Utility.PyX_Utility()
    connection_star_raw_dict = dict()

    # Draw phase belt
    if phase_belt == 60:
        u.pyx_draw_sector( [0,0], -30+BELT_BIAS,  30+BELT_BIAS, RADIUS+1)
        u.pyx_draw_sector( [0,0],  90+BELT_BIAS, 150+BELT_BIAS, RADIUS+1)
        u.pyx_draw_sector( [0,0], 210+BELT_BIAS, 270+BELT_BIAS, RADIUS+1)
        u.pyx_text(angular_location(  0, radius_bias=2), '$+u$', scale=1) # 'A'
        u.pyx_text(angular_location(120, radius_bias=2), '$+w$', scale=1) # 'B' NOTE THAT THE PHASE V and W are transposed!
        u.pyx_text(angular_location(240, radius_bias=2), '$+v$', scale=1) # 'C' NOTE THAT THE PHASE V and W are transposed!
        u.pyx_text(angular_location(180, radius_bias=2), '$-u$', scale=1) # 'X' 
        u.pyx_text(angular_location(300, radius_bias=2), '$-w$', scale=1) # 'Y' NOTE THAT THE PHASE V and W are transposed!
        u.pyx_text(angular_location( 60, radius_bias=2), '$-v$', scale=1) # 'Z' NOTE THAT THE PHASE V and W are transposed!
    elif phase_belt == 120:
        u.pyx_draw_sector( [0,0], -60+BELT_BIAS,  60+BELT_BIAS, RADIUS+1)
        u.pyx_draw_sector( [0,0],  60+BELT_BIAS, 180+BELT_BIAS, RADIUS+1)
        u.pyx_draw_sector( [0,0], 180+BELT_BIAS, 300+BELT_BIAS, RADIUS+1)
    else:
        for i in range(int(360 / (2*phase_belt))):
            u.pyx_draw_sector( [0,0], 2*i*phase_belt+-phase_belt/2+BELT_BIAS,  2*i*phase_belt+phase_belt/2+BELT_BIAS, RADIUS+1)
            u.pyx_text(angular_location(  2*i*phase_belt + BELT_BIAS, radius_bias=2), ABCDEFGHIJKLMNOPQRSTUVWXYZ[i], scale=1) # Add phase name labels to eaach zone
            # u.pyx_text(angular_location(180, radius_bias=2), '$-u$', scale=1)

    # Draw arrow and label
    for i in range(Q):
        # Draw arrow
        PHI = i*deg_elec_angle_between_adjacent_slots
        u.pyx_arrow(angular_location(PHI))

        # Draw label
        # distance_between_label_layers = 1.0
        u.pyx_text(angular_location(PHI+angle_between_arrow_and_label_star_of_slots, radius_bias=-(PHI//360)*distance_between_label_layers), r'\textbf{'+str(i+1)+'}', scale=0.8)

        # Draw 辅助线
        u.pyx_circle(RADIUS-(PHI//360)*1.0, bool_dashed=True, dash_list=[0,24], linewidth=0.020)
        # u.pyx_draw_sector([0,0], 0, 360, RADIUS-(PHI//360)*1.0, bool_stroke=True) # circle

        # Collect phase belt arrows for drawing connection star
        # PHI = limit_to_360_deg(PHI) // 这里的PHI不能归化到360°以内，否则会影响后面group a/c分组的判断！就是会致使dpnv_grouping_dict的生成结果是错误的。
        key = belong_to_which_phase_belt(PHI, phase_belt)
        val = (PHI, i+1)
        print(key, val)
        if key in connection_star_raw_dict:
            connection_star_raw_dict[key].append(val)
        else:
            connection_star_raw_dict[key] = [val]

        # print(i, PHI)
    # print(connection_star_raw_dict)
    # quit()
    return u, connection_star_raw_dict, phase_belt
def draw_connection_star(m, phase_belt, connection_star_raw_dict):

    u = PyX_Utility.PyX_Utility()

    if m == 3:
        u.pyx_text(angular_location(  0, radius_bias=2), '$u$', scale=1) # 'A'
        u.pyx_text(angular_location(120, radius_bias=2), '$w$', scale=1) # 'B' NOTE THAT THE PHASE V and W are transposed!
        u.pyx_text(angular_location(240, radius_bias=2), '$v$', scale=1) # 'C' NOTE THAT THE PHASE V and W are transposed!
    else:
        for i in range(int(360 / (2*phase_belt))):
            u.pyx_draw_sector( [0,0], 2*i*phase_belt+-phase_belt/2+BELT_BIAS,  2*i*phase_belt+phase_belt/2+BELT_BIAS, RADIUS+1)
            u.pyx_text(angular_location(  2*i*phase_belt + BELT_BIAS, radius_bias=2), ABCDEFGHIJKLMNOPQRSTUVWXYZ[i], scale=1) # Add phase name labels to eaach zone

    for key, val in connection_star_raw_dict.items():
        
        print('|', key, val)
        
        if key in 'abc':
            factor_reverse = -1
            phase_shift    = 180
        else: # in 'ABC'
            factor_reverse = 1
            phase_shift    = 0

        for el in val:
            # unpack
            PHI_ori = el[0]
            PHI, label = phase_shift+PHI_ori, str(factor_reverse*el[1])

            # Draw arrows
            u.pyx_arrow(angular_location(PHI))

            # The labels should drawn to the side a bit
            # PHI += angle_between_arrow_and_label * factor_reverse # switch side for reverse connected conductor
            PHI += angle_between_arrow_and_label

            # Draw labels
            if False:
                # Option 1 (looks crappier)
                u.pyx_text(angular_location(PHI+6*factor_reverse, radius_bias=-(PHI_ori//360)*1.0), label, scale=1)
            else:
                # Option 2 (looks better)
                # distance_between_label_layers = 1.25

                # radius = RADIUS-(PHI//360)*distance_between_label_layers

                if factor_reverse == -1:
                    radius = RADIUS - distance_between_label_layers
                else:
                    radius = RADIUS

                # X = (radius) * np.cos(4*factor_reverse/180*np.pi) # (radius, 0)旋转一个小角度（比如4度）
                # Y = (radius) * np.sin(4*factor_reverse/180*np.pi) # (radius, 0)旋转一个小角度（比如4度）
                X = (radius) * np.cos(4/180*np.pi) # (radius, 0)旋转一个小角度（比如4度）
                Y = (radius) * np.sin(4/180*np.pi) # (radius, 0)旋转一个小角度（比如4度）

                # X += -(PHI_ori//360)*distance_between_label_layers # 水平移动

                angle = PHI/180*np.pi
                X, Y = X*np.cos(angle) + Y*-np.sin(angle), X*np.sin(angle) + Y*np.cos(angle) # Park Trans.

                u.pyx_text((X,Y), r'\textbf{'+label+'}', scale=0.8) # T2 

                # Draw 辅助线
                u.pyx_circle(radius, bool_dashed=True, dash_list=[0,24], linewidth=0.020)
                print('draw_connection_star:::', label, PHI)
    return u
def draw_connection_star_at_another_frequency(connection_star_raw_dict, frequency_ratio, which_phase='Aa'):

    u = PyX_Utility.PyX_Utility()
    dpnv_grouping_dict = dict()
    dpnv_grouping_dict['GAC'] = []
    dpnv_grouping_dict['GBD'] = []

    # Draw 180e band (IT IS VERY IMPORTANT THAT WE USE -120 and -240 INSTEAD OF 120 and 240 HERE)
    if which_phase == 'Aa':
        band_phase_shift = -0
    elif which_phase == 'Bb':
        band_phase_shift = -120 # *0 == fix 180e band
    elif which_phase == 'Cc':
        band_phase_shift = -240 # *0 == fix 180e band
    else:
        raise Exception('Wrong which_phase: %s'%(which_phase))
    LB =     90-BELT_BIAS+band_phase_shift
    UB = 180+90-BELT_BIAS+band_phase_shift
    u.pyx_draw_sector( [0,0], LB, UB, RADIUS+1)

    # print(connection_star_raw_dict)
    # quit()

    # Draw the rest
    for key, val in connection_star_raw_dict.items():

        if key in which_phase:

            print('||||||', key, val)

            if key in 'abc':
                factor_reverse = -1
                phase_shift    = 180
            else: # in 'ABC'
                factor_reverse = 1
                phase_shift    = 0

            for el in val:
                # val: [(300.0, 11), (330.0, 12)]
                # el: (330.0, 12)

                # unpack
                el[0] # PHI at the original frequency
                PHI_ori = el[0] * frequency_ratio # PHI at the new frequency # at suspension frequency == 1/p*ps   or   at torque frequency == 1/ps*p
                PHI, label = phase_shift+PHI_ori, str(factor_reverse*el[1])


                # print('Group a/c:', label, LB, UB, PHI_ori, limit_to_360_deg(PHI), phase_shift)

                # Grouping a/c
                if belong_to_band(LB, UB, PHI):
                    print('Group a/c:', label, LB, UB, limit_to_360_deg(PHI))
                    dpnv_grouping_dict['GAC'] += [label]
                else:
                    print('Group b/d:', label, LB, UB, limit_to_360_deg(PHI))
                    dpnv_grouping_dict['GBD'] += [label]

                # Draw arrows
                u.pyx_arrow(angular_location(PHI))

                # The labels should drawn to the side a bit
                PHI += angle_between_arrow_and_label * factor_reverse

                # Draw labels
                if False:
                    # Option 1 (looks crappier)
                    u.pyx_text(angular_location(PHI+6*factor_reverse, radius_bias=-(PHI_ori//360)*1.0), label, scale=1)
                else:
                    # Option 2 (looks better)
                    # distance_between_label_layers = 1.25
                    X = (RADIUS) * np.cos(4*factor_reverse/180*np.pi) # (RADIUS, 0)旋转一个小角度（比如4度）
                    Y = (RADIUS) * np.sin(4*factor_reverse/180*np.pi) # (RADIUS, 0)旋转一个小角度（比如4度）
                    X += -(PHI_ori//360)*distance_between_label_layers # 水平移动
                    angle = PHI/180*np.pi
                    X, Y = X*np.cos(angle) + Y*-np.sin(angle), X*np.sin(angle) + Y*np.cos(angle) # Park Trans.
                    u.pyx_text((X,Y), r'\textbf{'+label+'}', scale=1)

                    # Draw 辅助线
                    u.pyx_circle(RADIUS-(PHI//360)*distance_between_label_layers, bool_dashed=True, dash_list=[0,24], linewidth=0.020)
    
    # backward compatable
    dpnv_grouping_dict[which_phase] = dpnv_grouping_dict['GAC']
    return u, dpnv_grouping_dict

def winding_distribution_factor_verPyrhonen(h, Q, p, m=3):

    # This is tested and shown to be wrong. cjh 2020-10-08
    # This is tested and shown to be wrong. cjh 2022-04-23 with h=1, Q=12, p=4, m=3

    # Option 2: Winding distribution factor Pyrhonen@(2.24)
    # alpha_u = 2*np.pi / Q * h*p # <- should be h*p or *p???
    alpha_u = 2*np.pi / Q * p # <- should be h*p or *p???
    # print('#debug', alpha_u/np.pi*180)
    q=Q/(2*p*m)
    if q%0 != 0:
        raise Exception('The formula is only valid for integral slot winding.')
    winding_distribution_factor_kw_at_h = np.sin(h*q*alpha_u/2) / (q*np.sin(h*alpha_u/2))
    return winding_distribution_factor_kw_at_h

def winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=None, Aa='Aa'):
    # Winding distribution factor Severson-2017-TIA

    # for key, val in connection_star_raw_dict.items():
    #     print(key, val)
    # quit()
    # A [(0.0, 1), (30.0, 2)]
    # Z [(60.0, 3), (90.0, 4)]
    # B [(120.0, 5), (150.0, 6)]
    # X [(180.0, 7), (210.0, 8)]
    # C [(240.0, 9), (270.0, 10)]
    # Y [(300.0, 11), (330.0, 12)]

    list_slots_of_a_phase = []
    list_phase_shift     = []
    for key, val in connection_star_raw_dict.items():
        if key in Aa:
            list_slots_of_a_phase += [el[-1] for el in val]
            if key in Aa[0]:
                list_phase_shift += [0]*len(val)
            elif key in Aa[1]:
                list_phase_shift += [np.pi]*len(val)

    # phasors in DPNV group a/c will be flipped (phase shifted 180^e)
    if phase_Aa_dpnv_grouping_dict is not None:
        phase_Aa_dpnv_grouping_list = phase_Aa_dpnv_grouping_dict[Aa]
        reversed_excitation_upper_layer = [abs(int(el)) for el in phase_Aa_dpnv_grouping_list]
        # print('reversed_excitation_upper_layer:', reversed_excitation_upper_layer)
        # print(val)
        for idx, slot_number in enumerate(list_slots_of_a_phase):
            if slot_number in reversed_excitation_upper_layer:
                # print(slot_number, reversed_excitation_upper_layer)
                list_phase_shift[idx] += np.pi

    # print(list_slots_of_a_phase)
    if bool_double_layer_winding == True:
        number_of_coils_per_phase = len(list_slots_of_a_phase)
        N = number_of_coils_per_phase # N = Q/m
    else:
        number_of_coil_sides_per_phase = len(list_slots_of_a_phase)
        N = number_of_coil_sides_per_phase # N = Q/m
    # print(N, list_phase_shift)
    winding_distribution_factor_kw_at_h = abs( 
                                               sum( 
                                                    [ np.exp( 
                                                              1j * (phase_angle_of_slot_i_at_frequency_h(slot_number, h, Q)+phase_shift)
                                                            ) for slot_number, phase_shift in zip(list_slots_of_a_phase, list_phase_shift) 
                                                    ]
                                                  )
                                             ) / N
    return winding_distribution_factor_kw_at_h
def draw_turn_function(connection_star_raw_dict, coil_pitch_y, Q, phase_Aa_dpnv_grouping_list=None, turn_func_bias=0, Aa='Aa', bool_double_layer_winding=True):

    def larger_than_Q(slot_number):
        if slot_number>Q:
            return slot_number - Q
        else:
            return slot_number

    if phase_Aa_dpnv_grouping_list is not None:
        reversed_excitation_upper_layer = [abs(int(el)) for el in phase_Aa_dpnv_grouping_list]
        reversed_excitation_lower_layer = [larger_than_Q(abs(int(el))+coil_pitch_y) for el in phase_Aa_dpnv_grouping_list]
        # print('|||||', reversed_excitation_upper_layer, reversed_excitation_lower_layer)
    else:
        reversed_excitation_upper_layer = reversed_excitation_lower_layer = None

    # initialize the canvas
    u = PyX_Utility.PyX_Utility()

    if bool_double_layer_winding == False:
        return u
        raise Exception('Not implemented for single layer.')

    phase_Aa_upper_layer_list = []
    phase_Aa_upper_layer_list_2 = []
    phase_Aa_lower_layer_list_2 = []
    connection_list = []
    for key, val in connection_star_raw_dict.items():
        if key not in Aa:
            continue
        print(key,val)
        phase_Aa_upper_layer_list   += [el[0] for el in val]
        phase_Aa_upper_layer_list_2 += [el[1]                               for el in val]
        phase_Aa_lower_layer_list_2 += [larger_than_Q(el[1] + coil_pitch_y) for el in val]
        if key in Aa[-1]:
            connection_list += [-1]*len(val)
        else:
            connection_list += [1]*len(val)

    # get accumulated_turns_min first to determine the location of the text "# 1, 2, 3, 4, ..., Q"
    turns_per_coil_side = 5
    accumulated_turns = 0 + turn_func_bias
    accumulated_turns_min = 0
    for i in range(Q):
        slot_number = i+1
        if slot_number in phase_Aa_upper_layer_list_2:
            sign = connection_list[phase_Aa_upper_layer_list_2.index(slot_number)]

            if reversed_excitation_upper_layer is not None and slot_number in reversed_excitation_upper_layer:
                sign *= -1
            accumulated_turns += sign*turns_per_coil_side
            if accumulated_turns < accumulated_turns_min:
                accumulated_turns_min = accumulated_turns

        if slot_number in phase_Aa_lower_layer_list_2:
            sign = -1 * connection_list[phase_Aa_lower_layer_list_2.index(slot_number)]

            if reversed_excitation_lower_layer is not None and slot_number in reversed_excitation_lower_layer:
                sign *= -1

            accumulated_turns += sign*turns_per_coil_side
            if accumulated_turns < accumulated_turns_min:
                accumulated_turns_min = accumulated_turns
    # print(accumulated_turns_min, '||||')
    # quit()
    slot_spacing = 5
    accumulated_turns = 0 + turn_func_bias
    last_location = [0, accumulated_turns]
    BIAS_LAYER = -RADIUS - 10
    # BIAS_LAYER = -RADIUS - 5 - 5
    示意用Marker位置 = - 10
    for i in range(Q):
        slot_number = i+1
        this_location = [slot_number*slot_spacing, accumulated_turns]

        # 1, 2, 3, 4, ..., Q
        u.pyx_text([this_location[0], 示意用Marker位置+1*accumulated_turns_min], str(i+1), scale=2)

        # Draw lines
        u.pyx_line(last_location, this_location)
        last_location = this_location

        # Draw marker from upper layer
        if slot_number in phase_Aa_upper_layer_list_2:
            sign = connection_list[phase_Aa_upper_layer_list_2.index(slot_number)]

            if reversed_excitation_upper_layer is not None and slot_number in reversed_excitation_upper_layer:
                sign *= -1

            accumulated_turns += sign*turns_per_coil_side
            # if accumulated_turns < accumulated_turns_min:
            #     accumulated_turns_min = accumulated_turns

            # Draw markers
            if sign == 1:
                u.pyx_marker_plus(this_location, size=0.8) # 波形中的Marker
                u.pyx_marker_plus([this_location[0], 示意用Marker位置+5 -0 + 0+1*accumulated_turns_min], size=0.8) # 示意用Marker
            else:
                u.pyx_marker_minus(this_location, size=0.8) # 波形中的Marker
                u.pyx_marker_minus([this_location[0], 示意用Marker位置+5 -0 + 0+1*accumulated_turns_min], size=0.8) # 示意用Marker

            this_location = [slot_number*slot_spacing, accumulated_turns]

        # Draw lines
        u.pyx_line(last_location, this_location)
        last_location = this_location

        # Draw marker from lower layer
        if slot_number in phase_Aa_lower_layer_list_2:
            sign = -1 * connection_list[phase_Aa_lower_layer_list_2.index(slot_number)]

            if reversed_excitation_lower_layer is not None and slot_number in reversed_excitation_lower_layer:
                sign *= -1

            accumulated_turns += sign*turns_per_coil_side
            # if accumulated_turns < accumulated_turns_min:
            #     accumulated_turns_min = accumulated_turns

            # Draw markers
            if sign == 1:
                u.pyx_marker_plus(this_location, size=0.8, rgb=[0,0,1]) # 波形中的Marker
                u.pyx_marker_plus([this_location[0], 示意用Marker位置+5 -0 - 2.5+1*accumulated_turns_min], size=0.8, rgb=[0,0,1]) # 示意用Marker
            else:
                u.pyx_marker_minus(this_location, size=0.8, rgb=[0,0,1]) # 波形中的Marker
                u.pyx_marker_minus([this_location[0], 示意用Marker位置+5 -0 - 2.5+1*accumulated_turns_min], size=0.8, rgb=[0,0,1]) # 示意用Marker

            this_location = [slot_number*slot_spacing, accumulated_turns]

        # Draw lines
        u.pyx_line(last_location, this_location)
        last_location = this_location

    # Draw upper and lower layer markers for ease of inspection
    u.pyx_line(                 [-0, 示意用Marker位置 +6.5+1*accumulated_turns_min], 
                [(Q+0)*slot_spacing, 示意用Marker位置 +6.5+1*accumulated_turns_min], arg=['dashed']) # , 'my-thick-line'
                                                                # pyx.style.dash(dash_list), pyx.style.linewidth(linewidth)
    # u.pyx_text([-1, -RADIUS-0 +3], 'Coil pitch is %d'%(coil_pitch_y), scale=1)

    # Reference Aaes
    u.pyx_arrow([-25-5+20, -1], [-25-3+20,  -1])
    u.pyx_arrow([-25-5+20, -1], [-25-5+20,  +1])
    u.pyx_text([-25+2+20, -2.5], 'Slot number', scale=2)
    u.pyx_text([-25-1+20,    2], 'Ampere-turns', scale=2)

    # print(phase_Aa_upper_layer_list_2)
    # print(phase_Aa_lower_layer_list_2)
    # print(connection_list)
    return u

def winding_short_pitch_factor(h, coil_pitch_y, Q, npp):
    y_Q = Q / (2*npp)

    # h: harmonic order
    # k_ph = np.sin(h*npp * coil_pitch_y/y_Q * np.pi*0.5) # if you use this, h*npp = v*npp = n, 此时h=3，是指极对数为3个npp的谐波。
                                                          # 此时，h的意义是相对于npp这个磁场为三次的谐波磁场，如果npp为2，那么h=3所对应的是气隙中6对极的谐波。

    k_ph = np.sin(h/npp * coil_pitch_y/y_Q * np.pi*0.5)   # if you use this, h/npp = n,         此时h=3，是指极对数为3的谐波。   <- 这是我们要的，
                                                          # 我们要计算气隙中某个极对数的谐波所对应分布绕组和短距绕组，并相乘，所以净极对数(h)要对得上。
                                                          # 换句话说，我们要的是h=3次谐波的短距系数，而不是相对于ps为3次的谐波的短距系数。
    # coil_pitch_y / (Q/2) * np.pi is the short pitch radian for 1 pole pair field
    # coil_pitch_y / (Q/4) * np.pi is the short pitch radian for 2 pole pair field
    # coil_pitch_y / (Q/(2*npp)) * np.pi is the short pitch radian for npp pole pair field
    # Bb "k_ph = np.sin(h/npp * coil_pitch_y/y_Q * np.pi*0.5)", you are trying to use the short pitch radian for npp pole pair field to calculate the pitch factor for a h pole pair field.
    # That is why you first need to convert the short pitch radian to mechanical radian first (divided Bb npp) and then convert it to h pole pair field's (multiplied Bb h).

    # Slack Conversation
        # Jiahao:
            # Hello Eric,
            # I am checking the windings.  Now the only thing seems not right is the pitch factor for the suspension winding.
            # P.S.: I am using (2.32) from Pyrhonen's book and it is correct for the torque winding.
            # For example, for the Q12, p1, ps2, y5 design, we have:
            # the distribution factor for harmonic h=2 is cos(30 deg)=0.866;
            # the pitch factor for harmonic h=2 is sin(2 * 5/ (12/4) * 90 deg)=-0.866;
            # so, the full winding factor should be kw = 0.866*0.866 = 0.75;
            # But your Example pdf file gives kw = 0.433.
            # Hypothetically, if one makes a mistake when calculating y_Q for suspension winding, the pitch factor for harmonic h=2 could become: h=2 is sin(2 * 5/ (12/2) * 90 deg)=0.5, such that kw = 0.866 * 0.5 = 0.433, which matches your results.
        # Eric:
            # Quick response: pitch factor is kpn = sin(n \gamma/2), where gamma is the coil pitch in (elec.) radians and n=vp.

    # 短距长距，角度差个符号，cos一下没了。
        # if coil_pitch_y > y_Q:
        #     # 长距
        #     # print('[Warning] Over-pitch winding detected!')
        #     the_angle = h* 0.5* (-np.pi + coil_pitch_y/y_Q*np.pi)
        #     k_ph = np.cos(the_angle)
        # elif abs(coil_pitch_y - y_Q)<1e-3: # EPS
        #     k_ph = 1.0
        # else:
        #     # 短距
        #     the_angle = h* 0.5* (np.pi - coil_pitch_y/y_Q*np.pi)
        #     k_ph = np.cos(the_angle)
        #     # k_ph = np.sin(h * coil_pitch_y/y_Q * np.pi*0.5) # this is only valid for short pitching
    return k_ph
def winding_short_pitch_factor_v2(h, coil_pitch_y, Q):
    y_Q = Q / (2) 
    k_ph = np.sin(h * coil_pitch_y/y_Q * np.pi*0.5)
    # Here, coil_pitch_y/y_Q * np.pi is the short pitch radian (elec.) for 1 pole pair field
    return k_ph

import pyx

class Winding_Derivation(object):
    """General Implementation of the Winding_Derivation Process."""
    def __init__(self, slot_pole_comb, bool_double_layer_winding=True):

        # General info
        self.m = m = slot_pole_comb[0]
        self.Q = Q = slot_pole_comb[1]
        self.p = p = slot_pole_comb[2]
        self.ps = ps = slot_pole_comb[3]
        self.coil_pitch_y = coil_pitch_y = slot_pole_comb[4]
        self.turn_func_bias = turn_func_bias = slot_pole_comb[5]

        self.t = t = gcd(Q,p) # For the fundamental frequency, there will be t phasors at each angular location and therefore t levels to the phasor star.
        self.ts= ts = gcd(Q,ps)
        self.q = q = Q/(2*p)/m
        self.qs= qs = Q/(2*ps)/m
        print(f'Q={Q} \np={p} \nps={ps} \nt={t} \nts={ts} \nq={q} \nqs={qs}.')
        if p%2==0: 
            print('DPNV Type 1.'); 
        elif p%2==1 and ps%2==0: 
            print('DPNV Type 2.') 
        else: 
            print('DPNV Type 3.')
        if q<=1:
            print('This winding has no distribution factor, right?')



        # Torque frequency as fundamental (h=p)
        print(f'Torque star of slots with Q={Q} and p={p}.')
        drawer_T1, connection_star_raw_dict, phase_belt = draw_star_of_slots(Q, p, m)
        self.connection_star_raw_dict = connection_star_raw_dict
        drawer_T2                           = draw_connection_star(m, phase_belt, connection_star_raw_dict)
        if m == 3: # suspension winding supports only 3 phase for now (it can be extended for multiphase)
            drawer_T3a, dpnv_grouping_dict_a    = draw_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*ps, which_phase='Aa') # p is original frequency, ps is the new frequency
            drawer_T3b, dpnv_grouping_dict_b    = draw_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*ps, which_phase='Bb') # (Unchanged) NOTE THAT THE PHASE V and W are transposed!
            drawer_T3c, dpnv_grouping_dict_c    = draw_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*ps, which_phase='Cc') # (Unchanged) NOTE THAT THE PHASE V and W are transposed!
            self.dpnv_grouping_dict_a = dpnv_grouping_dict_a
            self.dpnv_grouping_dict_b = dpnv_grouping_dict_b
            self.dpnv_grouping_dict_c = dpnv_grouping_dict_c
        # drawer_T33, _                       = draw_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*(3*p), which_phase='Aa') # 3次谐波
        # drawer_T35, _                       = draw_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*(5*p), which_phase='Aa') # 5次谐波
        # drawer_T37, _                       = draw_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*(7*p), which_phase='Aa') # 7次谐波
        # drawer_T1.pyx_text([0,RADIUS+5], '1) Torque phasor star ($Q=%d$, $p=%d$, $t=%d$)'%(Q,p,t))
        # drawer_T1.pyx_text([0,RADIUS+3], r'$q = %g$, $Q^\prime=%g$, $p^\prime=%g$' % (q, Q/t, p/t)); self.Q_prime = Q/t
        # drawer_T2.pyx_text([0,RADIUS+5], '2) Torque connection star (Pursue mAaimum emf)', scale=1)
        # drawer_T2.pyx_text([0,RADIUS+3], '2) Torque connection star', scale=1)

        # for ISMB 2021
        drawer_T1.pyx_text([0,RADIUS+5], 'Torque phasor star ($Q_s=%d$, $p=%d$, $t=%d$)'%(Q,p,t))
        drawer_T1.pyx_text([0,RADIUS+3], r'$q = %g$, $Q^\prime_s=%g$, $p^\prime=%g$' % (q, Q/t, p/t)); self.Q_prime = Q/t
        drawer_T2.pyx_text([0,RADIUS+5], 'Torque connection star (Pursue maximum emf)', scale=1)
        drawer_T2.pyx_text([0,RADIUS+3], 'Torque connection star', scale=1)



        def get_list_phase_slot_number(which_phase_belt_A='A', which_phase_belt_a='a'):
            list_phase_u_slot_number = [[el[-1] for el in val] for key, val in connection_star_raw_dict.items() if key in which_phase_belt_A]
            list_phase_u_slot_number += [[-el[-1] for el in val] for key, val in connection_star_raw_dict.items() if key in which_phase_belt_a]
            list_phase_u_slot_number = [item for sublist in list_phase_u_slot_number for item in sublist]
            return list_phase_u_slot_number
        if m == 3:
            self.list_phase_u_slot_number = list_phase_u_slot_number = get_list_phase_slot_number('A', 'a')
            drawer_T2.pyx_text([0,-RADIUS-2], 'Phase U: ' + ', '.join([str(el) for el in list_phase_u_slot_number]), scale=1)
            self.list_phase_v_slot_number = list_phase_v_slot_number = get_list_phase_slot_number('B', 'b')
            drawer_T2.pyx_text([0,-RADIUS-4], 'Phase V: ' + ', '.join([str(el) for el in list_phase_v_slot_number]), scale=1)
            self.list_phase_w_slot_number = list_phase_w_slot_number = get_list_phase_slot_number('C', 'c')
            drawer_T2.pyx_text([0,-RADIUS-6], 'Phase W: ' + ', '.join([str(el) for el in list_phase_w_slot_number]), scale=1)
        else:
            # solution for m>3
            self.list_slot_number_of_phase = dict()
            for _phaseNumber in range(m):
                _phaseName = ABCDEFGHIJKLMNOPQRSTUVWXYZ[_phaseNumber]
                self.list_slot_number_of_phase[_phaseName] = get_list_phase_slot_number(_phaseName, _phaseName.lower())
                drawer_T2.pyx_text([0,-RADIUS-3-3/(m/3)*_phaseNumber], f'Phase {_phaseName}: ' + ', '.join([str(el) for el in self.list_slot_number_of_phase[_phaseName]]), scale=3/(m/3))

        if m == 3:
            # for ISMB 2021
            # drawer_T3a.pyx_text([0,RADIUS+3], '3) Torque conn. star at sus. freq. U', scale=1); # drawer_T3a.pyx_text([0,RADIUS+2.5], 'Should have null emf but if you flip phasors in the $180^e$ band you have sus. emf')
            # drawer_T3b.pyx_text([0,RADIUS+3], '3) Torque conn. star at sus. freq. W', scale=1)
            # drawer_T3c.pyx_text([0,RADIUS+3], '3) Torque conn. star at sus. freq. V', scale=1)
            # drawer_T33.pyx_text([0,RADIUS+5], '3) Torque conn. star at 3rd harmonic (ph.U)', scale=1)
            # drawer_T35.pyx_text([0,RADIUS+5], '3) Torque conn. star at 5th harmonic (ph.U)', scale=1)
            # drawer_T37.pyx_text([0,RADIUS+5], '3) Torque conn. star at 7th harmonic (ph.U)', scale=1)
            drawer_T3a.pyx_text([0,-RADIUS-3], '*Group a/c: ' + ', '.join(dpnv_grouping_dict_a['Aa']) + ' (phase U)', scale=1)
            drawer_T3b.pyx_text([0,-RADIUS-3], '*Group a/c: ' + ', '.join(dpnv_grouping_dict_b['Bb']) + ' (phase W)', scale=1) # NOTE THAT THE PHASE V and W are transposed here!
            drawer_T3c.pyx_text([0,-RADIUS-3], '*Group a/c: ' + ', '.join(dpnv_grouping_dict_c['Cc']) + ' (phase V)', scale=1) # NOTE THAT THE PHASE V and W are transposed here!


        if bool_double_layer_winding:
            # Torque
            drawer_T4 = draw_turn_function(connection_star_raw_dict, coil_pitch_y, Q, turn_func_bias=turn_func_bias)
            # for ISMB 2021
            # drawer_T4.pyx_text([0, RADIUS+5], '4) Torque turn function', scale=1)

            if m == 3:
                # Suspension
                # print('|||', dpnv_grouping_dict_a)
                drawer_T4a = draw_turn_function(connection_star_raw_dict, coil_pitch_y, Q, phase_Aa_dpnv_grouping_list=dpnv_grouping_dict_a['Aa'], Aa='Aa', turn_func_bias=turn_func_bias)
                drawer_T4b = draw_turn_function(connection_star_raw_dict, coil_pitch_y, Q, phase_Aa_dpnv_grouping_list=dpnv_grouping_dict_b['Bb'], Aa='Bb', turn_func_bias=turn_func_bias)
                drawer_T4c = draw_turn_function(connection_star_raw_dict, coil_pitch_y, Q, phase_Aa_dpnv_grouping_list=dpnv_grouping_dict_c['Cc'], Aa='Cc', turn_func_bias=turn_func_bias)

                # for ISMB 2021
                # drawer_T4a.pyx_text([0, RADIUS+5], '4) Sus. turn function U', scale=1)
                # drawer_T4b.pyx_text([0, RADIUS+5], '4) Sus. turn function V', scale=1)
                # drawer_T4c.pyx_text([0, RADIUS+5], '4) Sus. turn function W', scale=1)
        else:
            drawer_T4 = draw_turn_function(connection_star_raw_dict, coil_pitch_y, Q, turn_func_bias=turn_func_bias)
            msg = 'Turn function for single layer winding is not implemented yet.'
            drawer_T4.pyx_text([0, RADIUS+5], msg, scale=1)
            print(msg)







        drawer_Text = PyX_Utility.PyX_Utility()

        kw_at_p  = winding_distribution_factor(Q, connection_star_raw_dict, p, bool_double_layer_winding)
        kw_at_ps = winding_distribution_factor(Q, connection_star_raw_dict, ps, bool_double_layer_winding)
        print('- Winding dristribution factor at p (=%d): %g.' % (p, kw_at_p))
        print('- Winding dristribution factor at ps (=%d): %g.' % (ps, kw_at_ps))
        if m==3:
            kw_at_ps  = winding_distribution_factor(Q, connection_star_raw_dict, ps, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_a, Aa='Aa')
            print('- Winding dristribution factor at ps (=%d) with flipped phasors in group a/c: %g.' % (ps, kw_at_ps))

        print('- Winding dristribution factor at h:');
        drawer_Text.pyx_text([-PLOT_SPACING*3, -5 -RADIUS-5], 'Torque winding dist. factor at h freq.:')
        drawer_Text.pyx_text([-PLOT_SPACING*2, -5 -RADIUS-5], 'Sus. winding dist. factor at h freq.:')
        count_TW = count_SW = 0
        self.torque_kd_at_h = dict()
        self.suspen_kd_at_h = dict()
        # import utility
        for h in range(25):
            temp_TW   = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding)
            # print(h, ': %g'%(temp_TW), end=' \n ')
            if m==3:
                temp_SW_A = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_a, Aa='Aa')
                temp_SW_B = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_b, Aa='Bb')
                temp_SW_C = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_c, Aa='Cc')
            # temp_TW   = winding_distribution_factor_verPyrhonen(h, Q, p)
            # print(temp_TW)
            # utility.blockPrint()
            # temp_SW_A = winding_distribution_factor_verPyrhonen(h, Q, ps)
            # temp_SW_B = winding_distribution_factor_verPyrhonen(h, Q, ps)
            # temp_SW_C = winding_distribution_factor_verPyrhonen(h, Q, ps)
            if abs(temp_TW) > 1e-4:
                msg = '\t%d: %.3f'%(h, temp_TW); print('\tT:', msg); count_TW += 1; self.torque_kd_at_h[h] = temp_TW
                drawer_Text.pyx_text([-PLOT_SPACING*3, -5 -RADIUS-5 -2.5*count_TW], msg)
            if m==3:
                if abs(temp_SW_A) > 1e-4:
                    msg = '\t%d: %.2f,\t%.2f,\t%.2f'%(h, temp_SW_A, temp_SW_B, temp_SW_C); print('\tS:', msg); count_SW += 1; self.suspen_kd_at_h[h] = temp_SW_A, temp_SW_B, temp_SW_C
                    drawer_Text.pyx_text([ -PLOT_SPACING*2, -5 -RADIUS-5 -2.5*count_SW], msg)
            # utility.enablePrint()

            # #debug
            # if h>1:
            #     quit()
        # quit()

        print('- Winding short pitch factor at h:') 
        drawer_Text.pyx_text([-PLOT_SPACING, -5 -RADIUS-5], 'Torque winding pitch factor at h freq.:')
        drawer_Text.pyx_text([            0, -5 -RADIUS-5], 'Sus. winding pitch factor at h freq.:')
        count_TW = count_SW = 0
        self.torque_kp_at_h = dict()
        self.suspen_kp_at_h = dict()
        for h in range(25):
            temp_TW = winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            temp_SW = winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            # print('\t'*5, h, temp_TW, temp_SW)
            if abs(temp_TW) > 1e-4:
                msg = '\t%d: %.3f'%(h, temp_TW); print('\t\tT:', msg); count_TW += 1; self.torque_kp_at_h[h] = temp_TW
                drawer_Text.pyx_text([-PLOT_SPACING, -5 -RADIUS-5 -2.5*count_TW], msg)
            if abs(temp_SW) > 1e-4:
                msg = '\t%d: %.3f'%(h, temp_SW); print('\t\tS:', msg); count_SW += 1; self.suspen_kp_at_h[h] = temp_SW
                drawer_Text.pyx_text([0, -5 -RADIUS-5 -2.5*count_SW], msg)

        print('- Winding factor at h:')
        drawer_Text.pyx_text([  PLOT_SPACING, -5 -RADIUS-5], 'Torque winding factor at h freq.:')
        drawer_Text.pyx_text([PLOT_SPACING*2, -5 -RADIUS-5], 'Sus. winding factor at h freq.:')
        count_TW = count_SW = 0
        self.torque_kw_at_h = dict()
        self.suspen_kw_at_h = dict()
        for h in range(25):
            temp_TW = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding)
            temp_TW *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            # print(h, temp_TW)
            if m==3:
                temp_SW_A = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_a, Aa='Aa')
                temp_SW_B = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_b, Aa='Bb')
                temp_SW_C = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_c, Aa='Cc')
                temp_SW_A *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
                temp_SW_B *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
                temp_SW_C *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            if abs(temp_TW) > 1e-4:
                msg = '\t%d: %.3f'%(h, temp_TW); print(msg); count_TW += 1; self.torque_kw_at_h[h] = temp_TW
                drawer_Text.pyx_text([PLOT_SPACING, -5 -RADIUS-5 -2.5*count_TW], msg)
            if m==3:
                if abs(temp_SW_A) > 1e-4:
                    msg = '\t%d: %.2f,\t%.2f,\t%.2f'%(h, temp_SW_A, temp_SW_B, temp_SW_C); print('\t\t\tS:', msg); count_SW += 1; self.suspen_kw_at_h[h] = temp_SW_A, temp_SW_B, temp_SW_C
                    drawer_Text.pyx_text([PLOT_SPACING*2, -5 -RADIUS-5 -2.5*count_SW], msg)






















        self.drawer_T1 = drawer_T1
        self.drawer_T2 = drawer_T2
        self.drawer_T4 = drawer_T4
        if m == 3:
            self.drawer_T3a = drawer_T3a
            self.drawer_T3b = drawer_T3b
            self.drawer_T3c = drawer_T3c
            self.drawer_T4a = drawer_T4a
            self.drawer_T4b = drawer_T4b
            self.drawer_T4c = drawer_T4c
        self.drawer_Text = drawer_Text

    def get_complex_number_winding_factor_of_coil_i(self, i, coil_pitch_y, Q, v, p):
        ''' In this formulation, the basic component is a coil rather than a coill side.
            absolute harminic index h = v*p, with v the relative harmonic index w.r.t. p.
        '''

        alpha_u = 2*np.pi / Q *p

        # print(v, p, alpha_u, coil_pitch_y)
        gamma =            coil_pitch_y*alpha_u/p
        radii = np.sin(v*p*coil_pitch_y*alpha_u/p/2)
        print(f'Coil span [mech.deg] = {gamma/np.pi*180} | [elec.deg] = {v*p*gamma / np.pi*180}', end=' | ')
        # print('\t radii:', radii)
        ELS_angles = -0.5*np.pi - v*p*alpha_u*(2*i+coil_pitch_y) / (2*p) # IMPORTANT: this is in elec.rad!!!
        CJH_angles =  0.5*np.pi - v*p*alpha_u*(2*i+coil_pitch_y) / (2*p) # IMPORTANT: this is in elec.rad!!!
        ELS_pitch_factor_per_coil = radii * np.exp(1j*ELS_angles)
        CJH_pitch_factor_per_coil = radii * np.exp(1j*CJH_angles)
        return ELS_pitch_factor_per_coil, CJH_pitch_factor_per_coil

    def get_complex_number_kw_per_phase(self, v, p, positive_connected_coils, negative_connected_coils):

        kp_els_list = []
        kp_cjh_list = []

        for i in positive_connected_coils:
            els, cjh = self.get_complex_number_winding_factor_of_coil_i(i, self.coil_pitch_y, Q=self.Q, v=v, p=p)
            # this is pitch factor of coil in positive zone
            kp_els_list.append(els)
            kp_cjh_list.append(cjh)

            print(f'\tkp@Coil+{i:02d}', end='\t|\t')
            # print('\tels = %g∠%.1f' % (np.abs(els), np.angle(els)/np.pi*180*1), end='\t|\t')
            print(  'cjh = %.3f∠%.1f' % (np.abs(cjh), np.angle(cjh)/np.pi*180*1))

        for i in negative_connected_coils:
            els, cjh = self.get_complex_number_winding_factor_of_coil_i(i, self.coil_pitch_y, Q=self.Q, v=v, p=p)
            # this is pitch factor of coil in negative zone
            els *= -1 # np.abs(els) * np.exp(1j*(np.angle(els)*1+np.pi)/p)
            cjh *= -1 # np.abs(cjh) * np.exp(1j*(np.angle(cjh)*1+np.pi)/p)
            kp_els_list.append(els)
            kp_cjh_list.append(cjh)

            print(f'\tkp@Coil-{i:02d}', end='\t|\t')
            # print('\tels = %g∠%.1f' % (np.abs(els), np.angle(els)/np.pi*180*1), end='\t|\t')
            print(  'cjh = %.3f∠%.1f' % (np.abs(cjh), np.angle(cjh)/np.pi*180*1))

        ''' When you sum up the vectors, they must be first converted to using elec.rad by multiplying the angle by p pole pairs.
        '''
        average = lambda x: np.sum(x)/len(x)
        # print([np.angle(el)/np.pi*180 for el in kp_els_list])
        # quit()
        _kw_els = average(kp_els_list)
        _kw_cjh = average(kp_cjh_list)
        return _kw_els, _kw_cjh

    def get_complex_number_kw(self, p_or_ps, v=1, bool_study_suspension_subharmonics=False):

        dict_kw_els = dict()
        dict_kw_cjh = dict()

        for ZONE in ['A', 'B', 'C']:
            pcc = [i for _angle, i in self.connection_star_raw_dict[ZONE]]
            try:
                ncc = [i for _angle, i in self.connection_star_raw_dict[ZONE.lower()]]
            except KeyError:
                print('负60度相带里没有槽电势矢量，所以abc不存在于connection_star_raw_dict.keys()。')
                ncc = [] # no negative connected coils

            if p_or_ps == self.ps or bool_study_suspension_subharmonics:
                # suspension winding has different connection patten from the torque winding
                if ZONE == 'A':
                    dpnv_grouping_AC = self.dpnv_grouping_dict_a['GAC']
                    dpnv_grouping_BD = self.dpnv_grouping_dict_a['GBD']
                if ZONE == 'B':
                    dpnv_grouping_AC = self.dpnv_grouping_dict_b['GAC']
                    dpnv_grouping_BD = self.dpnv_grouping_dict_b['GBD']
                if ZONE == 'C':
                    dpnv_grouping_AC = self.dpnv_grouping_dict_c['GAC']
                    dpnv_grouping_BD = self.dpnv_grouping_dict_c['GBD']
                connection_star_raw_results = [-int(el) for el in dpnv_grouping_AC] + [int(el) for el in dpnv_grouping_BD]
                pcc = [    el  for el in connection_star_raw_results if el>0]
                ncc = [abs(el) for el in connection_star_raw_results if el<0]
                print('sus-pcc:', pcc)
                print('sus-ncc:', ncc)
                # quit()

            _kw_els, _kw_cjh = self.get_complex_number_kw_per_phase(v=v, p=p_or_ps, positive_connected_coils=pcc, negative_connected_coils=ncc)
            dict_kw_els[f'{ZONE}'] = _kw_els
            dict_kw_cjh[f'{ZONE}'] = _kw_cjh
            dict_kw_els[f'{ZONE}_abs'] = np.abs(_kw_els)
            dict_kw_cjh[f'{ZONE}_abs'] = np.abs(_kw_cjh)
            dict_kw_els[f'{ZONE}_angle'] = np.angle(_kw_els)/np.pi*180 # no need to convert mech.deg to elec.deg, as it is already in elec.deg
            dict_kw_cjh[f'{ZONE}_angle'] = np.angle(_kw_cjh)/np.pi*180 # no need to convert mech.deg to elec.deg, as it is already in elec.deg
            # print('kw(els)=', _kw_els, '=', f'{np.abs(_kw_els)}∠{np.angle(_kw_els)/np.pi*180}')
            print('kw(cjh)=', _kw_cjh, '=', f'{np.abs(_kw_cjh):.3f}∠{np.angle(_kw_cjh)/np.pi*180:.1f}')

        self.dict_kw_els = dict_kw_els
        self.dict_kw_cjh = dict_kw_cjh
        return dict_kw_els, dict_kw_cjh

def main_derivation():

    # m, Q, p, ps, y, turn function bias (turn_func_bias)
    Slot_Pole_Combinations = [
                                # (15, 30, 2, 3, 10, 0),
                                # (3, 24, 8, 7, 1, 0), # 
                                # (3, 24, 8, 9, 1, 0), # 
                                (3, 24, 8, 11, 1, 0), # 
                                # (3, 6, 5, 4, 1, 0), # concentrated winding with working harmonics
                                # (3, 6, 5, 4, 2, 0), # concentrated winding with working harmonics
                                # (3, 6, 10, 11, 3, 0), # concentrated winding with working harmonics
                                # (3, 6, 10, 11, 1, 0), # concentrated winding with working harmonics
                                # (3, 81, 3, 2, 13, 0), # 哔哩哔哩 NJ伯乐在此
                                # (3, 6, 5, 4, 2, 0), # 6s/14pp slot FSPM y=2 （y=3 has no suspension field)
                                # (3, 6, 1, 2, 2, 0), # 6s/14pp slot FSPM y=2
                                # (3, 6, 4, 5, 2, 0), # 6s/14pp slot FSPM y=2
                                # (3, 12, 14, 13, 1, 0), # 12 slot FSPM y=1
                                # (3, 12, 20, 19, 2, 0), # 12 slot FSPM y=2
                                # (3, 12, 22, 21, 2, 0), # 12 slot FSPM y=2
                                # (3, 12, 22, 23, 2, 0), # 12 slot FSPM y=2
                                # (3, 24, 14, 13, 1, 0), # 24 slot FSPM
                                # (3, 24, 10, 11, 1, 0), # 24 slot FSPM
                                # (3, 24, 10, 11, 2, 0), # 24 slot FSPM (y=2  is gabbage)
                                # (3, 24, 10, 9, 1, 0), # 24 slot FSPM (ps=9)
                                # (3, 24, 9, 1, 1, 0), # Type 3 winding with odd p and ps=1
                                # (3, 12, 4, 1, 1, 0), # Eric's Homopolar Motor??? p is even number?
                                # (3, 18, 4, 1, 1, 0), # Eric's Homopolar Motor??? p is even number?
                                # (3, 12, 7, 8, 1, 0), # BLDC motor concentrated winding
                                # (3, 12, 4, 5, 1, 0), # FSPM-12/10
                                # (3, 12, 4, 5, 1, 0), # PEMD p4ps5
                                # (3, 24, 1, 2, 9, 0),
                                # (3, 18, 2, 1, 1, 0), # David Meeker SPMLoss example  https://femm.info/wiki/spmloss
                                # (3, 48, 2, 1, 12, 0), # XSQ
                                # (3, 9, 4, 5, 1, 0), # Slessinv TIE.R1
                                # (3, 18, 4, 5, 2, 0), # Slessinv TIE.R1 750W Servo
                                # (3, 27, 3, 2, 4, 0), # 0 <--- 这个如果作绕组你会发现W相的Group AC和Group BD的阴影刚好差了一点角度，导致三相AC/BD分组不对称。
                                # (3, 36, 3, 2, 5, 0), # 1              # ISMB 2021 Winding (Only Amplitude asymmetry)
                                # (3, 36, 3, 2, 6, 0), # 2
                                # (3, 36, 3, 2, 4, 0), # 3
                                # (3, 24, 2, 3, 5, 0), # 4            # ISMB 2021 Winding (Phase asymmetry in working harmonics)
                                # (3, 24, 2, 3, 6, 0), # 5
                                # (3, 24, 2, 3, 4, 0), # 6
                                # (3, 18, 2, 3, 4, 0), # 7            # ISMB 2021 Winding (No asymmetry in working harmonics)
                                # (3, 18, 2, 3, 5, 0), # 8
                                # (3, 18, 2, 3, 3, 0), # 9
                                #   (3, 36, 2, 3, 8, 0), # 10
                                #   (3, 36, 2, 3, 9, 0), # 11
                                #   (3, 36, 2, 3, 7, 0), # 12
                                #   (3, 27, 2, 3, 6, 0), # 13
                                #   (3, 36, 2, 3,    8,   0), # Dec. 15, 2020
                                #   (3, 36, 4, 5,    3,   0), # Jan. 19, 2021
                                #   (3, 18, 8, 7,    1,   0), # Mar. 25, 2021 哔哩哔哩：一介介一
                                #   (3, 18, 2, 3, 4, 0),
                                #   (3, 36, 3, 4, 5, 0),              # ISMB 2021 Winding design 7
                                #   (3, 18, 3, 4, 2, 0), # phase asymmetry
                                #   (3, 18, 3, 4, 3, 0), # phase asymmetry
                                #   (3, 24, 4, 5, 3, 0),
                                #   (3, 12, 5, 6, 1, 0), # SH-WuJi
                             ]
    bool_double_layer_winding = True

    for index, slot_pole_comb in enumerate(Slot_Pole_Combinations):

        # if index != 4 and index != 3 and index != 2:
        #     continue
        # if index != 0:
        #     continue

        wd = Winding_Derivation(slot_pole_comb, bool_double_layer_winding)

        ''' NEW ISMB 2021 Complex Number Winding Factor '''
        print('\n----------------------------------n=p')
        dict_kw_els, dict_kw_cjh = wd.get_complex_number_kw(wd.p)
        print('\n----------------------------------n=ps')
        dict_kw_els, dict_kw_cjh = wd.get_complex_number_kw(wd.ps, v=1)
        # dict_kw_els, dict_kw_cjh = wd.get_complex_number_kw(wd.ps, v=1/3)
        # dict_kw_els, dict_kw_cjh = wd.get_complex_number_kw(1, v=1, bool_study_suspension_subharmonics=True)
        phase_difference = [
                            dict_kw_els['A_angle'] - dict_kw_els['B_angle'],
                            dict_kw_els['B_angle'] - dict_kw_els['C_angle'],
                            dict_kw_els['C_angle'] - dict_kw_els['A_angle'],
                           ]
        print('phases of winding:', dict_kw_els['A_angle'], dict_kw_els['B_angle'], dict_kw_els['C_angle'])
        print('phase_difference:', phase_difference)
        for el in phase_difference:
            if abs(abs(el) - 240.0)>1e-5 and abs(abs(el) - 120.0)>1e-5:
                print('!!!Phase Asymmetry Detected')

        # continue








        fname = output_dir + 'wily_p%dps%dQ%dy%d'%(wd.p, wd.ps, wd.Q, wd.coil_pitch_y)

        if False:
            ''' ISMB 2021: Produce sub-figure for the paper
            '''
            wd.drawer_T1.cvs.writePDFfile(fname + '_T1')
            wd.drawer_T2.cvs.writePDFfile(fname + '_T2')
            wd.drawer_T3a.cvs.insert(wd.drawer_T3b.cvs, [pyx.trafo.translate(20*2,  0)]) # NOTE THAT THE PHASE V and W are transposed!
            wd.drawer_T3a.cvs.insert(wd.drawer_T3c.cvs, [pyx.trafo.translate(20*1,  0)]) # NOTE THAT THE PHASE V and W are transposed!
            wd.drawer_T3a.cvs.writePDFfile(fname + '_T3abc')
            wd.drawer_T4.cvs.writePDFfile(fname + '_T4')
            wd.drawer_T4a.cvs.writePDFfile(fname + '_T4a')
            wd.drawer_T4b.cvs.writePDFfile(fname + '_T4b')
            wd.drawer_T4c.cvs.writePDFfile(fname + '_T4c')
            print(f'Write pdf to {fname}')
            quit()

        # Collage
        wd.drawer_T1.cvs.insert(wd.drawer_T2.cvs,  [pyx.trafo.translate(PLOT_SPACING*1,  0)])
        if wd.m==3:
            wd.drawer_T1.cvs.insert(wd.drawer_T3a.cvs, [pyx.trafo.translate(PLOT_SPACING*2,  0)])
            wd.drawer_T1.cvs.insert(wd.drawer_T3b.cvs, [pyx.trafo.translate(PLOT_SPACING*3,  0)])
            wd.drawer_T1.cvs.insert(wd.drawer_T3c.cvs, [pyx.trafo.translate(PLOT_SPACING*4,  0)])
        # wd.drawer_T1.cvs.insert(wd.drawer_T33.cvs, [pyx.trafo.translate(PLOT_SPACING*5,  0)])
        # wd.drawer_T1.cvs.insert(wd.drawer_T35.cvs, [pyx.trafo.translate(PLOT_SPACING*6,  0)])
        # wd.drawer_T1.cvs.insert(wd.drawer_T37.cvs, [pyx.trafo.translate(PLOT_SPACING*7,  0)])
        wd.drawer_T1.cvs.insert(wd.drawer_T4.cvs,  [pyx.trafo.translate(             0, -50)])
        if wd.m==3:
            wd.drawer_T1.cvs.insert(wd.drawer_T4a.cvs, [pyx.trafo.translate(             0,-100)])
            wd.drawer_T1.cvs.insert(wd.drawer_T4b.cvs, [pyx.trafo.translate(             0,-150)])
            wd.drawer_T1.cvs.insert(wd.drawer_T4c.cvs, [pyx.trafo.translate(             0,-200)])
        # insert winding factor table
        wd.drawer_T1.cvs.insert(wd.drawer_Text.cvs, [pyx.trafo.translate(PLOT_SPACING*3, -220)])

        # Save collage as file
        wd.drawer_T1.cvs.writePDFfile(fname)
        print(f'save to {fname}')

        if wd.m!=3:
            quit()
        # continue

        # Q: 打出grouping和绕组函数波形，以及winding factor (including pitch factor)
        # Q: 绕组系数kw 和 绕组函数的DFT分量的幅值之间有什么关系？<-这个可以被用来检验你的谐波的短距系数有没有算对。
        # 短距系数计算？
        # 绕组连接图？
        # 按小博那样，直接顺便分析一下槽谐波啊！！！a voltage phasor diagram --Pyrhonen@p81



        # lazy me
        connection_star_raw_dict = wd.connection_star_raw_dict
        dpnv_grouping_dict_a = wd.dpnv_grouping_dict_a
        dpnv_grouping_dict_b = wd.dpnv_grouping_dict_b
        dpnv_grouping_dict_c = wd.dpnv_grouping_dict_c
        Q = wd.Q
        p = wd.p
        ps = wd.ps
        coil_pitch_y = wd.coil_pitch_y
        def reformat_wily_info(A,B,C,a,b,c, dpnv_grouping_dict_a, dpnv_grouping_dict_b, dpnv_grouping_dict_c):

            Q = len(A+B+C+a+b+c)
            layer_X_phases = [None]*Q
            layer_X_signs = [None]*Q

            for slot_number in A:
                layer_X_phases[slot_number-1] = 'U'
                layer_X_signs[slot_number-1] = '+'
            for slot_number in a:
                layer_X_phases[slot_number-1] = 'U'
                layer_X_signs[slot_number-1] = '-'
            for slot_number in B:
                layer_X_phases[slot_number-1] = 'V'
                layer_X_signs[slot_number-1] = '+'
            for slot_number in b:
                layer_X_phases[slot_number-1] = 'V'
                layer_X_signs[slot_number-1] = '-'
            for slot_number in C:
                layer_X_phases[slot_number-1] = 'W'
                layer_X_signs[slot_number-1] = '+'
            for slot_number in c:
                layer_X_phases[slot_number-1] = 'W'
                layer_X_signs[slot_number-1] = '-'

            grouping_AC = [0]*Q
            for slot_number in dpnv_grouping_dict_a['Aa']:
                grouping_AC[abs(int(slot_number))-1] = 1
            for slot_number in dpnv_grouping_dict_b['Bb']:
                grouping_AC[abs(int(slot_number))-1] = 1
            for slot_number in dpnv_grouping_dict_c['Cc']:
                grouping_AC[abs(int(slot_number))-1] = 1

            return layer_X_phases, layer_X_signs, grouping_AC

        print('\n---Winding definition in python:')
        # for key, val in connection_star_raw_dict.items():
        #     print('\t', key, [el[1] for el in val])
        # print('Group a/c:', dpnv_grouping_dict_a, dpnv_grouping_dict_b, dpnv_grouping_dict_c, sep='\t\n\t')
        
        if 'a' not in connection_star_raw_dict: # This is a special case when there is no emf vectors in zone X,Y,Z, e.g., Q6p2ps1
            layer_X_phases, layer_X_signs, grouping_AC = reformat_wily_info( 
                                [el[1] for el in connection_star_raw_dict['A']],
                                [el[1] for el in connection_star_raw_dict['B']],
                                [el[1] for el in connection_star_raw_dict['C']],
                                [],
                                [],
                                [],
                                dpnv_grouping_dict_a, dpnv_grouping_dict_b, dpnv_grouping_dict_c
                                )
        else:
            layer_X_phases, layer_X_signs, grouping_AC = reformat_wily_info( 
                                [el[1] for el in connection_star_raw_dict['A']],
                                [el[1] for el in connection_star_raw_dict['B']],
                                [el[1] for el in connection_star_raw_dict['C']],
                                [el[1] for el in connection_star_raw_dict['a']],
                                [el[1] for el in connection_star_raw_dict['b']],
                                [el[1] for el in connection_star_raw_dict['c']],
                                dpnv_grouping_dict_a, dpnv_grouping_dict_b, dpnv_grouping_dict_c
                                )
        print(r'''
            if DPNV_or_SEPA == True \
            and Qs == %d \
            and p == %d \
            and ps == %d \
            and coil_pitch_y == %d:
    '''%(Q,p,ps,coil_pitch_y), end='')
        print(r'''
                self.layer_X_phases = %s
                self.layer_X_signs  = %s
                self.coil_pitch_y   = coil_pitch_y
                self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
                self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)
    '''%(layer_X_phases, layer_X_signs), end='')
        print(r'''
                self.grouping_AC            = %s
                self.number_parallel_branch = %d
                self.number_winding_layer   = %d

                self.bool_3PhaseCurrentSource = False
                self.CommutatingSequenceD = 1
                self.CommutatingSequenceB = 0
    '''%(grouping_AC, 2, 2 if bool_double_layer_winding else 1), end='')

    return layer_X_phases, layer_X_signs, grouping_AC, coil_pitch_y

# https://pythoninformer.com/categories/category-pycairo/
# https://pythoninformer.com/python-libraries/pycairo/complex-shapes/
# https://pythoninformer.com/python-libraries/pycairo/drawing-shapes/
import cairo 
class winding_diagram:
    def __init__(self, layer_X_phases, layer_X_signs, coil_pitch_y, grouping_AC=None):
        self.grouping_AC = grouping_AC
        self.coil_pitch_y=coil_pitch_y

        self.layer_X_phases=layer_X_phases
        self.layer_X_signs=layer_X_signs

        self.layer_Y_phases = self.infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
        self.layer_Y_signs  = self.infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

        l_rightlayer1 = self.layer_X_phases
        l_rightlayer2 = self.layer_X_signs
        l_leftlayer1  = self.layer_Y_phases
        l_leftlayer2  = self.layer_Y_signs

        # dl: dict list
        self.dl_rightlayer = dl_rightlayer = {'U': [], 'V': [], 'W':[]}
        self.dl_leftlayer  = dl_leftlayer  = {'U': [], 'V': [], 'W':[]}
        count_slot = 0
        while True:
            try:
                count_slot += 1

                ABC     = l_rightlayer1.pop(0)
                up_down = l_rightlayer2.pop(0)
                for target in ['U', 'V', 'W']:
                    if ABC == target:
                        dl_rightlayer[target].append(up_down+str(count_slot))
                        break

                ABC     = l_leftlayer1.pop(0)
                up_down = l_leftlayer2.pop(0)
                for target in ['U', 'V', 'W']:
                    if ABC == target:
                        dl_leftlayer[target].append(up_down+str(count_slot))
                        break

            except IndexError as e: # nothing to pop
                # raise e
                count_slot -= 1
                break
        print('Qs:', count_slot)
        for k, v in dl_rightlayer.items():
            print('Right layer:', k, v)
        for k, v in dl_leftlayer.items():
            print('Left layer:', k, v)

        # dl_rightlayer['A'] = dl_rightlayer['U']
        # dl_leftlayer['A'] = dl_leftlayer['U']
        # dl_rightlayer['B'] = dl_rightlayer['V']
        # dl_leftlayer['B'] = dl_leftlayer['V']
        # dl_rightlayer['C'] = dl_rightlayer['W']
        # dl_leftlayer['C'] = dl_leftlayer['W']

        self.dl_grouping_AC = dict()
        self.dl_grouping_BD = dict()
        if self.grouping_AC is not None:
            for phase in ['U', 'V', 'W']:
                self.dl_grouping_AC[phase] = [el for el in dl_rightlayer[phase] if grouping_AC[abs(int(el))-1]==True ]
                self.dl_grouping_BD[phase] = [el for el in dl_rightlayer[phase] if grouping_AC[abs(int(el))-1]==False]
                print(f'Phase {phase}, Group A/C:', self.dl_grouping_AC[phase])
                print(f'         Group B/D:', self.dl_grouping_BD[phase])

    @staticmethod
    def infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(layer_X_phases, coil_pitch):
        return layer_X_phases[-coil_pitch:] + layer_X_phases[:-coil_pitch]

    @staticmethod
    def infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(layer_X_signs, coil_pitch):
        temp = layer_X_signs[-coil_pitch:] + layer_X_signs[:-coil_pitch]
        return [('-' if el == '+' else '+') for el in temp]

    def draw(self):

        slope = None
        for phase in ['U', 'V', 'W']:

            WIDTH       = 200
            HEIGHT      = 50
            PIXEL_SCALE = 10
            # with cairo.SVGSurface(f'{output_dir}tec-ismb-winding-diagram-{phase}.svg', WIDTH*PIXEL_SCALE, HEIGHT*PIXEL_SCALE) as surface:
            temp = self.dl_rightlayer['U']
            with cairo.SVGSurface(f'{output_dir}winding-diagram-{phase}-{"".join(temp)}.svg', WIDTH*PIXEL_SCALE, HEIGHT*PIXEL_SCALE) as surface:
                ctx = cairo.Context(surface)
                ctx.scale(PIXEL_SCALE, PIXEL_SCALE)
                ctx.set_source_rgba(0, 0, 0, 1)
                ctx.set_line_width(0.1)
                ctx.set_line_cap(cairo.LINE_CAP_ROUND)

                ctx.set_font_size(1)
                ctx.select_font_face("Times new roman", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

                Qs = len(self.dl_rightlayer['U'])*3
                # print(Qs, 'should be 36')

                COIL_LENGTH = 4 # 导体长度
                coil_spacing = 0.5 # 两个上层导体之间的间距
                coil_bias = 0.3 # 绘制双层绕组的第二层的偏移量
                THE_END   = (Qs  + 1 - 1/4.)*(1+coil_spacing)
                THE_SHIFT = (Qs) * (1+coil_spacing)

                CANVAS_OFFSET = COIL_LENGTH

                SLOPE_CONST = 2 # end-turn
                slope = SLOPE_CONST / (self.coil_pitch_y/2.*(1+coil_spacing) + coil_bias/2.)

                FOUR = COIL_LENGTH # 导体末端高度
                SIX = 2 + FOUR # 加上端部绕组后的高度

                ''' Draw coils in the slots as vertical lines for 3 phases
                '''
                for color, winding_layout_right, winding_layout_left, _ in zip( ['#5C9C31', '#E97675', '#6C8FAB'] , 
                                                                                [self.dl_rightlayer['U'], self.dl_rightlayer['V'], self.dl_rightlayer['W']],
                                                                                [self.dl_leftlayer['U'], self.dl_leftlayer['V'], self.dl_leftlayer['W']],
                                                                                list(range(3))):
                    if _ == 0:
                        ctx.set_source_rgba(92/255, 156/255, 49/255, 1.0)
                    elif _ == 1:
                        ctx.set_source_rgba(233/255, 118/255, 117/255, 1.0)
                    elif _ == 2:
                        ctx.set_source_rgba(108/255, 143/255, 171/255,  1.0)
                        # https://www.color-hex.com/color/6c8fab

                    count = 0
                    for el1, el2 in zip(winding_layout_right, winding_layout_left):
                        x1 = float(el1)
                        x2 = float(el2)
                        print(_, color, x1)

                        # def arrow_demo(ctx, x, y, width, height, a, b):
                        #     ctx.move_to(x, y + b)
                        #     ctx.line_to(x, y + height - b)
                        #     ctx.line_to(x + a, y + height - b)
                        #     ctx.line_to(x + a, y + height)
                        #     ctx.line_to(x + width, y + height/2)
                        #     ctx.line_to(x + a, y)
                        #     ctx.line_to(x + a, y + b)
                        #     ctx.close_path()
                        #     ctx.fill()

                        def varrow(ctx, P1, P2, arrow_width=0.2, arrow_height=0.5):
                            ctx.move_to(*P1)
                            ctx.line_to(*P2)
                            ctx.rel_line_to(-arrow_width, (-1 if (P2[1]>P1[1]) else 1) * arrow_height)
                            ctx.move_to(*P2)
                            ctx.rel_line_to(arrow_width,  (-1 if (P2[1]>P1[1]) else 1) * arrow_height)
                            ctx.set_line_join(cairo.LINE_JOIN_ROUND)
                            global_lw = ctx.get_line_width()
                            ctx.set_line_width(global_lw*0.8)
                            ctx.stroke()
                            ctx.set_line_width(global_lw)

                        # x means coil location (its meaning is equivalent to loc)
                        def draw_vertical_line_and_arrow_head(x, x2):

                            # 画上层导体
                            XU = abs(x)*(1+coil_spacing)
                            ctx.move_to(XU,           0+CANVAS_OFFSET)
                            ctx.line_to(XU, COIL_LENGTH+CANVAS_OFFSET)
                            ctx.set_dash([1, 0])
                            ctx.stroke()

                            # 画电流方向箭头
                            ARROW_BIAS = 1.0
                            annotate_bias = 0.5
                            if x>0:
                                xytext = [XU, 0.5*COIL_LENGTH + annotate_bias + CANVAS_OFFSET]
                                xy     = [XU, 0.5*COIL_LENGTH - annotate_bias + CANVAS_OFFSET]
                            else:
                                xy     = [XU, 0.5*COIL_LENGTH + annotate_bias + CANVAS_OFFSET]
                                xytext = [XU, 0.5*COIL_LENGTH - annotate_bias + CANVAS_OFFSET]
                            varrow(ctx, xy, xytext)

                            # 画下层导体
                            XL = abs(x2)*(1+coil_spacing) + coil_bias
                            ctx.move_to(XL,           0+CANVAS_OFFSET)
                            ctx.line_to(XL, COIL_LENGTH+CANVAS_OFFSET)
                            ctx.set_dash([0.5, 0.25])
                            ctx.stroke()

                            # 画数字123……
                            ctx.move_to(xy[0]+1.25*coil_bias, annotate_bias + 0.25 + CANVAS_OFFSET)
                            ctx.show_text(str(int(abs(x))))

                            # 画电流方向箭头
                            annotate_bias = 0.25
                            if x2>0:
                                xytext = [XL, 0.5*COIL_LENGTH + annotate_bias + ARROW_BIAS + CANVAS_OFFSET]
                                xy     = [XL, 0.5*COIL_LENGTH - annotate_bias + ARROW_BIAS + CANVAS_OFFSET]
                            else:
                                xy     = [XL, 0.5*COIL_LENGTH + annotate_bias + ARROW_BIAS + CANVAS_OFFSET]
                                xytext = [XL, 0.5*COIL_LENGTH - annotate_bias + ARROW_BIAS + CANVAS_OFFSET]
                            ctx.set_dash([1, 0])
                            varrow(ctx, xy, xytext)

                        def draw_end_turns(x, coil_pitch):

                            def mirror_plot(ctx, X, Y):
                                # 先画北极端部绕组
                                ctx.move_to(X[0], Y[0])
                                ctx.line_to(X[1], Y[1])
                                ctx.stroke()
                                # 再画南极端部绕组
                                ctx.move_to(X[0], Y[0] - 2*(Y[0] - COIL_LENGTH/2 - CANVAS_OFFSET) )
                                ctx.line_to(X[1], Y[1] - 2*(Y[1] - COIL_LENGTH/2 - CANVAS_OFFSET) )
                                ctx.stroke()
                                # print(Y[0], Y[1])

                            XU = abs(x)*(1+coil_spacing)

                            # going up
                            loc1 = (1+coil_spacing) * (abs(x) + coil_pitch/2) + coil_bias/2. # half pitch
                            # going down
                            loc2 = (1+coil_spacing) * (abs(x) + coil_pitch) + coil_bias # full pitch

                            if loc1 <= THE_END:
                                # going up
                                ctx.set_dash([1, 0])
                                mirror_plot(ctx,  [ XU, loc1], [ FOUR+CANVAS_OFFSET, 
                                                                SIX + CANVAS_OFFSET])
                                if loc2 < THE_END:
                                    # goibng down 
                                    ctx.set_dash([0.5, 0.3])
                                    mirror_plot(ctx,  [ loc1, loc2], [ SIX+CANVAS_OFFSET, 
                                                                    FOUR + CANVAS_OFFSET])
                                else:
                                    # going down (to be continued)
                                    ctx.set_dash([0.5, 0.3])
                                    mirror_plot(ctx,  [ loc1, THE_END], [ SIX+CANVAS_OFFSET, 
                                                   SIX-slope*(THE_END-loc1) + CANVAS_OFFSET])
                                    # going down (continued)
                                    mirror_plot(ctx,  [ THE_END-THE_SHIFT, loc2-THE_SHIFT], [ SIX-slope*(THE_END-loc1)+CANVAS_OFFSET, 
                                                                                                                FOUR + CANVAS_OFFSET])
                            else:

                                # going up (to be continued)
                                ctx.set_dash([1, 0])
                                mirror_plot(ctx,  [ XU, THE_END ], [ FOUR+CANVAS_OFFSET, 
                                                FOUR+slope*(THE_END-XU) + CANVAS_OFFSET])
                                # going up (continued)
                                mirror_plot(ctx,  [ THE_END-THE_SHIFT, loc1-THE_SHIFT], [ FOUR+slope*(THE_END-XU)+CANVAS_OFFSET, 
                                                                                                            SIX + CANVAS_OFFSET])

                                # going down 
                                ctx.set_dash([0.5, 0.3])
                                mirror_plot(ctx,  [ loc1-THE_SHIFT, loc2-THE_SHIFT], [ SIX+CANVAS_OFFSET, 
                                                                                    FOUR + CANVAS_OFFSET])

                        
                        draw_vertical_line_and_arrow_head(x1, x2) # right layer, left layer
                        draw_end_turns(x1, self.coil_pitch_y)

                        count += 1
                        # if count == 11:
                        #     break



                ''' Draw terminals for a phase
                '''
                if phase == 'U':
                    ctx.set_source_rgba(92/255, 156/255, 49/255, 1.0)
                elif phase == 'V':
                    ctx.set_source_rgba(233/255, 118/255, 117/255, 1.0)
                elif phase == 'W':
                    ctx.set_source_rgba(108/255, 143/255, 171/255,  1.0)

                if self.grouping_AC is not None:
                    def draw_terminals(slope, phase='U', coil_pitch=None):

                        terminal_bias = 0.5

                        self.dl_grouping_AC[phase] # ['-8', '+13', '+14', '-19', '-20', '+25']
                        self.dl_grouping_BD[phase] # ['+1', '+2', '-7', '+26', '-31', '-32']

                        def avoid_coil_distance_is_over_half_slots(grouping_number):
                            for index, str_number in enumerate(grouping_number):
                                if index+1 < len(grouping_number):
                                    if abs(int(grouping_number[index+1])) - abs(int(str_number)) > Qs/2:
                                        grouping_number = grouping_number[index+1:] + grouping_number[:index+1]
                                        print('Re-order as', grouping_number)
                            return grouping_number

                        self.dl_grouping_AC[phase] = avoid_coil_distance_is_over_half_slots(self.dl_grouping_AC[phase])
                        self.dl_grouping_BD[phase] = avoid_coil_distance_is_over_half_slots(self.dl_grouping_BD[phase])
                        print(self.dl_grouping_AC[phase])
                        print(self.dl_grouping_BD[phase])

                        for grouping_number, NINE_BIAS, grp in zip([self.dl_grouping_BD[phase], self.dl_grouping_AC[phase]],
                        # for grouping_number, NINE_BIAS, grp in zip(
                        #                                       (self.dl_grouping_BD[phase], self.dl_grouping_AC[phase][::-1]),
                                                              [1.0, 0.5],
                                                              ['b', 'a']):
                            HORIZONTAL_OFFSET = 0.0
                            THE_DIRECTION = None
                            NINE = SIX + NINE_BIAS

                            # Draw group A/C or B/D:
                            for index, str_number in enumerate(grouping_number):
                                x = int(str_number)
                                XU = abs(x)*(1+coil_spacing) # XU: location of upper layer

                                # Draw vertical coil terminal line
                                # Draw vertical coil terminal line
                                # Draw vertical coil terminal line
                                loc_triangle = XU+coil_pitch/2*(1+coil_spacing)+coil_bias/2 # top of the triangle
                                if x>0: 
                                    # solid line
                                    ctx.set_dash([1, 0])
                                    if loc_triangle - terminal_bias < THE_END: # 想画竖线，只要横坐标小于 THE_END，就不需要去掉 THE_SHIFT。
                                        端部连接处 = (loc_triangle - terminal_bias, NINE)
                                        ctx.move_to(端部连接处[0],            端部连接处[1] + CANVAS_OFFSET)
                                        ctx.line_to(端部连接处[0], SIX-terminal_bias*slope + CANVAS_OFFSET)
                                        ctx.stroke()
                                        # loc_triangle -= THE_SHIFT # 这里不要shift！
                                    else:
                                        端部连接处 = (loc_triangle - THE_SHIFT - terminal_bias, NINE)
                                        # y_span = [端部连接处[1]]+[SIX-terminal_bias*slope]
                                        # xy     = [端部连接处[0], sum(y_span)/2]
                                        # xytext = [端部连接处[0], sum(y_span)/2 + 0.25]
                                        ctx.move_to(端部连接处[0],            端部连接处[1] + CANVAS_OFFSET)
                                        ctx.line_to(端部连接处[0], SIX-terminal_bias*slope + CANVAS_OFFSET)
                                        ctx.stroke()

                                    # the other layer
                                    ctx.set_dash([0.0, 0.2])
                                    # ctx.move_to(端部连接处[0]+terminal_bias*2 - THE_SHIFT*(loc_triangle+terminal_bias>THE_END),            端部连接处[1] + CANVAS_OFFSET)
                                    # ctx.line_to(端部连接处[0]+terminal_bias*2 - THE_SHIFT*(loc_triangle+terminal_bias>THE_END), SIX-terminal_bias*slope + CANVAS_OFFSET)
                                    ctx.move_to(端部连接处[0]+terminal_bias*2,            端部连接处[1] + CANVAS_OFFSET)
                                    ctx.line_to(端部连接处[0]+terminal_bias*2, SIX-terminal_bias*slope + CANVAS_OFFSET)
                                    ctx.stroke()
                                    # ctx.set_dash([1, 0])

                                else: 
                                    # dashed line
                                    ctx.set_dash([0.0, 0.2])
                                    if loc_triangle + terminal_bias < THE_END:
                                        # loc_triangle -= THE_SHIFT
                                        端部连接处 = (loc_triangle + terminal_bias, NINE)
                                        ctx.move_to(端部连接处[0], 端部连接处[1] + CANVAS_OFFSET)
                                        ctx.line_to(端部连接处[0], SIX-terminal_bias*slope + CANVAS_OFFSET)
                                        ctx.stroke()
                                    else:
                                        端部连接处 = (loc_triangle - THE_SHIFT + terminal_bias, NINE)
                                        # y_span = [端部连接处[1]]+[SIX-terminal_bias*slope]
                                        # xy     = [端部连接处[0], sum(y_span)/2]
                                        # xytext = [端部连接处[0], sum(y_span)/2 + 0.25]
                                        ctx.move_to(端部连接处[0], 端部连接处[1] + CANVAS_OFFSET)
                                        ctx.line_to(端部连接处[0], SIX-terminal_bias*slope + CANVAS_OFFSET)
                                        ctx.stroke()

                                    # the other layer
                                    if loc_triangle - terminal_bias > THE_END:
                                        # 减了以后仍超过THE_END，简单，沿用之前已经shift过的端部连接处坐标即可。
                                        ctx.set_dash([1, 0])
                                        ctx.move_to(端部连接处[0]-terminal_bias*2,            端部连接处[1] + CANVAS_OFFSET)
                                        ctx.line_to(端部连接处[0]-terminal_bias*2, SIX-terminal_bias*slope + CANVAS_OFFSET)
                                        ctx.stroke()
                                        # ctx.set_dash([0.0, 0.2])
                                    else:
                                        # 减了以后没超THE_END，需要把shift给取消
                                        ctx.set_dash([1, 0])
                                        ctx.move_to(loc_triangle-terminal_bias,            端部连接处[1] + CANVAS_OFFSET)
                                        ctx.line_to(loc_triangle-terminal_bias, SIX-terminal_bias*slope + CANVAS_OFFSET)
                                        ctx.stroke()

                                # print('-----------', x)
                                # ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->", color=color, lw=0.5, alpha=1.0))


                                # 其实这里的代码默默保证了线圈的端部接线是保证电流从左向右流动的。
                                # 画线圈组的出头
                                # 画线圈组的出头
                                # 画线圈组的出头
                                if index == 0 or index+1 >= len(grouping_number):

                                    if index == 0:
                                        if x>0:
                                            ctx.set_dash([0.0, 0.2])
                                        else:
                                            ctx.set_dash([1, 0])

                                        XTerminal = loc_triangle + terminal_bias * (1 if x > 0 else -1)
                                        if XTerminal > THE_END:
                                            XTerminal -= THE_SHIFT
                                        ctx.move_to(XTerminal,           1.0 + NINE + CANVAS_OFFSET)
                                        ctx.line_to(XTerminal,                 NINE + CANVAS_OFFSET)
                                        ctx.stroke()

                                        ctx.set_dash([1, 0])
                                        ctx.arc(    XTerminal,           1.0 + NINE + CANVAS_OFFSET, 0.2, 0, 2*np.pi)
                                        ctx.rel_move_to(0, 0.5)
                                        ctx.show_text(f'{phase.lower()}{grp}{"+"}')
                                        ctx.stroke()
                                    else:
                                        if x<0:
                                            ctx.set_dash([0.0, 0.2])
                                        else:
                                            ctx.set_dash([1, 0])

                                        XTerminal = loc_triangle - terminal_bias * (1 if x > 0 else -1)
                                        if XTerminal > THE_END:
                                            XTerminal -= THE_SHIFT
                                        ctx.move_to(XTerminal,           1.0 + NINE + CANVAS_OFFSET)
                                        ctx.line_to(XTerminal,                 NINE + CANVAS_OFFSET)
                                        ctx.stroke()

                                        ctx.set_dash([1, 0])
                                        ctx.arc(    XTerminal,           1.0 + NINE + CANVAS_OFFSET, 0.2, 0, 2*np.pi)
                                        ctx.rel_move_to(0, 0.5)
                                        ctx.show_text(f'{phase.lower()}{grp}{"-"}')
                                        ctx.stroke()

                                # 其实这里的代码默默保证了线圈的端部接线是保证电流从左向右流动的。
                                # Draw horizontal line that connects coil to coil
                                # Draw horizontal line that connects coil to coil
                                # Draw horizontal line that connects coil to coil
                                if index+1 < len(grouping_number):
                                    x2 = int(grouping_number[index+1])

                                    if x2 < 0:
                                        # 如果x2是负的，就不需要垫高
                                        HORIZONTAL_OFFSET = 0.0
                                        ctx.move_to(端部连接处[0], NINE + CANVAS_OFFSET)
                                    else:
                                        # 如果x2是正的，则需要垫高
                                        HORIZONTAL_OFFSET += 0.5
                                        ctx.move_to(端部连接处[0], NINE + CANVAS_OFFSET)
                                        ctx.line_to(端部连接处[0], NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)

                                    # if abs(x2) != 35 and abs(x2) != 36:
                                    # if abs(x2) != 36:
                                    if True:
                                        loc_triangle2 = (abs(x2)+coil_pitch/2)*(1+coil_spacing)+coil_bias/2 # top of the triangle
                                        print(phase, x, x2)
                                        def draw_continued(shift):
                                            ctx.line_to(THE_END, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)

                                            ctx.move_to(THE_END - THE_SHIFT, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)

                                            if HORIZONTAL_OFFSET > 0:
                                                ctx.line_to(loc_triangle2 - shift + terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)
                                                ctx.line_to(loc_triangle2 - shift + terminal_bias, NINE + CANVAS_OFFSET)
                                            else:
                                                ctx.line_to(loc_triangle2 - shift - terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)

                                        # if abs(x2) > abs(x):
                                        if loc_triangle2 > loc_triangle:
                                            if abs(x2) == 36:
                                                print(x, x2, loc_triangle, loc_triangle2)

                                            if (loc_triangle > THE_END and loc_triangle2 > THE_END) \
                                            or (loc_triangle+terminal_bias > THE_END and loc_triangle2 > THE_END): # 为了解决样机绕组的V相的20线圈连21线圈时会出现的问题
                                                print("Both x's and x2's triangle locations are larger than THE_END.", x, x2)
                                                # 两个线圈的三角尖尖都超过了THE_END，比如x=-35，x2=-36
                                                # ctx.rel_move_to(-THE_SHIFT, 0)

                                                if HORIZONTAL_OFFSET > 0:
                                                    ctx.line_to(loc_triangle2-THE_SHIFT + terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)
                                                    ctx.line_to(loc_triangle2-THE_SHIFT + terminal_bias, NINE + CANVAS_OFFSET)
                                                else:
                                                    ctx.line_to(loc_triangle2-THE_SHIFT - terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)
                                            else:
                                                if loc_triangle2 < THE_END:
                                                    # 此时到了x2处，
                                                    if HORIZONTAL_OFFSET > 0:
                                                        ctx.line_to(loc_triangle2 + terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)
                                                        ctx.line_to(loc_triangle2 + terminal_bias, NINE + CANVAS_OFFSET)
                                                    else:
                                                        ctx.line_to(loc_triangle2 - terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)
                                                else:
                                                    # 需要向右延伸截断并从最左边续上
                                                    draw_continued(THE_SHIFT)
                                        else:
                                            #  x很大，比如24，而x很小，比如1。 所以 x2=1 < x=24
                                            if loc_triangle > THE_END and loc_triangle2 < THE_END:

                                                # 24的loc_triangle肯定超过THE_END了，但是1的没有loc_triangle2很小
                                                # ctx.rel_move_to(-THE_SHIFT, 0)

                                                if HORIZONTAL_OFFSET > 0:
                                                    ctx.line_to(loc_triangle2 + terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)
                                                    ctx.line_to(loc_triangle2 + terminal_bias, NINE + CANVAS_OFFSET)
                                                else:
                                                    ctx.line_to(loc_triangle2 - terminal_bias, NINE + CANVAS_OFFSET + HORIZONTAL_OFFSET)
                                            else:
                                                # 需要向右延伸截断并从最左边续上
                                                draw_continued(0)

                                    ctx.set_dash([0.0, 0.2])
                                    ctx.stroke()
                            # break # draw only group ac

                    draw_terminals(slope, phase=phase, coil_pitch=self.coil_pitch_y)
                    # break

if __name__ == '__main__':

    if True:
        phases, signs, grouping_AC, coil_pitch_y = main_derivation()
        quit()
    else:
        phases = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V']
        signs = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
        grouping_AC = [0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1]
        coil_pitch_y = 5

    diagram = winding_diagram(
        layer_X_phases = phases,
        layer_X_signs = signs,
        coil_pitch_y = coil_pitch_y,
        grouping_AC = grouping_AC,
    )
    diagram.draw()
