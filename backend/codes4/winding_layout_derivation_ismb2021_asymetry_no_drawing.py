"""
绕组布局推导脚本（无绘图版本）
仅计算绕组信息，不进行任何绘图操作
"""
ABCDEFGHIJKLMNOPQRSTUVWXYZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

from pylab import np
import math
from utility import gcd

# 仅保留计算所需的全局变量
BELT_BIAS = 5  # deg. elec. (用于相位带计算)

# 全局 verbose 标志（用于控制是否打印输出）
_verbose = False


def _print(*args, **kwargs):
    """条件打印函数：仅在 verbose=True 时打印"""
    if _verbose:
        print(*args, **kwargs)


def set_verbose(verbose=True):
    """设置全局 verbose 标志"""
    global _verbose
    _verbose = verbose


def limit_to_360_deg(PHI):
    while PHI < 0:
        PHI += 360
    while PHI >= 360:
        PHI -= 360
    return PHI


def belong_to_which_phase_belt(PHI, phase_belt):
    PHI = limit_to_360_deg(PHI)

    if phase_belt == 60:
        if PHI <= phase_belt*0.5 + BELT_BIAS or PHI > 360 - phase_belt*0.5 + BELT_BIAS:
            return 'A'  # 'u+'
        elif 180 - phase_belt*0.5 + BELT_BIAS < PHI <= 180 + phase_belt*0.5 + BELT_BIAS:
            return 'a'  # 'u-'
        elif 120 - phase_belt*0.5 + BELT_BIAS < PHI <= 120 + phase_belt*0.5 + BELT_BIAS:
            return 'B'  # 'v+'
        elif 300 - phase_belt*0.5 + BELT_BIAS < PHI <= 300 + phase_belt*0.5 + BELT_BIAS:
            return 'b'  # 'v-'
        elif 240 - phase_belt*0.5 + BELT_BIAS < PHI <= 240 + phase_belt*0.5 + BELT_BIAS:
            return 'C'  # 'w+'
        elif 60 - phase_belt*0.5 + BELT_BIAS < PHI <= 60 + phase_belt*0.5 + BELT_BIAS:
            return 'c'  # 'w-'
        else:
            raise Exception('Unexpected PHI=%g' % (PHI))
    elif phase_belt == 120:
        if PHI <= phase_belt*0.5 + BELT_BIAS or PHI > 360 - phase_belt*0.5 + BELT_BIAS:
            return 'A'  # 'u+'
        elif 120 - phase_belt*0.5 + BELT_BIAS < PHI <= 120 + phase_belt*0.5 + BELT_BIAS:
            return 'B'  # 'v+'
        elif 240 - phase_belt*0.5 + BELT_BIAS < PHI <= 240 + phase_belt*0.5 + BELT_BIAS:
            return 'C'  # 'w+'
        else:
            raise Exception('Unexpected PHI=%g' % (PHI))
    else:
        phase_name_positive_list = ABCDEFGHIJKLMNOPQRSTUVWXYZ
        phase_name_negative_list = ABCDEFGHIJKLMNOPQRSTUVWXYZ.lower()

        i = 0
        if PHI <= phase_belt*0.5 + BELT_BIAS or PHI > 360 - phase_belt*0.5 + BELT_BIAS:
            return phase_name_positive_list[0]  # 'u+'
        elif 180 - phase_belt*0.5 + BELT_BIAS < PHI <= 180 + phase_belt*0.5 + BELT_BIAS:
            return phase_name_negative_list[0]  # 'u-'
        else:
            while True:
                i += 1
                if i*2*phase_belt - phase_belt*0.5 + BELT_BIAS < PHI <= i*2*phase_belt + phase_belt*0.5 + BELT_BIAS:
                    return phase_name_positive_list[i]
                elif i*2*phase_belt+180 - phase_belt*0.5 + BELT_BIAS < PHI <= i*2*phase_belt+180 + phase_belt*0.5 + BELT_BIAS:
                    return phase_name_negative_list[i]
                if i > 100:
                    raise Exception('Dead loop.')


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


def phase_angle_of_slot_i_at_frequency_h(slot_number, h, Q):
    """计算槽 i 在频率 h 下的相位角"""
    deg_elec_angle_between_adjacent_slots = 2*math.pi * h / Q
    return deg_elec_angle_between_adjacent_slots * (slot_number-1)


def compute_star_of_slots(Q, p, m, verbose=None):
    """
    计算槽电势星形图数据（不绘图）
    返回 connection_star_raw_dict 和 phase_belt
    
    Args:
        verbose: 是否打印输出，None 时使用全局 _verbose 设置
    """
    if verbose is None:
        verbose = _verbose
    
    电角度 = electrical_degree = 360 * p
    槽距角 = deg_elec_angle_between_adjacent_slots = 360 * p / Q
    相带 = phase_belt = 360/(2*m)  # 60 or 120 deg. elec.
    
    if verbose:
        _print('\t Number of Zones = 2*p*m = %d' % (2*p*m))
        _print('\t 槽距角：\t', deg_elec_angle_between_adjacent_slots, 'deg. elec.')
        _print('\t 相带：\t', phase_belt, 'deg. elec.')

    connection_star_raw_dict = dict()

    # 计算每个槽的相位带归属
    for i in range(Q):
        PHI = i*deg_elec_angle_between_adjacent_slots
        # PHI 不能归化到360°以内，否则会影响后面group a/c分组的判断
        key = belong_to_which_phase_belt(PHI, phase_belt)
        val = (PHI, i+1)
        if verbose:
            _print(key, val)
        if key in connection_star_raw_dict:
            connection_star_raw_dict[key].append(val)
        else:
            connection_star_raw_dict[key] = [val]

    return connection_star_raw_dict, phase_belt


def compute_connection_star_at_another_frequency(connection_star_raw_dict, frequency_ratio, which_phase='Aa', verbose=None):
    """
    计算在另一个频率下的连接星形图分组（不绘图）
    返回 dpnv_grouping_dict
    
    Args:
        verbose: 是否打印输出，None 时使用全局 _verbose 设置
    """
    if verbose is None:
        verbose = _verbose
    
    dpnv_grouping_dict = dict()
    dpnv_grouping_dict['GAC'] = []
    dpnv_grouping_dict['GBD'] = []

    # 计算 180e 带
    if which_phase == 'Aa':
        band_phase_shift = -0
    elif which_phase == 'Bb':
        band_phase_shift = -120
    elif which_phase == 'Cc':
        band_phase_shift = -240
    else:
        raise Exception('Wrong which_phase: %s' % (which_phase))
    
    LB = 90 - BELT_BIAS + band_phase_shift
    UB = 180 + 90 - BELT_BIAS + band_phase_shift

    # 计算分组
    for key, val in connection_star_raw_dict.items():
        if key in which_phase:
            if verbose:
                _print('||||||', key, val)

            if key in 'abc':
                factor_reverse = -1
                phase_shift = 180
            else:  # in 'ABC'
                factor_reverse = 1
                phase_shift = 0

            for el in val:
                PHI_ori = el[0] * frequency_ratio  # PHI at the new frequency
                PHI, label = phase_shift+PHI_ori, str(factor_reverse*el[1])

                # Grouping a/c
                if belong_to_band(LB, UB, PHI):
                    if verbose:
                        _print('Group a/c:', label, LB, UB, limit_to_360_deg(PHI))
                    dpnv_grouping_dict['GAC'] += [label]
                else:
                    if verbose:
                        _print('Group b/d:', label, LB, UB, limit_to_360_deg(PHI))
                    dpnv_grouping_dict['GBD'] += [label]

    # backward compatible
    dpnv_grouping_dict[which_phase] = dpnv_grouping_dict['GAC']
    return dpnv_grouping_dict


def winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, phase_Aa_dpnv_grouping_dict=None, Aa='Aa'):
    """计算绕组分布因子"""
    list_slots_of_a_phase = []
    list_phase_shift = []
    
    for key, val in connection_star_raw_dict.items():
        if key in Aa:
            list_slots_of_a_phase += [el[-1] for el in val]
            if key in Aa[0]:
                list_phase_shift += [0]*len(val)
            elif key in Aa[1]:
                list_phase_shift += [math.pi]*len(val)

    # phasors in DPNV group a/c will be flipped (phase shifted 180^e)
    if phase_Aa_dpnv_grouping_dict is not None:
        phase_Aa_dpnv_grouping_list = phase_Aa_dpnv_grouping_dict[Aa]
        reversed_excitation_upper_layer = [abs(int(el)) for el in phase_Aa_dpnv_grouping_list]
        for idx, slot_number in enumerate(list_slots_of_a_phase):
            if slot_number in reversed_excitation_upper_layer:
                list_phase_shift[idx] += math.pi

    if bool_double_layer_winding == True:
        number_of_coils_per_phase = len(list_slots_of_a_phase)
        N = number_of_coils_per_phase  # N = Q/m
    else:
        number_of_coil_sides_per_phase = len(list_slots_of_a_phase)
        N = number_of_coil_sides_per_phase  # N = Q/m

    winding_distribution_factor_kw_at_h = abs(
        sum([
            np.exp(1j * (phase_angle_of_slot_i_at_frequency_h(slot_number, h, Q)+phase_shift))
            for slot_number, phase_shift in zip(list_slots_of_a_phase, list_phase_shift)
        ])
    ) / N
    return winding_distribution_factor_kw_at_h


def winding_short_pitch_factor_v2(h, coil_pitch_y, Q):
    """计算短距因子"""
    y_Q = Q / (2)
    k_ph = math.sin(h * coil_pitch_y/y_Q * math.pi*0.5)
    return k_ph


class Winding_Derivation(object):
    """绕组推导类（无绘图版本）"""
    
    def __init__(self, slot_pole_comb, bool_double_layer_winding=True, verbose=None):
        """
        Args:
            slot_pole_comb: (m, Q, p, ps, coil_pitch_y, turn_func_bias)
            bool_double_layer_winding: 是否为双层绕组
            verbose: 是否打印输出，None 时使用全局 _verbose 设置
        """
        if verbose is None:
            verbose = _verbose
        self.verbose = verbose
        
        # General info
        self.m = m = slot_pole_comb[0]
        self.Q = Q = slot_pole_comb[1]
        self.p = p = slot_pole_comb[2]
        self.ps = ps = slot_pole_comb[3]
        self.coil_pitch_y = coil_pitch_y = slot_pole_comb[4]
        self.turn_func_bias = turn_func_bias = slot_pole_comb[5]

        self.bool_double_layer_winding = bool_double_layer_winding

        self.t = t = gcd(Q, p)
        self.ts = ts = gcd(Q, ps)
        self.q = q = Q/(2*p)/m
        self.qs = qs = Q/(2*ps)/m
        
        if verbose:
            _print(f'Q={Q} \np={p} \nps={ps} \nt={t} \nts={ts} \nq={q} \nqs={qs}.')
        
        if p % 2 == 0:
            if verbose:
                _print('DPNV Type 1.')
        elif p % 2 == 1 and ps % 2 == 0:
            if verbose:
                _print('DPNV Type 2.')
        else:
            if verbose:
                _print('DPNV Type 3.')
        
        if q <= 1 and verbose:
            _print('This winding has no distribution factor, right?')

        # 计算槽电势星形图
        if verbose:
            _print(f'Torque star of slots with Q={Q} and p={p}.')
        connection_star_raw_dict, phase_belt = compute_star_of_slots(Q, p, m, verbose=verbose)
        self.connection_star_raw_dict = connection_star_raw_dict

        # 计算各相槽号
        def get_list_phase_slot_number(which_phase_belt_A='A', which_phase_belt_a='a'):
            list_phase_u_slot_number = [[el[-1] for el in val] for key, val in connection_star_raw_dict.items() if key in which_phase_belt_A]
            list_phase_u_slot_number += [[-el[-1] for el in val] for key, val in connection_star_raw_dict.items() if key in which_phase_belt_a]
            list_phase_u_slot_number = [item for sublist in list_phase_u_slot_number for item in sublist]
            return list_phase_u_slot_number

        if m == 3:
            self.list_phase_u_slot_number = list_phase_u_slot_number = get_list_phase_slot_number('A', 'a')
            if verbose:
                _print('Phase U:', ', '.join([str(el) for el in list_phase_u_slot_number]))
            self.list_phase_v_slot_number = list_phase_v_slot_number = get_list_phase_slot_number('B', 'b')
            if verbose:
                _print('Phase V:', ', '.join([str(el) for el in list_phase_v_slot_number]))
            self.list_phase_w_slot_number = list_phase_w_slot_number = get_list_phase_slot_number('C', 'c')
            if verbose:
                _print('Phase W:', ', '.join([str(el) for el in list_phase_w_slot_number]))
        else:
            # solution for m>3
            self.list_slot_number_of_phase = dict()
            for _phaseNumber in range(m):
                _phaseName = ABCDEFGHIJKLMNOPQRSTUVWXYZ[_phaseNumber]
                self.list_slot_number_of_phase[_phaseName] = get_list_phase_slot_number(_phaseName, _phaseName.lower())
                if verbose:
                    _print(f'Phase {_phaseName}:', ', '.join([str(el) for el in self.list_slot_number_of_phase[_phaseName]]))

        # 计算 DPNV 分组
        if m == 3:
            dpnv_grouping_dict_a = compute_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*ps, which_phase='Aa', verbose=verbose)
            dpnv_grouping_dict_b = compute_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*ps, which_phase='Bb', verbose=verbose)
            dpnv_grouping_dict_c = compute_connection_star_at_another_frequency(connection_star_raw_dict, 1/p*ps, which_phase='Cc', verbose=verbose)
            self.dpnv_grouping_dict_a = dpnv_grouping_dict_a
            self.dpnv_grouping_dict_b = dpnv_grouping_dict_b
            self.dpnv_grouping_dict_c = dpnv_grouping_dict_c
            if verbose:
                _print('*Group a/c (phase U):', ', '.join(dpnv_grouping_dict_a['Aa']))
                _print('*Group a/c (phase W):', ', '.join(dpnv_grouping_dict_b['Bb']))
                _print('*Group a/c (phase V):', ', '.join(dpnv_grouping_dict_c['Cc']))

        # 计算绕组因子
        kw_at_p = winding_distribution_factor(Q, connection_star_raw_dict, p, bool_double_layer_winding)
        kw_at_ps = winding_distribution_factor(Q, connection_star_raw_dict, ps, bool_double_layer_winding)
        if verbose:
            _print('- Winding distribution factor at p (=%d): %g.' % (p, kw_at_p))
            _print('- Winding distribution factor at ps (=%d): %g.' % (ps, kw_at_ps))
        
        if m == 3:
            kw_at_ps = winding_distribution_factor(Q, connection_star_raw_dict, ps, bool_double_layer_winding, 
                                                   phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_a, Aa='Aa')
            if verbose:
                _print('- Winding distribution factor at ps (=%d) with flipped phasors in group a/c: %g.' % (ps, kw_at_ps))

        if verbose:
            _print('- Winding distribution factor at h:')
        count_TW = count_SW = 0
        self.torque_kd_at_h = dict()
        self.suspen_kd_at_h = dict()
        
        for h in range(25):
            temp_TW = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding)
            if m == 3:
                temp_SW_A = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, 
                                                        phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_a, Aa='Aa')
                temp_SW_B = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, 
                                                        phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_b, Aa='Bb')
                temp_SW_C = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, 
                                                        phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_c, Aa='Cc')
            
            if abs(temp_TW) > 1e-4:
                msg = '\t%d: %.3f' % (h, temp_TW)
                if verbose:
                    _print('\tT:', msg)
                count_TW += 1
                self.torque_kd_at_h[h] = temp_TW
            
            if m == 3:
                if abs(temp_SW_A) > 1e-4:
                    msg = '\t%d: %.2f,\t%.2f,\t%.2f' % (h, temp_SW_A, temp_SW_B, temp_SW_C)
                    if verbose:
                        _print('\tS:', msg)
                    count_SW += 1
                    self.suspen_kd_at_h[h] = temp_SW_A, temp_SW_B, temp_SW_C

        if verbose:
            _print('- Winding short pitch factor at h:')
        count_TW = count_SW = 0
        self.torque_kp_at_h = dict()
        self.suspen_kp_at_h = dict()
        
        for h in range(25):
            temp_TW = winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            temp_SW = winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            
            if abs(temp_TW) > 1e-4:
                msg = '\t%d: %.3f' % (h, temp_TW)
                if verbose:
                    _print('\t\tT:', msg)
                count_TW += 1
                self.torque_kp_at_h[h] = temp_TW
            
            if abs(temp_SW) > 1e-4:
                msg = '\t%d: %.3f' % (h, temp_SW)
                if verbose:
                    _print('\t\tS:', msg)
                count_SW += 1
                self.suspen_kp_at_h[h] = temp_SW

        if verbose:
            _print('- Winding factor at h:')
        count_TW = count_SW = 0
        self.torque_kw_at_h = dict()
        self.suspen_kw_at_h = dict()
        
        for h in range(25):
            temp_TW = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding)
            temp_TW *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            
            if m == 3:
                temp_SW_A = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, 
                                                        phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_a, Aa='Aa')
                temp_SW_B = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, 
                                                        phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_b, Aa='Bb')
                temp_SW_C = winding_distribution_factor(Q, connection_star_raw_dict, h, bool_double_layer_winding, 
                                                        phase_Aa_dpnv_grouping_dict=dpnv_grouping_dict_c, Aa='Cc')
                temp_SW_A *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
                temp_SW_B *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
                temp_SW_C *= winding_short_pitch_factor_v2(h, coil_pitch_y, Q)
            
            if abs(temp_TW) > 1e-4:
                msg = '\t%d: %.3f' % (h, temp_TW)
                if verbose:
                    _print(msg)
                count_TW += 1
                self.torque_kw_at_h[h] = temp_TW
            
            if m == 3:
                if abs(temp_SW_A) > 1e-4:
                    msg = '\t%d: %.2f,\t%.2f,\t%.2f' % (h, temp_SW_A, temp_SW_B, temp_SW_C)
                    if verbose:
                        _print('\t\t\tS:', msg)
                    count_SW += 1
                    self.suspen_kw_at_h[h] = temp_SW_A, temp_SW_B, temp_SW_C

    def get_complex_number_winding_factor_of_coil_i(self, i, coil_pitch_y, Q, v, p):
        """计算线圈 i 的复数绕组因子"""
        alpha_u = 2*math.pi / Q * p
        gamma = coil_pitch_y*alpha_u/p
        radii = math.sin(v*p*coil_pitch_y*alpha_u/p/2)
        if self.verbose:
            _print(f'Coil span [mech.deg] = {gamma/math.pi*180} | [elec.deg] = {v*p*gamma / math.pi*180}', end=' | ')
        
        ELS_angles = -0.5*math.pi - v*p*alpha_u*(2*i+coil_pitch_y) / (2*p)
        CJH_angles = 0.5*math.pi - v*p*alpha_u*(2*i+coil_pitch_y) / (2*p)
        ELS_pitch_factor_per_coil = radii * np.exp(1j*ELS_angles)
        CJH_pitch_factor_per_coil = radii * np.exp(1j*CJH_angles)
        return ELS_pitch_factor_per_coil, CJH_pitch_factor_per_coil

    def get_complex_number_kw_per_phase(self, v, p, positive_connected_coils, negative_connected_coils):
        """计算每相的复数绕组因子"""
        kp_els_list = []
        kp_cjh_list = []

        for i in positive_connected_coils:
            els, cjh = self.get_complex_number_winding_factor_of_coil_i(i, self.coil_pitch_y, Q=self.Q, v=v, p=p)
            kp_els_list.append(els)
            kp_cjh_list.append(cjh)
            if self.verbose:
                _print(f'\tkp@Coil+{i:02d}', end='\t|\t')
                _print('cjh = %.3f∠%.1f' % (np.abs(cjh), np.angle(cjh)/math.pi*180*1))

        for i in negative_connected_coils:
            els, cjh = self.get_complex_number_winding_factor_of_coil_i(i, self.coil_pitch_y, Q=self.Q, v=v, p=p)
            els *= -1
            cjh *= -1
            kp_els_list.append(els)
            kp_cjh_list.append(cjh)
            if self.verbose:
                _print(f'\tkp@Coil-{i:02d}', end='\t|\t')
                _print('cjh = %.3f∠%.1f' % (np.abs(cjh), np.angle(cjh)/math.pi*180*1))

        average = lambda x: np.sum(x)/len(x)
        _kw_els = average(kp_els_list)
        _kw_cjh = average(kp_cjh_list)
        return _kw_els, _kw_cjh

    def get_complex_number_kw(self, p_or_ps, v=1, bool_study_suspension_subharmonics=False):
        """计算复数绕组因子"""
        dict_kw_els = dict()
        dict_kw_cjh = dict()

        for ZONE in ['A', 'B', 'C']:
            pcc = [i for _angle, i in self.connection_star_raw_dict[ZONE]]
            try:
                ncc = [i for _angle, i in self.connection_star_raw_dict[ZONE.lower()]]
            except KeyError:
                if self.verbose:
                    _print('负60度相带里没有槽电势矢量，所以abc不存在于connection_star_raw_dict.keys()。')
                ncc = []

            if p_or_ps == self.ps or bool_study_suspension_subharmonics:
                # suspension winding has different connection pattern from the torque winding
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
                pcc = [el for el in connection_star_raw_results if el > 0]
                ncc = [abs(el) for el in connection_star_raw_results if el < 0]
                if self.verbose:
                    _print('sus-pcc:', pcc)
                    _print('sus-ncc:', ncc)

            _kw_els, _kw_cjh = self.get_complex_number_kw_per_phase(v=v, p=p_or_ps, positive_connected_coils=pcc, negative_connected_coils=ncc)
            dict_kw_els[f'{ZONE}'] = _kw_els
            dict_kw_cjh[f'{ZONE}'] = _kw_cjh
            dict_kw_els[f'{ZONE}_abs'] = np.abs(_kw_els)
            dict_kw_cjh[f'{ZONE}_abs'] = np.abs(_kw_cjh)
            dict_kw_els[f'{ZONE}_angle'] = np.angle(_kw_els)/math.pi*180
            dict_kw_cjh[f'{ZONE}_angle'] = np.angle(_kw_cjh)/math.pi*180
            if self.verbose:
                _print('kw(cjh)=', _kw_cjh, '=', f'{np.abs(_kw_cjh):.3f}∠{np.angle(_kw_cjh)/math.pi*180:.1f}')

        return dict_kw_els, dict_kw_cjh


    def format_print_out_string(self):

        ''' NEW ISMB 2021 Complex Number Winding Factor '''
        if self.verbose:
            _print('\n----------------------------------n=p')
        self.dict_torque_kw_els, self.dict_torque_kw_cjh = self.get_complex_number_kw(self.p)
        if self.verbose:
            _print('\n----------------------------------n=ps')
        self.dict_suspension_kw_els, self.dict_suspension_kw_cjh = self.get_complex_number_kw(self.ps, v=1)
        
        phase_difference = [
            self.dict_suspension_kw_els['A_angle'] - self.dict_suspension_kw_els['B_angle'],
            self.dict_suspension_kw_els['B_angle'] - self.dict_suspension_kw_els['C_angle'],
            self.dict_suspension_kw_els['C_angle'] - self.dict_suspension_kw_els['A_angle'],
        ]
        if self.verbose:
            _print('phases of winding:', self.dict_suspension_kw_els['A_angle'], self.dict_suspension_kw_els['B_angle'], self.dict_suspension_kw_els['C_angle'])
            _print('phase_difference:', phase_difference)
        for el in phase_difference:
            if abs(abs(el) - 240.0) > 1e-5 and abs(abs(el) - 120.0) > 1e-5:
                if self.verbose:
                    _print('!!!Phase Asymmetry Detected')

        # 输出绕组定义信息
        connection_star_raw_dict = self.connection_star_raw_dict
        if self.m == 3:
            dpnv_grouping_dict_a = self.dpnv_grouping_dict_a
            dpnv_grouping_dict_b = self.dpnv_grouping_dict_b
            dpnv_grouping_dict_c = self.dpnv_grouping_dict_c
        else:
            dpnv_grouping_dict_a = None
            dpnv_grouping_dict_b = None
            dpnv_grouping_dict_c = None
        Q = self.Q
        p = self.p
        ps = self.ps
        coil_pitch_y = self.coil_pitch_y

        def reformat_wily_info(A, B, C, a, b, c, dpnv_grouping_dict_a, dpnv_grouping_dict_b, dpnv_grouping_dict_c):
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

        if self.verbose:
            _print('\n---Winding definition in python:')
        
        if 'a' not in connection_star_raw_dict:
            layer_X_phases, layer_X_signs, grouping_AC = reformat_wily_info(
                [el[1] for el in connection_star_raw_dict['A']],
                [el[1] for el in connection_star_raw_dict['B']],
                [el[1] for el in connection_star_raw_dict['C']],
                [], [], [],
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
        
        self.print_out_string = (
            'if DPNV_or_SEPA == True \\\n'
            'and Qs == %d \\\n'
            'and p == %d \\\n'
            'and ps == %d \\\n'
            'and coil_pitch_y == %d:\n'
            '    self.layer_X_phases = %s\n'
            '    self.layer_X_signs  = %s\n'
            '    self.coil_pitch_y   = coil_pitch_y\n'
            '    self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)\n'
            '    self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)\n'
            '    self.grouping_AC            = %s\n'
            '    self.number_parallel_branch = %d\n'
            '    self.number_winding_layer   = %d\n'
            '\n'
            '    self.bool_3PhaseCurrentSource = False\n'
            '    self.CommutatingSequenceD = 1\n'
            '    self.CommutatingSequenceB = 0\n'
            % (Q, p, ps, coil_pitch_y, layer_X_phases, layer_X_signs, grouping_AC, 2, 2 if self.bool_double_layer_winding else 1)
        )
        if self.verbose:
            _print(self.print_out_string, end='')

        self.layer_X_phases = layer_X_phases
        self.layer_X_signs = layer_X_signs
        self.grouping_AC = grouping_AC
        self.coil_pitch_y = coil_pitch_y

        return self.print_out_string

def main_derivation(m, Qs, p, ps, coil_pitch_y, verbose=None):
    if verbose is None:
        verbose = _verbose
    
    # m, Q, p, ps, y, turn function bias (turn_func_bias)
    Slot_Pole_Combinations = [
        (m, Qs, p, ps, coil_pitch_y, 0),
        # (3, 24, 4, 5, 1, 0),  # Q24p4ps5
        # 可以添加更多组合
    ]
    bool_double_layer_winding = True

    for index, slot_pole_comb in enumerate(Slot_Pole_Combinations):
        wd = Winding_Derivation(slot_pole_comb, bool_double_layer_winding, verbose=verbose)
        wd.format_print_out_string()
        return wd

if __name__ == '__main__':
    main_derivation()

