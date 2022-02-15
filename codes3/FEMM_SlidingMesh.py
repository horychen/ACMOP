#coding:utf-8
import femm, logging, os, sys, subprocess, winding_layout, utility
from collections import OrderedDict
from time import sleep
from time import time as clock_time
from pylab import plt, np

SELECT_ALL = 4
EPS = 1e-2 # unit mm

class Individual_Analyzer_FEMM_Edition(object):
    def __init__(self, p, ns, GroupSummary):

        self.p = p # number of pole pairs
        self.ns = ns # number of steps/points
        self.GroupSummary = GroupSummary

        self.femm_time = []
        self.femm_rotor_position_mech_deg = []

        self.femm_torque = []
        self.femm_forces = []
        self.femm_energy = []

        self.femm_circuitProperties = []
        self.femm_circuit_currents = []
        self.femm_circuit_voltages = []
        self.femm_circuit_fluxLinkages = []

        self.femm_motor_currents_d = []
        self.femm_motor_currents_q = []
        self.femm_bearing_currents_d = []
        self.femm_bearing_currents_q = []

        self.femm_motor_fluxLinkage_d = []
        self.femm_motor_fluxLinkage_q = []
        self.femm_bearing_fluxLinkage_d = []
        self.femm_bearing_fluxLinkage_q = []

        # get dict versions of the lists (to be used in collecting results parallel)
        self.list_of_attributes = []
        for attribute in dir(self):
            if 'femm' in attribute:
                self.list_of_attributes.append(attribute)
                exec(f'self.dict_{attribute} = dict()')

    def add(self, time, rotor_position_mech_deg, torque, forces, energy, circuitProperties, index=None):

        self.add_to_list_or_dict(time, self.femm_time, self.dict_femm_time, index=index)
        self.add_to_list_or_dict(rotor_position_mech_deg, self.femm_rotor_position_mech_deg, self.dict_femm_rotor_position_mech_deg, index=index)

        self.add_to_list_or_dict(torque, self.femm_torque, self.dict_femm_torque, index=index)
        self.add_to_list_or_dict(forces, self.femm_forces, self.dict_femm_forces, index=index)
        self.add_to_list_or_dict(energy, self.femm_energy, self.dict_femm_energy, index=index)
        # print('-'*100)
        # print(circuitProperties)
        # print(circuitProperties[0])
        # print(circuitProperties[0][0])
        # print(circuitProperties[0][0][0])
        circuit_currents     = list(map(lambda el: el[0], circuitProperties))
        circuit_voltages     = list(map(lambda el: el[1], circuitProperties))
        circuit_fluxLinkages = list(map(lambda el: el[2], circuitProperties))
        self.add_to_list_or_dict(circuit_currents, self.femm_circuit_currents, self.dict_femm_circuit_currents, index=index) # 6 coils UVW and AC/BD
        self.add_to_list_or_dict(circuit_voltages, self.femm_circuit_voltages, self.dict_femm_circuit_voltages, index=index) # 6 coils UVW and AC/BD
        self.add_to_list_or_dict(circuit_fluxLinkages, self.femm_circuit_fluxLinkages, self.dict_femm_circuit_fluxLinkages, index=index) # 6 coils UVW and AC/BD

        ''' Get UVW and dq quantities
        '''
        def get_motor_and_bearing_quantities_of_DPNV_winding(circuit_currents):
            circuit_current_coil_U_GrpA = circuit_currents[0]
            circuit_current_coil_U_GrpB = circuit_currents[3]
            circuit_current_coil_V_GrpA = circuit_currents[1]
            circuit_current_coil_V_GrpB = circuit_currents[4]
            circuit_current_coil_W_GrpA = circuit_currents[2]
            circuit_current_coil_W_GrpB = circuit_currents[5]
            motor_current_U   = np.array(circuit_current_coil_U_GrpA) + np.array(circuit_current_coil_U_GrpB)
            motor_current_V   = np.array(circuit_current_coil_V_GrpA) + np.array(circuit_current_coil_V_GrpB)
            motor_current_W   = np.array(circuit_current_coil_W_GrpA) + np.array(circuit_current_coil_W_GrpB)
            bearing_current_U = 0.5*(np.array(circuit_current_coil_U_GrpA) - np.array(circuit_current_coil_U_GrpB))
            bearing_current_V = 0.5*(np.array(circuit_current_coil_V_GrpA) - np.array(circuit_current_coil_V_GrpB))
            bearing_current_W = 0.5*(np.array(circuit_current_coil_W_GrpA) - np.array(circuit_current_coil_W_GrpB))
            return motor_current_U, motor_current_V, motor_current_W, bearing_current_U, bearing_current_V, bearing_current_W
        motor_current_U, motor_current_V, motor_current_W, bearing_current_U, bearing_current_V, bearing_current_W = get_motor_and_bearing_quantities_of_DPNV_winding(circuit_currents)
        motor_voltage_U, motor_voltage_V, motor_voltage_W, bearing_voltage_U, bearing_voltage_V, bearing_voltage_W = get_motor_and_bearing_quantities_of_DPNV_winding(circuit_voltages)
        motor_fluxLinkage_U, motor_fluxLinkage_V, motor_fluxLinkage_W, bearing_fluxLinkage_U, bearing_fluxLinkage_V, bearing_fluxLinkage_W = get_motor_and_bearing_quantities_of_DPNV_winding(circuit_fluxLinkages)

        # Amplitude invariant Clarke transform
        def amplitude_invariant_Clarke_transform(u,v,w):
            alpha = 2/3 * (    u - 0.5*v - 0.5*w)
            beta  = 2/3 * (np.sqrt(3)/2 * (v - w))
            gamma = 2/3 * (0.5*u + 0.5*v + 0.5*w)
            return alpha, beta, gamma

        # Currents
        motor_current_alpha, motor_current_beta, motor_current_gamma = amplitude_invariant_Clarke_transform(motor_current_U, motor_current_V, motor_current_W)
        bearing_current_alpha, bearing_current_beta, bearing_current_gamma = amplitude_invariant_Clarke_transform(bearing_current_U, bearing_current_V, bearing_current_W)
        # Fluxes
        motor_fluxLinkage_alpha, motor_fluxLinkage_beta, motor_fluxLinkage_gamma = amplitude_invariant_Clarke_transform(motor_fluxLinkage_U, motor_fluxLinkage_V, motor_fluxLinkage_W)
        bearing_fluxLinkage_alpha, bearing_fluxLinkage_beta, bearing_fluxLinkage_gamma = amplitude_invariant_Clarke_transform(bearing_fluxLinkage_U, bearing_fluxLinkage_V, bearing_fluxLinkage_W)

        def Park_transform(alpha, beta, theta):
            d =  np.cos(theta)*alpha + np.sin(theta)*beta
            q = -np.sin(theta)*alpha + np.cos(theta)*beta
            return d, q

        rotor_position_elec_rad = self.p * rotor_position_mech_deg/180*np.pi
        # Currents
        motor_current_d, motor_current_q     = Park_transform(motor_current_alpha, motor_current_beta, theta=rotor_position_elec_rad)
        bearing_current_d, bearing_current_q = Park_transform(bearing_current_alpha, bearing_current_beta, theta=rotor_position_elec_rad)
        self.add_to_list_or_dict(motor_current_d, self.femm_motor_currents_d, self.dict_femm_motor_currents_d, index=index)
        self.add_to_list_or_dict(motor_current_q, self.femm_motor_currents_q, self.dict_femm_motor_currents_q, index=index)
        self.add_to_list_or_dict(bearing_current_d, self.femm_bearing_currents_d, self.dict_femm_bearing_currents_d, index=index)
        self.add_to_list_or_dict(bearing_current_q, self.femm_bearing_currents_q, self.dict_femm_bearing_currents_q, index=index)
        # Fluxes
        motor_fluxLinkage_d, motor_fluxLinkage_q     = Park_transform(motor_fluxLinkage_alpha, motor_fluxLinkage_beta, theta=rotor_position_elec_rad)
        bearing_fluxLinkage_d, bearing_fluxLinkage_q = Park_transform(bearing_fluxLinkage_alpha, bearing_fluxLinkage_beta, theta=rotor_position_elec_rad)
        self.add_to_list_or_dict(motor_fluxLinkage_d, self.femm_motor_fluxLinkage_d, self.dict_femm_motor_fluxLinkage_d, index=index)
        self.add_to_list_or_dict(motor_fluxLinkage_q, self.femm_motor_fluxLinkage_q, self.dict_femm_motor_fluxLinkage_q, index=index)
        self.add_to_list_or_dict(bearing_fluxLinkage_d, self.femm_bearing_fluxLinkage_d, self.dict_femm_bearing_fluxLinkage_d, index=index)
        self.add_to_list_or_dict(bearing_fluxLinkage_q, self.femm_bearing_fluxLinkage_q, self.dict_femm_bearing_fluxLinkage_q, index=index)

    def add_to_list_or_dict(self, el, l, d, index=None):
        if index is None:
            l.append(el)
        else:
            d[index] = (el)

    def convert_dict_to_list(self):
        pass

    def get_ss_data(self):
        # Torque
        torque_average = self.torque_average = np.mean(self.femm_torque)
        torque_error = np.array(self.femm_torque) - torque_average
        ss_max_torque_error = (max(torque_error), min(torque_error))
        normalized_torque_ripple = (ss_max_torque_error[0] - ss_max_torque_error[1]) / torque_average

        # print('-'*100)
        # print(self.femm_circuit_currents)

        femm_currents_coil_U_GrpA = list(map(lambda el: el[0], self.femm_circuit_currents))
        femm_currents_coil_U_GrpB = list(map(lambda el: el[3], self.femm_circuit_currents))
        femm_currents_coil_V_GrpA = list(map(lambda el: el[1], self.femm_circuit_currents))
        femm_currents_coil_V_GrpB = list(map(lambda el: el[4], self.femm_circuit_currents))
        femm_currents_coil_W_GrpA = list(map(lambda el: el[2], self.femm_circuit_currents))
        femm_currents_coil_W_GrpB = list(map(lambda el: el[5], self.femm_circuit_currents))

        self.motor_current_U   = np.array(femm_currents_coil_U_GrpA) + np.array(femm_currents_coil_U_GrpB)
        self.motor_current_V   = np.array(femm_currents_coil_V_GrpA) + np.array(femm_currents_coil_V_GrpB)
        self.motor_current_W   = np.array(femm_currents_coil_W_GrpA) + np.array(femm_currents_coil_W_GrpB)
        self.bearing_current_U = np.array(femm_currents_coil_U_GrpA) - np.array(femm_currents_coil_U_GrpB)
        self.bearing_current_V = np.array(femm_currents_coil_V_GrpA) - np.array(femm_currents_coil_V_GrpB)
        self.bearing_current_W = np.array(femm_currents_coil_W_GrpA) - np.array(femm_currents_coil_W_GrpB)

        # Radial force
        force_x = list(map(lambda el: el[0], self.femm_forces))
        force_y = list(map(lambda el: el[1], self.femm_forces))
        self.sfv = sfv = utility.suspension_force_vector(force_x, force_y)
        print('[FEMM_SlidingMesh.py] SS data:')
        print('\t torque_average =', torque_average, 'Nm')
        print('\t ss_max_torque_error =', ss_max_torque_error , 'Nm')
        print('\t normalized_torque_ripple =', normalized_torque_ripple*100, '%')
        print('\t ss_avg_force_vector =', sfv.ss_avg_force_vector, 'N')
        print('\t ss_avg_force_magnitude =', sfv.ss_avg_force_magnitude, 'N')
        print('\t ss_avg_force_angle =', sfv.ss_avg_force_angle, 'deg')
        print('\t ss_max_force_err_abs =', sfv.ss_max_force_err_abs, 'N')
        print('\t ss_max_force_err_ang =', sfv.ss_max_force_err_ang, 'deg')

        # spec_performance_dict['Torque'] = None # this is dependent on stack length
        # spec_performance_dict['ForceAbs'] = None # this is dependent on stack length
        # spec_performance_dict['RotorWeight'] = None # this is dependent on stack length

        spec_performance_dict = dict()
        spec_performance_dict['cost_function'] = None
        spec_performance_dict['f1']  = None
        spec_performance_dict['f2']  = None
        spec_performance_dict['f3']  = None
        spec_performance_dict['FRW'] = None
        spec_performance_dict['normalized_torque_ripple'] = self.normalized_torque_ripple = normalized_torque_ripple
        spec_performance_dict['normalized_force_error_magnitude'] = self.normalized_force_error_magnitude = np.max(sfv.ss_max_force_err_abs) / sfv.ss_avg_force_magnitude
        spec_performance_dict['force_error_angle'] = self.force_error_angle = np.max(sfv.ss_max_force_err_ang)
        spec_performance_dict['project_name'] = None
        spec_performance_dict['individual_name'] = None
        spec_performance_dict['number_current_generation'] = None
        spec_performance_dict['individual_index'] = None
        spec_performance_dict['power_factor'] = None
        spec_performance_dict['rated_ratio'] = None
        spec_performance_dict['rated_stack_length_mm'] = None
        spec_performance_dict['rated_total_loss'] = None
        spec_performance_dict['rated_stator_copper_loss_along_stack'] = None
        spec_performance_dict['rated_rotor_copper_loss_along_stack'] = None
        spec_performance_dict['stator_copper_loss_in_end_turn'] = None
        spec_performance_dict['rotor_copper_loss_in_end_turn'] = None
        spec_performance_dict['rated_iron_loss'] = None
        spec_performance_dict['rated_windage_loss'] = None
        spec_performance_dict['str_results'] = None
        self.spec_performance_dict = spec_performance_dict

    def compute_objectives(self, acm_variant, select_fea_config_dict, toolFEA):
        # Power factor
        # power_factor = self.get_power_factor(ss_time_half_cycle=?, ss_voltage_half_cycle=?, ss_current_half_cycle=?, targetFreq=acm_variant.template.d['EX']['DriveW_Freq'])

        # Rated loss
        rated_total_loss, total_loss, \
                rated_stack_length_mm, rated_ratio, \
                required_torque, \
                Js, Jr, \
                Vol_Cu, \
                rated_stator_copper_loss_along_stack, \
                rated_rotor_copper_loss_along_stack, \
                stator_copper_loss_in_end_turn, \
                rotor_copper_loss_in_end_turn, \
                rated_iron_loss, \
                rated_windage_loss, \
                rated_magnet_Joule_loss = self.get_rated_loss(acm_variant, select_fea_config_dict, toolFEA)

        print(  f'''
                rated_total_loss, total_loss: {rated_total_loss}, {total_loss}, 
                rated_stack_length_mm, rated_ratio: {rated_stack_length_mm}, {rated_ratio}, 
                required_torque: {required_torque}, 
                Js, Jr: {Js}, {Jr}, 
                Vol_Cu: {Vol_Cu}, 
                rated_stator_copper_loss_along_stack: {rated_stator_copper_loss_along_stack}, 
                rated_rotor_copper_loss_along_stack: {rated_rotor_copper_loss_along_stack}, 
                stator_copper_loss_in_end_turn: {stator_copper_loss_in_end_turn}, 
                rotor_copper_loss_in_end_turn: {rotor_copper_loss_in_end_turn}, 
                rated_iron_loss: {rated_iron_loss}, 
                rated_windage_loss: {rated_windage_loss},
                rated_magnet_Joule_loss: {rated_magnet_Joule_loss}
                '''
            )

        # Mechnical constants        
        rotor_volume = acm_variant.template.get_rotor_volume() 
        rotor_weight = acm_variant.template.get_rotor_weight()
        shaft_power  = acm_variant.template.d['EX']['Omega'] * self.torque_average # make sure update_mechanical_parameters is called so that Omega corresponds to slip_freq_breakdown_torque
        efficiency   = shaft_power / (total_loss + shaft_power)  # 效率计算：机械功率/(损耗+机械功率)

        # THERMAL
        # if 'IM' in self.acm_variant.machine_type:
        #     stator_current_density = Js
        #     rotor_current_density  = Jr
        #     # print('Current density [Arms/m^2]:', stator_current_density, rotor_current_density, sep='\n')
        #     # if rotor_current_density > 8e6:
        #     #     print('rotor_current_density is over 8e6 Arms/m^2')
        # else:
        #                             # 基波电流幅值（在一根导体里的电流，六相逆变器中的GroupBDW相的电流，所以相当于已经考虑了并联支路数了）
        #     stator_current_density = dm.ui_info[2] / 1.4142135623730951 / (acm_variant.coils.mm2_slot_area*1e-6/acm_variant.template.d['EX']['DriveW_zQ'])
        #     print('[Loss-FEMMSlidingMesh.py] Data Magager: stator_current_density (GroupBDW) = %g Arms/m^2'%(stator_current_density))
        #     rotor_current_density = 0

        print('[Loss-FEMMSlidingMesh.py] Required torque: %g Nm'%(required_torque))
        print("[Loss-FEMMSlidingMesh.py] acm_variant.template.d['EX']['Omega']: %g rad/s"%(acm_variant.template.d['EX']['Omega']))
        rated_shaft_power  = acm_variant.template.d['EX']['Omega'] * required_torque
        rated_efficiency   = rated_shaft_power / (rated_total_loss + rated_shaft_power)  # 效率计算：机械功率/(损耗+机械功率)
        print('[FEMM_SlidingMesh.py] DEBUG efficiency =', efficiency)
        print('[FEMM_SlidingMesh.py] DEBUG rated_efficiency =', rated_efficiency)

        rated_rotor_volume = np.pi*(acm_variant.template.d['GP']['mm_r_or'].value*1e-3)**2 * (rated_stack_length_mm*1e-3)
        print('[Loss-FEMMSlidingMesh.py] rated_stack_length_mm, mm_template_stack_length:', rated_stack_length_mm, acm_variant.template.d['EX']['mm_template_stack_length'])

        # This weighted list suggests that peak-to-peak torque ripple of 5% is comparable with Em of 5% or Ea of 1 deg. Ref: Ye gu ECCE 2018
        # Eric suggests Ea is 1 deg. But I think this may be too much emphasis on Ea so large Trip does not matter anymore (not verified yet).
        list_weighted_ripples = [self.normalized_torque_ripple/0.05, self.sfv.normalized_force_error_magnitude/0.05, self.sfv.force_error_angle]

        if 'IM' in acm_variant.template.machine_type:
            # - Torque per Rotor Volume
            f1_IM = - required_torque / rated_rotor_volume
            f1 = f1_IM
        elif 'PMSM' in acm_variant.template.machine_type:
            # - Cost
            price_per_volume_steel    = 0.28  * 61023.744 # $/in^3 (M19 Gauge26) # 0.23 for low carbon, semi-processed 24 Gauge electrical steel
            price_per_volume_copper   = 1.2   * 61023.744 # $/in^3 wire or bar or end-ring
            price_per_volume_magnet   = 11.61 * 61023.744 # $/in^3 NdFeB PM
            # price_per_volume_aluminum = 0.88  / 16387.064 # $/in^3 wire or cast Al

            # Vol_Fe = (2*acm_variant.template.d['GP']['mm_r_os'].value*1e-3) ** 2 * (rated_stack_length_mm*1e-3) # 注意，硅钢片切掉的方形部分全部消耗了。# Option 1 (Jiahao)
            Vol_Fe = ( np.pi*(acm_variant.template.d['GP']['mm_r_os'].value*1e-3)**2 - np.pi*(acm_variant.template.d['GP']['mm_r_ri'].value*1e-3)**2 ) * (rated_stack_length_mm*1e-3) # Option 2 (Eric)

            Vol_PM = (acm_variant.rotorMagnet.mm2_magnet_area*1e-6) * (rated_stack_length_mm*1e-3)

            print('[Loss-FEMMSlidingMesh.py] Gross-Area_Fe', (acm_variant.template.d['GP']['mm_r_os'].value*1e-3) ** 2        *1e6, 'mm^2')
            print('[Loss-FEMMSlidingMesh.py] Gross-Area_Cu (est.)',   Vol_Cu/(rated_stack_length_mm*1e-3)                     *1e6, 'mm^2')
            print('[Loss-FEMMSlidingMesh.py] Gross-Area_PM', (acm_variant.rotorMagnet.mm2_magnet_area*1e-6)                   *1e6, 'mm^2')
            print('[Loss-FEMMSlidingMesh.py] Gross-Volume_Fe',    Vol_Fe, 'm^3')
            print('[Loss-FEMMSlidingMesh.py] Gross-Volume_Cu',    Vol_Cu, 'm^3')
            print('[Loss-FEMMSlidingMesh.py] Gross-Volume_PM',    Vol_PM, 'm^3')
            f1_PMSM =    Vol_Fe * price_per_volume_steel \
                    +    Vol_Cu * price_per_volume_copper\
                    +    Vol_PM * price_per_volume_magnet
            print('[Loss-FEMMSlidingMesh.py] Cost Fe: $', Vol_Fe * price_per_volume_steel )
            print('[Loss-FEMMSlidingMesh.py] Cost Cu: $', Vol_Cu * price_per_volume_copper)
            print('[Loss-FEMMSlidingMesh.py] Cost PM: $', Vol_PM * price_per_volume_magnet)
            f1 = f1_PMSM

        # - Efficiency @ Rated Power
        f2 = - rated_efficiency
        # Ripple Performance (Weighted Sum)
        f3 = sum(list_weighted_ripples)

        FRW = self.sfv.ss_avg_force_magnitude / rotor_weight
        print('[Loss-FEMMSlidingMesh.py] FRW:', FRW, ', Rotor weight:', rotor_weight, '[N], Stack length:', acm_variant.template.d['EX']['mm_template_stack_length'], 'mm, Rated stack length:', rated_stack_length_mm, 'mm')
        rated_rotor_volume = acm_variant.template.get_rotor_volume(stack_length=rated_stack_length_mm) 
        rated_rotor_weight = acm_variant.template.get_rotor_weight(stack_length=rated_stack_length_mm)
        print('[Loss-FEMMSlidingMesh.py] rated_rotor_volume:', rated_rotor_volume, 'rated_rotor_weight:', rated_rotor_weight)

        rated_results = [   rated_shaft_power, 
                            rated_efficiency,
                            rated_total_loss, 
                            rated_stator_copper_loss_along_stack, 
                            rated_rotor_copper_loss_along_stack, 
                            stator_copper_loss_in_end_turn, 
                            rotor_copper_loss_in_end_turn, 
                            rated_iron_loss, 
                            rated_windage_loss,
                            rated_rotor_volume,
                            rated_stack_length_mm,  # new!
                            acm_variant.template.d['EX']['mm_template_stack_length']]           # new! 在计算FRW的时候，我们只知道原来的叠长下的力，所以需要知道原来的叠长是多少。

        # machine_results = [power_factor, efficiency, self., normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]


    def get_rated_loss(self, acm_variant, select_fea_config_dict, toolFEA):

        ''' loss results with axial length: mm_template_stack_length
        '''
        # Copper loss
        stator_copper_loss, rotor_copper_loss, stator_copper_loss_along_stack, rotor_copper_loss_along_stack, Js, Jr, Vol_Cu = self.get_copper_loss(acm_variant, femm_solver=None)
        stator_copper_loss_in_end_turn = stator_copper_loss - stator_copper_loss_along_stack 
        rotor_copper_loss_in_end_turn  = rotor_copper_loss - rotor_copper_loss_along_stack
        # iron, coil proximity, magnet Joule losses
        iron_loss, prox_loss, magnet_Joule_loss = self.get_iron_loss_and_prox_loss_and_magnet_loss(acm_variant, select_fea_config_dict, toolFEA)
        windage_loss = utility.get_windage_loss(acm_variant, acm_variant.template.d['EX']['mm_template_stack_length'])
        total_loss   = stator_copper_loss_along_stack + magnet_Joule_loss + iron_loss + windage_loss + prox_loss

        ''' loss results with revised axial length: rated_stack_length_mm
        '''
        speed_rpm       = acm_variant.template.SI['ExcitationFreqSimulated'] * 60 / acm_variant.template.SI['p'] # rpm
        required_torque = acm_variant.template.SI['mec_power'] / (2*np.pi*speed_rpm)*60

        rated_ratio                          = required_torque / self.torque_average 
        rated_stack_length_mm                = rated_ratio * acm_variant.template.d['EX']['mm_template_stack_length']
        rated_stator_copper_loss_along_stack = rated_ratio * stator_copper_loss_along_stack
        rated_rotor_copper_loss_along_stack  = rated_ratio * rotor_copper_loss_along_stack
        rated_iron_loss                      = rated_ratio * iron_loss
        rated_windage_loss                   = utility.get_windage_loss(acm_variant, rated_stack_length_mm)
        rated_prox_loss                      = rated_ratio * prox_loss
        rated_magnet_Joule_loss              = rated_ratio * magnet_Joule_loss

        # total_loss   = copper_loss + iron_loss + windage_loss
        rated_total_loss = rated_stator_copper_loss_along_stack \
                        + rated_rotor_copper_loss_along_stack \
                        + stator_copper_loss_in_end_turn \
                        + rotor_copper_loss_in_end_turn \
                        + rated_iron_loss \
                        + rated_windage_loss \
                        + rated_prox_loss \
                        + rated_magnet_Joule_loss
        return  rated_total_loss, total_loss, \
                rated_stack_length_mm, rated_ratio, \
                required_torque, \
                Js, Jr, \
                Vol_Cu, \
                            rated_stator_copper_loss_along_stack, \
                            rated_rotor_copper_loss_along_stack, \
                            stator_copper_loss_in_end_turn, \
                            rotor_copper_loss_in_end_turn, \
                            rated_iron_loss, \
                            rated_windage_loss, \
                            rated_magnet_Joule_loss

    def get_copper_loss(self, acm_variant, femm_solver=None):
        if femm_solver is not None:
            try:
                # convert rotor current results (complex number) into its amplitude
                femm_solver.list_rotor_current_amp = [abs(el) for el in femm_solver.vals_results_rotor_current] # el is complex number

                _s, _r, _sAlongStack, _rAlongStack, _Js, _Jr = femm_solver.get_copper_loss_pyrhonen(acm_variant.slot_area_utilizing_ratio*femm_solver.stator_slot_area, 
                                                                                                                            femm_solver.rotor_slot_area, 
                                                                                                                            total_CurrentAmp=acm_variant.DriveW_CurrentAmp+acm_variant.BeariW_CurrentAmp)
                s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = femm_solver.get_copper_loss_Bolognani(acm_variant.slot_area_utilizing_ratio*femm_solver.stator_slot_area, 
                                                                                                                            femm_solver.rotor_slot_area, 
                                                                                                                            total_CurrentAmp=acm_variant.DriveW_CurrentAmp+acm_variant.BeariW_CurrentAmp,
                                                                                                                            STATOR_SLOT_FILL_FACTOR=acm_variant.template.SI['WindingFill'],
                                                                                                                            TEMPERATURE_OF_COIL=acm_variant.template.SI['Temperature'])

                msg1 = 'Pyrhonen : %g, %g | %g, %g | %g, %g ' % (_s, _r, _sAlongStack, _rAlongStack, _Js, _Jr) 
                msg2 = 'Bolognani: %g, %g | %g, %g | %g, %g ' % (s, r, sAlongStack, rAlongStack, Js, Jr) 
                print('[FEMM_SlidingMesh.py] copper loss', msg1)
                print('[FEMM_SlidingMesh.py] copper loss', msg2)
                logger = logging.getLogger(__name__)
                logger.debug(msg1)
                logger.debug(msg2)
            except Exception as e:
                raise e
        else:
            SI = acm_variant.template.SI
            GP = acm_variant.template.d['GP']
            EX = acm_variant.template.d['EX']
            wily = EX['wily']
            copper_loss_parameters = [GP['mm_d_sleeve'].value + GP['mm_d_fixed_air_gap'].value,
                                GP['mm_w_st'].value,
                                wily.number_parallel_branch,
                                EX['DriveW_zQ'],
                                wily.coil_pitch_y,
                                acm_variant.template.SI['Qs'],
                                EX['mm_template_stack_length'],
                                EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp'], # total current amplitude
                                GP['mm_r_or'].value,       # mm
                                GP['mm_r_os'].value*2*1e-3 # m, stator_yoke_diameter_Dsyi
                                ]
            s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = utility.get_copper_loss_Bolognani(EX['slot_current_utilizing_ratio']*acm_variant.coils.mm2_slot_area*1e-6, # slot_current_utilizing_ratio is 1 for combined winding and is less than 1 for separate winding
                                                                                                copper_loss_parameters=copper_loss_parameters,
                                                                                                STATOR_SLOT_FILL_FACTOR=acm_variant.template.SI['WindingFill'],
                                                                                                TEMPERATURE_OF_COIL=acm_variant.template.SI['Temperature'])

        return s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu

    def get_iron_loss_and_prox_loss_and_magnet_loss(self, acm_variant, select_fea_config_dict, toolFEA):
        if 'JMAG' in select_fea_config_dict:
            iron_loss = toolFEA.dm.jmag_loss_list[2]
            prox_loss = 0.0
            if 'IM' in acm_variant.template.machine_type:
                rotor_copper_loss = toolFEA.dm.jmag_loss_list[1]
            elif 'PM' in acm_variant.template.machine_type:
                magnet_Joule_loss = toolFEA.dm.jmag_loss_list[1]
        elif 'FEMM' in select_fea_config_dict:
            _, r, s, magnet_Joule_loss, _, prox_loss, _ = acm_variant.analyzer.iron_loss_calculation_by_dft(acm_variant)
            iron_loss = r+s
        return iron_loss, prox_loss, magnet_Joule_loss

    def get_power_factor(self, ss_time_half_cycle, ss_voltage_half_cycle, ss_current_half_cycle, targetFreq=1e3, numPeriodicalExtension=1000):
        # number_of_steps_at_steady_state: steps corresponding to half the period 
        # mytime  = ss_time_half_cycle
        # voltage = ss_voltage_half_cycle
        # current = ss_current_half_cycle

        # from pylab import *
        # print len(mytime), len(voltage), len(current)
        # figure()
        # plot(mytime, voltage)
        # plot(mytime, current)
        # show()
        power_factor, u, i, phase_diff_ui = utility.compute_power_factor_from_half_period(ss_voltage_half_cycle, ss_current_half_cycle, ss_time_half_cycle, targetFreq=targetFreq, numPeriodicalExtension=numPeriodicalExtension)
        self.ui_info = [power_factor, u, i, phase_diff_ui]
        print('[FEMM_SlidingMesh.py] ui_info', self.ui_info)
        return power_factor

    def iron_loss_preparing_data_matrices(self, index):
        if index==0:
            # Record the initial mesh elements if the first time through the loop
            self.nn = int(femm.mo_numelements())

            self.z  = np.zeros([self.nn,       1], dtype=np.complex64) # Location of the centroid of each element
            self.a  = np.zeros([self.nn,       1]) # Area of each element
            self.g  = np.zeros([self.nn,       1]) # Block label of each element
            for m in range(self.nn): # start from 0 for indexing but from 1 for counting 
                counting_element = m+1
                elm = femm.mo_getelement(counting_element)
                # z is a vector of complex numbers that represents the location of the centroid of each element.
                self.z[m] = elm[4-1] + 1j*elm[5-1]
                # element area in the length units used to draw the geometry
                self.a[m] = elm[6-1]
                # group number associated with the element
                self.g[m] = elm[7-1]

            # mo_getprobleminfo returns, among other things, the depth of the
            # machine in the into-the-page direction and the length units used to
            # draw the geometry. Both of these pieces of information will be needed
            # to integrate the losses over the volume of the machine.
            self.probinfo = femm.mo_getprobleminfo()

            self.A  = np.zeros([self.ns, self.nn]) # matrix that will hold the vector potential info
            self.b  = np.zeros([self.ns, self.nn], dtype=np.complex64) # matrix that will hold the flux density info

        # Store element flux densities B and magnetic potential A
        for m in range(self.nn):
            if self.g[m]==self.GroupSummary['magnet']: # Element is in a rotor magnet, marked with group numbers 11 and higher
                # Store vector potential at the element centroid for elements that are in PMs
                self.A[index,m] = femm.mo_geta( float(np.real(self.z[m])), 
                                                float(np.imag(self.z[m])) )
            elif self.g[m] == self.GroupSummary['stator_iron_core'] \
              or self.g[m] == self.GroupSummary['rotor_iron_core'] \
              or self.g[m] == self.GroupSummary['coils']: # Element is on the stator or rotor iron or coils
                # Store flux density at the element centroid for these elements
                b_temp = femm.mo_getb(  float(np.real(self.z[m])), 
                                        float(np.imag(self.z[m])) )
                self.b[index,m] = b_temp[0] + 1j*b_temp[1]

    def iron_loss_calculation_by_dft(self, acm_variant):
        def matlab_fft(matrix_input):
            return np.fft.fft(matrix_input.T).T

        # Compute the square of the amplitude of each harmonic at the centroid of
        # each element in the mesh. Matlab's built-in FFT function makes this easy.
        bxfft = np.abs(matlab_fft(np.real(self.b))) * (2/self.ns)
        byfft = np.abs(matlab_fft(np.imag(self.b))) * (2/self.ns)
        bsq  = (bxfft*bxfft) + (byfft*byfft)

        # Compute the volume of each element in units of meter**3
        stack_length_m = self.probinfo[3-1] # Length of the machine in the into-the-page direction
        lengthunits    = self.probinfo[4-1] # Length of drawing unit in meters
        self.v         = self.a*stack_length_m*lengthunits**2
        print('LOSS-DEBUG, stack_length_m, lengthunits:', stack_length_m, lengthunits)

        # compute fft of A at the center of each element
        Afft = matlab_fft(self.A)*(2/self.ns)
        # filter elements by the group number of magnet
        g3 = (self.g==self.GroupSummary['magnet'])
        # total volume of the magnet under consideration
        vmag = np.dot(self.v.T, g3) # [inner product]
        # average current in the magnet for each harmonic
        Jo = np.dot(Afft, self.v * g3) / vmag
        # subtract averages off of each element in the magnet
        Jm = Afft - np.dot(Jo, g3.T)

        if True:
            thisFrequency = acm_variant.template.d['EX']['Omega']/(2*np.pi) # thisSpeed/60 # mechanical speed in Hz

            # Make a vector representing the frequency associated with each harmonic
            # The last half of the entries are zeroed out so that we don't count each
            # harmonic twice--the upper half of the FFT a mirror of the lower half
            w_ = np.arange(self.ns)
            MyLowestHarmonic = self.p
            if acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS'] == 0.5:
                # the flux density can be extended to get the other half for calcuting iron losses 
                raise Exception('Iron loss cannot be calcultated because the DFT resolution is higher than flux density waveform fundamental frequency.')
            elif acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS'] < 1:
                raise Exception('Iron loss cannot be calcultated because the DFT resolution is higher than flux density waveform fundamental frequency.')
            elif acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS'] > 1:
                print('[Warning] More than 1 cycle is not beneficial for iron loss calculation as there is no anti-aliasing filter applied to flux density waveform and if high frequency component is present, then aliasing will ruins your results if you have many bins between 0 Hz and fundamental frequency.]')

            DFT_RESOLUTION_1 = MyLowestHarmonic*thisFrequency / acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS']
            w = DFT_RESOLUTION_1*w_ * (w_<(self.ns/2))
            DFT_RESOLUTION_2 = MyLowestHarmonic*thisFrequency # this is correct only when acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS'] = 1
            print('Freuqencies/Harmonics:', w, DFT_RESOLUTION_1, DFT_RESOLUTION_2)

            # iron loss (Dividing the result by cs corrects for the lamination stacking factor)
            if 'M19' in acm_variant.template.SI['Steel'] or 'M15' in acm_variant.template.SI['Steel']:
                # M-19 Steel
                ce = 0.530 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
                ch = 143.  # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
                cs = 0.95  # % Lamination stacking factor (nondimensional)
            elif acm_variant.template.SI['Steel'] == 'Arnon5':
                # %Arnon7
                ce = 0.07324 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
                ch = 187.6   # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
                cs = 0.96    # % Lamination stacking factor (nondimensional)
            g1=(self.g==self.GroupSummary['rotor_iron_core'])
            g2=(self.g==self.GroupSummary['stator_iron_core'])
            rotor_loss  = np.dot(np.dot((ch*w+ce*w*w), bsq), (self.v*g1)) / cs
            stator_loss = np.dot(np.dot((ch*w+ce*w*w), bsq), (self.v*g2)) / cs

            fig, axes = plt.subplots(3)
            x, y1, y2, y3 = [], [], [], []
            for jj in range(int(self.ns)):
                # _w = (w==jj*MyLowestHarmonic*thisFrequency)
                # if sum(_w) == 0: raise
                _w = (w==w[jj]) * w[jj]
                _rotor_loss  = np.dot(np.dot((ch*_w+ce*_w*_w), bsq), (self.v*g1)) / cs
                _stator_loss = np.dot(np.dot((ch*_w+ce*_w*_w), bsq), (self.v*g2)) / cs
                # plt.scatter(jj*MyLowestHarmonic*thisFrequency, _stator_loss)
                x.append(w[jj])
                y1.append(_rotor_loss)
                y2.append(_stator_loss)
            ax=axes[0]; ax.stem(x, y1); ax.set_ylabel('Rotor loss [W]')
            ax=axes[1]; ax.stem(x, y2); ax.set_ylabel('Stator loss [W]')

            # and prox losses can be totalled up in a similar way
            g4 = (self.g==self.GroupSummary['coils'])
            WindingFill = acm_variant.template.SI['WindingFill']
            TemperatureRise = acm_variant.template.SI['Temperature'] # temperature increase, degrees C
            AWG   = acm_variant.template.SI['AWG'] # Magnet wire gauge used in winding
            if AWG != 25: raise Exception('Not implemented for AWG other than 25.')
            dwire = 0.324861*0.0254*np.exp(-0.115942*AWG) # wire diameter in meters as a function of AWG
            owire = (58*10**6)/(1+TemperatureRise*0.004) # conductivity of the wire in S/m at prescribed deltaT
            cePhase = (np.pi**2/8)*dwire**2*WindingFill*owire
            prox_loss = np.dot(np.dot((cePhase*w*w), bsq), (self.v*g4))

            # copper loss (cannot be calculated becasue PhaseResistance is unknown)
                # PhaseResistance = 0.223*2.333 # phase resistance including end turns at 20 degC
                # Iphase  = np.sqrt(MyIdCurrent**2+MyIqCurrent**2)/np.sqrt(2)
                # PhaseOhmic = 3*(PhaseResistance*(1+TemperatureRise*0.004))*Iphase**2
            PhaseOhmic = 0.0

            # Add up eddy current losses in the magnets
            omag = 0.556*10**6 # conductivity of sintered NdFeB in S/m
            magnet_loss = 0.5 * np.dot((omag*(2*np.pi*w)**2), np.dot((np.abs(Jm)**2), self.v)) # g3 is already taken into account in variable Jm
                        # 0.5 = (1/sqrt(2)) ** 2, the peak value Jm should be converted to RMS value for calculating loss power


            for jj in range(int(self.ns)):
                _w = (w==w[jj]) * w[jj]
                _magnet_loss = 0.5 * np.dot((omag*(2*np.pi*_w)**2), np.dot((np.abs(Jm)**2), self.v))
                y3.append(_magnet_loss)
            ax=axes[2]; ax.stem(x, y3); ax.set_ylabel('Magnet loss [W]')
            _f = (w<w[8])
            _rotor_loss  = np.dot(np.dot((ch*_f*w+ce*_f*w*w), bsq), (self.v*g1)) / cs
            _stator_loss = np.dot(np.dot((ch*_f*w+ce*_f*w*w), bsq), (self.v*g2)) / cs
            _magnet_loss = 0.5 * np.dot((omag*(2*np.pi*_f*w)**2), np.dot((np.abs(Jm)**2), self.v))
            print(f'[FEMM_SlidingMesh.py] Conservative iron loss: {_rotor_loss}< {rotor_loss}; {_stator_loss} < {stator_loss}')
            print(f'[FEMM_SlidingMesh.py] Conservative magnet loss: {_magnet_loss}< {magnet_loss}')



            total_loss = rotor_loss + stator_loss + prox_loss + PhaseOhmic + magnet_loss
            print('LOSS-DEBUG: rotor_loss, stator_loss, prox_loss, PhaseOhmic, magnet_loss')
            print('LOSS-DEBUG:', rotor_loss, stator_loss, prox_loss, PhaseOhmic, magnet_loss)

            vsiron = np.dot(self.v.T, g1)
            vriron = np.dot(self.v.T, g2)
            vcoil = WindingFill * np.dot(self.v.T, g4)
            print('[FEMM_SlidingMesh.py] DEBUG vmag, vsiron, vriron, vcoil:', vmag, vsiron, vriron, vcoil, 'mm2')

            return [thisFrequency*60, rotor_loss, stator_loss, magnet_loss, PhaseOhmic, prox_loss, total_loss]

class FEMM_SlidingMesh(object):

    def __init__(self, acm_variant):
        self.acm_variant = acm_variant
        self.rotor_phase_name_list = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        # local references
        GP = self.acm_variant.template.d['GP']
        EX = self.acm_variant.template.d['EX']
        p  = self.acm_variant.template.SI['p']

        # figure out initial rotor angle
        deg_pole_span = 180/p
        self.initial_rotor_position_mech_deg = (deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + EX['wily'].deg_winding_U_phase_phase_axis_angle + deg_pole_span
        print('[FEMM_SlidingMesh.py] initial_rotor_position_mech_deg:', self.initial_rotor_position_mech_deg, 'deg || deg_winding_U_phase_phase_axis_angle:', '', EX['wily'].deg_winding_U_phase_phase_axis_angle)

        # figure out step size
        self.electrical_period = acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS']/EX['DriveW_Freq']
        self.number_of_steps   = acm_variant.template.fea_config_dict['femm.number_of_steps_2ndTSS']
        self.step_size_sec = self.electrical_period / (self.number_of_steps-1)  # <--- Explanation of -1 here:
                                                                                # In order to have full DFT resolution, we need to actually run one more point than the specified {fea_config_dict['femm.number_of_steps_2ndTSS']=}")
                                                                                # Or, we can make the step_size_sec larger such that self.step_size_sec * self.number_of_steps > one electrical period (i.e., self.step_size_sec * (self.number_of_steps-1) === one electrical period)
        self.step_size_mech_deg = EX['Omega'] * self.step_size_sec / np.pi * 180

    def open(self):
        femm.openfemm(True) # bHide # False for debug
        femm.newdocument(0) # magnetic

    def probdef(self, stack_length=100, time_harmonic_study_frequency=0):
        self.time_harmonic_study_frequency = time_harmonic_study_frequency
        # femm.smartmesh(False) <- This will not work due to bug of femm.__init__ 
        femm.callfemm_noeval('smartmesh(0)') # call this before probdef?
        femm.mi_probdef(time_harmonic_study_frequency, 'millimeters', 'planar', 1e-8, # must < 1e-8
                        stack_length, 18, 1) # The acsolver parameter (default: 0) specifies which solver is to be used for AC problems: 0 for successive approximation, 1 for Newton.
                                             # 1 for 'I intend to try the acsolver of Newton, as this is the default for JMAG@[Nonlinear Calculation] Setting Panel in the [Study Properties] Dialog Box'
        # femm.callfemm_noeval('mi_smartmesh(0)') # call this after probdef?
        self.bool_automesh = False # setting to false gives no effect?

        femm.smartmesh(-1) # let mi_smartmesh deside. You must turn it off in parasolver.py
        # self.bool_automesh = True

    def pre_process(self, project_file_name):
        # if self.deg_per_step == 0.0:
        #     print('Locked Rotor! Run 40 stEPS for one slip period.')
        #     self.im.update_mechanical_parameters(syn_freq=0.0)
        self.add_material()
        self.add_block_labels_static_solver()
        self.add_sliding_band_and_boundary(deg_innerAngle=0) # modifying parameter here has no effect
        self.update_circuit_excitation(time=0) # modifying parameter here has no effect
        femm.mi_saveas(project_file_name)
        print('[FEMM_SlidingMesh.py] Save as:', project_file_name)
        self.project_file_name = project_file_name

    def add_material(self):
        # mi_addmaterial('matname', mu x, mu y, H c, J, Cduct, Lam d, Phi hmax, lam fill, LamType, Phi hx, Phi hy, nstr, dwire)
        femm.mi_getmaterial('Air')
        femm.mi_getmaterial('Copper') # for coil
        # femm.mi_getmaterial('18 AWG') # for coil
        # femm.mi_getmaterial('Aluminum, 1100') # for bar?
        # femm.mi_getmaterial('304 Stainless Steel') # for shaft?

        # femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0);
        # femm.mi_addmaterial('Aluminum', 1, 1, 0, 0, 35, 0, 0, 1, 0, 0, 0)
        if 'Bar_Conductivity' in self.acm_variant.template.d['EX'].keys():
            # spec_derive_dict?
            femm.mi_addmaterial('Aluminum', 1, 1, 0, 0, self.acm_variant.template.d['EX']['Bar_Conductivity']*1e-6, 0, 0, 1, 0, 0, 0) # [MS/m]
        # femm.mi_addmaterial('Aluminum', 1, 1, 0, 0, 1/1.673e-2, 0, 0, 1, 0, 0, 0)

        # femm.mi_addmaterial('LinearIron', 2000, 2000, 0, 0, 0, 0, 0, 1, 0, 0, 0);

        if self.acm_variant.template.spec_input_dict['Steel'] == 'M19Gauge29':
            # femm.mi_getmaterial('M-19 Steel') # for Stator & Rotor Iron Cores (Nonlinear with B-H curve)
            femm.mi_addmaterial('M19Gauge29',0,0, 0,0, 0,0.3556,0, 0.95) # no lamination for testing consistency with JMAG
            hdata, bdata = np.loadtxt(self.acm_variant.template.fea_config_dict['dir.parent'] + 'BH/M-19-Steel-BH-Curve-afterJMAGsmooth.BH', unpack=True, usecols=(0,1))
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('M19Gauge29', bdata[n], hdata[n])

        elif self.acm_variant.template.spec_input_dict['Steel'] == 'Arnon5':
            # Arnon5 is 1/5 thick as M15, which is too thin to use and it is expensive as well
            femm.mi_addmaterial('Arnon5-final',0,0, 0,0, 0.0,0.127,0, 0.96)
            # BH = np.loadtxt(self.dir_parent + 'BH/Arnon5/Arnon5-final.txt', unpack=True, usecols=(0,1))
            BH = np.loadtxt('D:/DrH/Codes/c/Arnon5/Arnon5-final.txt', unpack=True, usecols=(0,1))

            bdata = BH[1][1:] # if not skip the first point, there will be two (0,0) in FEMM software, reason unknown.
            hdata = BH[0][1:] # if not skip the first point, there will be two (0,0) in FEMM software, reason unknown.
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('Arnon5-final', bdata[n], hdata[n])

        elif self.acm_variant.template.spec_input_dict['Steel'] == 'M15':
            femm.mi_addmaterial('My M-15 Steel',0,0, 0,0, 0,0.635,0, 0.98)
            BH = np.loadtxt(self.dir_codes + '../Arnon5/M-15-Steel-BH-Curve.txt', unpack=True, usecols=(0,1))
            bdata = BH[1]
            hdata = BH[0]
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('My M-15 Steel', bdata[n], hdata[n])


        if False:
            # A more interesting material to add is the iron with a nonlinear
            # BH curve.  First, we create a material in the same way as if we
            # were creating a linear material, except the values used for
            # permeability are merely placeholders.
            femm.mi_addmaterial('Arnon5',0,0,0,0, 0.0,0.127,0,0.96)
            # A set of points defining the BH curve is then specified.
            BH = np.loadtxt(self.dir_codes + 'Arnon5_Kang_after_JMAG_Smoothed.txt', unpack=True, usecols=(0,1))
            bdata = BH[1]
            hdata = BH[0]
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('Arnon5', bdata[n], hdata[n])

    @staticmethod
    def set_block_label(group_no, material_name, p, meshsize_if_no_automesh, incircuit='<None>', turns=0, automesh=True, magdir=0):
        femm.mi_addblocklabel(p[0],p[1])
        femm.mi_selectlabel(p[0],p[1])
        femm.mi_setblockprop(material_name, automesh, meshsize_if_no_automesh, incircuit, magdir, group_no, turns)
        femm.mi_clearselected()

    def add_block_labels_static_solver(self, fraction=1):
        # add_block_labels_static_solver is implemented with the new general DPNV winding implementation

        self.GroupSummary = {
            "air_gap"          : 9,
            "stator_iron_core" : 10,
            "coils"            : 11,
            "rotor_iron_core"  : 100,
            "magnet"           : 101,
            "shaft"            : 102,
        }

        SERIES_CONNECTED = 1
        PARALLEL_CONNECTED = 0

        if self.acm_variant.template.fea_config_dict['femm.Coarse_Mesh']==True: # Coarse mesh
            if self.acm_variant.template.fea_config_dict['femm.Coarse_Mesh_Level'] == 2:
                MESH_SIZE_ALUMINUM = 2 * 6    # 3
                MESH_SIZE_STEEL    = 2 * 6    # 4
                MESH_SIZE_AIR      = 2 * 0.75 # 0.5 
                MESH_SIZE_COPPER   = 2 * 10   # 8
            elif self.acm_variant.template.fea_config_dict['femm.Coarse_Mesh_Level'] == 3:
                MESH_SIZE_ALUMINUM = 1 * 6    # 3
                MESH_SIZE_STEEL    = 1 * 6    # 4
                MESH_SIZE_AIR      = 1 * 0.75 # 0.5 
                MESH_SIZE_COPPER   = 1 * 10   # 8
            else:
                raise Exception('Invalid femm.Coarse_Mesh_Level.')
        else:
            MESH_SIZE_ALUMINUM = 3
            MESH_SIZE_STEEL    = 4
            MESH_SIZE_AIR      = 0.5 
            MESH_SIZE_COPPER   = 8

        GP = self.acm_variant.template.d['GP']

        # Air region (Stator)
        X, Y = -(GP['mm_r_or'].value+5*EPS), 0
        self.set_block_label(9, 'Air', (X, Y), MESH_SIZE_AIR, automesh=self.bool_automesh)
        # Air region (Rotor)
        X, Y =   GP['mm_r_si'].value-5*EPS, 0
        self.set_block_label(9, 'Air', (X, Y), MESH_SIZE_AIR, automesh=self.bool_automesh)

        # # Air Gap Boundary for Rotor Motion #2
        # self.set_block_label(9, '<No Mesh>',   (0, im.Radius_OuterRotor+0.5*im.Length_AirGap), 5, automesh=self.bool_automesh)
        # self.set_block_label(9, 'Air',         (0, im.Radius_OuterRotor+0.7*im.Length_AirGap), 0.5, automesh=self.bool_automesh)
        # self.set_block_label(9, 'Air',         (0, im.Radius_OuterRotor+0.3*im.Length_AirGap), 0.5, automesh=self.bool_automesh)


        # shaft
        if fraction == 1:
            self.set_block_label(102, '<No Mesh>',         (0, 0),  20)
            # self.set_block_label(101, 'Air',         (0, 0),  10, automesh=True) # for deeply-saturated rotor yoke

        # Iron Core @225deg
        if 'M19' in self.acm_variant.template.spec_input_dict['Steel']:
            steel_name = 'M19Gauge29'
        elif self.acm_variant.template.spec_input_dict['Steel'] == 'M15':
            steel_name = 'My M-15 Steel' 
        elif self.acm_variant.template.spec_input_dict['Steel'] == 'Arnon5':
            steel_name = 'Arnon5-final'
        X, Y = -(GP['mm_r_ri'].value + 0.5*GP['mm_d_ri'].value), 0
        self.set_block_label(100, steel_name, (X, Y), MESH_SIZE_STEEL, automesh=self.bool_automesh)
        X, Y = (0.5*(GP['mm_r_si'].value + GP['mm_r_os'].value)), 0
        self.set_block_label(10, steel_name, (X, Y), MESH_SIZE_STEEL, automesh=self.bool_automesh)

        # Rotor Magnet
        THETA = GP['deg_alpha_rm'].value * 0.5 / 180 * np.pi
        R = GP['mm_r_or'].value - 0.5*GP['mm_d_rp'].value
        femm.mi_getmaterial('N40')
        for _ in range(2*self.acm_variant.template.SI['p']):
            X, Y = R*np.cos(THETA), R*np.sin(THETA)
            self.set_block_label(101, 'N40', (X, Y), MESH_SIZE_STEEL, automesh=self.bool_automesh, magdir=THETA/np.pi*180 + _%2*180)
            deg_alpha_rp = 360 / (2*self.acm_variant.template.SI['p'])
            THETA += deg_alpha_rp / 180 * np.pi

        # Circuit Configuration
        # Rotor Winding
        pass

        # Stator Winding
        EX = self.acm_variant.template.d['EX']
        wily = EX['wily']
        npb = wily.number_parallel_branch
        nwl = wily.number_winding_layer # number of windign layers 
        ampD = EX['DriveW_CurrentAmp']/npb
        ampB = EX['BeariW_CurrentAmp']
        turnsDPNV = EX['DriveW_zQ']/nwl #* 100 # DEBUG
        # turnsB = EX['BeariW_zQ']/nwl #* 100 # DEBUG
        print('[FEMM_SllidingMesh.py] DEBUG ampD, ampB, DriveW_CurrentAmp, BeariW_CurrentAmp:', ampD, ampB, EX['DriveW_CurrentAmp'], EX['BeariW_CurrentAmp'])
        print('[FEMM_SllidingMesh.py] circuit turns for combined winding:', turnsDPNV)
        # quit()

        # static solver
        omegaDrive = 2*np.pi*self.acm_variant.template.d['EX']['DriveW_Freq']
        omegaBeari = 2*np.pi*self.acm_variant.template.d['EX']['BeariW_Freq']
        phase_shift_drive = -120 if wily.CommutatingSequenceD == 1 else 120
        phase_shift_beari = -120 if wily.CommutatingSequenceB == 1 else 120
        self.dict_stator_current_function = []
        self.dict_stator_current_function.append(lambda t, ampD=ampD, ampB=ampB, phase_shift_drive=phase_shift_drive, phase_shift_beari=phase_shift_beari: ampD * np.sin(omegaDrive*t + 0*phase_shift_drive/180*np.pi) \
                                                                                                                                                         - ampB * np.sin(omegaBeari*t + 0*phase_shift_beari/180*np.pi))
        self.dict_stator_current_function.append(lambda t, ampD=ampD, ampB=ampB, phase_shift_drive=phase_shift_drive, phase_shift_beari=phase_shift_beari: ampD * np.sin(omegaDrive*t + 1*phase_shift_drive/180*np.pi) \
                                                                                                                                                         - ampB * np.sin(omegaBeari*t + 1*phase_shift_beari/180*np.pi))
        self.dict_stator_current_function.append(lambda t, ampD=ampD, ampB=ampB, phase_shift_drive=phase_shift_drive, phase_shift_beari=phase_shift_beari: ampD * np.sin(omegaDrive*t + 2*phase_shift_drive/180*np.pi) \
                                                                                                                                                         - ampB * np.sin(omegaBeari*t + 2*phase_shift_beari/180*np.pi))
        self.dict_stator_current_function.append(lambda t, ampD=ampD, ampB=ampB, phase_shift_drive=phase_shift_drive, phase_shift_beari=phase_shift_beari: ampD * np.sin(omegaDrive*t + 0*phase_shift_drive/180*np.pi) \
                                                                                                                                                         + ampB * np.sin(omegaBeari*t + 0*phase_shift_beari/180*np.pi))
        self.dict_stator_current_function.append(lambda t, ampD=ampD, ampB=ampB, phase_shift_drive=phase_shift_drive, phase_shift_beari=phase_shift_beari: ampD * np.sin(omegaDrive*t + 1*phase_shift_drive/180*np.pi) \
                                                                                                                                                         + ampB * np.sin(omegaBeari*t + 1*phase_shift_beari/180*np.pi))
        self.dict_stator_current_function.append(lambda t, ampD=ampD, ampB=ampB, phase_shift_drive=phase_shift_drive, phase_shift_beari=phase_shift_beari: ampD * np.sin(omegaDrive*t + 2*phase_shift_drive/180*np.pi) \
                                                                                                                                                         + ampB * np.sin(omegaBeari*t + 2*phase_shift_beari/180*np.pi))

        femm.mi_addcircprop('U-GrpAC', self.dict_stator_current_function[0](0.0), SERIES_CONNECTED)
        femm.mi_addcircprop('V-GrpAC', self.dict_stator_current_function[1](0.0), SERIES_CONNECTED)
        femm.mi_addcircprop('W-GrpAC', self.dict_stator_current_function[2](0.0), SERIES_CONNECTED)
        femm.mi_addcircprop('U-GrpBD', self.dict_stator_current_function[3](0.0), SERIES_CONNECTED)
        femm.mi_addcircprop('V-GrpBD', self.dict_stator_current_function[4](0.0), SERIES_CONNECTED)
        femm.mi_addcircprop('W-GrpBD', self.dict_stator_current_function[5](0.0), SERIES_CONNECTED)

        Qs = self.acm_variant.template.SI['Qs']

        # dict_dir = {'+':1, '-':-1} # wrong (not consistent with JMAG)
        dict_dir = {'+':-1, '-':1, 'o':0}
        R = GP['mm_r_si'].value + 0.8 * GP['mm_d_st'].value
        angle_per_slot = 2*np.pi/Qs

        # X layer winding's blocks
        ''' TODO Need to make sure the block is inside the coil object
        '''
        self.acm_variant.coils.innerCoord
        OFFSET = 2 / 180 * np.pi
        THETA = - angle_per_slot + 0.5*angle_per_slot - OFFSET # This 3 deg must be less than 360/Qs/2，取这么大是为了在GUI上看得清楚点。
        count = 0
        # for phase, up_or_down in zip(im.l_rightlayer1,im.l_rightlayer2):
        for phase, up_or_down, AC_or_BD in zip(wily.layer_X_phases, wily.layer_X_signs, wily.grouping_AC):
            # circuit_name = 'd' + phase
            circuit_name = phase + '-Grp' + ('AC' if AC_or_BD else 'BD')
            THETA += angle_per_slot
            X = R*np.cos(THETA); Y = R*np.sin(THETA)
            count += 1
            if fraction == 4:
                if not (count > Qs*0.5+EPS and count <= Qs*0.75+EPS): 
                    continue
            if fraction == 2:
                if not (count > Qs*0.5+EPS): 
                    continue
            # if phase == 'U':
            self.set_block_label(11, 'Copper', (X, Y), MESH_SIZE_COPPER, turns=turnsDPNV*dict_dir[up_or_down], automesh=self.bool_automesh, incircuit=circuit_name)

        # Y layer winding's blocks
        if fraction == 1:
            THETA = - angle_per_slot + 0.5*angle_per_slot + OFFSET

            grouping_AC_of_Y_layer = winding_layout.infer_Y_layer_grpAC_from_X_layer_and_coil_pitch_y(wily.grouping_AC, wily.coil_pitch_y)
            for phase, up_or_down, AC_or_BD in zip(wily.layer_Y_phases, wily.layer_Y_signs, grouping_AC_of_Y_layer):
                # circuit_name = 'b' + phase
                circuit_name = phase + '-Grp' + ('AC' if AC_or_BD else 'BD')
                THETA += angle_per_slot
                X = R*np.cos(THETA); Y = R*np.sin(THETA)

                self.set_block_label(11, 'Copper', (X, Y), MESH_SIZE_COPPER, turns=turnsDPNV*dict_dir[up_or_down], automesh=self.bool_automesh, incircuit=circuit_name)

        # 危险！FEMM默认把没有设置incircuit的导体都在无限远短接在一起——也就是说，你可能把定子悬浮绕组也短接到鼠笼上去了！
        # 所以，一定要设置好悬浮绕组，而且要用serial-connected，电流给定为 0 A。


    def add_sliding_band_and_boundary(self, deg_innerAngle, deg_OuterAngle=0, fraction=1):
        # if self.deg_per_step != 0.0:
        #     # 之前用的方法是打开上一个FEM文件，然后将其模型旋转deg_per_step，用不到rotor_position_in_deg的！
        #     # 当然，我们也试过AirGapBoundary（David Meeker推荐的），转动转子不需要重复剖分，但是计算出来的力不准（转矩是准的）
        #     # 现在，我们打开在0位置的fem文件，然后转动，saveas。这样，就不用不断地打开文件了
        #     femm.mi_selectgroup(100) # this only select the block labels
        #     femm.mi_selectgroup(101)
        #     femm.mi_selectcircle(0,0,self.acm_variant.template.d['GP']['mm_r_or'].value+EPS,SELECT_ALL) # this selects the nodes, segments, arcs
        #     femm.mi_moverotate(0,0, self.deg_per_step)
        #     # femm.mi_zoomnatural()

        GP = self.acm_variant.template.d['GP']
        Point_Stator = [GP['mm_r_si'].value - 0.3*GP['mm_d_fixed_air_gap'].value, 0]
        Point_Rotor  = [GP['mm_r_si'].value - 0.6*GP['mm_d_fixed_air_gap'].value, 0]
        self.drawArc([0,0], Point_Stator, [-Point_Stator[0], Point_Stator[1]])
        self.drawArc([0,0],               [-Point_Stator[0], Point_Stator[1]], Point_Stator)
        self.drawArc([0,0], Point_Rotor,  [-Point_Rotor[0], Point_Rotor[1]])
        self.drawArc([0,0],               [-Point_Rotor[0], Point_Rotor[1]], Point_Rotor)

        # sliding band is implemented as anti-... in FEMM
        femm.mi_addboundprop('WholeModelSlidingBand', 0,0,0, 0,0,0, 0,0, 6,deg_innerAngle,deg_OuterAngle) # Periodic Air Gap (see https://www.femm.info/wiki/slidingbandbenchmark)
        # 注意圆弧分成两段来画了
        femm.mi_selectarcsegment(0,  Point_Stator[0])
        femm.mi_setarcsegmentprop(1, "WholeModelSlidingBand", False, 137) # maxseg = 1
        femm.mi_clearselected()
        femm.mi_selectarcsegment(0, -Point_Stator[0])
        femm.mi_setarcsegmentprop(1, "WholeModelSlidingBand", False, 137) # maxseg = 1
        femm.mi_clearselected()
        femm.mi_selectarcsegment(0,  Point_Rotor[0])
        femm.mi_setarcsegmentprop(1, "WholeModelSlidingBand", False, 137) # maxseg = 1
        femm.mi_clearselected()
        femm.mi_selectarcsegment(0, -Point_Rotor[0])
        femm.mi_setarcsegmentprop(1, "WholeModelSlidingBand", False, 137) # maxseg = 1
        femm.mi_clearselected()
        self.set_block_label(102, '<No Mesh>',         (0, 0.5*Point_Stator[0]+0.5*Point_Rotor[0]),  20)

        # Boundary Conditions 
        # femm.mi_makeABC() # open boundary
        if fraction == 1:
            femm.mi_addboundprop('BC:A=0', 0,0,0, 0,0,0,0,0,0,0,0)

            femm.mi_selectarcsegment(0,-GP['mm_r_os'].value)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 10) # maxseg = 20 deg (only this is found effective)
            femm.mi_clearselected()
            femm.mi_selectarcsegment(0,GP['mm_r_os'].value)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 10)
            femm.mi_clearselected()

            femm.mi_selectarcsegment(0,-GP['mm_r_ri'].value)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 100)
            femm.mi_clearselected()
            femm.mi_selectarcsegment(0,GP['mm_r_ri'].value)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 100)
            femm.mi_clearselected()

        # Other arc-segment-specific mesh constraints are already done in draw_model()

    def run_transient_study(self):
        # local references
        # GP = self.acm_variant.template.d['GP']
        # EX = self.acm_variant.template.d['EX']
        # p  = self.acm_variant.template.SI['p']
        fea_config_dict = self.acm_variant.template.fea_config_dict

        # figure out initial rotor angle
        # deg_pole_span = 180/p
        # initial_rotor_position_mech_deg = (deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + EX['wily'].deg_winding_U_phase_phase_axis_angle + deg_pole_span
        # print('[FEMM_SlidingMesh.py] initial_rotor_position_mech_deg:', initial_rotor_position_mech_deg, 'deg || deg_winding_U_phase_phase_axis_angle:', '', EX['wily'].deg_winding_U_phase_phase_axis_angle)

        # figure out step size
        # electrical_period = fea_config_dict['femm.number_cycles_in_2ndTSS']/EX['DriveW_Freq']
        # number_of_steps   = fea_config_dict['femm.number_of_steps_2ndTSS']
        # step_size_sec = electrical_period / (number_of_steps-1)
        # step_size_mech_deg = EX['Omega'] * step_size_sec / np.pi * 180

        # transient FEA with sliding band (air gap boundary in FEMM)
        print('[FEMM_SlidingMesh.py] Run FEMM solver...')
        self.acm_variant.analyzer = Individual_Analyzer_FEMM_Edition(p=self.acm_variant.template.SI['p'], 
                                                                     ns=fea_config_dict['femm.number_of_steps_2ndTSS'],
                                                                     GroupSummary=self.GroupSummary)

        if fea_config_dict['femm.number_of_parallel_solve'] == 1:
            for index in range(self.number_of_steps):

                # update excitations
                time                         = index * self.step_size_sec
                RotorAngle_MechanicalDegrees = index * self.step_size_mech_deg
                current_sources = self.update_circuit_excitation(time)

                # update rotor position
                femm.mi_modifyboundprop('WholeModelSlidingBand', 10, self.initial_rotor_position_mech_deg+RotorAngle_MechanicalDegrees+self.acm_variant.template.fea_config_dict['femm.MechDeg_IdEqualToNonZeroAngle']) # change inner angle to 

                # if index > 0:
                #     femm.mi_setprevious(self.project_file_name[:-4] + f'-{index-1:03d}.ans', 1) # 1: must be incremental permeability (because frozen permability means the permeability is fixed).
                # Starting nonlinear field from previous solution is not possible with current FEMM's function


                femm.mi_saveas(self.project_file_name[:-4] + f'-{index:03d}.fem')
                # femm.mi_createmesh() # no need
                femm.mi_analyze(1)

                # collect results at: index
                femm.mi_loadsolution()
                self.collect_femm_results(index, time, RotorAngle_MechanicalDegrees, current_sources)

        else:
            for index in range(self.number_of_steps):
                time                         = index * self.step_size_sec
                RotorAngle_MechanicalDegrees = index * self.step_size_mech_deg
                current_sources = self.update_circuit_excitation(time)
                femm.mi_modifyboundprop('WholeModelSlidingBand', 10, self.initial_rotor_position_mech_deg+RotorAngle_MechanicalDegrees+self.acm_variant.template.fea_config_dict['femm.MechDeg_IdEqualToNonZeroAngle']) # change inner angle to 
                femm.mi_saveas(self.project_file_name[:-4] + f'-{index:03d}.fem')
                if os.path.exists(self.project_file_name[:-4] + f'-{index:03d}.ans'):
                    os.remove(self.project_file_name[:-4] + f'-{index:03d}.ans')
            self.parallel_solve_transient_FEA(self.step_size_sec, self.step_size_mech_deg)

        femm.mi_close()

    def collect_femm_results(self, index, time, RotorAngle_MechanicalDegrees, current_sources, bool_parallel=False):

        rotation_operator = np.exp(1j*RotorAngle_MechanicalDegrees/180*np.pi)

        if False:
            if index == 0:
                self.initialize_collection_of_field_data()

            # Field data @ anywhere
            for i, el in enumerate(self.list_of_femm_triangle_element):
                Bx, By = femm.mo_getb(el.x, el.y)
                self.list_of_femm_triangle_element[i].femm_Bx_list.append(Bx)
                self.list_of_femm_triangle_element[i].femm_By_list.append(By)

            # Show flux density plot by element
            plt.figure()
            B = np.array([np.sqrt(el.femm_Bx_list[index]**2+el.femm_By_list[index]**2) for el in self.list_of_femm_triangle_element if el.group==10 and el.x>0 and el.y>0]); scaled_B = (B - B.min()) / B.ptp()
            plt.scatter( [el.x                                                         for el in self.list_of_femm_triangle_element if el.group==10 and el.x>0 and el.y>0], 
                         [el.y                                                         for el in self.list_of_femm_triangle_element if el.group==10 and el.x>0 and el.y>0], 
                         edgecolors=plt.cm.copper(scaled_B), color='white', marker='.')

            B = np.array([np.sqrt(el.femm_Bx_list[index]**2+el.femm_By_list[index]**2) for el in self.list_of_femm_triangle_element if el.group==100 and el.x>0 and el.y>0]); scaled_B = (B - B.min()) / B.ptp()
            plt.scatter( [el.x                                                         for el in self.list_of_femm_triangle_element if el.group==100 and el.x>0 and el.y>0], 
                         [el.y                                                         for el in self.list_of_femm_triangle_element if el.group==100 and el.x>0 and el.y>0], 
                         edgecolors=plt.cm.copper(scaled_B), color='white', marker='.')

            B = np.array([np.sqrt(el.femm_Bx_list[index]**2+el.femm_By_list[index]**2) for el in self.list_of_femm_triangle_element if el.group==101 and el.x>0 and el.y>0]); scaled_B = (B - B.min()) / B.ptp()
            plt.scatter( [el.x                                                         for el in self.list_of_femm_triangle_element if el.group==101 and el.x>0 and el.y>0], 
                         [el.y                                                         for el in self.list_of_femm_triangle_element if el.group==101 and el.x>0 and el.y>0], 
                         edgecolors=plt.cm.copper(scaled_B), color='white', marker='.')
            plt.gca().set_aspect('equal')
            plt.show()

        # Iron loss: field (B) and potential (A) data
        self.acm_variant.analyzer.iron_loss_preparing_data_matrices(index)

        # Torque, force, fluxes
        torque = femm.mo_gapintegral('WholeModelSlidingBand',0)
        forces = femm.mo_gapintegral('WholeModelSlidingBand',1)
        energy = femm.mo_gapintegral('WholeModelSlidingBand',2)
        circuitProperties = ( femm.mo_getcircuitproperties('U-GrpAC'),
                              femm.mo_getcircuitproperties('V-GrpAC'),
                              femm.mo_getcircuitproperties('W-GrpAC'),
                              femm.mo_getcircuitproperties('U-GrpBD'),
                              femm.mo_getcircuitproperties('V-GrpBD'),
                              femm.mo_getcircuitproperties('W-GrpBD') )
        if bool_parallel:
            self.acm_variant.analyzer.add(time, RotorAngle_MechanicalDegrees, torque, forces, energy, circuitProperties, index=index)
        else:
            self.acm_variant.analyzer.add(time, RotorAngle_MechanicalDegrees, torque, forces, energy, circuitProperties, index=None)
        femm.mo_close()

        # print to screen
        print(f'\t {index}, {time*1000} ms, {RotorAngle_MechanicalDegrees} deg : {torque:.2f} Nm, [{forces[0]:.2f} N, {forces[1]:.2f}] N, {energy:.2f} J', end=' | ')
        print(' A, '.join(f'{current:.1f}' for current in current_sources), 'A')
        for el in circuitProperties:
            print('\t\t', el)

    def compute_objectives(self, select_fea_config_dict):

        self.acm_variant.analyzer.get_ss_data()
        self.acm_variant.analyzer.compute_objectives(self.acm_variant, select_fea_config_dict, self)

        # Bx = list(map(lambda el: el[0], self.stator_Bx_data))
        # By = list(map(lambda el: el[0], self.stator_By_data))
        # plt.plot(Bx, By, 'o')
        # plt.scatter(x, y, marker='+', s=150, linewidths=4, c=y, cmap=plt.cm.coolwarm)
        # plt.show()

    def save_results_to_disk(self):
        pass

    def update_circuit_excitation(self, time):
        # rotor current
        # for i in range(self.rotor_slot_per_pole):
        #     circuit_name = 'r%s'%(self.rotor_phase_name_list[i])
        #     femm.mi_modifycircprop(circuit_name, 1, self.dict_rotor_current_function[i](time))

        # stator current
        femm.mi_modifycircprop('U-GrpAC', 1, self.dict_stator_current_function[0](time))
        femm.mi_modifycircprop('V-GrpAC', 1, self.dict_stator_current_function[1](time))
        femm.mi_modifycircprop('W-GrpAC', 1, self.dict_stator_current_function[2](time))
        femm.mi_modifycircprop('U-GrpBD', 1, self.dict_stator_current_function[3](time))
        femm.mi_modifycircprop('V-GrpBD', 1, self.dict_stator_current_function[4](time))
        femm.mi_modifycircprop('W-GrpBD', 1, self.dict_stator_current_function[5](time))

        return [
            self.dict_stator_current_function[0](time),
            self.dict_stator_current_function[1](time),
            self.dict_stator_current_function[2](time),
            self.dict_stator_current_function[3](time),
            self.dict_stator_current_function[4](time),
            self.dict_stator_current_function[5](time),
        ]

    def parallel_solve_transient_FEA(self, step_size_sec, step_size_mech_deg):

        fea_config_dict = self.acm_variant.template.fea_config_dict
        number_of_parallel_solve = fea_config_dict['femm.number_of_parallel_solve']
        number_of_points = fea_config_dict['femm.number_of_steps_2ndTSS']
        number_of_points_per_solve = number_of_points // number_of_parallel_solve
        procs = []
        for i in range(number_of_parallel_solve):
            proc = subprocess.Popen([sys.executable, 'parasolveSlidingMesh.py', str(i), str(number_of_parallel_solve), str(number_of_points), fea_config_dict['output_dir'], self.project_file_name], bufsize=-1)
            procs.append(proc)

        for proc in procs:
            code = proc.wait() # return exit code
            print('[FEMM_SlidingMesh.py] DEBUG process return code:', code)

        ''' Collecting results after all .fem files are solved (easy to code but apparently this takes more time for colelcting)
        '''
        list_of_completed_ans_file = [0] * number_of_points
        flag_run_while = True
        print('[FEMM_SlidingMesh.py] Start collecting .ans results one by one...')
        while flag_run_while:
            for index, _bool in enumerate(list_of_completed_ans_file):
                if _bool == 0:
                    # try to load .ans file
                    try:
                        femm.opendocument(self.project_file_name[:-4] + f'-{index:03d}.ans')
                        time                         = index * step_size_sec
                        RotorAngle_MechanicalDegrees = index * step_size_mech_deg
                        self.collect_femm_results(index, time, RotorAngle_MechanicalDegrees, current_sources=[])
                        list_of_completed_ans_file[index] = True
                    except Exception as e:
                        raise e
            flag_run_while = sum(list_of_completed_ans_file) != number_of_points
            print('\t\t', list_of_completed_ans_file)

    def initialize_collection_of_field_data(self):

        # get slot area for copper loss calculation
        # femm.mo_groupselectblock(11) # fraction is 1 
        # self.stator_slot_area = femm.mo_blockintegral(5) / self.im.Qs # unit: m^2 (verified by GUI operation)
        # femm.mo_clearblock()

        # femm.mo_groupselectblock(101)
        # self.rotor_slot_area = femm.mo_blockintegral(5) / self.im.Qr
        # femm.mo_clearblock()

        self.number_of_elements = femm.mo_numelements()
        if True:
            self.list_of_femm_triangle_element = []
            self.stator_Area_data = []
            self.stator_xy_complex_data = []
            self.rotor_Area_data = []
            self.rotor_xy_complex_data = []
            print('[FEMM_SlidingMesh.py] DEBUG 1/4 model for computing iron loss...')
            for id_element in range(1, self.number_of_elements+1):
                indexNode1, indexNode2, indexNode3, x, y, area, group = femm.mo_getelement(id_element)
                self.list_of_femm_triangle_element.append(utility.femm_triangle_element(indexNode1, indexNode2, indexNode3, x, y, area, group, [], []))
            print('[FEMM_SlidingMesh.py] elements:', id_element, self.number_of_elements)

            # debug
            colors = {
                9: None,
                10: 'red',
                11: None,
                100: 'black',
                101: 'blue',
                102: None,
            }
            from pylab import plt
            plt.figure()
            plt.scatter( [el.x for el in self.list_of_femm_triangle_element if el.group==10 and el.x>0 and el.y>0], 
                         [el.y for el in self.list_of_femm_triangle_element if el.group==10 and el.x>0 and el.y>0], 
                                                                            color=colors[10], marker='.')
            plt.scatter( [el.x for el in self.list_of_femm_triangle_element if el.group==100 and el.x>0 and el.y>0], 
                         [el.y for el in self.list_of_femm_triangle_element if el.group==100 and el.x>0 and el.y>0], 
                                                                            color=colors[100], marker='.')
            plt.scatter( [el.x for el in self.list_of_femm_triangle_element if el.group==101 and el.x>0 and el.y>0], 
                         [el.y for el in self.list_of_femm_triangle_element if el.group==101 and el.x>0 and el.y>0], 
                                                                            color=colors[101], marker='.')
            plt.gca().set_aspect('equal')
            plt.show()
        else:
            self.stator_Area_data = []
            self.stator_xy_complex_data = []
            self.rotor_Area_data = []
            self.rotor_xy_complex_data = []
            print('[FEMM_SlidingMesh.py] DEBUG 1/4 model for computing iron loss...')
            for id_element in range(1, self.number_of_elements+1):
                _, _, _, x, y, area, group = femm.mo_getelement(id_element)

                # consider 1/4 model for loss (this is valid if we presume the suspension two pole field is weak)
                # Use the mesh info of the initial rotor position 
                if group == 10: # stator iron
                    new_xy_complex = (x+1j*y) * np.exp(1j* np.pi/self.acm_variant.template.SI['Qs']) # 反正是选1/4模型，我们选使得分割线经过槽的1/4模型。
                    if new_xy_complex.real>0 and new_xy_complex.imag>0:
                        self.stator_Area_data.append(area)
                        self.stator_xy_complex_data.append(x+1j*y)
                elif group == 100: # rotor iron
                    if y>0 and x>0:
                        self.rotor_Area_data.append(area)
                        self.rotor_xy_complex_data.append(x+1j*y)
                elif group == 101: # magnet
                    if y>0 and x>0:
                        self.rotor_Area_data.append(area)
                        self.rotor_xy_complex_data.append(x+1j*y)

            self.stator_Bx_data = []
            self.stator_By_data = []
            for i in range(len(self.stator_xy_complex_data)):
                self.stator_Bx_data.append([])
                self.stator_By_data.append([])
            self.rotor_Bx_data = []
            self.rotor_By_data = []
            for i in range(len(self.rotor_xy_complex_data)):
                self.rotor_Bx_data.append([])
                self.rotor_By_data.append([])

            # debug
            from pylab import plt
            plt.figure()
            for ind, complex_number in enumerate(self.stator_xy_complex_data):
                plt.scatter( complex_number.real, complex_number.imag)
            plt.gca().set_aspect('equal')

            plt.figure()
            for ind, complex_number in enumerate(self.rotor_xy_complex_data):
                plt.scatter( complex_number.real, complex_number.imag)        
            plt.gca().set_aspect('equal')

            plt.show()

    def getSketch(self, sketchName, color=None):
        pass

    ''' Drawer
    '''
    def draw_spmsm(self, acm_variant):
        # Rotor Core
        list_regions_1 = acm_variant.rotorCore.draw(self, bool_draw_whole_model=True)
        self.bMirror = False
        self.iRotateCopy = acm_variant.rotorCore.p*2
        # region1 = self.prepareSection(list_regions_1, color=color_rgb_A)

        # Shaft
        # list_regions = acm_variant.shaft.draw(self)
        # self.bMirror = False
        # self.iRotateCopy = 1
        # region0 = self.prepareSection(list_regions)

        # Rotor Magnet
        list_regions = acm_variant.rotorMagnet.draw(self, bool_draw_whole_model=True)
        self.bMirror = False
        self.iRotateCopy = acm_variant.rotorMagnet.notched_rotor.p*2
        # region2 = self.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

        # Sleeve
        # list_regions = acm_variant.sleeve.draw(self)
        # self.bMirror = False
        # self.iRotateCopy = acm_variant.rotorMagnet.notched_rotor.p*2
        # regionS = self.prepareSection(list_regions)

        # Stator Core
        list_regions = acm_variant.stator_core.draw(self, bool_draw_whole_model=True)
        self.bMirror = True
        self.iRotateCopy = acm_variant.stator_core.Q
        # region3 = self.prepareSection(list_regions, color=color_rgb_A)

        # Stator Winding
        list_regions = acm_variant.coils.draw(self, bool_draw_whole_model=True)
        self.bMirror = False
        self.iRotateCopy = acm_variant.coils.stator_core.Q
        # region4 = self.prepareSection(list_regions)

        # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
        EX = acm_variant.template.d['EX']
        CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * EX['WindingFill'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
        CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
        CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。

        # Maybe there is a bug here... regarding the excitation for suspension winding...
        variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
        variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
        EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
        EX['DriveW_CurrentAmp'] = acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
        EX['BeariW_CurrentAmp'] = acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

        slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
        print('[FEMM_SlidingMesh.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

        return True

    @staticmethod
    def mirror_and_copyrotate(Q, Radius, fraction):
        # Mirror
        femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL) # this EPS is sometime necessary to selece the arc at Radius.
        femm.mi_mirror2(0,0,-Radius,0, SELECT_ALL)

        # Rotate
        femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL)
        femm.mi_copyrotate2(0, 0, 360./Q, int(Q)/fraction, SELECT_ALL)

    @staticmethod
    def draw_arc(p1, p2, angle, maxseg=1, center=None, **kwarg):
        femm.mi_drawarc(p1[0],p1[1],p2[0],p2[1],angle/np.pi*180,maxseg) # [deg]

    @staticmethod
    def add_arc(p1, p2, angle, maxseg=1, center=None, **kwarg):
        femm.mi_addarc(p1[0],p1[1],p2[0],p2[1],angle/np.pi*180,maxseg) # [deg]

    @staticmethod
    def draw_line(p1, p2):
        femm.mi_drawline(p1[0],p1[1],p2[0],p2[1])

    @staticmethod
    def add_line(p1, p2):
        femm.mi_addsegment(p1[0],p1[1],p2[0],p2[1])

    @staticmethod
    def drawLine(p1, p2):
        femm.mi_drawline(p1[0],p1[1],p2[0],p2[1])
        return []

    @staticmethod
    def drawArc(centerxy, startxy, endxy, maxseg=1):
        
        v1 = np.array([startxy[0] - centerxy[0], startxy[1] - centerxy[1]])
        v2 = np.array([endxy[0]   - centerxy[0], endxy[1]   - centerxy[1]])

        cos夹角 = (v1[0]*v2[0] + v1[1]*v2[1]) / (np.sqrt(v1.dot(v1))*np.sqrt(v2.dot(v2)))
        angle = np.arccos(cos夹角)
        if angle == 0:
            angle = 360

        p1, p2 = startxy, endxy
        femm.mi_drawarc(p1[0],p1[1],p2[0],p2[1],angle/np.pi*180,maxseg)
        return []


    def some_solver_related_operations_rotor_before_mirror_rotation(self, im, P6, P8):

        if im.use_drop_shape_rotor_bar == True:
            # constraint to reduce element number @rotor-P6
            femm.mi_selectarcsegment(P6[0], P6[1])
            femm.mi_setarcsegmentprop(8, "<None>", False, 100)
            femm.mi_clearselected()

            # constraint to reduce element number @rotor-P8
            femm.mi_selectarcsegment(P8[0], P8[1])
            femm.mi_setarcsegmentprop(8, "<None>", False, 100)
            femm.mi_clearselected()
        else:
            # constraint to reduce element number @rotor-P8
            femm.mi_selectarcsegment(P8[0], P8[1])
            femm.mi_setarcsegmentprop(8, "<None>", False, 100)
            femm.mi_clearselected()

    def some_solver_related_operations_fraction(self, im, fraction):
        # Boundary
        if fraction == 1:
            femm.mi_drawarc(im.Radius_Shaft,0, -im.Radius_Shaft,0, 180, 20) # 边界不要用太小的segment咯！避免剖分过细（这里设置无效）
            femm.mi_drawarc(-im.Radius_Shaft,0, im.Radius_Shaft,0, 180, 20)
            femm.mi_drawarc(im.Radius_OuterStatorYoke,0, -im.Radius_OuterStatorYoke,0, 180, 20)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, im.Radius_OuterStatorYoke,0, 180, 20)
        elif fraction == 4:
            femm.mi_drawarc(-im.Radius_Shaft,0, 0, -im.Radius_Shaft, 90, 10)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, 0, -im.Radius_OuterStatorYoke, 90, 10)
            femm.mi_selectrectangle(-EPS-im.Radius_Shaft,EPS,EPS-im.Radius_OuterStatorYoke,im.Radius_OuterStatorYoke,SELECT_ALL)
            femm.mi_selectrectangle(EPS,-EPS-im.Radius_Shaft,im.Radius_OuterStatorYoke,EPS-im.Radius_OuterStatorYoke,SELECT_ALL)
            femm.mi_deleteselected()

            # between 2rd and 3th quarters
            p1 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
            p2 = (-im.Radius_Shaft, 0)
            self.add_line(p1, p2)
            p2 = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
            self.add_line(p1, p2)
            p1 = (-im.Radius_OuterRotor-0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            self.draw_line(p1, p2)
            p2 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)
            self.draw_line(p1, p2)
            p1 = (-im.Radius_OuterStatorYoke, 0)
            self.add_line(p1, p2)

            # between 3rd and 4th quarters
            p1 = (0, -im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2)
            p2 = (0, -im.Radius_Shaft)
            self.add_line(p1, p2)
            p2 = (0, -im.Location_RotorBarCenter-im.Radius_of_RotorSlot)
            self.add_line(p1, p2)
            p1 = (0, -im.Radius_OuterRotor-0.5*im.Length_AirGap)
            self.draw_line(p1, p2)
            p2 = (0, -im.Radius_OuterRotor-im.Length_AirGap)
            self.draw_line(p1, p2)
            p1 = (0, -im.Radius_OuterStatorYoke)
            self.add_line(p1, p2)
        elif fraction == 2:
            femm.mi_drawarc(-im.Radius_Shaft,0, im.Radius_Shaft,0, 180, 15)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, im.Radius_OuterStatorYoke,0, 180, 15)
            femm.mi_selectrectangle(EPS-im.Radius_OuterStatorYoke,EPS, -EPS+im.Radius_OuterStatorYoke,EPS+im.Radius_OuterStatorYoke, SELECT_ALL)
            femm.mi_deleteselected()

            # between 2rd and 3th quarters
            p1 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
            p2 = (-im.Radius_Shaft, 0)
            self.add_line(p1, p2)
            p2 = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
            self.add_line(p1, p2)
            p1 = (-im.Radius_OuterRotor-0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            self.draw_line(p1, p2)
            p2 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)
            self.draw_line(p1, p2)
            p1 = (-im.Radius_OuterStatorYoke, 0)
            self.add_line(p1, p2)

            # between 1rd and 4th quarters
            p1 = (+im.Location_RotorBarCenter2-im.Radius_of_RotorSlot2, 0)
            p2 = (+im.Radius_Shaft, 0)
            self.add_line(p1, p2)
            p2 = (+im.Location_RotorBarCenter+im.Radius_of_RotorSlot, 0)
            self.add_line(p1, p2)
            p1 = (+im.Radius_OuterRotor+0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            self.draw_line(p1, p2)
            p2 = (+im.Radius_OuterRotor+im.Length_AirGap, 0)
            self.draw_line(p1, p2)
            p1 = (+im.Radius_OuterStatorYoke, 0)
            self.add_line(p1, p2)
        else:
            raise Exception('not supported fraction = %d' % (fraction))
        # Air Gap Boundary for Rotor Motion #1
        # R = im.Radius_OuterRotor+0.6*im.Length_AirGap
        # femm.mi_drawarc(R,0, -R,0, 180, 5)
        # femm.mi_drawarc(-R,0, R,0, 180, 5)
        # R = im.Radius_OuterRotor+0.4*im.Length_AirGap
        # femm.mi_drawarc(R,0, -R,0, 180, 5)
        # femm.mi_drawarc(-R,0, R,0, 180, 5)






    def get_air_gap_B(self, number_of_points=360):
        im = self.im
        femm.opendocument(self.output_file_name + '.fem')
        femm.mi_loadsolution()

        list_B_magitude = []
        R = im.Radius_OuterRotor + 0.25*im.Length_AirGap
        for i in range(number_of_points):
            THETA = i / 180.0 * np.pi
            X = R*np.cos(THETA)
            Y = R*np.sin(THETA)
            B_vector_complex = femm.mo_getb(X, Y)
            B_X_complex = B_vector_complex[0]
            B_Y_complex = B_vector_complex[1]
            B_X_real = np.real(B_X_complex)
            B_Y_real = np.real(B_Y_complex)
            # Assume the magnitude is all due to radial component
            B_magitude = np.sqrt(B_X_real**2 + B_Y_real**2)
            inner_product = B_X_real * X + B_Y_real *Y
            list_B_magitude.append( B_magitude * utility.copysign(1, inner_product) )
        return list_B_magitude

    def femm_integrate_4_current(self, fname, fraction, dir_output=None, returnData=False):
        '''Make sure femm is opened
        Returns:
            [type] -- [list of complex number of rotor currents from FEMM]
        '''

        # get corresponding rotor current conditions for later static FEA
        femm.opendocument(fname)

        if True:
            # physical amount of Cage
            im = self.im
            vals_results_rotor_current = []
            # R = 0.5*(im.Location_RotorBarCenter + im.Location_RotorBarCenter2)
            R = im.Location_RotorBarCenter # Since 5/23/2019
            angle_per_slot = 2*np.pi/im.Qr
            THETA_BAR = np.pi - angle_per_slot + EPS # add EPS for the half bar
            # print 'number of rotor_slot per partial model', self.rotor_slot_per_pole * int(4/fraction)
            for i in range(self.rotor_slot_per_pole * int(4/fraction)):
                THETA_BAR += angle_per_slot
                THETA = THETA_BAR
                X = R*np.cos(THETA); Y = R*np.sin(THETA)
                femm.mo_selectblock(X, Y) # or you can select circuit rA rB ...
                vals_results_rotor_current.append(femm.mo_blockintegral(7)) # integrate for current
                femm.mo_clearblock()
            # the other half bar of rA
            THETA_BAR += angle_per_slot
            THETA = THETA_BAR - 2*EPS
            X = R*np.cos(THETA); Y = R*np.sin(THETA)
            femm.mo_selectblock(X, Y)
            vals_results_rotor_current.append(femm.mo_blockintegral(7)) # integrate for current
            femm.mo_clearblock()

            ################################################################
            # Also collect slot area information for loss evaluation in JMAG optimization 
            ################################################################
            if True:
                # get stator slot area for copper loss calculation
                femm.mo_groupselectblock(11)
                stator_slot_area = femm.mo_blockintegral(5) / (im.Qs/fraction) # unit: m^2 (verified by GUI operation)
                femm.mo_clearblock()

                # get rotor slot area for copper loss calculation
                femm.mo_groupselectblock(101)
                rotor_slot_area = femm.mo_blockintegral(5) / (im.Qs/fraction)
                femm.mo_clearblock()

            femm.mo_close()
            # return [-el for el in vals_results_rotor_current[self.rotor_slot_per_pole:2*self.rotor_slot_per_pole]] # 用第四象限的转子电流，因为第三象限的被切了一半，麻烦！
            # vals_results_rotor_current[self.rotor_slot_per_pole:2*self.rotor_slot_per_pole]这里用的都是第四象限的转子电流了，我们后面默认用的是第三象限的转子电流，即rA1 rB1 ...，所以要反相一下(-el)
            vals_results_rotor_current = [-el for el in vals_results_rotor_current[self.rotor_slot_per_pole:2*self.rotor_slot_per_pole]]
        # vals_results_rotor_current = self.femm_integrate_4_current(self.fraction)

        if dir_output is None:
            dir_output = self.dir_run_sweeping

        if returnData == False: # no return then write to file
            with open(dir_output + "femm_rotor_current_conditions.txt", "w") as stream:
                str_results = ''
                for el in vals_results_rotor_current:
                    stream.write("%g %g \n" % (el.real, el.imag))
            print('done. append to eddycurrent_results.txt.')
            return None
        else:
            return vals_results_rotor_current, stator_slot_area, rotor_slot_area

    def read_Torque_and_B_data(self, ans_file, rotation_operator, handle_torque):
        # (str_rotor_position, rotation_operator):
        femm.opendocument(ans_file)

        # Physical Amount on the Rotor
        femm.mo_groupselectblock(100) # rotor iron
        femm.mo_groupselectblock(101) # rotor bars
        Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
        Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
        torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
        femm.mo_clearblock()
        # write results to a data file (write to partial files to avoid compete between parallel instances)
        handle_torque.write("%s %g %g %g\n" % ( ans_file[-8:-4], torque, Fx, Fy ))    

        # stator iron (group==10)
        for ind, stator_xy_complex_number in enumerate(self.stator_xy_complex_data):
            # 1. What we need for iron loss evaluation is the B waveform at a fixed point (x,y). 
            #    For example, (x,y) is the centeroid of element in stator tooth.
            Bx, By = femm.mo_getb( stator_xy_complex_number.real,
                                   stator_xy_complex_number.imag)
            self.stator_Bx_data[ind].append(Bx)
            self.stator_By_data[ind].append(By)

        # rotor iron (group==100)
        for ind, rotor_xy_complex_number in enumerate(self.rotor_xy_complex_data):
            # 2. The element at (x,y) is no longer the same element from last rotor position.
            #    To find the exact element from last rotor position,
            #    we rotate the (x,y) forward as we rotate the model (rotor), get the B value there: (x,y)*rotation_operator, and correct the (Bx,By)/rotation_operator
            new_xy_complex = rotor_xy_complex_number * rotation_operator
            Bx, By = femm.mo_getb( new_xy_complex.real, 
                                   new_xy_complex.imag )
            new_BxBy_complex = (Bx + 1j*By) / rotation_operator
            self.rotor_Bx_data[ind].append(new_BxBy_complex.real)
            self.rotor_By_data[ind].append(new_BxBy_complex.imag)

        femm.mo_close()

    def get_copper_loss_Bolognani(self, stator_slot_area, rotor_slot_area, STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1.0, TEMPERATURE_OF_COIL=75, total_CurrentAmp=None):
        # make sure these two values 
        # space_factor_kCu = SLOT_FILL_FACTOR in Pyrhonen09 design
        # space_factor_kAl = 1 in Pyrhonen09 design
        # as well as temperature
        #   TEMPERATURE_OF_COIL: temperature increase, degrees C
        #   Here, by conductor it means the parallel branch (=1) is considered, and number_of_coil_per_slot = 8 (and will always be 8, see example below).
        #   Just for example, if parallel branch is 2, number_of_coil_per_slot will still be 8, even though we cut the motor in half and count there will be 16 wire cross-sections,
        #   because the reduction in resistance due to parallel branch will be considered into the variable resistance_per_coil.
        #   Remember that the number_of_coil_per_slot (in series) is limited by the high back EMF rather than slot area.
        rho_Copper = (3.76*TEMPERATURE_OF_COIL+873)*1e-9/55. # resistivity

        air_gap_length_delta     = self.im.design_parameters[0]*1e-3 # m

        # http://127.0.0.1:4000/tech/ECCE-2019-Documentation/

        ################################################################
        # Stator Copper Loss 
        ################################################################
        tooth_width_w_t          = self.im.design_parameters[1]*1e-3 # m
        Area_S_slot              = stator_slot_area
        area_copper_S_Cu         = STATOR_SLOT_FILL_FACTOR * Area_S_slot
        a                        = self.im.wily.number_parallel_branch
        zQ                       = self.im.DriveW_zQ
        coil_pitch_yq            = self.im.wily.coil_pitch_y
        Q                        = self.im.Qs
        # the_radius_m             = 1e-3*(0.5*(self.acm_variant.template.d['GP']['mm_r_or'].value + self.im.Length_AirGap + self.im.Radius_InnerStatorYoke))
        stack_length_m           = 1e-3*self.im.stack_length
        # number_of_phase          = 3
        # Ns                       = zQ * self.im.Qs / (2 * number_of_phase * a) # 3 phase winding
        # density_of_copper        = 8960 
        # k_R                      = 1 # AC resistance factor
        current_rms_value        = total_CurrentAmp / 1.4142135623730951 # for one phase
        # Area_conductor_Sc        = Area_S_slot * STATOR_SLOT_FILL_FACTOR / zQ

        Js = (current_rms_value/a) * zQ / area_copper_S_Cu # 逆变器电流current_rms_value在流入电机时，
                                                           # 受到并联支路分流，(current_rms_value/a)才是实际导体中流动的电流值，
                                                           # 这样的电流在一个槽内有zQ个，所以Islot=(current_rms_value/a) * zQ
                                                           # 槽电流除以槽内铜的面积，就是电流密度

        stator_inner_diameter_D = 2*(air_gap_length_delta + self.acm_variant.template.d['GP']['mm_r_or'].value*1e-3)
        slot_height_h_t = 0.5*(self.im.stator_yoke_diameter_Dsyi - stator_inner_diameter_D)
        slot_pitch_pps = np.pi * (stator_inner_diameter_D + slot_height_h_t) / Q
        kov = 1.8 # \in [1.6, 2.0]
        end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + tooth_width_w_t) + slot_pitch_pps*kov * (coil_pitch_yq - 1)

        Vol_Cu = area_copper_S_Cu * (stack_length_m + end_winding_length_Lew) * Q
        stator_copper_loss = rho_Copper * Vol_Cu * Js**2

        Vol_Cu_along_stack = area_copper_S_Cu * (end_winding_length_Lew) * Q
        stator_copper_loss_along_stack = rho_Copper * Vol_Cu_along_stack * Js**2

        print('Stator current [Arms]:', current_rms_value, 'Js:', Js)


        ################################################################
        # Rotor Copper Loss
        ################################################################
        tooth_width_w_t          = self.im.design_parameters[2]*1e-3 # m
        Area_S_slot              = rotor_slot_area
        area_copper_S_Cu         = ROTOR_SLOT_FILL_FACTOR * Area_S_slot
        a                        = 1
        zQ                       = 1
        coil_pitch_yq            = self.im.Qr/self.im.DriveW_poles
        Q                        = self.im.Qr
        # the_radius_m             = 1e-3*(self.acm_variant.template.d['GP']['mm_r_or'].value - self.im.Radius_of_RotorSlot)
        stack_length_m           = 1e-3*self.im.stack_length
        # number_of_phase          = self.im.Qr/self.im.DriveW_poles
        # Ns                       = zQ * self.im.Qr / (2 * number_of_phase * a) # turns in series = 2 for 4 pole winding; = 1 for 2 pole winding
        # density_of_copper        = 8960 
        # k_R                      = 1 # AC resistance factor
        current_rms_value = sum(self.list_rotor_current_amp) / ( 1.4142135623730951 * len(self.list_rotor_current_amp) )
        # Area_conductor_Sc        = Area_S_slot * ROTOR_SLOT_FILL_FACTOR / zQ        
        Jr = (current_rms_value/a) * zQ / area_copper_S_Cu # 逆变器电流current_rms_value在流入电机时，


        rotor_outer_diameter_Dor = 2*(self.acm_variant.template.d['GP']['mm_r_or'].value*1e-3)
        slot_height_h_t = self.im.rotor_slot_height_h_sr
        slot_pitch_pps = np.pi * (rotor_outer_diameter_Dor - slot_height_h_t) / Q
        kov = 1.6 
        end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + tooth_width_w_t) + slot_pitch_pps*kov * (coil_pitch_yq - 1)

        Vol_Cu = area_copper_S_Cu * (stack_length_m + end_winding_length_Lew) * Q
        rotor_copper_loss = rho_Copper * Vol_Cu * Jr**2

        Vol_Cu_along_stack = area_copper_S_Cu * (end_winding_length_Lew) * Q
        rotor_copper_loss_along_stack = rho_Copper * Vol_Cu_along_stack * Jr**2

        print('Rotor current [Arms]:', current_rms_value, 'Jr:', Jr)

        return stator_copper_loss, rotor_copper_loss, stator_copper_loss_along_stack, rotor_copper_loss_along_stack, Js, Jr, Vol_Cu

    def get_copper_loss_pyrhonen(self, stator_slot_area, rotor_slot_area, STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1.0, TEMPERATURE_OF_COIL=75, total_CurrentAmp=None):

        # 这里有一个问题：就是用的几何参数是im_template的值，而不是来自 self.im.design_parameter
        # 这里有一个问题：就是用的几何参数是im_template的值，而不是来自 self.im.design_parameter
        # 这里有一个问题：就是用的几何参数是im_template的值，而不是来自 self.im.design_parameter

        # make sure these two values 
        # space_factor_kCu = SLOT_FILL_FACTOR in Pyrhonen09 design
        # space_factor_kAl = 1 in Pyrhonen09 design
        # as well as temperature
        #   TEMPERATURE_OF_COIL: temperature increase, degrees C
        #   Here, by conductor it means the parallel branch (=1) is considered, and number_of_coil_per_slot = 8 (and will always be 8, see example below).
        #   Just for example, if parallel branch is 2, number_of_coil_per_slot will still be 8, even though we cut the motor in half and count there will be 16 wire cross-sections,
        #   because the reduction in resistance due to parallel branch will be considered into the variable resistance_per_coil.
        #   Remember that the number_of_coil_per_slot (in series) is limited by the high back EMF rather than slot area.
        im = self.im
        rho_Copper = (3.76*TEMPERATURE_OF_COIL+873)*1e-9/55. # resistivity


        # Stator Copper Loss 
        Area_slot                = stator_slot_area
        a                        = im.wily.number_parallel_branch
        zQ                       = im.DriveW_zQ
        coil_pitch_by_slot_count = im.wily.coil_pitch_y
        Q                        = im.Qs
        the_radius_m             = 1e-3*(0.5*(im.Radius_OuterRotor + im.Length_AirGap + im.Radius_InnerStatorYoke))
        stack_length_m           = 1e-3*im.stack_length
        number_of_phase          = 3
        Ns                       = zQ * im.Qs / (2 * number_of_phase * a) # 3 phase winding
        density_of_copper        = 8960 
        k_R                      = 1 # AC resistance factor
        current_rms_value        = total_CurrentAmp / 1.4142135623730951 # for one phase

        Area_conductor_Sc        = Area_slot * STATOR_SLOT_FILL_FACTOR / zQ
        Js = current_rms_value / (a * Area_conductor_Sc)

        coil_span_W = coil_pitch_by_slot_count / Q * (the_radius_m * 2*np.pi)
        lav = 2*stack_length_m + 2.4*coil_span_W + 0.1
        mass_Cu             = density_of_copper * Ns * lav * Area_conductor_Sc
        mass_Cu_along_stack = density_of_copper * Ns * (2*stack_length_m) * Area_conductor_Sc

        stator_copper_loss             = number_of_phase * rho_Copper / density_of_copper * k_R * Js**2 * mass_Cu
        stator_copper_loss_along_stack = number_of_phase * rho_Copper / density_of_copper * k_R * Js**2 * mass_Cu_along_stack # this part of loss will scale as the depth of the 2D model is scaled


        # Rotor Copper Loss
        Area_slot                = rotor_slot_area
        a                        = 1
        zQ                       = 1
        coil_pitch_by_slot_count = im.Qr/im.DriveW_poles
        Q                        = im.Qr
        the_radius_m             = 1e-3*(im.Radius_OuterRotor - im.Radius_of_RotorSlot)
        stack_length_m           = 1e-3*im.stack_length
        number_of_phase          = im.Qr/im.DriveW_poles
        Ns                       = zQ * im.Qr / (2 * number_of_phase * a) # turns in series = 2 for 4 pole winding; = 1 for 2 pole winding
        density_of_copper        = 8960 
        k_R                      = 1 # AC resistance factor
        current_rms_value        = None

        Area_conductor_Sc        = Area_slot * ROTOR_SLOT_FILL_FACTOR / zQ        
        # print('list_rotor_current_amp', self.list_rotor_current_amp) # self.list_rotor_current_amp is defined in population.py
        rotor_copper_loss             = 0.0
        rotor_copper_loss_along_stack = 0.0
        # sum_rotor_current_density     = 0.0
        list_Jr = []
        for amp in self.list_rotor_current_amp:
            current_rms_value              = amp / 1.4142135623730951
            list_Jr.append(current_rms_value / (a * Area_conductor_Sc))
        Jr = sum(list_Jr) / len(list_Jr) # take average for Jr
        # print('Jr=%g Arms/m^2'%(Jr))

        coil_span_W = coil_pitch_by_slot_count / Q * (the_radius_m * 2*np.pi)
        lav = 2*stack_length_m + 2.4*coil_span_W + 0.1
        mass_Cu             = density_of_copper * Ns * lav * Area_conductor_Sc
        mass_Cu_along_stack = density_of_copper * Ns * (2*stack_length_m) * Area_conductor_Sc

        rotor_copper_loss             = number_of_phase * rho_Copper / density_of_copper * k_R * Jr**2 * mass_Cu
        rotor_copper_loss_along_stack = number_of_phase * rho_Copper / density_of_copper * k_R * Jr**2 * mass_Cu_along_stack # this part of loss will scale as the depth of the 2D model is scaled

        # print('stator slot area', stator_slot_area, 'm^2')
        # print('rotor slot area', rotor_slot_area, 'm^2')
        return stator_copper_loss, rotor_copper_loss, stator_copper_loss_along_stack, rotor_copper_loss_along_stack, Js, Jr

    # def compute_iron_loss(self, MAX_FREQUENCY=50e3, SLOT_FILL_FACTOR=0.5, TEMPERATURE_OF_COIL=75):
    #     # http://www.femm.info/wiki/SPMLoss
    #     # % Now, total core loss can be computed in one fell swoop...

    #     # Iron Loss
    #     # % Dividing the result by cs corrects for the lamination stacking factor
    #     if 'M19' in self.acm_variant.template.SI['Steel'] or 'M15' in self.acm_variant.template.SI['Steel']:
    #         # M-19 Steel
    #         ce = 0.530 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
    #         ch = 143.  # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
    #         cs = 0.95  # % Lamination stacking factor (nondimensional)
    #     elif self.acm_variant.template.SI['Steel'] == 'Arnon5':
    #         # %Arnon7
    #         ce = 0.07324 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
    #         ch = 187.6   # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
    #         cs = 0.96    # % Lamination stacking factor (nondimensional)

    #     # % Get parameters for proximity effect loss computation for phase windings
    #     # AWG     = 25       # % Magnet wire gauge used in winding
    #     # dwire   = 0.324861*0.0254*exp(-0.115942*AWG)   # % wire diameter in meters as a function of AWG
    #     # owire   = (58.*1e6) / (1+TEMPERATURE_OF_COIL*0.004) # % conductivity of the wire in S/m at prescribed deltaT
    #     # cePhase = SLOT_FILL_FACTOR * (pi**2/8.) * dwire**2 *owire

    #     # dff = MyLowestHarmonic*thisFrequency*w.*(w<(ns/2));
    #     # try:
    #     #     NFFT = self.number_ans
    #     # except Exception as e:
    #     #     NFFT = len(np.arange(0, 180, self.deg_per_step))
    #     #     raise(e)

    #     if False:
    #         def test(Bx_data, base_freq):
    #             print('There are in total', len(Bx_data), 'elements per step.')
    #             print('There are in total', len(Bx_data[0]), 'steps.')

    #             fig_dft, axes_dft = subplots(2,1)
    #             fig, ax = subplots()
    #             for id_element in range(0, len(Bx_data), 200): # typical id
    #                 Bx_waveform = Bx_data[id_element]
    #                 ax.plot(np.arange(0, 180, self.deg_per_step), Bx_waveform, label=id_element, alpha=0.3)

    #                 # DFT
    #                 samp_freq = 1. / (self.deg_per_step/180.*pi / self.im.Omega)
    #                 utility.basefreqDFT(Bx_waveform, samp_freq, ax_time_domain=axes_dft[0], ax_freq_domain=axes_dft[1], base_freq=base_freq)
    #             ax.legend()
    #         test(self.stator_Bx_data, 500)
    #         test(self.stator_By_data, 500)
    #         test(self.rotor_Bx_data, 1)
    #         test(self.rotor_By_data, 1)
    #         show()

    #     global threshold 
    #     def remove_noises(bxfft, threshold_tuner=0.3): # square window in freq domain
    #         # remove noises in frequency domain to estimate correct power spectrum
    #         # https://dsp.stackexchange.com/questions/9054/removing-noise-from-audio-using-fourier-transform-in-matlab
    #         # https://dsp.stackexchange.com/questions/6220/why-is-it-a-bad-idea-to-filter-by-zeroing-out-fft-bins
    #         # https://www.mathworks.com/help/matlab/math/fourier-transforms.html
    #         global threshold 
    #         threshold = threshold_tuner * np.mean(bxfft)
    #         noises_places = np.where(bxfft<threshold, 0, 1)
    #         bxfft *= noises_places
    #         # print 'threshold', threshold
    #         return bxfft

    #     if self.acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS']<1:
    #         raise Exception('At least one cycle is needed to be solved for calculating machine losses, but number_cycles_in_2ndTSS is', self.acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS'])
    #     NFFT = self.acm_variant.template.fea_config_dict['femm.number_of_steps_2ndTSS'] # length of the flux density component (Bx or By) waveform
    #     samp_freq = 1 / self.step_size_sec
    #     dft_freq = 0.5*samp_freq*np.linspace(0,1,int(NFFT/2+1))
    #     print('[FEMM_SlidingMesh.py] DFT resolution is:', dft_freq[1], 'Hz, and dft_freq is', dft_freq)

    #     stator_eddycurrent_loss_harmonics = np.zeros(len(dft_freq))
    #     stator_hysteresis_loss_harmonics = np.zeros(len(dft_freq))
    #     stator_volume = 0.0
    #     rotor_eddycurrent_loss_harmonics = np.zeros(len(dft_freq))
    #     rotor_hysteresis_loss_harmonics = np.zeros(len(dft_freq))
    #     rotor_volume = 0.0
    #     magnet_eddycurrent_loss_harmonics = np.zeros(len(dft_freq))
    #     magnet_volume = 0.0
    #     coil_proximity_loss_harmonics = np.zeros(len(dft_freq))

    #     plt.figure()
    #     plt.plot(self.list_of_femm_triangle_element[0].femm_Bx_list, self.list_of_femm_triangle_element[0].femm_By_list, color='k')
    #     plt.plot(self.list_of_femm_triangle_element[100].femm_Bx_list, self.list_of_femm_triangle_element[100].femm_By_list, color='r')
    #     plt.plot(self.list_of_femm_triangle_element[1000].femm_Bx_list, self.list_of_femm_triangle_element[1000].femm_By_list, color='b')
    #     plt.plot(self.list_of_femm_triangle_element[10000].femm_Bx_list, self.list_of_femm_triangle_element[10000].femm_By_list, color='g')
    #     plt.show()

    #     for el in self.list_of_femm_triangle_element:

    #         bxfft = utility.singleSidedDFT(el.femm_Bx_list, samp_freq)
    #         byfft = utility.singleSidedDFT(el.femm_By_list, samp_freq)

    #             # # test remove noises in spectrum
    #             # index_enough = -1
    #             # fig_dft, axes_dft = plt.subplots(4,1, sharex=True)
    #             # axes_dft[3].plot(dft_freq[:index_enough], bxfft[:index_enough], '>',alpha=0.4)
    #             # _bxfft = remove_noises(bxfft)
    #             # _byfft = remove_noises(byfft)
    #             # # test remove noises in spectrum
    #             # axes_dft[0].plot(dft_freq[:index_enough], stator_eddycurrent_loss_harmonics[:index_enough], '+',alpha=0.4)
    #             # axes_dft[0].set_xlabel('Frequency [Hz]')
    #             # axes_dft[0].set_ylabel('Stator\nEddy\nCurrent\nLoss [W]')
    #             # axes_dft[1].plot(dft_freq[:index_enough], stator_hysteresis_loss_harmonics[:index_enough], '^',alpha=0.4)
    #             # axes_dft[1].set_xlabel('Frequency [Hz]')
    #             # axes_dft[1].set_ylabel('Stator\nHysteresis\nLoss [W]')
    #             # axes_dft[2].plot(dft_freq[:index_enough], _bxfft[:index_enough], '>',alpha=0.4)
    #             # axes_dft[2].plot(dft_freq[:index_enough], np.ones(len(dft_freq[:index_enough]))*threshold)
    #             # axes_dft[3].plot(dft_freq[:index_enough], np.ones(len(dft_freq[:index_enough]))*threshold)
    #             # plt.show()

    #         # bxfft = remove_noises(bxfft)
    #         # byfft = remove_noises(byfft)
    #         bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density

    #         volume_element = (el.area*1e-6) * (self.acm_variant.template.d['EX']['mm_template_stack_length']*1e-3) # Compute the volume of each element in units of meter^3

    #         # iron loss (stator)
    #         if el.group == self.GroupSummary['stator_iron_core']:
    #             stator_eddycurrent_loss_harmonics += (ce*dft_freq**2 * volume_element/cs ) * bsq
    #             stator_hysteresis_loss_harmonics  += (ch*dft_freq    * volume_element/cs ) * bsq
    #             stator_volume += volume_element

    #         # iron loss (rotor)
    #         if el.group == self.GroupSummary['rotor_iron_core']:
    #             rotor_eddycurrent_loss_harmonics += (ce*dft_freq**2 * volume_element/cs ) * bsq
    #             rotor_hysteresis_loss_harmonics  += (ch*dft_freq    * volume_element/cs ) * bsq
    #             rotor_volume += volume_element

    #         # proximity loss
    #         # if el.group == self.GroupSummary['coils']:
    #         #     coil_proximity_loss_harmonics  += 0

    #         # magnet (eddy current) loss
    #         # if el.group == self.GroupSummary['magnet']:
    #         #     magnet_eddycurrent_loss_harmonics  += 0

    #     print('[FEMM_SlidingMesh.py] DEBUG signal lengths:', NFFT, len(dft_freq), len(bxfft), len(el.femm_Bx_list))

    #     try:
    #         index_enough = next(ind for ind, el in enumerate(dft_freq) if el > MAX_FREQUENCY) # 50e3 Hz is enough 10 times base frequency 500 Hz
    #     except StopIteration:
    #         index_enough = None

    #     stator_eddycurrent_loss = sum( stator_eddycurrent_loss_harmonics[:index_enough] ) 
    #     stator_hysteresis_loss  = sum( stator_hysteresis_loss_harmonics[:index_enough] )
    #     rotor_eddycurrent_loss  = sum( rotor_eddycurrent_loss_harmonics[:index_enough] ) 
    #     rotor_hysteresis_loss   = sum( rotor_hysteresis_loss_harmonics[:index_enough]  )

    #     print(  stator_eddycurrent_loss,
    #             stator_hysteresis_loss,
    #             stator_volume,
    #             rotor_eddycurrent_loss,
    #             rotor_hysteresis_loss,
    #             rotor_volume)

    #     '''
    #     # find the index of dft_freq that corresponds to 50e3 Hz, because higher results are contaminated by noises in fourier analysis

    #     fig_dft, axes_dft = plt.subplots(4,1, sharex=True)
    #     axes_dft[0].plot(dft_freq[:index_enough], stator_eddycurrent_loss_harmonics[:index_enough], '+',alpha=0.4)
    #     axes_dft[0].set_xlabel('Frequency [Hz]')
    #     axes_dft[0].set_ylabel('Stator\nEddy\nCurrent\nLoss [W]')
    #     axes_dft[1].plot(dft_freq[:index_enough], stator_hysteresis_loss_harmonics[:index_enough], '^',alpha=0.4)
    #     axes_dft[1].set_xlabel('Frequency [Hz]')
    #     axes_dft[1].set_ylabel('Stator\nHysteresis\nLoss [W]')

    #     # Rotor iron loss
    #     rotor_eddycurrent_loss_harmonics = np.zeros(len(dft_freq))
    #     rotor_hysteresis_loss_harmonics = np.zeros(len(dft_freq))
    #     for id_element, area_element in enumerate(self.rotor_Area_data):
    #         rotor_Bx_waveform = self.rotor_Bx_data[id_element]
    #         rotor_By_waveform = self.rotor_By_data[id_element]
    #         bxfft = utility.singleSidedDFT(rotor_Bx_waveform, samp_freq)
    #         byfft = utility.singleSidedDFT(rotor_By_waveform, samp_freq)

    #         bxfft = remove_noises(bxfft)
    #         byfft = remove_noises(byfft)

    #         bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
    #         volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3

    #         rotor_eddycurrent_loss_harmonics += (ce*dft_freq**2 * volume_element/cs ) * bsq
    #         rotor_hysteresis_loss_harmonics  += (ch*dft_freq    * volume_element/cs ) * bsq
    #         rotor_volume += volume_element
    #     rotor_eddycurrent_loss = sum( rotor_eddycurrent_loss_harmonics[:index_enough] ) 
    #     rotor_hysteresis_loss  = sum( rotor_hysteresis_loss_harmonics[:index_enough]  )

    #     axes_dft[2].plot(dft_freq[:index_enough], rotor_eddycurrent_loss_harmonics[:index_enough], '+',alpha=0.4)
    #     axes_dft[2].set_xlabel('Frequency [Hz]')
    #     axes_dft[2].set_ylabel('\nRotor\nEddy Current\nLoss [W]')
    #     axes_dft[3].plot(dft_freq[:index_enough], rotor_hysteresis_loss_harmonics[:index_enough], '^',alpha=0.4)
    #     axes_dft[3].set_xlabel('Frequency [Hz]')
    #     axes_dft[3].set_ylabel('\nRotor\nHysteresis\nLoss [W]')
    #     '''
    #         # Copper Loss - Stator Winding Proximity Effect (Rotor side is neglected because the slip frequency is low)
    #         # this should be done with g==2, i.e., field data on coil area
    #         # % and prox losses can be totalled up in a similar way as iron loss
    #         # prox_loss += np.dot(cePhase * dft_freq**2, bsq) * volume_element

    #     # Did you use 1/4 model for the loss calculation?
    #     number_of_fraction = 1
    #     return ( number_of_fraction*stator_eddycurrent_loss, 
    #              number_of_fraction*stator_hysteresis_loss, 
    #              number_of_fraction*rotor_eddycurrent_loss, 
    #              number_of_fraction*rotor_hysteresis_loss, 
    #              number_of_fraction*stator_volume, 
    #              number_of_fraction*rotor_volume )
