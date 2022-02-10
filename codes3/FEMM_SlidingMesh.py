#coding:utf-8

from this import d
import femm
import numpy as np
import logging
import os
from collections import OrderedDict
import utility
from time import sleep
from time import time as clock_time
import winding_layout

SELECT_ALL = 4
EPS = 1e-2 # unit mm

class Individual_Analyzer_FEMM_Edition(object):
    def __init__(self, p):

        self.p = p

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

    def add(self, time, rotor_position_mech_deg, torque, forces, energy, circuitProperties):

        self.femm_time . append( time )
        self.femm_rotor_position_mech_deg. append( rotor_position_mech_deg )

        self.femm_torque . append( torque )
        self.femm_forces . append( forces )
        self.femm_energy . append( energy )
        # print('-'*100)
        # print(circuitProperties)
        # print(circuitProperties[0])
        # print(circuitProperties[0][0])
        # print(circuitProperties[0][0][0])
        circuit_currents     = list(map(lambda el: el[0], circuitProperties))
        circuit_voltages     = list(map(lambda el: el[1], circuitProperties))
        circuit_fluxLinkages = list(map(lambda el: el[2], circuitProperties))
        self.femm_circuit_currents . append( circuit_currents ) # 6 coils UVW and AC/BD
        self.femm_circuit_voltages . append( circuit_voltages ) # 6 coils UVW and AC/BD
        self.femm_circuit_fluxLinkages . append( circuit_fluxLinkages ) # 6 coils UVW and AC/BD

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
            bearing_current_U = np.array(circuit_current_coil_U_GrpA) - np.array(circuit_current_coil_U_GrpB)
            bearing_current_V = np.array(circuit_current_coil_V_GrpA) - np.array(circuit_current_coil_V_GrpB)
            bearing_current_W = np.array(circuit_current_coil_W_GrpA) - np.array(circuit_current_coil_W_GrpB)
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
        self.femm_motor_currents_d. append( motor_current_d )
        self.femm_motor_currents_q. append( motor_current_q )
        self.femm_bearing_currents_d. append( bearing_current_d )
        self.femm_bearing_currents_q. append( bearing_current_q )
        # Fluxes
        motor_fluxLinkage_d, motor_fluxLinkage_q     = Park_transform(motor_fluxLinkage_alpha, motor_fluxLinkage_beta, theta=rotor_position_elec_rad)
        bearing_fluxLinkage_d, bearing_fluxLinkage_q = Park_transform(bearing_fluxLinkage_alpha, bearing_fluxLinkage_beta, theta=rotor_position_elec_rad)
        self.femm_motor_fluxLinkage_d. append( motor_fluxLinkage_d )
        self.femm_motor_fluxLinkage_q. append( motor_fluxLinkage_q )
        self.femm_bearing_fluxLinkage_d. append( bearing_fluxLinkage_d )
        self.femm_bearing_fluxLinkage_q. append( bearing_fluxLinkage_q )

    def get_ss_data(self):
        # Torque
        torque_average = np.mean(self.femm_torque)
        torque_error = np.array(self.femm_torque) - torque_average
        ss_max_torque_error = max(torque_error), min(torque_error)
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
        print('\t normalized_torque_ripple =', normalized_torque_ripple, 'Nm')
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
        spec_performance_dict['normalized_torque_ripple'] = normalized_torque_ripple
        spec_performance_dict['normalized_force_error_magnitude'] = np.max(sfv.ss_max_force_err_abs) / sfv.ss_avg_force_magnitude
        spec_performance_dict['force_error_angle'] = np.max(sfv.ss_max_force_err_ang)
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

class FEMM_SlidingMesh(object):

    def __init__(self, acm_variant):
        self.acm_variant = acm_variant
        # self.vangogh = FEMM_Solver.VanGogh_FEMM(acm_variant)

        self.deg_per_step = acm_variant.template.fea_config_dict['femm.deg_per_step'] # deg, we need this for show_results
        self.dir_parent    = acm_variant.template.fea_config_dict['dir.parent']

        self.rotor_phase_name_list = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    def open(self):
        femm.openfemm(True) # bHide # False for debug
        femm.newdocument(0) # magnetic

    def probdef(self, stack_length=100, time_harmonic_study_frequency=0):
        self.time_harmonic_study_frequency = time_harmonic_study_frequency
        # femm.smartmesh(False) <- This will not work due to bug of femm.__init__  # let mi_smartmesh deside. You must turn it off in parasolver.py
        femm.callfemm_noeval('smartmesh(0)') # call this before probdef?
        femm.mi_probdef(time_harmonic_study_frequency, 'millimeters', 'planar', 1e-8, # must < 1e-8
                        stack_length, 18, 1) # The acsolver parameter (default: 0) specifies which solver is to be used for AC problems: 0 for successive approximation, 1 for Newton.
                                             # 1 for 'I intend to try the acsolver of Newton, as this is the default for JMAG@[Nonlinear Calculation] Setting Panel in the [Study Properties] Dialog Box'
        # femm.callfemm_noeval('mi_smartmesh(0)') # call this after probdef?
        self.bool_automesh = False # setting to false gives no effect?

        # femm.smartmesh(True) # let mi_smartmesh deside. You must turn it off in parasolver.py
        # self.bool_automesh = True # setting to false gives no effect?

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
            hdata, bdata = np.loadtxt(self.dir_parent + 'BH/M-19-Steel-BH-Curve-afterJMAGsmooth.BH', unpack=True, usecols=(0,1))
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
        print('[FEMM_SllidingMesh.py] DEBUG ampD, ampB:', ampD, ampB)
        print('[FEMM_SllidingMesh.py] circuit turns for combined winding:', turnsDPNV)

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
        GP = self.acm_variant.template.d['GP']
        EX = self.acm_variant.template.d['EX']
        p  = self.acm_variant.template.SI['p']

        # figure out initial rotor angle
        deg_pole_span = 180/p
        initial_rotor_position_mech_deg = (deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + EX['wily'].deg_winding_U_phase_phase_axis_angle + deg_pole_span
        print('[FEMM_SlidingMesh.py] initial_rotor_position_mech_deg:', initial_rotor_position_mech_deg, 'deg || deg_winding_U_phase_phase_axis_angle:', '', EX['wily'].deg_winding_U_phase_phase_axis_angle)

        # figure out step size
        electrical_period = self.acm_variant.template.fea_config_dict['femm.number_cycles_in_2ndTSS']/EX['DriveW_Freq']
        number_of_steps   = self.acm_variant.template.fea_config_dict['femm.number_of_steps_2ndTSS']
        step_size_sec = electrical_period / number_of_steps
        step_size_mech_deg = EX['Omega'] * step_size_sec / np.pi * 180

        # transient FEA with sliding band (air gap boundary in FEMM)
        print('[FEMM_SlidingMesh.py] Run FEMM solver...')
        self.acm_variant.analyzer = Individual_Analyzer_FEMM_Edition(p=p)
        for index in range(number_of_steps):
            time                         = index * step_size_sec
            RotorAngle_MechanicalDegrees = index * step_size_mech_deg
            current_sources = self.update_circuit_excitation(time)
            femm.mi_modifyboundprop('WholeModelSlidingBand', 10, initial_rotor_position_mech_deg+RotorAngle_MechanicalDegrees+self.acm_variant.template.fea_config_dict['femm.MechDeg_IdEqualToNonZeroAngle']) # change inner angle to 

            # if index > 0:
            #     femm.mi_setprevious(self.project_file_name[:-4] + f'-{index-1:03d}.ans', 1) # 1: must be incremental permeability (because frozen permability means the permeability is fixed).
            # Starting nonlinear field from previous solution is not possible with current FEMM's function
            femm.mi_saveas(self.project_file_name[:-4] + f'-{index:03d}.fem')

            # femm.mi_createmesh() # no need
            femm.mi_analyze(1)

            # collect results at: index
            femm.mi_loadsolution()
            torque = femm.mo_gapintegral('WholeModelSlidingBand',0)
            forces = femm.mo_gapintegral('WholeModelSlidingBand',1)
            energy = femm.mo_gapintegral('WholeModelSlidingBand',2)
            circuitProperties = ( femm.mo_getcircuitproperties('U-GrpAC'),
                                  femm.mo_getcircuitproperties('V-GrpAC'),
                                  femm.mo_getcircuitproperties('W-GrpAC'),
                                  femm.mo_getcircuitproperties('U-GrpBD'),
                                  femm.mo_getcircuitproperties('V-GrpBD'),
                                  femm.mo_getcircuitproperties('W-GrpBD') )
            self.acm_variant.analyzer.add(time, RotorAngle_MechanicalDegrees, torque, forces, energy, circuitProperties)
            print(f'\t {index}, {time*1000} ms, {RotorAngle_MechanicalDegrees} deg : {torque:.2f} Nm, [{forces[0]:.2f} N, {forces[1]:.2f}] N, {energy:.2f} J', end=' | ')
            print(' A, '.join(f'{current:.1f}' for current in current_sources), 'A')
            for el in circuitProperties:
                print('\t\t', el)
            femm.mo_close()

            # save field results for incremental study
            # femm.mi_saveas('ptemp-acmop.fem')

        femm.mi_close()

    def compute_objectives(self):
        pass
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

    # def save(self):
        # print('[FEMM_SlidingMesh.py] save to:', self.output_file_name + '.fem')
        # femm.mi_saveas(self.output_file_name + '.fem')
        # self.model_rotor_rotate(time)
        # femm.mi_saveas(fem_file)

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
        CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * EX['fill_factor'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
        CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
        CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。

        # Maybe there is a bug here... regarding the excitation for suspension winding...
        variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
        variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
        EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
        EX['DriveW_CurrentAmp'] = acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
        EX['BeariW_CurrentAmp'] = acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

        slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
        print('[inner_rotor_motor.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

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
        femm.mi_drawarc(p1[0],p1[1],p2[0],p2[1],angle/pi*180,maxseg) # [deg]

    @staticmethod
    def add_arc(p1, p2, angle, maxseg=1, center=None, **kwarg):
        femm.mi_addarc(p1[0],p1[1],p2[0],p2[1],angle/pi*180,maxseg) # [deg]

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
            THETA = i / 180.0 * pi
            X = R*cos(THETA)
            Y = R*sin(THETA)
            B_vector_complex = femm.mo_getb(X, Y)
            B_X_complex = B_vector_complex[0]
            B_Y_complex = B_vector_complex[1]
            B_X_real = np.real(B_X_complex)
            B_Y_real = np.real(B_Y_complex)
            # Assume the magnitude is all due to radial component
            B_magitude = sqrt(B_X_real**2 + B_Y_real**2)
            inner_product = B_X_real * X + B_Y_real *Y
            list_B_magitude.append( B_magitude * copysign(1, inner_product) )
        return list_B_magitude

    def has_results(self, dir_run=None):
        # print 'self.freq', self.freq
        if dir_run == None:
            dir_run = self.dir_run

        a = [f for f in os.listdir(dir_run) if '.ans' in f].__len__()
        b = [f for f in os.listdir(dir_run) if '.fem' in f].__len__()
        if a == 0:
            if 'no_duplicates.txt' in os.listdir(dir_run):
                return True # 直接把femm的结果从服务器上拷过来也行
            else:
                return False
        print('[FEMM.has_results] ans count: %d. fem count: %d.' % (a, b))
        return a == b

    def show_results(self, bool_plot=True):

        if self.flag_eddycurrent_solver == True:
            print('show results for eddy current solver')
            return self.show_results_eddycurrent(bool_plot=bool_plot)

        if self.flag_static_solver == True:
            print('show results for static solver')
            return self.show_results_static(bool_plot=bool_plot)

        return None

    def show_results_eddycurrent(self, bool_plot):
        if self.fraction == 1:
            self.fraction = 2 # this fixes a bug calling femm_integrate_4_current without calling run_frequency_sweeping before
            raise Exception('You should initialize FEMM Solver with freq!=0')
    
        ans_file_list = os.listdir(self.dir_run_sweeping)
        ans_file_list = [f for f in ans_file_list if '.ans' in f]

        femm.openfemm(True)
        list_ans_file = []
        list_torque = []
        # write results to a data file
        str_results = ''
        for ind, f in enumerate( ans_file_list ):
            # femm.opendocument(self.dir_run_sweeping + f[:-4] + '.fem')
            # femm.mi_loadsolution()
            femm.opendocument(self.dir_run_sweeping + f)

            # physical amount on rotor
            femm.mo_groupselectblock(100)
            femm.mo_groupselectblock(101)
            Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
            Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
            torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
            femm.mo_clearblock()

            # rotor current
            # _ = self.femm_integrate_4_current()
            # # TODO: this is for testing phase
            # if float(f[3:-6]) == 3:
            #     print '\n', f[3:-6], 'Hz'
            #     for el in vals_results_rotor_current:
            #         print abs(el)

            str_results += "%s %g %g %g\n" % ( f[3:-6], torque, Fx, Fy ) 
            list_ans_file.append(f)
            list_torque.append(torque)
            femm.mo_close() 
        with open(self.dir_run_sweeping + "eddycurrent_results.txt", "w") as stream:
            stream.write(str_results)

        # find breakdown torque and slip frequency that we are interested
        index, breakdown_torque = utility.get_max_and_index(list_torque)
        slip_freq_breakdown_torque = list_ans_file[index][3:-6]
        print("FEMM's breakdown data: %s Hz, %g Nm" % (slip_freq_breakdown_torque, breakdown_torque))
        self.im.update_mechanical_parameters(float(slip_freq_breakdown_torque))
        # if self.im.slip_freq_breakdown_torque != float(slip_freq_breakdown_torque):
        #     raise Exception('[DEBUG] JMAG disagrees with FEMM (!= 3 Hz).')

        # write rotor currents to file
        self.femm_integrate_4_current(self.dir_run_sweeping + list_ans_file[index], self.fraction)

        femm.closefemm()

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
            angle_per_slot = 2*pi/im.Qr
            THETA_BAR = pi - angle_per_slot + EPS # add EPS for the half bar
            # print 'number of rotor_slot per partial model', self.rotor_slot_per_pole * int(4/fraction)
            for i in range(self.rotor_slot_per_pole * int(4/fraction)):
                THETA_BAR += angle_per_slot
                THETA = THETA_BAR
                X = R*cos(THETA); Y = R*sin(THETA)
                femm.mo_selectblock(X, Y) # or you can select circuit rA rB ...
                vals_results_rotor_current.append(femm.mo_blockintegral(7)) # integrate for current
                femm.mo_clearblock()
            # the other half bar of rA
            THETA_BAR += angle_per_slot
            THETA = THETA_BAR - 2*EPS
            X = R*cos(THETA); Y = R*sin(THETA)
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


    def show_results_static(self, bool_plot=True):
        # Collect results from all the .ans file in the dir_run folder of FEMM.

        # recall that for static FEA, you call show_results once when half .ans files are generated from watchdog
        self.freq = 0 # needed for 

        # TODO 判断，如果文件存在，且不是空的！
        # if exists .txt file, then load it
        missed_ans_file_list = []
        if os.path.exists(self.dir_run + "static_results.txt"):
            data = np.loadtxt(self.dir_run + "static_results.txt", unpack=True, usecols=(0,1,2,3))

            # use dict to eliminate duplicates
            results_dict = {}
            for i in range(len(data[0])):
                results_dict[data[0][i]] = (data[1][i], data[2][i], data[3][i]) 
            keys_without_duplicates = list(OrderedDict.fromkeys(data[0])) # remove duplicated item because it is a dict now!
            keys_without_duplicates.sort()

            # check for missed .ans files
            if len( np.arange(0, 180, self.deg_per_step) ) == len(keys_without_duplicates):
                pass
            else:
                for self.rotor_position_in_deg in np.arange(0, 180, self.deg_per_step):
                    flag_missed = True
                    for key in keys_without_duplicates:
                        if int('%04d'%(10*self.rotor_position_in_deg)) == key:
                            flag_missed = False
                            break
                    if flag_missed == True:
                        missed_ans_file_list.append( self.get_output_file_name(booL_dir=False) + '%04d'%(10*self.rotor_position_in_deg) + '.ans' )
            print('missed:', missed_ans_file_list)

            # typical print gives: 5 1813 1795 1799.0
            # print len(missed_ans_file_list), len(data[0]), len(keys_without_duplicates), keys_without_duplicates[-1]
            # quit()

            # write data without duplicates to file
            with open(self.dir_run + "static_results_no_duplicates.txt", 'w') as f:
                for key in keys_without_duplicates:
                    f.writelines('%g %g %g %g\n' % (key, results_dict[key][0], results_dict[key][1], results_dict[key][2]))
                print('[FEMM.show_results_static] the last key is', max(keys_without_duplicates), '[begin from 0]. the length of keys is', len(keys_without_duplicates))

            data = np.loadtxt(self.dir_run + "static_results_no_duplicates.txt", unpack=True, usecols=(0,1,2,3))

            last_index = len(data[0])
        else:
            last_index = 0

        ans_file_list = os.listdir(self.dir_run)
        ans_file_list = [f for f in ans_file_list if '.ans' in f]

        # last_index == 0 则说明是第一次运行后处理
        if last_index > 0:
            if len(ans_file_list) <= last_index:
                if bool_plot == True:
                    self.plot_results(data)
                return data
            else:
                print('There are new .ans files. Now append them')

        # iter ans_file_list and missed_ans_file_list, and write to .txt 
        femm.openfemm(True) # bHide
        print('there are total %d .ans files'%(len(ans_file_list)))
        print('I am going to append the rest %d ones.'%(len(ans_file_list) - last_index))
        for ind, f in enumerate( ans_file_list[last_index:] + missed_ans_file_list ):
            if ind >= len(ans_file_list[last_index:]):
                print('...open missed .ans files')
                if os.path.exists(self.dir_run + f) == False:
                    print('run mi_analyze for %s' % (f))
                    femm.opendocument(self.dir_run + f[:-4] + '.fem')
                    femm.mi_analyze(1)
                else:
                    femm.opendocument(self.dir_run + f[:-4] + '.fem')                    
            else:
                print(last_index + ind, end=' ')
                femm.opendocument(self.dir_run + f[:-4] + '.fem')

            # load solution (if corrupted, re-run)
            try:
                femm.mi_loadsolution()
            except Exception as e:
                logger = logging.getLogger(__name__)
                logger.error('The .ans file to this .fem file is corrupted. re-run the .fem file %s'%(f), exc_info=True)
                femm.opendocument(self.dir_run + f[:-4] + '.fem')
                femm.mi_analyze(1)
                femm.mi_loadsolution()

            # get the physical amounts on the rotor
            try:
                femm.mo_groupselectblock(100)
                femm.mo_groupselectblock(101)
                Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
                Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
                torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
                femm.mo_clearblock()

                # Air Gap Boundary for Rotor Motion #5
                # gap_torque = femm.mo_gapintegral("AGB4RM", 0)
                # gap_force = femm.mo_gapintegral("AGB4RM", 1)
                # print gap_force, gap_torque, torque, Fx, Fy

                # write results to a data file
                with open(self.dir_run + "static_results.txt", "a") as stream:
                    stream.write("%s %g %g %g\n" % ( f[-8:-4], torque, Fx, Fy ))
            except Exception as e:
                logger = logging.getLogger(__name__)
                logger.error('Encounter error while post-processing (integrating, etc.).', exc_info=True)
                raise e

            # avoid run out of RAM when there are a thousand of ans files loaded into femm...
            # if ind % 10 == 0:
            #     femm.closefemm() 
            #     femm.openfemm(True)
            femm.mo_close()  # use mo_ to close .ans file
            femm.mi_close()  # use mi_ to close .fem file

        print('done. append to static_results.txt.')
        femm.closefemm()

        try:
            data
        except:
            print('call this method again to plot...')
            return None
        if bool_plot == True:
            self.plot_results(data)
        return data

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

        coil_span_W = coil_pitch_by_slot_count / Q * (the_radius_m * 2*pi)
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

        coil_span_W = coil_pitch_by_slot_count / Q * (the_radius_m * 2*pi)
        lav = 2*stack_length_m + 2.4*coil_span_W + 0.1
        mass_Cu             = density_of_copper * Ns * lav * Area_conductor_Sc
        mass_Cu_along_stack = density_of_copper * Ns * (2*stack_length_m) * Area_conductor_Sc

        rotor_copper_loss             = number_of_phase * rho_Copper / density_of_copper * k_R * Jr**2 * mass_Cu
        rotor_copper_loss_along_stack = number_of_phase * rho_Copper / density_of_copper * k_R * Jr**2 * mass_Cu_along_stack # this part of loss will scale as the depth of the 2D model is scaled

        # print('stator slot area', stator_slot_area, 'm^2')
        # print('rotor slot area', rotor_slot_area, 'm^2')
        return stator_copper_loss, rotor_copper_loss, stator_copper_loss_along_stack, rotor_copper_loss_along_stack, Js, Jr

    def get_iron_loss(self, MAX_FREQUENCY=50e3, SLOT_FILL_FACTOR=0.5, TEMPERATURE_OF_COIL=75):
        # http://www.femm.info/wiki/SPMLoss
        # % Now, total core loss can be computed in one fell swoop...
        im = self.im
        samp_freq = 1. / (self.deg_per_step/180.*pi / self.im.Omega)
        print('Sampling frequency is', samp_freq, 'Hz', '<=> deg_per_step is', self.deg_per_step)

        # Iron Loss
        # % Dividing the result by cs corrects for the lamination stacking factor
        if 'M19' in self.im.spec_input_dict['Steel'] or 'M15' in self.im.spec_input_dict['Steel']:
            # M-19 Steel
            ce = 0.530 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
            ch = 143.  # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
            cs = 0.95  # % Lamination stacking factor (nondimensional)
        elif self.im.spec_input_dict['Steel'] == 'Arnon5':
            # %Arnon7
            ce = 0.07324 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
            ch = 187.6   # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
            cs = 0.96    # % Lamination stacking factor (nondimensional)

        # % Get parameters for proximity effect loss computation for phase windings
        # AWG     = 25       # % Magnet wire gauge used in winding
        # dwire   = 0.324861*0.0254*exp(-0.115942*AWG)   # % wire diameter in meters as a function of AWG
        # owire   = (58.*1e6) / (1+TEMPERATURE_OF_COIL*0.004) # % conductivity of the wire in S/m at prescribed deltaT
        # cePhase = SLOT_FILL_FACTOR * (pi**2/8.) * dwire**2 *owire

        # dff = MyLowestHarmonic*thisFrequency*w.*(w<(ns/2));  
        try:
            NFFT = self.number_ans
        except:
            NFFT = len(np.arange(0, 180, self.deg_per_step))

        dft_freq = 0.5*samp_freq*np.linspace(0,1,NFFT/2+1)
        print(dft_freq[1])
        dft_freq = dft_freq + im.slip_freq_breakdown_torque
        print(dft_freq[1])

        stator_eddycurrent_loss = 0.0
        stator_hysteresis_loss = 0.0
        rotor_eddycurrent_loss = 0.0
        rotor_hysteresis_loss = 0.0
        stator_volume = 0.0
        rotor_volume = 0.0
        # prox_loss = 0.0

        # DEBUG:
        # for id_element in range(self.number_ans):
        # for id_element in range(self.number_of_elements):

        # Element-wise calculation
        try:
            self.stator_Area_data
        except:
            # def what_loadtxt_should_be_doing(fname):
            #     data = []
            #     with open(fname, 'r') as f:
            #         read_iterator = csv_reader(f, delimiter =' ')
            #         for ind, row in enumerate(self.whole_row_reader(read_iterator)):
            #             print row
            #             data.append([float(el) for el in row])
            #     return data
            # self.stator_Bx_data   = what_loadtxt_should_be_doing(self.dir_run + 'stator_Bx_data.txt')
            # self.stator_By_data   = what_loadtxt_should_be_doing(self.dir_run + 'stator_By_data.txt')
            # self.stator_Area_data = what_loadtxt_should_be_doing(self.dir_run + 'stator_Area_data.txt')
            # self.rotor_Bx_data    = what_loadtxt_should_be_doing(self.dir_run + 'rotor_Bx_data.txt')
            # self.rotor_By_data    = what_loadtxt_should_be_doing(self.dir_run + 'rotor_By_data.txt')
            # self.rotor_Area_data  = what_loadtxt_should_be_doing(self.dir_run + 'rotor_Area_data.txt')

            self.stator_Bx_data   = np.loadtxt(self.dir_run + 'stator_Bx_data.txt'  )
            self.stator_By_data   = np.loadtxt(self.dir_run + 'stator_By_data.txt'  )
            self.stator_Area_data = np.loadtxt(self.dir_run + 'stator_Area_data.txt')
            self.rotor_Bx_data    = np.loadtxt(self.dir_run + 'rotor_Bx_data.txt'   )
            self.rotor_By_data    = np.loadtxt(self.dir_run + 'rotor_By_data.txt'   )
            self.rotor_Area_data  = np.loadtxt(self.dir_run + 'rotor_Area_data.txt' )
    
        # print np.shape(self.stator_Bx_data)
        # print np.shape(self.stator_Bx_data[0])
        # quit()

        from pylab import subplots, show
        if False:
            def test(Bx_data, base_freq):
                print('There are in total', len(Bx_data), 'elements per step.')
                print('There are in total', len(Bx_data[0]), 'steps.')

                fig_dft, axes_dft = subplots(2,1)
                fig, ax = subplots()
                for id_element in range(0, len(Bx_data), 200): # typical id
                    Bx_waveform = Bx_data[id_element]
                    ax.plot(np.arange(0, 180, self.deg_per_step), Bx_waveform, label=id_element, alpha=0.3)

                    # DFT
                    samp_freq = 1. / (self.deg_per_step/180.*pi / self.im.Omega)
                    utility.basefreqDFT(Bx_waveform, samp_freq, ax_time_domain=axes_dft[0], ax_freq_domain=axes_dft[1], base_freq=base_freq)
                ax.legend()
            test(self.stator_Bx_data, 500)
            test(self.stator_By_data, 500)
            test(self.rotor_Bx_data, 1)
            test(self.rotor_By_data, 1)
            show()


        # remove noises in frequency domain to estimate correct power spectrum
        # https://dsp.stackexchange.com/questions/9054/removing-noise-from-audio-using-fourier-transform-in-matlab
        # https://dsp.stackexchange.com/questions/6220/why-is-it-a-bad-idea-to-filter-by-zeroing-out-fft-bins
        # https://www.mathworks.com/help/matlab/math/fourier-transforms.html

        global threshold 
        def remove_noises(bxfft, threshold_tunner=0.3): # square window in freq domain
            global threshold 
            threshold = threshold_tunner * np.mean(bxfft)
            noises_places = np.where(bxfft<threshold, 0, 1)
            bxfft *= noises_places
            # print 'threshold', threshold
            return bxfft

        # Stator iron loss
        stator_eddycurrent_loss_dft = np.zeros(len(dft_freq))
        stator_hysteresis_loss_dft = np.zeros(len(dft_freq))
        for id_element, area_element in enumerate(self.stator_Area_data):
            stator_Bx_waveform = self.stator_Bx_data[id_element]
            stator_By_waveform = self.stator_By_data[id_element]
            bxfft = utility.singleSidedDFT(stator_Bx_waveform, samp_freq)
            byfft = utility.singleSidedDFT(stator_By_waveform, samp_freq)

            # # test remove noises in spectrum
            # index_enough = -1
            # fig_dft, axes_dft = subplots(4,1, sharex=True)
            # axes_dft[3].plot(dft_freq[:index_enough], bxfft[:index_enough], '>',alpha=0.4)

            bxfft = remove_noises(bxfft)
            byfft = remove_noises(byfft)

            bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
            volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3            

            stator_eddycurrent_loss_dft += (ce*dft_freq**2 * volume_element/cs ) * bsq
            stator_hysteresis_loss_dft  += (ch*dft_freq    * volume_element/cs ) * bsq
            stator_volume += volume_element

            # # test remove noises in spectrum
            # axes_dft[0].plot(dft_freq[:index_enough], stator_eddycurrent_loss_dft[:index_enough], '+',alpha=0.4)
            # axes_dft[0].set_xlabel('Frequency [Hz]')
            # axes_dft[0].set_ylabel('Stator\nEddy\nCurrent\nLoss [W]')
            # axes_dft[1].plot(dft_freq[:index_enough], stator_hysteresis_loss_dft[:index_enough], '^',alpha=0.4)
            # axes_dft[1].set_xlabel('Frequency [Hz]')
            # axes_dft[1].set_ylabel('Stator\nHysteresis\nLoss [W]')
            # axes_dft[2].plot(dft_freq[:index_enough], bxfft[:index_enough], '>',alpha=0.4)
            # axes_dft[2].plot(dft_freq[:index_enough], np.ones(len(dft_freq[:index_enough]))*threshold)
            # axes_dft[3].plot(dft_freq[:index_enough], np.ones(len(dft_freq[:index_enough]))*threshold)
            # show()

        # find the index of dft_freq that corresponds to 50e3 Hz, because higher results are contaminated by noises in fourier analysis
        index_enough = next(ind for ind, el in enumerate(dft_freq) if el > MAX_FREQUENCY) # 50e3 Hz is enough 10 times base frequency 500 Hz
        stator_eddycurrent_loss = sum( stator_eddycurrent_loss_dft[:index_enough] ) 
        stator_hysteresis_loss  = sum( stator_hysteresis_loss_dft[:index_enough] )

        fig_dft, axes_dft = subplots(4,1, sharex=True)
        axes_dft[0].plot(dft_freq[:index_enough], stator_eddycurrent_loss_dft[:index_enough], '+',alpha=0.4)
        axes_dft[0].set_xlabel('Frequency [Hz]')
        axes_dft[0].set_ylabel('Stator\nEddy\nCurrent\nLoss [W]')
        axes_dft[1].plot(dft_freq[:index_enough], stator_hysteresis_loss_dft[:index_enough], '^',alpha=0.4)
        axes_dft[1].set_xlabel('Frequency [Hz]')
        axes_dft[1].set_ylabel('Stator\nHysteresis\nLoss [W]')

        # Rotor iron loss
        rotor_eddycurrent_loss_dft = np.zeros(len(dft_freq))
        rotor_hysteresis_loss_dft = np.zeros(len(dft_freq))
        for id_element, area_element in enumerate(self.rotor_Area_data):
            rotor_Bx_waveform = self.rotor_Bx_data[id_element]
            rotor_By_waveform = self.rotor_By_data[id_element]
            bxfft = utility.singleSidedDFT(rotor_Bx_waveform, samp_freq)
            byfft = utility.singleSidedDFT(rotor_By_waveform, samp_freq)

            bxfft = remove_noises(bxfft)
            byfft = remove_noises(byfft)

            bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
            volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3

            rotor_eddycurrent_loss_dft += (ce*dft_freq**2 * volume_element/cs ) * bsq
            rotor_hysteresis_loss_dft  += (ch*dft_freq    * volume_element/cs ) * bsq
            rotor_volume += volume_element
        rotor_eddycurrent_loss = sum( rotor_eddycurrent_loss_dft[:index_enough] ) 
        rotor_hysteresis_loss  = sum( rotor_hysteresis_loss_dft[:index_enough]  )

        axes_dft[2].plot(dft_freq[:index_enough], rotor_eddycurrent_loss_dft[:index_enough], '+',alpha=0.4)
        axes_dft[2].set_xlabel('Frequency [Hz]')
        axes_dft[2].set_ylabel('\nRotor\nEddy Current\nLoss [W]')
        axes_dft[3].plot(dft_freq[:index_enough], rotor_hysteresis_loss_dft[:index_enough], '^',alpha=0.4)
        axes_dft[3].set_xlabel('Frequency [Hz]')
        axes_dft[3].set_ylabel('\nRotor\nHysteresis\nLoss [W]')


        # # Stator iron loss
        # for id_element, area_element in enumerate(self.stator_Area_data):             
        #     stator_Bx_waveform = self.stator_Bx_data[id_element]
        #     stator_By_waveform = self.stator_By_data[id_element]
        #     bxfft = utility.singleSidedDFT(stator_Bx_waveform, samp_freq)
        #     byfft = utility.singleSidedDFT(stator_By_waveform, samp_freq)
        #     bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
        #     volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3
        #     stator_eddycurrent_loss += ( ce*np.dot(dft_freq**2, bsq) * volume_element ) / cs
        #     stator_hysteresis_loss  += ( ch*np.dot(dft_freq, bsq) * volume_element ) / cs
        #     stator_volume += volume_element

        # # Rotor iron loss
        # for id_element, area_element in enumerate(self.rotor_Area_data):
        #     rotor_Bx_waveform = self.rotor_Bx_data[id_element]
        #     rotor_By_waveform = self.rotor_By_data[id_element]
        #     bxfft = utility.singleSidedDFT(rotor_Bx_waveform, samp_freq)
        #     byfft = utility.singleSidedDFT(rotor_By_waveform, samp_freq)
        #     bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
        #     volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3
        #     rotor_eddycurrent_loss += ( ce*np.dot(dft_freq**2, bsq) * volume_element ) / cs
        #     rotor_hysteresis_loss  += ( ch*np.dot(dft_freq, bsq) * volume_element ) / cs
        #     rotor_volume += volume_element

            # Copper Loss - Stator Winding Proximity Effect (Rotor side is neglected because the slip frequency is low)
            # this should be done with g==2, i.e., field data on coil area
            # % and prox losses can be totalled up in a similar way as iron loss
            # prox_loss += np.dot(cePhase * dft_freq**2, bsq) * volume_element

        # we use only 1/4 model for the loss calculation
        number_of_fraction = 4
        return ( number_of_fraction*stator_eddycurrent_loss, 
                 number_of_fraction*stator_hysteresis_loss, 
                 number_of_fraction*rotor_eddycurrent_loss, 
                 number_of_fraction*rotor_hysteresis_loss, 
                 number_of_fraction*stator_volume, 
                 number_of_fraction*rotor_volume )

# def get_magnet_loss(self):
    # pass
    # % Add up eddy current losses in the magnets
    # % Magnet properties
    # RotorMagnets = Num_pole;
    # omag = 0.556*10^6;                  % conductivity of sintered NdFeB in S/m
    # % compute fft of A at the center of each element
    # Jm=fft(A)*(2/ns);
    # for k=1:RotorMagnets
    #     g3=(g==(10+k));
    #     % total volume of the magnet under consideration;
    #     vmag=v'*g3;
    #     % average current in the magnet for each harmonic
    #     Jo=(Jm*(v.*g3))/vmag;
    #     % subtract averages off of each each element in the magnet
    #     Jm = Jm - Jo*g3';
    # end
    # %     magnet_loss = (1/2)*((omag*(2*pi*w).^2)*(abs(Jm).^2)*v);
    # magnet_loss = (1/2)*((omag*(w).^2)*(abs(Jm).^2)*v);

    # total_loss = rotor_loss + stator_loss + prox_loss + SlotOhmic + magnet_loss + Air_friction_loss;
    # results = [results; thisSpeed, rotor_loss, stator_loss, magnet_loss, SlotOhmic + prox_loss, Air_friction_loss, total_loss];

