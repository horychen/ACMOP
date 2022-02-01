import inner_rotor_motor, pyrhonen_procedure_as_function
import logging
from collections import OrderedDict
from utility import acmop_parameter
from pylab import np
from pprint import pprint

import CrossSectInnerNotchedRotor
import CrossSectStator
import Location2D

class bearingless_induction_template(inner_rotor_motor.template_machine_as_numbers):
    ''' This is bearingless induction motor template
    '''
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(bearingless_induction_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'BLIM'
        self.name = '__BLIM'

        # 初始化搜索空间
        GP = self.d['GP']
        OP = self.d['OP']
        SD = self.SD
        # childGP = OrderedDict({
        #     # IM Peculiar
        #     "deg_alpha_rm"      : acmop_parameter("free",      "magnet_pole_span_angle",        None, [None, None], lambda GP,SD:None),
        #     "mm_d_rp"           : acmop_parameter("free",      "inter_polar_iron_thickness",    None, [None, None], lambda GP,SD:None),
        #     "deg_alpha_rs"      : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "magnet_segment_span_angle",     None, [None, None], lambda GP,SD:None),
        #     "mm_d_rs"           : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "inter_segment_iron_thickness",  None, [None, None], lambda GP,SD:None),
        # })
        # GP.update(childGP)

       # Get Analytical Design
        self.Pyrhonen2009Book(fea_config_dict, SD, GP, OP) # I Stopped here.

        # 定义搜索空间，determine bounds
        original_template_neighbor_bounds = self.get_template_neighbor_bounds(GP, SD)        
        self.bounds_denorm = self.define_search_space(GP, original_template_neighbor_bounds)

        # Template's Other Properties (Shared by the swarm)
        OP = self.get_other_properties_after_geometric_parameters_are_initialized(GP, SD)
        # BEARING Winding Excitation Properties
        if True:
            OP['BeariW_zQ']         = OP['DriveW_zQ']
            OP['BeariW_CurrentAmp'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] * (OP['DriveW_CurrentAmp'] / fea_config_dict['TORQUE_CURRENT_RATIO'])
            OP['BeariW_Freq']       = OP['DriveW_Freq']
            OP['BeariW_Rs']         = OP['DriveW_Rs'] * OP['BeariW_zQ'] / OP['DriveW_zQ']
            OP['BeariW_poles']      = SD['ps']*2
            OP['slot_current_utilizing_ratio'] = fea_config_dict['SUSPENSION_CURRENT_RATIO'] + fea_config_dict['TORQUE_CURRENT_RATIO'] # will be less than 1 for separate winding

    def Pyrhonen2009Book(self, fea_config_dict, SD, GP, OP):

        air_gap_flux_density_B = SD['guess_air_gap_flux_density_B']
        # air_gap_flux_density_B = 0.9
        stator_tooth_flux_density_B_ds = SD['guess_stator_tooth_flux_density_B_ds']
        # stator_tooth_flux_density_B_ds = 1.5

        stator_yoke_flux_density_Bys = SD['guess_stator_yoke_flux_density_Bys']

        if SD['p'] >= 2:
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.75
            alpha_rm_over_alpha_rp = 1.0
            # stator_yoke_flux_density_Bys = 1.2
        else:
            # penalty for p=1 motor, i.e., large yoke height
            ROTOR_STATOR_YOKE_HEIGHT_RATIO = 0.5
            alpha_rm_over_alpha_rp = 0.75
            # stator_yoke_flux_density_Bys = 1.5

        stator_outer_diameter_Dse = 0.225 # m

        speed_rpm = SD['ExcitationFreqSimulated'] * 60 / SD['p'] # rpm

        rotor_outer_radius_r_or = pyrhonen_procedure_as_function.eric_specify_tip_speed_get_radius(SD['tip_speed'], speed_rpm)
        rotor_outer_diameter_Dr = rotor_outer_radius_r_or*2
        sleeve_length = 3
        stator_inner_radius_r_is  = rotor_outer_radius_r_or + (sleeve_length+SD['minimum_mechanical_air_gap_length_mm'])*1e-3 # m (sleeve 3 mm, air gap 0.75 mm)
        stator_inner_diameter_Dis = stator_inner_radius_r_is*2
        split_ratio = stator_inner_diameter_Dis / stator_outer_diameter_Dse

        stator_yoke_height_h_ys = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis * alpha_rm_over_alpha_rp / (2*stator_yoke_flux_density_Bys * 2*SD['p'])
        stator_tooth_height_h_ds = (stator_outer_diameter_Dse - stator_inner_diameter_Dis) / 2 - stator_yoke_height_h_ys
        stator_slot_height_h_ss = stator_tooth_height_h_ds
        stator_tooth_width_b_ds = air_gap_flux_density_B * np.pi * stator_inner_diameter_Dis / (stator_tooth_flux_density_B_ds* SD['Qs'])

        OP['stator_slot_area'] = stator_slot_area = np.pi/(4*SD['Qs']) * ((stator_outer_diameter_Dse - 2*stator_yoke_height_h_ys)**2 - stator_inner_diameter_Dis**2) - stator_tooth_width_b_ds * stator_tooth_height_h_ds

        slot_pitch_pps = np.pi * (stator_inner_diameter_Dis + stator_slot_height_h_ss) / SD['Qs']
        kov = 1.8 # \in [1.6, 2.0]
        OP['end_winding_length_Lew'] = end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + stator_tooth_width_b_ds) + slot_pitch_pps*kov * (SD['coil_pitch_y'] - 1)

        Q = SD['Qs']
        p = SD['p']
        # STATOR
        GP['deg_alpha_st'].value         = 360/Q - 2 # deg
        GP['deg_alpha_so'].value         = GP['deg_alpha_st'].value/2 # im_template uses alpha_so as 0.
        GP['mm_r_si'].value              = 1e3*stator_inner_radius_r_is # mm
        GP['mm_r_os'].value              = 1e3*stator_outer_diameter_Dse/2 # mm
        GP['mm_d_so'].value              = 1 # mm
        GP['mm_d_sp'].value              = 1.5*GP['mm_d_so'].value
        GP['mm_d_st'].value              = 1e3*(0.5*stator_outer_diameter_Dse - stator_yoke_height_h_ys) - GP['mm_r_si'].value - GP['mm_d_sp'].value  # mm
        GP['mm_d_sy'].value              = 1e3*stator_yoke_height_h_ys # mm
        GP['mm_w_st'].value              = 1e3*stator_tooth_width_b_ds # mm
        # ROTOR
        GP['mm_d_sleeve'].value          = sleeve_length
        GP['mm_d_fixed_air_gap'].value   = SD['minimum_mechanical_air_gap_length_mm']
        GP['split_ratio'].value          = split_ratio
        GP['mm_d_pm'].value              = 4  # mm
        GP['mm_d_ri'].value              = 1e3*ROTOR_STATOR_YOKE_HEIGHT_RATIO*stator_yoke_height_h_ys # TODO：This ratio (0.75) is randomly specified
        GP['mm_r_or'].value              = 1e3*rotor_outer_radius_r_or
        GP['mm_r_ri'].value              = 1e3*stator_inner_radius_r_is - GP['mm_d_pm'].value - GP['mm_d_ri'].value - GP['mm_d_sleeve'].value - GP['mm_d_fixed_air_gap'].value
        # SPMSM specific
        GP['deg_alpha_rm'].value         = 0.95*360/(2*p) # deg
        GP['mm_d_rp'].value              = 3  # mm
        GP['deg_alpha_rs'].value         = 0.975*GP['deg_alpha_rm'].value / SD['no_segmented_magnets']
        GP['mm_d_rs'].value              = 0.20*GP['mm_d_rp'].value # d_pm > d_rp and d_pm > d_rs

        # Those are some obsolete variables that are convenient to have.
        # template.Radius_OuterStatorYoke = spec_geometry_dict['Radius_OuterStatorYoke'] = 1e3*0.5*stator_outer_diameter_Dse # mm
        # template.Radius_OuterRotor      = spec_geometry_dict['Radius_OuterRotor'] = 1e3*rotor_outer_radius_r_or # mm

        # Those variables are for PMSM convenience
        # template.rotor_steel_outer_radius = spec_geometry_dict['rotor_steel_outer_radius'] = template.Radius_OuterRotor

        # required_torque = SD['mec_power']/(2*np.pi*speed_rpm)*60
        # rotor_volume_Vr = required_torque/(2*SD['TangentialStress'])
        # template.required_torque = required_torque

    def get_template_neighbor_bounds(self, GP, SD):
        ''' The bounds are determined around the template design.
        '''

        Q = self.SD['Qs']
        p = self.SD['p']
        s = self.SD['no_segmented_magnets']

        GP = self.d['GP']

        original_template_neighbor_bounds = {
            # STATOR
            "deg_alpha_st": [ 0.35*360/Q, 0.9*360/Q],
            "mm_d_so":      [  0.5,   5],                                                       
            "mm_d_st":      [0.8*GP['mm_d_st'].value, 1.2*GP['mm_d_st'].value],                
            "mm_r_os":      [1.0*GP['mm_r_os'].value, 1.2*GP['mm_r_os'].value], 
            "mm_w_st":      [0.8*GP['mm_w_st'].value, 1.2*GP['mm_w_st'].value],                
            # ROTOR
            "mm_d_sleeve":  [3,   6],
            "split_ratio":  [0.4, 0.6], # Binder-2020-MLMS-0953@Fig.7
            "mm_d_pm":      [2.5, 7],
            "mm_d_ri":      [0.8*GP['mm_d_ri'].value,  1.2*GP['mm_d_ri'].value],                              
            # SPMSM specific
            "deg_alpha_rm": [0.6*360/(2*p),          1.0*360/(2*p)],                                     
            "mm_d_rp":      [2.5,   6],                                                         
            "deg_alpha_rs": [0.8*360/(2*p)/s,        0.975*360/(2*p)/s],                               
            "mm_d_rs":      [2.5,   6]
        }
        return original_template_neighbor_bounds

    def build_design_parameters_list(self):
        GP = self.d['GP']
        SD = self.SD
        # obsolete feature
        design_parameters = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0, #'mm_r_st',
            0.0, #'mm_r_sf',
            0.0, #'mm_r_sb',
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
            ]
        return design_parameters

class bearingless_spmsm_design_variant(inner_rotor_motor.variant_machine_as_objects):

    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        if x_denorm is None:
            raise

        # 初始化父类
        super(bearingless_spmsm_design_variant, self).__init__(template, x_denorm, counter, counter_loop)

        # 检查几何变量之间是否有冲突
        GP = self.template.d['GP']
        SD = self.template.SD
        self.check_invalid_design(GP, SD)

        # Parts
        self.rotorCore = CrossSectInnerNotchedRotor.CrossSectInnerNotchedRotor(
                            name = 'NotchedRotor',
                            mm_d_pm      = GP['mm_d_pm'].value,
                            deg_alpha_rm = GP['deg_alpha_rm'].value, # angular span of the pole: class type DimAngular
                            deg_alpha_rs = GP['deg_alpha_rs'].value, # segment span: class type DimAngular
                            mm_d_ri      = GP['mm_d_ri'].value, # rotor iron thickness: class type DimLinear
                            mm_r_ri      = GP['mm_r_ri'].value, # inner radius of rotor: class type DimLinear
                            mm_d_rp      = GP['mm_d_rp'].value, # interpolar iron thickness: class type DimLinear
                            mm_d_rs      = GP['mm_d_rs'].value, # inter segment iron thickness: class type DimLinear
                            p = template.SD['p'], # Set pole-pairs to 2
                            s = template.SD['no_segmented_magnets'], # Set magnet segments/pole to 4
                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0))

        self.shaft = CrossSectInnerNotchedRotor.CrossSectShaft(name = 'Shaft',
                                                      notched_rotor = self.rotorCore
                                                    )

        self.rotorMagnet = CrossSectInnerNotchedRotor.CrossSectInnerNotchedMagnet( name = 'RotorMagnet',
                                                      notched_rotor = self.rotorCore
                                                    )

        self.stator_core = CrossSectStator.CrossSectInnerRotorStator( name = 'StatorCore',
                                            deg_alpha_st = GP['deg_alpha_st'].value, #40,
                                            deg_alpha_so = GP['deg_alpha_so'].value, #20,
                                            mm_r_si = GP['mm_r_si'].value,
                                            mm_d_so = GP['mm_d_so'].value,
                                            mm_d_sp = GP['mm_d_sp'].value,
                                            mm_d_st = GP['mm_d_st'].value,
                                            mm_d_sy = GP['mm_d_sy'].value,
                                            mm_w_st = GP['mm_w_st'].value,
                                            mm_r_st = 0.0, # =0
                                            mm_r_sf = 0.0, # =0
                                            mm_r_sb = 0.0, # =0
                                            Q = template.SD['Qs'],
                                            location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                            )

        self.coils = CrossSectStator.CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                    stator_core = self.stator_core)

        self.sleeve = CrossSectInnerNotchedRotor.CrossSectSleeve(
                            name = 'Sleeve',
                            notched_magnet = self.rotorMagnet,
                            d_sleeve = GP['mm_d_sleeve'].value
                            )

        #03 Mechanical Parameters
        self.update_mechanical_parameters()

    def check_invalid_design(self, GP, SD):
        pass

class test:
    def fea_bearingless_induction(self, im_template, x_denorm, counter, counter_loop):
        logger = logging.getLogger(__name__)
        print('Run FEA for individual #%d'%(counter))

        # get local design variant
        im_variant = population.bearingless_induction_motor_design.local_design_variant(im_template, 0, counter, x_denorm)

        # print('::', im_template.Radius_OuterRotor, im_template.Width_RotorSlotOpen)
        # print('::', im_variant.Radius_OuterRotor, im_variant.Width_RotorSlotOpen)
        # quit()

        # TODO: Change indivudal name to be more useful
        if counter_loop == 1:
            im_variant.name = 'ind%d'%(counter)
        else:
            im_variant.name = 'ind%d-redo%d'%(counter, counter_loop)
        # im_variant.spec = im_template.spec
        self.im_variant = im_variant
        self.femm_solver = FEMM_Solver.FEMM_Solver(self.im_variant, flag_read_from_jmag=False, freq=50) # eddy+static
        im = None



        if counter_loop == 1:
            self.project_name          = 'proj%d'%(counter)
        else:
            self.project_name          = 'proj%d-redo%d'%(counter, counter_loop)
        self.expected_project_file = self.output_dir + "%s.jproj"%(self.project_name)

        original_study_name = im_variant.name + "Freq"
        tran2tss_study_name = im_variant.name + 'Tran2TSS'

        self.dir_femm_temp         = self.output_dir + 'femm_temp/'
        self.femm_output_file_path = self.dir_femm_temp + original_study_name + '.csv'

        # self.jmag_control_state = False

        # local scripts
        def open_jmag(expected_project_file_path):
            if self.app is None:
                # app = win32com.client.Dispatch('designer.Application.181')
                app = win32com.client.Dispatch('designer.Application.171')
                # app = win32com.client.gencache.EnsureDispatch('designer.Application.171') # https://stackoverflow.com/questions/50127959/win32-dispatch-vs-win32-gencache-in-python-what-are-the-pros-and-cons

                if self.fea_config_dict['designer.Show'] == True:
                    app.Show()
                else:
                    app.Hide()
                # app.Quit()
                self.app = app # means that the JMAG Designer is turned ON now.
                self.bool_run_in_JMAG_Script_Editor = False

                def add_steel(self):
                    print('[First run on this computer detected]', im_template.spec_input_dict['Steel'], 'is added to jmag material library.')
                    import population
                    if 'M15' in im_template.spec_input_dict['Steel']:
                        population.add_M1xSteel(self.app, self.fea_config_dict['dir.parent'], steel_name="M-15 Steel")
                    elif 'M19' in im_template.spec_input_dict['Steel']:
                        population.add_M1xSteel(self.app, self.fea_config_dict['dir.parent'])
                    elif 'Arnon5' == im_template.spec_input_dict['Steel']:
                        population.add_Arnon5(self.app, self.fea_config_dict['dir.parent'])

                # too avoid tons of the same material in JAMG's material library
                fname = self.fea_config_dict['dir.parent'] + '.jmag_state.txt'
                if not os.path.exists(fname):
                    with open(fname, 'w') as f:
                        f.write(self.fea_config_dict['pc_name'] + '/' + im_template.spec_input_dict['Steel'] + '\n')
                    add_steel(self)
                else:
                    with open(fname, 'r') as f:
                        for line in f.readlines():
                            if self.fea_config_dict['pc_name'] + '/' + im_template.spec_input_dict['Steel'] not in line:
                                add_steel(self)

            else:
                app = self.app

            print('[acm_designer.py]', expected_project_file_path)
            if os.path.exists(expected_project_file_path):
                print('[acm_designer.py] JMAG project exists already. I learned my lessions. I will NOT delete it but create a new one with a different name instead.')
                # os.remove(expected_project_file_path)
                attempts = 2
                temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)
                while os.path.exists(temp_path):
                    attempts += 1
                    temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)

                expected_project_file_path = temp_path

            app.NewProject("Untitled")
            app.SaveAs(expected_project_file_path)
            logger.debug('Create JMAG project file: %s'%(expected_project_file_path))

            return app

        def draw_jmag(app):
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Draw the model in JMAG Designer
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            DRAW_SUCCESS = self.draw_jmag_induction(app,
                                                counter, 
                                                im_variant,
                                                im_variant.name)
            if DRAW_SUCCESS == 0:
                # TODO: skip this model and its evaluation
                cost_function = 99999 # penalty
                logging.getLogger(__name__).warn('Draw Failed for %s-%s\nCost function penalty = %g.%s', self.project_name, im_variant.name, cost_function, self.im_variant.show(toString=True))
                raise Exception('Draw Failed: Are you working on the PC? Sometime you by mistake operate in the JMAG Geometry Editor, then it fails to draw.')
                return None
            elif DRAW_SUCCESS == -1:
                raise Exception(' DRAW_SUCCESS == -1:')

            # JMAG
            if app.NumModels()>=1:
                model = app.GetModel(im_variant.name)
            else:
                logger.error('there is no model yet for %s'%(im_variant.name))
                raise Exception('why is there no model yet? %s'%(im_variant.name))
            return model

        def rotating_static_FEA():

            # wait for femm to finish, and get your slip of breakdown
            new_fname = self.dir_femm_temp + original_study_name + '.csv'
            with open(new_fname, 'r') as f:
                data = f.readlines()
                freq = float(data[0][:-1])
                torque = float(data[1][:-1])
            slip_freq_breakdown_torque, breakdown_torque, breakdown_force = freq, torque, None
            # Now we have the slip, set it up!
            im_variant.update_mechanical_parameters(slip_freq_breakdown_torque) # do this for records only

            # Must run this script after slip_freq_breakdown_torque is known to get new results, but if results are present, it is okay to use this script for post-processing.
            if not self.femm_solver.has_results():
                print('run_rotating_static_FEA')
                # utility.blockPrint()
                self.femm_solver.run_rotating_static_FEA()
                self.femm_solver.parallel_solve()
                # utility.enablePrint()

            # collecting parasolve with post-process
            # wait for .ans files
            # data_femm_solver = self.femm_solver.show_results_static(bool_plot=False) # this will wait as well?
            while not self.femm_solver.has_results():
                print(clock_time())
                sleep(3)
            results_dict = {}
            for f in [f for f in os.listdir(self.femm_solver.dir_run) if 'static_results' in f]:
                data = np.loadtxt(self.femm_solver.dir_run + f, unpack=True, usecols=(0,1,2,3))
                for i in range(len(data[0])):
                    results_dict[data[0][i]] = (data[1][i], data[2][i], data[3][i]) 
            keys_without_duplicates = [key for key, item in results_dict.items()]
            keys_without_duplicates.sort()
            with open(self.femm_solver.dir_run + "no_duplicates.txt", 'w') as fw:
                for key in keys_without_duplicates:
                    fw.writelines('%g %g %g %g\n' % (key, results_dict[key][0], results_dict[key][1], results_dict[key][2]))
            data_femm_solver = np.array([ keys_without_duplicates, 
                                             [results_dict[key][0] for key in keys_without_duplicates], 
                                             [results_dict[key][1] for key in keys_without_duplicates], 
                                             [results_dict[key][2] for key in keys_without_duplicates]])
            # print(data_femm_solver)
            # from pylab import plt
            # plt.figure(); plt.plot(data_femm_solver[0])
            # plt.figure(); plt.plot(data_femm_solver[1])
            # plt.figure(); plt.plot(data_femm_solver[2])
            # plt.figure(); plt.plot(data_femm_solver[3])
            # plt.show()
            # quit()
            return data_femm_solver

        def rotating_eddy_current_FEA(im, app, model):

            # Freq Sweeping for break-down Torque Slip
            # remember to export the B data using subroutine 
            # and check export table results only
            study = im.add_study(app, model, self.dir_csv_output_folder, choose_study_type='frequency')

            # Freq Study: you can choose to not use JMAG to find the breakdown slip.
            # Option 1: you can set im.slip_freq_breakdown_torque by FEMM Solver
            # Option 2: Use JMAG to sweeping the frequency
            # Does study has results?
            if study.AnyCaseHasResult():
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.check_csv_results(study.GetName())
            else:
                # mesh
                im.add_mesh(study, model)

                # Export Image
                    # for i in range(app.NumModels()):
                    #     app.SetCurrentModel(i)
                    #     model = app.GetCurrentModel()
                    #     app.ExportImage(r'D:\Users\horyc\OneDrive - UW-Madison\pop\run#10/' + model.GetName() + '.png')
                app.View().ShowAllAirRegions()
                # app.View().ShowMeshGeometry() # 2nd btn
                app.View().ShowMesh() # 3rn btn
                app.View().Zoom(3)
                app.View().Pan(-im.Radius_OuterRotor, 0)
                app.ExportImageWithSize('./' + model.GetName() + '.png', 2000, 2000)
                app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted because only ouput table results are selected.

                # run
                study.RunAllCases()
                app.Save()

                def check_csv_results(dir_csv_output_folder, study_name, returnBoolean=False, file_suffix='_torque.csv'): # '_iron_loss_loss.csv'
                    # print self.dir_csv_output_folder + study_name + '_torque.csv'
                    if not os.path.exists(dir_csv_output_folder + study_name + file_suffix):
                        if returnBoolean == False:
                            print('Nothing is found when looking into:', dir_csv_output_folder + study_name + file_suffix)
                            return None
                        else:
                            return False
                    else:
                        if returnBoolean == True:
                            return True

                    try:
                        # check csv results 
                        l_slip_freq = []
                        l_TorCon    = []
                        l_ForCon_X  = []
                        l_ForCon_Y  = []

                        with open(dir_csv_output_folder + study_name + '_torque.csv', 'r') as f: 
                            for ind, row in enumerate(utility.csv_row_reader(f)):
                                if ind >= 5:
                                    try:
                                        float(row[0])
                                    except:
                                        continue
                                    l_slip_freq.append(float(row[0]))
                                    l_TorCon.append(float(row[1]))

                        with open(dir_csv_output_folder + study_name + '_force.csv', 'r') as f: 
                            for ind, row in enumerate(utility.csv_row_reader(f)):
                                if ind >= 5:
                                    try:
                                        float(row[0])
                                    except:
                                        continue
                                    # l_slip_freq.append(float(row[0]))
                                    l_ForCon_X.append(float(row[1]))
                                    l_ForCon_Y.append(float(row[2]))

                        breakdown_force = max(np.sqrt(np.array(l_ForCon_X)**2 + np.array(l_ForCon_Y)**2))

                        index, breakdown_torque = utility.get_index_and_max(l_TorCon)
                        slip_freq_breakdown_torque = l_slip_freq[index]
                        return slip_freq_breakdown_torque, breakdown_torque, breakdown_force
                    except NameError as e:
                        logger = logging.getLogger(__name__)
                        logger.error('No CSV File Found.', exc_info=True)
                        raise e

                # evaluation based on the csv results
                print(':::ZZZ', self.dir_csv_output_folder)
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = check_csv_results(self.dir_csv_output_folder, study.GetName())

            # this will be used for other duplicated studies
            original_study_name = study.GetName()
            im.csv_previous_solve = self.dir_csv_output_folder + original_study_name + '_circuit_current.csv'
            im.update_mechanical_parameters(slip_freq_breakdown_torque, syn_freq=im.DriveW_Freq)


            # EC Rotate: Rotate the rotor to find the ripples in force and torque # 不关掉这些云图，跑第二个study的时候，JMAG就挂了：app.View().SetVectorView(False); app.View().SetFluxLineView(False); app.View().SetContourView(False)
            ecrot_study_name = original_study_name + "-FFVRC"
            casearray = [0 for i in range(1)]
            casearray[0] = 1
            model.DuplicateStudyWithCases(original_study_name, ecrot_study_name, casearray)

            app.SetCurrentStudy(ecrot_study_name)
            study = app.GetCurrentStudy()

            divisions_per_slot_pitch = 24 # self.fea_config_dict['ec_rotate_divisions_per_slot_pitch']  # 24
            study.GetStep().SetValue("Step", divisions_per_slot_pitch) 
            study.GetStep().SetValue("StepType", 0)
            study.GetStep().SetValue("FrequencyStep", 0)
            study.GetStep().SetValue("Initialfrequency", slip_freq_breakdown_torque)

                # study.GetCondition(u"RotCon").SetValue(u"MotionGroupType", 1)
            study.GetCondition("RotCon").SetValue("Displacement", + 360.0/im.Qr/divisions_per_slot_pitch)

            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

            # model.RestoreCadLink()
            study.Run()
            app.Save()
            # model.CloseCadLink()

        def transient_FEA_as_reference(im, slip_freq_breakdown_torque):
            # Transient Reference
            tranRef_study_name = "TranRef"
            if model.NumStudies()<4:
                model.DuplicateStudyWithType(tran2tss_study_name, "Transient2D", tranRef_study_name)
                app.SetCurrentStudy(tranRef_study_name)
                study = app.GetCurrentStudy()

                # 将一个滑差周期和十个同步周期，分成 400 * end_point / (1.0/im.DriveW_Freq) 份。
                end_point = 0.5/slip_freq_breakdown_torque + 10.0/im.DriveW_Freq
                # Pavel Ponomarev 推荐每个电周期400~600个点来捕捉槽效应。
                division = self.fea_config_dict['designer.TranRef-StepPerCycle'] * end_point / (1.0/im.DriveW_Freq)  # int(end_point * 1e4)
                                                                        # end_point = division * 1e-4
                study.GetStep().SetValue("Step", division + 1) 
                study.GetStep().SetValue("StepType", 1) # regular inverval
                study.GetStep().SetValue("StepDivision", division)
                study.GetStep().SetValue("EndPoint", end_point)

                # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
                study.GetStudyProperties().SetValue("DirectSolverType", 1)

                study.RunAllCases()
                app.Save()

        # debug for tia-iemdc-ecce-2019
        # data_femm_solver = rotating_static_FEA()
        # from show_results_iemdc19 import show_results_iemdc19
        # show_results_iemdc19(   self.dir_csv_output_folder, 
        #                         im_variant, 
        #                         femm_solver_data=data_femm_solver, 
        #                         femm_rotor_current_function=self.femm_solver.get_rotor_current_function()
        #                     )
        # quit()

        if self.fea_config_dict['bool_re_evaluate']==False:

            # this should be summoned even before initializing femm, and it will decide whether the femm results are reliable
            app = open_jmag(self.expected_project_file) # will set self.jmag_control_state to True

            ################################################################
            # Begin from where left: Frequency Study
            ################################################################
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Eddy Current Solver for Breakdown Torque and Slip
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # check for existing results
            if os.path.exists(self.femm_output_file_path):
                # for file in os.listdir(self.dir_femm_temp):
                #     if original_study_name in file:
                #         print('----------', original_study_name, file)
                print('Remove legacy femm output files @ %s'%(self.femm_output_file_path))
                os.remove(self.femm_output_file_path)
                os.remove(self.femm_output_file_path[:-4]+'.fem')
                # quit()
                # quit()

            # delete me
            # model = draw_jmag(app)

            # At this point, no results exist from femm.
            print('Run greedy_search_for_breakdown_slip...')
            femm_tic = clock_time()
            # self.femm_solver.__init__(im_variant, flag_read_from_jmag=False, freq=50.0)
            if im_variant.DriveW_poles == 4 and self.fea_config_dict['femm.use_fraction'] == True:
                print('FEMM model only solves for a fraction of 2.\n'*3)
                self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                    bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=2)
            else:
                # p >= 3 is not tested so do not use fraction for now
                self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                    bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=1) # 转子导条必须形成通路

            ################################################################
            # Begin from where left: Transient Study
            ################################################################
            model = draw_jmag(app)

            # EC-Rotate
            # rotating_eddy_current_FEA(im_variant, app, model)

            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # TranFEAwi2TSS for ripples and iron loss
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # add or duplicate study for transient FEA denpending on jmag_run_list
            # FEMM+JMAG (注意，这里我们用50Hz作为滑差频率先设置起来，等拿到breakdown slip freq的时候，再更新变量slip和study properties的时间。)
            study = im_variant.add_TranFEAwi2TSS_study( 50.0, app, model, self.dir_csv_output_folder, tran2tss_study_name, logger)
            app.SetCurrentStudy(tran2tss_study_name)
            study = app.GetCurrentStudy()
            self.mesh_study(im_variant, app, model, study)

            # wait for femm to finish, and get your slip of breakdown
            slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.femm_solver.wait_greedy_search(femm_tic)

            # Now we have the slip, set it up!
            im_variant.update_mechanical_parameters(slip_freq_breakdown_torque) # do this for records only
            if im_variant.the_slip != slip_freq_breakdown_torque / im_variant.DriveW_Freq:
                raise Exception('Check update_mechanical_parameters().')
            study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im_variant.the_slip))
            if True:
                number_of_steps_2ndTSS = self.fea_config_dict['designer.number_of_steps_2ndTSS'] 
                DM = app.GetDataManager()
                DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
                refarray = [[0 for i in range(3)] for j in range(3)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = 0.5/slip_freq_breakdown_torque #0.5 for 17.1.03l # 1 for 17.1.02y
                refarray[1][1] =    number_of_steps_2ndTSS                          # 16 for 17.1.03l #32 for 17.1.02y
                refarray[1][2] =        50
                refarray[2][0] = refarray[1][0] + 0.5/im_variant.DriveW_Freq #0.5 for 17.1.03l 
                refarray[2][1] =    number_of_steps_2ndTSS  # also modify range_ss! # don't forget to modify below!
                refarray[2][2] =        50
                DM.GetDataSet("SectionStepTable").SetTable(refarray)
                number_of_total_steps = 1 + 2 * number_of_steps_2ndTSS # [Double Check] don't forget to modify here!
                study.GetStep().SetValue("Step", number_of_total_steps)
                study.GetStep().SetValue("StepType", 3)
                study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

            # static FEA solver with FEMM (need eddy current FEA results)
            # print('::', im_variant.Omega, im_variant.Omega/2/np.pi)
            # print('::', femm_solver.im.Omega, femm_solver.im.Omega/2/np.pi)
            # quit()
            # rotating_static_FEA()

            # debug JMAG circuit
            # app.Save()
            # quit()

            # Run JMAG study
            self.run_study(im_variant, app, study, clock_time())

            # export Voltage if field data exists.
            if self.fea_config_dict['delete_results_after_calculation'] == False:
                # Export Circuit Voltage
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                # app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(self.dir_csv_output_folder + im_variant.name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(self.dir_csv_output_folder + tran2tss_study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")

            # TranRef
            # transient_FEA_as_reference(im_variant, slip_freq_breakdown_torque)

        else:
            # FEMM的转子电流，哎，是个麻烦事儿。。。

            self.femm_solver.vals_results_rotor_current = []

            new_fname = self.dir_femm_temp + original_study_name + '.csv'
            try:
                with open(new_fname, 'r') as f:
                    buf = f.readlines()
                    for idx, el in enumerate(buf):
                        if idx == 0:
                            slip_freq_breakdown_torque = float(el)
                            continue
                        if idx == 1:
                            breakdown_torque = float(el)
                            continue
                        if idx == 2:
                            self.femm_solver.stator_slot_area = float(el)
                            continue
                        if idx == 3:
                            self.femm_solver.rotor_slot_area = float(el)
                            continue
                        # print(el)
                        # print(el.split(','))
                        temp = el.split(',')
                        self.femm_solver.vals_results_rotor_current.append( float(temp[0])+ 1j*float(temp[1]) )
                        # print(self.femm_solver.vals_results_rotor_current)
                self.dirty_backup_stator_slot_area = self.femm_solver.stator_slot_area
                self.dirty_backup_rotor_slot_area = self.femm_solver.rotor_slot_area
                self.dirty_backup_vals_results_rotor_current = self.femm_solver.vals_results_rotor_current
            except FileNotFoundError as error:
                print(error)
                print('Use dirty_backup to continue...') # 有些时候，不知道为什么femm的结果文件（.csv）没了，这时候曲线救国，凑活一下吧
                self.femm_solver.stator_slot_area = self.dirty_backup_stator_slot_area           
                self.femm_solver.rotor_slot_area = self.dirty_backup_rotor_slot_area            
                self.femm_solver.vals_results_rotor_current = self.dirty_backup_vals_results_rotor_current 




            # 电机的电流值取决于槽的面积。。。。
            THE_mm2_slot_area = self.femm_solver.stator_slot_area*1e6

            if 'VariableStatorSlotDepth' in self.fea_config_dict['which_filter']:
                # set DriveW_CurrentAmp using the calculated stator slot area.
                print('[A]: DriveW_CurrentAmp is updated.')

                # 槽深变化，电密不变，所以电流也会变化。
                CurrentAmp_in_the_slot = THE_mm2_slot_area * im_variant.fill_factor * im_variant.Js*1e-6 * np.sqrt(2)
                CurrentAmp_per_conductor = CurrentAmp_in_the_slot / im_variant.DriveW_zQ
                CurrentAmp_per_phase = CurrentAmp_per_conductor * im_variant.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor

                im_variant.CurrentAmp_per_phase = CurrentAmp_per_phase

                im_variant.DriveW_CurrentAmp = self.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                im_variant.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

                slot_area_utilizing_ratio = (im_variant.DriveW_CurrentAmp + im_variant.BeariW_CurrentAmp) / im_variant.CurrentAmp_per_phase
                print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)
                
                print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                print('---im_variant.DriveW_CurrentAmp =', im_variant.DriveW_CurrentAmp)
                print('---im_variant.BeariW_CurrentAmp =', im_variant.BeariW_CurrentAmp)
                print('---TORQUE_CURRENT_RATIO:', self.fea_config_dict['TORQUE_CURRENT_RATIO'])
                print('---SUSPENSION_CURRENT_RATIO:', self.fea_config_dict['SUSPENSION_CURRENT_RATIO'])





        ################################################################
        # Load data for cost function evaluation
        ################################################################
        im_variant.results_to_be_unpacked = results_to_be_unpacked = utility.build_str_results(self.axeses, im_variant, self.project_name, tran2tss_study_name, self.dir_csv_output_folder, self.fea_config_dict, self.femm_solver)
        if results_to_be_unpacked is not None:
            if self.fig_main is not None:
                try:
                    self.fig_main.savefig(self.output_dir + im_variant.name + 'results.png', dpi=150)
                except Exception as e:
                    print('Directory exists?', self.output_dir + im_variant.name + 'results.png')
                    print(e)
                    print('\n\n\nIgnore error and continue.')
                finally:
                    utility.pyplot_clear(self.axeses)
            # show()
            return im_variant
        else:
            raise Exception('results_to_be_unpacked is None.')
        # winding analysis? 之前的python代码利用起来啊
        # 希望的效果是：设定好一个设计，马上进行运行求解，把我要看的数据都以latex报告的形式呈现出来。
        # OP_PS_Qr36_M19Gauge29_DPNV_NoEndRing.jproj

    def draw_jmag_induction(self, app, individual_index, im_variant, model_name, bool_trimDrawer_or_vanGogh=True, doNotRotateCopy=False):

        if individual_index == -1: # 后处理是-1
            print('Draw model for post-processing')
            if individual_index+1 + 1 <= app.NumModels():
                logger = logging.getLogger(__name__)
                logger.debug('The model already exists for individual with index=%d. Skip it.', individual_index)
                return -1 # the model is already drawn

        elif individual_index+1 <= app.NumModels(): # 一般是从零起步
            logger = logging.getLogger(__name__)
            logger.debug('The model already exists for individual with index=%d. Skip it.', individual_index)
            return -1 # the model is already drawn

        # open JMAG Geometry Editor
        app.LaunchGeometryEditor()
        geomApp = app.CreateGeometryEditor()
        # geomApp.Show()
        geomApp.NewDocument()
        doc = geomApp.GetDocument()
        ass = doc.GetAssembly()

        # draw parts
        try:
            if bool_trimDrawer_or_vanGogh:
                d = population.TrimDrawer(im_variant) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.plot_shaft("Shaft")

                d.plot_rotorCore("Rotor Core")
                d.plot_cage("Cage")

                d.plot_statorCore("Stator Core")

                d.plot_coil("Coil")
                # d.plot_airWithinRotorSlots(u"Air Within Rotor Slots")

                if 'VariableStatorSlotDepth' in self.fea_config_dict['which_filter']:
                    # set DriveW_CurrentAmp using the calculated stator slot area.
                    print('[A]: DriveW_CurrentAmp is updated.')

                    # 槽深变化，电密不变，所以电流也会变化。
                    CurrentAmp_in_the_slot = d.mm2_slot_area * im_variant.fill_factor * im_variant.Js*1e-6 * np.sqrt(2)
                    CurrentAmp_per_conductor = CurrentAmp_in_the_slot / im_variant.DriveW_zQ
                    CurrentAmp_per_phase = CurrentAmp_per_conductor * im_variant.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                    variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor

                    im_variant.CurrentAmp_per_phase = CurrentAmp_per_phase

                    im_variant.DriveW_CurrentAmp = self.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                    im_variant.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

                    slot_area_utilizing_ratio = (im_variant.DriveW_CurrentAmp + im_variant.BeariW_CurrentAmp) / im_variant.CurrentAmp_per_phase
                    print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)

                    print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                    print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                    print('---im_variant.DriveW_CurrentAmp =', im_variant.DriveW_CurrentAmp)
                    print('---im_variant.BeariW_CurrentAmp =', im_variant.BeariW_CurrentAmp)
                    print('---TORQUE_CURRENT_RATIO:', self.fea_config_dict['TORQUE_CURRENT_RATIO'])
                    print('---SUSPENSION_CURRENT_RATIO:', self.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            else:
                d = VanGogh_JMAG(im_variant, doNotRotateCopy=doNotRotateCopy) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.draw_model()
            self.d = d
        except Exception as e:
            print('See log file to plotting error.')
            logger = logging.getLogger(__name__)
            logger.error('The drawing is terminated. Please check whether the specified bounds are proper.', exc_info=True)

            raise e

            # print 'Draw Failed'
            # if self.pc_name == 'Y730':
            #     # and send the email to hory chen
            #     raise e

            # or you can skip this model and continue the optimization!
            return False # indicating the model cannot be drawn with the script.

        # Import Model into Designer
        doc.SaveModel(True) # True=on : Project is also saved. 
        model = app.GetCurrentModel() # model = app.GetModel(u"IM_DEMO_1")
        model.SetName(model_name)
        model.SetDescription(im_variant.model_name_prefix + '\n' + im_variant.show(toString=True))

        if doNotRotateCopy:
            im_variant.pre_process_structural(app, d.listKeyPoints)
        else:
            im_variant.pre_process(app)

        model.CloseCadLink() # this is essential if you want to create a series of models
        return True

    def run_study(self, im_variant, app, study, toc):
        logger = logging.getLogger(__name__)
        if self.fea_config_dict['designer.JMAG_Scheduler'] == False:
            print('[acm_designer.py] Run jam.exe...')
            # if run_list[1] == True:
            try:
                study.RunAllCases()
            except Exception as error:
                raise error
            msg = '[acm_designer.py] Time spent on %s is %g s.'%(study.GetName() , clock_time() - toc)
            logger.debug(msg)
            print(msg)
        else:
            print('[acm_designer.py] Submit to JMAG_Scheduler...')
            job = study.CreateJob()
            job.SetValue("Title", study.GetName())
            job.SetValue("Queued", True)
            job.Submit(False) # Fallse:CurrentCase, True:AllCases
            logger.debug('Submit %s to queue (Tran2TSS).'%(im_variant.individual_name))
            # wait and check
            # study.CheckForCaseResults()
        app.Save()

            # if the jcf file already exists, it pops a msg window
            # study.WriteAllSolidJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Solid', True) # True : Outputs cases that do not have results 
            # study.WriteAllMeshJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Mesh', True)

            # # run
            # if self.fea_config_dict['JMAG_Scheduler'] == False:
            #     study.RunAllCases()
            #     app.Save()
            # else:
            #     job = study.CreateJob()
            #     job.SetValue(u"Title", study.GetName())
            #     job.SetValue(u"Queued", True)
            #     job.Submit(True)
            #     logger.debug('Submit %s to queue (Freq).'%(im_variant.individual_name))
            #     # wait and check
            #     # study.CheckForCaseResults()

    def mesh_study(self, im_variant, app, model, study):

        # this `if' judgment is effective only if JMAG-DeleteResultFiles is False 
        # if not study.AnyCaseHasResult(): 
        # mesh
        im_variant.add_mesh(study, model)

        # Export Image
        app.View().ShowAllAirRegions()
        # app.View().ShowMeshGeometry() # 2nd btn
        app.View().ShowMesh() # 3rn btn
        app.View().Zoom(3)
        app.View().Pan(-im_variant.Radius_OuterRotor, 0)
        app.ExportImageWithSize(self.output_dir + model.GetName() + '.png', 2000, 2000)
        app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.
  