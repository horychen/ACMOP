import win32com.client
import os, logging
import numpy as np
EPS=0.01 # mm
class JMAG(object): #< ToolBase & DrawerBase & MakerExtrudeBase & MakerRevolveBase
    # JMAG Encapsulation for the JMAG Designer of JSOL Corporation.    
    def __init__(self, fea_config_dict, spec_input_dict=None):
        self.jd = None       # The activexserver selfect for JMAG Designer
        self.app = None      # app = jd
        self.projName = None # The name of JMAG Designer project (a string)
        self.geomApp = None  # The Geometry Editor selfect
        self.doc = None      # The document selfect in Geometry Editor
        self.ass = None      # The assemble selfect in Geometry Editor
        self.sketch = None   # The sketch selfect in Geometry Editor
        self.model = None    # The model selfect in JMAG Designer
        self.study = None    # The study selfect in JMAG Designer
        self.view = None     # The view selfect in JMAG Designer
        self.workDir = './'
        self.sketchNameList = []
        self.bMirror = True
        self.edge4Ref = None
        self.iRotateCopy = 0    # this is an integer
        self.consts      = None # Program constants (not used)
        self.defaultUnit = 'Millimeter' # Default length unit is mm (not used)

        self.fea_config_dict = fea_config_dict
        self.spec_input_dict = spec_input_dict

        # self.output_dir = self.fea_config_dict['dir.parent'] + self.fea_config_dict['run_folder']
        # self.dir_csv_output_folder = self.output_dir + 'csv/'
        # if not os.path.isdir(self.output_dir):
        #     os.makedirs(self.output_dir)
        # if not os.path.isdir(self.dir_csv_output_folder):
        #     os.makedirs(self.dir_csv_output_folder)

        # # post-process feature
        # self.fig_main, self.axeses = plt.subplots(2, 2, sharex=True, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        # utility.pyplot_clear(self.axeses)

        # self.folder_to_be_deleted = None

        # if os.path.exists(self.output_dir+'swarm_MOO_log.txt'):
        #     os.rename(self.output_dir+'swarm_MOO_log.txt', self.output_dir+'swarm_MOO_log_backup.txt')
        # open(self.output_dir+'swarm_MOO_log.txt', 'a').close()


        #   File "C:\Users\horyc\Anaconda3\lib\site-packages\win32com\client\dynamic.py", line 527, in __getattr__
        #     raise AttributeError("%s.%s" % (self._username_, attr))
        # AttributeError: designer.Application.171.Hide
        # 减少对app.Hide的调用，初始默认是Hide
        self.hide_or_show = True

    def open(self, expected_project_file_path):
        if self.app is None:
            try:
                app = win32com.client.Dispatch('designer.Application.171')
                # app = win32com.client.gencache.EnsureDispatch('designer.Application.171') # https://stackoverflow.com/questions/50127959/win32-dispatch-vs-win32-gencache-in-python-what-are-the-pros-and-cons
            except:
                try:
                    app = win32com.client.Dispatch('designer.Application.181')
                    # app = win32com.client.gencache.EnsureDispatch('designer.Application.171')
                except:
                    raise Exception('No JMAG Designer 17 or 18 is found in this PC.')
            if self.fea_config_dict['designer.Show'] == True:
                if self.hide_or_show == True:
                    app.Show()
                    self.hide_or_show = False
            else:
                if self.hide_or_show == False:
                    app.Hide()
                    self.hide_or_show = True
            # app.Quit()
            self.app = app # means that the JMAG Designer is turned ON now.

            def add_steel(self):
                print('[JMAG.py] [First run on %s detected]'%(self.fea_config_dict['pc_name']), self.spec_input_dict['Steel'], 'is added to jmag material library.')
                import population
                if 'M15' in self.spec_input_dict['Steel']:
                    population.add_M1xSteel(self.app, self.fea_config_dict['dir.parent'], steel_name="M-15 Steel")
                elif 'M19' in self.spec_input_dict['Steel']:
                    population.add_M1xSteel(self.app, self.fea_config_dict['dir.parent'])
                elif 'Arnon5' == self.spec_input_dict['Steel']:
                    population.add_Arnon5(self.app, self.fea_config_dict['dir.parent'])        

            # too avoid tons of the same material in JAMG's material library
            fname = self.fea_config_dict['dir.parent'] + '.jmag_state.txt'
            if not os.path.exists(fname):
                with open(fname, 'w') as f:
                    f.write(self.fea_config_dict['pc_name'] + '/' + self.spec_input_dict['Steel'] + '\n')
                add_steel(self)
            else:
                with open(fname, 'r') as f:
                    flag_already_there = False
                    for line in f.readlines():
                        print('[JMAG.py]', self.fea_config_dict['pc_name'], self.spec_input_dict['Steel'])
                        if self.fea_config_dict['pc_name'] + '/' + self.spec_input_dict['Steel'] in line:
                            flag_already_there = True
                            break
                    if flag_already_there == False:
                        add_steel(self)
        else:
            app = self.app

        print('[JMAG.py] expected_project_file_path:', expected_project_file_path)
        if os.path.exists(expected_project_file_path):
            print('[JMAG.py] JMAG project exists already. I learned my lessions. I will NOT delete it but create a new one with a different name instead.')
            # os.remove(expected_project_file_path)
            attempts = 2
            temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)
            while os.path.exists(temp_path):
                attempts += 1
                temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)

            expected_project_file_path = temp_path

        app.NewProject("Untitled")
        app.SaveAs(expected_project_file_path)
        logger = logging.getLogger(__name__)
        logger.debug('Create JMAG project file: %s'%(expected_project_file_path))
        return app

    def close(self):
        self.app.Quit()

    def save(self, name, description):
        self.doc.SaveModel(False) # True: Project is also saved. 
        model = self.app.GetCurrentModel()
        model.SetName(name)
        model.SetDescription(description)
    
    ''' PM Motor
    '''
    def pre_process_PMSM(self, app, model, acm_variant):
        # pre-process : you can select part by coordinate!
        ''' Group '''
        def group(name, id_list):
            model.GetGroupList().CreateGroup(name)
            for the_id in id_list:
                model.GetGroupList().AddPartToGroup(name, the_id)
                # model.GetGroupList().AddPartToGroup(name, name) #<- this also works

        part_ID_list = model.GetPartIDs()
        # print(part_ID_list)
        # quit()

        # view = app.View()
        # view.ClearSelect()
        # sel = view.GetCurrentSelection()
        # sel.SelectPart(123)
        # sel.SetBlockUpdateView(False)
        SI = acm_variant.template.spec_input_dict
        p = SI['p']
        s = SI['no_segmented_magnets']
        Q = SI['Qs']
                                #   轴 转子 永磁体  护套 定子 绕组
        if len(part_ID_list) != int(1 + 1 + p*2*s + 1 + 1 + Q*2):
            msg = 'Number of Parts is unexpected. Should be %d but get %d.\n'%(int(1 + 1 + p*2*s + 1 + 1 + Q*2), len(part_ID_list)) + self.show(toString=True)
            logger = logging.getLogger(__name__)
            logger.error(msg)
            raise utility.ExceptionBadNumberOfParts(msg)

        self.id_backiron = id_backiron = part_ID_list[0]
        id_shaft = part_ID_list[1]
        partIDRange_Magnet = part_ID_list[2:int(2+p*s*2)]
        id_sleeve = part_ID_list[int(2+p*s*2)]
        id_statorCore = part_ID_list[int(2+p*s*2)+1]
        partIDRange_Coil = part_ID_list[int(2+p*s*2)+2 : int(2+p*s*2)+2 + int(Q*2)]

        # debug
        # print(id_backiron)
        # print(id_shaft)
        # print(partIDRange_Magnet)
        # print(id_sleeve)
        # print(id_statorCore)
        # print(partIDRange_Coil)

        model.SuppressPart(id_sleeve, 1)

        group("Magnet", partIDRange_Magnet)
        group("Coils", partIDRange_Coil)

        ''' Add Part to Set for later references '''
        def add_part_to_set(name, x, y, ID=None):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection()
            if ID is None:
                # print x,y
                sel.SelectPartByPosition(x,y,0) # z=0 for 2D
            else:
                sel.SelectPart(ID)
            model.GetSetList().GetSet(name).AddSelected(sel)

        # def edge_set(name,x,y):
        #     model.GetSetList().CreateEdgeSet(name)
        #     model.GetSetList().GetSet(name).SetMatcherType(u"Selection")
        #     model.GetSetList().GetSet(name).ClearParts()
        #     sel = model.GetSetList().GetSet(name).GetSelection()
        #     sel.SelectEdgeByPosition(x,y,0) # sel.SelectEdge(741)
        #     model.GetSetList().GetSet(name).AddSelected(sel)
        # edge_set(u"AirGapCoast", 0, self.template.d['GP']['mm_r_or'].value+0.5*self.Length_AirGap)

        # Shaft
        add_part_to_set('ShaftSet', 0.0, 0.0, ID=id_shaft) # 坐标没用，不知道为什么，而且都给了浮点数了

        # Create Set for 4 poles Winding
        Angle_StatorSlotSpan = 360/Q
        # R = self.mm_r_si + self.mm_d_sp + self.mm_d_st *0.5 # this is not generally working (JMAG selects stator core instead.)
        # THETA = 0.25*(Angle_StatorSlotSpan)/180.*np.pi
        R = np.sqrt(acm_variant.coils.PCoil[0]**2 + acm_variant.coils.PCoil[1]**2)
        THETA = np.arctan(acm_variant.coils.PCoil[1]/acm_variant.coils.PCoil[0])
        X = R*np.cos(THETA)
        Y = R*np.sin(THETA)
        countXL = 0
        wily = acm_variant.template.d['EX']['wily']
        # try:
        #     wily.layer_X_phases
        # except AttributeError:
        #     print("[inner_rotor_motor.py] Reproduce design using jsonpickle will encounter error here: 'dict' object has no attribute 'layer_X_phases', implying that the object wily has become a dict after jsonpickle.")
        #     WILY = namedtuple('WILY', acm_variant.template.d['EX']['wily'])
        #     wily_as_obj =      WILY(**acm_variant.template.d['EX']['wily']) # https://stackoverflow.com/questions/43921240/pythonic-way-to-convert-a-dictionary-into-namedtuple-or-another-hashable-dict-li
        #     wily = wily_as_obj

        for UVW, UpDown in zip(wily.layer_X_phases,wily.layer_X_signs):
            countXL += 1 
            add_part_to_set("CoilLX%s%s %d"%(UVW,UpDown,countXL), X, Y)

            # print(X, Y, THETA)
            THETA += Angle_StatorSlotSpan/180.*np.pi
            X = R*np.cos(THETA)
            Y = R*np.sin(THETA)

        # Create Set for 2 poles Winding
        # THETA = 0.75*(Angle_StatorSlotSpan)/180.*np.pi # 这里这个角度的选择，决定了悬浮绕组产生悬浮力的方向！！！！！
        THETA = np.arctan(-acm_variant.coils.PCoil[1]/acm_variant.coils.PCoil[0]) + (2*np.pi)/Q
        X = R*np.cos(THETA)
        Y = R*np.sin(THETA)
        countYL = 0
        for UVW, UpDown in zip(wily.layer_Y_phases,wily.layer_Y_signs):
            countYL += 1 
            add_part_to_set("CoilLY%s%s %d"%(UVW,UpDown,countYL), X, Y)

            THETA += Angle_StatorSlotSpan/180.*np.pi
            X = R*np.cos(THETA)
            Y = R*np.sin(THETA)

        # Create Set for Magnets
        GP = acm_variant.template.d['GP']
        R = GP['mm_r_si'].value - GP['mm_d_sleeve'].value - GP['mm_d_fixed_air_gap'].value - 0.5*GP['mm_d_pm'].value
        alpha_rs = GP['deg_alpha_rs'].value /180*np.pi
        deg_pole_span = 360 / (p*2)

        if s>1:
            deg_alpha_notch  = (GP['deg_alpha_rm'].value - s*GP['deg_alpha_rs'].value) / (s-1) # inter-segment notch占的角度
            alpha_notch = deg_alpha_notch /180*np.pi

        list_xy_magnets = []
        # list_xy_airWithinRotorSlot = []
        for ind in range(int(p*2)):
            natural_ind = ind + 1

            if s==1:
                      # v---This negative sign means we walk CCW to assign sets.
                THETA = - (180/p-GP['deg_alpha_rm'].value + 0.5*GP['deg_alpha_rm'].value + deg_pole_span*ind) /180.*np.pi
                X = R*np.cos(THETA)
                Y = R*np.sin(THETA)

                add_part_to_set("Magnet %d"%(natural_ind), X, Y)
                list_xy_magnets.append([X,Y])
            else:     # v---This negative sign means we walk CCW to assign sets.
                THETA = - ( 180/p-GP['deg_alpha_rm'].value + 0.5*GP['deg_alpha_rs'].value + deg_pole_span*ind ) /180*np.pi # initial position
                # THETA = ( 0.5*self.deg_alpha_rs + deg_pole_span*ind ) /180*np.pi # initial position
                for s in range(s):
                    X = R*np.cos(THETA)
                    Y = R*np.sin(THETA)
                    add_part_to_set("Magnet %d s%d"%(natural_ind, s), X, Y)
                    list_xy_magnets.append([X,Y])
                    THETA -= alpha_notch + alpha_rs
                        # ^---This negative sign means we walk CCW to assign sets.

        # Create Set for Motion Region
        def part_list_set(name, list_xy, list_part_id=None, prefix=None):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection() 
            for xy in list_xy:
                sel.SelectPartByPosition(xy[0],xy[1],0) # z=0 for 2D
            if list_part_id is not None:
                for ID in list_part_id:
                    sel.SelectPart(ID)
            model.GetSetList().GetSet(name).AddSelected(sel)
        part_list_set('Motion_Region', list_xy_magnets, list_part_id=[id_backiron, id_shaft])

        part_list_set('MagnetSet', list_xy_magnets)
        return True
    def add_magnetic_transient_study(self, app, model, dir_csv_output_folder, study_name, acm_variant):
        logger = logging.getLogger(__name__)
        # spmsm_variant = self

        model.CreateStudy("Transient2D", study_name)
        app.SetCurrentStudy(study_name)
        study = model.GetStudy(study_name)

        # SS-ATA
        # study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
        # study.GetStudyProperties().SetValue("SpecifySlip", 0)
        # study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)
        # study.GetStudyProperties().SetValue(u"TimePeriodicType", 2) # This is for TP-EEC but is not effective

        # misc
        GP = acm_variant.template.d['GP']
        EX = acm_variant.template.d['EX']
        wily = EX['wily']
        # try:
        #     wily.deg_winding_U_phase_phase_axis_angle
        # except AttributeError:
        #     print("[inner_rotor_motor.py] Reproduce design using jsonpickle will encounter error here: 'dict' object has no attribute 'deg_winding_U_phase_phase_axis_angle', implying that the object wily has become a dict after jsonpickle.")
        #     WILY = namedtuple('WILY', acm_variant.template.d['EX']['wily'])
        #     wily_as_obj =      WILY(**acm_variant.template.d['EX']['wily']) # https://stackoverflow.com/questions/43921240/pythonic-way-to-convert-a-dictionary-into-namedtuple-or-another-hashable-dict-li
        #     wily = wily_as_obj
        
        study.GetStudyProperties().SetValue("ConversionType", 0)
        study.GetStudyProperties().SetValue("NonlinearMaxIteration", acm_variant.template.fea_config_dict['designer.max_nonlinear_iteration'])
        study.GetStudyProperties().SetValue("ModelThickness", EX['mm_template_stack_length']) # [mm] Stack Length

        # Material
        self.add_material(study, acm_template=acm_variant.template)

        # Conditions - Motion
        study.CreateCondition("RotationMotion", "RotCon") # study.GetCondition(u"RotCon").SetXYZPoint(u"", 0, 0, 1) # megbox warning
        # print('the_speed:', acm_variant.the_speed)
        study.GetCondition("RotCon").SetValue("AngularVelocity", int(EX['the_speed']))
        study.GetCondition("RotCon").ClearParts()
        study.GetCondition("RotCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)
        # Implementation of id=0 control:
        #   After rotate the rotor by half the inter-pole notch span, The d-axis initial position is at pole pitch angle divided by 2.
        #   The U-phase current is sin(omega_syn*t) = 0 at t=0 and requires the d-axis to be at the winding phase axis (to obtain id=0 control)
        deg_pole_span = 180/acm_variant.template.SI['p']
        #                                                              inter-pole notch (0.5 for half)         rotate to x-axis    winding placing bias (half adjacent slot angle)      reverse north and south pole to make torque positive.
        print('[inner_rotor_motor.py] [PMSM JMAG] InitialRotationAngle :',(deg_pole_span-GP['deg_alpha_rm'].value)*0.5, - deg_pole_span*0.5, + wily.deg_winding_U_phase_phase_axis_angle,     + deg_pole_span)
        print('[inner_rotor_motor.py] [PMSM JMAG] InitialRotationAngle =',(deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + wily.deg_winding_U_phase_phase_axis_angle     + deg_pole_span, 'deg')
        study.GetCondition("RotCon").SetValue(u"InitialRotationAngle",    (deg_pole_span-GP['deg_alpha_rm'].value)*0.5 - deg_pole_span*0.5 + wily.deg_winding_U_phase_phase_axis_angle     + deg_pole_span) 


        study.CreateCondition("Torque", "TorCon") # study.GetCondition(u"TorCon").SetXYZPoint(u"", 0, 0, 0) # megbox warning
        study.GetCondition("TorCon").SetValue("TargetType", 1)
        study.GetCondition("TorCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("TorCon").ClearParts()

        study.CreateCondition("Force", "ForCon")
        study.GetCondition("ForCon").SetValue("TargetType", 1)
        study.GetCondition("ForCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("ForCon").ClearParts()


        # Conditions - FEM Coils & Conductors (i.e. stator/rotor winding)
        self.add_circuit(app, model, study, acm_variant, bool_3PhaseCurrentSource=wily.bool_3PhaseCurrentSource)


        # True: no mesh or field results are needed
        study.GetStudyProperties().SetValue("OnlyTableResults", acm_variant.template.fea_config_dict['designer.OnlyTableResults'])

        # Linear Solver
        if False:
            # sometime nonlinear iteration is reported to fail and recommend to increase the accerlation rate of ICCG solver
            study.GetStudyProperties().SetValue("IccgAccel", 1.2) 
            study.GetStudyProperties().SetValue("AutoAccel", 0)
        else:
            # this can be said to be super fast over ICCG solver.
            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

        if acm_variant.template.fea_config_dict['designer.MultipleCPUs'] == True:
            # This SMP(shared memory process) is effective only if there are tons of elements. e.g., over 100,000.
            # too many threads will in turn make them compete with each other and slow down the solve. 2 is good enough for eddy current solve. 6~8 is enough for transient solve.
            study.GetStudyProperties().SetValue("UseMultiCPU", True)
            study.GetStudyProperties().SetValue("MultiCPU", 2) 

        # 上一步的铁磁材料的状态作为下一步的初值，挺好，但是如果每一个转子的位置转过很大的话，反而会减慢非线性迭代。
        # 我们的情况是：0.33 sec 分成了32步，每步的时间大概在0.01秒，0.01秒乘以0.5*497 Hz = 2.485 revolution...
        # study.GetStudyProperties().SetValue(u"NonlinearSpeedup", 0) # JMAG17.1以后默认使用。现在后面密集的步长还多一点（32步），前面16步慢一点就慢一点呗！

        # two sections of different time step size
        if True:
            number_cycles_in_1stTSS = acm_variant.template.fea_config_dict['designer.number_cycles_in_1stTSS']
            number_cycles_in_2ndTSS = acm_variant.template.fea_config_dict['designer.number_cycles_in_2ndTSS']
            number_cycles_in_3rdTSS = acm_variant.template.fea_config_dict['designer.number_cycles_in_3rdTSS']
            number_cycles_prolonged = acm_variant.template.fea_config_dict['designer.number_cycles_prolonged']
            number_of_steps_1stTSS = acm_variant.template.fea_config_dict['designer.number_of_steps_1stTSS'] 
            number_of_steps_2ndTSS = acm_variant.template.fea_config_dict['designer.number_of_steps_2ndTSS'] 
            number_of_steps_3rdTSS = number_cycles_in_3rdTSS*acm_variant.template.fea_config_dict['designer.StepPerCycle_3rdTSS']
            DM = app.GetDataManager()
            DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
            if number_cycles_prolonged == 0:
                if number_cycles_in_3rdTSS == 0:
                    refarray = [[0 for i in range(3)] for j in range(3)]
                    refarray[0][0] = 0
                    refarray[0][1] =    1
                    refarray[0][2] =        50
                    refarray[1][0] = number_cycles_in_1stTSS/EX['DriveW_Freq']
                    refarray[1][1] =    number_of_steps_1stTSS
                    refarray[1][2] =        50
                    refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/EX['DriveW_Freq']
                    refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                    refarray[2][2] =        50
                else:
                    refarray = [[0 for i in range(3)] for j in range(4)]
                    refarray[0][0] = 0
                    refarray[0][1] =    1
                    refarray[0][2] =        50
                    refarray[1][0] = number_cycles_in_1stTSS/EX['DriveW_Freq']
                    refarray[1][1] =    number_of_steps_1stTSS
                    refarray[1][2] =        50
                    refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/EX['DriveW_Freq']
                    refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                    refarray[2][2] =        50
                    refarray[3][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS+number_cycles_in_3rdTSS)/EX['DriveW_Freq']
                    refarray[3][1] =    number_of_steps_3rdTSS
                    refarray[3][2] =        50
            else:
                refarray = [[0 for i in range(3)] for j in range(4)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = number_cycles_in_1stTSS/EX['DriveW_Freq']
                refarray[1][1] =    number_of_steps_1stTSS
                refarray[1][2] =        50
                refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/EX['DriveW_Freq']
                refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                refarray[2][2] =        50
                refarray[3][0] = refarray[2][0] + number_cycles_prolonged/EX['DriveW_Freq'] 
                refarray[3][1] =    number_cycles_prolonged*acm_variant.template.fea_config_dict['designer.TranRef-StepPerCycle'] 
                refarray[3][2] =        50
                # refarray[4][0] = refarray[3][0] + 0.5/EX['DriveW_Freq'] # 最后来一个超密的半周期400步
                # refarray[4][1] =    400
                # refarray[4][2] =        50
            number_of_total_steps = 1 + number_of_steps_1stTSS + number_of_steps_2ndTSS + number_of_steps_3rdTSS + number_cycles_prolonged*acm_variant.template.fea_config_dict['designer.TranRef-StepPerCycle'] # [Double Check] don't forget to modify here!
            print('[inner_rotor_motor.py]: refarray:', refarray)
            DM.GetDataSet("SectionStepTable").SetTable(refarray)
            study.GetStep().SetValue("Step", number_of_total_steps)
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

        # add equations
        study.GetDesignTable().AddEquation("freq")
        study.GetDesignTable().AddEquation("speed")
        study.GetDesignTable().GetEquation("freq").SetType(0)
        study.GetDesignTable().GetEquation("freq").SetExpression("%g"%((EX['DriveW_Freq'])))
        study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency")
        study.GetDesignTable().GetEquation("speed").SetType(1)
        study.GetDesignTable().GetEquation("speed").SetExpression("freq * %d"%(60/(EX['DriveW_poles']/2)))
        study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed of four pole")

        # speed, freq, slip
        study.GetCondition("RotCon").SetValue("AngularVelocity", 'speed')
        if acm_variant.template.spec_input_dict['DPNV_or_SEPA']==False:
            app.ShowCircuitGrid(True)
            study.GetCircuit().GetComponent("CS4").SetValue("Frequency", "freq")
            study.GetCircuit().GetComponent("CS2").SetValue("Frequency", "freq")

        # max_nonlinear_iteration = 50
        # study.GetStudyProperties().SetValue(u"NonlinearMaxIteration", max_nonlinear_iteration)
        # study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
        # study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)

        # # add other excitation frequencies other than 500 Hz as cases
        # for case_no, DriveW_Freq in enumerate([50.0, slip_freq_breakdown_torque]):
        #     slip = slip_freq_breakdown_torque / DriveW_Freq
        #     study.GetDesignTable().AddCase()
        #     study.GetDesignTable().SetValue(case_no+1, 0, DriveW_Freq)
        #     study.GetDesignTable().SetValue(case_no+1, 1, slip)

        # 你把Tran2TSS计算周期减半！
        # 也要在计算铁耗的时候选择1/4或1/2的数据！（建议1/4）
        # 然后，手动添加end step 和 start step，这样靠谱！2019-01-09：注意设置铁耗条件（iron loss condition）的Reference Start Step和End Step。

        # Iron Loss Calculation Condition
        # Stator 
        if True:
            cond = study.CreateCondition("Ironloss", "IronLossConStator")
            cond.SetValue("RevolutionSpeed", "freq*60/%d"%(0.5*EX['DriveW_poles']))
            cond.ClearParts()
            sel = cond.GetSelection()
            sel.SelectPartByPosition(acm_variant.template.d['GP']['mm_r_si'].value + EPS, 0 ,0) # 2022-02-04 这里发现代码有点歧义：注意，实际上acm_variant.template.d['GP']已经被修改了，acm_variant.template.d['GP'] = acm_variant.GP。
            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results and to have a FFT plot
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3) # 3:Custom
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            cond.SetValue("StartReferenceStep", number_of_total_steps+1 - 0.5*int(number_of_steps_2ndTSS/number_cycles_in_2ndTSS)) # 1/2 period <=> number_of_steps_2ndTSS
                # cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTSS*0.5) # 1/4 period <=> number_of_steps_2ndTSS*0.5
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
        # Check CSV reults for iron loss (You cannot check this for Freq study) # CSV and save space
        study.GetStudyProperties().SetValue("CsvOutputPath", dir_csv_output_folder) # it's folder rather than file!
        study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss")
        study.GetStudyProperties().SetValue("DeleteResultFiles", acm_variant.template.fea_config_dict['delete_results_after_calculation'])
        # Terminal Voltage/Circuit Voltage: Check for outputing CSV results 
        study.GetCircuit().CreateTerminalLabel("TerminalGroupACU", 8, -13)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupACV", 8, -11)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupACW", 8, -9)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupBDU", 23, -13)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupBDV", 23, -11)
        study.GetCircuit().CreateTerminalLabel("TerminalGroupBDW", 23, -9)
        # Export Stator Core's field results only for iron loss calculation (the csv file of iron loss will be clean with this setting)
            # study.GetMaterial(u"Rotor Core").SetValue(u"OutputResult", 0) # at least one part on the rotor should be output or else a warning "the jplot file does not contains displacement results when you try to calc. iron loss on the moving part." will pop up, even though I don't add iron loss condition on the rotor.
        # study.GetMeshControl().SetValue(u"AirRegionOutputResult", 0)
        # study.GetMaterial("Shaft").SetValue("OutputResult", 0)
        # study.GetMaterial("Cage").SetValue("OutputResult", 0)
        # study.GetMaterial("Coil").SetValue("OutputResult", 0)
        # Rotor
        if True:
            cond = study.CreateCondition("Ironloss", "IronLossConRotor")
            cond.SetValue("BasicFrequencyType", 2)
            cond.SetValue("BasicFrequency", "freq")
                # cond.SetValue(u"BasicFrequency", u"slip*freq") # this require the signal length to be at least 1/4 of slip period, that's too long!
            cond.ClearParts()
            sel = cond.GetSelection()
            # sel.SelectPartByPosition(acm_variant.mm_r_ri + EPS, 0 ,0) # Why this is not working??? Because it is integer.... you must use 0.0 instead of 0!!!
            sel.SelectPart(self.id_backiron)

            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3)
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            cond.SetValue("StartReferenceStep", number_of_total_steps+1 - 0.5*int(number_of_steps_2ndTSS/number_cycles_in_2ndTSS)) # 1/2 period <=> number_of_steps_2ndTSS
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 2) # specify reference steps for 1/2 or 1/4 period and extend it to whole period (2 means 1/2 peirodicity; 4 means 1/4 peirodicity)
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
        self.study_name = study_name
        return study
    def add_structural_static_study(self):
        pass
    def add_mesh(self, study, model):
        pass
    # TranFEAwi2TSS
    def add_material(self, study, acm_template):
        if 'M19' in acm_template.spec_input_dict['Steel']:
            study.SetMaterialByName("StatorCore", "M-19 Steel Gauge-29")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 95)
                # study.GetMaterial(u"Stator Core").SetValue(u"UserConductivityValue", 1900000)

            study.SetMaterialByName("NotchedRotor", "M-19 Steel Gauge-29")
            study.GetMaterial("NotchedRotor").SetValue("Laminated", 1)
            study.GetMaterial("NotchedRotor").SetValue("LaminationFactor", 95)

        elif 'M15' in acm_template.spec_input_dict['Steel']:
            study.SetMaterialByName("StatorCore", "M-15 Steel")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 98)

            study.SetMaterialByName("NotchedRotor", "M-15 Steel")
            study.GetMaterial("NotchedRotor").SetValue("Laminated", 1)
            study.GetMaterial("NotchedRotor").SetValue("LaminationFactor", 98)

        elif acm_template.spec_input_dict['Steel'] == 'Arnon5':
            study.SetMaterialByName("StatorCore", "Arnon5-final")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 96)

            study.SetMaterialByName("NotchedRotor", "Arnon5-final")
            study.GetMaterial("NotchedRotor").SetValue("Laminated", 1)
            study.GetMaterial("NotchedRotor").SetValue("LaminationFactor", 96)

        else:
            msg = 'Warning: default material is used: DCMagnetic Type/50A1000.'
            print(msg)
            logging.getLogger(__name__).warn(msg)
            study.SetMaterialByName("StatorCore", "DCMagnetic Type/50A1000")
            study.GetMaterial("StatorCore").SetValue("UserConductivityType", 1)
            study.SetMaterialByName("NotchedRotor", "DCMagnetic Type/50A1000")
            study.GetMaterial("NotchedRotor").SetValue("UserConductivityType", 1)

        study.SetMaterialByName("Coils", "Copper")
        study.GetMaterial("Coils").SetValue("UserConductivityType", 1)

        # study.SetMaterialByName("Cage", "Aluminium")
        # study.GetMaterial("Cage").SetValue("EddyCurrentCalculation", 1)
        # study.GetMaterial("Cage").SetValue("UserConductivityType", 1)
        # study.GetMaterial("Cage").SetValue("UserConductivityValue", self.Bar_Conductivity)

        # N40H Reversible
        study.SetMaterialByName(u"Magnet", u"Arnold/Reversible/N40H")
        study.GetMaterial(u"Magnet").SetValue(u"EddyCurrentCalculation", 1)
        available_temperature_list = [-40, 20, 60, 80, 100, 120, 150, 180, 200, 220] # according to JMAG
        magnet_temperature = min(available_temperature_list, key=lambda x:abs(x-acm_template.spec_input_dict['Temperature']))        
        print('[inner_rotor_motor.py] magnet_temperature is', magnet_temperature, 'deg C.')
        study.GetMaterial(u"Magnet").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)
        study.GetMaterial(u"Magnet").SetValue(u"Poles", acm_template.d['EX']['DriveW_poles'])

        study.GetMaterial(u"Magnet").SetDirectionXYZ(1, 0, 0)
        study.GetMaterial(u"Magnet").SetAxisXYZ(0, 0, 1)
        study.GetMaterial(u"Magnet").SetOriginXYZ(0, 0, 0)
        study.GetMaterial(u"Magnet").SetPattern(u"RadialCircular")
        study.GetMaterial(u"Magnet").SetValue(u"StartAngle", 0.5* 360/(2*acm_template.SI['p']) ) # 半个极距


        # add_carbon_fiber_material(app)
    def add_circuit(self, app, model, study, acm_variant, bool_3PhaseCurrentSource=True):
        # Circuit - Current Source
        app.ShowCircuitGrid(True)
        study.CreateCircuit()

        # 4 pole motor Qs=24 dpnv implemented by two layer winding (6 coils). In this case, drive winding has the same slot turns as bearing winding
        def circuit(Grouping,turns,Rs,ampD,ampB,freq,phase=0, CommutatingSequenceD=0, CommutatingSequenceB=0, x=10,y=10, bool_3PhaseCurrentSource=True):
            study.GetCircuit().CreateSubCircuit("Star Connection", "Star Connection %s"%(Grouping), x, y) # è¿™äº›æ•°å­—æŒ‡çš„æ˜¯gridçš„ä¸ªæ•°ï¼Œç¬¬å‡ è¡Œç¬¬å‡ åˆ—çš„æ ¼ç‚¹å¤„
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil1").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil1").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil2").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil2").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil3").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil3").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil1").SetName("CircuitCoil%sU"%(Grouping))
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil2").SetName("CircuitCoil%sV"%(Grouping))
            study.GetCircuit().GetSubCircuit("Star Connection %s"%(Grouping)).GetComponent("Coil3").SetName("CircuitCoil%sW"%(Grouping))
            # Star Connection_2 is GroupAC
            # Star Connection_4 is GroupBD

            if bool_3PhaseCurrentSource == True: # must use this for frequency analysis

                study.GetCircuit().CreateComponent("3PhaseCurrentSource", "CS%s"%(Grouping))
                study.GetCircuit().CreateInstance("CS%s"%(Grouping), x-4, y+1)
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("Amplitude", ampD+ampB)
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("Frequency", "freq") # this is not needed for freq analysis # "freq" is a variable
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("PhaseU", phase)
                # Commutating sequence is essencial for the direction of the field to be consistent with speed: UVW rather than UWV
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("CommutatingSequence", CommutatingSequenceD) 
            else: 
                I1 = "CS%s-1"%(Grouping)
                I2 = "CS%s-2"%(Grouping)
                I3 = "CS%s-3"%(Grouping)
                study.GetCircuit().CreateComponent("CurrentSource", I1)
                study.GetCircuit().CreateInstance(                   I1, x-4, y+3)
                study.GetCircuit().CreateComponent("CurrentSource", I2)
                study.GetCircuit().CreateInstance(                   I2, x-4, y+1)
                study.GetCircuit().CreateComponent("CurrentSource", I3)
                study.GetCircuit().CreateInstance(                   I3, x-4, y-1)

                phase_shift_drive = -120 if CommutatingSequenceD == 1 else 120
                phase_shift_beari = -120 if CommutatingSequenceB == 1 else 120

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 0*phase_shift_drive) # "freq" variable cannot be used here. So pay extra attension here when you create new case of a different freq.
                f2 = app.FunctionFactory().Sin(ampB, freq, 0*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I1).SetFunction(func)

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 1*phase_shift_drive)
                f2 = app.FunctionFactory().Sin(ampB, freq, 1*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I2).SetFunction(func)

                func = app.FunctionFactory().Composite()
                f1 = app.FunctionFactory().Sin(ampD, freq, 2*phase_shift_drive)
                f2 = app.FunctionFactory().Sin(ampB, freq, 2*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I3).SetFunction(func)

            study.GetCircuit().CreateComponent("Ground", "Ground")
            study.GetCircuit().CreateInstance("Ground", x+2, y+1)
        # 这里电流幅值中的0.5因子源自DPNV导致的等于2的平行支路数。没有考虑到这一点，是否会对initial design的有效性产生影响？
        # 仔细看DPNV的接线，对于转矩逆变器，绕组的并联支路数为2，而对于悬浮逆变器，绕组的并联支路数为1。

        EX = acm_variant.template.d['EX']
        wily = EX['wily']

        npb = wily.number_parallel_branch
        nwl = wily.number_winding_layer # number of windign layers 
        # if acm_variant.template.fea_config_dict['DPNV_separate_winding_implementation'] == True or acm_variant.template.spec_input_dict['DPNV_or_SEPA'] == False:
        if acm_variant.template.spec_input_dict['DPNV_or_SEPA'] == False:
            # either a separate winding or a DPNV winding implemented as a separate winding
            ampD =  0.5 * (EX['DriveW_CurrentAmp']/npb + acm_variant.BeariW_CurrentAmp) # 为了代码能被四极电机和二极电机通用，代入看看就知道啦。
            ampB = -0.5 * (EX['DriveW_CurrentAmp']/npb - acm_variant.BeariW_CurrentAmp) # 关于符号，注意下面的DriveW对应的circuit调用时的ampB前还有个负号！
            if bool_3PhaseCurrentSource != True:
                raise Exception('Logic Error Detected.')
        else:
            # case: DPNV as an actual two layer winding
            ampD = EX['DriveW_CurrentAmp']/npb
            ampB = EX['BeariW_CurrentAmp']
            if bool_3PhaseCurrentSource != False:
                raise Exception('Logic Error Detected.')

        circuit('GroupAC',  EX['DriveW_zQ']/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=EX['DriveW_Rs'],ampD= ampD,
                              ampB=-ampB, freq=acm_variant.template.d['EX']['DriveW_Freq'], phase=0,
                              CommutatingSequenceD=wily.CommutatingSequenceD,
                              CommutatingSequenceB=wily.CommutatingSequenceB)
        circuit('GroupBD',  EX['BeariW_zQ']/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=EX['BeariW_Rs'],ampD= ampD,
                              ampB=+ampB, freq=acm_variant.template.d['EX']['BeariW_Freq'], phase=0,
                              CommutatingSequenceD=wily.CommutatingSequenceD,
                              CommutatingSequenceB=wily.CommutatingSequenceB,x=25) # CS4 corresponds to uauc (conflict with following codes but it does not matter.)

        # Link FEM Coils to Coil Set     
        # if acm_variant.template.fea_config_dict['DPNV_separate_winding_implementation'] == True or acm_variant.template.spec_input_dict['DPNV_or_SEPA'] == False:
        if acm_variant.template.spec_input_dict['DPNV_or_SEPA'] == False:
            def link_FEMCoils_2_CoilSet(Grouping,l1,l2):
                # link between FEM Coil Condition and Circuit FEM Coil
                for UVW in ['U','V','W']:
                    which_phase = "%s%s-Phase"%(Grouping,UVW)
                    study.CreateCondition("FEMCoil", which_phase)
                    condition = study.GetCondition(which_phase)
                    condition.SetLink("CircuitCoil%s%s"%(Grouping,UVW))
                    condition.GetSubCondition("untitled").SetName("Coil Set 1")
                    condition.GetSubCondition("Coil Set 1").SetName("delete")
                count = 0
                dict_dir = {'+':1, '-':0, 'o':None}
                # select the part to assign the FEM Coil condition
                for UVW, UpDown in zip(l1,l2):
                    count += 1 
                    if dict_dir[UpDown] is None:
                        # print 'Skip', UVW, UpDown
                        continue
                    which_phase = "%s%s-Phase"%(Grouping,UVW)
                    condition = study.GetCondition(which_phase)
                    condition.CreateSubCondition("FEMCoilData", "Coil Set %d"%(count))
                    subcondition = condition.GetSubCondition("Coil Set %d"%(count))
                    subcondition.ClearParts()
                    subcondition.AddSet(model.GetSetList().GetSet("Coil%s%s%s %d"%(Grouping,UVW,UpDown,count)), 0) # 未修改，这里的和Set对不上的，Grouping要换成LX或LY
                    subcondition.SetValue("Direction2D", dict_dir[UpDown])
                # clean up
                for UVW in ['U','V','W']:
                    which_phase = "%s%s-Phase"%(Grouping,UVW)
                    condition = study.GetCondition(which_phase)
                    condition.RemoveSubCondition("delete")
            link_FEMCoils_2_CoilSet('GroupAC', # 这边还是古老的，还未修改啊。。。
                                    wily.dict_coil_connection['layer X phases'], # 40 for 4 poles, +1 for UVW, 
                                    wily.dict_coil_connection['layer X signs'])                 # +2 for up or down, 
            link_FEMCoils_2_CoilSet('GroupBD', 
                                    wily.dict_coil_connection['layer Y phases'], # 20 for 2 poles, +1 for UVW, .
                                    wily.dict_coil_connection['layer Y signs'])                 # +2 for up or down,  这里的2和4等价于leftlayer和rightlayer。
        else:
            # 两个改变，一个是激励大小的改变（本来是200A 和 5A，现在是205A和195A），
            # 另一个绕组分组的改变，现在的A相是上层加下层为一相，以前是用俩单层绕组等效的。

            # Link FEM Coils to Coil Set as double layer short pitched winding
            # Create FEM Coil Condition
            # here we map circuit component `Coil2A' to FEM Coil Condition 'phaseAuauc
            # here we map circuit component `Coil4A' to FEM Coil Condition 'phaseAubud
            for suffix, poles in zip(['GroupAC', 'GroupBD'], [EX['DriveW_poles'], EX['BeariW_poles']]): # 仍然需要考虑poles，是因为为Coil设置Set那里的代码还没有更新。这里的2(acm_variant.DriveW_poles)和4(acm_variant.BeariW_poles)等价于leftlayer和rightlayer。
                for UVW in ['U','V','W']:
                    study.CreateCondition("FEMCoil", 'phase'+UVW+suffix)
                    # link between FEM Coil Condition and Circuit FEM Coil
                    condition = study.GetCondition('phase'+UVW+suffix)
                    condition.SetLink("CircuitCoil%s%s"%(suffix,UVW))
                    condition.GetSubCondition("untitled").SetName("delete")
            countXL = 0 # countXL indicates which slot the current rightlayer is in.
            index = 0
            dict_dir = {'+':1, '-':0}
            coil_pitch = wily.coil_pitch_y #wily.dict_coil_connection[0]
            # select the part (via `Set') to assign the FEM Coil condition
            for UVW, UpDown in zip(wily.layer_X_phases, wily.layer_X_signs):

                countXL += 1 
                if wily.grouping_AC[index] == 1:
                    suffix = 'GroupAC'
                else:
                    suffix = 'GroupBD'
                condition = study.GetCondition('phase'+UVW+suffix)

                # right layer
                # print (countXL, "Coil Set %d"%(countXL), end=' ')
                condition.CreateSubCondition("FEMCoilData", "Coil Set Layer X %d"%(countXL))
                subcondition = condition.GetSubCondition("Coil Set Layer X %d"%(countXL))
                subcondition.ClearParts()
                subcondition.AddSet(model.GetSetList().GetSet("CoilLX%s%s %d"%(UVW,UpDown,countXL)), 0) # poles=4 means right layer, rather than actual poles
                subcondition.SetValue("Direction2D", dict_dir[UpDown])

                # left layer
                Q = acm_variant.template.SI['Qs']
                if coil_pitch > 0:
                    if countXL+coil_pitch <= Q:
                        count_leftlayer = countXL+coil_pitch
                        index_leftlayer = index+coil_pitch
                    else:
                        count_leftlayer = int(countXL+coil_pitch - Q)
                        index_leftlayer = int(index+coil_pitch - Q)
                else:
                    if countXL+coil_pitch > 0:
                        count_leftlayer = countXL+coil_pitch
                        index_leftlayer = index+coil_pitch
                    else:
                        count_leftlayer = int(countXL+coil_pitch + Q)
                        index_leftlayer = int(index+coil_pitch + Q)

                # Check if it is a distributed windg???
                if wily.distributed_or_concentrated == False:
                    print('[inner_rotor_motor.py] Concentrated winding!')
                    UVW    = wily.layer_Y_phases[index_leftlayer]
                    UpDown = wily.layer_Y_signs[index_leftlayer]
                else:
                    # print('Distributed winding.')
                    if wily.layer_Y_phases[index_leftlayer] != UVW:
                        print('[inner_rotor_motor.py] [Warning] Potential bug in your winding layout detected.')
                        raise Exception('Bug in winding layout detected.')
                    # 右层导体的电流方向是正，那么与其串联的一个coil_pitch之处的左层导体就是负！不需要再检查l_leftlayer2了~
                    if UpDown == '+': 
                        UpDown = '-'
                    else:
                        UpDown = '+'
                # print (count_leftlayer, "Coil Set %d"%(count_leftlayer))
                condition.CreateSubCondition("FEMCoilData", "Coil Set Layer Y %d"%(count_leftlayer))
                subcondition = condition.GetSubCondition("Coil Set Layer Y %d"%(count_leftlayer))
                subcondition.ClearParts()
                subcondition.AddSet(model.GetSetList().GetSet("CoilLY%s%s %d"%(UVW,UpDown,count_leftlayer)), 0) # poles=2 means left layer, rather than actual poles
                subcondition.SetValue("Direction2D", dict_dir[UpDown])
                # print 'coil_pitch=', coil_pitch
                # print layer_X_phases[index], UVW
                # print l_leftlayer1[index_leftlayer]
                # print layer_X_phases
                # print l_leftlayer1
                index += 1
            # clean up
            for suffix in ['GroupAC', 'GroupBD']:
                for UVW in ['U','V','W']:
                    condition = study.GetCondition('phase'+UVW+suffix)
                    condition.RemoveSubCondition("delete")
            # raise Exception('Test DPNV PE.')


    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 画图
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def addConstraintCocentricity(self, vA, vB):
        print(vA.GetName(), vB.GetName())
        ref1 = self.sketch.GetItem(vA.GetName())
        ref2 = self.doc.CreateReferenceFromItem(ref1)
        ref3 = self.sketch.GetItem(vB.GetName())
        ref4 = self.doc.CreateReferenceFromItem(ref3)
        self.sketch.CreateBiConstraint(u"concentricity", ref2, ref4)

        # ref1 = geomApp.GetDocument().GetAssembly().GetItem(u"StatorCore").GetItem(u"Vertex.7")
        # ref2 = geomApp.GetDocument().CreateReferenceFromItem(ref1)
        # ref3 = geomApp.GetDocument().GetAssembly().GetItem(u"StatorCore").GetItem(u"Vertex.9")
        # ref4 = geomApp.GetDocument().CreateReferenceFromItem(ref3)
        # geomApp.GetDocument().GetAssembly().GetItem(u"StatorCore").CreateBiConstraint(u"concentricity", ref2, ref4)
        # geomApp.GetDocument().GetAssembly().GetItem(u"StatorCore").GetItem(u"Vertex.9").SetProperty(u"X", 56.2082073145118)
        # geomApp.GetDocument().GetAssembly().GetItem(u"StatorCore").GetItem(u"Vertex.9").SetProperty(u"Y", -32.4518236236997)

    def drawLine(self, startxy, endxy, returnVertexName=False):
        # DRAWLINE Draw a line.
        #    drawLine([start_x, _y], [end_x, _y]) draws a line

        if self.sketch is None:
            self.sketch = self.getSketch(0)
            self.sketch.OpenSketch()
        
        # A = self.sketch.CreateVertex(startxy[0], startxy[1])
        # print(A.GetX(), A.GetY())
        # B = self.sketch.CreateVertex(endxy[0], endxy[1])
        # print(A.GetX(), A.GetY())
        line = self.sketch.CreateLine(startxy[0],startxy[1],endxy[0],endxy[1])
        if returnVertexName==False:
            return [line]
        else:
            return [line], [A,B]

    def drawArc(self, centerxy, startxy, endxy, returnVertexName=False):
        
        if self.sketch is None:
            self.sketch = self.getSketch(0)
            self.sketch.OpenSketch()
        
        # A = self.sketch.CreateVertex(startxy[0], startxy[1])
        # B = self.sketch.CreateVertex(endxy[0], endxy[1])
        # C = self.sketch.CreateVertex(centerxy[0], centerxy[1])
        arc = self.sketch.CreateArc(centerxy[0], centerxy[1],
                                    startxy[0], startxy[1],
                                    endxy[0], endxy[1])
        if returnVertexName==False:
            return [arc]
        else:
            return [arc], [A,B,C]

    def drawCircle(self, centerxy, radius, returnVertexName=False):
        
        if self.sketch is None:
            self.sketch = self.getSketch(0)
            self.sketch.OpenSketch()
        
        # A = self.sketch.CreateVertex(centerxy[0], centerxy[1])
        arc = self.sketch.CreateCircle(centerxy[0], centerxy[1], radius)

        if returnVertexName==False:
            return [arc]
        else:
            return [arc], [A]

    def checkGeomApp(self):
        if self.geomApp is None:
            self.app.LaunchGeometryEditor()
            self.geomApp = self.app.CreateGeometryEditor(True)
            self.doc = self.geomApp.NewDocument()                
        geomApp = self.geomApp
        return geomApp
    
    def getSketch(self, sketchName, color=None):

        if sketchName in self.sketchNameList:
            self.sketch = self.ass.GetItem(sketchName)
            # open sketch for drawing (must be closed before switch to another sketch)
            self.sketch.OpenSketch()
            return self.sketch
        else:
            self.sketchNameList.append(sketchName)

        self.geomApp = self.checkGeomApp()
        self.doc = self.geomApp.GetDocument()
        self.ass = self.doc.GetAssembly()
        ref1 = self.ass.GetItem('XY Plane')
        ref2 = self.doc.CreateReferenceFromItem(ref1)
        self.sketch = self.ass.CreateSketch(ref2)
        self.sketch.SetProperty('Name', sketchName)

        if color is not None:
            self.sketch.SetProperty('Color', color)
        
        # open sketch for drawing (must be closed before switch to another sketch)
        self.sketch.OpenSketch()
        return self.sketch
    
    def prepareSection(self, token, bMirrorMerge=True, bRotateMerge=True, **kwarg): # csToken is a list of cross section's token

        list_regions = token['list_regions']

        list_region_objects = []
        for idx, list_segments in enumerate(list_regions):
            # print(list_segments)

            # Region
            self.doc.GetSelection().Clear()
            for segment in list_segments:
                # print(segment)
                # debugging = list_segments(i).GetName()
                self.doc.GetSelection().Add(self.sketch.GetItem(segment.GetName()))

            self.sketch.CreateRegions()
            # self.sketch.CreateRegionsWithCleanup(EPS, True) # StatorCore will fail
            if idx == 0:
                region_object = self.sketch.GetItem('Region') # This is how you get access to the region you create.
            else:
                region_object = self.sketch.GetItem('Region.%d'%(idx+1)) # This is how you get access to the region you create.
            list_region_objects.append(region_object)

        # remove region
        if 'list_regions_to_remove' in token.keys():
            list_regions_to_remove = token['list_regions_to_remove']        
            for region_object, boo in zip(list_region_objects, list_regions_to_remove):
                if boo:
                    self.doc.GetSelection().Clear()
                    self.doc.GetSelection().Add(region_object)
                    self.doc.GetSelection().Delete()

        for idx, region_object in enumerate(list_region_objects):
            # Mirror
            if self.bMirror == True:
                if self.edge4Ref is None:
                    self.regionMirrorCopy(region_object, edge4Ref=None, symmetryType=2, bMerge=bMirrorMerge) # symmetryType=2 means x-axis as ref
                else:
                    self.regionMirrorCopy(region_object, edge4Ref=self.edge4ref, symmetryType=None, bMerge=bMirrorMerge) # symmetryType=2 means x-axis as ref

            # RotateCopy
            if self.iRotateCopy != 0:
                # print('Copy', self.iRotateCopy)
                self.regionCircularPattern360Origin(region_object, self.iRotateCopy, bMerge=bRotateMerge)


        self.sketch.CloseSketch()
        return list_region_objects
    
    def regionMirrorCopy(self, region, edge4Ref=None, symmetryType=None, bMerge=True):
        # Default: edge4ref=None, symmetry_type=None, bMerge=True

        mirror = self.sketch.CreateRegionMirrorCopy()
        mirror.SetProperty('Merge', bMerge)
        ref2 = self.doc.CreateReferenceFromItem(region)
        mirror.SetPropertyByReference('Region', ref2)
        
        if edge4Ref is None:
            if symmetryType is None:
                raise Exception('At least give one of edge4ref and symmetry_type')
            else:
                mirror.SetProperty('SymmetryType', symmetryType)
        else:
            ref1 = self.sketch.GetItem(edge4Ref.GetName()) # e.g., u"Line"
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            mirror.SetPropertyByReference('Symmetry', ref2)

        if bMerge == False and region.GetName() == 'Region':
            new_region = self.ass.GetItem('Region.1')
        # return new_region 

    def regionCircularPattern360Origin(self, region, Q_float, bMerge=True):
        # index is used to define name of region

        Q_float = float(Q_float) # don't ask me, ask JSOL

        circular_pattern = self.sketch.CreateRegionCircularPattern()
        circular_pattern.SetProperty('Merge', bMerge)

        ref2 = self.doc.CreateReferenceFromItem(region)
        circular_pattern.SetPropertyByReference('Region', ref2)
        face_region_string = circular_pattern.GetProperty('Region')

        # else:
        #     face_region_string = circular_pattern.GetProperty('Region.%d'%(index+1))
        # %face_region_string = face_region_string[0]
        
        # 想办法避免调用这个函数，比如你可以把绕组变成两个part，一个是上层，一个是下层。
        # if do_you_have_region_in_the_mirror == true
        
        # if True:
        #     # origin_is = origin.GetName()
        #     # ref1 = self.ass.GetItem(self.sketch.GetName()).GetItem('Vertex.3')
        #     origin = self.sketch.CreateVertex(0,0)
        #     ref1 = self.ass.GetItem(self.sketch.GetName()).GetItem(origin.GetName())
        #     ref2 = self.doc.CreateReferenceFromItem(ref1)
        #     circular_pattern.SetPropertyByReference('Center', ref2)
        # elif True:
        #     # Matlab's actxserver cannot pass integer to JMAG (the following 1)
        #     circular_pattern.SetProperty('CenterType', 1)
        #     circular_pattern.SetProperty('CenterPosX', 2.0)
        #     circular_pattern.SetProperty('CenterPosY', 5.0)
        # else:
        # Matlab's actxserver cannot pass integer to JMAG (the following 2)
        circular_pattern.SetProperty('CenterType', 2) # origin I guess

        # print('Copy', Q_float)
        circular_pattern.SetProperty('Angle', '360/%d'% Q_float)
        circular_pattern.SetProperty('Instance', str(Q_float))

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 分析
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

    def draw_jmag_model(self, app, individual_index, im_variant, model_name, bool_trimDrawer_or_vanGogh=True, doNotRotateCopy=False):

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

    @staticmethod
    def run_study(acm_variant, app, study, fea_config_dict, toc):
        logger = logging.getLogger(__name__)
        if fea_config_dict['designer.JMAG_Scheduler'] == False:
            print('[JMAG.py] Run jam.exe...')
            # if run_list[1] == True:
            try:
                study.RunAllCases()
            except Exception as error:
                raise error
            from time import time as clock_time
            msg = 'Time spent on %s is %g s.'%(study.GetName() , clock_time() - toc)
            logger.debug(msg)
            print(msg)
        else:
            print('Submit to JMAG_Scheduler...')
            job = study.CreateJob()
            job.SetValue("Title", study.GetName())
            job.SetValue("Queued", True)
            job.Submit(False) # Fallse:CurrentCase, True:AllCases
            logger.debug('Submit %s to queue (Tran2TSS).'%(acm_variant.individual_name))
            # wait and check
            # study.CheckForCaseResults()
        app.Save()

    @staticmethod
    def mesh_study(acm_variant, app, model, study, output_dir):

        # this `if' judgment is effective only if JMAG-DeleteResultFiles is False 
        # if not study.AnyCaseHasResult(): 
        # mesh
        # acm_variant.add_mesh(study, model)
        if True:
            # this is for multi slide planes, which we will not be using
            refarray = [[0 for i in range(2)] for j in range(1)]
            refarray[0][0] = 3
            refarray[0][1] = 1
            study.GetMeshControl().GetTable("SlideTable2D").SetTable(refarray) 

            study.GetMeshControl().SetValue("MeshType", 1) # make sure this has been exe'd: study.GetCondition(u"RotCon").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)
            study.GetMeshControl().SetValue("RadialDivision", 4) # for air region near which motion occurs
            study.GetMeshControl().SetValue("CircumferentialDivision", 720) #1440) # for air region near which motion occurs 这个数足够大，sliding mesh才准确。
            study.GetMeshControl().SetValue("AirRegionScale", 1.05) # [Model Length]: Specify a value within the following area. (1.05 <= value < 1000)
            study.GetMeshControl().SetValue("MeshSize", 4) # mm
            study.GetMeshControl().SetValue("AutoAirMeshSize", 0)
            study.GetMeshControl().SetValue("AirMeshSize", 4) # mm
            study.GetMeshControl().SetValue("Adaptive", 0)

            # This is not neccessary for whole model FEA. In fact, for BPMSM simulation, it causes mesh error "The copy target region is not found".
            # study.GetMeshControl().CreateCondition("RotationPeriodicMeshAutomatic", "autoRotMesh") # with this you can choose to set CircumferentialDivision automatically

            study.GetMeshControl().CreateCondition("Part", "MagnetMeshCtrl")
            study.GetMeshControl().GetCondition("MagnetMeshCtrl").SetValue("Size", acm_variant.template.fea_config_dict['designer.meshSize_Magnet'])
            study.GetMeshControl().GetCondition("MagnetMeshCtrl").ClearParts()
            study.GetMeshControl().GetCondition("MagnetMeshCtrl").AddSet(model.GetSetList().GetSet("MagnetSet"), 0)

            study.GetMeshControl().CreateCondition("Part", "ShaftMeshCtrl")
            study.GetMeshControl().GetCondition("ShaftMeshCtrl").SetValue("Size", 10) # 10 mm
            study.GetMeshControl().GetCondition("ShaftMeshCtrl").ClearParts()
            study.GetMeshControl().GetCondition("ShaftMeshCtrl").AddSet(model.GetSetList().GetSet("ShaftSet"), 0)

            def mesh_all_cases(study):
                numCase = study.GetDesignTable().NumCases()
                for case in range(0, numCase):
                    study.SetCurrentCase(case)
                    if study.HasMesh() == False:
                        study.CreateMesh()
                    # if case == 0:
                    #     app.View().ShowAllAirRegions()
                    #     app.View().ShowMeshGeometry()
                    #     app.View().ShowMesh()
            mesh_all_cases(study)

        # Export Image
        app.View().ShowAllAirRegions()
        # app.View().ShowMeshGeometry() # 2nd btn
        app.View().ShowMesh() # 3rn btn
        app.View().Zoom(3)
        app.View().Pan(-acm_variant.template.d['GP']['mm_r_or'].value, 0)
        app.ExportImageWithSize(output_dir + model.GetName() + '.png', 2000, 2000)
        app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.

    ''' BELOW is for JMAG Designer
    '''
    def draw_spmsm(self, acm_variant, bool_pyx=False):

        # blue
        # color_rgb_A = np.array([113, 142, 164])/255
        # color_rgb_B = np.array([73, 109, 137])/255

        # yellow
        # color_rgb_A = np.array([255, 252, 170])/255
        # color_rgb_B = np.array([212, 208, 166])/255

        # gray
        color_rgb_A = np.array([236,236,236])/255
        color_rgb_B = np.array([226,226,226])/255

        # Rotor Core
        list_regions_1 = acm_variant.rotorCore.draw(self)
        self.bMirror = False
        self.iRotateCopy = acm_variant.rotorCore.p*2
        region1 = self.prepareSection(list_regions_1, color=color_rgb_A)

        # Shaft
        if not bool_pyx:
            list_regions = acm_variant.shaft.draw(self)
            self.bMirror = False
            self.iRotateCopy = 1
            region0 = self.prepareSection(list_regions)

        # Rotor Magnet
        list_regions = acm_variant.rotorMagnet.draw(self)
        self.bMirror = False
        self.iRotateCopy = acm_variant.rotorMagnet.notched_rotor.p*2
        region2 = self.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

        # This is only for post-processing and it is for handle a un-fixable filling bug with PyX.
        if bool_pyx:
            region1 = self.prepareSection(list_regions_1, color=color_rgb_A)


        # Sleeve
        if not bool_pyx:
            list_regions = acm_variant.sleeve.draw(self)
            self.bMirror = False
            self.iRotateCopy = acm_variant.rotorMagnet.notched_rotor.p*2
            regionS = self.prepareSection(list_regions)

        # Stator Core
        list_regions = acm_variant.stator_core.draw(self)
        self.bMirror = True
        self.iRotateCopy = acm_variant.stator_core.Q
        region3 = self.prepareSection(list_regions, color=color_rgb_A)

        if not bool_pyx:
            # Stator Winding
            list_regions = acm_variant.coils.draw(self)
            self.bMirror = False
            self.iRotateCopy = acm_variant.coils.stator_core.Q
            region4 = self.prepareSection(list_regions)

            # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
            EX = acm_variant.template.d['EX']
            CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * EX['fill_factor'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
            CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                # try:
                #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                # except AttributeError:
                #     # print(EX['wily'])
                #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily']['number_parallel_branch']
                #     print("[inner_rotor_motor.py] Reproduce design using jsonpickle will encounter error here: 'dict' object has no attribute 'number_parallel_branch', implying that the object wily has become a dict after jsonpickle.")
                #     # quit() 

            # Maybe there is a bug here... regarding the excitation for suspension winding...
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
            EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
            EX['DriveW_CurrentAmp'] = acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            EX['BeariW_CurrentAmp'] = acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            # acm_variant.spec_geometry_dict['DriveW_CurrentAmp'] = acm_variant.DriveW_CurrentAmp

            slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
            print('[inner_rotor_motor.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

            # print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
            # print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
            # print('---acm_variant.DriveW_CurrentAmp =', acm_variant.DriveW_CurrentAmp)
            # print('---acm_variant.BeariW_CurrentAmp =', acm_variant.BeariW_CurrentAmp)
            # print('---TORQUE_CURRENT_RATIO:', acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'])
            # print('---SUSPENSION_CURRENT_RATIO:', acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            # Import Model into Designer
            self.save(acm_variant.name, self.show(acm_variant, toString=True))

        return True

    ''' JMAG Description
    '''
    def show(self, acm_variant, toString=False):
        def get_tuple_list(object):
            attrs = list(vars(object).items())
            key_list = [el[0] for el in attrs]
            val_list = [el[1] for el in attrs]
            the_dict = dict(list(zip(key_list, val_list)))
            sorted_key = sorted(key_list, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item)) # this is also useful for string beginning with digiterations '15 Steel'.
            tuple_list = [(key, the_dict[key]) for key in sorted_key]
            return tuple_list
        variant_tuple_list = get_tuple_list(acm_variant)
        template_tuple_list = get_tuple_list(acm_variant.template)
        if toString==False:
            print('- Bearingless PMSM Individual #%s\n\t' % (acm_variant.name), end=' ')
            print(', \n\t'.join("%s = %s" % item for item in variant_tuple_list))
            print(', \n\t'.join("%s = %s" % item for item in template_tuple_list))
            return ''
        else:
            return '\n- Bearingless PMSM Individual #%s\n\t' % (acm_variant.name) \
                    + ', \n\t'.join("%s = %s" % item for item in variant_tuple_list) \
                    + '\n' + '--'*10 + '\n- Template:\n\t' \
                    + ', \n\t'.join("%s = %s" % item for item in template_tuple_list)


if __name__ == '__main__':
    app = win32com.client.Dispatch('designer.Application.171')
    app.Show()
    app.NewProject("Untitled")
    quit()


if __name__ == '__main__':
    from utility import my_execfile
    my_execfile('./default_setting.py', g=globals(), l=locals())
    fea_config_dict

    toolJd = JMAG(fea_config_dict)

    project_name          = 'proj%d'%(0)
    expected_project_file_path = './' + "%s.jproj"%(project_name)

    toolJd.open(expected_project_file_path)

    # toolJd.getSketch('RotorCore', '#FE840E')
    # toolJd.iRotateCopy = 0 # rotorCore.Qr
    # comp1.make(toolJd,toolJd)

    # toolJd.getSketch('RotorBar', '#0E001E')
    # toolJd.iRotateCopy = 0
    # makeToken = comp2.make(toolJd,toolJd)
    
    # # Import Model into Designer
    # toolJd.doc.SaveModel(false) # True: Project is also saved. 
    # model = toolJd.app.GetCurrentModel()
    # model.SetName('model_temp')
    # model.SetDescription('eMach IM Tutorial')

    # # Pre-process
    # toolJd.preProcess(makeToken)
    # # model.CloseCadLink()

