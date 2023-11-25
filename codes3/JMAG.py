import win32com.client, os, logging, utility, numpy
from pylab import np, plt, mpl
print('The mpl backend is', mpl.rcParams['backend'])
print('The mpl backend is', mpl.rcParams['backend'])
print('The mpl backend is', mpl.rcParams['backend'])
mpl.use('Agg') # ('pdf') #   # https://github.com/matplotlib/matplotlib/issues/21950
EPS=0.01 # mm

class data_manager(object):

    def __init__(self):
        self.basic_info = []
        self.time_list = []
        self.TorCon_list = []
        self.ForConX_list = []
        self.ForConY_list = []
        self.ForConAbs_list = []

        self.jmag_loss_list = None
        self.femm_loss_list = None

    def unpack(self, bool_more_info=False):
        if bool_more_info:
            return self.basic_info, self.time_list, self.TorCon_list, self.ForConX_list, self.ForConY_list, self.ForConAbs_list, \
                self.DisplacementAngle_list, \
                self.circuit_current(which='GroupACU'), \
                self.circuit_current(which='GroupACV'), \
                self.circuit_current(which='GroupACW'), \
                self.circuit_current(which='GroupBDU'), \
                self.circuit_current(which='GroupBDV'), \
                self.circuit_current(which='GroupBDW'), \
                self.terminal_voltage(which='GroupACU'), \
                self.terminal_voltage(which='GroupACV'), \
                self.terminal_voltage(which='GroupACW'), \
                self.terminal_voltage(which='GroupBDU'), \
                self.terminal_voltage(which='GroupBDV'), \
                self.terminal_voltage(which='GroupBDW'), \
                self.coil_fluxLinkage(which='GroupACU'), \
                self.coil_fluxLinkage(which='GroupACV'), \
                self.coil_fluxLinkage(which='GroupACW'), \
                self.coil_fluxLinkage(which='GroupBDU'), \
                self.coil_fluxLinkage(which='GroupBDV'), \
                self.coil_fluxLinkage(which='GroupBDW')
        else:
            return self.basic_info, self.time_list, self.TorCon_list, self.ForConX_list, self.ForConY_list, self.ForConAbs_list

    def terminal_voltage(self, which='GroupBDW'): # 2A 2B 2C 4A 4B 4C
        return self.Current_dict['Terminal%s [Case 1]'%(which)]
        # 端点电压是相电压吗？应该是，我们在中性点设置了地电位

    def coil_fluxLinkage(self, which='GroupBDW'):
        return self.FluxLinkage_dict['CircuitCoil%s'%(which)]

    def circuit_current(self, which='GroupBDW'): # 2A 2B 2C 4A 4B 4C
        return self.Current_dict['CircuitCoil%s'%(which)]

    def get_voltage_and_current(self, number_of_steps_at_steady_state):

        # key = 'GroupBDW' if DPNV_or_SEPA == True else 'Torque'
        key = 'Default'
        try:
            # 4C <- the C-phase of the 4 pole winding
            mytime  = self.Current_dict['Time(s)'][-number_of_steps_at_steady_state:]
            voltage =      self.terminal_voltage()[-number_of_steps_at_steady_state:]
            current =       self.circuit_current(which=key)[-number_of_steps_at_steady_state:]
        except KeyError as error:
            raise error

        # if len(mytime) > len(voltage):
        #     mytime = mytime[:len(voltage)]

        # print len(mytime), len(voltage), number_of_steps_at_steady_state
        # print len(mytime), len(voltage)
        # print len(mytime), len(voltage)

        # for access to plot
        self.myvoltage = voltage
        self.mycurrent = current
        self.mytime    = mytime

    def power_factor(self, number_of_steps_at_steady_state, targetFreq=1e3, numPeriodicalExtension=1000):
        # number_of_steps_at_steady_state: steps corresponding to half the period 

        # for key, val in self.Current_dict.iteritems():
        #     if 'Terminal' in key:
        #         print key, val
        # quit()

        self.get_voltage_and_current(number_of_steps_at_steady_state)
        mytime  = self.mytime
        voltage = self.myvoltage
        current = self.mycurrent

        # from pylab import *
        # print len(mytime), len(voltage), len(current)
        # figure()
        # plot(mytime, voltage)
        # plot(mytime, current)
        # show()
        power_factor, u, i, phase_diff_ui = utility.compute_power_factor_from_half_period(voltage, current, mytime, targetFreq=targetFreq, numPeriodicalExtension=numPeriodicalExtension)
        self.ui_info = [power_factor, u, i, phase_diff_ui]
        return power_factor


class JMAG(object): #< ToolBase & DrawerBase & MakerExtrnudeBase & MakerRevolveBase
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
        # self.hide_or_show = None

    def open(self, expected_project_file_path):
        if self.app is None:
            try:
                # app = win32com.client.Dispatch('designer.Application.200')
                app = win32com.client.Dispatch('designer.Application.171')
                print('JMAG 17.1 is not found. Will use any other JMAG version avaiilable.')
                # app = win32com.client.gencache.EnsureDispatch('designer.Application.171') # https://stackoverflow.com/questions/50127959/win32-dispatch-vs-win32-gencache-in-python-what-are-the-pros-and-cons
            except:
                try:
                    app = win32com.client.Dispatch('designer.Application')
                    # app = win32com.client.Dispatch('designer.Application.181')
                    # app = win32com.client.gencache.EnsureDispatch('designer.Application.171')
                except:
                    raise Exception('No JMAG Designer 17 is found in this PC.')

            self.JMAG_version_string = app.VersionString(0)
            self.JMAG_version_number = float(app.VersionString(0)[:2])

            if self.fea_config_dict['designer.Show'] == True:
                app.Show()
            else:
                app.Hide()
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
        print('[JMAG.py] expected_project_file_path (to abs path):', os.path.abspath(expected_project_file_path))
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
        app.SaveAs(os.path.abspath(expected_project_file_path)) # must be absolute path!
        logger = logging.getLogger(__name__)
        logger.info(r'Create JMAG project file: %s'%(expected_project_file_path))
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
    def pre_process_CPPM(self, app, model, acm_variant):
        # pre-process : you can select part by coordinate!
        ''' Group '''
        def group(name, id_list):
            model.GetGroupList().CreateGroup(name)
            for the_id in id_list:
                model.GetGroupList().AddPartToGroup(name, the_id)
                # model.GetGroupList().AddPartToGroup(name, name) #<- this also works

        part_ID_list = model.GetPartIDs()

        # view = app.View()
        # view.ClearSelect()
        # sel = view.GetCurrentSelection()
        # sel.SelectPart(123)
        # sel.SetBlockUpdateView(False)
        SI = acm_variant.template.spec_input_dict
        p = SI['p']
        s = SI['no_segmented_magnets']
        Q = SI['Qs']
        print(Q)
        print(int(1 + 1 + p*1*s + 1 + Q*2))
        # quit()
                                #   轴 转子 永磁体  定子 绕组
        if len(part_ID_list) != int(1 + 1 + p*1*s + 1 + Q*2):
        # if len(part_ID_list) != int(1 + 1 + p*1*s +      1 + Q*2):
        # if len(part_ID_list) != int(1 + Q*2 + 1):
            msg = 'Number of Parts is unexpected. Should be %d but get %d.\n'%(int(1 + 1 + p*2*s + 1 + 1 + Q*2), len(part_ID_list)) + self.show(toString=True)
            logger = logging.getLogger(__name__)
            logger.error(msg)
            raise utility.ExceptionBadNumberOfParts(msg)

        self.id_rotorCore = id_rotorCore = part_ID_list[0]
        id_shaft = part_ID_list[1]
        partIDRange_Magnet = part_ID_list[2:int(2+p*s*1)]                           # 此处C.P.应为1*p
        # id_sleeve = part_ID_list[int(2+p*s*2)]
        self.id_statorCore = id_statorCore = part_ID_list[int(2+p*s*1)]
        partIDRange_Coil = part_ID_list[int(2+p*s*1) + 1 : int(2+p*s*1)+ 1 + int(Q*2)]

        # debug
        # print(p)
        # print(int(s))
        # print(id_rotorCore)
        # print(id_shaft)
        # print(partIDRange_Magnet)
        # print(len(partIDRange_Magnet))
        # # print(id_sleeve)
        # print(id_statorCore)
        # print(partIDRange_Coil)
        # print(len(partIDRange_Coil))
        # quit()

        self.bool_suppressShaft = False
        # model.SuppressPart(id_sleeve, 1)

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
        # edge_set(u"AirGapCoast", 0, self.template.d['GP']['mm_r_ro'].value+0.5*self.Length_AirGap)

        # Shaft
        add_part_to_set('ShaftSet', 0.0, 0.0, ID=id_shaft) # 坐标没用，不知道为什么，而且都给了浮点数了

        # Create Set for 8 poles Winding
        Angle_StatorSlotSpan = 360/Q
        # R = self.mm_r_si + self.mm_d_stt + self.mm_d_st *0.5 # this is not generally working (JMAG selects stator core instead.)
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
        R = GP['mm_r_si'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value - 0.5*GP['mm_d_pm'].value
        # alpha_rs = GP['deg_alpha_rs'].value /180*np.pi
        deg_pole_span = 360 / (p*2)

        # 这里暂时进不去
        if s>1:
            deg_alpha_notch  = (GP['deg_alpha_rm'].value - s*GP['deg_alpha_rs'].value) / (s-1) # inter-segment notch占的角度
            alpha_notch = deg_alpha_notch /180*np.pi
        # 从这开始
        list_xy_magnets = []
        # list_xy_airWithinRotorSlot = []
        for ind in range(int(p)):
            natural_ind = ind + 1
            
            if s==1:
                      # v---This negative sign means we walk CCW to assign sets.
                THETA = - (180/p-GP['deg_alpha_rm'].value + 0.5*GP['deg_alpha_rm'].value + deg_pole_span*2*ind) /180.*np.pi
                X = R*np.cos(THETA)
                Y = R*np.sin(THETA)

                add_part_to_set("Magnet %d"%(natural_ind), X, Y)
                list_xy_magnets.append([X,Y])


            ### 这里进不去（暂时）

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

        part_list_set('Motion_Region', list_xy_magnets, list_part_id=[id_rotorCore, id_shaft])

        part_list_set('MagnetSet', list_xy_magnets)

        
        # debug
        # print(p)
        # print(int(s))
        # print(id_rotorCore)
        # print(id_shaft)
        # print(partIDRange_Magnet)
        # print(len(partIDRange_Magnet))
        # # print(id_sleeve)
        # print(id_statorCore)
        # print(partIDRange_Coil)
        # print(part_ID_list)
        # print(part_ID_list)
        # print(part_ID_list)
        # print(part_ID_list)
        # print(part_ID_list)
        # print(len(partIDRange_Coil))
        # quit()
        
        return True
    def pre_process_fluxAlternator(self, app, model, acm_variant):
        # pre-process : you can select part by coordinate!
        ''' Group '''
        def group(name, id_list):
            model.GetGroupList().CreateGroup(name)
            for the_id in id_list:
                model.GetGroupList().AddPartToGroup(name, the_id)
                # model.GetGroupList().AddPartToGroup(name, name) #<- this also works

        part_ID_list = model.GetPartIDs()
        print(part_ID_list)

        SI = acm_variant.template.spec_input_dict
        p  = SI['p']
        pe = SI['pe']
        Q  = SI['Qs']
                                #   转子 轴 永磁体 定子 绕组
        if len(part_ID_list) != int(1 + 1 + pe*2 + 1 + Q*4):
            msg = 'Number of Parts is unexpected. Should be %d but get %d.\n'%(int(1 + 1 + pe*2 + 1 + Q*4), len(part_ID_list)) + self.show(toString=True)
            logger = logging.getLogger(__name__)
            logger.error(msg)
            raise utility.ExceptionBadNumberOfParts(msg)


        self.id_rotorCore = id_rotorCore = part_ID_list[0]
        id_shaft           = part_ID_list[1]
        partIDRange_Magnet = part_ID_list[2:int(2+pe*2)]
        self.id_statorCore = id_statorCore = part_ID_list[  int(2+pe*2)]
        partIDRange_Coil   = part_ID_list[  int(2+pe*2)+1 : int(2+pe*2)+1 + int(Q*4)]

        # debug
        # print(id_rotorCore)
        # print(id_shaft)
        # print(partIDRange_Magnet)
        # print(id_statorCore)
        # print(partIDRange_Coil)

        # model.SuppressPart(id_shaft, 1) # cause internal error if you suppress parts that relate to condition settings
        self.bool_suppressShaft = False

        # group("Magnet", partIDRange_Magnet)
        group("Magnet-CW",  [partIDRange_Magnet[0]])
        group("Magnet-CCW", [partIDRange_Magnet[1]])
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
        # edge_set(u"AirGapCoast", 0, self.template.d['GP']['mm_r_ro'].value+0.5*self.Length_AirGap)

        # Shaft
        add_part_to_set('ShaftSet', 0.0, 0.0, ID=id_shaft) # 坐标没用，不知道为什么，而且都给了浮点数了


        wily_not_used = acm_variant.template.d['EX']['wily']

        # Create Set for 4 poles Winding
        countXL = 0
        for coil in acm_variant.coils.wily:
            countXL += 1 
            UVW = coil['LayerX-Phase']
            UpDown = coil['LayerX-Direction']
            X = coil['LayerX-X']
            Y = coil['LayerX-Y']
            add_part_to_set("CoilLX%s%s %d"%(UVW,UpDown,countXL), X, Y)

        # Create Set for 2 poles Winding
        countYL = 0
        for coil in acm_variant.coils.wily:
            countYL += 1 
            UVW = coil['LayerY-Phase']
            UpDown = coil['LayerY-Direction']
            X = coil['LayerY-X']
            Y = coil['LayerY-Y']
            add_part_to_set("CoilLY%s%s %d"%(UVW,UpDown,countYL), X, Y)

        # Create Set for Magnets
        GP = acm_variant.template.d['GP']

        list_xy_magnets = []
        # for ind in range(int(p*2)):
        counterMagnet = 0
        for magnet in acm_variant.statorMagnet.maly:
            counterMagnet += 1
                    # v---This negative sign means we walk CCW to assign sets.
            # THETA = - (180/p-GP['deg_alpha_rm'].value + 0.5*GP['deg_alpha_rm'].value + deg_pole_span*ind) /180.*np.pi
            # X = R*np.cos(THETA)
            # Y = R*np.sin(THETA)
            X = magnet['X']
            Y = magnet['Y']
            add_part_to_set("Magnet %d"%(counterMagnet), X, Y)
            list_xy_magnets.append([X,Y])

        # Create Set for Motion Region
        def add_parts_to_set(name, list_xy, list_part_id=None, prefix=None):
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

        if 'PMSM' in acm_variant.template.name:
            add_parts_to_set('Motion_Region', list_xy_magnets, list_part_id=[id_rotorCore, id_shaft])
        elif 'Flux_Alternator' in acm_variant.template.name:
            add_parts_to_set('Motion_Region', [], list_part_id=[id_rotorCore, id_shaft])

        add_parts_to_set('MagnetSet', list_xy_magnets)
        return True
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

        self.id_rotorCore = id_rotorCore = part_ID_list[0]
        id_shaft = part_ID_list[1]
        partIDRange_Magnet = part_ID_list[2:int(2+p*s*2)]
        id_sleeve = part_ID_list[int(2+p*s*2)]
        self.id_statorCore = id_statorCore = part_ID_list[int(2+p*s*2)+1]
        partIDRange_Coil = part_ID_list[int(2+p*s*2)+2 : int(2+p*s*2)+2 + int(Q*2)]

        # debug
        # print(id_rotorCore)
        # print(id_shaft)
        # print(partIDRange_Magnet)
        # print(id_sleeve)
        # print(id_statorCore)
        # print(partIDRange_Coil)

        self.bool_suppressShaft = False
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
        # edge_set(u"AirGapCoast", 0, self.template.d['GP']['mm_r_ro'].value+0.5*self.Length_AirGap)

        # Shaft
        add_part_to_set('ShaftSet', 0.0, 0.0, ID=id_shaft) # 坐标没用，不知道为什么，而且都给了浮点数了

        # Create Set for 4 poles Winding
        Angle_StatorSlotSpan = 360/Q
        # R = self.mm_r_si + self.mm_d_stt + self.mm_d_st *0.5 # this is not generally working (JMAG selects stator core instead.)
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
        R = GP['mm_r_si'].value - GP['mm_d_sleeve'].value - GP['mm_d_mech_air_gap'].value - 0.5*GP['mm_d_pm'].value
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
        part_list_set('Motion_Region', list_xy_magnets, list_part_id=[id_rotorCore, id_shaft])

        part_list_set('MagnetSet', list_xy_magnets)
        return True
    def pre_process_FSPM(self, app, model, acm_variant):

        ''' 1. Group: make multiple parts that have the same name as a group'''
        def group(name, id_list):
            model.GetGroupList().CreateGroup(name)
            for the_id in id_list:
                # model.GetGroupList().AddPartToGroup(name, name) #<- this also works
                model.GetGroupList().AddPartToGroup(name, the_id)

        part_ID_list = model.GetPartIDs()
        SI = acm_variant.template.spec_input_dict
        p = SI['p']
        pe = SI['pe']
        Qs = SI['Qs']
                                   #   轴 转子 永磁体  定子 绕组
        self.number_of_parts_in_JMAG = 1 + 1 + pe*2 + 1 + Qs*2

        if len(part_ID_list) != self.number_of_parts_in_JMAG:
            msg = 'Number of Parts is unexpected. Should be %d but get %d.\n'%(int(self.number_of_parts_in_JMAG), len(part_ID_list))
            #  + self.show(toString=True)
            logger = logging.getLogger(__name__)
            logger.error(msg)
            raise utility.ExceptionBadNumberOfParts(msg)

        self.id_rotorCore = id_rotorCore = part_ID_list[0]
        id_shaft = part_ID_list[1]
        partIDRange_Magnet = part_ID_list[2:int(2+pe*2)]
        self.id_statorCore = id_statorCore = part_ID_list[int(2+pe*2)]
        partIDRange_Coil = part_ID_list[int(2+pe*2)+1 : int(2+pe*2)+1 + int(Qs*2)]

        # group("Magnet", partIDRange_Magnet)
        # group("Magnet-CW",  partIDRange_Magnet[0::2])
        # group("Magnet-CCW", partIDRange_Magnet[1::2])
        group("Magnet-CCW",  partIDRange_Magnet[0::2]) # for reversing torque sign
        group("Magnet-CW", partIDRange_Magnet[1::2])
        group("Coils", partIDRange_Coil)


        ''' 2. Add Part to Set for later references '''
        def add_part_to_set(name, x, y, ID=None):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection()
            if ID is None:
                sel.SelectPartByPosition(x,y,0) # z=0 for 2D
            else:
                sel.SelectPart(ID)
            model.GetSetList().GetSet(name).AddSelected(sel)

        # Shaft Set
        # add_part_to_set('ShaftSet', 0.0, 0.0, ID=id_shaft) # 坐标没用，不知道为什么，而且都给了浮点数了
        self.bool_suppressShaft = True
        if self.bool_suppressShaft:
            model.SuppressPart(id_shaft, 1)

        # Create Sets for Coils
        wily = acm_variant.template.d['EX']['wily']
        Angle_StatorSlotSpan = 360/Qs
        R = np.sqrt(acm_variant.coils.PCoil[0]**2 + acm_variant.coils.PCoil[1]**2) # switch coil naming for intuitive concentric winding (cf. PMSM)

        # Coils belonging to Layer X
        THETA = np.arctan(-acm_variant.coils.PCoil[1]/acm_variant.coils.PCoil[0]) + (2*np.pi)/Qs
        X = R*np.cos(THETA)
        Y = R*np.sin(THETA)
        countXL = 0
        for UVW, UpDown in zip(wily.layer_X_phases,wily.layer_X_signs):
            countXL += 1 
            add_part_to_set("CoilLX%s%s %d"%(UVW,UpDown,countXL), X, Y)

            # print(X, Y, THETA)
            THETA += Angle_StatorSlotSpan/180.*np.pi
            X = R*np.cos(THETA)
            Y = R*np.sin(THETA)

        # Coils belonging to Layer Y
        THETA = np.arctan(acm_variant.coils.PCoil[1]/acm_variant.coils.PCoil[0]) # switch coil naming for intuitive concentric winding (cf. PMSM)
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
        R = 0.5*(GP['mm_r_si'].value + GP['mm_r_so'].value)
        deg_pole_span = 360 / (pe*2)
        list_xy_magnets = []
        for ind in range(int(pe*2)):
            natural_ind = ind + 1
                  # v---This negative sign means we walk CCW to assign sets.
            THETA = -(deg_pole_span*ind) /180*np.pi
            X = R*np.cos(THETA)
            Y = R*np.sin(THETA)
            add_part_to_set("Magnet %d"%(natural_ind), X, Y)
            list_xy_magnets.append([X,Y])

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

        part_list_set('MagnetSet', list_xy_magnets)

        if self.bool_suppressShaft == False:
            add_part_to_set('ShaftSet', 0.0, 0.0, ID=id_shaft)
            part_list_set('Motion_Region', [], list_part_id=[id_rotorCore, id_shaft])
        else:
            part_list_set('Motion_Region', [], list_part_id=[id_rotorCore, ])

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
        SI = acm_variant.template.SI
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
        self.add_material(study, acm_variant)

        # Conditions - Motion
        study.CreateCondition("RotationMotion", "RotCon") # study.GetCondition(u"RotCon").SetXYZPoint(u"", 0, 0, 1) # megbox warning
        # print('the_speed:', acm_variant.speed_rpm)
        study.GetCondition("RotCon").SetValue("AngularVelocity", int(EX['the_speed']))
        study.GetCondition("RotCon").ClearParts()
        study.GetCondition("RotCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)

        study.GetCondition("RotCon").SetValue(u"InitialRotationAngle", acm_variant.InitialRotationAngle)
        # if acm_variant.Rotation_Axis == -1:
        # study.GetCondition("RotCon").SetValue("Rotation Axis", "DownWard")


        study.CreateCondition("Torque", "TorCon") # study.GetCondition(u"TorCon").SetXYZPoint(u"", 0, 0, 0) # megbox warning
        study.GetCondition("TorCon").SetValue("TargetType", 1)
        study.GetCondition("TorCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("TorCon").ClearParts()

        study.CreateCondition("Force", "ForCon")
        study.GetCondition("ForCon").SetValue("TargetType", 1)
        study.GetCondition("ForCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("ForCon").ClearParts()


        # Conditions - FEM Coils & Conductors (i.e. stator/rotor winding)
        if acm_variant.boolCustomizedCircuit == True:
            acm_variant.add_circuit_customized(app, model, study)
        else:
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

            # FREQUENCY = EX['DriveW_Freq'] # PMSM
            FREQUENCY = SI['ExcitationFreqSimulated'] # FSPM

            if number_cycles_prolonged == 0:
                if number_cycles_in_3rdTSS == 0:
                    refarray = [[0 for i in range(3)] for j in range(3)]
                    refarray[0][0] = 0
                    refarray[0][1] =    1
                    refarray[0][2] =        50
                    refarray[1][0] = number_cycles_in_1stTSS/FREQUENCY
                    refarray[1][1] =    number_of_steps_1stTSS
                    refarray[1][2] =        50
                    refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/FREQUENCY
                    refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                    refarray[2][2] =        50
                else:
                    refarray = [[0 for i in range(3)] for j in range(4)]
                    refarray[0][0] = 0
                    refarray[0][1] =    1
                    refarray[0][2] =        50
                    refarray[1][0] = number_cycles_in_1stTSS/FREQUENCY
                    refarray[1][1] =    number_of_steps_1stTSS
                    refarray[1][2] =        50
                    refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/FREQUENCY
                    refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                    refarray[2][2] =        50
                    refarray[3][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS+number_cycles_in_3rdTSS)/FREQUENCY
                    refarray[3][1] =    number_of_steps_3rdTSS
                    refarray[3][2] =        50
            else:
                refarray = [[0 for i in range(3)] for j in range(4)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = number_cycles_in_1stTSS/FREQUENCY
                refarray[1][1] =    number_of_steps_1stTSS
                refarray[1][2] =        50
                refarray[2][0] = (number_cycles_in_1stTSS+number_cycles_in_2ndTSS)/FREQUENCY
                refarray[2][1] =    number_of_steps_2ndTSS # 最后的number_of_steps_2ndTSS（32）步，必须对应半个周期，从而和后面的铁耗计算相对应。
                refarray[2][2] =        50
                refarray[3][0] = refarray[2][0] + number_cycles_prolonged/FREQUENCY 
                refarray[3][1] =    number_cycles_prolonged*acm_variant.template.fea_config_dict['designer.TranRef-StepPerCycle'] 
                refarray[3][2] =        50
                # refarray[4][0] = refarray[3][0] + 0.5/FREQUENCY # 最后来一个超密的半周期400步
                # refarray[4][1] =    400
                # refarray[4][2] =        50
            number_of_total_steps = 1 + number_of_steps_1stTSS + number_of_steps_2ndTSS + number_of_steps_3rdTSS + number_cycles_prolonged*acm_variant.template.fea_config_dict['designer.TranRef-StepPerCycle'] # [Double Check] don't forget to modify here!
            # print('[inner_rotor_motor.py]: refarray:', refarray)
            DM.GetDataSet("SectionStepTable").SetTable(refarray)
            study.GetStep().SetValue("Step", number_of_total_steps)
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

        # add equations
        study.GetDesignTable().AddEquation("freq")
        study.GetDesignTable().AddEquation("speed")
        study.GetDesignTable().GetEquation("freq").SetType(0)
        study.GetDesignTable().GetEquation("freq").SetExpression("%g"%(FREQUENCY))
        study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency in Hz")
        study.GetDesignTable().GetEquation("speed").SetType(1)
        # study.GetDesignTable().GetEquation("speed").SetExpression("freq * %d"%(60/(0.5*EX['RotorPoleNumber'])))
        if 'number_of_rotor_pole_pairs' not in SI.keys():
            number_of_rotor_pole_pairs = SI['p']
        else:
            number_of_rotor_pole_pairs = SI['number_of_rotor_pole_pairs']
        study.GetDesignTable().GetEquation("speed").SetExpression("freq * %f"%(60 / number_of_rotor_pole_pairs ))
        study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed in r/min")

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
        if self.fea_config_dict['designer.AddIronLossCondition']:
            cond = study.CreateCondition("Ironloss", "IronLossConStator")
            cond.SetValue("RevolutionSpeed", "freq*60/%d"%(0.5*EX['DriveW_poles']))
            cond.SetValue(u"Poles", EX['DriveW_poles'])
            cond.ClearParts()
            sel = cond.GetSelection()
            # sel.SelectPartByPosition(acm_variant.template.d['GP']['mm_r_si'].value + EPS, 0 ,0) # 2022-02-04 这里发现代码有点歧义：注意，实际上acm_variant.template.d['GP']已经被修改了，acm_variant.template.d['GP'] = acm_variant.GP。 # btw, this works!
            sel.SelectPart(self.id_statorCore)
            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results and to have a FFT plot
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3) # 3:Custom
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            # cond.SetValue("StartReferenceStep", number_of_total_steps+1 - 0.5*int(number_of_steps_2ndTSS/number_cycles_in_2ndTSS)) # 1/4 period <=> 0.5*number_of_steps_2ndTSS
            if number_cycles_in_2ndTSS < 0.5:
                raise Exception('Invalid number_cycles_in_2ndTSS:', number_cycles_in_2ndTSS)
            cond.SetValue("StartReferenceStep", number_of_total_steps+1 - int(number_of_steps_2ndTSS*0.5/number_cycles_in_2ndTSS)) # 1/2 period <=> number_of_steps_2ndTSS
                # cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTSS*0.5) # 1/4 period <=> number_of_steps_2ndTSS*0.5
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 2) # specify reference steps for 1/2 period and extend it to whole period. Don't use 1/4 period (for reasons, see JMAG help on this topic)
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
        # Rotor
        if self.fea_config_dict['designer.AddIronLossCondition']:
            cond = study.CreateCondition("Ironloss", "IronLossConRotor")
            cond.SetValue("BasicFrequencyType", 2)
            cond.SetValue("BasicFrequency", "freq")
                # cond.SetValue(u"BasicFrequency", u"slip*freq") # this require the signal length to be at least 1/4 of slip period, that's too long!
            cond.ClearParts()
            sel = cond.GetSelection()
            # sel.SelectPartByPosition(acm_variant.mm_r_ri + EPS, 0 ,0) # Why this is not working??? Because it is integer.... you must use 0.0 instead of 0!!!
            sel.SelectPart(self.id_rotorCore)

            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3)
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            cond.SetValue("StartReferenceStep", number_of_total_steps+1 - int(number_of_steps_2ndTSS*0.5/number_cycles_in_2ndTSS)) # 1/2 period <=> number_of_steps_2ndTSS
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 2) # specify reference steps for 1/2 or 1/4 period and extend it to whole period (2 means 1/2 peirodicity; 4 means 1/4 peirodicity)
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 



        # Check CSV reults for iron loss (You cannot check this for Freq study) # CSV and save space
        study.GetStudyProperties().SetValue("CsvOutputPath", dir_csv_output_folder) # it's folder rather than file!
        if acm_variant.template.fea_config_dict["designer.AddIronLossCondition"]:
            # study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss") # old
            # study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;TerminalVoltage;JouleLoss;StoredEnergy;TotalDisplacementAngle;Inductance;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss") # new since 2022
            study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;TerminalVoltage;JouleLoss;StoredEnergy;TotalDisplacementAngle;Inductance;FEMCoilInductance;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss") # new since 2022
        else:
            study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;TerminalVoltage;JouleLoss;StoredEnergy;TotalDisplacementAngle;Inductance") # no iron loss condition
        study.GetStudyProperties().SetValue("DeleteResultFiles", acm_variant.template.fea_config_dict['delete_results_after_calculation'])

        # Terminal Voltage/Circuit Voltage: Check for outputing CSV results 
        if 'Flux_Alternator' in acm_variant.template.name:
            study.GetCircuit().CreateTerminalLabel("TerminalLabel"+acm_variant.circuit_coil_names[0], 9,  1)
            study.GetCircuit().CreateTerminalLabel("TerminalLabel"+acm_variant.circuit_coil_names[1], 9,  6) # seek WriteTable
            study.GetCircuit().CreateTerminalLabel("TerminalLabel"+acm_variant.circuit_coil_names[2], 9, -4)
            study.GetCircuit().CreateTerminalLabel("TerminalLabel"+acm_variant.circuit_coil_names[3], 9, -9)
            print('[JMAG.py] Cannot create terminal label:', "TerminalLabel"+acm_variant.circuit_coil_names[0])
            print('[JMAG.py] Cannot create terminal label:', "TerminalLabel"+acm_variant.circuit_coil_names[1])
            print('[JMAG.py] Cannot create terminal label:', "TerminalLabel"+acm_variant.circuit_coil_names[2])
            print('[JMAG.py] Cannot create terminal label:', "TerminalLabel"+acm_variant.circuit_coil_names[3])
        else:
            if self.JMAG_version_number >= 20:
                # 新版JMAG把Y轴反过来了
                study.GetCircuit().CreateTerminalLabel("TerminalGroupACU", 8, 13) # seek WriteTable
                study.GetCircuit().CreateTerminalLabel("TerminalGroupACV", 8, 11)
                study.GetCircuit().CreateTerminalLabel("TerminalGroupACW", 8, 9)
                study.GetCircuit().CreateTerminalLabel("TerminalGroupBDU", 23, 13)
                study.GetCircuit().CreateTerminalLabel("TerminalGroupBDV", 23, 11)
                study.GetCircuit().CreateTerminalLabel("TerminalGroupBDW", 23, 9)
            else:
                study.GetCircuit().CreateTerminalLabel("TerminalGroupACU", 8, -13) # seek WriteTable
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




        self.study_name = study_name
        return study
    def add_structural_static_study(self):
        pass
    def add_mesh(self, study, model):
        pass
    # TranFEAwi2TSS
    def add_material(self, study, acm_variant):
        acm_template = acm_variant.template

        # if 'PMSM' in acm_variant.template.name:
        #     rotorCoreName = "NotchedRotor"
        # elif 'Flux_Alternator' in acm_variant.template.name:
        #     rotorCoreName = "SalientPoleRotor"
        rotorCoreName = acm_variant.rotorCore.name

        if 'M19' in acm_template.spec_input_dict['Steel']:
            study.SetMaterialByName("StatorCore", "M-19 Steel Gauge-29")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 95)
                # study.GetMaterial(u"Stator Core").SetValue(u"UserConductivityValue", 1900000)

            study.SetMaterialByName(rotorCoreName, "M-19 Steel Gauge-29")
            study.GetMaterial(rotorCoreName).SetValue("Laminated", 1)
            study.GetMaterial(rotorCoreName).SetValue("LaminationFactor", 98)

        elif 'M15' in acm_template.spec_input_dict['Steel']:
            study.SetMaterialByName("StatorCore", "M-15 Steel")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 98)

            study.SetMaterialByName(rotorCoreName, "M-15 Steel")
            study.GetMaterial(rotorCoreName).SetValue("Laminated", 1)
            study.GetMaterial(rotorCoreName).SetValue("LaminationFactor", 98)

        elif acm_template.spec_input_dict['Steel'] == 'Arnon5':
            study.SetMaterialByName("StatorCore", "Arnon5-final")
            study.GetMaterial("StatorCore").SetValue("Laminated", 1)
            study.GetMaterial("StatorCore").SetValue("LaminationFactor", 96)

            study.SetMaterialByName(rotorCoreName, "Arnon5-final")
            study.GetMaterial(rotorCoreName).SetValue("Laminated", 1)
            study.GetMaterial(rotorCoreName).SetValue("LaminationFactor", 96)

        elif acm_template.spec_input_dict['Steel'] == '35CS250':
            study.SetMaterialByName(u"StatorCore", u"35CS250")
            study.SetMaterialByName(rotorCoreName, u"35CS250")

        else:
            msg = 'Warning: default material is used: DCMagnetic Type/50A1000.'
            print(msg)
            logging.getLogger(__name__).warn(msg)
            study.SetMaterialByName("StatorCore", "DCMagnetic Type/50A1000")
            study.GetMaterial("StatorCore").SetValue("UserConductivityType", 1)
            study.SetMaterialByName(rotorCoreName, "DCMagnetic Type/50A1000")
            study.GetMaterial(rotorCoreName).SetValue("UserConductivityType", 1)

        study.SetMaterialByName("Coils", "Copper")
        study.GetMaterial("Coils").SetValue("UserConductivityType", 1)

        # study.SetMaterialByName("Cage", "Aluminium")
        # study.GetMaterial("Cage").SetValue("EddyCurrentCalculation", 1)
        # study.GetMaterial("Cage").SetValue("UserConductivityType", 1)
        # study.GetMaterial("Cage").SetValue("UserConductivityValue", self.Bar_Conductivity)

        # N40H Reversible
        available_temperature_list = [-40, 20, 60, 80, 100, 120, 150, 180, 200, 220] # according to JMAG
        magnet_temperature = min(available_temperature_list, key=lambda x:abs(x-acm_template.spec_input_dict['Temperature']))        
        print('[inner_rotor_motor.py] magnet_temperature is', magnet_temperature, 'deg C.')
        if 'PMSM' in acm_variant.template.name:
            study.SetMaterialByName(u"Magnet", u"Arnold/Reversible/N40H")
            study.GetMaterial(u"Magnet").SetValue(u"EddyCurrentCalculation", 1)
            study.GetMaterial(u"Magnet").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)

            study.GetMaterial(u"Magnet").SetValue(u"Poles", acm_template.d['EX']['DriveW_poles'])
            study.GetMaterial(u"Magnet").SetDirectionXYZ(1, 0, 0)
            study.GetMaterial(u"Magnet").SetAxisXYZ(0, 0, 1)
            study.GetMaterial(u"Magnet").SetOriginXYZ(0, 0, 0)
            study.GetMaterial(u"Magnet").SetPattern(u"RadialCircular")
            study.GetMaterial(u"Magnet").SetValue(u"StartAngle", 0.5* 360/(2*acm_template.SI['p']) ) # 半个极距

        elif 'FSPM' in acm_variant.template.name:
            study.SetMaterialByName(u"Magnet-CW", u"Arnold/Reversible/N40H")
            study.GetMaterial(      u"Magnet-CW").SetValue(u"EddyCurrentCalculation", 1)
            study.GetMaterial(      u"Magnet-CW").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)
            study.GetMaterial(      u"Magnet-CW").SetOrientation(False)
            study.GetMaterial(      u"Magnet-CW").SetDirectionXYZ(0, 0, 1)
            study.GetMaterial(      u"Magnet-CW").SetAxisXYZ(0, 0, 1)
            study.GetMaterial(      u"Magnet-CW").SetOriginXYZ(0, 0, 0)
            study.GetMaterial(      u"Magnet-CW").SetPattern(u"Circular")

            study.SetMaterialByName(u"Magnet-CCW", u"Arnold/Reversible/N40H")
            study.GetMaterial(      u"Magnet-CCW").SetValue(u"EddyCurrentCalculation", 1)
            study.GetMaterial(      u"Magnet-CCW").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)
            study.GetMaterial(      u"Magnet-CCW").SetOriginXYZ(0, 0, 0)
            study.GetMaterial(      u"Magnet-CCW").SetPattern(u"Circular")

        elif 'Flux_Alternator' in acm_variant.template.name:
            study.SetMaterialByName(u"Magnet-CW", u"Arnold/Reversible/N40H")
            study.GetMaterial(      u"Magnet-CW").SetValue(u"EddyCurrentCalculation", 1)
            study.GetMaterial(      u"Magnet-CW").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)
            study.GetMaterial(      u"Magnet-CW").SetOrientation(False)
            study.GetMaterial(      u"Magnet-CW").SetDirectionXYZ(0, 0, 1)
            study.GetMaterial(      u"Magnet-CW").SetAxisXYZ(0, 0, 1)
            study.GetMaterial(      u"Magnet-CW").SetOriginXYZ(0, 0, 0)
            study.GetMaterial(      u"Magnet-CW").SetPattern(u"Circular")

            study.SetMaterialByName(u"Magnet-CCW", u"Arnold/Reversible/N40H")
            study.GetMaterial(      u"Magnet-CCW").SetValue(u"EddyCurrentCalculation", 1)
            study.GetMaterial(      u"Magnet-CCW").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)
            study.GetMaterial(      u"Magnet-CCW").SetOriginXYZ(0, 0, 0)
            study.GetMaterial(      u"Magnet-CCW").SetPattern(u"Circular")
        
        
        elif 'CPPM' in acm_variant.template.name:

            study.SetMaterialByName(u"Magnet", u"Arnold/Reversible/N40H")
            study.GetMaterial(      u"Magnet").SetValue(u"EddyCurrentCalculation", 1)
            study.GetMaterial(      u"Magnet").SetValue(u"Temperature", magnet_temperature) # 80 deg TEMPERATURE (There is no 75 deg C option)

            study.GetMaterial(      u"Magnet").SetValue(u"Poles", acm_template.d['EX']['DriveW_poles'])
            study.GetMaterial(      u"Magnet").SetDirectionXYZ(1, 1, 0)
            study.GetMaterial(      u"Magnet").SetAxisXYZ(0, 0, 1)
            study.GetMaterial(      u"Magnet").SetOriginXYZ(0, 0, 0)
            study.GetMaterial(      u"Magnet").SetPattern(u"RadialCircular")
            study.GetMaterial(      u"Magnet").SetValue(u"StartAngle", 0.5*360/(2*acm_template.SI['p']) ) # 半个极距

            # study.GetDesignTable().AddParameterVariableName(u"StatorPM: Direction")
            # study.GetDesignTable().AddParameterVariableName(u"StatorPM: Inward/Outward")
            # study.GetDesignTable().AddCase()
            # study.GetDesignTable().SetValue(1, 2, u"(1, 0, 0)")
            # study.GetDesignTable().SetValue(1, 3, u"Inward")

            #debug
            # study.SetMaterialByName(u"Magnet", u"Ambient/Air/Air")

        # add_carbon_fiber_material(app)
    def add_circuit(self, app, model, study, acm_variant, bool_3PhaseCurrentSource=True):
        # Circuit - Current Source
        app.ShowCircuitGrid(True)
        study.CreateCircuit()

        # 4 pole motor Qs=24 dpnv implemented by two layer winding (6 coils). In this case, drive winding has the same slot turns as bearing winding
        def circuit(Grouping,turns,Rs,ampD,ampB,freq,phase=0, CommutatingSequenceD=0, CommutatingSequenceB=0, x=10,y=10, bool_3PhaseCurrentSource=True):
            study.GetCircuit().CreateSubCircuit("Star Connection", "Star Connection %s"%(Grouping), x, y)
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
            elif 'CPPM' in acm_variant.template.name: 
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
                dcB = ampB/np.sqrt(2)

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 0*phase_shift_drive) # "freq" variable cannot be used here. So pay extra attension here when you create new case of a different freq.
                f2 = app.FunctionFactory().Constant(dcB)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I1).SetFunction(func)

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 1*phase_shift_drive)
                f2 = app.FunctionFactory().Constant(dcB)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I2).SetFunction(func)

                func = app.FunctionFactory().Composite()
                f1 = app.FunctionFactory().Sin(ampD, freq, 2*phase_shift_drive)
                f2 = app.FunctionFactory().Constant(dcB)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I3).SetFunction(func)
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
        # for k,v in EX.items():
        #     print(k,v)
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
                    # print('[JMAG.py] Concentrated winding!')
                    UVW    = wily.layer_Y_phases[index_leftlayer]
                    UpDown = wily.layer_Y_signs[index_leftlayer]
                else:
                    # print('Distributed winding.')
                    if wily.layer_Y_phases[index_leftlayer] != UVW:
                        print('[JMAG.py] [Warning] Potential bug in your winding layout detected.')
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
            for region_object, boo in zip(list_region_objects, token['list_regions_to_remove']):
                if boo:
                    self.doc.GetSelection().Clear()
                    self.doc.GetSelection().Add(region_object)
                    self.doc.GetSelection().Delete()
        if 'inner_or_outer_region_to_remove' in token.keys():
            for REGION_NAME, boo in zip(['Region', 'Region.2'], token['inner_or_outer_region_to_remove']):
                if boo:
                    self.doc.GetSelection().Clear()
                    self.doc.GetSelection().Add(self.sketch.GetItem(REGION_NAME))
                    self.doc.GetSelection().Delete()

        for idx, region_object in enumerate(list_region_objects):
            # Mirror
            if self.bMirror == True:
                if self.edge4Ref is None:
                    self.regionMirrorCopy(region_object, edge4Ref=None, symmetryType=2, bMerge=bMirrorMerge) # symmetryType=2 means x-axis as ref
                else:
                    self.regionMirrorCopy(region_object, edge4Ref=self.edge4ref, symmetryType=None, bMerge=bMirrorMerge) # symmetryType=2 means x-axis as ref

            # RotateCopy
            if self.iRotateCopy >= 2:
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
                logger.debug('The model already exists for individual with index=%d. Skip it.' % individual_index)
                return -1 # the model is already drawn

        elif individual_index+1 <= app.NumModels(): # 一般是从零起步
            logger = logging.getLogger(__name__)
            logger.debug('The model already exists for individual with index=%d. Skip it.' % individual_index)
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
            logger.info(msg)
            # print(msg)
        else:
            print('Submit to JMAG_Scheduler...')
            job = study.CreateJob()
            job.SetValue("Title", study.GetName())
            job.SetValue("Queued", True)
            job.Submit(False) # Fallse:CurrentCase, True:AllCases
            logger.info('Submit %s to queue (Tran2TSS).'%(acm_variant.individual_name))
            # wait and check
            # study.CheckForCaseResults()
        app.Save()

    # @staticmethod
    def mesh_study(self, acm_variant, app, model, study, output_dir):

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

            if self.bool_suppressShaft == False:
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

        if False:
            # Export Image
            app.View().ShowAllAirRegions()
            # app.View().ShowMeshGeometry() # 2nd btn
            app.View().ShowMesh() # 3rn btn
            app.View().Zoom(3)
            app.View().Pan(-acm_variant.template.d['GP']['mm_r_ro'].value, 0)
            app.ExportImageWithSize(output_dir + model.GetName() + '.png', 2000, 2000)
            app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.

    ''' BELOW is for JMAG Designer
    '''
    def draw_CPPM(self, acm_variant, bool_pyx=False, bool_draw_whole_model=False):
        if bool_draw_whole_model:
            self.bMirror = False
            self.iRotateCopy = 1
        # blue
        # color_rgb_A = np.array([113, 142, 164])/255
        # color_rgb_B = np.array([73, 109, 137])/255

        # yellow
        # color_rgb_A = np.array([255, 252, 170])/255
        # color_rgb_B = np.array([212, 208, 166])/255

        # gray
        color_rgb_A = np.array([236,236,236])/255
        color_rgb_B = np.array([226,226,226])/255

        # Rotor Core 1
        if 1:
            list_regions_1 = acm_variant.rotorCore.draw(self, bool_draw_whole_model=bool_draw_whole_model)
            self.bMirror = False
            self.iRotateCopy = acm_variant.rotorCore.p
            region1 = self.prepareSection(list_regions_1, color=color_rgb_A)

        # Shaft
        if 1:
            list_regions = acm_variant.shaft.draw(self, bool_draw_whole_model=bool_draw_whole_model)
            region0 = self.prepareSection(list_regions)

        # Rotor Magnet
        if 1:
            list_regions = acm_variant.rotorMagnet.draw(self, bool_draw_whole_model=bool_draw_whole_model)
            self.bMirror = False
            self.iRotateCopy = acm_variant.rotorCore.p
            region2 = self.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

        # Sleeve
        # list_regions = acm_variant.sleeve.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        # regionS = self.prepareSection(list_regions)

        # This is only for post-processing and it is for handle a un-fixable filling bug with PyX.
        # if bool_pyx:
        #     region1 = self.prepareSection(list_regions_1, color=color_rgb_A)

        # Stator Core
        if 1:
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

            self.calculate_excitation_current(acm_variant)

            # # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
            # EX = acm_variant.template.d['EX']
            # CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * EX['WindingFill'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
            # CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
            # CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
            #     # try:
            #     #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
            #     # except AttributeError:
            #     #     # print(EX['wily'])
            #     #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily']['number_parallel_branch']
            #     #     print("[inner_rotor_motor.py] Reproduce design using jsonpickle will encounter error here: 'dict' object has no attribute 'number_parallel_branch', implying that the object wily has become a dict after jsonpickle.")
            #     #     # quit() 

            # # Maybe there is a bug here... regarding the excitation for suspension winding...
            # variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            # variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
            # EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
            # EX['DriveW_CurrentAmp'] = acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            # EX['BeariW_CurrentAmp'] = acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp
            # print('[inner_rotor_motor.py] Excitations have been over-written by the constraint on Js! Total, DriveW, BeariW [A]:', 
            #                                                                                             EX['CurrentAmp_per_phase'],
            #                                                                                             EX['DriveW_CurrentAmp'],
            #                                                                                             EX['BeariW_CurrentAmp'])

            # # acm_variant.spec_geometry_dict['DriveW_CurrentAmp'] = acm_variant.DriveW_CurrentAmp

            # slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
            # print('[JMAG.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

            # # print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
            # # print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
            # # print('---acm_variant.DriveW_CurrentAmp =', acm_variant.DriveW_CurrentAmp)
            # # print('---acm_variant.BeariW_CurrentAmp =', acm_variant.BeariW_CurrentAmp)
            # # print('---TORQUE_CURRENT_RATIO:', acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'])
            # # print('---SUSPENSION_CURRENT_RATIO:', acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            # Import Model into Designer
            self.save(acm_variant.name, self.show(acm_variant, toString=True))

        return True

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

            self.calculate_excitation_current(acm_variant)

                # # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
                # EX = acm_variant.template.d['EX']
                # CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * EX['WindingFill'] * EX['Js']*1e-6 * np.sqrt(2) #/2.2*2.8
                # CurrentAmp_per_conductor = CurrentAmp_in_the_slot / EX['DriveW_zQ']
                # CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                #     # try:
                #     #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily'].number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                #     # except AttributeError:
                #     #     # print(EX['wily'])
                #     #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wily']['number_parallel_branch']
                #     #     print("[inner_rotor_motor.py] Reproduce design using jsonpickle will encounter error here: 'dict' object has no attribute 'number_parallel_branch', implying that the object wily has become a dict after jsonpickle.")
                #     #     # quit() 

                # # Maybe there is a bug here... regarding the excitation for suspension winding...
                # variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
                # variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_parallel_branch is 1 for suspension winding
                # EX['CurrentAmp_per_phase'] = CurrentAmp_per_phase
                # EX['DriveW_CurrentAmp'] = acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                # EX['BeariW_CurrentAmp'] = acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp
                # print('[inner_rotor_motor.py] Excitations have been over-written by the constraint on Js! Total, DriveW, BeariW [A]:', 
                #                                                                                             EX['CurrentAmp_per_phase'],
                #                                                                                             EX['DriveW_CurrentAmp'],
                #                                                                                             EX['BeariW_CurrentAmp'])

                # # acm_variant.spec_geometry_dict['DriveW_CurrentAmp'] = acm_variant.DriveW_CurrentAmp

                # slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
                # print('[JMAG.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')

                # # print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                # # print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                # # print('---acm_variant.DriveW_CurrentAmp =', acm_variant.DriveW_CurrentAmp)
                # # print('---acm_variant.BeariW_CurrentAmp =', acm_variant.BeariW_CurrentAmp)
                # # print('---TORQUE_CURRENT_RATIO:', acm_variant.template.fea_config_dict['TORQUE_CURRENT_RATIO'])
                # # print('---SUSPENSION_CURRENT_RATIO:', acm_variant.template.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            # Import Model into Designer
            self.save(acm_variant.name, self.show(acm_variant, toString=True))

        return True
    
    def draw_doublySalient(self, acm_variant, bool_draw_whole_model=True):
        if bool_draw_whole_model:
            self.bMirror = False
            self.iRotateCopy = 1

        # gray
        color_rgb_A = np.array([236,236,236])/255
        color_rgb_B = np.array([226,226,226])/255

        # Rotor Core
        list_regions_1 = acm_variant.rotorCore.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        region1 = self.prepareSection(list_regions_1, color=color_rgb_A)

        # Shaft
        list_regions = acm_variant.shaft.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        # Stator Magnet
        list_regions = acm_variant.statorMagnet.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        region2 = self.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

        # Sleeve
        # list_regions = acm_variant.sleeve.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        # regionS = self.prepareSection(list_regions)

        # Stator Core
        list_regions = acm_variant.stator_core.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        region3 = self.prepareSection(list_regions, color=color_rgb_A)

        # Stator Winding
        list_regions = acm_variant.coils.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        region4 = self.prepareSection(list_regions)

        self.calculate_excitation_current(acm_variant)

        # Import Model into Designer
        self.save(acm_variant.name, self.show(acm_variant, toString=True))
        return True
    def draw_FSPM(self, acm_variant, bool_draw_whole_model=True):
        if bool_draw_whole_model:
            self.bMirror = False
            self.iRotateCopy = 1

        # gray
        color_rgb_A = np.array([236,236,236])/255
        color_rgb_B = np.array([226,226,226])/255

        SI = acm_variant.template.SI

        # Rotor Core
        list_regions_1 = acm_variant.rotorCore.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        self.bMirror = False
        self.iRotateCopy = SI['pm']
        region1 = self.prepareSection(list_regions_1, color=color_rgb_A)

        # Shaft
        list_regions = acm_variant.shaft.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        region0 = self.prepareSection(list_regions)

        # Stator Magnet
        list_regions = acm_variant.statorMagnet.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        self.bMirror = False
        self.iRotateCopy = SI['pe']*2
        region2 = self.prepareSection(list_regions, bRotateMerge=False, color=color_rgb_B)

        # Sleeve
        # list_regions = acm_variant.sleeve.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        # regionS = self.prepareSection(list_regions)

        # Stator Core
        list_regions = acm_variant.stator_core.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        self.bMirror = False
        self.iRotateCopy = SI['Qs']
        region3 = self.prepareSection(list_regions, color=color_rgb_A)

        # Stator Winding
        list_regions = acm_variant.coils.draw(self, bool_draw_whole_model=bool_draw_whole_model)
        region4 = self.prepareSection(list_regions)


        self.calculate_excitation_current(acm_variant)

        # Import Model into Designer
        self.save(acm_variant.name, self.show(acm_variant, toString=True))
        return True

    @staticmethod
    def calculate_excitation_current(acm_variant):
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
        # print('[inner_rotor_motor.py] Excitations have been over-written by the constraint on Js! Total, DriveW, BeariW [A]:', 
                                                                                                    # EX['CurrentAmp_per_phase'],
                                                                                                    # EX['DriveW_CurrentAmp'],
                                                                                                    # EX['BeariW_CurrentAmp'])

        slot_current_utilizing_ratio = (EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp']) / EX['CurrentAmp_per_phase']
        # print('[JMAG.py]---Heads up! slot_current_utilizing_ratio is', slot_current_utilizing_ratio, '  (PS: =1 means it is combined winding)')
        pass

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

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 从硬盘提取数据
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    @staticmethod
    def add_plots(axeses, dm, title=None, label=None, zorder=None, time_list=None, sfv=None, torque=None, range_ss=None, alpha=0.7):

        info = '%s' % (title)
        torque_average = sum(torque[-range_ss:])/len(torque[-range_ss:])
        info += '\nAverage Torque: %g Nm' % (torque_average)
        # torque error = torque - avg. torque
        torque_error = np.array(torque) - torque_average
        ss_max_torque_error = max(torque_error[-range_ss:]), min(torque_error[-range_ss:])
        # we use  half of peak-to-peak value to compute error rather than use peak-to-peak value
        normalized_torque_ripple   = 1.0*(ss_max_torque_error[0] - ss_max_torque_error[1]) / torque_average
        info += '\nNormalized Torque Ripple: %g %%' % (normalized_torque_ripple*100)

        info += '\nAverage Force Mag: %g N'% (sfv.ss_avg_force_magnitude)
        # we use half of peak-to-peak value to compute error rather than use peak-to-peak value
        normalized_force_error_magnitude = sfv.normalized_force_error_magnitude
        info += '\nNormalized Force Error Mag: %g%%, (+)%g%% (-)%g%%' % (normalized_force_error_magnitude*100,
                                                                    sfv.ss_max_force_err_abs[0]/sfv.ss_avg_force_magnitude*100,
                                                                    sfv.ss_max_force_err_abs[1]/sfv.ss_avg_force_magnitude*100)
        # we use peak value to compute error rather than use peak-to-peak value
        # 跟Eric讨论过后，确定了悬浮力的角度误差不能用峰峰值的一半，而是要用最大值和最小值中绝对值更大的那一个。
        force_error_angle = sfv.force_error_angle
        info += '\nMaximum Force Error Angle: %g [deg], (+)%g deg (-)%g deg' % (force_error_angle,
                                                                    sfv.ss_max_force_err_ang[0],
                                                                    sfv.ss_max_force_err_ang[1])
        info += '\nExtra Info:'
        info += '\n\tAverage Force Vector: (%g, %g) N' % (sfv.ss_avg_force_vector[0], sfv.ss_avg_force_vector[1])
        info += '\n\tTorque Ripple (Peak-to-Peak): %g Nm'% ( max(torque[-range_ss:]) - min(torque[-range_ss:]))
        info += '\n\tForce Mag Ripple (Peak-to-Peak): %g N'% (sfv.ss_max_force_err_abs[0] - sfv.ss_max_force_err_abs[1])

        if axeses is not None:
            # plot for torque and force
            ax = axeses[0][0]; ax.plot(time_list, torque,                                           alpha=alpha, label=label, zorder=zorder)
            ax.plot(time_list, np.ones(len(time_list)) * torque_average, 'k-')
            ax = axeses[0][1]; ax.plot(time_list, sfv.force_abs,                                    alpha=alpha, label=label, zorder=zorder)
            ax.plot(time_list, sfv.force_x,                                      alpha=alpha, label=label, zorder=zorder)
            ax.plot(time_list, sfv.force_y,                                      alpha=alpha, label=label, zorder=zorder)
            ax.plot(time_list, np.ones(len(time_list)) * sfv.ss_avg_force_magnitude, 'k-')
            ax = axeses[1][0]; ax.plot(time_list, 100*sfv.force_err_abs/sfv.ss_avg_force_magnitude, label=label, alpha=alpha, zorder=zorder)
            ax = axeses[1][1]; 
            ax.plot(time_list, sfv.force_err_ang_old_way,                        label=label, alpha=alpha, zorder=zorder)
            ax.plot(time_list, sfv.force_err_ang_new_way,                        label=label, alpha=alpha, zorder=zorder)
            ax.plot(time_list, sfv.force_ang,                        label=label, alpha=alpha, zorder=zorder)
            ax.plot(time_list, np.ones(len(time_list)) * sfv.ss_avg_force_angle, 'k-')

        # plot for visialization of power factor 
        # dm.get_voltage_and_current(range_ss)
        # ax = axeses[2][0]; ax.plot(dm.mytime, dm.myvoltage, label=label, alpha=alpha, zorder=zorder)
        # ax = axeses[2][0]; ax.plot(dm.mytime, dm.mycurrent, label=label, alpha=alpha, zorder=zorder)

        return info, torque_average, normalized_torque_ripple, sfv.ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle

    @staticmethod
    def read_csv_results_4_general_purpose(study_name, path_prefix, fea_config_dict, femm_solver, acm_variant=None):

        machine_type = acm_variant.template.machine_type

        # Read TranFEAwi2TSS results

        # logging.getLogger(__name__).debug('Look into: ' + path_prefix)

        # Torque
        basic_info = []
        time_list = []
        TorCon_list = []
        with open(path_prefix + study_name + '_torque.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if count<=8:
                    try:
                        float(row[1])
                    except:
                        continue
                    else:
                        basic_info.append((row[0], float(row[1])))
                else:
                    try:
                        time_list.append(float(row[0]))
                        TorCon_list.append(float(row[1]))
                    except ValueError as e:
                        print('You need to manually delete the csv files as there might be cases results in them. (It shuold has only one case.)')
                        raise e

        # Force
        basic_info = []
        # time_list = []
        ForConX_list = []
        ForConY_list = []
        with open(path_prefix + study_name + '_force.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if count<=8:
                    try:
                        float(row[1])
                    except:
                        continue
                    else:
                        basic_info.append((row[0], float(row[1])))
                else:
                    # time_list.append(float(row[0]))
                    ForConX_list.append(float(row[1]))
                    ForConY_list.append(float(row[2]))
        ForConAbs_list = np.sqrt(np.array(ForConX_list)**2 + np.array(ForConY_list)**2 )

        # Current
        key_list = []
        Current_dict = dict()
        with open(path_prefix + study_name + '_circuit_current.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if count<=8:
                    if 'Time' in row[0]: # Time(s)
                        for key in row:
                            key_list.append(key)
                            Current_dict[key] = []
                    else:
                        continue
                else:
                    for ind, val in enumerate(row):
                        Current_dict[key_list[ind]].append(float(val))
        Current_dict['CircuitCoilDefault'] = Current_dict[key_list[-1]]
        # print(key_list)

        # FluxLinkage (2022)
        key_list = []
        FluxLinkage_dict = dict()
        with open(path_prefix + study_name + '_flux_of_fem_coil.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if count<=8:
                    if 'Time' in row[0]: # Time(s)
                        for key in row:
                            key_list.append(key)
                            FluxLinkage_dict[key] = []
                    else:
                        continue
                else:
                    for ind, val in enumerate(row):
                        FluxLinkage_dict[key_list[ind]].append(float(val))
        FluxLinkage_dict['CircuitCoilDefault'] = FluxLinkage_dict[key_list[-1]]

        # Displacement angle (2022)
        DisplacementAngle_list = []
        with open(path_prefix + study_name + '_total_rotational_displacement.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if count<=8:
                    pass
                else:
                    DisplacementAngle_list.append(float(row[1]))

        # Terminal Voltage 
        new_key_list = []
        if fea_config_dict['delete_results_after_calculation'] == False:
            # file name is by individual_name like ID32-2-4_EXPORT_CIRCUIT_VOLTAGE.csv rather than ID32-2-4Tran2TSS_circuit_current.csv
            fname = path_prefix + study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv"
            # print 'Terminal Voltage - look into:', fname
            with open(fname, 'r') as f:
                count = 0
                for row in utility.csv_row_reader(f):
                    count +=1
                    if count==1: # Time | Terminal1 | Terminal2 | ... | Termial6
                        if 'Time' in row[0]: # Time, s
                            for key in row:
                                new_key_list.append(key) # Yes, you have to use a new key list, because the ind below bgeins at 0.
                                Current_dict[key] = []
                        else:
                            raise Exception('Problem with csv file for terminal voltage.')
                    else:
                        for ind, val in enumerate(row):
                            Current_dict[new_key_list[ind]].append(float(val))
        key_list += new_key_list

        # Loss
        # Iron Loss
        with open(path_prefix + study_name + '_iron_loss_loss.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if 'IM' in machine_type:
                    if count>8:
                        rotor_iron_loss = float(row[2]) # Rotor Core
                        stator_iron_loss = float(row[3]) # Stator Core
                        print('[utility.py] Iron loss:', stator_iron_loss, rotor_iron_loss)
                        break
                elif 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type:
                    if count>7:
                        print('[JMAG.py] This should be 0:', float(row[0]))
                        rotor_iron_loss = float(row[1]) # Rotor Core
                        stator_iron_loss = float(row[4]) # Stator Core
                        print('[JMAG.py] Iron loss:', stator_iron_loss, rotor_iron_loss)
                        break
        with open(path_prefix + study_name + '_joule_loss_loss.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if 'IM' in machine_type:
                    if count>8:
                        rotor_eddycurrent_loss  = float(row[2]) # Rotor Core
                        stator_eddycurrent_loss = float(row[3]) # Stator Core
                        print('[utility.py] Eddy current loss:', stator_eddycurrent_loss, rotor_eddycurrent_loss)
                        break
                elif 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type:
                    if count>7:
                        rotor_eddycurrent_loss  = float(row[1]) # Rotor Core
                        stator_eddycurrent_loss = float(row[4]) # Stator Core
                        print('[utility.py] Eddy current loss:', stator_eddycurrent_loss, rotor_eddycurrent_loss)
                        break
        with open(path_prefix + study_name + '_hysteresis_loss_loss.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if 'IM' in machine_type:
                    if count>8:
                        rotor_hysteresis_loss = float(row[2]) # Rotor Core
                        stator_hysteresis_loss = float(row[3]) # Stator Core
                        print('[utility.py] Hysteresis loss:', stator_hysteresis_loss, rotor_hysteresis_loss)
                        break
                elif 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type:
                    if count>7:
                        rotor_hysteresis_loss  = float(row[1]) # Rotor Core
                        stator_hysteresis_loss = float(row[4]) # Stator Core
                        print('[utility.py] Hysteresis loss:', stator_hysteresis_loss, rotor_hysteresis_loss)
                        break

        # Joule Loss (Copper and Magnet)
        rotor_Joule_loss_list = []
        with open(path_prefix + study_name + '_joule_loss.csv', 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if 'IM' in machine_type:
                    if count == 8:
                        headers = row
                        for idx_coil, h in enumerate(headers): # on server there are 3 air regions... while on PC there are 2...
                            if 'Coil' in h:
                                break

                    if count>8:
                        if count==8+1:
                            if 'Coil' not in headers[idx_coil]:
                                print('[utility.py]', headers)
                                raise Exception('Error when load csv data for Coil.')
                            stator_copper_loss = float(row[idx_coil]) # Coil # it is the same over time, this value does not account for end coil

                        if 'Cage' not in headers[idx_coil-1]:
                            print('[utility.py]', headers)
                            raise Exception('Error when load csv data for Cage.')
                        rotor_Joule_loss_list.append(float(row[idx_coil-1])) # Cage

                elif 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type:
                    if count == 7: # 少一个slip变量，所以不是8，是7。
                        headers = row
                        for idx_coil, h in enumerate(headers):
                            if 'Coil' in h:
                                break

                    if count>7:
                        if count==7+1:
                            if 'Coil' not in headers[idx_coil]:
                                print('[utility.py]', headers)
                                raise Exception('Error when load csv data for Coil.')
                            stator_copper_loss = float(row[idx_coil]) # Coil # it is the same over time, this value does not account for end coil

                        if 'Magnet' not in headers[idx_coil-1]:
                            print('[utility.py]', headers)
                            raise Exception('Error when load csv data for Magnet.')
                        rotor_Joule_loss_list.append(float(row[idx_coil-1])) # Magnet

        # use the last 1/4 period data to compute average copper loss of Tran2TSS rather than use that of Freq study

        # effective_part = rotor_Joule_loss_list[-int(0.5*fea_config_dict['designer.number_of_steps_2ndTSS']):] # number_of_steps_2ndTSS = steps for half peirod
        effective_part = rotor_Joule_loss_list[-int(fea_config_dict['designer.number_of_steps_2ndTSS']):]
        if len(effective_part) == 0: # there is no rotor eddy current (e.g., FSPM motor)
            rotor_Joule_loss = 0.0
            print('[JMAG.py] csv results: effective_part is []')
            print('[JMAG.py] csv results: effective_part is []')
            print('[JMAG.py] csv results: effective_part is []')
            # print(rotor_Joule_loss_list)
            # print(effective_part)
        else:
            rotor_Joule_loss = sum(effective_part) / len(effective_part)
        if 'PMSM' in machine_type or 'CPPM' in machine_type:
            print('[utility.py] Magnet Joule loss:', rotor_Joule_loss)

        if femm_solver is not None:
            # blockPrint()
            try:
                # convert rotor current results (complex number) into its amplitude
                femm_solver.list_rotor_current_amp = [abs(el) for el in femm_solver.vals_results_rotor_current] # el is complex number
                # settings not necessarily be consistent with Pyrhonen09's design: , STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1., TEMPERATURE_OF_COIL=75

                # slot_area_utilizing_ratio = (acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp) / acm_variant.CurrentAmp_per_phase
                # if slot_area_utilizing_ratio < 1:
                #     print('Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio, 'which means you are simulating a separate winding? If not, contrats--you found a bug...')
                #     print('DW, BW, Total:', acm_variant.DriveW_CurrentAmp, acm_variant.BeariW_CurrentAmp, acm_variant.CurrentAmp_per_phase)

                _s, _r, _sAlongStack, _rAlongStack, _Js, _Jr = femm_solver.get_copper_loss_pyrhonen(acm_variant.slot_area_utilizing_ratio*femm_solver.stator_slot_area, 
                                                                                                                            femm_solver.rotor_slot_area, 
                                                                                                                            total_CurrentAmp=acm_variant.DriveW_CurrentAmp+acm_variant.BeariW_CurrentAmp)
                s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = femm_solver.get_copper_loss_Bolognani(acm_variant.slot_area_utilizing_ratio*femm_solver.stator_slot_area, 
                                                                                                                            femm_solver.rotor_slot_area, 
                                                                                                                            total_CurrentAmp=acm_variant.DriveW_CurrentAmp+acm_variant.BeariW_CurrentAmp)

                msg1 = 'Pyrhonen : %g, %g | %g, %g | %g, %g ' % (_s, _r, _sAlongStack, _rAlongStack, _Js, _Jr) 
                msg2 = 'Bolognani: %g, %g | %g, %g | %g, %g ' % (s, r, sAlongStack, rAlongStack, Js, Jr) 
                logger = logging.getLogger(__name__)
                logger.debug(msg1)
                logger.debug(msg2)
            except Exception as e:
                raise e
            # enablePrint()
        else:
            SI = acm_variant.template.SI
            GP = acm_variant.template.d['GP']
            EX = acm_variant.template.d['EX']
            wily = EX['wily']
            copper_loss_parameters = [GP['mm_d_sleeve'].value + GP['mm_d_mech_air_gap'].value,
                                GP['mm_w_st'].value,
                                wily.number_parallel_branch,
                                EX['DriveW_zQ'],
                                wily.coil_pitch_y,
                                acm_variant.template.SI['Qs'],
                                EX['mm_template_stack_length'],
                                EX['DriveW_CurrentAmp'] + EX['BeariW_CurrentAmp'], # total current amplitude
                                GP['mm_r_ro'].value,       # mm
                                GP['mm_r_so'].value*2*1e-3 # m, stator_yoke_diameter_Dsyi
                                ]
            # slot_area_utilizing_ratio = (acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp) / acm_variant.CurrentAmp_per_phase
            # if slot_area_utilizing_ratio < 1:
            #     print('Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio, 'which means you are simulating a separate winding? If not, contrats--you found a bug...')
            #     print('DW, BW, Total:', acm_variant.DriveW_CurrentAmp, acm_variant.BeariW_CurrentAmp, acm_variant.CurrentAmp_per_phase)
            s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = utility.get_copper_loss_Bolognani(
                EX['slot_current_utilizing_ratio']*acm_variant.coils.mm2_slot_area*1e-6, 
                copper_loss_parameters=copper_loss_parameters, 
                STATOR_SLOT_FILL_FACTOR=acm_variant.template.SI['WindingFill'],
                TEMPERATURE_OF_COIL=acm_variant.template.SI['Temperature'])
            # s, r, sAlongStack, rAlongStack, Js, Jr = 0, 0, 0, 0, 0, 0

        dm = data_manager()
        dm.basic_info     = basic_info
        dm.time_list      = time_list
        dm.TorCon_list    = TorCon_list
        dm.ForConX_list   = ForConX_list
        dm.ForConY_list   = ForConY_list
        dm.ForConAbs_list = ForConAbs_list
        dm.Current_dict   = Current_dict
        dm.key_list       = key_list
        dm.jmag_loss_list = [   stator_copper_loss, 
                                rotor_Joule_loss, 
                                stator_iron_loss+rotor_iron_loss, 
                                stator_eddycurrent_loss+rotor_eddycurrent_loss, 
                                stator_hysteresis_loss+rotor_hysteresis_loss ]
        dm.femm_loss_list = [s, r, sAlongStack, rAlongStack, Js, Jr ]
        dm.Vol_Cu         = Vol_Cu
        dm.FluxLinkage_dict = FluxLinkage_dict
        dm.DisplacementAngle_list = DisplacementAngle_list
        return dm

    # def build_str_results(self, axeses, acm_variant, project_name, tran_study_name, dir_csv_output_folder, fea_config_dict, femm_solver=None):
    def build_str_results(self, acm_variant, project_name, tran_study_name, dir_csv_output_folder, fea_config_dict, femm_solver=None):
        # originate from fobj

        if fea_config_dict['bool_post_processing'] == False:
            self.fig_main, self.axeses = plt.subplots(2, 2, sharex=True, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
            utility.pyplot_clear(self.axeses)
        else:
            self.fig_main, self.axeses = None, None

        machine_type = acm_variant.template.machine_type

        try:
            self.dm = dm = self.read_csv_results_4_general_purpose(tran_study_name, dir_csv_output_folder, fea_config_dict, femm_solver, acm_variant=acm_variant)

            # Get peak to peak coil flux linkage
            coil_flux_linkage_peak2peak_value_results = []
            for k, v in dm.FluxLinkage_dict.items():
                coil_flux_linkage_peak2peak_value_results.append( max(v) - min(v) )
                print(k, max(v), min(v))
            coil_flux_linkage_peak2peak_value = np.average(coil_flux_linkage_peak2peak_value_results)
            print(coil_flux_linkage_peak2peak_value_results)

        except Exception as e:
            print(e)
            logging.getLogger(__name__).error('Error when loading csv results for Tran2TSS. Check the Report of JMAG Designer. (Maybe Material is not added.)', exc_info=True)
            msg = 'CSV results are not found. Will re-build and re-run the JMAG project...' 
            raise utility.ExceptionReTry(msg)
            # return None
            # raise e

        if fea_config_dict['designer.number_cycles_in_3rdTSS'] == 0:
            number_of_steps_at_steady_state = fea_config_dict['designer.number_of_steps_2ndTSS']
        else:
            number_of_steps_3rdTSS = fea_config_dict['designer.number_cycles_in_3rdTSS']*self.fea_config_dict['designer.StepPerCycle_3rdTSS']
            number_of_steps_at_steady_state = number_of_steps_3rdTSS
        dm.number_of_steps_at_steady_state = number_of_steps_at_steady_state

        basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
        sfv = utility.suspension_force_vector(ForConX_list, ForConY_list, range_ss=number_of_steps_at_steady_state) # samples in the tail that are in steady state
        str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle = \
            self.add_plots( self.axeses, dm,
                        title=tran_study_name,
                        label='Transient FEA w/ 2 Time Step Sections',
                        zorder=8,
                        time_list=time_list,
                        sfv=sfv,
                        torque=TorCon_list,
                        range_ss=sfv.range_ss)
        str_results += '\n\tbasic info:' +   ''.join(  [str(el) for el in basic_info])

        if dm.jmag_loss_list is None:
            raise Exception('Loss data is not loaded?')
        else:
            str_results += '\n\tjmag loss info: '  + ', '.join(['%g'%(el) for el in dm.jmag_loss_list]) # dm.jmag_loss_list = [stator_copper_loss, rotor_copper_loss, stator_iron_loss, stator_eddycurrent_loss, stator_hysteresis_loss]

        str_results += '\n\tfemm loss info: '  + ', '.join(['%g'%(el) for el in dm.femm_loss_list])

        if fea_config_dict['delete_results_after_calculation'] == False:
            power_factor = dm.power_factor(number_of_steps_at_steady_state, targetFreq=acm_variant.template.d['EX']['DriveW_Freq'])
            str_results += '\n\tPF: %g' % (power_factor)

        # compute the fitness 
        rotor_volume = acm_variant.template.get_rotor_volume() 
        rotor_weight = acm_variant.template.get_rotor_weight()
        shaft_power  = acm_variant.template.d['EX']['Omega'] * torque_average # make sure update_mechanical_parameters is called so that Omega corresponds to slip_freq_breakdown_torque

        if 'IM' in acm_variant.template.machine_type:
            if False: # fea_config_dict['jmag_run_list'][0] == 0
                # by JMAG only
                copper_loss  = dm.jmag_loss_list[0] + dm.jmag_loss_list[1] 
                iron_loss    = dm.jmag_loss_list[2] 
            else:
                # by JMAG for iron loss and FEMM for copper loss
                if dm.femm_loss_list[0] is None: # this will happen for running release_design.py
                    copper_loss  = dm.jmag_loss_list[0] + dm.jmag_loss_list[1]
                else:
                    copper_loss  = dm.femm_loss_list[0] + dm.femm_loss_list[1]
                iron_loss = dm.jmag_loss_list[2] 
        elif 'PM' in acm_variant.template.machine_type:
            # Rotor magnet loss by JMAG
            magnet_Joule_loss = dm.jmag_loss_list[1]
            # Stator copper loss by Binder and Bolognani 2006
            copper_loss = dm.femm_loss_list[0] + magnet_Joule_loss
            iron_loss = dm.jmag_loss_list[2] 
        else:
            raise Exception('Unknown machine type:', acm_variant.template.machine_type)

        windage_loss = utility.get_windage_loss(acm_variant, acm_variant.template.d['EX']['mm_template_stack_length'])

        # 这样计算效率，输出转矩大的，铁耗大一倍也没关系了，总之就是气隙变得最小。。。要不就不要优化气隙了。。。
        total_loss   = copper_loss + iron_loss + windage_loss
        efficiency   = shaft_power / (total_loss + shaft_power)  # 效率计算：机械功率/(损耗+机械功率)
        str_results  += '\n\teta, windage, total_loss: %g, %g, %g' % (efficiency, windage_loss, total_loss)

        # for easy access to codes
        machine_results = [power_factor, efficiency, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
        machine_results.extend(dm.jmag_loss_list)
        if dm.femm_loss_list is None:
            raise
        machine_results.extend(dm.femm_loss_list)
        machine_results.extend([windage_loss, total_loss])

        str_machine_results = ','.join('%g'%(el) for el in machine_results if el is not None) # note that femm_loss_list can be None called by release_design.py

        cost_function_O1, list_cost_O1 = utility.compute_list_cost(utility.use_weights('O1'), rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss)
        cost_function_O2, list_cost_O2 = utility.compute_list_cost(utility.use_weights('O2'), rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss)

        ################################################################
        # NEW CODES for rated performance
        ################################################################
        # caculate the fitness
        print('[utility.py]', '-'*40)
        print('[utility.py] Calculate the fitness for', acm_variant.name)

        # LOSS
        if 'IM' in machine_type:
            stator_copper_loss_along_stack = dm.femm_loss_list[2]
            magnet_Joule_loss = 0.0
            rotor_copper_loss_along_stack  = dm.femm_loss_list[3]

            stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack
            rotor_copper_loss_in_end_turn  = dm.femm_loss_list[1] - rotor_copper_loss_along_stack

        elif 'PMSM' in machine_type:
            stator_copper_loss_along_stack = dm.femm_loss_list[2]
            magnet_Joule_loss
            rotor_copper_loss_along_stack = 0.0

            stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack
            rotor_copper_loss_in_end_turn  = 0

        elif 'FSPM' in machine_type:
            stator_copper_loss_along_stack = dm.femm_loss_list[2]
            magnet_Joule_loss
            rotor_copper_loss_along_stack = 0.0

            stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack
            rotor_copper_loss_in_end_turn  = 0

        elif 'CPPM' in machine_type:
            stator_copper_loss_along_stack = dm.femm_loss_list[2]
            magnet_Joule_loss
            rotor_copper_loss_along_stack = 0.0

            stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack
            rotor_copper_loss_in_end_turn  = 0

        speed_rpm       = acm_variant.template.SI['ExcitationFreqSimulated'] * 60 / acm_variant.template.SI['p'] # rpm
        required_torque = acm_variant.template.SI['mec_power'] / (2*np.pi*speed_rpm)*60

        rated_ratio                          = required_torque / torque_average 
        rated_stack_length_mm                = rated_ratio * acm_variant.template.d['EX']['mm_template_stack_length']
        rated_stator_copper_loss_along_stack = rated_ratio * stator_copper_loss_along_stack
        rated_magnet_Joule_loss              = rated_ratio * magnet_Joule_loss
        rated_rotor_copper_loss_along_stack  = rated_ratio * rotor_copper_loss_along_stack
        rated_iron_loss                      = rated_ratio * dm.jmag_loss_list[2]
        rated_windage_loss                   = utility.get_windage_loss(acm_variant, rated_stack_length_mm)

        # total_loss   = copper_loss + iron_loss + windage_loss
        rated_total_loss =  rated_stator_copper_loss_along_stack \
                        + rated_magnet_Joule_loss \
                        + rated_rotor_copper_loss_along_stack \
                        + stator_copper_loss_in_end_turn \
                        + rotor_copper_loss_in_end_turn \
                        + rated_iron_loss \
                        + rated_windage_loss


        # THERMAL
        if 'IM' in machine_type:
            stator_current_density = dm.femm_loss_list[4]
            rotor_current_density  = dm.femm_loss_list[5]
            # print('Current density [Arms/m^2]:', stator_current_density, rotor_current_density, sep='\n')
            # if rotor_current_density > 8e6:
            #     print('rotor_current_density is over 8e6 Arms/m^2')
        else:
                                    # 基波电流幅值（在一根导体里的电流，六相逆变器中的GroupBDW相的电流，所以相当于已经考虑了并联支路数了）
            stator_current_density = dm.ui_info[2] / 1.4142135623730951 / (acm_variant.coils.mm2_slot_area*1e-6/acm_variant.template.d['EX']['DriveW_zQ'])
            print('[utility.py] Data Magager: stator_current_density (GroupBDW) = %g Arms/m^2'%(stator_current_density))
            rotor_current_density = 0

        print('[utility.py] Required torque: %g Nm'%(required_torque))
        print("[utility.py] acm_variant.template.d['EX']['Omega']: %g rad/s"%(acm_variant.template.d['EX']['Omega']))
        rated_shaft_power  = acm_variant.template.d['EX']['Omega'] * required_torque
        rated_efficiency   = rated_shaft_power / (rated_total_loss + rated_shaft_power)  # 效率计算：机械功率/(损耗+机械功率)

        rated_rotor_volume = np.pi*(acm_variant.template.d['GP']['mm_r_ro'].value*1e-3)**2 * (rated_stack_length_mm*1e-3)
        print('[utility.py] rated_stack_length_mm =', rated_stack_length_mm)

        # This weighted list suggests that peak-to-peak torque ripple of 5% is comparable with Em of 5% or Ea of 1 deg. Ref: Ye gu ECCE 2018
        # Eric suggests Ea is 1 deg. But I think this may be too much emphasis on Ea so large Trip does not matter anymore (not verified yet).
        list_weighted_ripples = [normalized_torque_ripple/0.05, normalized_force_error_magnitude/0.05, force_error_angle]


        # Torque per Rotor Volume
        TRV = required_torque / rated_rotor_volume
        Cost = 0


        # - Cost # Note 1/ 1.6387e-5 = 61023.744 # Note 1 cubic inch is 1.6387e-5 cubic meter
        price_per_volume_steel    = 0.28  * 61023.744 # $/in^3 (M19 Gauge26) # 0.23 for low carbon, semi-processed 24 Gauge electrical steel
        price_per_volume_copper   = 1.2   * 61023.744 # $/in^3 wire or bar or end-ring
        price_per_volume_magnet   = 11.61 * 61023.744 # $/in^3 NdFeB PM
        # price_per_volume_aluminum = 0.88  / 16387.064 # $/in^3 wire or cast Al
        # Vol_Fe = (2*acm_variant.template.d['GP']['mm_r_so'].value*1e-3) ** 2 * (rated_stack_length_mm*1e-3) # 注意，硅钢片切掉的方形部分全部消耗了。# Option 1 (Jiahao)
        Vol_Fe = ( np.pi*(acm_variant.template.d['GP']['mm_r_so'].value*1e-3)**2 - np.pi*(acm_variant.template.d['GP']['mm_r_ri'].value*1e-3)**2 ) * (rated_stack_length_mm*1e-3) # Option 2 (Eric)
        if 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type:
            if 'FSPM' in machine_type:
                Vol_PM = (acm_variant.statorMagnet.mm2_magnet_area*1e-6) * (rated_stack_length_mm*1e-3)
                print('[utility.py] Area_PM', (acm_variant.statorMagnet.mm2_magnet_area*1e-6))
            else:
                Vol_PM = (acm_variant.rotorMagnet.mm2_magnet_area*1e-6) * (rated_stack_length_mm*1e-3)
                print('[utility.py] Area_PM', (acm_variant.rotorMagnet.mm2_magnet_area*1e-6))
        else:
            Vol_PM = 0.0
        # print('[utility.py] Area_Fe', (acm_variant.template.d['GP']['mm_r_so'].value*1e-3) ** 2)
        # print('[utility.py] Area_Cu (est.)', dm.Vol_Cu/(rated_stack_length_mm*1e-3))
        # print('[utility.py] Volume_Fe',    Vol_Fe)
        # print('[utility.py] Volume_Cu', dm.Vol_Cu)
        # print('[utility.py] Volume_PM',    Vol_PM)
        Cost =    Vol_Fe * price_per_volume_steel \
                + dm.Vol_Cu * price_per_volume_copper\
                +    Vol_PM * price_per_volume_magnet
        Cost_Fe =    Vol_Fe * price_per_volume_steel 
        Cost_Cu = dm.Vol_Cu * price_per_volume_copper
        Cost_PM =    Vol_PM * price_per_volume_magnet
        print(f'[utility.py] Cost_Fe: {Cost_Fe}')
        print(f'[utility.py] Cost_Cu: {Cost_Cu}')
        print(f'[utility.py] Cost_PM: {Cost_PM}')

        if acm_variant.template.fea_config_dict["moo.fitness_OA"] == 'TorqueDensity':
            f1 = -TRV
        elif acm_variant.template.fea_config_dict["moo.fitness_OA"] == 'Cost':
            f1 = Cost
        else:
            raise 

        if acm_variant.template.fea_config_dict["moo.fitness_OB"] == 'Efficiency':
            # - Efficiency @ Rated Power
            f2 = - rated_efficiency
        elif acm_variant.template.fea_config_dict["moo.fitness_OB"] == 'TorqueRipple':
            f2 = normalized_torque_ripple
        elif acm_variant.template.fea_config_dict["moo.fitness_OB"] == 'ForceErrorMagnitude':
            f2 = normalized_force_error_magnitude
        else:
            raise 

        if acm_variant.template.fea_config_dict["moo.fitness_OC"] == 'BearinglessRippleSum':
            # Ripple Performance (Weighted Sum)
            f3 = sum(list_weighted_ripples)
        elif acm_variant.template.fea_config_dict["moo.fitness_OC"] == 'ForceErrorMagnitude':
            f3 = normalized_force_error_magnitude
        elif acm_variant.template.fea_config_dict["moo.fitness_OC"] == 'ForceErrorAngle':
            f3 = force_error_angle
        else:
            raise 

        if acm_variant.template.fea_config_dict["moo.fitness_OD"] is not None:
            f4 = 0

        FRW = ss_avg_force_magnitude / rotor_weight
        print('[utility.py] FRW:', FRW, 'Rotor weight:', rotor_weight, 'Stack length:', acm_variant.template.d['EX']['mm_template_stack_length'], 'Rated stack length:', rated_stack_length_mm)
        rated_rotor_volume = acm_variant.template.get_rotor_volume(stack_length=rated_stack_length_mm) 
        rated_rotor_weight = acm_variant.template.get_rotor_weight(stack_length=rated_stack_length_mm)
        print('[utility.py] rated_rotor_volume:', rated_rotor_volume, 'rated_rotor_weight:', rated_rotor_weight)

        rated_results = [   rated_shaft_power, 
                            rated_efficiency,
                            rated_total_loss, 
                            rated_stator_copper_loss_along_stack, 
                            rated_magnet_Joule_loss,
                            rated_rotor_copper_loss_along_stack, 
                            stator_copper_loss_in_end_turn, 
                            rotor_copper_loss_in_end_turn, 
                            rated_iron_loss, 
                            rated_windage_loss,
                            rated_rotor_volume,
                            rated_stack_length_mm,  # new!
                            acm_variant.template.d['EX']['mm_template_stack_length']]           # new! 在计算FRW的时候，我们只知道原来的叠长下的力，所以需要知道原来的叠长是多少。

        # print(type(acm_variant.counter)==type(''))
        # print(type(acm_variant.counter)==type(''))
        # print(type(acm_variant.counter)==type(''))

        str_results = '\n-------\n%s-%s\n%d,%d,O1=%g,O2=%g,f1=%g,f2=%g,f3=%g\n%s\n%s\n%s\n' % (
                        project_name, acm_variant.get_individual_name(), 
                        -1 if type(acm_variant.counter)==type('') else int(acm_variant.counter//acm_variant.template.fea_config_dict["moo.popsize"]), # generation count
                        -1 if type(acm_variant.counter)==type('') else acm_variant.counter, # individual count
                        cost_function_O1, cost_function_O2, f1, f2, f3,
                        str_machine_results,
                        ','.join(['%g'%(el) for el in rated_results]), # 改为输出 rated_results
                        ','.join(['%g'%(el) for el in acm_variant.template.build_design_parameters_list()]) ) + str_results

        # str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss, cost_function = results_to_be_unpacked

        # write design evaluation data to file
        with open(dir_csv_output_folder[:-4] + 'swarm_data.txt', 'a') as f:
            f.write(str_results)

        # if fea_config_dict['use_weights'] == 'O1':
        #     cost_function = cost_function_O1
        # elif fea_config_dict['use_weights'] == 'O2':
        #     cost_function = cost_function_O2
        # else:
        #     raise Exception('Not implemented error.')

        return (cost_function_O1, cost_function_O2), f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle, \
                project_name, acm_variant.get_individual_name(), \
                -1 if type(acm_variant.counter)==type('') else int(acm_variant.counter//acm_variant.template.fea_config_dict["moo.popsize"]), \
                acm_variant.counter,\
                power_factor, \
                rated_ratio, \
                rated_stack_length_mm, \
                rated_total_loss, \
                rated_stator_copper_loss_along_stack, \
                rated_magnet_Joule_loss, \
                rated_rotor_copper_loss_along_stack, \
                stator_copper_loss_in_end_turn, \
                rotor_copper_loss_in_end_turn, \
                rated_iron_loss, \
                rated_windage_loss, \
                str_results, \
                acm_variant.coils.mm2_slot_area, \
                coil_flux_linkage_peak2peak_value, \
                TRV, Cost, Cost_Fe, Cost_Cu, Cost_PM, \
                ss_avg_force_magnitude, rotor_weight, torque_average
    # str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss, cost_function


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

