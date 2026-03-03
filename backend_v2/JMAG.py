try:
    import win32com.client
except ImportError:
    import sys
    from unittest.mock import MagicMock
    sys.modules["win32com"] = MagicMock()
    sys.modules["win32com.client"] = MagicMock()

import os, logging, utility, numpy, builtins

try:
    import pythoncom  # 用于 COM 初始化
except ImportError:
    import sys
    from unittest.mock import MagicMock
    sys.modules["pythoncom"] = MagicMock()

import matplotlib
matplotlib.use('Agg')
from pylab import np, plt, mpl; import math
# logger = logging.getLogger(__name__)
# logger.debug('The mpl backend is %s', mpl.rcParams['backend'])
# mpl.use('Agg') # ('pdf') #   # https://github.com/matplotlib/matplotlib/issues/21950
EPS=0.01 # mm
from time import time as clock_time

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '_default', 'main_test'))

class JMAG(object): #< ToolBase & DrawerBase & MakerExtrnudeBase & MakerRevolveBase
    # JMAG Encapsulation for the JMAG Designer of JSOL Corporation.    
    def __init__(self, fea_config_dict=None):
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
        self.bool_suppressShaft = False
        self.JMAG_version_number = 20

        # Default fea_config_dict
        self.fea_config_dict = fea_config_dict if fea_config_dict is not None else {}

        self.flag_material_already_loaded = False

        self.verbose_drawing = False

    def open(self, Steel_name: str, expected_project_file_path: str, pc_name: str, dir_parent: str):
        if self.app is None:
            # 在 Streamlit 等多线程环境中，需要显式初始化 COM
            # 使用 COINIT_APARTMENTTHREADED 模式（适合单线程单元模型）
            try:
                pythoncom.CoInitializeEx(pythoncom.COINIT_APARTMENTTHREADED)
            except Exception:
                # COM 已经初始化或初始化失败，继续执行（可能已在主线程初始化）
                # 捕获所有异常，因为 CoInitializeEx 可能在不同情况下抛出不同的异常
                pass

            # 打开JMAG Designer
            try:
                app = win32com.client.Dispatch('designer.Application.200')
                # app = win32com.client.Dispatch('designer.Application.171')
                # app = win32com.client.gencache.EnsureDispatch('designer.Application.171') # https://stackoverflow.com/questions/50127959/win32-dispatch-vs-win32-gencache-in-python-what-are-the-pros-and-cons
            except:
                # print('JMAG 20.0 is not found. Will use any other JMAG version avaiilable.')
                try:
                    app = win32com.client.Dispatch('designer.Application')
                    # app = win32com.client.Dispatch('designer.Application.181')
                    # app = win32com.client.gencache.EnsureDispatch('designer.Application.171')
                except:
                    raise Exception('COM call to JMAG Designer is not successful.')

            self.JMAG_version_string = app.VersionString(0)
            self.JMAG_version_number = float(app.VersionString(0)[:2])

            if self.fea_config_dict['designer.show'] == True:
                app.Show()
            else:
                app.Hide()
            # app.Quit()
            self.app = app # means that the JMAG Designer is turned ON now.

            # Check if Steel_name does not contain 'M-15', 'M-19', or 'Arnon'
            if not any(substring in Steel_name for substring in ['M-15', 'M-19', 'Arnon']):
                print('No custom steel is added to JMAG Designer.')

            def add_steel(app, dir_parent, Steel_name: str):

                def add_M1xSteel(app, dir_parent, Steel_name="M-19 Steel Gauge-29"):

                    if '19' in Steel_name:
                        try:
                            BH = np.loadtxt(dir_parent + '../BH/M-19-Steel-BH-Curve-afterJMAGsmooth.BH', unpack=True, usecols=(0,1)) # after JMAG smooth, it beomces HB rather than BH
                        except OSError:
                            BH = np.loadtxt(dir_parent + './BH/M-19-Steel-BH-Curve-afterJMAGsmooth.BH', unpack=True, usecols=(0,1)) # after JMAG smooth, it beomces HB rather than BH
                    elif '15' in Steel_name:
                        BH = np.loadtxt(dir_parent + '../BH/M-15-Steel-BH-Curve.txt', unpack=True, usecols=(0,1))


                    app.GetMaterialLibrary().CreateCustomMaterial(Steel_name, "Custom Materials")
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("Density", 7.85)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("MagneticSteelPermeabilityType", 2)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("CoerciveForce", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).GetTable("BhTable").SetName("Untitled")

                    refarray = BH.T.tolist()

                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).GetTable("BhTable").SetTable(refarray)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("DemagnetizationCoerciveForce", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("MagnetizationSaturated", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("MagnetizationSaturated2", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("YoungModulus", 210000)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("ShearModulus", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("YoungModulusX", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("YoungModulusY", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("YoungModulusZ", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("ShearModulusXY", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("ShearModulusYZ", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("ShearModulusZX", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G11", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G12", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G13", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G14", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G15", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G16", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G22", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G23", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G24", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G25", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G26", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G33", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G34", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G35", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G36", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G44", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G45", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G46", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G55", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G56", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("G66", 0)

                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("MagnetizationSaturated2", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("MagnetizationSaturatedMakerValue", 0)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("Loss_Type", 1)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("LossConstantKhX", 143.)
                    app.GetMaterialLibrary().GetUserMaterial(Steel_name).SetValue("LossConstantKeX", 0.530)

                def add_Arnon5(app, dir_parent):
                    app.GetMaterialLibrary().CreateCustomMaterial("Arnon5-final", "Custom Materials")
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("Density", 7.85)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagneticSteelPermeabilityType", 2)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("CoerciveForce", 0)
                    # app.GetMaterialLibrary().GetUserMaterial(u"Arnon5-final").GetTable("BhTable").SetName(u"SmoothZeroPointOne")

                    BH = np.loadtxt(dir_parent + 'Arnon5/Arnon5-final.txt', unpack=True, usecols=(0,1))
                    refarray = BH.T.tolist()

                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").GetTable("BhTable").SetTable(refarray)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("DemagnetizationCoerciveForce", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturated", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturated2", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ExtrapolationMethod", 1)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulus", 210000)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulus", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulusX", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulusY", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulusZ", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulusXY", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulusYZ", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulusZX", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G11", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G12", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G13", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G14", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G15", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G16", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G22", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G23", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G24", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G25", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G26", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G33", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G34", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G35", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G36", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G44", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G45", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G46", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G55", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G56", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G66", 0)

                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturated2", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturatedMakerValue", 0)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("Loss_Type", 1)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("LossConstantKhX", 186.6)
                    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("LossConstantKeX", 0.07324)

                def add_Arnon7(app, dir_parent):
                    pass



                if 'M-15' in Steel_name:
                    add_M1xSteel(app, dir_parent, Steel_name="M-15 Steel")
                elif 'M-19' in Steel_name:
                    add_M1xSteel(app, dir_parent)
                elif 'Arnon5' == Steel_name:
                    add_Arnon5(app, dir_parent)       
                else:
                    return 'not custom steel'

            # to avoid tons of the same material in JMAG's material library
            fname = dir_parent + 'BH/.jmag_state.txt'
            pc_name = self.fea_config_dict['pc_name']
            steel_material_entry = pc_name + '/' + Steel_name

            # Check if file exists; if not, act as if material is not yet loaded
            if os.path.exists(fname):
                with open(fname, 'r') as f:
                    for line in f:
                        if steel_material_entry in line:
                            self.flag_material_already_loaded = True
                            # print('[JMAG.py] steel material already there:', steel_material_entry)
                            break

            if not os.path.exists(fname):
                if self.flag_material_already_loaded == False:
                    add_steel(app, dir_parent, Steel_name=Steel_name)
                    self.flag_material_already_loaded = True
                    with open(fname, 'a') as f:
                        f.write(steel_material_entry + '\n')
        else:
            app = self.app

        logger = logging.getLogger(__name__)
        logger.info('expected_project_file_path: %s', expected_project_file_path)
        if os.path.exists(expected_project_file_path):
            try:
                os.remove(expected_project_file_path)
                logger.info('Deleted existing JMAG project file: %s', expected_project_file_path)
            except Exception as e:
                logger.warning('Could not delete existing JMAG project file (maybe locked): %s. Error: %s', expected_project_file_path, e)
                attempts = 2
                temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)
                while os.path.exists(temp_path):
                    attempts += 1
                    temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)
                expected_project_file_path = temp_path
                logger.info('Will use alternative project file: %s', expected_project_file_path)

        app.NewProject("Untitled")
        app.SaveAs(os.path.abspath(expected_project_file_path)) 
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
    def pre_process_PMSM(self, app, model, acm_variant):
        wp = acm_variant.user_input['winding']

        logger = logging.getLogger(__name__)
        print(f"DEBUG JMAG FILE: {__file__}")
        start_time = clock_time()
        # logger.info('pre-process PMSM... Starting time: %g s.'%start_time)
        
        def verbose_print(*args):
            if getattr(builtins, 'verbose', False):
                print(*args)

        # view = app.View()
        # view.ClearSelect()
        # sel = view.GetCurrentSelection()
        # sel.SelectPart(123)
        # sel.SetBlockUpdateView(False)

        def export_image(app, model, path2SwarmData, suffix='.png'):
            wp = acm_variant.user_input['winding']
            app.View().ShowAllAirRegions()
            # app.View().ShowMeshGeometry() # 2nd btn
            app.View().ShowMesh() # 3rn btn
            app.View().Zoom(3)
            app.View().Pan(-acm_variant.user_input['geometry']['r_rotor_outer'], 0)
            app.ExportImageWithSize(path2SwarmData + '/jmag_screenshots/'  + suffix, 2000, 2000)

        # 2022-04-12: this is a temporary fix for JMAG Designer 21.0
        # self.jd.GetProject().GetModel(acm_variant.user_input['evaluation']['project_name']).GetStudy(acm_variant.user_input['evaluation']['project_name']).GetMeshControl().GetCondition(u"MagnetMeshCtrl").SetValue(u"Size", 0.002)
        if self.fea_config_dict['designer.show']:
            app.View().Pan(-acm_variant.user_input['geometry']['r_rotor_outer'], 0)

        # Use winding and target attributes
        
        p = wp['pole_count'] // 2
        s = 1
        Q = wp['slot_count']
        project_name = acm_variant.user_input['evaluation']['project_name']
        path2SwarmData = acm_variant.user_input['evaluation']['project_loc']

        if len(acm_variant.geometry.parts) < 3: # Rotor, Stator, Magnet(s)
            #+ self.show(acm_variant,toString=False)
            export_image(app, model, path2SwarmData + '/', suffix='-BadNumberOfParts.png')
            raise Exception('Too few parts: ' + str(len(acm_variant.geometry.parts)))

        export_image(app, model, path2SwarmData, suffix=project_name+'.png')

        # pre-process : you can select part by coordinate!
        ''' Group '''
        def group(name, id_list):
            wp = acm_variant.user_input['winding']
            model.GetGroupList().CreateGroup(name)
            for the_id in id_list:
                model.GetGroupList().AddPartToGroup(name, the_id)
                # model.GetGroupList().AddPartToGroup(name, name) #<- this also works

        part_ID_list = model.GetPartIDs()
        if isinstance(part_ID_list, int) or isinstance(part_ID_list, str):
            part_ID_list = (part_ID_list,)

        pole_count = wp['pole_count']
        self.id_rotorCore = id_rotorCore = part_ID_list[0]
        partIDRange_Magnet = list(part_ID_list[1 : 1 + pole_count]) if len(part_ID_list) > 1 else []
        self.id_statorCore = id_statorCore = part_ID_list[1 + pole_count] if len(part_ID_list) > 1 + pole_count else None
        partIDRange_Coil = list(part_ID_list[2 + pole_count :]) if len(part_ID_list) > 2 + pole_count else []

        # debug
        verbose_print(f"id_rotorCore={id_rotorCore}")
        verbose_print(f"partIDRange_Magnet={partIDRange_Magnet}")
        verbose_print(f"id_statorCore={id_statorCore}")
        verbose_print(f"partIDRange_Coil={partIDRange_Coil}")
        # model.SuppressPart(id_sleeve, 1)

        group("Magnet", partIDRange_Magnet)
        verbose_print(f"Group 'Magnet' created with IDs: {partIDRange_Magnet}")
        group("Coils", partIDRange_Coil)
        verbose_print(f"Group 'Coils' created with IDs: {partIDRange_Coil}")

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

        # Create Set for layer_X_phases
        Angle_StatorSlotSpan = 360/Q

        def get_PCoil(acm_variant):
            wp = acm_variant.user_input['winding']
            # Use the inner_coords centroids of the Coil part
            for part in acm_variant.geometry.parts:
                if 'coil' in part['name'].lower():
                    if 'inner_coords' in part and part['inner_coords']:
                        # Pick the region with the most negative Y for JMAG compatibility
                        print(f"[coil inner_coords] {part['inner_coords']}")
                        return min(part['inner_coords'].values(), key=lambda p: p[1])
            # Rough fallback based on GP
            gp = acm_variant.user_input['geometry']
            return [gp['r_rotor_outer'] + gp['d_air_gap'] + gp['d_tooth']*0.9, -gp['w_tooth']*0.5*1.1]         
        PCoil = get_PCoil(acm_variant)
        print(f"PCoil={PCoil}, Angle_StatorSlotSpan={Angle_StatorSlotSpan}")
        
        # Start JMAG logic
        self.doc.GetSelection().Clear()
        
        
        R = math.sqrt(PCoil[0]**2 + PCoil[1]**2)
        THETA = math.atan2(PCoil[1], PCoil[0])
        X = R*math.cos(THETA)
        Y = R*math.sin(THETA)
        countXL = 0

        for UVW, UpDown in zip(wp['layer_X_phases'],wp['layer_X_signs']):
            countXL += 1 
            print(f"CoilLX: count={countXL}, Phase={UVW}, Sign={UpDown}, Location=({X:.4f}, {Y:.4f}), Angle={math.degrees(THETA):.4f} deg")
            add_part_to_set("CoilLX%s%s %d"%(UVW,UpDown,countXL), X, Y)

            # print(X, Y, THETA)
            THETA += Angle_StatorSlotSpan/180.*math.pi
            X = R*math.cos(THETA)
            Y = R*math.sin(THETA)

        # Create Set for layer_Y_phases
        if PCoil[1] >= 0:
            raise Exception(f'Pay attention to the coil setup. The codes are written assuming the first coil is in the 12th slot. In other words PCoil[1] should negative, but {PCoil[1]=}')
        else:
            THETA = math.atan2(-PCoil[1], PCoil[0]) - (2*np.pi) / wp['slot_count']
        X = R*math.cos(THETA)
        Y = R*math.sin(THETA)
        countYL = 0
        for UVW, UpDown in zip(wp['layer_Y_phases'],wp['layer_Y_signs']):
            countYL += 1 
            print(f"CoilLY: count={countYL}, Phase={UVW}, Sign={UpDown}, Location=({X:.4f}, {Y:.4f}), Angle={math.degrees(THETA):.4f} deg")
            add_part_to_set("CoilLY%s%s %d"%(UVW,UpDown,countYL), X, Y)

            THETA += Angle_StatorSlotSpan/180.*math.pi
            X = R*math.cos(THETA)
            Y = R*math.sin(THETA)

        # Create Set for Magnets
        gp = acm_variant.user_input['geometry']
        R = gp['r_rotor_outer'] - 0.5 * gp['d_magnet']
        alpha_slot_pitch = 360.0 / wp['slot_count']
        alpha_pole_pitch = 360.0 / wp['pole_count']

        list_xy_magnets = []
        # list_xy_airWithinRotorSlot = []
        for ind in range(int(p*2)):
            if s==1:
                THETA = alpha_pole_pitch*ind /180.*math.pi
                X = R*math.cos(THETA)
                Y = R*math.sin(THETA)
                natural_ind = ind + 1
                print(f"Magnet: ind={ind}, natural_ind={natural_ind}, Location=({X:.4f}, {Y:.4f}), Angle={math.degrees(THETA):.4f} deg")
                add_part_to_set("Magnet %d"%(natural_ind), X, Y)
                list_xy_magnets.append([X,Y])
            else:     # v---This negative sign means we walk CCW to assign sets.
                raise NotImplementedError

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
        part_list_set('Motion_Region', list_xy_magnets, list_part_id=[id_rotorCore, ])
        part_list_set('MagnetSet', list_xy_magnets)
        msg = 'Time spent on pre-process PMSM is %g s.'%(clock_time() - start_time)
        logger.info(msg)

        return True

    def add_magnetic_transient_study(self, app, model, path2FEACsv, study_name, acm_variant):
        wp = acm_variant.user_input['winding']
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
        
        fea_config = acm_variant.user_input['fea_config_dict']
        
        study.GetStudyProperties().SetValue("ConversionType", 0)
        study.GetStudyProperties().SetValue("NonlinearMaxIteration", self.fea_config_dict['designer.max_nonlinear_iteration'])
        study.GetStudyProperties().SetValue("ModelThickness", wp['l_stack']) # [mm] Stack Length

        # Material
        self.add_material(study, acm_variant)

        # Conditions - Motion
        study.CreateCondition("RotationMotion", "RotCon") # study.GetCondition(u"RotCon").SetXYZPoint(u"", 0, 0, 1) # megbox warning
        study.GetCondition("RotCon").SetValue("AngularVelocity", int(wp['rated_speed']))
        study.GetCondition("RotCon").ClearParts()
        study.GetCondition("RotCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)

        study.GetCondition("RotCon").SetValue(u"InitialRotationAngle", acm_variant.user_input['target']['initial_rotation_angle'])
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
        # if wp['bool_CustomizedCircuit'] == True:
        #     acm_variant.add_circuit_customized(app, model, study)
        # else:

        self.add_circuit(app, model, study, acm_variant, bool_3PhaseCurrentSource=wp['bool_3PhaseCurrentSource'])

        # True: no mesh or field results are needed
        study.GetStudyProperties().SetValue("OnlyTableResults", fea_config['designer.OnlyTableResults'])

        # Linear Solver
        if False:
            # sometime nonlinear iteration is reported to fail and recommend to increase the accerlation rate of ICCG solver
            study.GetStudyProperties().SetValue("IccgAccel", 1.2) 
            study.GetStudyProperties().SetValue("AutoAccel", 0)
        else:
            # this can be said to be super fast over ICCG solver.
            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

        if fea_config['designer.MultipleCPUs']:
            # This SMP(shared memory process) is effective only if there are tons of elements. e.g., over 100,000.
            # too many threads will in turn make them compete with each other and slow down the solve. 2 is good enough for eddy current solve. 6~8 is enough for transient solve.
            study.GetStudyProperties().SetValue("UseMultiCPU", True)
            study.GetStudyProperties().SetValue("MultiCPU", 2) 

        # two sections of different time step size
        if True:
            number_cycles_in_1stTSS = fea_config['designer.number_cycles_in_1stTSS']
            number_cycles_in_2ndTSS = fea_config['designer.number_cycles_in_2ndTSS']
            number_cycles_in_3rdTSS = fea_config['designer.number_cycles_in_3rdTSS']
            number_cycles_prolonged = fea_config['designer.number_cycles_prolonged']
            number_of_steps_1stTSS = fea_config['designer.number_of_steps_1stTSS']
            number_of_steps_2ndTSS = fea_config['designer.number_of_steps_2ndTSS']
            number_of_steps_3rdTSS = number_cycles_in_3rdTSS * fea_config['designer.StepPerCycle_3rdTSS']
            DM = app.GetDataManager()
            DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")

            FREQUENCY = abs(wp['rated_speed'] / 60.0 * (wp['pole_count'] / 2.0))
            print(f"DEBUG JMAG: FREQUENCY = {FREQUENCY}, number_cycles_in_1stTSS = {number_cycles_in_1stTSS}")

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
                refarray[3][1] =    number_cycles_prolonged*acm_variant.user_input['fea_config_dict']['designer.TranRef-StepPerCycle'] 
                refarray[3][2] =        50
                # refarray[4][0] = refarray[3][0] + 0.5/FREQUENCY # 最后来一个超密的半周期400步
                # refarray[4][1] =    400
                # refarray[4][2] =        50
            number_of_total_steps = 1 + number_of_steps_1stTSS + number_of_steps_2ndTSS + number_of_steps_3rdTSS + number_cycles_prolonged*acm_variant.user_input['fea_config_dict']['designer.TranRef-StepPerCycle'] # [Double Check] don't forget to modify here!
            # print('[inner_rotor_motor.py]: refarray:', refarray)
            DM.GetDataSet("SectionStepTable").SetTable(refarray)
            study.GetStep().SetValue("Step", number_of_total_steps)
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

        # add equations
        study.GetDesignTable().AddEquation("freq")
        study.GetDesignTable().GetEquation("freq").SetType(0)
        study.GetDesignTable().GetEquation("freq").SetExpression("%g"%(FREQUENCY))
        study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency in Hz")

        study.GetDesignTable().AddEquation("speed")
        study.GetDesignTable().GetEquation("speed").SetType(1)
        study.GetDesignTable().GetEquation("speed").SetExpression("freq * %f" % (60 / (wp['pole_count'] // 2) ))
        study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed in r/min")

        # speed, freq, slip
        study.GetCondition("RotCon").SetValue("AngularVelocity", 'speed')

        # if wp['bool_DPNVorSEPA'] == False:
        #     app.ShowCircuitGrid(True)
        #     study.GetCircuit().GetComponent("CS4").SetValue("Frequency", FREQUENCY)
        #     study.GetCircuit().GetComponent("CS2").SetValue("Frequency", FREQUENCY)

        # max_nonlinear_iteration = 50
        # study.GetStudyProperties().SetValue(u"NonlinearMaxIteration", max_nonlinear_iteration)
        # study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
        # study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)

        # # add other excitation frequencies other than 500 Hz as cases
        # for case_no, FreqS in enumerate([50.0, slip_freq_breakdown_torque]):
        #     slip = slip_freq_breakdown_torque / FreqS
        #     study.GetDesignTable().AddCase()
        #     study.GetDesignTable().SetValue(case_no+1, 0, FreqS)
        #     study.GetDesignTable().SetValue(case_no+1, 1, slip)

        # 你把Tran2TSS计算周期减半！
        # 也要在计算铁耗的时候选择1/4或1/2的数据！（建议1/4）
        # 然后，手动添加end step 和 start step，这样靠谱！2019-01-09：注意设置铁耗条件（iron loss condition）的Reference Start Step和End Step。
        # print("TODO: 铁耗如果选了四分之一周期，JMAG会自动把结果扩展到全周期，需要保证磁场是从零开始的，详见JMAG帮助文档说明。如果是二分之一周期，则不会有这个问题。")

        # Iron Loss Calculation Condition
        # Stator 
        if fea_config['designer.AddIronLossCondition']:
            cond = study.CreateCondition("Ironloss", "IronLossConStator")
            cond.SetValue("RevolutionSpeed", int(wp['rated_speed']))
            cond.SetValue(u"Poles", wp['pole_count'])
            cond.ClearParts()
            sel = cond.GetSelection()
            # sel.SelectPartByPosition(acm_variant.template.SI['GP']['mm_r_si'].value + EPS, 0 ,0) # 2022-02-04 这里发现代码有点歧义：注意，实际上acm_variant.template.SI['GP']已经被修改了，acm_variant.template.SI['GP'] = acm_variant.GP。 # btw, this works!
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
        if acm_variant.user_input['fea_config_dict']['designer.AddIronLossCondition']:
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
        # Convert path2FEACsv to a full absolute path string
        csv_path_str = os.path.abspath(path2FEACsv).replace('\\', '/')
        if not csv_path_str.endswith('/'):
            csv_path_str += '/'
        study.GetStudyProperties().SetValue("CsvOutputPath", csv_path_str) # it's folder rather than file!
        if acm_variant.user_input['fea_config_dict']["designer.AddIronLossCondition"]:
            if self.JMAG_version_number >= 21:
                study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;TerminalVoltage;JouleLoss;StoredEnergy;TotalDisplacementAngle;Inductance;FEMCoilInductance;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss") # new since 2022
            else:
                study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss") # old
        else:
            if self.JMAG_version_number >= 21:
                study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;TerminalVoltage;JouleLoss;StoredEnergy;TotalDisplacementAngle;Inductance") # no iron loss condition
            else:
                study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle")
        study.GetStudyProperties().SetValue("DeleteResultFiles", acm_variant.user_input['fea_config_dict']['delete_results_after_calculation'])

        # Terminal Voltage/Circuit Voltage: Check for outputing CSV results 
        if 'Flux_Alternator' in acm_variant.user_input['target']['machine_class']:
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
            else: # JMAG_version_number == 17
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
        wp = acm_variant.user_input['winding']

        def safe_set_material(study, part_name, material_name):
            wp = acm_variant.user_input['winding']
            study.SetMaterialByName(part_name, material_name)
            # Check if material was actually set (some versions don't throw exception)
            if study.GetMaterial(part_name).GetName() == "":
                raise Exception("Material not found")

        # Find rotorCore part name
        rotorCoreName = "rotorCore" # Default
        for part in acm_variant.geometry.parts:
            if 'rotor' in part['name'].lower() and 'core' in part['name'].lower():
                rotorCoreName = part['name']
                break

        
        # Steel
        mat_dict = acm_variant.user_input['material']
        
        safe_set_material(study, rotorCoreName, mat_dict['rotor_core_steel_name'])
        study.GetMaterial(rotorCoreName).SetValue("Laminated", 1)
        study.GetMaterial(rotorCoreName).SetValue("LaminationFactor", mat_dict['lamination_factor'])

        safe_set_material(study, "statorCore", mat_dict['stator_core_steel_name'])
        study.GetMaterial("statorCore").SetValue("Laminated", 1)
        study.GetMaterial("statorCore").SetValue("LaminationFactor", mat_dict['lamination_factor'])

        # Copper
        study.SetMaterialByName("Coils", "Copper")

        # Magnet
        target_dict = acm_variant.user_input['target']
        if 'PMSM' in target_dict['machine_class']:
            safe_set_material(study, u"Magnet", mat_dict['magnet_material_name'])
            study.GetMaterial(u"Magnet").SetValue(u"EddyCurrentCalculation", 1)
            study.GetMaterial(u"Magnet").SetValue(u"Temperature", mat_dict['magnet_temperature']) 

            study.GetMaterial(u"Magnet").SetValue(u"Poles", wp['pole_count'])
            study.GetMaterial(u"Magnet").SetDirectionXYZ(1, 0, 0)
            study.GetMaterial(u"Magnet").SetAxisXYZ(0, 0, -1)
            study.GetMaterial(u"Magnet").SetOriginXYZ(0, 0, 0)
            study.GetMaterial(u"Magnet").SetPattern(u"RadialCircular")
            study.GetMaterial(u"Magnet").SetOrientation(False) # False: 南极朝右，北极朝左，True: 南极朝左，北极朝右
            study.GetMaterial(u"Magnet").SetValue(u"StartAngle", mat_dict['magnet_start_angle']) 
            study.GetMaterial(u"Magnet").SetValue(u"UseAnisotropicMagnet", 0)

        # add_carbon_fiber_material(app)

    def add_circuit(self, app, model, study, acm_variant, bool_3PhaseCurrentSource=True):
        wp = acm_variant.user_input['winding']
        # Circuit - Current Source
        app.ShowCircuitGrid(True)
        study.CreateCircuit()

        # 4 pole motor Qs=24 dpnv implemented by two layer winding (6 coils). In this case, drive winding has the same slot turns as bearing winding
        def circuit(Grouping,turns,Rs,ampD,ampB,freq,phase=0, CommutatingSequenceD=0, CommutatingSequenceB=0, x=10,y=10, bool_3PhaseCurrentSource=True):
            wp = acm_variant.user_input['winding']
            print(f'Add {Grouping} circuit')
            if bool_3PhaseCurrentSource:
                study.GetCircuit().CreateComponent("3-Phase Current Source", Grouping)
            else:
                study.GetCircuit().CreateComponent("3-Phase Voltage Source", Grouping)
            study.GetCircuit().CreateComponent("3-Phase Coil", Grouping)
            study.GetCircuit().CreateComponent("3-Phase Ground", Grouping)

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
                print('Use 3-Phase Current Source')
                # quit()

                study.GetCircuit().CreateComponent("3PhaseCurrentSource", "CS%s"%(Grouping))
                study.GetCircuit().CreateInstance("CS%s"%(Grouping), x-4, y+1)
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("Amplitude", ampD+ampB)
                study.GetCircuit().GetComponent("CS%s"%(Grouping)).SetValue("Frequency", freq) # bypassed JMAG Equation string
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
                f1 = app.FunctionFactory().Sin(ampD, freq, 0*phase_shift_drive) # JMAG Equation "freq" variable cannot be used here (Known API bug). So pay extra attension here when you create new case of a different freq.
                target_dict = acm_variant.user_input['target']
                mc = target_dict['machine_class']
                if 'CPPM' in mc or 'CSPPM' in mc: 
                    dcB = ampB/math.sqrt(2)
                    f2 = app.FunctionFactory().Constant(dcB)
                else:
                    f2 = app.FunctionFactory().Sin(ampB, freq, 0*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I1).SetFunction(func)

                func = app.FunctionFactory().Composite()
                f1 = app.FunctionFactory().Sin(ampD, freq, 1*phase_shift_drive)
                if 'CPPM' in mc or 'CSPPM' in mc: 
                    dcB = -0.5*ampB/math.sqrt(2)
                    f2 = app.FunctionFactory().Constant(dcB)
                else:
                    f2 = app.FunctionFactory().Sin(ampB, freq, 1*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I2).SetFunction(func)

                func = app.FunctionFactory().Composite()
                f1 = app.FunctionFactory().Sin(ampD, freq, 2*phase_shift_drive)
                if 'CPPM' in mc or 'CSPPM' in mc: 
                    dcB = -0.5*ampB/math.sqrt(2)
                    f2 = app.FunctionFactory().Constant(dcB)
                else:
                    f2 = app.FunctionFactory().Sin(ampB, freq, 2*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I3).SetFunction(func)

            study.GetCircuit().CreateComponent("Ground", "Ground")
            study.GetCircuit().CreateInstance("Ground", x+2, y+1)
        # 这里电流幅值中的0.5因子源自DPNV导致的等于2的平行支路数。没有考虑到这一点，是否会对initial design的有效性产生影响？
        # 仔细看DPNV的接线，对于转矩逆变器，绕组的并联支路数为2，而对于悬浮逆变器，绕组的并联支路数为1。

        
        npb = wp['number_parallel_branch']
        nwl = wp['number_winding_layer']
        # if acm_variant.user_input['fea_config_dict']['DPNV_separate_winding_implementation'] == True or acm_variant.template.spec_input_dict['bool_DPNVorSEPA'] == False:
        # if wp['bool_DPNVorSEPA'] == False:
        #     # either a separate winding or a DPNV winding implemented as a separate winding
        #     ampD =  0.5 * (wp['drive_winding_current']/npb + wp['bearing_winding_current']) # 为了代码能被四极电机和二极电机通用，代入看看就知道啦。
        #     ampB = -0.5 * (wp['drive_winding_current']/npb - wp['bearing_winding_current']) # 关于符号，注意下面的DriveW对应的circuit调用时的ampB前还有个负号！
        #     if bool_3PhaseCurrentSource != True:
        #         raise Exception('Logic Error Detected.')
        # else:

        # case: DPNV as an actual two layer winding
        ampD = wp['drive_winding_current']/npb
        ampB = wp['bearing_winding_current']
        if bool_3PhaseCurrentSource != False:
            raise Exception('Logic Error Detected.')

        circuit('GroupAC',  wp['wires_per_slot']/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=wp['phase_resistance'],ampD= ampD,
                            ampB=-ampB, freq=wp['excitation_frequency_simulated'], phase=0,
                            CommutatingSequenceD=wp['CommutatingSequenceD'],
                            CommutatingSequenceB=wp['CommutatingSequenceB'])
                            
        circuit('GroupBD',  wp['wires_per_slot']/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=wp['phase_resistance'],ampD= ampD,
                                ampB=+ampB, freq=wp['excitation_frequency_simulated'], phase=0,
                                CommutatingSequenceD=wp['CommutatingSequenceD'],
                                CommutatingSequenceB=wp['CommutatingSequenceB'],x=25) # CS4 corresponds to uauc (conflict with following codes but it does not matter.)

        # Link FEM Coils to Coil Set     
        # if acm_variant.user_input['fea_config_dict']['DPNV_separate_winding_implementation'] == True or acm_variant.template.spec_input_dict['bool_DPNVorSEPA'] == False:
        # if wp['bool_DPNVorSEPA'] == False:
        #     def link_FEMCoils_2_CoilSet(SetPrefix, Grouping, l1, l2):
        #         wp = acm_variant.user_input['winding']
        #         # link between FEM Coil Condition and Circuit FEM Coil
        #         for UVW in ['U','V','W']:
        #             which_phase = "%s%s-Phase"%(Grouping,UVW)
        #             study.CreateCondition("FEMCoil", which_phase)
        #             condition = study.GetCondition(which_phase)
        #             condition.SetLink("CircuitCoil%s%s"%(Grouping,UVW))
        #             condition.GetSubCondition("untitled").SetName("Coil Set 1")
        #             condition.GetSubCondition("Coil Set 1").SetName("delete")
        #         count = 0
        #         dict_dir = {'+':1, '-':0, 'o':None}
        #         # select the part to assign the FEM Coil condition
        #         for UVW, UpDown in zip(l1,l2):
        #             count += 1 
        #             if dict_dir[UpDown] is None:
        #                 # print 'Skip', UVW, UpDown
        #                 continue
        #             which_phase = "%s%s-Phase"%(Grouping,UVW)
        #             condition = study.GetCondition(which_phase)
        #             condition.CreateSubCondition("FEMCoilData", "Coil Set %d"%(count))
        #             subcondition = condition.GetSubCondition("Coil Set %d"%(count))
        #             subcondition.ClearParts()
        #             subcondition.AddSet(model.GetSetList().GetSet("Coil%s%s%s %d"%(SetPrefix,UVW,UpDown,count)), 0)
        #             subcondition.SetValue("Direction2D", dict_dir[UpDown])
        #         # clean up
        #         for UVW in ['U','V','W']:
        #             which_phase = "%s%s-Phase"%(Grouping,UVW)
        #             condition = study.GetCondition(which_phase)
        #             condition.RemoveSubCondition("delete")
        #     link_FEMCoils_2_CoilSet('LX', 'GroupAC', 
        #                             wp['dict_coil_connection']['layer X phases'], 
        #                             wp['dict_coil_connection']['layer X signs']) 
        #     link_FEMCoils_2_CoilSet('LY', 'GroupBD', 
        #                             wp['dict_coil_connection']['layer Y phases'], 
        #                             wp['dict_coil_connection']['layer Y signs'])
        # else:
        if True:
            # 两个改变，一个是激励大小的改变（本来是200A 和 5A，现在是205A和195A），
            # 另一个绕组分组的改变，现在的A相是上层加下层为一相，以前是用俩单层绕组等效的。

            # Link FEM Coils to Coil Set as double layer short pitched winding
            # Create FEM Coil Condition
            # here we map circuit component `Coil2A' to FEM Coil Condition 'phaseAuauc
            # here we map circuit component `Coil4A' to FEM Coil Condition 'phaseAubud
            for suffix in ['GroupAC', 'GroupBD']: # 仍然需要考虑poles，是因为为Coil设置Set那里的代码还没有更新。这里的2(acm_variant.DriveW_poles)和4(acm_variant.BeariW_poles)等价于leftlayer和rightlayer。
                for UVW in ['U','V','W']:
                    study.CreateCondition("FEMCoil", 'phase'+UVW+suffix)
                    # link between FEM Coil Condition and Circuit FEM Coil
                    condition = study.GetCondition('phase'+UVW+suffix)
                    condition.SetLink("CircuitCoil%s%s"%(suffix,UVW))
                    condition.GetSubCondition("untitled").SetName("delete")
            countXL = 0 # countXL indicates which slot the current rightlayer is in.
            index = 0
            dict_dir = {'+':1, '-':0}
            coil_pitch = wp['coil_pitch_y'] #wp['dict_coil_connection'][0]
            # select the part (via `Set') to assign the FEM Coil condition
            for UVW, UpDown in zip(wp['layer_X_phases'], wp['layer_X_signs']):


                countXL += 1 
                if wp['grouping_AC'][index] == 1:
                    suffix = 'GroupAC'
                else:
                    suffix = 'GroupBD'
                condition = study.GetCondition('phase'+UVW+suffix)
                condition = study.GetCondition('phase'+UVW+suffix)

                # right layer
                # print (countXL, "Coil Set %d"%(countXL), end=' ')
                condition.CreateSubCondition("FEMCoilData", "Coil Set Layer X %d"%(countXL))
                subcondition = condition.GetSubCondition("Coil Set Layer X %d"%(countXL))
                subcondition.ClearParts()
                subcondition.AddSet(model.GetSetList().GetSet("CoilLX%s%s %d"%(UVW,UpDown,countXL)), 0) # poles=4 means right layer, rather than actual poles
                subcondition.SetValue("Direction2D", dict_dir[UpDown])

                # left layer
                Q = wp['slot_count']
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
                if wp['bool_distributed_or_concentrated'] == False:
                    # print('[JMAG.py] Concentrated winding!')
                    UVW    = wp['layer_Y_phases'][index_leftlayer]
                    UpDown = wp['layer_Y_signs'][index_leftlayer]
                else:
                    # print('Distributed winding.')
                    if wp['layer_Y_phases'][index_leftlayer] != UVW:
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

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 画图
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

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
            return [line] # startxy[0],startxy[1],endxy[0],endxy[1]
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
            return [arc] # centerxy[0], centerxy[1], startxy[0], startxy[1], endxy[0], endxy[1]
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
            # quit()
            if idx == 0:
                region_object = self.sketch.GetItem('Region') # This is how you get access to the region you create.
                # quit()
            else:
                region_object = self.sketch.GetItem('Region.%d'%(idx+1)) # This is how you get access to the region you create.
            list_region_objects.append(region_object)
            # print(list_region_objects)
            # quit()
        # remove region
        if 'list_regions_to_remove' in token.keys():
            for region_object, boo in zip(list_region_objects, token['list_regions_to_remove']):
                if boo:
                    self.doc.GetSelection().Clear()
                    self.doc.GetSelection().Add(region_object)
                    self.doc.GetSelection().Delete()
        # quit()
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
                self.regionCircularPattern360Origin(idx, region_object, self.iRotateCopy, bMerge=bRotateMerge)
            # quit()
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

    def regionCircularPattern360Origin(self, idx, region, Q_float, bMerge=True):
        # index is used to define name of region

        Q_float = float(Q_float) # don't ask me, ask JSOL

        circular_pattern = self.sketch.CreateRegionCircularPattern()
        circular_pattern.SetProperty('Merge', bMerge)

        ref2 = self.doc.CreateReferenceFromItem(region)
        # ref2 = 
        circular_pattern.SetPropertyByReference('Region', ref2)
        # circular_pattern.SetPropertyByReference('Region.2', ref2)
        # if idx == 0:
        face_region_string = circular_pattern.GetProperty('Region')

        # else:
            # face_region_string = circular_pattern.GetProperty('Region.%d'%(index+1))
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
            self.SI = d
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
        model.SetDescription(im_variant.model_name_prefix + '\n' + im_variant.show(toString=False))

        if doNotRotateCopy:
            im_variant.pre_process_structural(app, d.listKeyPoints)
        else:
            im_variant.pre_process(app)

        model.CloseCadLink() # this is essential if you want to create a series of models
        return True

    @staticmethod
    def run_study(acm_variant, app, study, fea_config_dict, toc):
        logger = logging.getLogger(__name__)
        print(f"DEBUG JMAG: run_study called for study: {study.GetName()}")
        if fea_config_dict['designer.JMAG_Scheduler'] == False:
            logger.info('Run jam.exe...')
            print("DEBUG JMAG: Running jam.exe (study.RunAllCases)...")
            # if run_list[1] == True:
            try:
                study.RunAllCases()
            except Exception as error:
                raise error
            msg = 'Time spent on %s is %g s.'%(study.GetName() , clock_time() - toc)
            logger.info(msg)
            # print(msg)
        else:
            print('Submit to JMAG_Scheduler...')
            job = study.CreateJob()
            job.SetValue("Title", study.GetName())
            job.SetValue("Queued", True)
            job.Submit(False) # Fallse:CurrentCase, True:AllCases
            logger.info('Submit %s to queue (Tran2TSS).'%(acm_variant.user_input['evaluation']['project_name']))
            # wait and check
            # study.CheckForCaseResults()
        app.Save()

    def mesh_study(self, acm_variant, app, model, study, output_dir):

        # this `if' judgment is effective only if JMAG-DeleteResultFiles is False 
        # if not study.AnyCaseHasResult(): 

        fea_config = acm_variant.user_input['fea_config_dict']
        meshctrl = study.GetMeshControl()

        # Air Mesh Control
        # this is for multi slide planes, which we will not be using
        refarray = [[0 for i in range(2)] for j in range(1)]
        refarray[0][0] = 3
        refarray[0][1] = 1
        meshctrl.GetTable("SlideTable2D").SetTable(refarray) 

        meshctrl.SetValue("MeshType", 1) # make sure this has been exe'd: study.GetCondition(u"RotCon").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)
        meshctrl.SetValue("RadialDivision", 8) # for air region near which motion occurs
        meshctrl.SetValue("CircumferentialDivision", fea_config['designer.CircumferentialDivision']) #1440) # for air region near which motion occurs 这个数足够大，sliding mesh才准确。
        meshctrl.SetValue("AirRegionScale", 1.05) # [Model Length]: Specify a value within the following area. (1.05 <= value < 1000)

        meshctrl.SetValue("AirMeshSize", fea_config['designer.meshSizeAir']) # mm
        meshctrl.SetValue("AutoAirMeshSize", 0)
        meshctrl.SetValue("Adaptive", 0)

        meshctrl.SetValue("MeshSize", fea_config['designer.meshSize_General']) # mm

        # This is not neccessary for whole model FEA. In fact, for BPMSM simulation, it causes mesh error "The copy target region is not found".
        # meshctrl.CreateCondition("RotationPeriodicMeshAutomatic", "autoRotMesh") # with this you can choose to set CircumferentialDivision automatically

        # Add mesh for Magnet Set
        meshctrl.CreateCondition("Part", "MagnetMeshCtrl")
        meshctrl.GetCondition("MagnetMeshCtrl").SetValue("Size", fea_config['designer.meshSize_Magnet'])
        meshctrl.GetCondition("MagnetMeshCtrl").ClearParts()
        meshctrl.GetCondition("MagnetMeshCtrl").AddSet(model.GetSetList().GetSet("MagnetSet"), 0)

        # Add mesh for Stator Part
        meshctrl.CreateCondition(u"Part", u"StatorMeshCtrl")
        meshctrl.GetCondition(u"StatorMeshCtrl").SetValue(u"Size", fea_config['designer.meshSize_Stator'])
        meshctrl.GetCondition(u"StatorMeshCtrl").ClearParts()
        sel = meshctrl.GetCondition(u"StatorMeshCtrl").GetSelection()
        sel.SelectPart(self.id_statorCore)
        meshctrl.GetCondition(u"StatorMeshCtrl").AddSelected(sel)

        # Add mesh for Rotor Part
        meshctrl.CreateCondition(u"Part", u"RotorMeshCtrl")
        meshctrl.GetCondition(u"RotorMeshCtrl").SetValue(u"Size", fea_config['designer.meshSize_Rotor'])
        meshctrl.GetCondition(u"RotorMeshCtrl").ClearParts()
        sel = meshctrl.GetCondition(u"RotorMeshCtrl").GetSelection()
        sel.SelectPart(self.id_rotorCore)
        meshctrl.GetCondition(u"RotorMeshCtrl").AddSelected(sel)

        # if self.bool_suppressShaft == False:
        #     meshctrl.CreateCondition("Part", "ShaftMeshCtrl")
        #     meshctrl.GetCondition("ShaftMeshCtrl").SetValue("Size", meshSize_Shaft) 
        #     meshctrl.GetCondition("ShaftMeshCtrl").ClearParts()
        #     meshctrl.GetCondition("ShaftMeshCtrl").AddSet(model.GetSetList().GetSet("ShaftSet"), 0)

        def mesh_all_cases(study):
            numCase = study.GetDesignTable().NumCases()
            print(f"DEBUG JMAG: mesh_all_cases for {numCase} cases")
            for case in range(0, numCase):
                study.SetCurrentCase(case)
                if study.HasMesh() == False:
                    print(f"DEBUG JMAG: Creating mesh for case {case}...")
                    try:
                        study.CreateMesh()
                        print(f"DEBUG JMAG: Mesh created for case {case}")
                    except Exception as e:
                        print(f"DEBUG JMAG: Mesh creation failed for case {case}: {e}")
                        raise e
                # if case == 0:
                #     app.View().ShowAllAirRegions()
                #     app.View().ShowMeshGeometry()
                #     app.View().ShowMesh()
        mesh_all_cases(study)

        # Export Image of the Meshed Regions
        if False:
            app.View().ShowAllAirRegions()
            # app.View().ShowMeshGeometry() # 2nd btn
            app.View().ShowMesh() # 3rn btn
            app.View().Zoom(3)
            app.View().Pan(-acm_variant.geometry.r_rotor_outer.value, 0)
            app.ExportImageWithSize(output_dir + model.GetName() + '.png', 2000, 2000)
            app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.

    ''' BELOW is for JMAG Designer
    '''
    def draw_spmsm(self, acm_variant, bool_pyx=False):
        wp = acm_variant.user_input['winding']

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
        list_regions = acm_variant.statorCore.draw(self)
        self.bMirror = True
        self.iRotateCopy = acm_variant.statorCore.Q
        region3 = self.prepareSection(list_regions, color=color_rgb_A)

        if not bool_pyx:
            # Stator Winding
            list_regions = acm_variant.coils.draw(self)
            self.bMirror = False
            self.iRotateCopy = acm_variant.coils.statorCore.Q
            region4 = self.prepareSection(list_regions)

            self.calculate_excitation_current(acm_variant)

                # # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
                # EX = acm_variant.template.SIEX
                # CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * acm_variant.winding.fill_factor * EX['Js']*1e-6 * math.sqrt(2) #/2.2*2.8
                # CurrentAmp_per_conductor = CurrentAmp_in_the_slot / acm_variant.winding.wires_per_slot
                # CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wp'].number_of_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                #     # try:
                #     #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wp'].number_of_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                #     # except AttributeError:
                #     #     # print(EX['wp'])
                #     #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wp']['number_of_parallel_branch']
                #     #     print("[inner_rotor_motor.py] Reproduce design using jsonpickle will encounter error here: 'dict' object has no attribute 'number_of_parallel_branch', implying that the object wp has become a dict after jsonpickle.")
                #     #     # quit() 

                # # Maybe there is a bug here... regarding the excitation for suspension winding...
                # variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
                # variant_BeariW_CurrentAmp =  CurrentAmp_per_conductor * 1 # number_of_parallel_branch is 1 for suspension winding
                # wp['pole_count'] // 2hase_current_amplitude = CurrentAmp_per_phase
                # acm_variant.winding.drive_winding_current = acm_variant.user_input['fea_config_dict']['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                # acm_variant.winding.bearing_winding_current = acm_variant.user_input['fea_config_dict']['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp
                # print('[inner_rotor_motor.py] Excitations have been over-written by the constraint on Js! Total, DriveW, BeariW [A]:', 
                #                                                                                             wp['pole_count'] // 2hase_current_amplitude,
                #                                                                                             acm_variant.winding.drive_winding_current,
                #                                                                                             acm_variant.winding.bearing_winding_current)

                # # acm_variant.spec_geometry_dict['DriveW_CurrentAmp'] = acm_variant.DriveW_CurrentAmp

                # slot_current_utilizing_ratio_for_torque = (acm_variant.winding.drive_winding_current + acm_variant.winding.bearing_winding_current) / wp['pole_count'] // 2hase_current_amplitude
                # print('[JMAG.py]---Heads up! slot_current_utilizing_ratio_for_torque is', slot_current_utilizing_ratio_for_torque, '  (PS: =1 means it is combined winding)')

                # # print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                # # print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                # # print('---acm_variant.DriveW_CurrentAmp =', acm_variant.DriveW_CurrentAmp)
                # # print('---acm_variant.BeariW_CurrentAmp =', acm_variant.BeariW_CurrentAmp)
                # # print('---TORQUE_CURRENT_RATIO:', acm_variant.user_input['fea_config_dict']['TORQUE_CURRENT_RATIO'])
                # # print('---SUSPENSION_CURRENT_RATIO:', acm_variant.user_input['fea_config_dict']['SUSPENSION_CURRENT_RATIO'])

            # Import Model into Designer
            self.save(acm_variant.user_input['target']['machine_class'], self.show(acm_variant, toString=False))

        # import builtins
        # builtins.ad.visualize_dict['GeometricComponentsObjects']['rotorCore'] = acm_variant.rotorCore
        # builtins.ad.visualize_dict['GeometricComponentsObjects']['shaft'] = acm_variant.shaft
        # builtins.ad.visualize_dict['GeometricComponentsObjects']['rotorMagnet'] = acm_variant.rotorMagnet
        # builtins.ad.visualize_dict['GeometricComponentsObjects']['sleeve'] = acm_variant.sleeve
        # builtins.ad.visualize_dict['GeometricComponentsObjects']['statorCore'] = acm_variant.statorCore
        # builtins.ad.visualize_dict['GeometricComponentsObjects']['coils'] = acm_variant.coils

        return True

    # @staticmethod
    # def calculate_excitation_current(acm_variant):
    #     # 根据绕组的形状去计算可以放铜导线的面积，然后根据电流密度计算定子电流
    #     EX = acm_variant.template.SIEX
    #     CurrentAmp_in_the_slot = acm_variant.coils.mm2_slot_area * acm_variant.winding.fill_factor * EX['Js']*1e-6 * math.sqrt(2) #/2.2*2.8
    #     CurrentAmp_per_conductor = CurrentAmp_in_the_slot / acm_variant.winding.wires_per_slot
    #     CurrentAmp_per_phase = CurrentAmp_per_conductor * EX['wp'].number_of_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。

    #     # Maybe there is a bug here... regarding the excitation for suspension winding...
    #     # variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
    #     # variant_BeariW_CurrentAmp =  CurrentAmp_per_phase * 1 # number_of_parallel_branch is 1 for suspension winding
    #     wp['pole_count'] // 2hase_current_amplitude = CurrentAmp_per_phase
    #     variant_DriveW_CurrentAmp = acm_variant.winding.drive_winding_current = acm_variant.user_input['fea_config_dict']['circuit.TORQUE_CURRENT_RATIO'] * CurrentAmp_per_phase
    #     variant_BeariW_CurrentAmp = acm_variant.winding.bearing_winding_current = acm_variant.user_input['fea_config_dict']['circuit.SUSPENSION_CURRENT_RATIO'] * CurrentAmp_per_phase
    #     # print('[inner_rotor_motor.py] Excitations have been over-written by the constraint on Js! Total, DriveW, BeariW [A]:', 
    #                                                                                                 # wp['pole_count'] // 2hase_current_amplitude,
    #                                                                                                 # acm_variant.winding.drive_winding_current,
    #                                                                                                 # acm_variant.winding.bearing_winding_current)

    #     slot_current_utilizing_ratio_for_torque = (acm_variant.winding.drive_winding_current + acm_variant.winding.bearing_winding_current) / wp['pole_count'] // 2hase_current_amplitude
    #     print('[JMAG.py]---Heads up! slot_current_utilizing_ratio_for_torque is', slot_current_utilizing_ratio_for_torque, '  (PS: =1 means it is combined winding)')
    #     print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
    #     print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
    #     print('---acm_variant.DriveW_CurrentAmp =', variant_DriveW_CurrentAmp)
    #     print('---acm_variant.BeariW_CurrentAmp =', variant_BeariW_CurrentAmp)
    #     print('---TORQUE_CURRENT_RATIO:', acm_variant.user_input['fea_config_dict']['circuit.TORQUE_CURRENT_RATIO'])
    #     print('---SUSPENSION_CURRENT_RATIO:', acm_variant.user_input['fea_config_dict']['circuit.SUSPENSION_CURRENT_RATIO'])
    #     print('---zQ:', acm_variant.winding.wires_per_slot)


    ''' JMAG Description
    '''
    def show(self, acm_variant, toString=False):
        pass
        # def get_tuple_list(object):
        #     attrs = list(vars(object).items())
        #     key_list = [el[0] for el in attrs]
        #     val_list = [el[1] for el in attrs]
        #     the_dict = dict(list(zip(key_list, val_list)))
        #     sorted_key = sorted(key_list, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item)) # this is also useful for string beginning with digiterations '15 Steel'.
        #     tuple_list = [(key, the_dict[key]) for key in sorted_key]
        #     return tuple_list
        # variant_tuple_list = get_tuple_list(acm_variant)
        # # template_tuple_list = get_tuple_list(acm_variant.template)
        # if toString==False:
        #     logger = logging.getLogger(__name__)
        #     logger.info('- Bearingless PMSM Individual #%s', acm_variant.user_input['target']['machine_class'])
        #     logger.info('\t%s', ', \n\t'.join("%s = %s" % item for item in variant_tuple_list))
        #     # print(', \n\t'.join("%s = %s" % item for item in template_tuple_list))
        #     return ''
        # else:
        #     return '\n- Bearingless PMSM Individual #%s\n\t' % (acm_variant.user_input['target']['machine_class']) \
        #             + ', \n\t'.join("%s = %s" % item for item in variant_tuple_list) \
        #             + '\n' + '--'*10 + '\n- Template:\n\t' \
        #             + ', \n\t'.join("%s = %s" % item for item in template_tuple_list)

    @staticmethod
    def add_plots(axeses, dm, title=None, label=None, zorder=None, time_list=None, sfv=None, torque=None, range_ss=None, alpha=0.7):

        info = '%s' % (title)
        torque_average = sum(torque[-range_ss:])/len(torque[-range_ss:])
        if torque_average < 0.0:
            logger = logging.getLogger(__name__)
            logger.warning('Average Torque: %g Nm which is negative, it is converted to positive to continue.' % (torque_average))
            torque_average = abs(torque_average)

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
        # print(normalized_force_error_magnitude)
        # print(sfv.ss_max_force_err_abs[1])
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

        # if axeses is not None:
        #     # plot for torque and force
        #     ax = axeses[0][0]; ax.plot(time_list, torque,                                           alpha=alpha, label=label, zorder=zorder)
        #     ax.plot(time_list, np.ones(len(time_list)) * torque_average, 'k-')
        #     ax = axeses[0][1]; ax.plot(time_list, sfv.force_abs,                                    alpha=alpha, label=label, zorder=zorder)
        #     ax.plot(time_list, sfv.force_x,                                      alpha=alpha, label=label, zorder=zorder)
        #     ax.plot(time_list, sfv.force_y,                                      alpha=alpha, label=label, zorder=zorder)
        #     ax.plot(time_list, np.ones(len(time_list)) * sfv.ss_avg_force_magnitude, 'k-')
        #     ax = axeses[1][0]; ax.plot(time_list, 100*sfv.force_err_abs/sfv.ss_avg_force_magnitude, label=label, alpha=alpha, zorder=zorder)
        #     ax = axeses[1][1]; 
        #     ax.plot(time_list, sfv.force_err_ang_old_way,                        label=label, alpha=alpha, zorder=zorder)
        #     ax.plot(time_list, sfv.force_err_ang_new_way,                        label=label, alpha=alpha, zorder=zorder)
        #     ax.plot(time_list, sfv.force_ang,                        label=label, alpha=alpha, zorder=zorder)
        #     ax.plot(time_list, np.ones(len(time_list)) * sfv.ss_avg_force_angle, 'k-')

        # plot for visialization of power factor 
        # dm.get_voltage_and_current(range_ss)
        # ax = axeses[2][0]; ax.plot(dm.mytime, dm.myvoltage, label=label, alpha=alpha, zorder=zorder)
        # ax = axeses[2][0]; ax.plot(dm.mytime, dm.mycurrent, label=label, alpha=alpha, zorder=zorder)

        return info, torque_average, normalized_torque_ripple, sfv.ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle

    @staticmethod
    def read_csv_results_4_general_purpose(study_name, path_prefix, fea_config_dict, femm_solver, acm_variant=None):
        wp = acm_variant.user_input['winding']

        machine_type = acm_variant.user_input['target']['machine_class']

        # Read TranFEAwi2TSS results

        # logging.getLogger(__name__).debug('Look into: ' + path_prefix)

        # Torque
        basic_info = []
        time_list = []
        TorCon_list = []
        with open(os.path.join(path_prefix, study_name + '_torque.csv'), 'r') as f:
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
        with open(os.path.join(path_prefix, study_name + '_force.csv'), 'r') as f:
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
        with open(os.path.join(path_prefix, study_name + '_circuit_current.csv'), 'r') as f:
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
        flux_csv_path = os.path.join(path_prefix, study_name + '_flux_of_fem_coil.csv')
        if os.path.exists(flux_csv_path):
            with open(flux_csv_path, 'r') as f:
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
        else:
            print(f"Warning: FEM Coil Flux CSV not found at {flux_csv_path}, skipping flux linkage processing.")

        # Displacement angle (2022)
        DisplacementAngle_list = []
        with open(os.path.join(path_prefix, study_name + '_total_rotational_displacement.csv'), 'r') as f:
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
            fname = os.path.join(path_prefix, study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
            # print 'Terminal Voltage - look into:', fname
            if os.path.exists(fname):
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
        with open(os.path.join(path_prefix, study_name + '_iron_loss_loss.csv'), 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if 'IM' in machine_type:
                    if count == 7:  
                        if 'rotorCore' not in row[2] or 'statorCore' not in row[3]:
                            raise
                    if count>8:
                        rotor_iron_loss = float(row[2]) # Rotor Core
                        stator_iron_loss = float(row[3]) # Stator Core
                        logger = logging.getLogger(__name__)
                        logger.info('Iron loss: %s %s', stator_iron_loss, rotor_iron_loss)
                        break
                elif 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type or 'CSPPM' in machine_type:
                    if count == 7:  # Validate header row
                        if 'rotorCore' not in row[1] or 'statorCore' not in row[2]:
                            if 'rotorCore' not in row[1]:
                                raise ValueError(f'Expected "rotorCore" at column index 1 in header, but found: {row[1] if len(row) > 1 else "missing"}')
                            if 'statorCore' not in row[2]:
                                raise ValueError(f'Expected "statorCore" at column index 2 in header, but found: {row[2] if len(row) > 2 else "missing"}')
                    if count>7:
                        if float(row[0]) != 0.0:
                            logger = logging.getLogger(__name__)
                            logger.debug('This should be 0: %s. This is likely due to you set up cases in your JMAG project. This automation supports case number of 1 only.', float(row[0]))
                        rotor_iron_loss = float(row[1]) # Rotor Core
                        stator_iron_loss = float(row[2]) # Stator Core
                        logger = logging.getLogger(__name__)
                        logger.info('Iron loss: %s W (Stator) and %s (Rotor)', stator_iron_loss, rotor_iron_loss)
                        break
        with open(os.path.join(path_prefix, study_name + '_joule_loss_loss.csv'), 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if 'IM' in machine_type:
                    if count>8:
                        rotor_eddycurrent_loss  = float(row[2]) # Rotor Core
                        stator_eddycurrent_loss = float(row[3]) # Stator Core
                        logger = logging.getLogger(__name__)
                        logger.info('Eddy current loss: %s %s', stator_eddycurrent_loss, rotor_eddycurrent_loss)
                        break
                elif 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type or 'CSPPM' in machine_type:
                    if count == 7:  # Validate header row
                        if 'rotorCore' not in row[1] or 'statorCore' not in row[2]:
                            if 'rotorCore' not in row[1]:
                                raise ValueError(f'Expected "rotorCore" at column index 1 in header, but found: {row[1] if len(row) > 1 else "missing"}')
                            if 'statorCore' not in row[2]:
                                raise ValueError(f'Expected "statorCore" at column index 2 in header, but found: {row[2] if len(row) > 2 else "missing"}')
                    if count>7:
                        rotor_eddycurrent_loss  = float(row[1]) # Rotor Core
                        stator_eddycurrent_loss = float(row[2]) # Stator Core
                        logger = logging.getLogger(__name__)
                        logger.info('Eddy current loss: %s %s', stator_eddycurrent_loss, rotor_eddycurrent_loss)
                        break
        with open(os.path.join(path_prefix, study_name + '_hysteresis_loss_loss.csv'), 'r') as f:
            count = 0
            for row in utility.csv_row_reader(f):
                count +=1
                if 'IM' in machine_type:
                    if count>8:
                        rotor_hysteresis_loss = float(row[2]) # Rotor Core
                        stator_hysteresis_loss = float(row[3]) # Stator Core
                        logger = logging.getLogger(__name__)
                        logger.info('Hysteresis loss: %s %s', stator_hysteresis_loss, rotor_hysteresis_loss)
                        break
                elif 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type or 'CSPPM' in machine_type:
                    if count == 7:  # Validate header row
                        if 'rotorCore' not in row[1] or 'statorCore' not in row[2]:
                            if 'rotorCore' not in row[1]:
                                raise ValueError(f'Expected "rotorCore" at column index 1 in header, but found: {row[1] if len(row) > 1 else "missing"}')
                            if 'statorCore' not in row[2]:
                                raise ValueError(f'Expected "statorCore" at column index 2 in header, but found: {row[2] if len(row) > 2 else "missing"}')
                    if count>7:
                        rotor_hysteresis_loss  = float(row[1]) # Rotor Core
                        stator_hysteresis_loss = float(row[2]) # Stator Core
                        logger = logging.getLogger(__name__)
                        logger.info('Hysteresis loss: %s %s', stator_hysteresis_loss, rotor_hysteresis_loss)
                        break

        # Joule Loss (Copper and Magnet)
        rotor_Joule_loss_list = []
        with open(os.path.join(path_prefix, study_name + '_joule_loss.csv'), 'r') as f:
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

                elif 'PM' in machine_type:
                    if count == 7: # 少一个slip变量，所以不是8，是7。
                        header_row = row
                        for idx_coil, h in enumerate(header_row):
                            if 'Coil' in h:
                                break

                    if count>7:
                        if count==7+1:
                            if 'Coil' not in header_row[idx_coil]:
                                print('[utility.py]', header_row)
                                raise Exception('Error when load csv data for Coil.')
                            stator_copper_loss = float(row[idx_coil]) # Coil # it is the same over time, this value does not account for end coil

                        if 'Magnet' not in header_row[idx_coil-1]:
                            print('[utility.py]', header_row)
                            raise Exception('Error when load csv data for Magnet.')
                        rotor_Joule_loss_list.append(float(row[idx_coil-1])) # Magnet

        # use the last 1/4 period data to compute average copper loss of Tran2TSS rather than use that of Freq study

        # 是指时域波形中已经让涡流（转子导体包括铜棒和永磁体）达到稳态的部分。
        # effective_part = rotor_Joule_loss_list[-int(0.5*fea_config_dict['designer.number_of_steps_2ndTSS']):] # number_of_steps_2ndTSS = steps for half peirod
        effective_part = rotor_Joule_loss_list[-int(fea_config_dict['designer.number_of_steps_2ndTSS']):]
        if len(effective_part) == 0: # there is no rotor eddy current (e.g., FSPM motor)
            rotor_Joule_loss = 0.0
            print('[JMAG.py] csv results: effective_part is []')
            # print(rotor_Joule_loss_list)
            # print(effective_part)
        else:
            rotor_Joule_loss = sum(effective_part) / len(effective_part)
        if 'PM' in machine_type:
            logger = logging.getLogger(__name__)
            logger.info('Magnet Joule loss: %s', rotor_Joule_loss)

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


            geo = acm_variant.user_input['geometry']
            win = wp
            copper_loss_parameters = [
                geo['d_air_gap'],
                geo['w_tooth'],
                wp['number_parallel_branch'],
                win['wires_per_slot'],
                wp['coil_pitch_y'],
                win['slot_count'],
                win['l_stack'],
                win['drive_winding_current'] + win['bearing_winding_current'], # total current amplitude
                geo['r_rotor_outer'],       # mm
                geo['r_stator_outer'] * 2 * 1e-3 # m, stator_yoke_diameter_Dsyi
            ]
            # slot_area_utilizing_ratio = (acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp) / acm_variant.CurrentAmp_per_phase
            # if slot_area_utilizing_ratio < 1:
            #     print('Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio, 'which means you are simulating a separate winding? If not, contrats--you found a bug...')
            #     print('DW, BW, Total:', acm_variant.DriveW_CurrentAmp, acm_variant.BeariW_CurrentAmp, acm_variant.CurrentAmp_per_phase)
            s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = utility.get_copper_loss_Bolognani(
                1.0 * acm_variant.geometry.slot_area * 1e-6, # win['separate_winding_utilization_ratio_for_torque'] = 1.0
                copper_loss_parameters=copper_loss_parameters, 
                STATOR_SLOT_FILL_FACTOR=win['fill_factor'],
                TEMPERATURE_OF_COIL=acm_variant.user_input['material']['magnet_temperature'])
            # s, r, sAlongStack, rAlongStack, Js, Jr = 0, 0, 0, 0, 0, 0

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

            def terminal_voltage(self, which='GroupBDW'):
                for w in [which, 'GroupACW']:
                    for suffix in [' [Case 1]', '']:
                        k = 'Terminal%s%s' % (w, suffix)
                        if k in self.Current_dict:
                            return self.Current_dict[k]
                print(f"Warning: KeyError bypassed for terminal_voltage. Keys like 'Terminal{which}/GroupACW' not found. Using default zeros array.")
                return [0.0] * len(self.Current_dict.get('Time(s)', []))

            def coil_fluxLinkage(self, which='GroupBDW'):
                for w in [which, 'GroupACW']:
                    for suffix in [' [Case 1]', '']:
                        k = 'CircuitCoil%s%s' % (w, suffix)
                        if k in self.FluxLinkage_dict:
                            return self.FluxLinkage_dict[k]
                print(f"Warning: KeyError bypassed for coil_fluxLinkage. Keys like 'CircuitCoil{which}/GroupACW' not found. Using default zeros array.")
                return [0.0] * len(self.Current_dict.get('Time(s)', []))

            def circuit_current(self, which='GroupBDW'):
                for w in [which, 'GroupACW']:
                    for suffix in [' [Case 1]', '']:
                        k = 'CircuitCoil%s%s' % (w, suffix)
                        if k in self.Current_dict:
                            return self.Current_dict[k]
                # Fallback to Default
                if 'CircuitCoilDefault' in self.Current_dict:
                    return self.Current_dict['CircuitCoilDefault']
                print(f"Warning: KeyError bypassed for circuit_current. Keys like 'CircuitCoil{which}/GroupACW/Default' not found. Using default zeros array.")
                return [0.0] * len(self.Current_dict.get('Time(s)', []))

            def get_voltage_and_current(self, number_of_steps_at_steady_state):

                # key = 'GroupBDW' if bool_DPNVorSEPA == True else 'Torque'
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

    @staticmethod
    def read_csv_results_multiple_cases(study_name, path_prefix, fea_config_dict, femm_solver, acm_variant=None):
        import copy
        import numpy as np
        import os
        import logging
        import utility
        import math
        wp = acm_variant.user_input['winding']
        machine_type = acm_variant.user_input['target']['machine_class']

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
                self.Current_dict = {}
                self.FluxLinkage_dict = {}
                self.key_list = []
                
            def set_derived(self, jmag_loss, femm_loss, Vol_Cu):
                self.jmag_loss_list = jmag_loss
                self.femm_loss_list = femm_loss
                self.Vol_Cu = Vol_Cu

            def unpack(self, bool_more_info=False):
                if bool_more_info:
                    return self.basic_info, self.time_list, self.TorCon_list, self.ForConX_list, self.ForConY_list, self.ForConAbs_list, \
                        self.DisplacementAngle_list, \
                        self.circuit_current(which='GroupACU'), self.circuit_current(which='GroupACV'), self.circuit_current(which='GroupACW'), \
                        self.circuit_current(which='GroupBDU'), self.circuit_current(which='GroupBDV'), self.circuit_current(which='GroupBDW'), \
                        self.terminal_voltage(which='GroupACU'), self.terminal_voltage(which='GroupACV'), self.terminal_voltage(which='GroupACW'), \
                        self.terminal_voltage(which='GroupBDU'), self.terminal_voltage(which='GroupBDV'), self.terminal_voltage(which='GroupBDW'), \
                        self.coil_fluxLinkage(which='GroupACU'), self.coil_fluxLinkage(which='GroupACV'), self.coil_fluxLinkage(which='GroupACW'), \
                        self.coil_fluxLinkage(which='GroupBDU'), self.coil_fluxLinkage(which='GroupBDV'), self.coil_fluxLinkage(which='GroupBDW')
                else:
                    return self.basic_info, self.time_list, self.TorCon_list, self.ForConX_list, self.ForConY_list, self.ForConAbs_list

            def terminal_voltage(self, which='GroupBDW'):
                for w in [which, 'GroupACW']:
                    for suffix in [' [Case 1]', '']:
                        k = 'Terminal%s%s' % (w, suffix)
                        if k in self.Current_dict:
                            return self.Current_dict[k]
                return [0.0] * len(self.Current_dict.get('Time(s)', []))

            def coil_fluxLinkage(self, which='GroupBDW'):
                for w in [which, 'GroupACW']:
                    for suffix in [' [Case 1]', '']:
                        k = 'CircuitCoil%s%s' % (w, suffix)
                        if k in self.FluxLinkage_dict:
                            return self.FluxLinkage_dict[k]
                return [0.0] * len(self.Current_dict.get('Time(s)', []))

            def circuit_current(self, which='GroupBDW'):
                for w in [which, 'GroupACW']:
                    for suffix in [' [Case 1]', '']:
                        k = 'CircuitCoil%s%s' % (w, suffix)
                        if k in self.Current_dict:
                            return self.Current_dict[k]
                if 'CircuitCoilDefault' in self.Current_dict:
                    return self.Current_dict['CircuitCoilDefault']
                return [0.0] * len(self.Current_dict.get('Time(s)', []))

            def get_voltage_and_current(self, number_of_steps_at_steady_state):
                key = 'Default'
                mytime  = self.Current_dict['Time(s)'][-number_of_steps_at_steady_state:]
                voltage = self.terminal_voltage()[-number_of_steps_at_steady_state:]
                current = self.circuit_current(which=key)[-number_of_steps_at_steady_state:]
                self.myvoltage = voltage
                self.mycurrent = current
                self.mytime    = mytime

            def power_factor(self, number_of_steps_at_steady_state, targetFreq=1e3, numPeriodicalExtension=1000):
                self.get_voltage_and_current(number_of_steps_at_steady_state)
                power_factor, u, i, phase_diff_ui = utility.compute_power_factor_from_half_period(self.myvoltage, self.mycurrent, self.mytime, targetFreq=targetFreq, numPeriodicalExtension=numPeriodicalExtension)
                self.ui_info = [power_factor, u, i, phase_diff_ui]
                return power_factor


        def get_col_case_map(filepath):
            col_to_case = {}
            case_ids = []
            with open(filepath, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count == 3:
                        for i, val in enumerate(row):
                            if i == 0: continue
                            if val:
                                col_to_case[i] = val
                                if val not in case_ids:
                                    case_ids.append(val)
                        break
            return case_ids, col_to_case

        torque_csv_path = os.path.join(path_prefix, study_name + '_torque.csv')
        case_ids, col_map = get_col_case_map(torque_csv_path)
        
        dm_dict = {cid: data_manager() for cid in case_ids}

        # 1. Torque
        with open(torque_csv_path, 'r') as f:
            for count, row in enumerate(utility.csv_row_reader(f), 1):
                if count <= 8:
                    pass
                else:
                    t = float(row[0])
                    for cid in case_ids:
                        if len(dm_dict[cid].time_list) == 0 or dm_dict[cid].time_list[-1] != t:
                            dm_dict[cid].time_list.append(t)
                    for i, val in enumerate(row):
                        if i == 0: continue
                        if i in col_map:
                            dm_dict[col_map[i]].TorCon_list.append(float(val))

        # 2. Force
        force_csv_path = os.path.join(path_prefix, study_name + '_force.csv')
        if os.path.exists(force_csv_path):
            _, col_map_force = get_col_case_map(force_csv_path)
            header_map_force = {}
            with open(force_csv_path, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count <= 8:
                        if count == 8:
                            for i, val in enumerate(row): header_map_force[i] = val
                    else:
                        for i, val in enumerate(row):
                            if i == 0: continue
                            if i in col_map_force:
                                cid = col_map_force[i]
                                if '1st' in header_map_force.get(i, ''):
                                    dm_dict[cid].ForConX_list.append(float(val))
                                elif '2nd' in header_map_force.get(i, ''):
                                    dm_dict[cid].ForConY_list.append(float(val))
            for cid in case_ids:
                if len(dm_dict[cid].ForConX_list) > 0:
                    dm_dict[cid].ForConAbs_list = np.sqrt(np.array(dm_dict[cid].ForConX_list)**2 + np.array(dm_dict[cid].ForConY_list)**2).tolist()

        # 3. Current
        current_csv_path = os.path.join(path_prefix, study_name + '_circuit_current.csv')
        if os.path.exists(current_csv_path):
            _, col_map_current = get_col_case_map(current_csv_path)
            header_map_current = {}
            with open(current_csv_path, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count <= 8:
                        if count == 8:
                            for i, val in enumerate(row):
                                header_map_current[i] = val
                                if i == 0:
                                    for cid in case_ids: dm_dict[cid].Current_dict[val] = []
                                elif i in col_map_current:
                                    cid = col_map_current[i]
                                    dm_dict[cid].Current_dict[val] = []
                                    dm_dict[cid].key_list.append(val)
                    else:
                        for i, val in enumerate(row):
                            if i == 0:
                                for cid in case_ids: dm_dict[cid].Current_dict[header_map_current[0]].append(float(val))
                            elif i in col_map_current:
                                cid = col_map_current[i]
                                dm_dict[cid].Current_dict[header_map_current[i]].append(float(val))
            for cid in case_ids:
                if dm_dict[cid].key_list:
                    dm_dict[cid].Current_dict['CircuitCoilDefault'] = dm_dict[cid].Current_dict[dm_dict[cid].key_list[-1]]

        # 4. Flux Linkage
        flux_csv_path = os.path.join(path_prefix, study_name + '_flux_of_fem_coil.csv')
        if os.path.exists(flux_csv_path):
            _, col_map_flux = get_col_case_map(flux_csv_path)
            header_map_flux = {}
            with open(flux_csv_path, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count <= 8:
                        if count == 8:
                            for i, val in enumerate(row):
                                header_map_flux[i] = val
                                if i == 0:
                                    for cid in case_ids: dm_dict[cid].FluxLinkage_dict[val] = []
                                elif i in col_map_flux:
                                    cid = col_map_flux[i]
                                    dm_dict[cid].FluxLinkage_dict[val] = []
                    else:
                        for i, val in enumerate(row):
                            if i == 0:
                                for cid in case_ids: dm_dict[cid].FluxLinkage_dict[header_map_flux[0]].append(float(val))
                            elif i in col_map_flux:
                                cid = col_map_flux[i]
                                dm_dict[cid].FluxLinkage_dict[header_map_flux[i]].append(float(val))
            for cid in case_ids:
                keys = list(dm_dict[cid].FluxLinkage_dict.keys())
                if len(keys) > 1:
                    dm_dict[cid].FluxLinkage_dict['CircuitCoilDefault'] = dm_dict[cid].FluxLinkage_dict[keys[-1]]

        # 5. Displacement angle
        disp_csv_path = os.path.join(path_prefix, study_name + '_total_rotational_displacement.csv')
        if os.path.exists(disp_csv_path):
            _, col_map_disp = get_col_case_map(disp_csv_path)
            for cid in case_ids: dm_dict[cid].DisplacementAngle_list = []
            with open(disp_csv_path, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count > 8:
                        for i, val in enumerate(row):
                            if i == 0: continue
                            if i in col_map_disp:
                                cid = col_map_disp[i]
                                dm_dict[cid].DisplacementAngle_list.append(float(val))

        # 6. Terminal Voltage
        voltage_fname = os.path.join(path_prefix, study_name + "_circuit_voltage.csv")
        if not os.path.exists(voltage_fname):
            voltage_fname = os.path.join(path_prefix, study_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
        if os.path.exists(voltage_fname):
            _, col_map_volt = get_col_case_map(voltage_fname)
            header_map_volt = {}
            with open(voltage_fname, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count <= 8:
                        if 'Time' in row[0]:
                            header_line = count
                            for i, val in enumerate(row):
                                header_map_volt[i] = val
                                if i == 0: pass
                                elif i in col_map_volt:
                                    cid = col_map_volt[i]
                                    dm_dict[cid].Current_dict[val] = []
                    elif 'header_line' in locals() and count > header_line:
                        for i, val in enumerate(row):
                            if i == 0: pass
                            elif i in col_map_volt:
                                cid = col_map_volt[i]
                                dm_dict[cid].Current_dict[header_map_volt[i]].append(float(val))

        # 7. Losses
        for cid in case_ids:
            dm_dict[cid].stator_iron_loss = 0.0
            dm_dict[cid].rotor_iron_loss = 0.0
            dm_dict[cid].stator_eddycurrent_loss = 0.0
            dm_dict[cid].rotor_eddycurrent_loss = 0.0
            dm_dict[cid].stator_hysteresis_loss = 0.0
            dm_dict[cid].rotor_hysteresis_loss = 0.0
            
        def parse_loss(filename, stator_key, rotor_key):
            fname = os.path.join(path_prefix, study_name + filename)
            if not os.path.exists(fname): return
            _, col_map_loss = get_col_case_map(fname)
            header_map_loss = {}
            with open(fname, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count <= 8:
                        if 'Frequency' in row[0] or 'Time' in row[0]:
                            for i, val in enumerate(row): header_map_loss[i] = val
                    else:
                        for i, val in enumerate(row):
                            if i == 0: continue
                            if i in col_map_loss:
                                cid = col_map_loss[i]
                                if 'statorCore' in header_map_loss[i]:
                                    val_stator = getattr(dm_dict[cid], stator_key)
                                    setattr(dm_dict[cid], stator_key, val_stator + float(val))
                                elif 'rotorCore' in header_map_loss[i]:
                                    val_rotor = getattr(dm_dict[cid], rotor_key)
                                    setattr(dm_dict[cid], rotor_key, val_rotor + float(val))

        parse_loss('_iron_loss_loss.csv', 'stator_iron_loss', 'rotor_iron_loss')
        parse_loss('_joule_loss_loss.csv', 'stator_eddycurrent_loss', 'rotor_eddycurrent_loss')
        parse_loss('_hysteresis_loss_loss.csv', 'stator_hysteresis_loss', 'rotor_hysteresis_loss')

        # 8. Joule Loss
        joule_csv_path = os.path.join(path_prefix, study_name + '_joule_loss.csv')
        if os.path.exists(joule_csv_path):
            _, col_map_joule = get_col_case_map(joule_csv_path)
            header_map_joule = {}
            for cid in case_ids:
                dm_dict[cid].stator_copper_loss = 0.0
                dm_dict[cid].rotor_Joule_loss_list = []
            with open(joule_csv_path, 'r') as f:
                for count, row in enumerate(utility.csv_row_reader(f), 1):
                    if count <= 8:
                        if count == 8:
                            for i, val in enumerate(row): header_map_joule[i] = val
                    else:
                        for i, val in enumerate(row):
                            if i == 0: continue
                            if i in col_map_joule:
                                cid = col_map_joule[i]
                                if 'Coil' in header_map_joule.get(i, '') or 'Coils' in header_map_joule.get(i, '') or 'statorCore' in header_map_joule.get(i, ''):
                                    if 'Coil' in header_map_joule.get(i, '') or 'Coils' in header_map_joule.get(i, ''):
                                        if dm_dict[cid].stator_copper_loss == 0.0:
                                            dm_dict[cid].stator_copper_loss = float(val)
                                elif 'Magnet' in header_map_joule.get(i, '') or 'Cage' in header_map_joule.get(i, ''):
                                    dm_dict[cid].rotor_Joule_loss_list.append(float(val))

        # Combine JMAG losses
        for cid in case_ids:
            try:
                effective_part = dm_dict[cid].rotor_Joule_loss_list[-int(fea_config_dict['designer.number_of_steps_2ndTSS']):]
                rotor_Joule_loss = sum(effective_part) / len(effective_part) if effective_part else 0.0
            except:
                rotor_Joule_loss = 0.0
            dm_dict[cid].jmag_loss_list = [
                dm_dict[cid].stator_copper_loss,
                rotor_Joule_loss,
                dm_dict[cid].stator_iron_loss + dm_dict[cid].rotor_iron_loss,
                dm_dict[cid].stator_eddycurrent_loss + dm_dict[cid].rotor_eddycurrent_loss,
                dm_dict[cid].stator_hysteresis_loss + dm_dict[cid].rotor_hysteresis_loss
            ]

            if femm_solver is not None:
                try:
                    s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = femm_solver.get_copper_loss_Bolognani(
                        acm_variant.slot_area_utilizing_ratio * femm_solver.stator_slot_area,
                        femm_solver.rotor_slot_area,
                        total_CurrentAmp=acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp)
                except Exception as e:
                    raise e
            else:
                geo = acm_variant.user_input['geometry']
                win = wp
                copper_loss_parameters = [
                    geo['d_air_gap'],
                    geo['w_tooth'],
                    wp['number_parallel_branch'],
                    win['wires_per_slot'],
                    wp['coil_pitch_y'],
                    win['slot_count'],
                    win['l_stack'],
                    win['drive_winding_current'] + win['bearing_winding_current'],
                    geo['r_rotor_outer'],
                    geo['r_stator_outer'] * 2 * 1e-3
                ]
                s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = utility.get_copper_loss_Bolognani(
                    1.0 * acm_variant.geometry.slot_area * 1e-6,
                    copper_loss_parameters=copper_loss_parameters,
                    STATOR_SLOT_FILL_FACTOR=win['fill_factor'],
                    TEMPERATURE_OF_COIL=acm_variant.user_input['material']['magnet_temperature'])
            
            dm_dict[cid].femm_loss_list = [s, r, sAlongStack, rAlongStack, Js, Jr]
            dm_dict[cid].Vol_Cu = Vol_Cu
        
        ordered_dms = [dm_dict[str(x)] for x in sorted([int(c) for c in case_ids])]
        return ordered_dms

    def build_str_results_for_multiple_cases(self, acm_variant, project_name, tran_study_name, path2FEACsv, fea_config_dict, femm_solver=None):
        import numpy as np
        import math
        import utility
        import logging
        wp = acm_variant.user_input['winding']
        machine_type = acm_variant.user_input['target']['machine_class']

        dm_list = self.read_csv_results_multiple_cases(tran_study_name, path2FEACsv, fea_config_dict, femm_solver, acm_variant)
        
        out_f1, out_f2, out_f3, out_FRW = [], [], [], []
        out_normalized_torque_ripple, out_normalized_force_error_magnitude, out_force_error_angle = [], [], []
        out_power_factor, out_rated_ratio, out_rated_stack_length_mm, out_rated_total_loss = [], [], [], []
        out_rated_stator_copper_loss_along_stack, out_rated_magnet_Joule_loss, out_rated_rotor_copper_loss_along_stack = [], [], []
        out_stator_copper_loss_in_end_turn, out_rotor_copper_loss_in_end_turn, out_rated_iron_loss, out_rated_windage_loss = [], [], [], []
        out_str_results = []
        out_coil_flux_linkage_peak2peak_value = []
        out_TRV, out_Cost, out_Cost_Fe, out_Cost_Cu, out_Cost_PM = [], [], [], [], []
        out_ss_avg_force_magnitude, out_rotor_weight, out_torque_average = [], [], []
        cost_function_O1_list, cost_function_O2_list = [], []

        for dm in dm_list:
            self.dm = dm
            coil_flux_linkage_peak2peak_value_results = []
            for k, v in dm.FluxLinkage_dict.items():
                if v: coil_flux_linkage_peak2peak_value_results.append(max(v) - min(v))
            if coil_flux_linkage_peak2peak_value_results:
                coil_flux_linkage_peak2peak_value = np.average(coil_flux_linkage_peak2peak_value_results)
            else:
                coil_flux_linkage_peak2peak_value = 0.0
                
            if fea_config_dict['designer.number_cycles_in_3rdTSS'] == 0:
                number_of_steps_at_steady_state = fea_config_dict['designer.number_of_steps_2ndTSS']
            else:
                number_of_steps_3rdTSS = fea_config_dict['designer.number_cycles_in_3rdTSS']*fea_config_dict['designer.StepPerCycle_3rdTSS']
                number_of_steps_at_steady_state = number_of_steps_3rdTSS
            dm.number_of_steps_at_steady_state = number_of_steps_at_steady_state

            basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
            sfv = utility.suspension_force_vector(ForConX_list, ForConY_list, range_ss=number_of_steps_at_steady_state)
            
            str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle = self.add_plots(None, dm, title=tran_study_name, label='Transient FEA', zorder=8, time_list=time_list, sfv=sfv, torque=TorCon_list, range_ss=sfv.range_ss)
            
            if fea_config_dict['delete_results_after_calculation'] == False:
                try: power_factor = dm.power_factor(number_of_steps_at_steady_state, targetFreq=wp['excitation_frequency_simulated'])
                except: power_factor = 0.0
            else:
                power_factor = 0.0

            geom = acm_variant.user_input['geometry']
            wind = wp
            
            r_ro = geom.get('r_rotor_outer', 4.0)
            r_ri = geom.get('r_shaft', 1.0)
            l_st = wind.get('l_stack', 16.0)
            rated_speed = wind.get('rated_speed', 20000.0)
            
            rotor_volume = math.pi * ((r_ro*1e-3)**2 - (r_ri*1e-3)**2) * (l_st*1e-3)
            rotor_weight = rotor_volume * 7650.0  
            shaft_power  = rated_speed/60. * 2*math.pi * torque_average

            if 'PMSM' in machine_type or 'SPMSM' in machine_type:
                magnet_Joule_loss = dm.jmag_loss_list[1]
                if len(dm.femm_loss_list) > 0 and dm.femm_loss_list[0] is not None:
                    copper_loss = dm.femm_loss_list[0] + magnet_Joule_loss
                else:
                    copper_loss = magnet_Joule_loss
                iron_loss = dm.jmag_loss_list[2]
            else:
                copper_loss = dm.jmag_loss_list[0] + dm.jmag_loss_list[1]
                iron_loss = dm.jmag_loss_list[2]

            windage_loss = 0.0
            total_loss   = copper_loss + iron_loss + windage_loss
            efficiency   = shaft_power / (total_loss + shaft_power) if total_loss + shaft_power > 0 else 0

            if 'PM' in machine_type or 'SPMSM' in machine_type:
                stator_copper_loss_along_stack = dm.femm_loss_list[2]
                magnet_Joule_loss = dm.jmag_loss_list[1]
                rotor_copper_loss_along_stack = 0.0
                stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack
                rotor_copper_loss_in_end_turn = 0
            else:
                stator_copper_loss_along_stack = dm.femm_loss_list[2]
                magnet_Joule_loss = 0.0
                rotor_copper_loss_along_stack  = dm.femm_loss_list[3]
                stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack
                rotor_copper_loss_in_end_turn  = dm.femm_loss_list[1] - rotor_copper_loss_along_stack

            rated_power = acm_variant.user_input['target']['rated_power']
            if rated_power is None:
                rated_power = wind.get('rated_power', 250.0)

            required_torque = rated_power / (2*math.pi*rated_speed)*60
            
            rated_ratio = required_torque / (torque_average if torque_average != 0 else 1e-9)
            rated_stack_length_mm = rated_ratio * l_st
            rated_stator_copper_loss_along_stack = rated_ratio * stator_copper_loss_along_stack
            rated_magnet_Joule_loss = rated_ratio * magnet_Joule_loss
            rated_rotor_copper_loss_along_stack = rated_ratio * rotor_copper_loss_along_stack
            rated_iron_loss = rated_ratio * dm.jmag_loss_list[2]
            rated_windage_loss = 0.0

            rated_total_loss = rated_stator_copper_loss_along_stack + rated_magnet_Joule_loss + rated_rotor_copper_loss_along_stack + stator_copper_loss_in_end_turn + rotor_copper_loss_in_end_turn + rated_iron_loss + rated_windage_loss

            rated_shaft_power  = rated_power
            rated_efficiency   = rated_shaft_power / (rated_total_loss + rated_shaft_power) if rated_total_loss + rated_shaft_power > 0 else 0
            rated_rotor_volume = math.pi * ((r_ro*1e-3)**2 - (r_ri*1e-3)**2) * (rated_stack_length_mm*1e-3)
            TRV = required_torque / rated_rotor_volume

            r_so = geom.get('r_stator_outer', 40.0)
            Vol_Fe = ( math.pi*(r_so*1e-3)**2 - math.pi*(r_ri*1e-3)**2 ) * (rated_stack_length_mm*1e-3)
            
            price_per_volume_steel    = 0.28  * 61023.744
            price_per_volume_copper   = 1.2   * 61023.744
            price_per_volume_magnet   = 11.61 * 61023.744

            if 'PM' in machine_type or 'SPMSM' in machine_type:
                magnet_area = getattr(acm_variant.geometry, 'magnet_area', geom.get('magnet_area', 100.0))
                Vol_PM = (magnet_area*1e-6) * (rated_stack_length_mm*1e-3)
            else:
                Vol_PM = 0.0

            Cost = Vol_Fe * price_per_volume_steel + dm.Vol_Cu * price_per_volume_copper + Vol_PM * price_per_volume_magnet
            Cost_Fe = Vol_Fe * price_per_volume_steel
            Cost_Cu = dm.Vol_Cu * price_per_volume_copper
            Cost_PM = Vol_PM * price_per_volume_magnet

            fitness_mapping = {
                'TorqueDensity': -TRV,
                'Cost': Cost,
                'TorqueDensityOverSquireRootCopperLoss': -TRV / math.sqrt(rated_stator_copper_loss_along_stack + stator_copper_loss_in_end_turn),
                'Efficiency': -rated_efficiency,
                'TorqueRipple': normalized_torque_ripple,
                'ForceErrorMagnitude': normalized_force_error_magnitude,
                'ForceErrorAngle': force_error_angle,
                'IronLoss': rated_iron_loss,
            }
            eval_config = acm_variant.user_input['eval_config']
            f1 = fitness_mapping.get(eval_config.get("moo.fitness_OA"), 0.0)
            f2 = fitness_mapping.get(eval_config.get("moo.fitness_OB"), 0.0)
            f3 = fitness_mapping.get(eval_config.get("moo.fitness_OC"), 0.0)
            FRW = ss_avg_force_magnitude / (rotor_weight if rotor_weight != 0 else 1)
            
            cost_function_O1, list_cost_O1 = utility.compute_list_cost(utility.use_weights('O1'), rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss)
            cost_function_O2, list_cost_O2 = utility.compute_list_cost(utility.use_weights('O2'), rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss)
            
            out_f1.append(f1)
            out_f2.append(f2)
            out_f3.append(f3)
            out_FRW.append(FRW)
            out_normalized_torque_ripple.append(normalized_torque_ripple)
            out_normalized_force_error_magnitude.append(normalized_force_error_magnitude)
            out_force_error_angle.append(force_error_angle)
            out_power_factor.append(power_factor)
            out_rated_ratio.append(rated_ratio)
            out_rated_stack_length_mm.append(rated_stack_length_mm)
            out_rated_total_loss.append(rated_total_loss)
            out_rated_stator_copper_loss_along_stack.append(rated_stator_copper_loss_along_stack)
            out_rated_magnet_Joule_loss.append(rated_magnet_Joule_loss)
            out_rated_rotor_copper_loss_along_stack.append(rated_rotor_copper_loss_along_stack)
            out_stator_copper_loss_in_end_turn.append(stator_copper_loss_in_end_turn)
            out_rotor_copper_loss_in_end_turn.append(rotor_copper_loss_in_end_turn)
            out_rated_iron_loss.append(rated_iron_loss)
            out_rated_windage_loss.append(rated_windage_loss)
            out_str_results.append(str_results)
            out_coil_flux_linkage_peak2peak_value.append(coil_flux_linkage_peak2peak_value)
            out_TRV.append(TRV)
            out_Cost.append(Cost)
            out_Cost_Fe.append(Cost_Fe)
            out_Cost_Cu.append(Cost_Cu)
            out_Cost_PM.append(Cost_PM)
            out_ss_avg_force_magnitude.append(ss_avg_force_magnitude)
            out_rotor_weight.append(rotor_weight)
            out_torque_average.append(torque_average)
            cost_function_O1_list.append(cost_function_O1)
            cost_function_O2_list.append(cost_function_O2)

        gen = -1
        ind = -1

        return (cost_function_O1_list, cost_function_O2_list), out_f1, out_f2, out_f3, out_FRW, \
               out_normalized_torque_ripple, out_normalized_force_error_magnitude, out_force_error_angle, \
               project_name, machine_type, \
               gen, ind, \
               out_power_factor, out_rated_ratio, out_rated_stack_length_mm, out_rated_total_loss, \
               out_rated_stator_copper_loss_along_stack, out_rated_magnet_Joule_loss, out_rated_rotor_copper_loss_along_stack, \
               out_stator_copper_loss_in_end_turn, out_rotor_copper_loss_in_end_turn, out_rated_iron_loss, out_rated_windage_loss, \
               out_str_results, acm_variant.geometry.slot_area, out_coil_flux_linkage_peak2peak_value, \
               out_TRV, out_Cost, out_Cost_Fe, out_Cost_Cu, out_Cost_PM, \
               out_ss_avg_force_magnitude, out_rotor_weight, out_torque_average

    # def build_str_results(self, axeses, acm_variant, project_name, tran_study_name, path2FEACsv, fea_config_dict, femm_solver=None):
    def build_str_results_for_single_case(self, acm_variant, project_name, tran_study_name, path2FEACsv, fea_config_dict, femm_solver=None):
        wp = acm_variant.user_input['winding']
        print(f"DEBUG JMAG: build_str_results called with acm_variant type: {type(acm_variant)}")
        print(f"DEBUG JMAG: project_name: {project_name}")
        machine_type = acm_variant.user_input['target']['machine_class']

        try:
            self.dm = dm = self.read_csv_results_4_general_purpose(tran_study_name, path2FEACsv, fea_config_dict, femm_solver, acm_variant=acm_variant)

            # Get peak to peak coil flux linkage
            coil_flux_linkage_peak2peak_value_results = []
            for k, v in dm.FluxLinkage_dict.items():
                coil_flux_linkage_peak2peak_value_results.append( max(v) - min(v) )
                logger = logging.getLogger(__name__)
                logger.debug('%s max: %s, min: %s', k, max(v), min(v))
            if len(coil_flux_linkage_peak2peak_value_results) > 0:
                coil_flux_linkage_peak2peak_value = np.average(coil_flux_linkage_peak2peak_value_results)
            else:
                coil_flux_linkage_peak2peak_value = 0.0
            logger = logging.getLogger(__name__)
            logger.debug('coil_flux_linkage_peak2peak_value_results: %s', coil_flux_linkage_peak2peak_value_results)

        except Exception as e:
            logger = logging.getLogger(__name__)
            logger.error('Exception in build_str_results: %s', e)
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
            self.add_plots( None, dm,
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
            power_factor = dm.power_factor(number_of_steps_at_steady_state, targetFreq=wp['excitation_frequency_simulated'])
            str_results += '\n\tPF: %g' % (power_factor)


        # compute the fitness 
        geom = acm_variant.user_input['geometry']
        wind = wp
        
        r_ro = geom['r_rotor_outer']
        if r_ro is None:
            r_ro = 4.0
            logging.getLogger(__name__).warning("KeyError bypassed for 'r_rotor_outer'. Using default value: %g", r_ro)

        r_ri = geom['r_shaft']
        if r_ri is None:
            r_ri = 1.0
            logging.getLogger(__name__).warning("KeyError bypassed for 'r_shaft'/'r_rotor_inner'. Using default value: %g", r_ri)

        l_st = wind['l_stack']
        if l_st is None:
            l_st = 16.0
            logging.getLogger(__name__).warning("KeyError bypassed for 'l_stack'. Using default value: %g", l_st)

        rated_speed = wind['rated_speed']
        if rated_speed is None:
            rated_speed = 20000.0
            logging.getLogger(__name__).warning("KeyError bypassed for 'rated_speed'/'speed_rpm'. Using default value: %g", rated_speed)

        rotor_volume = math.pi * ((r_ro*1e-3)**2 - (r_ri*1e-3)**2) * (l_st*1e-3)
        rotor_weight = rotor_volume * 7650.0  # approximate steel density kg/m^3
        shaft_power  = rated_speed/60. * 2*math.pi * torque_average

        if 'IM' in acm_variant.user_input['target']['machine_class']:
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
        elif 'PM' in acm_variant.user_input['target']['machine_class'] or 'SPMSM' in acm_variant.user_input['target']['machine_class']:
            # Rotor magnet loss by JMAG
            magnet_Joule_loss = dm.jmag_loss_list[1]
            # Stator copper loss by Binder and Bolognani 2006
            if len(dm.femm_loss_list) > 0 and dm.femm_loss_list[0] is not None:
                copper_loss = dm.femm_loss_list[0] + magnet_Joule_loss
            else:
                copper_loss = magnet_Joule_loss # fallback
            iron_loss = dm.jmag_loss_list[2] 
        else:
            raise Exception('Unknown machine type:', acm_variant.user_input['target']['machine_class'])

        # TODO: Fix windage loss calculation with new dictionary
        # windage_loss = utility.get_windage_loss(acm_variant, wp['stack_length_specified'])
        windage_loss = 0.0

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
        logger = logging.getLogger(__name__)
        logger.info('Calculate the fitness for %s', acm_variant.user_input['target']['machine_class'])
        logger.info('with x_denorm_dict: %s', acm_variant.user_input)

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

        wind_dict = wp
        rated_power = acm_variant.user_input['target']['rated_power']
        if rated_power is None:
            rated_power = wind_dict['rated_power']
            if rated_power is None:
                rated_power = 250.0
                logging.getLogger(__name__).warning("KeyError bypassed for 'rated_power'. Using default value: %g", rated_power)

        rated_speed = wind_dict['rated_speed']
        if rated_speed is None:
            rated_speed = 20000.0
            logging.getLogger(__name__).warning("KeyError bypassed for 'rated_speed'/'speed_rpm'. Using default value: %g", rated_speed)

        stack_length = wp['l_stack']
        if stack_length is None:
            stack_length = 16.0
            logging.getLogger(__name__).warning("KeyError bypassed for 'l_stack'. Using default value: %g", stack_length)

        required_torque = rated_power / (2*math.pi*rated_speed)*60

        rated_ratio                          = required_torque / (torque_average if torque_average != 0 else 1e-9)
        rated_stack_length_mm                = rated_ratio * stack_length
        rated_stator_copper_loss_along_stack = rated_ratio * stator_copper_loss_along_stack
        rated_magnet_Joule_loss              = rated_ratio * magnet_Joule_loss
        rated_rotor_copper_loss_along_stack  = rated_ratio * rotor_copper_loss_along_stack
        rated_iron_loss                      = rated_ratio * dm.jmag_loss_list[2]
        rated_windage_loss                   = 0.0

        # print(acm_variant.template.SIacm_variant.winding.stack_length_specified)
        # print(rated_ratio)
        # print(torque_average)
        # print(required_torque)
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
        elif 'PMSM' in machine_type:
                                    # 基波电流幅值（在一根导体里的电流，六相逆变器中的GroupBDW相的电流，所以相当于已经考虑了并联支路数了）
            stator_current_density = dm.ui_info[2] / 1.4142135623730951 / (acm_variant.geometry.slot_area*1e-6/wp['wires_per_slot'])
            logger = logging.getLogger(__name__)
            logger.info('Data Magager: stator_current_density (GroupBDW) = %g Arms/m^2', stator_current_density)
            rotor_current_density = 0

        rated_shaft_power  = rated_power
        rated_efficiency   = rated_shaft_power / (rated_total_loss + rated_shaft_power)  # 效率计算：机械功率/(损耗+机械功率)

        rated_rotor_volume = math.pi * ((r_ro*1e-3)**2 - (r_ri*1e-3)**2) * (rated_stack_length_mm*1e-3)

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
        
        r_so = geom['r_stator_outer']
        if r_so is None:
            r_so = 40.0
            logging.getLogger(__name__).warning("KeyError bypassed for 'r_stator_outer'. Using default value: %g", r_so)
        Vol_Fe = ( math.pi*(r_so*1e-3)**2 - math.pi*(r_ri*1e-3)**2 ) * (rated_stack_length_mm*1e-3) # Option 2 (Eric)
        
        if 'PMSM' in machine_type or 'FSPM' in machine_type or 'CPPM' in machine_type or 'CSPPM' in machine_type:
            if 'FSPM' in machine_type:
                # Find statorMagnet part
                statorMagnetArea = geom['stator_magnet_area']
                if statorMagnetArea is None:
                    statorMagnetArea = 100.0
                    logging.getLogger(__name__).warning("KeyError bypassed for 'stator_magnet_area'. Using default value: %g", statorMagnetArea)
                Vol_PM = (statorMagnetArea*1e-6) * (rated_stack_length_mm*1e-3)
                logger = logging.getLogger(__name__)
                logger.info('Area_PM %s', (statorMagnetArea*1e-6))
            else:
                magnet_area = getattr(acm_variant.geometry, 'magnet_area', geom['magnet_area'])
                if magnet_area is None:
                    magnet_area = 100.0
                    logging.getLogger(__name__).warning("KeyError bypassed for 'magnet_area'. Using default value: %g", magnet_area)
                Vol_PM = (magnet_area*1e-6) * (rated_stack_length_mm*1e-3)
        else:
            Vol_PM = 0.0
        # print('[utility.py] Area_Fe', (acm_variant.template.SI['GP']['mm_r_so'].value*1e-3) ** 2)
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
        # print('[utility.py] Vol_Fe',    Vol_Fe)
        # print('[utility.py] Volume_Cu',    dm.Vol_Cu)
        # print('[utility.py] Volume_PM',    Vol_PM)
        # print(f'[utility.py] Cost_Fe: {Cost_Fe}')
        # print(f'[utility.py] Cost_Cu: {Cost_Cu}')
        # print(f'[utility.py] Cost_PM: {Cost_PM}')

        fitness_mapping = {
            'TorqueDensity': -TRV,
            'Cost': Cost,
            'TorqueDensityOverSquireRootCopperLoss': -TRV / math.sqrt(rated_stator_copper_loss_along_stack + stator_copper_loss_in_end_turn),
            'Efficiency': - rated_efficiency,
            'TorqueRipple': normalized_torque_ripple,
            'ForceErrorMagnitude': normalized_force_error_magnitude,
            'ForceErrorAngle': force_error_angle,
            'IronLoss': rated_iron_loss,
            # 'CoggingTorque': raise Exception("CoggingTorque is not implemented")
        }
        eval_config = acm_variant.user_input['eval_config']
        
        fit_OA_key = eval_config["moo.fitness_OA"]
        f1 = fitness_mapping[fit_OA_key]

        fit_OB_key = eval_config["moo.fitness_OB"]
        f2 = fitness_mapping[fit_OB_key]

        fit_OC_key = eval_config["moo.fitness_OC"]
        f3 = fitness_mapping[fit_OC_key]

        FRW = ss_avg_force_magnitude / (rotor_weight if rotor_weight != 0 else 1)
        logger = logging.getLogger(__name__)
        logger.info('FRW: %s, Rotor weight: %s, Stack length: %s, Rated stack length: %s', FRW, rotor_weight, stack_length, rated_stack_length_mm)
        rated_rotor_volume = math.pi * ((r_ro*1e-3)**2 - (r_ri*1e-3)**2) * (rated_stack_length_mm*1e-3)
        rated_rotor_weight = rated_rotor_volume * 7650.0  # approximate steel density kg/m^3
        logger.info('rated_rotor_volume: %s, rated_rotor_weight: %s', rated_rotor_volume, rated_rotor_weight)

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
                            stack_length]           # new! 在计算FRW的时候，我们只知道原来的叠长下的力，所以需要知道原来的叠长是多少。

        counter = eval_config['counter']
        popsize = eval_config['moo.popsize']
        gen = int(counter // popsize) if isinstance(counter, (int, float)) and popsize > 0 else -1
        ind = counter if isinstance(counter, (int, float)) else -1

        str_results = '\n-------\n%s-%s\n%d,%d,O1=%g,O2=%g,f1=%g,f2=%g,f3=%g\n%s\n%s\n' % (
                        project_name, acm_variant.user_input['target']['machine_class'],
                        gen, 
                        ind, 
                        cost_function_O1, cost_function_O2, f1, f2, f3,
                        str_machine_results,
                        ','.join(['%g'%(el) for el in rated_results]), # 改为输出 rated_results
                        ) + str_results

        # str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss, cost_function = results_to_be_unpacked

        # write design evaluation data to file
        with open(path2FEACsv[:-4] + 'swarm_data.txt', 'a') as f:
            f.write(str_results)

        # if fea_config_dict['use_weights'] == 'O1':
        #     cost_function = cost_function_O1
        # elif fea_config_dict['use_weights'] == 'O2':
        #     cost_function = cost_function_O2
        # else:
        #     raise Exception('Not implemented error.')

        return (cost_function_O1, cost_function_O2), f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle, \
                project_name, acm_variant.user_input['target']['machine_class'], \
                gen, \
                ind,\
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
                acm_variant.geometry.slot_area, \
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

