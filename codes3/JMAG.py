import win32com.client
import os
import logging
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
                    raise Exception('No JMAG Designer is found.')
            if self.fea_config_dict['designer.Show'] == True:
                app.Show()
            else:
                app.Hide()
            # app.Quit()
            self.app = app # means that the JMAG Designer is turned ON now.

            def add_steel(self):
                print('[First run on %s detected]'%(self.fea_config_dict['pc_name']), self.spec_input_dict['Steel'], 'is added to jmag material library.')
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
                        print(self.fea_config_dict['pc_name'], self.spec_input_dict['Steel'])
                        if self.fea_config_dict['pc_name'] + '/' + self.spec_input_dict['Steel'] in line:
                            flag_already_there = True
                            break
                    if flag_already_there == False:
                        add_steel(self)
        else:
            app = self.app

        print(expected_project_file_path)
        if os.path.exists(expected_project_file_path):
            print('JMAG project exists already. I learned my lessions. I will NOT delete it but create a new one with a different name instead.')
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
                # debugging = list_segments(i).GetName()
                self.doc.GetSelection().Add(self.sketch.GetItem(segment.GetName()))

            self.sketch.CreateRegions()
            # self.sketch.CreateRegionsWithCleanup(EPS, True) # StatorCore will fail

            if idx == 0:
                region_object = self.sketch.GetItem('Region') # This is how you get access to the region you create.
            else:
                region_object = self.sketch.GetItem('Region.%d'%(idx+1)) # This is how you get access to the region you create.
            list_region_objects.append(region_object)
        # raise

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

    def run_study(self, im_variant, app, study, toc):
        logger = logging.getLogger(__name__)
        if self.fea_config_dict['JMAG_Scheduler'] == False:
            print('Run jam.exe...')
            # if run_list[1] == True:
            study.RunAllCases()
            msg = 'Time spent on %s is %g s.'%(study.GetName() , clock_time() - toc)
            logger.debug(msg)
            print(msg)
        else:
            print('Submit to JMAG_Scheduler...')
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

