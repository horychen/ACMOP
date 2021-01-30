import inner_rotor_motor
from collections import OrderedDict

class bearingless_spmsm_template(inner_rotor_motor.template_machine_as_numbers):
    ''' This is a surface mounted PM motor but it might have saliency if alpha_rm is less than 180/p.
        就是说，允许永磁体陷入转子铁芯，只是永磁体外没有铁包裹防止永磁体飞出，而是需要额外增加碳纤维套。
    '''
    def __init__(self, fea_config_dict, spec_input_dict):
        # 初始化父类
        super(bearingless_spmsm_template, self).__init__(fea_config_dict, spec_input_dict)

        # 基本信息
        self.machine_type = 'SPMSM'
        self.name = '__SPMSM'

        # 初始化搜索空间
        GP = self.d['GP']
        SD = self.SD
        childGP = OrderedDict({
            # SPMSM Peculiar
            "deg_alpha_rm"      : acmop_parameter("free",      "magnet_pole_span_angle",        None, [None, None]),
            "mm_d_rp"           : acmop_parameter("free",      "inter_polar_iron_thickness",    None, [None, None]),
            "deg_alpha_rs"      : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "magnet_segment_span_angle",     None, [None, None]),
            "mm_d_rs"           : acmop_parameter("free" if SD['no_segmented_magnets']!=1 else "fixed",   "inter_segment_iron_thickness",  None, [None, None]),
        })
        self.d.update( {"GP": childGP} )



class bearingless_spmsm_design_variant(inner_rotor_motor.variant_machine_as_objects):

    def __init__(self, spmsm_template=None, x_denorm=None, counter=None, counter_loop=None):
        # 初始化父类
        super(bearingless_spmsm_template, self).__init__(spmsm_template, x_denorm, counter, counter_loop)


































# circumferential segmented rotor 
if __name__ == '__main__':
    import JMAG
    import Location2D
    import CrossSectInnerNotchedRotor
    import CrossSectStator
    # from pylab import np

    if True:
        from utility import my_execfile
        my_execfile('./default_setting.py', g=globals(), l=locals())
        fea_config_dict
        toolJd = JMAG.JMAG(fea_config_dict)

        project_name          = 'proj%d'%(0)
        expected_project_file_path = './' + "%s.jproj"%(project_name)
        toolJd.open(expected_project_file_path)

    spmsm_template = bearingless_spmsm_template()
    spmsm_template.fea_config_dict = fea_config_dict
    Q = 6

    spmsm_template.deg_alpha_st = 360/Q*0.8   # deg_alpha_st # span angle of tooth: class type DimAngular
    spmsm_template.deg_alpha_so = 0                          # deg_alpha_so # angle of tooth edge: class type DimAngular
    spmsm_template.mm_r_si      = 50     # mm_r_si           # inner radius of stator teeth: class type DimLinear
    spmsm_template.mm_d_so      = 5      # mm_d_so           # tooth edge length: class type DimLinear
    spmsm_template.mm_d_sp      = 1.5*spmsm_template.mm_d_so # mm_d_sp      # tooth tip length: class type DimLinear
    spmsm_template.mm_d_st      = 15     # mm_d_st      # tooth base length: class type DimLinear
    spmsm_template.mm_d_sy      = 15     # mm_d_sy      # back iron thickness: class type DimLinear
    spmsm_template.mm_w_st      = 13     # mm_w_st      # tooth base width: class type DimLinear
    spmsm_template.mm_r_st      = 0         # mm_r_st      # fillet on outter tooth: class type DimLinear
    spmsm_template.mm_r_sf      = 0         # mm_r_sf      # fillet between tooth tip and base: class type DimLinear
    spmsm_template.mm_r_sb      = 0         # mm_r_sb      # fillet at tooth base: class type DimLinear
    spmsm_template.Q            = 6      # number of stator slots (integer)
    spmsm_template.sleeve_length        = 2 # mm
    spmsm_template.fixed_air_gap_length = 0.75 # mm
    spmsm_template.mm_d_pm      = 6      # mm_d_pm          # manget depth
    spmsm_template.deg_alpha_rm = 60     # deg_alpha_rm     # angular span of the pole: class type DimAngular
    spmsm_template.deg_alpha_rs = 10     # deg_alpha_rs     # segment span: class type DimAngular
    spmsm_template.mm_d_ri      = 8      # mm_d_ri          # rotor iron thickness: class type DimLinear
    spmsm_template.mm_r_ri      = 40     # mm_r_ri          # inner radius of rotor: class type DimLinear
    spmsm_template.mm_d_rp      = 5      # mm_d_rp          # interpolar iron thickness: class type DimLinear
    spmsm_template.mm_d_rs      = 3      # mm_d_rs          # inter segment iron thickness: class type DimLinear
    spmsm_template.p = 2     # p     # number of pole pairs
    spmsm_template.s = 3     # s     # number of segments  

    build_design_parameters_list(spmsm_template)

    spmsm_template.DriveWinding_Freq       = 1000
    spmsm_template.DriveWinding_Rs         = 0.1 # TODO
    spmsm_template.DriveWinding_zQ         = 1
    spmsm_template.DriveWinding_CurrentAmp = None # this means it depends on the slot area
    spmsm_template.DriveWinding_poles = 2*spmsm_template.p

    spmsm_template.Js = 4e6 # Arms/m^2
    spmsm_template.fill_factor = 0.45

    spmsm_template.stack_length = 100 # mm

    # logger = logging.getLogger(__name__) 
    # logger.info('spmsm_variant ID %s is initialized.', self.name)

    spmsm = bearingless_spmsm_design(   spmsm_template=spmsm_template,
                                        free_variables=None,
                                        counter=None, 
                                        counter_loop=None
                                        )
    # Rotor Core
    list_segments = spmsm.rotorCore.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = 0 #spmsm.rotorCore.p*2
    region1 = toolJd.prepareSection(list_segments)

    # Rotor Magnet    
    list_regions = spmsm.rotorMagnet.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = 0 #spmsm.rotorMagnet.notched_rotor.p*2
    region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')


# notched rotor 
if __name__ == '__main__':
    import JMAG
    import Location2D
    import CrossSectInnerNotchedRotor
    import CrossSectStator
    # from pylab import np

    if True:
        from utility import my_execfile
        my_execfile('./default_setting.py', g=globals(), l=locals())
        fea_config_dict
        toolJd = JMAG.JMAG(fea_config_dict, spec_input_dict=None)

        project_name          = 'proj%d'%(0)
        expected_project_file_path = './' + "%s.jproj"%(project_name)
        toolJd.open(expected_project_file_path)

    spmsm_template = bearingless_spmsm_template()
    spmsm_template.fea_config_dict = fea_config_dict
    Q = 6

    spmsm_template.deg_alpha_st = 360/Q*0.8   # deg_alpha_st # span angle of tooth: class type DimAngular
    spmsm_template.deg_alpha_so = 0                          # deg_alpha_so # angle of tooth edge: class type DimAngular
    spmsm_template.mm_r_si      = 50     # mm_r_si           # inner radius of stator teeth: class type DimLinear
    spmsm_template.mm_d_so      = 5      # mm_d_so           # tooth edge length: class type DimLinear
    spmsm_template.mm_d_sp      = 1.5*spmsm_template.mm_d_so # mm_d_sp      # tooth tip length: class type DimLinear
    spmsm_template.mm_d_st      = 15     # mm_d_st      # tooth base length: class type DimLinear
    spmsm_template.mm_d_sy      = 15     # mm_d_sy      # back iron thickness: class type DimLinear
    spmsm_template.mm_w_st      = 13     # mm_w_st      # tooth base width: class type DimLinear
    spmsm_template.mm_r_st      = 0         # mm_r_st      # fillet on outter tooth: class type DimLinear
    spmsm_template.mm_r_sf      = 0         # mm_r_sf      # fillet between tooth tip and base: class type DimLinear
    spmsm_template.mm_r_sb      = 0         # mm_r_sb      # fillet at tooth base: class type DimLinear
    spmsm_template.Q            = 6      # number of stator slots (integer)
    spmsm_template.sleeve_length        = 2 # mm
    spmsm_template.fixed_air_gap_length = 0.75 # mm
    spmsm_template.mm_d_pm      = 6      # mm_d_pm          # manget depth
    spmsm_template.deg_alpha_rm = 60     # deg_alpha_rm     # angular span of the pole: class type DimAngular
    spmsm_template.deg_alpha_rs = 60     # deg_alpha_rs     # segment span: class type DimAngular
    spmsm_template.mm_d_ri      = 8      # mm_d_ri          # inner radius of rotor: class type DimLinear
    spmsm_template.mm_r_ri      = 40     # mm_r_ri          # rotor iron thickness: class type DimLinear
    spmsm_template.mm_d_rp      = 5      # mm_d_rp          # interpolar iron thickness: class type DimLinear
    spmsm_template.mm_d_rs      = 0*3      # mm_d_rs          # inter segment iron thickness: class type DimLinear
    spmsm_template.p = 2     # p     # number of pole pairs
    spmsm_template.s = 1     # s     # number of segments  

    build_design_parameters_list(spmsm_template)

    spmsm_template.DriveWinding_Freq       = 1000
    spmsm_template.DriveWinding_Rs         = 0.1 # TODO
    spmsm_template.DriveWinding_zQ         = 1
    spmsm_template.DriveWinding_CurrentAmp = None # this means it depends on the slot area
    spmsm_template.DriveWinding_poles = 2*spmsm_template.p

    spmsm_template.Js = 4e6 # Arms/m^2
    spmsm_template.fill_factor = 0.45

    spmsm_template.stack_length = 100 # mm

    # logger = logging.getLogger(__name__) 
    # logger.info('spmsm_variant ID %s is initialized.', self.name)

    spmsm = bearingless_spmsm_design(   spmsm_template=spmsm_template,
                                        free_variables=None,
                                        counter=None, 
                                        counter_loop=None
                                        )
    # Rotor Core
    list_segments = spmsm.rotorCore.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = spmsm.rotorCore.p*2
    region1 = toolJd.prepareSection(list_segments)

    # Rotor Magnet    
    list_regions = spmsm.rotorMagnet.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = spmsm.rotorMagnet.notched_rotor.p*2
    region2 = toolJd.prepareSection(list_regions)

    # Rotor Magnet
    sleeve = CrossSectInnerNotchedRotor.CrossSectSleeve(
                    name = 'Sleeve',
                    notched_magnet = spmsm.rotorMagnet,
                    d_sleeve = spmsm_template.sleeve_length
                    )

    list_regions = sleeve.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = spmsm.rotorMagnet.notched_rotor.p*2
    regionS = toolJd.prepareSection(list_regions)


    # # Stator Core
    # list_regions = spmsm.stator_core.draw(toolJd)
    # toolJd.bMirror = True
    # toolJd.iRotateCopy = spmsm.stator_core.Q
    # region1 = toolJd.prepareSection(list_regions)

    # # Stator Winding
    # list_regions = spmsm.coils.draw(toolJd)
    # toolJd.bMirror = False
    # toolJd.iRotateCopy = spmsm.coils.stator_core.Q
    # region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')


def add_carbon_fiber_material(app):
    app.GetMaterialLibrary().CreateCustomMaterial(u"CarbonFiber", u"Custom Materials")
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"Density", 1.6)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"CoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"DemagnetizationCoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"MagnetizationSaturated", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulus", 110000)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulus", 5000)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"PoissonRatio", 0.1)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"Thermal Expansion", 8.4e-06)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulusX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulusY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"YoungModulusZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulusXY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulusYZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"ShearModulusZX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G11", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G12", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G13", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G14", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G15", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G16", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G22", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G23", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G24", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G25", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G26", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G33", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G34", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G35", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G36", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G44", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G45", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G46", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G55", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G56", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"CarbonFiber").SetValue(u"G66", 0)





    # -*- coding: utf-8 -*-
    app = designer.GetApplication()
    app.GetMaterialLibrary().CopyMaterial(u"Arnold/Reversible/N40H")
    app.GetMaterialLibrary().GetUserMaterial(u"N40H(reversible) copy").SetValue(u"Name", u"MyN40H(reversible)")
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"Density", 7.5)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"CoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"DemagnetizationCoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"MagnetizationSaturated", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"MagnetizationSaturatedMakerValue", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulus", 160000)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"PoissonRatio", 0.24)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulusX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulusY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"YoungModulusZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"ShearModulusXY", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"ShearModulusYZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"ShearModulusZX", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G11", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G12", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G13", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G14", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G15", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G16", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G22", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G23", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G24", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G25", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G26", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G33", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G34", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G35", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G36", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G44", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G45", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G46", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G55", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G56", 0)
    app.GetMaterialLibrary().GetUserMaterial(u"MyN40H(reversible)").SetValue(u"G66", 0)



