import acmop

# demo: PMSM Q18 p4
mop = acmop.AC_Machine_Optiomization_Wrapper(
    select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting',
    # select_fea_config_dict = '#04 FEMM PMSM Evaluation Setting',
    select_spec            = 'PMSM Q12p4y1 A', #'PMSM Q18p4y2 Beijing ShiDaiChaoQun',
    project_loc            = fr'D:/DrH/Codes/acmop/_default/',
    bool_show_GUI          = True
)

# quit()

#########################
# Call the five modules
#########################
# mop.part_winding() # Module 1
# mop.acm_template # Module 2 (the execution code has been moved to the end of __post_init__ of AC_Machine_Optiomization_Wrapper)
if False:
    mop.part_evaluation() # Module 3
    # mop.part_optimization(acm_template) # Module 4
else:
    if True:
        motor_design_variant = mop.reproduce_design_from_jsonpickle('p4ps5-Q12y1-0999') # Module 5 - reproduction of any design variant object
    else:
        mop.part_post_optimization_analysis(project_name='proj12-SPMSM_IDQ12p4s1') # Module 5 - visualize swarm data
        # mop.part_post_optimization_analysis(project_name='proj212-SPMSM_IDQ12p1s1') # Module 5
