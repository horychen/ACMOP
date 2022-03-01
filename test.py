
''' Install Anaconda3-2021.05-Windows-x86_64 Python 3.8.8 for PYGMO to work '''

import os
os.system('cd codes3 && python acmop.py')
quit()


import os, sys
try:sys.path.insert(0, os.path.dirname(__file__)+'/codes3/')
except:sys.path.insert(0, 'D:/DrH/Codes/acmop/codes3/')
finally:import acmop

mop = acmop.AC_Machine_Optiomization_Wrapper(
    # select_spec='IM Q24p1y9 Qr32 Round Bar',
    # select_fea_config_dict = '#019 JMAG IM Nine Variables',

    select_spec            = 'PMSM Q12p4y1 PEMD-2020', #'PMSM Q18p4y2 Beijing ShiDaiChaoQun',
    select_fea_config_dict = '#02 JMAG PMSM Evaluation Setting',
    # select_fea_config_dict = '#04 FEMM PMSM Evaluation Setting',

    project_loc            = fr'../_default/',
    bool_show_GUI          = True
)
mop.part_evaluation() # Module 3


