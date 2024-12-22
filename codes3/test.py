
''' Install Anaconda3-2021.05-Windows-x86_64 Python 3.8.8 for PYGMO to work '''

# import os
# os.system('cd codes3 && python acmop.py')
import sys
import os
sys.path.append(os.path.abspath('../codes3'))

import JMAG

import matplotlib
matplotlib.use('TkAgg')  # 设置 Matplotlib 后端
import matplotlib.pyplot as plt
import numpy as np

def main():
    dm = 0 
    fig, axeses = plt.subplots(2, 2)
    time_list = np.linspace(0, 10, 100)
    torque = np.sin(time_list)
    force_x = np.cos(time_list)
    force_y = np.sin(time_list)
    force_abs = np.sqrt(force_x**2 + force_y**2)
    force_err_abs = np.zeros(100),  # 确保 force_err_abs 和 time_list 长度相同
    sfv = type('sfv', (object,), {
        'force_abs': force_abs,
        'force_x': force_x,
        'force_y': force_y,
        'ss_avg_force_magnitude': np.mean(force_abs),
        'normalized_force_error_magnitude': 0.1,
        'force_err_abs': force_err_abs,  # 确保 force_err_abs 和 time_list 长度相同
        'ss_max_force_err_abs': [np.max(force_err_abs), np.min(force_err_abs)],
        'force_error_angle': 5,
        'ss_max_force_err_ang': [3, -3],
        'ss_avg_force_vector': [np.mean(force_x), np.mean(force_y)],

        'force_err_ang_old_way': np.random.uniform(-5, 5, len(time_list)),
        'force_err_ang_new_way': np.random.uniform(-5, 5, len(time_list)),
        'force_ang': np.random.uniform(-5, 5, len(time_list)),
        'ss_avg_force_angle': np.mean(np.random.uniform(-5, 5, len(time_list)))
    })()
    # print(max(force_err_abs))
    JMAG.JMAG.add_plots(axeses, dm, title='Test Plot', label='Test', zorder=1, time_list=time_list, sfv=sfv, torque=torque, range_ss=20)

    plt.show()

if __name__ == '__main__':
    main()
    if force_x == 1:
        print('force_x == 1')

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


