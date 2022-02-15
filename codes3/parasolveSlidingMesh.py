# coding:u8
# https://stackoverflow.com/questions/19156467/run-multiple-instances-of-python-script-simultaneously
# https://docs.python.org/2/library/subprocess.html#subprocess.Popen

# wait 
# https://stackoverflow.com/questions/100624/python-on-windows-how-to-wait-for-multiple-child-processes

import os, sys, femm
from time import time

id_solver = int(sys.argv[1])
number_of_parallel_solve = int(sys.argv[2])
number_of_points = int(sys.argv[3])
output_dir = sys.argv[4]
project_file_name = sys.argv[5]
print('[parasolveSlidingMesh.py] ParaSolve instance ID:', id_solver)

femm.openfemm(True) # bHide
# femm.smartmesh(0)
femm.callfemm_noeval('smartmesh(0)')
# this is essential to reduce elements counts from >50000 to ~20000.
# print('mi_smartmesh is off')

for index in range(id_solver, number_of_points, number_of_parallel_solve):

    if os.path.exists(project_file_name[:-4] + f'-{index:03d}.ans'):
        raise Exception('Previous .ans is not deleted!!!')

    tic = time()
    femm.opendocument(project_file_name[:-4] + f'-{index:03d}.fem')
    # femm.mi_smartmesh(0)
    # femm.mi_saveas(f'temp-{id_solver}.fem')
    try:
        femm.mi_analyze(1) # None for inherited. 1 for a minimized window,
        # femm.mi_loadsolution()
        # femm.mo_smooth('off') # flux smoothing algorithm is off
        # nn = femm.mo_numelements()
        # print('number of element is', nn)
        # femm.mo_close()
    except Exception as error:
        print(error.args)
        raise error
        # print 'Is it: Material properties have not been defined for all regions? Check the following file:'
        # print i, fem_file_list[i]
    femm.mi_close()
    toc = time()
    print(index, project_file_name[:-4] + f'-{index:03d}.fem', toc - tic, 's')
femm.closefemm()
