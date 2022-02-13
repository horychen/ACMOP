# coding:u8
# https://stackoverflow.com/questions/19156467/run-multiple-instances-of-python-script-simultaneously
# https://docs.python.org/2/library/subprocess.html#subprocess.Popen

# wait 
# https://stackoverflow.com/questions/100624/python-on-windows-how-to-wait-for-multiple-child-processes

import os
import sys
import femm
from time import time

from numpy import exp, pi # sqrt
# from numpy import savetxt, c_

id_solver = int(sys.argv[1])
number_of_parallel_solve = int(sys.argv[2])
number_of_points = int(sys.argv[3])
output_dir = sys.argv[4]
project_file_name = sys.argv[5]
print('[parasolveSlidingMesh.py] ParaSolve instance ID:', id_solver)

femm.openfemm(True) # bHide
femm.callfemm_noeval('smartmesh(0)')
# this is essential to reduce elements counts from >50000 to ~20000.
print('mi_smartmesh is off')

for index in range(id_solver, number_of_points, number_of_parallel_solve):

    if os.path.exists(project_file_name[:-4] + f'-{index:03d}.ans'):
        raise Exception('Previous .ans is not deleted!!!')

    tic = time()
    femm.opendocument(project_file_name[:-4] + f'-{index:03d}.fem')
    try:
        femm.mi_analyze(1) # None for inherited. 1 for a minimized window,
        # write_Torque_and_B_data_to_file(output_file_name[-4:], exp(1j*rotor_position)) # this function is moved to FEMM_Solver.py as keep...

        # if True:
        #     # call this after mi_analyze
        #     femm.mi_loadsolution()
        #     # Physical Amount on the Rotor
        #     femm.mo_groupselectblock(100) # rotor iron
        #     femm.mo_groupselectblock(101) # rotor bars
        #     Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
        #     Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
        #     torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
        #     femm.mo_clearblock()
        #     # write results to a data file (write to partial files to avoid compete between parallel instances)
        #     handle_torque.write("%s %g %g %g\n" % ( output_file_name[-4:], torque, Fx, Fy )) # output_file_name[-4:] = str_rotor_position
        #     # close post-process
        #     femm.mo_close()
    except Exception as error:
            print(error.args)
            raise error
        # print 'Is it: Material properties have not been defined for all regions? Check the following file:'
        # print i, fem_file_list[i]
    femm.mi_close()
    toc = time()
    print(index, project_file_name[:-4] + f'-{index:03d}.fem', toc - tic, 's')
femm.closefemm()

# handle_torque.close()
# handle_B_data.close()



def write_Torque_and_B_data_to_file(str_rotor_position, rotation_operator):
    # call this after mi_analyze
    femm.mi_loadsolution()

    # Physical Amount on the Rotor
    femm.mo_groupselectblock(100) # rotor iron
    femm.mo_groupselectblock(101) # rotor bars
    Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
    Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
    torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
    femm.mo_clearblock()
    # write results to a data file (write to partial files to avoid compete between parallel instances)
    handle_torque.write("%s %g %g %g\n" % ( str_rotor_position, torque, Fx, Fy ))    

    # Field Amount of 1/4 model (this is valid if we presume the suspension two pole field is weak)
    number_of_elements = femm.mo_numelements()
    stator_Bx_data = []
    stator_By_data = []
    stator_Area_data = []
    rotor_Bx_data = []
    rotor_By_data = []
    rotor_Area_data = []
    # one_list = []
    for id_element in range(1, number_of_elements+1):
        _, _, _, x, y, area, group = femm.mo_getelement(id_element)
        if y>0 and x>0:
            if group==10: # stator iron
                # 1. What we need for iron loss evaluation is the B waveform at a fixed point (x,y). 
                #    For example, (x,y) is the centeroid of element in stator tooth.
                Bx, By = femm.mo_getb(x, y)
                stator_Bx_data.append(Bx)
                stator_By_data.append(By)
                stator_Area_data.append(area)

            if group==100: # rotor iron
                # 2. The element at (x,y) is no longer the same element from last rotor position.
                #    To find the exact element from last rotor position,
                #    we rotate the (x,y) forward as we rotate the model (rotor), get the B value there: (x,y)*rotation_operator, and correct the (Bx,By)/rotation_operator
                complex_new_xy = (x + 1j*y) * rotation_operator
                Bx, By = femm.mo_getb( complex_new_xy.real, 
                                       complex_new_xy.imag )
                complex_new_BxBy = (Bx + 1j*By) * rotation_operator
                rotor_Bx_data.append(complex_new_BxBy.real)
                rotor_By_data.append(complex_new_BxBy.imag)
                rotor_Area_data.append(area)

            # one_list.append(sqrt(Bx**2 + By**2))
            # one_list.append(area)
    # option 1
    handle_stator_B_data.write(str_rotor_position + ',' + ','.join(['%g,%g,%g'%(Bx,By,A) for Bx,By,A in zip(stator_Bx_data, stator_By_data, stator_Area_data) ]) + '\n')
    handle_rotor_B_data.write(str_rotor_position  + ',' + ','.join(['%g,%g,%g'%(Bx,By,A) for Bx,By,A in zip(rotor_Bx_data, rotor_By_data, rotor_Area_data) ]) + '\n')

        # option 2: one_list
        # handle_B_data.write(str_rotor_position + ',' + ','.join(['%g'%(B) for B in B_data ]) + ','.join(['%g'%(A) for A in Area_data ]) + '\n')

        # numpy is slower than open().write!!!
        # tic = time()
        # # savetxt(handle_B_data, c_[one_list])
        # savetxt(handle_B_data, one_list)
        # toc = time()
        # print toc - tic, 's\n\n'

    femm.mo_close()
