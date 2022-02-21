# coding:u8
# https://stackoverflow.com/questions/19156467/run-multiple-instances-of-python-script-simultaneously
# https://docs.python.org/2/library/subprocess.html#subprocess.Popen

# wait 
# https://stackoverflow.com/questions/100624/python-on-windows-how-to-wait-for-multiple-child-processes

import os, sys, femm
from time import time
from pylab import np

id_solver = int(sys.argv[1])
number_of_parallel_solve = int(sys.argv[2])
ns = int(sys.argv[3]) # number of steps
output_dir = sys.argv[4]
project_file_name = sys.argv[5]
GroupSummary = {
    "stator_iron_core" : int(sys.argv[6]),
    "coils"            : int(sys.argv[7]),
    "rotor_iron_core"  : int(sys.argv[8]),
    "magnet"           : int(sys.argv[9]),
}
print('[parasolveSlidingMesh.py] ParaSolve instance ID:', id_solver)

femm.openfemm(True) # bHide
# femm.smartmesh(0)
femm.callfemm_noeval('smartmesh(0)')
# this is essential to reduce elements counts from >50000 to ~20000.
# print('mi_smartmesh is off')

bool_initialized = False
paraResults = {}
for index in range(id_solver, ns, number_of_parallel_solve):

    if os.path.exists(project_file_name[:-4] + f'-{index:03d}.ans'):
        raise Exception('Previous .ans is not deleted!!!')

    tic = time()
    femm.opendocument(project_file_name[:-4] + f'-{index:03d}.fem')
    # femm.mi_smartmesh(0)
    # femm.mi_saveas(f'temp-{id_solver}.fem')
    femm.mi_analyze(1) # None for inherited. 1 for a minimized window,
    print(f'{index:02d}', project_file_name[:-4] + f'-{index:03d}.fem | Solving time: {time() - tic:.1f} s ', end='')
    femm.mi_loadsolution()
    # femm.mo_smooth('off') # flux smoothing algorithm is off

    if bool_initialized==False:
        bool_initialized = True
        # Record the initial mesh elements if the first time through the loop
        nn = int(femm.mo_numelements())

        M  = np.zeros([ns,       22]) # 9 columns of performance data matrix

        z  = np.zeros([nn,       1], dtype=np.complex64) # Location of the centroid of each element
        a  = np.zeros([nn,       1]) # Area of each element
        g  = np.zeros([nn,       1]) # Block label of each element
        for m in range(nn): # start from 0 for indexing but from 1 for counting 
            counting_element = m+1
            elm = femm.mo_getelement(counting_element)
            # z is a vector of complex numbers that represents the location of the centroid of each element.
            z[m] = elm[4-1] + 1j*elm[5-1]
            # element area in the length units used to draw the geometry
            a[m] = elm[6-1]
            # group number associated with the element
            g[m] = elm[7-1]

        # mo_getprobleminfo returns, among other things, the depth of the
        # machine in the into-the-page direction and the length units used to
        # draw the geometry. Both of these pieces of information will be needed
        # to integrate the losses over the volume of the machine.
        if id_solver == 0:
            probinfo = femm.mo_getprobleminfo()

        A  = np.zeros([ns, nn]) # matrix that will hold the vector potential info
        b  = np.zeros([ns, nn], dtype=np.complex64) # matrix that will hold the flux density info

    # Store element flux densities B and magnetic potential A
    for m in range(nn):
        if g[m]==GroupSummary['magnet']: # Element is in a rotor magnet, marked with group numbers 11 and higher
            # Store vector potential at the element centroid for elements that are in PMs
            A[index,m] = femm.mo_geta( float(np.real(z[m])), 
                                            float(np.imag(z[m])) )
        elif g[m] == GroupSummary['stator_iron_core'] \
            or g[m] == GroupSummary['rotor_iron_core'] \
            or g[m] == GroupSummary['coils']: # Element is on the stator or rotor iron or coils
            # Store flux density at the element centroid for these elements
            b_temp = femm.mo_getb(  float(np.real(z[m])), 
                                    float(np.imag(z[m])) )
            b[index,m] = b_temp[0] + 1j*b_temp[1]

    # Torque, force, fluxes
    torque = femm.mo_gapintegral('WholeModelSlidingBand',0)
    forces = femm.mo_gapintegral('WholeModelSlidingBand',1)
    energy = femm.mo_gapintegral('WholeModelSlidingBand',2)

    cP_Uac = femm.mo_getcircuitproperties('U-GrpAC')
    cP_Vac = femm.mo_getcircuitproperties('V-GrpAC')
    cP_Wac = femm.mo_getcircuitproperties('W-GrpAC')
    cP_Ubd = femm.mo_getcircuitproperties('U-GrpBD')
    cP_Vbd = femm.mo_getcircuitproperties('V-GrpBD')
    cP_Wbd = femm.mo_getcircuitproperties('W-GrpBD')
    M[index,0] = torque
    M[index,1] = forces[0]
    M[index,2] = forces[1]
    M[index,3] = energy
    M[index,4] = cP_Uac[0]
    M[index,5] = cP_Vac[0]
    M[index,6] = cP_Wac[0]
    M[index,7] = cP_Ubd[0]
    M[index,8] = cP_Vbd[0]
    M[index,9] = cP_Wbd[0]
    M[index,10] = cP_Uac[1]
    M[index,11] = cP_Vac[1]
    M[index,12] = cP_Wac[1]
    M[index,13] = cP_Ubd[1]
    M[index,14] = cP_Vbd[1]
    M[index,15] = cP_Wbd[1]
    M[index,16] = cP_Uac[2]
    M[index,17] = cP_Vac[2]
    M[index,18] = cP_Wac[2]
    M[index,19] = cP_Ubd[2]
    M[index,20] = cP_Vbd[2]
    M[index,21] = cP_Wbd[2]
    femm.mo_close()
    femm.mi_close()
    toc = time()
    print(f'Total time: {toc - tic:.1f} s.')
femm.closefemm()

paraResults['b'] = b
paraResults['A'] = A
paraResults['M'] = M
# print(id_solver, M)
if id_solver == 0:
    paraResults['z'] = z
    paraResults['a'] = a
    paraResults['g'] = g
    paraResults['probinfo'] = probinfo
import pickle
with open(f'paraResults{id_solver}.pkl', 'wb') as pickle_file:
    pickle.dump(paraResults, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
