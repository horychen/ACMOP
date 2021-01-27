from pylab import np
import matplotlib.pyplot as plt

def Get_FRW(raw):

    loc1 = raw[2].find('f1')
    loc2 = raw[2].find('f2')
    loc3 = raw[2].find('f3')
    f1 = float(raw[2][loc1+3:loc2-1])
    f2 = float(raw[2][loc2+3:loc3-1])
    f3 = float(raw[2][loc3+3:])
    Cost = f1
    Efficiency = -f2
    RipplePerformanceSum = f3

    machine_data = [float(x) for x in raw[3].split(',')]
    rated_data   = [float(x) for x in raw[4].split(',')]

    # print(raw)
    individual_ss_avg_force_magnitude               = machine_data[4]
    individual_rated_stack_length                   = rated_data[10]
    individual_original_stack_length                = rated_data[11]
    individual_rated_rotor_volume                   = rated_data[9]
    individual_original_rotor_volume                = individual_rated_rotor_volume/individual_rated_stack_length*individual_original_stack_length
    individual_rated_rotor_weight                   = (individual_rated_rotor_volume*8050*9.8)
    individual_original_rotor_weight                = individual_rated_rotor_weight/individual_rated_stack_length*individual_original_stack_length
    individual_FRW = individual_ss_avg_force_magnitude/individual_original_rotor_weight
    # print(  individual_ss_avg_force_magnitude,
    #         individual_rated_rotor_volume    ,
    #         individual_rated_rotor_weight    ,
    #         individual_rated_stack_length    ,
    #         individual_original_stack_length ,
    #         individual_original_rotor_weight ,
    #         individual_FRW, sep='\n')

    individual_torque_average = machine_data[2]
    individual_TRV = individual_torque_average / individual_original_rotor_volume
    Trip = machine_data[3]
    Em = machine_data[5]
    Ea = machine_data[6]
    eta = rated_data[1]
    # print(  individual_torque_average, 
    #         individual_TRV, sep='\n')

    return individual_TRV, individual_FRW, Trip, Em, Ea, eta, Cost, Efficiency, RipplePerformanceSum

def the_script():

    # print(sum(sizes)*total_loss)

    raw = raw_string.split('\n')
    DATA = Get_FRW(raw)
    print('------------Q%dp%d'% (Q, p))
    print(' '.join(['%g'%(el) for el in DATA]))
    figname = './loss_donut_chart_Q%dp%d%s.png' % (Q,p,proj_name)

    #colors # https://www.schemecolor.com/color/green
    # colors = ['#C766A1','#F49762','#FFEC8A','#A1D47B', '#32D081', '#CCEFAB', '#F66867', '#F7DD7D', '#5BC5EA', '#3D93DD']
    colors = ['#C766A1','#F49762','#FFEC8A','#A1D47B']
    labels = ['Iron', 'Magnet', 'Copper', 'Windage']

    #explsion
    explode = tuple([0.025 for i in range(len(labels))]) # (0.025,0.025,0.025,0.025,0.025,0.025)


    import matplotlib as mpl
    mpl.rcParams['font.size'] = 15.0
    mpl.rcParams['font.family'] = ['Times New Roman']
    # font = {'family' : 'Times New Roman', #'serif',
    #     'color' : 'darkblue',
    #     'weight' : 'normal',
    #     'size' : 14,}
    fig = plt.figure()
    ax1 = plt.gca()
    ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.1f%%', startangle=0, pctdistance=0.50, explode = explode)
    #draw circle
    centre_circle = plt.Circle((0,0),0.70,fc='white')
    ax1.add_artist(centre_circle)
    # Equal aspect ratio ensures that pie is drawn as a circle
    ax1.axis('equal')  
    # plt.tight_layout()
    fig.savefig(figname)


    # https://medium.com/@kvnamipara/a-better-visualisation-of-pie-charts-by-matplotlib-935b7667d77f

    # 如果转速可以变化，可以画成Stack Plots
    # https://www.youtube.com/watch?v=xN-Supd4H38
    # ax1.legend(loc=(0.05, 0.07))


r'run#62399/'  # spec_ECCE_PMSM_ (Q6p2)
best_idx=3228
proj_name='proj3280'
Q = 6
p = 2
raw_string = '''-------
proj3280-IDp2s1
0,3280,O1=4.82883,O2=12.2886,f1=187.023,f2=-0.922349,f3=18.8713
0.933699,0.918724,7.92742,0.0753712,13.9166,0.206778,13.2284,15696.7,1373.7,236.81,198.839,37.9711,328.271,0,126.56,0,4e+06,0,264.455,2203.24
50000,0.922349,4209.41,254.089,2757.91,201.71,0,475.431,520.267,0.000861708,185.894,92.5926
29.8532,29,44.6469,1.91563,2.87344,34.9151,35.2504,28.1269,0,0,0,6,5.48432,0.75,2.50013,84.876,84.876,18.3845,17.5279,2.50013,0,2,1
ind3280TranPMSM
Average Torque: 7.92742 Nm
Normalized Torque Ripple: 7.53712 %
Average Force Mag: 13.9166 N
Normalized Force Error Mag: 20.6778%, (+)10.6734% (-)-20.6778%
Maximum Force Error Angle: 13.2284 [deg], (+)7.86956 deg (-)-13.2284 deg
Extra Info:
    Average Force Vector: (13.4082, -3.72717) N
    Torque Ripple (Peak-to-Peak): 0.597499 Nm
    Force Mag Ripple (Peak-to-Peak): 4.36301 N
    basic info:('Case', 1.0)('freq', 1000.0)('0', 11.2488688449)
    jmag loss info: 15696.7, 1373.7, 236.81, 198.839, 37.9711
    femm loss info: 328.271, 0, 126.56, 0, 4e+06, 0
    PF: 0.933699
    eta, windage, total_loss: 0.918724, 264.455, 2203.24
'''
total_loss = 2203.24
sizes = np.array( [198.839+37.9711, 1373.7, 328.271, 264.455]) / total_loss
the_script()

r'run#62499/'  # spec_PEMD_BPMSM_Q12p2
best_idx=4235
proj_name='proj4240'
Q = 12
p = 2
raw_string = '''-------
proj4240-IDp2s1
0,4240,O1=3.34728,O2=11.9919,f1=199.968,f2=-0.940189,f3=18.8043
0.619281,0.938846,10.289,0.118452,34.3622,0.246812,11.4991,51615,454.917,1070.57,896.367,174.199,393.48,0,266.609,0,4e+06,0,186.535,2105.5
50000,0.940189,3180.78,412.404,703.688,126.87,0,1656,281.811,0.000715861,143.227,92.5926
21.088,14,46.0652,0.650198,0.975297,34.2786,30.7038,13.6394,0,0,0,12,5.42858,0.75,4.31275,72.1857,72.1857,26.7356,8.83826,2.64537,0,2,1
ind4240TranPMSM
Average Torque: 10.289 Nm
Normalized Torque Ripple: 11.8452 %
Average Force Mag: 34.3622 N
Normalized Force Error Mag: 24.6812%, (+)24.6812% (-)-22.42%
Maximum Force Error Angle: 11.4991 [deg], (+)8.90082 deg (-)-11.4991 deg
Extra Info:
    Average Force Vector: (-28.6381, 18.9901) N
    Torque Ripple (Peak-to-Peak): 1.21875 Nm
    Force Mag Ripple (Peak-to-Peak): 16.185 N
    basic info:('Case', 1.0)('freq', 1000.0)('0', 7.72796064466)
    jmag loss info: 51615, 454.917, 1070.57, 896.367, 174.199
    femm loss info: 393.48, 0, 266.609, 0, 4e+06, 0
    PF: 0.619281
    eta, windage, total_loss: 0.938846, 186.535, 2105.5
'''
total_loss = 2105.5
sizes = np.array( [896.367+174.199, 454.917, 393.48, 186.535]) / total_loss
the_script()


r'run#62599/'  # spec_PEMD_BPMSM_Q6p1)
best_idx=3538
proj_name='proj3543'
Q = 6
p = 1
raw_string = '''-------
proj3543-IDp1s1
0,3543,O1=2.85563,O2=9.52411,f1=189.139,f2=-0.944572,f3=12.6163
0.9658,0.945917,25.4134,0.179109,63.7517,0.126277,6.50855,37407.2,2432.1,901.458,659.753,241.705,520.904,0,405.183,0,4e+06,0,710.329,4564.79
50000,0.944572,2934.05,253.751,1523.14,115.721,0,564.551,476.889,0.000554417,57.9875,92.5926
53.6424,29,61.771,4.63907,6.95861,26.0559,24.4846,36.6188,0,0,0,6,5.85443,0.75,6.89133,122.818,122.818,20.9926,27.2826,3.86948,0,1,1
ind3543TranPMSM
Average Torque: 25.4134 Nm
Normalized Torque Ripple: 17.9109 %
Average Force Mag: 63.7517 N
Normalized Force Error Mag: 12.6277%, (+)12.4033% (-)-12.6277%
Maximum Force Error Angle: 6.50855 [deg], (+)6.50855 deg (-)-6.14227 deg
Extra Info:
    Average Force Vector: (-63.5519, 5.04301) N
    Torque Ripple (Peak-to-Peak): 4.55175 Nm
    Force Mag Ripple (Peak-to-Peak): 15.9577 N
    basic info:('Case', 1.0)('freq', 500.0)('0', -81.2077860673)
    jmag loss info: 37407.2, 2432.1, 901.458, 659.753, 241.705
    femm loss info: 520.904, 0, 405.183, 0, 4e+06, 0
    PF: 0.9658
    eta, windage, total_loss: 0.945917, 710.329, 4564.79
'''
total_loss = 4564.79
sizes = np.array( [659.753+241.705, 2432.1, 520.904, 710.329]) / total_loss
the_script()

r'run#62699/'  # spec_PEMD_BPMSM_Q12p4)
best_idx=2623
proj_name='proj2761'
Q = 12
p = 4
raw_string = '''-------
proj2761-IDp4s1
0,2761,O1=1.67271,O2=2.56622,f1=111.871,f2=-0.956511,f3=2.05655
0.877466,0.956982,19.4827,0.0173638,35.5549,0.0189782,1.32971,83124.1,741.741,1466.4,1333.78,132.618,271.946,0,141.623,0,4e+06,0,271.25,2751.34
50000,0.956511,2273.34,115.692,605.933,130.323,0,1197.91,223.482,0.000346803,75.6394,92.5926
10.5687,14,44.9525,4.06447,6.09671,41.6358,29.1343,19.4131,0,0,0,12,6,0.75,2.81223,44.9999,44.9999,3.69991,31.6904,2.81223,0,4,1
ind2761TranPMSM
Average Torque: 19.4827 Nm
Normalized Torque Ripple: 1.73638 %
Average Force Mag: 35.5549 N
Normalized Force Error Mag: 1.89782%, (+)1.71668% (-)-1.89782%
Maximum Force Error Angle: 1.32971 [deg], (+)1.11393 deg (-)-1.32971 deg
Extra Info:
    Average Force Vector: (26.6952, 23.4844) N
    Torque Ripple (Peak-to-Peak): 0.338293 Nm
    Force Mag Ripple (Peak-to-Peak): 1.28513 N
    basic info:('Case', 1.0)('freq', 2000.0)('0', 27.2342524488)
    jmag loss info: 83124.1, 741.741, 1466.4, 1333.78, 132.618
    femm loss info: 271.946, 0, 141.623, 0, 4e+06, 0
    PF: 0.877466
    eta, windage, total_loss: 0.956982, 271.25, 2751.34
'''
total_loss = 2751.34
sizes = np.array( [1333.78+132.618, 741.741, 271.946, 271.25]) / total_loss
the_script()

r'run#62799/'  # spec_PEMD_BPMSM_Q24p1
best_idx=17
proj_name='proj18'
Q = 24
p = 1
raw_string = '''-------
proj18-IDp1s1
0,18,O1=2.77849,O2=6.2196,f1=163.306,f2=-0.972844,f3=5.62524
0.952466,0.972669,14.5441,0.123586,30.2746,0.0450643,2.25223,187.53,61.755,669.859,499.238,170.621,317.567,0,235.191,0,4e+06,0,234.685,1283.87
50000,0.972844,1395.7,257.369,67.5781,82.3761,0,733.024,255.349,0.000569152,101.324,92.5926
7.21255,6,46.0347,1.70947,2.5642,26.8597,14.04,8.0177,0,0,0,24,3,0.75,6.66985,123.744,123.744,22.8053,12.8096,5.86961,0,1,1
ind18TranPMSM
Average Torque: 14.5441 Nm
Normalized Torque Ripple: 12.3586 %
Average Force Mag: 30.2746 N
Normalized Force Error Mag: 4.50643%, (+)3.66859% (-)-4.50643%
Maximum Force Error Angle: 2.25223 [deg], (+)2.09437 deg (-)-2.25223 deg
Extra Info:
    Average Force Vector: (-30.2712, -0.45051) N
    Torque Ripple (Peak-to-Peak): 1.79744 Nm
    Force Mag Ripple (Peak-to-Peak): 2.47495 N
    basic info:('Case', 1.0)('freq', 500.0)('0', -35.5308485767)
    jmag loss info: 187.53, 61.755, 669.859, 499.238, 170.621
    femm loss info: 317.567, 0, 235.191, 0, 4e+06, 0
    PF: 0.952466
    eta, windage, total_loss: 0.972669, 234.685, 1283.87
'''
total_loss = 1283.87
sizes = np.array( [669.859, 61.755, 317.567, 234.685]) / total_loss
the_script()


r'run#62899/'  # spec_PEMD_BPMSM_Q12p1
best_idx=3680
proj_name='proj3681'
Q = 12
p = 1
raw_string = '''-------
proj3681-IDp1s1
0,3681,O1=3.16206,O2=9.40385,f1=181.832,f2=-0.981832,f3=13.1549
0.999,0.982089,18.0316,0.112617,24.7578,0.156405,7.77445,31055.8,85.9826,279.705,208.901,70.8042,435.51,0,337.099,0,4e+06,0,231.907,1033.1
50000,0.981832,925.217,297.539,75.8922,98.411,0,246.881,206.494,0.0004563,81.7265,92.5926
20.4055,14,48.9056,4.97929,7.46894,27.1287,55.2312,17.0857,0,0,0,12,5.99869,0.75,6.75944,117.6,117.6,20.4682,14.9293,2.87997,0,1,1
ind3681TranPMSM
Average Torque: 18.0316 Nm
Normalized Torque Ripple: 11.2617 %
Average Force Mag: 24.7578 N
Normalized Force Error Mag: 15.6405%, (+)15.6405% (-)-9.49383%
Maximum Force Error Angle: 7.77445 [deg], (+)7.66744 deg (-)-7.77445 deg
Extra Info:
    Average Force Vector: (-23.4532, -7.93071) N
    Torque Ripple (Peak-to-Peak): 2.03065 Nm
    Force Mag Ripple (Peak-to-Peak): 6.22272 N
    basic info:('Case', 1.0)('freq', 500.0)('0', -23.7178610371)
    jmag loss info: 31055.8, 85.9826, 279.705, 208.901, 70.8042
    femm loss info: 435.51, 0, 337.099, 0, 4e+06, 0
    PF: 0.999
    eta, windage, total_loss: 0.982089, 231.907, 1033.1
'''
total_loss = 1033.1
sizes = np.array( [279.705, 85.9826, 435.51,  231.907]) / total_loss
the_script()




plt.show()


