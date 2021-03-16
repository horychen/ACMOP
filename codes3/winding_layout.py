    #禁止在cache时打印
def print(*arg, **kwarg):
    pass

def infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(layer_X_phases, coil_pitch):
    return layer_X_phases[-coil_pitch:] + layer_X_phases[:-coil_pitch]
def infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(layer_X_signs, coil_pitch):
    temp = layer_X_signs[-coil_pitch:] + layer_X_signs[:-coil_pitch]
    return [('-' if el == '+' else '+') for el in temp]
def infer_Y_layer_grpAC_from_X_layer_and_coil_pitch_y(grouping_AC, coil_pitch):
    return grouping_AC[-coil_pitch:] + grouping_AC[:-coil_pitch]
    # >>> a = [1,2,3,4,5,6,7,8,0]
    # >>> a[-7:] + a[:-7]
    # [3, 4, 5, 6, 7, 8, 0, 1, 2]

class winding_layout_v2(object):
    def __init__(self, DPNV_or_SEPA, Qs, p, ps=None, coil_pitch_y=None, pr=None):
        ''' Naming convention:
        # right layer = 1st layer = X layer = torque layer for separate winding
        # left layer  = 2nd layer = Y layer = suspension layer for separate winding
        '''

        # number of phase
        m = 3
        if p % m == 0 or (ps is not None and ps % m == 0):
            if DPNV_or_SEPA is not None:
                print('Warning: asymmetric suspension winding for DPNV.')

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Regular Motor Winding
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        if DPNV_or_SEPA is None \
        and Qs == 12 \
        and p == 2 \
        and coil_pitch_y == 3:

            self.layer_X_phases = ['U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V']
            self.layer_X_signs  = ['+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.number_parallel_branch = 1
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = True

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Separate Winding for Bearingless Motor
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        if DPNV_or_SEPA == False \
        and Qs == 24 \
        and p == 2 \
        and ps == 1:
            self.layer_X_phases = [ 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V' ]*2
            self.layer_X_signs  = [ '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-' ]*2
            self.layer_Y_phases = [ 'U', 'U', 'V', 'V', 'V', 'V', 'W', 'W', 'W', 'W', 'U', 'U', 'U', 'U', 'V', 'V', 'V', 'V', 'W', 'W', 'W', 'W', 'U', 'U' ]
            self.layer_Y_signs  = [ '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-' ]
            if coil_pitch_y is None:
                self.coil_pitch_y = int(Qs/(2*p)) # = Qs / poles, for single layer
            else:
                raise Exception('Do not specify coil_pitch_y for single layer winding.')

            self.number_parallel_branch = 1.
            self.number_winding_layer = 1 # for torque winding

            # Excitation options
            self.bool_3PhaseCurrentSource = True
            # Commutating sequence is predefined in codes (population.py) for separate winding, which is correct when rotor rotates CCW. This is still, however, used in FEMM_Solver.py, which should also be avoided.
            self.CommutatingSequenceD = 0 # D stands for Drive winding (i.e., torque winding)
            self.CommutatingSequenceB = 0 # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field

        if DPNV_or_SEPA == False \
        and Qs == 36 \
        and p == 3 \
        and ps == 2:
            self.layer_X_phases = [ 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V' ]*3
            self.layer_X_signs  = [ '+', '+', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+' ]*3
            self.layer_Y_phases = [ 'U', 'U', 'U', 'W', 'W', 'W', 'V', 'V', 'V', 'U', 'U', 'U', 'W', 'W', 'W', 'V', 'V', 'V' ]*2
            self.layer_Y_signs  = [ '+', '+', '+', '-', '-', '-', '+', '+', '+', '-', '-', '-', '+', '+', '+', '-', '-', '-' ]*2
            if coil_pitch_y is None:
                self.coil_pitch_y = int(Qs/(2*p)) # = Qs / poles, for single layer
            else:
                raise Exception('Do not specify coil_pitch_y for single layer winding.')

            self.number_parallel_branch = 1.
            self.number_winding_layer = 1 # for torque winding

            # Excitation options
            self.bool_3PhaseCurrentSource = True
            # Commutating sequence is predefined in codes (population.py) for separate winding, which is correct when rotor rotates CCW. This is still, however, used in FEMM_Solver.py, which should also be avoided.
            self.CommutatingSequenceD = 0 # D stands for Drive winding (i.e., torque winding)
            self.CommutatingSequenceB = 0 # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Combined Winding for Bearingless Motor
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

        if DPNV_or_SEPA == True \
        and Qs == 36 \
        and p == 3 \
        and ps == 2 \
        and coil_pitch_y == 5:

            self.layer_X_phases = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        # This winding is wrong
        if DPNV_or_SEPA == True \
        and Qs == 27 \
        and p == 3 \
        and ps == 2 \
        and coil_pitch_y == 4:

            self.layer_X_phases = ['U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V', 'U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V', 'U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if DPNV_or_SEPA == True \
        and Qs == 24 \
        and p == 2 \
        and ps == 1:
            self.layer_X_phases = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V'] # ExampleQ24p2m3ps1: torque winding outer layer
            self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.coil_pitch_y = int(Qs/(2*p)) if coil_pitch_y is None else coil_pitch_y # Y layer can be inferred from coil pitch and X layer diagram
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            # grouping AC is valid for the X layer signs, the Y layer can always be inferred from the X layer signs and the coil pitch.
            self.grouping_AC    = [  0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0] # 只取决于1s tlayer/outer layer/right layer的反相情况
            self.number_parallel_branch = 2.
            self.number_winding_layer = 2 # for torque winding and this means there could be a short pitch

            # Excitation options
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.CommutatingSequenceD = 1 # D stands for Drive winding (i.e., torque winding)
            self.CommutatingSequenceB = 0 # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field

        if DPNV_or_SEPA == True \
        and Qs == 24 \
        and p == 2 \
        and ps == 3 \
        and coil_pitch_y == 5:

            self.layer_X_phases = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if DPNV_or_SEPA == True \
        and Qs == 24 \
        and p == 1 \
        and ps == 2 \
        and coil_pitch_y == 9:

            self.layer_X_phases = ['U', 'U', 'U', 'W', 'W', 'W', 'W', 'V', 'V', 'V', 'V', 'U', 'U', 'U', 'U', 'W', 'W', 'W', 'W', 'V', 'V', 'V', 'V', 'U']
            self.layer_X_signs  = ['+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        if DPNV_or_SEPA == True \
        and Qs == 6 \
        and p == 1 \
        and ps == 2 \
        and coil_pitch_y == 2:

            self.layer_X_phases = ['U', 'W', 'V', 'U', 'W', 'V']
            self.layer_X_signs  = ['+', '-', '+', '-', '+', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 1, 0, 1, 0, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0


        if DPNV_or_SEPA == True \
        and Qs == 6 \
        and p == 2 \
        and ps == 1 \
        and coil_pitch_y == 1:

            self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W']
            self.layer_X_signs  = ['+', '+', '+', '+', '+', '+']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 1, 0, 1, 0, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0


        if DPNV_or_SEPA == True \
        and Qs == 12 \
        and p == 1 \
        and ps == 2 \
        and coil_pitch_y == 5:

            self.layer_X_phases = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0


        if DPNV_or_SEPA == True \
        and Qs == 12 \
        and p == 2 \
        and ps == 1 \
        and coil_pitch_y == 3:

            self.layer_X_phases = ['U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V']
            self.layer_X_signs  = ['+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0


        if DPNV_or_SEPA == True \
        and Qs == 12 \
        and p == 4 \
        and ps == 5 \
        and coil_pitch_y == 1:

            self.layer_X_phases = ['U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W']
            self.layer_X_signs  = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0


        if DPNV_or_SEPA == True \
        and Qs == 18 \
        and p == 2 \
        and ps == 1 \
        and coil_pitch_y == 4:

            self.layer_X_phases = ['U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V', 'U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        # TIA-ISMB 2020 stator winding for the double layer rotor winding
        if DPNV_or_SEPA == True \
        and Qs == 18 \
        and p == 2 \
        and ps == 3 \
        and coil_pitch_y == 4:

            self.layer_X_phases = ['U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V', 'U', 'W', 'W', 'V', 'U', 'U', 'W', 'V', 'V']
            self.layer_X_signs  = ['+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '-']
            self.coil_pitch_y   = coil_pitch_y
            self.layer_Y_phases = infer_Y_layer_phases_from_X_layer_and_coil_pitch_y(self.layer_X_phases, self.coil_pitch_y)
            self.layer_Y_signs  = infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(self.layer_X_signs, self.coil_pitch_y)

            self.grouping_AC            = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1]
            self.number_parallel_branch = 2
            self.number_winding_layer   = 2

            self.bool_3PhaseCurrentSource = False
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # End of winding definition
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        try: 
            self.coil_pitch_y
            self.distributed_or_concentrated = False if abs(self.coil_pitch_y) == 1 else True

            # below is valid for PMSM only

            q = SPP = Qs / (2*p * m)
            # print('[wily] q =', q)
            # if q%1 > 0:
            #     print('[wily] Fractional slot winding with an SPP of %g.'%(q))
            #     if q>1:
            #         raise Exception('This case is not thought thorough.')

            # Example: Q24p1
            # + + + + 0 0 0 0 0 - - - -
            # 0 1 2 3 4 5 6 7 8 9
            # 电流激励是sin(t)，0时刻的时候，电流导致的磁场刚好和U相相轴正交，所以转子的初始位置就应该是U轴相轴的位置。

            # Valid for PMSM                                                这个1是由于1号槽从x轴便宜了半个槽距角所导致的，那么为什么1号槽这么特别？因为我们默认了u相是从1号槽开始的，但是有时候U相可能从24号槽开始，那么就要进一步进入这种情况
            # self.initial_excitation_bias_compensation_deg = 360/Qs*0.5 * (1 + self.coil_pitch_y)

            # 寻找初始d轴激励角度
            if q%1 == 0:
                # integral slot

                phase_U_belt = self.layer_X_phases[:int(q)]
                number_of_U = sum([1 for el in phase_U_belt if el =='U'])
                if number_of_U<q:
                    # print('目前只有%d个字母U，需要寻找一共q(=%d)个字母U。'%(number_of_U, q))
                    for phase_U_starting_slot_number in range(-1, -int(q)-1, -1):
                        new_phase_U_belt = self.layer_X_phases[phase_U_starting_slot_number:] + phase_U_belt
                        number_of_U = sum([1 for el in new_phase_U_belt if el =='U'])
                        # print(phase_U_starting_slot_number, number_of_U)
                        if number_of_U < q:
                            continue
                        else:
                            break
                else:
                    phase_U_starting_slot_number = 1 # 没有等于0的哈，等于0你得以槽的中线作图绘制定子铁芯；现在的代码是以齿的中线作图绘制定子铁芯的哦。

                self.deg_winding_U_phase_phase_axis_angle = 360/Qs*0.5 * (phase_U_starting_slot_number + self.coil_pitch_y+(SPP-1) )
            else:
                # This clause includes fractional slot winding with an SPP value below 1.
                phase_U_starting_slot_number = 1

                # print('[wily] Fractional slot winding with an SPP of %g.'%(q))

                if q<1:
                    # fractional slot and q<1
                    self.deg_winding_U_phase_phase_axis_angle = 360/Qs*0.5 * (phase_U_starting_slot_number + self.coil_pitch_y )
                else:
                    # fractional slot and q>1
                    self.deg_winding_U_phase_phase_axis_angle = 360/Qs*0.5 * (phase_U_starting_slot_number + self.coil_pitch_y )
                    msg = '[wily][Warning] This case (q=%g) is not thought thorough, so you must inspect the initial excitation angle and initial rotor position manually to make sure it is id=0 control.'%(q)
                    print(msg)
                    if q>2:
                        raise Exception(msg)

                    # 2019/12/25: Q18, p2, ps1 case is verified to be okay with q=1.5.

            print('[wily] phase_U_starting_slot_number', phase_U_starting_slot_number)
            print('[wily] self.deg_winding_U_phase_phase_axis_angle:', self.deg_winding_U_phase_phase_axis_angle)
            # quit()
        except:
            print(Qs,p,ps,coil_pitch_y)
            raise Exception('Error: This winding is not implemented.')


        # 这是实际在pre_procee中调用的字典
        self.dict_coil_connection = {'layer X phases': self.layer_X_phases, 'layer X signs':self.layer_X_signs,   # 这里的命名规则是按照seprate winding的情况来的。
                                     'layer Y phases': self.layer_Y_phases, 'layer Y signs':self.layer_Y_signs}   # 这里的命名规则是按照seprate winding的情况来的。


        # ACMDM: motor mode and suspension mode (ACMDM is the initiative my first try to clean up the code, which has been aborted)
        self.DPNV_or_SEPA = DPNV_or_SEPA
        self.Qs = Qs
        self.p = p
        self.m = 3
        self.SPP = self.Qs / (2*self.p * self.m)
        self.ps = ps
        self.pr = pr
        if DPNV_or_SEPA is not None:
            if self.number_winding_layer == 1: # this is equivalent to separate winding for now
                self.list_layer_motor_phases      = [self.layer_X_phases]
                self.list_layer_motor_signs       = [self.layer_X_signs]
                self.list_layer_suspension_phases = [self.layer_Y_phases]
                self.list_layer_suspension_signs  = [self.layer_Y_signs]

            elif self.number_winding_layer == 2: # this is equivalent to combined winding for now
                self.list_layer_motor_phases      = [self.layer_X_phases, self.layer_Y_phases] #??? 这有啥用？？？
                self.list_layer_motor_signs       = [self.layer_X_signs, self.layer_Y_signs]
                self.list_layer_suspension_phases = self.list_layer_motor_phases[::]
                suspension_layer_X_signs = [(el if ac==0 else ('-' if el == '+' else '+')) for el, ac in zip(self.layer_X_signs, self.grouping_AC)]
                self.list_layer_suspension_signs  = [ suspension_layer_X_signs, 
                                                      infer_Y_layer_signs_from_X_layer_and_coil_pitch_y(suspension_layer_X_signs, self.coil_pitch_y)]
            else:
                raise Exception('Invalid number_winding_layer.')

            # save for winding function analysis
            self.ox_distribution_three_phase = []
            for phases, signs in zip(self.list_layer_motor_phases, self.list_layer_motor_signs):
                temp = [phase+sign for phase, sign in zip(phases, signs)]
                self.ox_distribution_three_phase = temp if self.ox_distribution_three_phase == [] else [a+b for a,b in zip(self.ox_distribution_three_phase, temp)]
            def replace_uvw_with_ox(string, UVW):
                string = string.replace(UVW+'+', 'x')
                string = string.replace(UVW+'-', 'o')
                return string
            def replace_uvw_with_empty_string(string, UVW):
                string = string.replace(UVW+'+', '')
                string = string.replace(UVW+'-', '')
                return string
            self.ox_distribution_phase_U = [replace_uvw_with_ox(el, 'U') for el in self.ox_distribution_three_phase]
            self.ox_distribution_phase_U = [replace_uvw_with_empty_string(el, 'V') for el in self.ox_distribution_phase_U]
            self.ox_distribution_phase_U = [replace_uvw_with_empty_string(el, 'W') for el in self.ox_distribution_phase_U]
            # print(self.ox_distribution_three_phase)
            # print(self.ox_distribution_phase_U)

# if __name__ == '__main__':
#     wily = winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=2, ps=1, coil_pitch=6)
#     quit()

class pole_specific_winding_with_neutral():
    def __init__(self, Qr, p, ps, coil_pitch_y=None):

        # Since p=2, this is a double layer windidng reduced to single layer implementation with neutral plate
        if Qr == 30 \
        and p == 2 \
        and ps == 3 \
        and coil_pitch_y == 10:
            self.pairs = [(1, 11, 21), (2, 12, 22), (3, 13, 23), (4, 14, 24), (5, 15, 25), (6, 16, 26), (7, 17, 27), (8, 18, 28), (9, 19, 29), (10, 20, 30)]

        # Since p=2, this is a double layer windidng reduced to single layer implementation with neutral plate
        if Qr == 18 \
        and p == 2 \
        and ps == 3 \
        and coil_pitch_y == 6:
            self.pairs = [(1, 7, 13), (2, 8, 14), (3, 9, 15), (4, 10, 16), (5, 11, 17), (6, 12, 18)]

        # this is for test purpose
        # to test the Qr30ps3 rotor circuit and see if it induces current under field without one pole pair component.
        # I only swap p and ps, and leave others the same with the rotor winding above
        if Qr == 30 \
        and p == 3 \
        and ps == 2 \
        and coil_pitch_y == 10:
            self.pairs = [(1, 11, 21), (2, 12, 22), (3, 13, 23), (4, 14, 24), (5, 15, 25), (6, 16, 26), (7, 17, 27), (8, 18, 28), (9, 19, 29), (10, 20, 30)]

        # this is for test purpose
        # to test the Qr30ps3 rotor circuit and see if it induces current under field WITH one pole pair component.
        # I only swap p and ps, and leave others the same with the rotor winding above
        if Qr == 30 \
        and p == 1 \
        and ps == 2 \
        and coil_pitch_y == 10:
            self.pairs = [(1, 11, 21), (2, 12, 22), (3, 13, 23), (4, 14, 24), (5, 15, 25), (6, 16, 26), (7, 17, 27), (8, 18, 28), (9, 19, 29), (10, 20, 30)]


        if Qr == 20 \
        and p == 3 \
        and ps == 2 \
        and coil_pitch_y == 10:
            self.pairs = [(1, 11), (2, 12), (3, 13), (4, 14), (5, 15), (6, 16), (7, 17), (8, 18), (9, 19), (10, 20)]

        if Qr == 24 \
        and p == 3 \
        and ps == 2 \
        and coil_pitch_y == 12:
            self.pairs = [(1, 13), (2, 14), (3, 15), (4, 16), (5, 17), (6, 18), (7, 19), (8, 20), (9, 21), (10, 22), (11, 23), (12, 24)]

        if Qr == 28 \
        and p == 3 \
        and ps == 2 \
        and coil_pitch_y == 14:
            '''
            For Qs=27, find valid set of pole specific winding with neutral plate.
            ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
            c/d = 1/2
            Qr  = 28
            p   = 3
            ps  = 2
            k1  = 1
            k   = 14
            t   = 2
            nl  = 15
            y   =  14
            m/n =  14/3 (this is consistent with Pyrhonen@(2.83): m cannot be multiples of n.)
            FSW
            q=z/n=1/3 = 28/(2*3*14)
            p*gamma =180
            kw_h = sin(h*gamma/2) =1.00, 0.00, -1.00, -0.00, 1.00, 0.00, -1.00, -0.00, 1.00, ...
            number_of_coils = 14.0
            number of coils per phase = 1.0
            ----
            Winding Layout (m=14 phases):
             | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |10 |11 |12 |13 |14 |15 |16 |17 |18 |19 |20 |21 |22 |23 |24 |25 |26 |27 |28
             | a | b | c | d | e | f | g | h | i | j | k | l | m | n | a | b | c | d | e | f | g | h | i | j | k | l | m | n
             [(0, 14), (1, 15), (2, 16), (3, 17), (4, 18), (5, 19), (6, 20), (7, 21), (8, 22), (9, 23), (10, 24), (11, 25), (12, 26), (13, 27)]
            [Finished in 0.1s]
            '''
            # self.pairs=[(0, 14), (1, 15), (2, 16), (3, 17), (4, 18), (5, 19), (6, 20), (7, 21), (8, 22), (9, 23), (10, 24), (11, 25), (12, 26), (13, 27)]
            self.pairs=[(1, 15), (2, 16), (3, 17), (4, 18), (5, 19), (6, 20), (7, 21), (8, 22), (9, 23), (10, 24), (11, 25), (12, 26), (13, 27), (14, 28)]

from pylab import np, plt, fft, linspace
import scipy.integrate as integrate

def nextpow2(L):
    n = 0
    while 2**n < L:
        n += 1
    return n

def periodic2pi(x):
    # make x periodic positive
    while(x<0):
        x += 2*np.pi

    # make it within 2*pi
    if int(x / (2*np.pi)) != 0:
        x = x % (2*np.pi)

    return x

def segmented_func(x, lst_x, lst_y):
    ''' The value in lst_x must increase monotonously.
        根据判断，x在lst_x中的位置，返回对应的阶梯函数值。
    '''
    x = periodic2pi(x)
    
    if len(lst_x) != len(lst_y):
        raise Exception('Error: lst_x & y must has same length.')

    for i in range(0,len(lst_x)):
        if x > lst_x[i]: # x 超过了lst_x[i]的位置，
            if i == len(lst_x)-1: # last i。比最后一个位置lst_x[-1]大，比2pi小。
                return lst_y[-1]

        else: # x 比lst_x[i]小，但是比lst_x[i-1]大。
            if i-1 < 0: # first i。比第一个位置lst_x[0]小，比0大。
                return lst_y[-1] # 注意lst_y[0]的值不是0，但是lst_y[-1]是0.
            else: # 普通情况：
                return lst_y[i-1]

class PhaseWinding:
    '''[PhaseWinding]
    
    [The winding function and turn function are consistent with Lipo's 2012 book.]
    '''
    def __init__(self, Qs, m, turns_per_slot, ox_distribution_phase_U, desc_type = 'STEP_DISTRIBUTED'):
        if desc_type == 'STEP_DISTRIBUTED':
            pass
        else:
            raise Exception('Not implemented for other type of descripsion.')

        self.slot_per_phase = Qs / m # Only upper layer counts for one slot.
        self.turns_per_slot = turns_per_slot
        self.degree_between_slots = 360.0 / Qs # mechanical deg.
        self.radian_between_slots = self.degree_between_slots / 180.0 * np.pi # mechanical rad.
        self.ox_distribution_phase_U = ox_distribution_phase_U
        print(self.ox_distribution_phase_U)

        # turn function 
        self.setTurnFuncObject(self.ox_distribution_phase_U) # define self.turn_func
        # <turn function> defined by Lipo 2012
        self.avg_val_of_turn_func = integrate.quad(self.turn_func, 0, 2*np.pi)[0] / (2*np.pi) # [1] is error of integration
        # winding function
        self.winding_func = lambda x: self.turn_func(x)-self.avg_val_of_turn_func

        # get self.sym_begin_pos 
        self.setSymPos()
        # symmetric turn function
        self.sym_turn_func = lambda x: self.turn_func(x + self.sym_begin_pos)
        # symmetric winding function
        self.sym_winding_func = lambda x: self.winding_func(x + self.sym_begin_pos)

    def setTurnFuncObject(self, ox_distribution_phase_U):
        '''  ['x','x','x','x',
              'n','n','n','n',
              'o','o','o','o',
              'n','n','n','n']
              非n地方才记录，x增，o减。
        '''
        x_val, y_val, lst_x, lst_y = 0, 0, [], []
        for el in ox_distribution_phase_U:
            x_val += self.radian_between_slots
            if el == 'x':
                y_val += self.turns_per_slot
                lst_y.append(y_val)
                lst_x.append(x_val)
            elif el == 'xx':
                y_val += 2* self.turns_per_slot                
                lst_y.append(y_val)
                lst_x.append(x_val)
            elif el == 'o':
                y_val -= self.turns_per_slot                
                lst_y.append(y_val)
                lst_x.append(x_val)
            elif el == 'oo':
                y_val -= 2 * self.turns_per_slot                
                lst_y.append(y_val)
                lst_x.append(x_val)
            else:
                pass
        self.lst_x = lst_x
        self.lst_y = lst_y
        self.turn_func = lambda x: segmented_func(x, self.lst_x, self.lst_y)

    def setSymPos(self, index=2):
        '''
            可选择两种不同的对称位置。用作原始turn_func和wind_func的偏移。举例：
            和lst_x一一对应的lst_y:
                50
                100
                150
                200      <- pos2 @ max y 
                150x
                100x
                50x
                0        <- pos1 @ min y
        '''
        # symmetrical pos_2 = pos_200 + (pos_150x - pos_200) / 2
        index = self.lst_y.index(max(self.lst_y))
        self.sym_begin_pos_2 = self.lst_x[index] + (self.lst_x[index+1] - self.lst_x[index]) / 2.

        # symmetrical pos_1
        self.sym_begin_pos_1 = self.sym_begin_pos_2 + np.pi

        if index == 1:
            self.sym_begin_pos = self.sym_begin_pos_1
        else:
            self.sym_begin_pos = self.sym_begin_pos_2

    def plot2piFft(self, func, Fs, L):
        ''' Fs is the sampling freq. 
            L is length of signal list.
            This plot is for a func that has period of 2pi.

            If you found the time domain wave is not very accurate,
            that is because you set too small Fs, which leads to
            to big step Ts.
        '''
        base_freq = 1.0/(2*np.pi) #频域横坐标除以基频，即以基频为单位，此处的基频为 2*pi rad/s
        Ts = 1.0/Fs
        t = [el*Ts for el in range(0,L)]
        x = [func(el) for el in t]

        # https://www.ritchievink.com/blog/2017/04/23/understanding-the-fourier-transform-by-example/

        # 小明给的代码：
        # sampleF = Fs
        # print('小明：')
        # for f, Y in zip(
        #                 np.arange(0, len(x)*sampleF,1) * 1/len(x) * sampleF, 
        #                 np.log10(np.abs(np.fft.fft(x) / len(x))) 
        #              ):
            # print('\t', f, Y)


        L_4pi = int(4*np.pi / Ts) +1 # 画前两个周期的
        
        self.fig_plot2piFft = plt.figure(7)
        plt.subplot(211)
        plt.plot(t[:L_4pi], x[:L_4pi])
        #title('Signal in Time Domain')
        #xlabel('Time / s')
        #ylabel('x(t)')
        plt.title('Winding Function')
        plt.xlabel('Angular location along air gap [mech. rad.]')
        plt.ylabel('Current Linkage by unit current [Ampere]')

        NFFT = 2**nextpow2(L)
        print('NFFT =', NFFT, '= 2^%g' % (nextpow2(L)), '>= L =', L)
        y = fft(x,NFFT) # y is a COMPLEX defined in numpy
        Y = [2 * el.__abs__() / L for el in y] # /L for spectrum aplitude consistent with actual signal. 2* for single-sided. abs for amplitude.
        f = Fs/2.0/base_freq*linspace(0,1,int(NFFT/2+1)) # unit is base_freq Hz
        #f = Fs/2.0*linspace(0,1,NFFT/2+1) # unit is Hz

        plt.subplot(212)
        plt.plot(f, Y[0:int(NFFT/2+1)])
        plt.title('Single-Sided Amplitude Spectrum of x(t)')
        plt.xlabel('Frequency divided by base_freq [base freq * Hz]')
        #plt.ylabel('|Y(f)|')
        plt.ylabel('Amplitude [1]')
        plt.xlim([0,50])
        # plt.show()

    def plotFuncObj(self, func):
        x = np.arange(0, 2*np.pi, 0.5/180*np.pi)
        y = [func(el) for el in x]
        x = [el/np.pi for el in x] # x for plot. unit is pi
        
        self.fig_plotFuncObj = plt.figure()
        ax = plt.subplot(111) #注意:一般都在ax中设置,不在plot中设置
        ax.plot(x, y)

        xmajorLocator   = plt.MultipleLocator(0.25) #将x主刻度标签设置为20的倍数
        xmajorFormatter = plt.FormatStrFormatter('%.2fπ') #设置x轴标签文本的格式
        ax.xaxis.set_major_locator(xmajorLocator)  
        ax.xaxis.set_major_formatter(xmajorFormatter)

        ##matplotlib.pyplot.minorticks_on()
        ##xminorLocator   = MultipleLocator(0.25)
        ##xminorFormatter = FormatStrFormatter(u'%.2fπ')
        ##ax.xaxis.set_minor_locator(xminorLocator)  
        ##ax.xaxis.set_minor_formatter(xminorFormatter)

        plt.xlabel('Angular location along the gap [mech. rad.]')
        plt.ylabel('Turns of winding [1]')
        # plt.title('Turn Function or Winding Function')
        plt.grid(True) # or ax.grid(True)
        # plt.gcf().savefig("turn_function.png")
        # plt.show()

if __name__ == '__main__':

    # wily = winding_layout_v2(DPNV_or_SEPA=False, Qs=24, p=2, ps=1)
    # wily = winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=2, ps=1, coil_pitch_y=6)
    wily = winding_layout_v2(DPNV_or_SEPA=True, Qs=24, p=1, ps=2, coil_pitch_y=9)

    zQ = 100 # number of conductors/turns per slot
    turns_per_layer = zQ / wily.number_winding_layer
    U_phase = PhaseWinding(wily.Qs, wily.m, turns_per_layer, wily.ox_distribution_phase_U)
    U_phase.plotFuncObj(U_phase.winding_func)
    U_phase.plot2piFft(U_phase.winding_func, Fs=1/(2*np.pi/3600), L=65536*2**4) # 采样频率：在2pi的周期内取720个点


    plt.show()

    if False:
        # PMLSM AAAAA
        # full teeth despite static end effect
        zQ = 50
        turns_per_layer = zQ / wily.number_winding_layer
        ox_distribution_one_phase = ['o','xx','oo','xx',
                                     'o','n','n','n',
                                     'n','n','n','n',
                                     'x','oo','xx','oo',
                                     'x','n','n','n',
                                     'n','n','n','n']
        AA_ph = PhaseWinding(wily.Qs, wily.m, turns_per_layer, ox_distribution_one_phase)
        AA_ph.plot2piFft(AA_ph.winding_func, 10, 10000)







# Obsolte after 20191126
class winding_layout(object):
    def __init__(self, DPNV_or_SEPA, Qs, p, ps=None):

        # separate winding
        if DPNV_or_SEPA == False \
        and Qs == 24 \
        and p == 2:
            self.l41=[ 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V', ]
            self.l42=[ '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', ]
            # separate style for one phase: ---- ++++
            self.l21=[ 'U', 'U', 'V', 'V', 'V', 'V', 'W', 'W', 'W', 'W', 'U', 'U', 'U', 'U', 'V', 'V', 'V', 'V', 'W', 'W', 'W', 'W', 'U', 'U', ]
            self.l22=[ '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', ]
            self.coil_pitch = 6 # = Qs / poles for single layer
            self.CommutatingSequenceD = 0
            self.CommutatingSequenceB = 0
            self.number_parallel_branch = 1.
            self.bool_3PhaseCurrentSource = True
            self.no_winding_layer = 1 # for torque winding

        # combined winding
        if DPNV_or_SEPA == True \
        and Qs == 24 \
        and p == 2:
            # DPNV winding implemented as DPNV winding (GroupAC means it experiences flip phasor excitation from suspension inverter, while GroupBD does not.)
            #                     U-GrBD                        U-GrBD    W-GrBD                        W-GrBD    V-GrBD                        V-GrBD
            #                               W-GrAC    V-GrAC                        V-GrAC    U-GrAC                        U-GrAC    W-GrAC             : flip phases 19-14??? slot of phase U??? (这个例子的这句话看不懂)
            self.l_rightlayer1 = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V'] # ExampleQ24p2m3ps1: torque winding outer layer
            self.l_rightlayer2 = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.l_leftlayer1  = self.l_rightlayer1[::] # ExampleQ24p2m3ps1: torque winding inner layer
            self.l_leftlayer2  = self.l_rightlayer2[::]
            self.grouping_AC   = [  0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0] # 只取决于outerlayer/rightlayer的反相情况
            self.coil_pitch    = 6 # left layer can be inferred from coil pitch and right layer diagram
            self.CommutatingSequenceD = 1 # D stands for Drive winding (i.e., torque winding)
            self.CommutatingSequenceB = 0 # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field
            self.number_parallel_branch = 2.
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.no_winding_layer = 2 # for torque winding and this means there could be a short pitch

            # backward compatibility
            self.l41 = self.l_rightlayer1
            self.l42 = self.l_rightlayer2
            self.l21 = self.l_leftlayer1
            self.l22 = self.l_leftlayer2

        # combined winding
        if DPNV_or_SEPA == True \
        and Qs == 24 \
        and p == 1:
            # DPNV winding implemented as DPNV winding (GroupAC means it experiences flip phasor excitation from suspension inverter, while GroupBD does not.)
            #                     U-GroupBD                               V-GroupBD                               W-GroupBD
            #                                         W-GroupAC                               U-GroupAC                               V-GroupAC           : flip phases 13-16 slot of phase U
            self.l_rightlayer1 = ['U', 'U', 'U', 'U', 'W', 'W', 'W', 'W', 'V', 'V', 'V', 'V', 'U', 'U', 'U', 'U', 'W', 'W', 'W', 'W', 'V', 'V', 'V', 'V'] # ExampleQ24p1m3ps2: torque winding outer layer
            self.l_rightlayer2 = ['+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-'] 
            self.l_leftlayer1  = ['U', 'W', 'W', 'W', 'W', 'V', 'V', 'V', 'V', 'U', 'U', 'U', 'U', 'W', 'W', 'W', 'W', 'V', 'V', 'V', 'V', 'U', 'U', 'U'] # ExampleQ24p1m3ps2: torque winding inner layer
            self.l_leftlayer2  = ['+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+']
            self.grouping_AC   = [  0,   0,   0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0,   0,   0,   1,   1,   1,   1] # 只取决于rightlayer的反相情况
            self.coil_pitch    = 9 # left layer can be inferred from coil pitch and right layer diagram
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0
            self.number_parallel_branch = 2.
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.no_winding_layer = 2 # for torque winding and this means there could be a short pitch
            self.initial_excitation_bias_compensation_deg = 360/24*0.5 # for torque winding # Note that the initial excitation direction is biased (not aligned with x-axis) due to the fact that the u-phase winding is not aligned with x-axis

            # backward compatibility
            self.l41 = self.l_rightlayer1
            self.l42 = self.l_rightlayer2
            self.l21 = self.l_leftlayer1
            self.l22 = self.l_leftlayer2


        # PMSM ONLY

        # combined winding
        # concentrated winding
        if DPNV_or_SEPA == True \
        and Qs == 6 \
        and p == 2:
            # DPNV winding implemented as DPNV winding (GroupAC means it experiences flip phasor excitation from suspension inverter, while GroupBD does not.)
            self.l_rightlayer1 = ['U', 'V', 'W', 'U', 'V', 'W'] # torque winding right layer
            self.l_rightlayer2 = ['+', '+', '+', '+', '+', '+']
            self.l_leftlayer1  = ['W', 'U', 'V', 'W', 'U', 'V']
            self.l_leftlayer2  = ['-', '-', '-', '-', '-', '-']
            # self.grouping_AC   = [  0,   0,   0,   1,   1,   1] # 只取决于outerlayer/rightlayer的反相情况，AC是在悬浮逆变器激励下会反相的
            self.CommutatingSequenceB = 1 # [???]

            self.grouping_AC   = [  0,   1,   0,   1,   0,   1] # Jingwei's layout
            self.CommutatingSequenceB = 0 # [CHECKED]

            self.coil_pitch    = -1 # left layer can be inferred from coil pitch and right layer diagram
                                    # We use negative coil_pitch to indicate concentrated winding
            self.CommutatingSequenceD = 1
            self.number_parallel_branch = 2.
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.no_winding_layer = 2 # for torque winding and this means there could be a short pitch

            self.initial_excitation_bias_compensation_deg = 0 # for torque winding

            # backward compatibility
            self.l41 = self.l_rightlayer1
            self.l42 = self.l_rightlayer2
            self.l21 = self.l_leftlayer1
            self.l22 = self.l_leftlayer2

        # combined winding
        # distrubuted winding
        if DPNV_or_SEPA == True \
        and Qs == 6 \
        and p == 1:
            # DPNV winding implemented as DPNV winding (GroupAC means it experiences flip phasor excitation from suspension inverter, while GroupBD does not.)
            self.l_rightlayer1 = ['U', 'W', 'V', 'U', 'W', 'V'] # torque winding right layer
            self.l_rightlayer2 = ['+', '-', '+', '-', '+', '-']
            self.l_leftlayer1  = ['W', 'V', 'U', 'W', 'V', 'U']
            self.l_leftlayer2  = ['-', '+', '-', '+', '-', '+']
            self.grouping_AC   = [  0,   1,   0,   1,   0,   1] # 只取决于outerlayer/rightlayer的反相情况，AC是在悬浮逆变器激励下会反相的
                                                                # Same with Jingwei's layout
            self.coil_pitch    = 2 # left layer can be inferred from coil pitch and right layer diagram
            self.CommutatingSequenceD = 1
            self.CommutatingSequenceB = 0 # [CHECKED]
            self.number_parallel_branch = 2.
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.no_winding_layer = 2 # for torque winding and this means there could be a short pitch

            self.initial_excitation_bias_compensation_deg = 0 # for torque winding

            # backward compatibility
            self.l41 = self.l_rightlayer1
            self.l42 = self.l_rightlayer2
            self.l21 = self.l_leftlayer1
            self.l22 = self.l_leftlayer2



        # combined winding
        # distrubuted winding
        if DPNV_or_SEPA == True \
        and Qs == 12 \
        and p == 2:
            # DPNV winding implemented as DPNV winding (GroupAC means it experiences flip phasor excitation from suspension inverter, while GroupBD does not.)
            self.l_rightlayer1 = ['U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V'] # torque winding right layer
            self.l_rightlayer2 = ['+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-']
            self.l_leftlayer1  = ['U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V']
            self.l_leftlayer2  = ['+', '-', '+', '-', '+', '-', '+', '-', '+', '-', '+', '-']

            # self.grouping_AC   = [  0,   1,   0,   0,   0,   0,   1,   0,   1,   1,   1,   1] # 只取决于outerlayer/rightlayer的反相情况，AC是在悬浮逆变器激励下会反相的
            # self.CommutatingSequenceB = 1 # 0 # [CHECKED] # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field
            self.grouping_AC   = [  0,   1,   1,   0,   0,   1,   1,   0,   0,   1,   1,   0] # Jingwei's layout
            self.CommutatingSequenceB = 0 # 0 # [CHECKED] # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field
            
            self.CommutatingSequenceD = 1 # D stands for Drive winding (i.e., torque winding)
            self.coil_pitch    = 3 # left layer can be inferred from coil pitch and right layer diagram
            self.number_parallel_branch = 2.
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.no_winding_layer = 2 # for torque winding and this means there could be a short pitch

            self.initial_excitation_bias_compensation_deg = 0 # for u phase torque winding

            # backward compatibility
            self.l41 = self.l_rightlayer1
            self.l42 = self.l_rightlayer2
            self.l21 = self.l_leftlayer1
            self.l22 = self.l_leftlayer2

        # combined winding
        # concentrated winding
        ps = 5 # ps should be specified in future revision
        if DPNV_or_SEPA == True \
        and Qs == 12 \
        and p == 4 \
        and ps == 5:
            # DPNV winding implemented as DPNV winding (GroupAC means it experiences flip phasor excitation from suspension inverter, while GroupBD does not.)
            # self.l_rightlayer1 = ['U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W'] # torque winding right layer
            # self.l_rightlayer2 = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+'] # This configuration gives negative torque
            # self.l_leftlayer1  = ['W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V']
            # self.l_leftlayer2  = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
            # self.l_rightlayer1 = ['U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V'] # Want to change the commutating sequence so the torque is positive? No, it is not that intuitive.
            # self.l_rightlayer2 = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
            # self.l_leftlayer1  = ['V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W', 'V', 'U', 'W']
            # self.l_leftlayer2  = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
            self.l_rightlayer1 = ['U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W'] # Want to change torque sign? Yes, this works
            self.l_rightlayer2 = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
            self.l_leftlayer1  = ['W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V', 'W', 'U', 'V']
            self.l_leftlayer2  = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
            # self.grouping_AC   = [  0,   0,   0,   1,   1,   1,   1,   1,   1,   0,   0,   0] # 只取决于outerlayer/rightlayer的反相情况，AC是在悬浮逆变器激励下会反相的
            self.grouping_AC   = [  0,   0,   1,   1,   0,   0,   1,   1,   0,   0,   1,   1] # Jingwei's layout
            self.coil_pitch    = -1 # left layer can be inferred from coil pitch and right layer diagram
                                    # We use negative coil_pitch to indicate concentrated winding
            self.CommutatingSequenceD = 1 # D stands for Drive winding (i.e., torque winding)
            self.CommutatingSequenceB = 0 # [CHECKED] # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field
            self.number_parallel_branch = 2.
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.no_winding_layer = 2 # for torque winding and this means there could be a short pitch

            self.initial_excitation_bias_compensation_deg = 0 # for torque winding

            # backward compatibility
            self.l41 = self.l_rightlayer1
            self.l42 = self.l_rightlayer2
            self.l21 = self.l_leftlayer1
            self.l22 = self.l_leftlayer2



        # combined winding
        # distrubuted winding
        if DPNV_or_SEPA == True \
        and Qs == 12 \
        and p == 1:
            # DPNV winding implemented as DPNV winding (GroupAC means it experiences flip phasor excitation from suspension inverter, while GroupBD does not.)
            self.l_rightlayer1 = ['U', 'U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V'] # torque winding right layer
            self.l_rightlayer2 = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
            self.l_leftlayer1  = ['U', 'W', 'W', 'V', 'V', 'U', 'U', 'W', 'W', 'V', 'V', 'U']
            self.l_leftlayer2  = ['+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+']

            self.grouping_AC   = [  0,   0,   1,   1,   0,   0,   1,   1,   0,   0,   1,   1]
            self.CommutatingSequenceB = 0 # 0 # [CHECKED] # B stands for Bearing winding (i.e., suspension winding), commutating sequence decides the direction of the rotating field
            
            self.CommutatingSequenceD = 1 # D stands for Drive winding (i.e., torque winding)
            self.coil_pitch    = 5 # left layer can be inferred from coil pitch and right layer diagram
            self.number_parallel_branch = 2.
            self.bool_3PhaseCurrentSource = False # 3PhaseCurrentSource is a macro in circuit setup of JMAG
            self.no_winding_layer = 2 # for torque winding and this means there could be a short pitch

            self.initial_excitation_bias_compensation_deg = 0 # for u phase torque winding

            # backward compatibility
            self.l41 = self.l_rightlayer1
            self.l42 = self.l_rightlayer2
            self.l21 = self.l_leftlayer1
            self.l22 = self.l_leftlayer2


        try: 
            self.coil_pitch
            self.distributed_or_concentrated = False if abs(self.coil_pitch) == 1 else True
        except:
            raise Exception('Error: Not implemented for this winding.')


        # new names
        self.layer_A1 = self.l41
        self.layer_A2 = self.l42
        self.layer_B1 = self.l21
        self.layer_B2 = self.l22
        self.Qs = Qs
        self.p = p


        # # combined winding (special case in which the list grouping_AC is not needed)
        # if DPNV_or_SEPA == True \
        # and Qs == 24 \
        # and p == 2:
        #     # DPNV winding implemented as separate winding
        #     # if self.fea_config_dict['DPNV_separate_winding_implementation'] == True or self.fea_config_dict['DPNV'] == False: 
        #         # You may see this msg because there are more than one designs in the initial_design.txt file.
        #         # msg = 'Not implemented error. In fact, this equivalent implementation works for 4 pole motor only.'
        #         # logging.getLogger(__name__).warn(msg)

        #     # this is legacy codes for easy implementation in FEMM
        #     self.l41=[ 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V', 'W', 'W', 'U', 'U', 'V', 'V']
        #     self.l42=[ '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
        #     # DPNV style for one phase: -- oo ++ oo
        #     self.l21=[  'U', 'U', 'W', 'W', 'V', 'V', 
        #                 'U', 'U', 'W', 'W', 'V', 'V', 
        #                 'U', 'U', 'W', 'W', 'V', 'V', 
        #                 'U', 'U', 'W', 'W', 'V', 'V']
        #     self.l22=[  '-', '-', 'o', 'o', '+', '+', # 横着读和竖着读都是负零正零。 
        #                 'o', 'o', '-', '-', 'o', 'o', 
        #                 '+', '+', 'o', 'o', '-', '-', 
        #                 'o', 'o', '+', '+', 'o', 'o']
        #     self.coil_pitch = 6
        #     self.CommutatingSequenceD = 0
        #     self.CommutatingSequenceB = 0
        #     self.number_parallel_branch = 1.
        #     self.bool_3PhaseCurrentSource = True
        #     self.no_winding_layer = 1 # for torque winding

