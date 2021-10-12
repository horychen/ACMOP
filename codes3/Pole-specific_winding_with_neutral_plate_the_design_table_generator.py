class ValidSet_of_PoleSpecificWindingWithNeutralPlate(object):
    def __init__(self, layers):
        print('_________________________\nlayers =', layers) # =1 or =2

        # ps in N_Even for single layer winding and k in N for double layer winding
        for ps in range(2, MAX_PS, 3-layers):
            print('\t' + 'ps =', ps)

            # p = ps +- 1
            for p in [ps-1, ps+1]:
                print('\t'*2 + f'p = {p:d}', end=' ')
                valid_Qr = self.get_ValidSet_of_Qr(QS, p)
                print('| Valid Qr set:', valid_Qr)

                # k1 in N < ps
                for k1 in range(1, ps):
                    print('\t'*3 + 'k1 =', k1)

                    # c/d = reduced( k1/ps )
                    c = Fraction(k1, ps).numerator
                    d = Fraction(k1, ps).denominator
                    print('\t'*3 + 'c/d = %d/%d'%(c, d))

                    # k in N_Even for single layer winding and k in N for double layer winding
                    # list_of_k = [k for k in range(1, MAX_K) if k%2==0] if layers==1 else range(1, MAX_K) # this turns out to be wrong
                    list_of_k = [k for k in range(1, MAX_K)]
                    for k in list_of_k:
                        print('\t'*4 + 'k=%3d'%(k), end='')

                        Qr = k*d
                        print('   ' + 'Qr=%3d'%(Qr), end='')
                        if Qr not in valid_Qr:
                            print()
                            continue
                        
                        nl = Qr/2 + 1 if layers == 1 else 0
                        print('   ' + 'nl=%2d'%(nl), end='')

                        y = k*c
                        print('   ' + 'y=%3d'%(y), end='')

                        # This step is core. Pick an m to enforce z=1 (q=z/n).
                        m = Fraction(Qr, 2*p).numerator
                        n = Fraction(Qr, 2*p).denominator
                        print('   ' + 'm/n=%3d/%d'%(m,n), end='')

                        # fractional slot winding?
                        print('   ISW ' if n ==1 else '   _FSW', end='')

                        # SPP
                        print('   ' + 'q=z/n=%d/%d'%(Fraction(Qr,2*p*m).numerator, Fraction(Qr,2*p*m).denominator), end='')

                        # coil pitch
                        p_aster_gamma = p*y/Qr*360
                        while p_aster_gamma > 360: p_aster_gamma -= 360
                        if p_aster_gamma > 180: p_aster_gamma = 360 - p_aster_gamma
                        print('   ' + 'p*gamma=%3g'%(p_aster_gamma), end='')

                        # winding factor = pitch factor because z=1, i.e., no distribution, i.e., one coil per coil group
                        gamma = y/Qr * 2*math.pi
                        print('   ' + 'kw_h = sin(h*gamma/2) =', end='')
                        for h in range(1, 10):
                            kw_h = math.sin(h*gamma/2)
                            print(f'{kw_h:.2f}', end=', ')
                        print('...')

    def get_ValidSet_of_Qr(self, Qs, p):
        valid_Qr = []
        for Qr in range(2, 2*Qs, 2):
            if abs(Qs - Qr) == 2 \
            or abs(Qs - Qr) == 2*p+1 \
            or abs(Qs - Qr) == 2*p-1 \
            or Qr == Qs \
            or Qr == 0.5*Qs \
            or Qr == 2*Qs:
            # or abs(Qs - Qr) == 2*p+2 \
            # or abs(Qs - Qr) == 2*p-2 \
            # or Qr > 1.25*Qs \
            # or Qr       % (6*p) == 0 \
            # or (Qr+2*p) % (6*p) == 0 \
            # or (Qr-2*p) % (6*p) == 0 \
            # or abs(Qr-Qs) == 2*p \
            # or abs(Qr-2*Qs) == 2*p \
            # or abs(Qr-Qs) == p \
            # or abs(Qr-0.5*Qs) == p:
                continue
            else:
                valid_Qr.append(Qr)
        return valid_Qr

    def get_design(self, layers, ps, p, k1, k):
        # c/d = reduced( k1/ps )
        self.c = c = Fraction(k1, ps).numerator
        self.d = d = Fraction(k1, ps).denominator
        print('c/d\t= %d/%d'%(c, d))

        self.Qr = Qr = k*d
        print(f'Qr\t= {Qr}')
        print(f'p\t= {p}')
        print(f'ps\t= {ps}')
        print(f'k1\t= {k1}')
        print(f'k\t= {k}')
        print(f't\t= {math.gcd(Qr,p)}')

        self.nl = Qr/2 + 1 if layers == 1 else 0
        print('nl\t= %2d'%(self.nl))

        y = k*c
        self.coil_pitch_y = coil_pitch_y = y
        print('y\t= %3d'%(y))

        # This step is core. Pick an m to enforce z=1 (q=z/n).
        self.m = m = Fraction(Qr, 2*p).numerator
        self.n = n = Fraction(Qr, 2*p).denominator
        print('m/n\t= %3d/%d'%(m,n), '(this is consistent with Pyrhonen@(2.83): m cannot be multiples of n.)')

        # fractional slot winding?
        print('ISW' if n==1 else 'FSW')

        # SPP
        self.q = Qr/(2*p*m)
        self.z = z = Fraction(Qr,2*p*m).numerator
        self.n = n = Fraction(Qr,2*p*m).denominator
        print(f'q=z/n={z:d}/{n:d} = {Qr:d}/(2*{p}*{m})')

        # coil pitch
        p_aster_gamma = p*y/Qr*360
        while p_aster_gamma > 360: p_aster_gamma -= 360
        if p_aster_gamma > 180: p_aster_gamma = 360 - p_aster_gamma
        print('p*gamma\t=%3g'%(p_aster_gamma))

        # winding factor = pitch factor because z=1, i.e., no distribution, i.e., one coil per coil group
        self.gamma = gamma = y/Qr * 2*math.pi
        print('kw_h = sin(h*gamma/2) =', end='')
        for h in range(1, 10):
            kw_h = math.sin(h*gamma/2)
            print(f'{kw_h:.2f}', end=', ')
        print('...')


        if layers == 1:
            number_of_coils = Qr/2 # each coil takes up two slots.
        elif layers == 2:
            number_of_coils = Qr
        self.number_of_coils = number_of_coils
        print('number_of_coils =', number_of_coils)
        print('number of coils per phase =', number_of_coils / m)


        # --------------------------------------------------- Tentative winding layout
        char_bias = 97
        sU = ''
        sL = ''
        self.pairs = []
        if y*2 == Qr and layers==1:
            for i in range(y):
                sU += ' | ' + chr(char_bias+i)
                self.pairs.append( (i+1, i+y+1) )
            for ind, _ in enumerate(range(y, Qr)):
                sU += ' | ' + chr(char_bias+ind)
        elif y*3 == Qr and layers==2:
            for i in range(Qr-y):
                sU += ' | ' + chr(char_bias+i)
                # sL += ' | ' + chr(char_bias+Qr+i-y)
                self.pairs.append( (i+1, i+y+1) )
            for i in range(Qr-y, Qr):
                sU += ' | ' + chr(char_bias+i)
                # sL += ' | ' + chr(char_bias+i+y)
                self.pairs.append( (i+1, i+y-Qr+1) )
        else:
            raise Exception('Not implemented.')
        print(f'----\nWinding Layout (m={m:d} phases):')
        print("".join([f' |{i+1:2d}' for i in range(Qr)]))
        print(sU)
        print(sL)
        print('\n'*3)
        print('----------Winding layout transcript used in winding_layout.py is as follows:')
        print(r'''
        if Qr == %d \
        and p == %d \
        and ps == %d \
        and coil_pitch_y == %d:'''%(Qr,p,ps,coil_pitch_y), end='\n\t\t\tself.pairs = ')
        if layers == 2:
            print(design.reduce_to_single_layer())
        else:
            print(self.pairs)

    def reduce_to_single_layer(self):
        list_groups = []
        # print('\n|||', self.pairs)
        # quit()
        for _, target in enumerate(self.pairs):
            bool_continue = False
            for group in list_groups:
                # print(target[0], 'in', group, '?') # 11 in (1, 11, 21) ?
                if target[0] in group:
                    bool_continue = True # avoid to include (1, 11, 21) and (11, 21, 1) at the same time
            if bool_continue == True:
                continue

            # print(target, end=' | ')
            bool_appended = False
            group = list(target) # [1, 11]  |  [11, 21]
            for i in range(_+1, len(self.pairs)):
                pair = self.pairs[i] # (2, 12), (11, 21)  |  (1, 11)
                if pair[0] in group: # 2, 11  |  1
                    if pair[1] not in group:
                        group.append(pair[1])
                        bool_appended = True
                if pair[1] in group: # 12, 21  |  11
                    if pair[0] not in group:
                        group.append(pair[0])
                        bool_appended = True
            # print(group)
            if bool_appended:
                list_groups.append(tuple(group))

        # print('Before:', self.pairs)
        # print('After:', list_groups)
        return list_groups

if __name__ == '__main__':
    layers = 1; QS = 36
    layers = 1; QS = 24
    layers = 2; QS = 24
    # layers = 2; QS = 18

    print(f'For Qs={QS:d}, find valid set of pole specific winding with neutral plate.')
    from fractions import Fraction
    import math
    MAX_K = 20+1
    MAX_PS = 4+1




    import os, sys
    # sys.stdout = open(os.devnull, "w")
    design = ValidSet_of_PoleSpecificWindingWithNeutralPlate(layers=layers)
    sys.stdout = sys.__stdout__

    # print('\n'+'~*'*40)
    # design.get_design(layers=1, ps=2, p=3, k1=1, k=12)
    # print('\n'+'~*'*40)
    # design.get_design(layers=1, ps=2, p=3, k1=1, k=10)

    # print('\n'+'~*'*40)
    # design.get_design(layers=2, ps=3, p=2, k1=1, k=10)

    print('\n'+'~*'*40)
    design.get_design(layers=2, ps=3, p=2, k1=1, k=6)


