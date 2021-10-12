import PyX_classes, PyX_Utility
import os, pyx
from pylab import np

''' TEC-ISMB-2021 三维鼠笼转子 示意图
'''
if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THIck
    def RL_Circuit(x, y, voltage_source_label=None, inductance_label=None, resistance_label=None, current_label=None):
        ## 反电势作为电压源
        U1 = PyX_classes.voltage_source(center=(x, y), bool_vertical=True)
        U1.draw(pu, label=voltage_source_label)
        ## 电感
        L1 = PyX_classes.inductance(location=(x, y+3), bool_vertical=True)
        L1.draw(pu, label=inductance_label)
        ## 电感
        R1 = PyX_classes.resistance(location=(x, y+7), bool_vertical=True)
        R1.draw(pu, label=resistance_label)
        ## 连接点
        PBottom = PyX_classes.connection_point(location=(x, y-3))
        PBottom.draw(pu)
        PTop = PyX_classes.connection_point(location=(x, y+11))
        PTop.draw(pu)
        ## 连起来
        PyX_classes.easy_connect(pu, U1, L1)
        PyX_classes.easy_connect(pu,     L1, R1)
        PyX_classes.easy_connect(pu,     U1, PBottom)
        PyX_classes.easy_connect(pu,     R1, PTop)

        return U1, L1, R1, PTop, PBottom

    x = 0
    y = 0
    U1, L1, R1, PTop, PBottom = dict(), dict(), dict(), dict(), dict()
    U1[0], L1[0], R1[0], PTop[0], PBottom[0] = RL_Circuit(x-6,  y+2, voltage_source_label=r'$\bar U_{r1,n}$', inductance_label=r'$L_{r\sigma,n}}$', resistance_label=r'$R_{r,n}}$', current_label=r'\bar I_{r1,n}')
    U1[1], L1[1], R1[1], PTop[1], PBottom[1] = RL_Circuit(x+0,  y)
    U1[2], L1[2], R1[2], PTop[2], PBottom[2] = RL_Circuit(x+5,  y)
    U1[3], L1[3], R1[3], PTop[3], PBottom[3] = RL_Circuit(x+10, y+2)
    PyX_classes.easy_connect(pu, PBottom[0], PBottom[1])
    PyX_classes.easy_connect(pu, PBottom[2], PBottom[1])
    PyX_classes.easy_connect(pu, PBottom[2], PBottom[3])
    PyX_classes.easy_connect(pu, PTop[0], PTop[1])
    PyX_classes.easy_connect(pu, PTop[2], PTop[1])
    PyX_classes.easy_connect(pu, PTop[2], PTop[3])

    PTop[4] = PyX_classes.connection_point(location=(x-6, y+2+11+3)); PTop[4].draw(pu)
    PTop[5] = PyX_classes.connection_point(location=(x+0, y  +11+6)); PTop[5].draw(pu)
    PTop[6] = PyX_classes.connection_point(location=(x+5, y  +11+6)); PTop[6].draw(pu)
    PTop[7] = PyX_classes.connection_point(location=(x+10,y+2+11+3)); PTop[7].draw(pu)



    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/{os.path.basename(__file__)[:-3]}'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    # os.system(f'sumatraPDF2 {fname}.pdf')
    # os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')
