import sys
print(sys.version)
import PyX_classes, PyX_Utility
import os, pyx
from pylab import np

''' TEC-ISMB-2021 鼠笼磁场时域图
'''
if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    scale = 1.25
    x = 17*scale
    y = 8*scale

    ## 画坐标图
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.thin
    PyX_classes.quadrature_coordinate(ylength=y/2+0.5, xlength=x).draw(pu, ymirror=True)

    ## 画波形
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
    graph = PyX_classes.function_as_graph()
    hat_B_delta_ps = 2
    hat_B_delta_p  = 4
    freq = 0.2
    period = 1/freq
    def fblue(t):
        return  - hat_B_delta_ps*np.sin(2*np.pi*(freq) * t )
    def fgreen(t):
        return  - hat_B_delta_p*np.sin(2*np.pi*(freq*3/4) * t )
    graph.draw( pu, [(t, fblue(t)) for t in np.arange( 0, 4*period+.01, period/40)], settings=['blue'] )
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashdotted
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THick
    graph.draw( pu, [(t, fgreen(t)) for t in np.arange( 0, 4*period+.01, period/20)], settings=['darkgreen'] )


    ## 画导条
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THIck
    radius = 0.5
    pu.pyx_circle(radius, center=[period*0.5+0*period, radius], settings=['RawSienna'])
    pu.pyx_circle(radius, center=[period*0.5+1*period, radius], settings=['RawSienna'])
    pu.pyx_circle(radius, center=[period*0.5+2*period, radius], settings=['RawSienna'])
    pu.pyx_circle(radius, center=[period*0.5+3*period, radius], settings=['RawSienna'])


    ## 标注
    # 纵轴
    pu.pyx_text([2, y/2], r'$B_\delta(\alpha)$', size=3)
    pu.pyx_text([x-0.5, 1.0], r'$\alpha$', size=3)
    def xAnnotation(left, right, yline, ytext, label):
        pu.pyx_line([left,  ytext],
                    [left,  yline])
        pu.pyx_line([right, ytext],
                    [right, yline])
        pu.pyx_arrow_both_ends( [left, ytext],
                                [right, ytext])
        pu.pyx_text( [(left+right)/2, ytext+1], label, size=3)
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashed
    xAnnotation(0.5*period, 1.5*period,
                yline = 0,
                ytext = -5,
                label=r'$\alpha_c$')


    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/{os.path.basename(__file__)[:-3]}'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    print('save to:', fname)

    ## 裁切、打开文件
    os.system(f'sumatraPDF {fname}.pdf')
    # os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')
