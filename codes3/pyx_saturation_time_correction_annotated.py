import PyX_classes, PyX_Utility
import os, pyx
from pylab import np

''' IMIFE saturation time based correction 示意图
'''
if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    scale = 1.25
    x = 20*scale
    y = 10*scale
    origin = [0, 0]
    pt0 = origin
    pt1 = [x/4, 0]
    pt2 = [x*3/4, 0]

    ell_limit = y/2-0.9

    def text_bias(p, xbias=-1, ybias=1):
        return [p[0]+xbias, p[1]+ybias]

    ## 画坐标图
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.thin
    PyX_classes.quadrature_coordinate(ylength=y/2+0.5, xlength=x).draw(pu, ymirror=True)

    ## 画限幅
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.Thick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashed
    pu.pyx_line([0,  ell_limit], [x-8,  ell_limit], settings=['RawSienna'])
    pu.pyx_line([0, -ell_limit], [x-8, -ell_limit], settings=['RawSienna'])
    pu.pyx_text([x-5,  ell_limit], r'$\ell=\psi^*$', settings=['RawSienna'])
    pu.pyx_text([x-5, -ell_limit], r'$-\ell=-\psi^*$', settings=['RawSienna'])

    ## 波形最大值
    psi_mu_max =  (y/2-2)
    psi_mu_min = -(y/2)
    psi_2_max  = psi_mu_max
    psi_2_min  = -ell_limit
    pu.pyx_line([0,  psi_mu_max], [(pt1[0]-pt0[0])/2, psi_mu_max])
    pu.pyx_text(  [-2.2,  psi_mu_max],   r'$\psi_{\alpha2,\max}$')
    pu.pyx_marker([   0,  psi_mu_max]);
    ## 波形的最小值
    pu.pyx_text(  [-2.2, -ell_limit], r'$\psi_{\alpha2,\min}$') #(pt2[0]-pt1[0])/2+pt1[0]
    pu.pyx_marker([   0, -ell_limit]);


    ## 画波形
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.Thick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
    graph = PyX_classes.function_as_graph()
    def f1(t):
        return  psi_mu_max*np.sin(2*np.pi*(2/x) * t )
    list_of_saturated_time = []
    def f2(t):
        val  =  psi_mu_min*np.sin(2*np.pi*(1/x) *(t-pt1[0]))
        if abs(val) > ell_limit: 
            list_of_saturated_time.append(t)
            return np.sign(val)*ell_limit
        else:
            return val
    graph.draw( pu, [(t, f1(t)) for t in np.arange(     0, pt1[0]+.1, 0.25)] )
    graph.draw( pu, [(t, f2(t)) for t in np.arange(pt1[0], pt2[0]+.1, 0.25)] )

    ## 画饱和时间
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.Thick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashed
    t_alpha_sat_max = 0
    t_alpha_sat_min = max(list_of_saturated_time) - min(list_of_saturated_time)
    tSatMinLeft  = min(list_of_saturated_time)
    tSatMinRight = max(list_of_saturated_time)
    def xAnnotation(left, right, yline, ytext, label):
        pu.pyx_line([left,  ytext],
                    [left,  yline])
        pu.pyx_line([right, ytext],
                    [right, yline])
        pu.pyx_arrow_both_ends( [left, ytext],
                                [right, ytext])
        pu.pyx_text( [(left+right)/2, ytext+1], 
                     label )
    xAnnotation(tSatMinLeft, tSatMinRight, 
                yline = -ell_limit,
                ytext = -ell_limit/2,
                label=r'$t_{\alpha,\mathrm{sat,min}}$')


    ## 标注
    # 纵轴
    pu.pyx_text([0, y/2+1.5], r'$\psi_{\alpha2}$')

    # 时间轴
    pu.pyx_marker(pt0); pu.pyx_text(text_bias(pt0), r'$t_0$')
    pu.pyx_marker(pt1); pu.pyx_text(text_bias(pt1), r'$t_1$')
    pu.pyx_marker(pt2); pu.pyx_text(text_bias(pt2,xbias=1), r'$t_2$')
    xAnnotation(pt1[0], pt2[0],
                yline = 0,
                ytext = 2,
                label=r'$\Delta t$')

    # 上饱和时间
    # pu.pyx_text([pt1[0]/2+3, ell_limit+1.5], r'$t_{\alpha,\mathrm{sat,max}}=0$')
    pu.pyx_text(       [x-5, ell_limit+1.5], r'$t_{\alpha,\mathrm{sat,max}}=0$')

    # 重心标注
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
    pu.pyx_marker([ 0,    (psi_2_max+psi_2_min)/2]);
    pu.pyx_text(  [-2, -2+(psi_2_max+psi_2_min)/2], 
                  r'$\frac{\psi_{\alpha2,\max}+\psi_{\alpha2,\min}}{2}$')
    pu.pyx_arrow( [-2.5, -1+(psi_2_max+psi_2_min)/2],
                  [-0.2,-0.1+(psi_2_max+psi_2_min)/2] )
    # 中分线
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashdotted
    # 竖线
    pu.pyx_line([pt1[0]/2, psi_2_max], 
                [pt1[0]/2, psi_2_min], settings=['tint-red'])
    # 横线
    pu.pyx_line([0,        (psi_2_max+psi_2_min)/2], 
                [pt1[0]/2, (psi_2_max+psi_2_min)/2], settings=['tint-red'])
    # 直角横线
    size = 0.5
    pu.pyx_line([pt1[0]/2-size, -size+(psi_2_max+psi_2_min)/2], 
                [pt1[0]/2,      -size+(psi_2_max+psi_2_min)/2], settings=['tint-red'])
    # 直角竖线
    pu.pyx_line([pt1[0]/2-size, -size+(psi_2_max+psi_2_min)/2], 
                [pt1[0]/2-size,      +(psi_2_max+psi_2_min)/2], settings=['tint-red'])






    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/pyx_sat_time_corr_annotated'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    os.system(f'sumatraPDF2 {fname}.pdf')
    os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')
