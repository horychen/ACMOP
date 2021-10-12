import PyX_classes, PyX_Utility
import os, pyx
from pylab import np

''' SlessInv saturation time based correction 示意图 (Ideal Situation Only)
'''
if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    scale = 1.25
    x = 20*scale
    y = 10*scale
    origin = [0, 0]

    freq1 = 1.75 # Change here
    freq2 = 1.75 # Change here

    pt0 = origin
    pt1 = [x/freq1/2, 0]
    pt2 = [pt1[0]+x/freq2/2, 0]
    pt3 = [pt2[0]+x/freq2/2, 0]

    ell_limit = y/2-0.9

    def text_bias(p, xbias=-1, ybias=1):
        return [p[0]+xbias, p[1]+ybias]

    ## 画坐标图
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.thin
    PyX_classes.quadrature_coordinate(ylength=y/2+0.5, xlength=x).draw(pu, ymirror=True, settings=['purple-blue'])

    ## 画限幅
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.Thick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashed
    pu.pyx_line([0,  ell_limit], [x-5,  ell_limit], settings=['RawSienna'])
    pu.pyx_line([0, -ell_limit], [x-8, -ell_limit], settings=['RawSienna'])
    pu.pyx_text([x-5+3,  ell_limit], r'$ \ell= K_{\rm Active}$', settings=['RawSienna'])
    pu.pyx_text([x-6+3, -ell_limit], r'$-\ell=-K_{\rm Active}$', settings=['RawSienna'])

    ## 波形最大值
    BIAS = 2.0
    psi_mu_max =  (y/2)
    psi_mu_min = -(y/2)
    psi_2_max  = psi_mu_max -BIAS
    psi_2_min  = -ell_limit
    pu.pyx_line([0,  psi_2_max], [(pt1[0]-pt0[0])/2, psi_2_max])
    pu.pyx_text(  [-2.2,  psi_2_max],   r'$\psi_{\alpha2,\max}$')
    pu.pyx_marker([   0,  psi_2_max]);
    ## 波形的最小值
    pu.pyx_text(  [-2.2, -ell_limit], r'$\psi_{\alpha2,\min}$') #(pt2[0]-pt1[0])/2+pt1[0]
    pu.pyx_marker([   0, -ell_limit]);


    ## 画波形
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.Thick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
    graph = PyX_classes.function_as_graph()
    def f1(t):
        return  psi_mu_max*np.sin(2*np.pi*(freq1/x) * t ) - BIAS
    list_of_saturated_time_before_t2 = []
    list_of_saturated_time_after_t2 = []
    def f2(t):
        val  =  psi_mu_min*np.sin(2*np.pi*(freq2/x) *(t-pt1[0])) - BIAS
        if abs(val) > ell_limit: 
            if t<pt2[0]:
                list_of_saturated_time_before_t2.append(t)
            else:
                list_of_saturated_time_after_t2.append(t)
            return np.sign(val)*ell_limit
        else:
            return val
    graph.draw( pu, [(t, f1(t)) for t in np.arange(     0, pt1[0]+.1, 0.05)] )
    # graph.draw( pu, [(t, f2(t)) for t in np.arange(pt1[0], pt3[0]+.1, 0.25)] ) # Change here
    graph.draw( pu, [(t, f2(t)) for t in np.arange(pt1[0], pt2[0]+.1, 0.05)] )
    graph.draw( pu, [(t, f1(t)) for t in np.arange(pt2[0], pt3[0]+.1, 0.05)] )

    ## 画饱和时间
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.Thick
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashed
    t_alpha_sat_max = 0
    t_alpha_sat_min = max(list_of_saturated_time_before_t2) - min(list_of_saturated_time_before_t2)
    def xAnnotation(left, right, yline, ytext, label, ytextbias=1):
        pu.pyx_line([left,  ytext],
                    [left,  yline])
        pu.pyx_line([right, ytext],
                    [right, yline])
        pu.pyx_arrow_both_ends( [left, ytext],
                                [right, ytext])
        pu.pyx_text( [(left+right)/2, ytext+ytextbias], 
                     label )
    # 下饱和时间
    tSatMinLeft  = min(list_of_saturated_time_before_t2)
    tSatMinRight = max(list_of_saturated_time_before_t2)
    xAnnotation(tSatMinLeft, tSatMinRight, 
                yline = -ell_limit,
                ytext = -ell_limit/2-1,
                label=r'$t_{\alpha,\mathrm{sat,min}}$')
    # 第二个上饱和时间
    # tSatMinLeft  = min(list_of_saturated_time_after_t2)
    # tSatMinRight = max(list_of_saturated_time_after_t2)
    # xAnnotation(tSatMinLeft, tSatMinRight, 
    #             yline = ell_limit,
    #             ytext = ell_limit+1,
    #             label=r'$t_{\alpha,\mathrm{sat,max}}$',
    #             ytextbias=1)

    # 第一个上饱和时间和它的箭头，我还给text加了一个tbox
    pu.pyx_line(        [pt1[0]/2+6, ell_limit+1.5], [(pt1[0]-pt0[0])/2, psi_2_max]) # draw annotation for r'$t_{\alpha,\mathrm{sat,max}}=0$'
    pu.pyx_line(        [pt1[0]/2+6, ell_limit+1.5], [(pt2[0]+pt3[0])/2, psi_2_max]) # draw annotation for r'$t_{\alpha,\mathrm{sat,max}}=0$'
    tbox = pu.pyx_text( [pt1[0]/2+6, ell_limit+1.5], 
                        r'$t_{\alpha,\mathrm{sat,max}}=0$',
                        BoxColor='warm-grey')
    # pu.pyx_text(       [x-5, ell_limit+1.5], r'$t_{\alpha,\mathrm{sat,max}}=0$')
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THick
    # pu.pyx_arrow( [pt1[0]/2, -0.2+ell_limit+1.0],
    #               [pt1[0]/2,  0.2+psi_2_max],
    #               settings=['warm-grey'] )


    ## 标注
    # 纵轴
    pu.pyx_text([0, y/2+1.5], r'$\psi_{\alpha2}(t)$', settings=['purple-blue'])

    # 时间轴
    pu.pyx_text([x-1, -1], r'$t$', settings=['purple-blue'])
    pu.pyx_marker((pt0[0]+0.75, pt0[1])); pu.pyx_text(text_bias(pt0,xbias=-0.5), r'$t_0$')
    pu.pyx_marker((pt1[0]-0.75, pt1[1])); pu.pyx_text(text_bias(pt1,xbias=0.5), r'$t_1$')
    pu.pyx_marker((pt2[0]+0.75, pt2[1])); pu.pyx_text(text_bias(pt2,xbias=-0.5), r'$t_2$')
    pu.pyx_marker((pt3[0]-0.75, pt3[1])); pu.pyx_text(text_bias(pt3,xbias=0.5), r'$t_3$')
    xAnnotation(pt1[0]-0.75, pt2[0]+0.75,
                yline = 0,
                ytext = 2,
                label=r'$\Delta t$')

    # 重心标注
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
    pu.pyx_marker([ 0,       (psi_2_max+psi_2_min)/2]);
    pu.pyx_text(  [-2, -2.25+(psi_2_max+psi_2_min)/2], 
                  r'$\frac{\psi_{\alpha2,\max}+\psi_{\alpha2,\min}}{2}$',
                  BoxColor='warm-grey')
    pu.pyx_arrow( [-2.5, -1+(psi_2_max+psi_2_min)/2],
                  [-0.2,-0.1+(psi_2_max+psi_2_min)/2] )
    # 中分线
    PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.dashdotted
    # PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.Thick
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
    # 线段相等符号
    def rotate_by_pivot(p, pivot, angle):
        result = [None]*2
        result[0] = (p[0]-pivot[0])* np.cos(angle) + (p[1]-pivot[1])*np.sin(angle) + pivot[0]
        result[1] = (p[0]-pivot[0])*-np.sin(angle) + (p[1]-pivot[1])*np.cos(angle) + pivot[1]
        return result
    if True:
        PyX_Utility.global_settings['linestyle'] = pyx.style.linestyle.solid
        def draw_equal_sign(pivot):
            pA = [pivot[0]+0.25, pivot[1]+0.15]
            pB = [pivot[0]-0.25, pivot[1]+0.15]
            pu.pyx_line(rotate_by_pivot(pA, pivot, 15/180*np.pi), 
                        rotate_by_pivot(pB, pivot, 15/180*np.pi), 
                        settings=['tint-red'])
            pC = [pivot[0]+0.25, pivot[1]-0.15]
            pD = [pivot[0]-0.25, pivot[1]-0.15]
            pu.pyx_line(rotate_by_pivot(pC, pivot, 15/180*np.pi), 
                        rotate_by_pivot(pD, pivot, 15/180*np.pi), 
                        settings=['tint-red'])

        # 上
        pivot = [pt1[0]/2, psi_2_max/2]
        draw_equal_sign(pivot)
        # 下
        pivot = [pt1[0]/2, psi_2_min/2]
        draw_equal_sign(pivot)


    ## 理想情况和非理想情况区间标注
    # pu.pyx_horizental_arrow_both_ends_with_text_inside(
    #     'ideal situation',
    #     [pt0[0], -y/2-0.5],
    #     [pt2[0], -y/2-0.5],
    #     BoxColor='warm-grey')
    # pu.pyx_horizental_arrow_both_ends_with_text_inside( 
    #     'non-ideal situation',
    #     [pt1[0], -y/2-2.0],
    #     [pt3[0], -y/2-2.0])
        # BoxColor='tint-yellow')



    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/pyx_sat_time_corr_annotated_ideal'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    # os.system(f'sumatraPDF2 {fname}.pdf')
    os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')
