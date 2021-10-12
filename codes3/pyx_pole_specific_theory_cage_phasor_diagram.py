import PyX_classes, PyX_Utility
import os, pyx
from pylab import np

''' TEC-ISMB-2021 鼠笼相量图
'''
if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THIck

    Q_r_prime = 8
    alpha_c = 2*np.pi / Q_r_prime
    for i in range(Q_r_prime):
        radius = 10
        x = radius *  np.cos(i*alpha_c)
        y = radius * -np.sin(i*alpha_c)
        pu.pyx_arrow((x,y))
        radius = 9
        x = radius *  np.cos(i*alpha_c+15/180*np.pi)
        y = radius * -np.sin(i*alpha_c+15/180*np.pi)
        pu.pyx_text((x,y), r'$\bar U_{r%d,n}$'%(i+1), settings=['blue'], scale=1.5)
        radius = 11
        x = radius *  np.cos(i*alpha_c+0/180*np.pi)
        y = radius * -np.sin(i*alpha_c+0/180*np.pi)
        pu.pyx_text((x,y), str(i+1), settings=[], scale=1.5)

    radius = 3
    x = radius *  np.cos(i*alpha_c)
    y = radius * -np.sin(i*alpha_c)
    pu.pyx_arc((x,y), (radius,0), settings=['earrow'])
    radius = 5
    x = radius * np.cos(0.5*alpha_c)
    y = radius * np.sin(0.5*alpha_c)
    pu.pyx_text((x,y), r'$n\alpha_c$', scale=1.5)


    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/{os.path.basename(__file__)[:-3]}'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    os.system(f'sumatraPDF2 {fname}.pdf')
    # os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')
