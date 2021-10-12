import PyX_classes, PyX_Utility
import os, pyx
from pylab import np

''' TEC-ISMB-2021 鼠笼相量图
'''
if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THIck

    pu.pyx_marker((0,0))

    radius = 10
    x = radius
    y = 0
    pu.pyx_arrow((x,y))
    p, ps = 3, 4
    Qr = 20
    alpha_c = 360 / Qr
    pu.pyx_text((x-3+0,   y+1), r'$1,$', scale=1.5)
    print(f'alpha={ps*0*alpha_c}')
    pu.pyx_text((x-3+1.5, y+1), r'$6,$', scale=1.5)
    print(f'alpha={ps*5*alpha_c}')
    pu.pyx_text((x-3+0,   y-1), r'$11,$', scale=1.5)
    print(f'alpha={ps*10*alpha_c}')
    pu.pyx_text((x-3+1.5,   y-1), r'$16$', scale=1.5)
    print(f'alpha={ps*15*alpha_c}')


    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/{os.path.basename(__file__)[:-3]}a'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    # os.system(f'sumatraPDF2 {fname}.pdf')
    # os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')

    print('save to:', fname)

if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THIck

    pu.pyx_marker((0,0))

    radius = 10
    x_ori = radius
    y_ori = 0
    p, ps = 3, 4
    Qr = 20
    alpha_c = 360 / Qr
    print(f'alpha_c={alpha_c}')

    pu.pyx_arrow((x_ori,y_ori))
    pu.pyx_text((x_ori-1.5,   y_ori-1), r'$1$', scale=1.5)
    print(f'alpha={p*0*alpha_c}')

    def rotate_arrow_and_text(BarNo):
        alpha = p*(BarNo-1)*alpha_c / 180 * np.pi
        x = x_ori* np.cos(alpha) + y_ori*-np.sin(alpha)
        y = x_ori* np.sin(alpha) + y_ori* np.cos(alpha)
        pu.pyx_arrow((x,y))
        pu.pyx_text(( np.sign(x)*(abs(x)-1), 
                      np.sign(y)*(abs(y)-1)
                     ), rf'${BarNo}$', scale=1.5)
        print(f'alpha={alpha/np.pi*180}')

    rotate_arrow_and_text(6)
    rotate_arrow_and_text(11)
    rotate_arrow_and_text(16)


    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/{os.path.basename(__file__)[:-3]}b'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    # os.system(f'sumatraPDF2 {fname}.pdf')
    # os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')

    print('save to:', fname)
