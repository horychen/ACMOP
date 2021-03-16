import PyX_classes, PyX_Utility
import os, pyx
from pylab import np

''' TEC-ISMB-2021 鼠笼相量图
'''
if __name__ == '!__main__':
    pu = PyX_Utility.PyX_Utility()

    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THIck

    pu.pyx_marker((0,0))

    radius = 10
    x = radius
    y = 0
    pu.pyx_arrow((x,y))
    pu.pyx_text((x-3+0,   y+1), r'$1,$', scale=1.5)
    pu.pyx_text((x-3+1.5, y+1), r'$2,$', scale=1.5)
    pu.pyx_text((x-3+0,   y-1), r'$3,$', scale=1.5)
    pu.pyx_text((x-3+3,   y-1), r'$4,\dots$', scale=1.5)


    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/{os.path.basename(__file__)[:-3]}a'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    os.system(f'sumatraPDF2 {fname}.pdf')
    # os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')



if __name__ == '__main__':
    pu = PyX_Utility.PyX_Utility()

    PyX_Utility.global_settings['linewidth'] = pyx.style.linewidth.THIck

    pu.pyx_marker((0,0))

    radius = 10
    x = radius
    y = 0
    pu.pyx_arrow((x,y))
    pu.pyx_text((x-3+1.5,   y-1), r'$1,3,\dots$', scale=1.5)
    pu.pyx_arrow((-x,y))
    pu.pyx_text((-x+4.5,   y-1), r'$2,4,\dots$', scale=1.5)


    ## 输出文件
    fname=fr'{os.path.dirname(__file__)}/../release/{os.path.basename(__file__)[:-3]}b'
    pu.cvs.writePDFfile(fname)
    # pu.cvs.writePDFfile(fname)
    # pu.cvs.writeSVGfile(fname)

    ## 裁切、打开文件
    os.system(f'sumatraPDF2 {fname}.pdf')
    # os.system(f'pdfcrop {fname}.pdf {fname}-crop.pdf')
