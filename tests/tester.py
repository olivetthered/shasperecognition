sys.path.append('/usr/share/inkscape/extensions')
import inkex
import pathmodifier
import numpy
import shapereco

#import cubicsuperpath, simplepath, simpletransform
#e=pathmodifier.PathModifier()
shapereco.set_debug_on()
#e = shapereco.Shapereco()
print sys.argv


def test(n):
    if isinstance(n,str):
        n = [n]
    e.defPathName = n
    e.effect()
    print n
    e.dumpSegs()


def multiTest():
    for p in ['path2987',
              'path2985',
              'path3013',
              'path3045',              
              'path3043',
              ] :
        test(p)
        
        #v = e.seg.variance()
        #n = len(e.seg.points)
        print ' ****************** '
        #print p,' v=%.2f / %.2f  r=%.2f  r*sqrt(n)=%.2f  npoints=%d'%(v, e.seg.length, v/e.seg.length,v/e.seg.length*numpy.sqrt(n) ,n  )


def generateCircle(n=30):
    d = numpy.pi*2/n
    c = numpy.zeros((30,2))
    for i in range(n):
        c[i] = [numpy.cos(d*i), numpy.sin(d*i), ]
    return c
#e.affect(output=False)

ar = numpy.array
s1 = shapereco.Segment.from2Points(ar([-159.2  ,197.3]) , ar([-135.6  ,213.9]) )
s2 = shapereco.Segment.from2Points( ar([-135.2,  213.8]), ar( [-123.5,  221.9]) )

p=  numpy.array( [ [1,1] ,[2,2] ,[4,4],
                   [4.5, 4.2],
                   [5, 4.6],
                   [6, 4.8],
                   [7, 5],
                   ] , dtype=float)
e = shapereco.ShapeReco()
tan = shapereco.buildTangeants( p )

from numpy import array

angle =  array([ 3.545,  3.902,  4.261,  4.633,  5.022,  5.429,  5.78 ,  6.283,
        6.283,  2.775,  2.731,  2.683,  2.634,  2.543,  2.45 ,  2.342,
        2.207,  2.062,  1.906,  1.728,  1.529,  1.323,  1.124,  0.943,
        0.763,  0.59 ,  0.433,  0.3  ,  0.175,  0.054, -0.076, -0.217,
       -0.365, -0.506, -0.651, -0.793, -0.929, -1.045, -1.146, -1.231,
       -1.314, -1.401, -1.495, -1.588, -1.685, -1.781, -1.877, -1.978,
       -2.071, -2.149, -2.213, -2.263, -2.307, -2.349, -2.381, -2.407,
       -2.43 , -2.454, -2.479, -2.497, -2.518, -2.55 , -2.597, -2.664])
deltasD= array([ 0.017,  0.04 ,  0.056,  0.091,  0.103,  0.115,  0.138,  0.173,
        0.185,  0.185,  0.281,  0.341,  0.55 ,  0.667,  0.719,  0.81 ,
        0.868,  1.011,  1.069,  1.129,  1.246,  1.41 ,  1.51 ,  1.699,
        1.881,  1.991,  2.066,  2.188,  2.259,  2.306,  2.376,  2.54 ,
        2.634,  2.797,  2.871,  2.996,  3.048,  3.165,  3.239,  3.301,
        3.354,  3.4  ,  3.449,  3.495,  3.542,  3.59 ,  3.661,  3.703,
        3.74 ,  3.782,  3.825,  3.874,  3.973,  4.023,  4.065,  4.107,
        4.149,  4.201,  4.234,  4.261,  4.294,  4.32 ,  4.353,  4.379])
dist_to_fit =  array([ 0.004,  0.004,  0.012,  0.017,  0.02 ,  0.024,  0.027,  0.03 ,
        0.035,  0.038,  0.041,  0.052,  0.056,  0.063,  0.07 ,  0.081,
        0.084,  0.088,  0.093,  0.093,  0.102,  0.14 ,  0.183,  0.192,
        0.194,  0.197,  0.204,  0.229,  0.265,  0.273,  0.291,  0.299,
        0.306,  0.317,  0.32 ,  0.332,  0.333,  0.345,  0.373,  0.378,
        0.389,  0.396,  0.427,  0.438,  0.488,  0.498,  0.52 ,  0.523,
        0.528,  0.531,  0.556,  0.578,  0.584,  0.589,  0.621,  0.624,
        0.918,  0.969,  1.029,  1.083,  1.456,  1.845,  2.406,  2.425])
