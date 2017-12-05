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
e = shapereco.Shapereco()
tan = shapereco.buildTangeants( p )
