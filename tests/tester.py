sys.path.append('/usr/share/inkscape/extensions')
import inkex
import pathmodifier
import simplepath
import numpy
import shapereco

#import cubicsuperpath, simplepath, simpletransform
#e=pathmodifier.PathModifier()
#shapereco.set_debug_on()
e = shapereco.ShapeReco()
print sys.argv

testA = numpy.array(range(5))

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
tan = shapereco.buildTangents( p )

from numpy import array
e.parse('tests/failureL.svg')
#e.parse('tests/picture.svg')
#e.parse('doc/demo.svg')



def readTangents( nid , aSize=0.1, doPrint=False):
    el = e.getElementById(nid)
    if el is None:
        print "Cant find ", nid
        return
    points, svgL=shapereco.toArray(simplepath.parsePath(el.get('d')))
    #print points
    tg=shapereco.buildTangents( points )
    res = e.checkForCircle(points, tg)
    e.points = points
    e.tg =tg
    (deltasD,angles, fits , dist_to_fit) = e.temp
    clList=shapereco.clusterValues(angles,aSize,refScaleAbs='abs')
    clList.sort(key=lambda cl:len(cl))    
    #print [len(cl) for cl in clList]
    #print (len(clList[-1])+len(clList[-2]))/float(len(angles))
    #print res

    deltaA = angles[1:] - angles[:-1]
    deltasDD =  (deltasD[1:] -deltasD[:-1])
    deltasDD[numpy.where(deltasDD==0.)] = 1e-5*deltasD[0]
    dAdD = abs(deltaA/deltasDD)
    belowT, count = True,0
    x,y,w,h = shapereco.computeBox(points)
    for i,v in enumerate(dAdD):
        if v>6 and belowT:
            count+=1
            belowT = False
        belowT= (v<6)

    if doPrint:
        print 'cl Lsit =',clList
        print 'dist to fit=', dist_to_fit
    #return dist_to_fit[int(len(angles)*0.9)], (len(clList[-1])+len(clList[-2]))/float(len(angles)) , len(clList)/float(len(angles))
    #return dist_to_fit[int(len(angles)*0.6)], len([c for c in clList if len(c)<3])/float(len(angles)), len(clList)/float(len(angles))
    return deltasD[-1]/3.14, count, numpy.sum(deltasDD[numpy.where(dAdD<0.4)])/(deltasD[-1]-deltasD[0])
    
def testMany(aSize=0.1):


    rect=[ "path4725" ,"path5186" ,"path4494" ,"path6787", "path16095"]
    
    cir=["path4511", "path4637" ,"path6181"    ,"path5062"    ,"path7415" ,"path9773" , "path9066", "path5531", "path5535", "path9291"]
    noconv = [   "path9362",  "path4661"] 
    tri = ["path9285"]
    
    lrec = [ readTangents(o,aSize) for o in rect ]
    lcir =[ readTangents(o,aSize) for o in cir ]
    lnoconv =[ readTangents(o,aSize) for o in noconv ]
    ltri =[ readTangents(o,aSize) for o in tri ]
    print "rec = ", [ "(%.2f , %.2f, %.2f)"%r for r in lrec]
    print "cir = ", [ "(%.2f , %.2f, %.2f)"%r for r in lcir]
    print "ncv = ", [ "(%.2f , %.2f, %.2f)"%r for r in lnoconv]
    print "ntri = ", [ "(%.2f , %.2f, %.2f)"%r for r in ltri]

def plotTangents( nid ):
    from  matplotlib import pyplot as plt
    readTangents(nid)
    (deltasD,angles, a  , dist_to_fit) = e.temp
    angles_vs_d = numpy.stack([ deltasD, angles] ,1)
    seg =shapereco.regLin( angles_vs_d )
    fits = -(seg.a*angles_vs_d[:,0]+seg.c)/seg.b
    
    deltaA = angles[1:] - angles[:-1]
    deltasDD =  (deltasD[1:] -deltasD[:-1])
    deltasDD[numpy.where(deltasDD==0.)] = 1e-5*deltasD[0]
    dAdD = abs(deltaA/deltasDD)
    plt.plot(deltasD, angles)
    plt.plot(deltasD, fits)
    #plt.plot(deltasD, dist_to_fit)
    plt.plot(deltasD[:-1], dAdD)
    axes = plt.gca()
    axes.set_ylim([-3.2, 6.2])
    plt.grid()
    plt.show()



def plotSegs(shape, *indices):
    from  matplotlib import pyplot as plt
    print indices
    if indices==(): indices=range(len(shape.listOfPaths))
    lines=[]
    for i in indices:
        s = shape.listOfPaths[i]
        print [s.point1[0], s.pointN[0] ], [s.point1[1], s.pointN[1]  ]
        lines +=[ plt.plot( [s.point1[0], s.pointN[0] ], [s.point1[1], s.pointN[1] ] ) ]
    plt.grid()
    plt.show()

def plotSegPointVsLineEq(s):
    from  matplotlib import pyplot as plt
    y = numpy.array([s.point1[1], s.pointN[1] ])
    x = numpy.array([s.point1[0], s.pointN[0] ])
    plt.plot( x,  y )
    if s.b != 0:
        plt.plot(x, (-x*s.a-s.c)/s.b )
    else:
        print ' !!!! b=0'
        plt.plot(x, [0,0] )
    plt.grid()
    plt.show()
