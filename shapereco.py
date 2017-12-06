#!/usr/bin/env python
'''
Copyright (C) 2017 , Pierre-Antoine Delsart

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
barraud@math.univ-lille1.fr

Quick description:
This extension uses all selected path, ignoring all other selected objects.
It tries to regularize hand drawn paths BY :
 - evaluating if the path is a full circle or ellipse
 - else finding sequences of aligned points and replacing them by a simple segment.
 - changing the segments angles to the closest remarkable angle (pi/2, pi/3, pi/6, etc...)
 - eqalizing all segments lengths which are close to each other
 - replacing 4 segments paths by a rectangle object if this makes sens (giving the correct rotation to the rectangle). 

Requires numpy.

'''

import sys
sys.path.append('/usr/share/inkscape/extensions')
import inkex
import simplepath
import gettext
_ = gettext.gettext




import numpy
numpy.set_printoptions(precision=3)
# *************************************************************
# a list of geometric helper functions 
def toArray(parsedList):
    """Interprets a list of [(command, args),...]
    where command is a letter coding for a svg path command
          args are the argument of the above command
    """
    interpretCommand = {
        'C' : lambda x, prevL : x[-2:], # bezier curve. Ignore the curve.
        'L' : lambda x, prevL : x[0:2],
        'M' : lambda x, prevL : x[0:2],
        'Z' : lambda x, prevL : prevL[0],
        }

    points =[]
    for i,(c, arg) in enumerate(parsedList):
        #debug('toArray ', i, c , arg)
        newp = interpretCommand[c](arg, points)
        points.append( newp)
    a=numpy.array( points )

    # Some times we have points *very* close to each other
    # these do not bring any meaning full info, so we remove them
    #
    x,y, w,h = computeBox(a)
    sizeC = 0.5*(w+h)
    deltas =  a[1:] - a[:-1] 
    rel_norms = numpy.sqrt(numpy.sum( deltas**2, 1 )) / sizeC
    keep = numpy.concatenate([numpy.where( rel_norms >0.005 )[0],numpy.array([len(a)-1])])

    return a[keep] , [ parsedList[i] for i in keep]

rotMat = numpy.matrix( [[1,-1],[1,1]] )/numpy.sqrt(2)
unrotMat = numpy.matrix( [[1,1],[-1,1]] )/numpy.sqrt(2)

def setupKnownAngles():
    pi = numpy.pi
    l = [ i*pi/8 for i in range(0, 9)] +[ i*pi/6 for i in [1,2,4,5,] ]
    knownAngle = numpy.array( l )
    return numpy.concatenate( [-knownAngle[:0:-1], knownAngle ])
knownAngle = setupKnownAngles()

_twopi =  2*numpy.pi
_pi = numpy.pi

def deltaAngle(a1,a2):
    d = a1 - a2 
    return d if d > -_pi else d+_twopi

def closeAngleAbs(a1,a2):
    d = abs(a1 - a2 )
    return min( abs(d-_pi), abs( _twopi - d) , d)

def deltaAngleAbs(a1,a2):
    return abs(in_mPi_pPi(a1 - a2 ))

def in_mPi_pPi(a):
    if(a>_pi): return a-_twopi
    if(a<-_pi): return a+_twopi
    return a
vec_in_mPi_pPi = numpy.vectorize(in_mPi_pPi)
from numpy import sqrt

def D2(p1, p2):
    return ((p1-p2)**2).sum()

def D(p1, p2):
    return sqrt(D2(p1,p2) )

def computeBox(a):
    """returns the bounding box enclosing the array of points a
    in the form (x,y, width, height) """
    xmin , ymin = a[:,0].min(), a[:,1].min()
    xmax , ymax = a[:,0].max(), a[:,1].max()

    return xmin, ymin, xmax-xmin, ymax-ymin

def dirAndLength(p1,p2):
    #l = max(D(p1, p2) ,1e-4)
    l = D(p1,p2)
    uv = (p1-p2)/l
    return l,uv

def length(p1,p2):
    return sqrt( D2(p1,p2) )


# *************************************************************
# debugging 
def void(*l):
    pass
def debug_on(*l):
    #inkex.errormsg(' '.join(str(i) for i in l) ) 
    sys.stderr.write(' '.join(str(i) for i in l) +'\n') 
debug = void

def set_debug_on():
    global debug
    debug = debug_on

#set_debug_on()

# *************************************************************
# Internal Objects
class Path(object):
    """Private representation of a sequence of points.
    A SVG node of type 'path' is splitted in several of these Path objects.
    """
    next = None # next Path in the sequence of path corresponding to a SVG node
    prev = None # previous Path in the sequence of path corresponding to a SVG node
    startIndexSource = 0 # position of first point in the full original array of point
    sourcepoints = None  # the full list of points from which this path is a subset
    
    def __init__(self, points):
        """points an array of points """
        self.points = points
        self.init()

    def init(self):
        self.effectiveNPoints = len(self.points)
        if self.effectiveNPoints>1:
            self.length , self.univ = dirAndLength(self.points[0], self.points[-1])
        else:
            self.length , self.univ = 0, numpy.array([0,0])
        if self.effectiveNPoints>0:
            self.pointN=self.points[-1]
            self.point1=self.points[0]
            
    def isSegment(self):
        return False
    def quality(self):
        return 1000

        

    def dump(self):
        n = len(self.points)
        if n>0:
            return 'path at '+str(self.points[0])+ ' to '+ str(self.points[-1])+'    npoints=%d / %d (eff)'%(n,self.effectiveNPoints)
        else:
            return 'path Void !'
        

    def removeLastPoints(self,n):
        self.points = self.points[:-n]
        self.init()
    def removeFirstPoints(self,n):
        self.points = self.points[n:]
        self.startIndexSource += n
        self.init()

    def costheta(self,seg):
        return self.unitv.dot(seg.unitv)

    def translate(self, tr):
        """Translate this path by tr"""
        #self.c = -self.c -self.a*tr[0] -self.b*tr[0]
        self.points = self.points + tr

    def formatedSegment(self, firstP=False):
        if firstP:            
            segment = [ ['M',[self.point1[0],self.point1[1] ] ],
                        ['L',[self.pointN[0],self.pointN[1] ] ]
                        ]
        else:
            segment = [ ['L',[self.pointN[0],self.pointN[1] ] ] ] 
        return segment

    def setIntersectWithNext(self, next=None):
        pass

    def mergedWithNext(self, newPath=None):
        """ Returns the combination of self and self.next.
        sourcepoints has to be set
        """
        if newPath is None: newPath = Path( self.sourcepoints[self.startIndexSource:self.startIndexSource+len(self.points)+len(self.next.points)] )

        newPath.sourcepoints = self.sourcepoints
        newPath.startIndexSource = self.startIndexSource
        newPath.prev = self.prev
        if self.prev : newPath.prev.next = newPath
        newPath.next = self.next.next
        if newPath.next:
            newPath.next.prev = newPath
        return newPath

# *************************************************************
#     
class Segment(Path):
    """ A segment. Defined by its line equation ax+by+c=0 and the points from orignal paths
    it is ensured that a**2+b**2 = 1
    """
    QUALITYCUT = 0.9
    
    newAngle    = None # temporary angle set during the "parralelization" step
    newLength = None   # temporary lenght set during the "parralelization" step

    # Segment Builders
    @staticmethod
    def from2Points( p1, p2, refPoints = None):
        dirV = p2-p1
        center = 0.5*(p2+p1)
        return Segment.fromCenterAndDir(center, dirV, refPoints)

    @staticmethod
    def fromCenterAndDir( center, dirV, refPoints=None):
        b = dirV[0]
        a = -dirV[1]
        c = - (a*center[0]+b*center[1])

        if refPoints is None:
            refPoints = numpy.array([ center-0.5*dirV, center+0.5*dirV] )
        s = Segment( a, b, c,  refPoints)
        return s

    
    def __init__(self, a,b,c, points, doinit=True):
        """a,b,c: the line parameters.
        points : the array of 2D points represented by this Segment
        doinit : if true will compute additionnal parameters to this Segment (first/last points, unit vector,...)
        """
        self.a = a
        self.b = b
        self.c = c
        
        self.points = points
        d = numpy.sqrt(a**2+b**2)
        if d != 1. :
            self.a /= d
            self.b /= d
            self.c /= d

        if doinit :
            self.init()

    def init(self):
        a,b,c = self.a, self.b, self.c
        x,y = self.points[0]
        self.point1 = numpy.array( [ b*(x*b-y*a) - c*a, a*(y*a-x*b) - c*b ] )
        x,y = self.points[-1]
        self.pointN = numpy.array( [ b*(x*b-y*a) - c*a, a*(y*a-x*b) - c*b ] )
        uv = self.computeDirLength()
        self.distancesToLine =  self.computeDistancesToLine(self.points)
        self.normalv = numpy.array( [ a, b ])

        self.angle = numpy.arccos( uv[0] )*numpy.sign(uv[1] )


    def computeDirLength(self):
        """re-compute and set unit vector and length """
        self.length , uv = dirAndLength(self.pointN, self.point1)
        self.unitv = uv
        return uv

    def isSegment(self):
        return True

    def recomputeEndPoints(self):
        a,b,c = self.a, self.b, self.c
        x,y = self.points[0]
        self.point1 = numpy.array( [ b*(x*b-y*a) - c*a, a*(y*a-x*b) - c*b ] )
        x,y = self.points[-1]
        self.length = numpy.sqrt( D2(self.pointN, self.point1) )

    def projectPoint(self,p):
        """ return the point projection of p onto this segment"""
        a,b,c = self.a, self.b, self.c
        x,y = p
        return numpy.array( [ b*(x*b-y*a) - c*a, a*(y*a-x*b) - c*b ] )        
        

    def intersect(self, seg):
        """Returns the intersection of this line with the line seg"""
        nu, nv = self.normalv, seg.normalv
        u = numpy.array([[-self.c],[-seg.c]])
        doRotation = min(nu.min(),nv.min()) <1e-4
        if doRotation:
            # rotate to avoid numerical issues
            nu = numpy.array(rotMat.dot(nu))[0]
            nv = numpy.array(rotMat.dot(nv))[0]
        debug('  intersection ' ,nu, nv, self.angle, seg.angle)
        m = numpy.matrix( (nu, nv) )        

        i =  (m**-1*u) 
        if doRotation:
            i = unrotMat*i
            
        i=numpy.array( i).swapaxes(0,1)[0]
        return i

    def setIntersectWithNext(self, next=None):
        """Modify self such as self.pointN is the intersection with next segment """
        if next is None:
            next = self.next
        if next and next.isSegment():
            if self.normalv.dot(next.unitv) < 1e-3:
                return
            debug(' Intersect',self, next,  ' from ', self.point1, self.pointN, ' to ' ,next.point1, next.pointN,)
            inter = self.intersect(next)
            debug('  --> ', inter, '  d=', D(self.pointN, inter) )
            next.point1 = inter
            self.pointN = inter
            self.computeDirLength()
            next.computeDirLength()
            
    def computeDistancesToLine(self, points):
        """points: array of points.
        returns the array of distances to this segment"""
        return abs(self.a*points[:,0]+self.b*points[:,1]+self.c)


    def distanceTo(self,point):
        return abs(self.a*point[0]+self.b*point[1]+self.c)        

    def inverse(self):
        """swap all x and y values.  """
        def inv(v):
            v[0], v[1] = v[1] , v[0]
        for v in [self.point1 , self.pointN , self.unitv, self.normalv]:
            inv(v)

        self.points = numpy.roll(self.points,1,axis=1)
        self.a, self.b = self.b, self.a
        self.angle = numpy.arccos( self.unitv[0] )*numpy.sign(self.unitv[1] )
        return

    def dumpShort(self):
        return 'seg  '+str(self.startIndexSource)+'  '+str(self.point1 )+'to '+str(self.pointN)+ ' npoints=%d | angle,offset=(%.2f,%.2f )'%(len(self.points),self.angle, self.c)+'  ',self.normalv

    def dump(self):
        v = self.variance()
        n = len(self.points)
        return 'seg  '+str(self.point1 )+' , '+str(self.pointN)+ '  v/l=%.2f / %.2f = %.2f  r*sqrt(n)=%.2f  npoints=%d | angle,offset=(%.2f,%.2f )'%(v, self.length, v/self.length,v/self.length*numpy.sqrt(n) ,n  , self.angle, self.c)
        
    def variance(self):
        d = self.distancesToLine
        return numpy.sqrt( (d**2).sum()/len(d) )

    def quality(self):
        n = len(self.points)
        return min(self.variance()/self.length*numpy.sqrt(n) , 1000)

    def formatedSegment(self, firstP=False):
        if firstP:            
            segment = [ ['M',[self.point1[0],self.point1[1] ] ],
                        ['L',[self.pointN[0],self.pointN[1] ] ]
                        ]
        else:
            segment = [ ['L',[self.pointN[0],self.pointN[1] ] ] ]
        #debug("Segment, format : ", segment)
        return segment
        
    def replaceInList(self, startPos, fullList):
        code0 = fullList[startPos][0]
        segment = [ [code0,[self.point1[0],self.point1[1] ] ],
                     ['L',[self.pointN[0],self.pointN[1] ] ]
                    ]
        l = fullList[:startPos]+segment+fullList[startPos+len(self.points):]
        return l




    def mergedWithNext(self):
        """ Returns the combination of self and self.next.
        sourcepoints has to be set
        """
        #spoints = self.sourcepoints[self.startIndexSource:self.startIndexSource+len(self.points)+len(self.next.points)]
        spoints = numpy.concatenate([self.points,self.next.points])
        newSeg = fitSingleSegment(spoints)
        
        newSeg = Path.mergedWithNext(self, newSeg)
        if newSeg.next:
            if newSeg.next.isSegment():
                newSeg.deltaAngle = abs(newSeg.angle - newSeg.next.angle)
            else:
                newSeg.deltaAngle = 10
        return newSeg

    
    def refit(self):
        xmin,ymin,w,h = computeBox(self.points)
        inverse = w<h
        if inverse:
            self.points = numpy.roll(self.points,1,axis=1) 

        s , o = regLin(self.points, returnOnlyPars=True)
        self.a=s
        self.b=-1
        self.c = o
        self.init()
        if inverse:
            self.inverse()

    def barycenter(self):
        return 0.5*(self.point1+self.pointN)
    def barycenter_origin(self):
        return self.points.sum(axis=0)/len(self.points)

    def box_origin(self):
        return computeBox(self.points)


    def translate(self, tr):
        """Translate this segment by tr """
        self.c = -self.c -self.a*tr[0] -self.b*tr[0]
        self.pointN = self.pointN+tr
        self.point1 = self.point1+tr
        
    def adjustToNewAngle(self):        
        """reset all parameters so that angle is change to self.newAngle """
        self.a,self.b,self.c = parametersFromPointAngle( 0.5*(self.point1+self.pointN), self.newAngle)
        self.angle = self.newAngle
        self.normalv = numpy.array( [ self.a, self.b ])
        self.unitv = numpy.array( [ self.b, -self.a ])
        if abs(self.angle) > numpy.pi/2 :
            if self.b > 0: self.unitv *= -1
        elif self.b<0 : self.unitv  *= -1

        self.point1 = self.projectPoint(self.point1) # reset point1 
        if self.next is None or not self.next.isSegment():
                # move the last point (no intersect with next)

                pN = self.projectPoint(self.pointN)
                dirN = pN - self.point1                
                lN = length(pN, self.point1)
                self.pointN = dirN/lN*self.length + self.point1
                #print ' ... adjusting last seg angle ',p.dump() , ' normalv=', p.normalv, 'unitv ', p.unitv
        else:
            self.setIntersectWithNext()


    def adjustToNewDistance(self):
        self.pointN = self.newLength* self.unitv + self.point1
        self.length = self.newLength
        #pass

    def tempLength(self):
        if self.newLength : return self.newLength
        else : return self.length
    def tempAngle(self):
        if self.newAngle: return self.newAngle
        return self.angle




# *************************************************************
# *************************************************************
# Groups of Path
#
class PathGroup(object):
    """A group of Path representing one SVG node.
     - a list of Path
     - a list of SVG commands describe the full node (=SVG path element)
     - a reference to the inkscape node object
     
    """
    def __init__(self, listOfPaths, refSVGPathList, refNode=None):
        self.refNode = refNode
        self.listOfPaths = listOfPaths
        self.refSVGPathList = refSVGPathList
        
    def addToNode(self, node):
        newList = reformatList( self.refSVGPathList, self.listOfPaths)        
        ele = addPath( newList , node)
        debug("PathGroup ", newList)
        return ele

class TangeantEnvelop(PathGroup):
    """Specialization where the Path objects are all Segments and represent tangeants to a curve """
    def addToNode(self, node):
        newList = [ ]
        for s in self.listOfPaths:
            newList += s.formatedSegment(firstP=True)
        debug("TangeantEnvelop ", newList)
        ele = addPath( newList , node)
        return ele

class Circle(PathGroup):
    """Specialization where the list of Path objects
    is to be replaced by a Circle specified by a center and a radius.

    If an other radius 'rmax' is given than this represent an ellipse.
    """

    def __init__(self, center, rad, listOfPaths, refNode=None, rmax=None, angle=0.):
        self.listOfPaths = listOfPaths
        self.refNode = refNode
        self.center = center
        self.radius = rad
        if rmax:
            self.type ='ellipse'
        else:
            self.type = 'circle'
        self.rmax = rmax
        self.angle = angle
        
    def addToNode(self, refnode):
        """Add a node in the xml structure corresponding to this rect
        refnode : xml node used as a reference, new point will be inserted a same level"""
        ele = inkex.etree.Element('{http://www.w3.org/2000/svg}'+self.type)

        ele.set('cx',str(self.center[0]))
        ele.set('cy',str(self.center[1]))
        if self.rmax:
            ele.set('ry',str(self.radius))
            ele.set('rx',str(self.rmax))
            ele.set('transform', 'rotate(%3.2f,%f,%f)'%(numpy.degrees(self.angle),self.center[0],self.center[1]))
        else:
            ele.set('r',str(self.radius))
        refnode.xpath('..')[0].append(ele)
        return ele

    
class Rectangle(PathGroup):
    """Specialization where the list of Path objects
    is to be replaced by a Rectangle specified by a center and size (w,h) and a rotation angle.

    """
    def __init__(self, center, size, angle, listOfPaths, refNode=None):
        self.listOfPaths = listOfPaths
        self.refNode = refNode
        self.center = center
        self.size = size
        self.bbox = size
        self.angle = angle
        pos = self.center - numpy.array( size )/2
        if angle != 0 :
            cosa = numpy.cos(angle)
            sina = numpy.sin(angle)            
            self.rotMat = numpy.matrix( [ [ cosa, sina], [-sina, cosa] ] )
            self.rotMatstr = 'matrix(%1.7f,%1.7f,%1.7f,%1.7f,0,0)'%(cosa, sina, -sina, cosa)


            #debug(' !!!!! Rotated rectangle !!', self.size, self.bbox,  ' angles ', a, self.angle ,' center',self.center)
        else :
            self.rotMatstr = None
        self.pos = pos
        debug(' !!!!! Rectangle !!', self.size, self.bbox,  ' angles ', self.angle ,' center',self.center)

    def addToNode(self, refnode):
        """Add a node in the xml structure corresponding to this rect
        refnode : xml node used as a reference, new point will be inserted a same level"""
        ele = inkex.etree.Element('{http://www.w3.org/2000/svg}rect')
        self.fill(ele)
        refnode.xpath('..')[0].append(ele)
        return ele
        
    def fill(self,ele):
        w, h = self.size
        ele.set('width',str(w))
        ele.set('height',str(h))
        w, h = self.bbox
        ele.set('x',str(self.pos[0]))
        ele.set('y',str(self.pos[1]))
        if self.rotMatstr:
            ele.set('transform', 'rotate(%3.2f,%f,%f)'%(numpy.degrees(self.angle),self.center[0],self.center[1]))
            #ele.set('transform', self.rotMatstr)

    @staticmethod
    def isRectangle( pathGroup):
        """Check if the segments in pathGroups can form a rectangle.
        Returns a Rectangle or None"""
        if isinstance(pathGroup, Circle ): return None
        segmentList = [p for p in pathGroup.listOfPaths if p.isSegment() ]#or p.effectiveNPoints >0]
        if len(segmentList) != 4:
            debug( 'rectangle Failed at length ', len(segmentList))
            return None
        a,b,c,d = segmentList

        if length(a.point1, d.pointN)> 0.2*(a.length+d.length)*0.5:
            debug('rectangle test failed closing ', length(a.point1, d.pointN), a.length, d.length)
            return None
        
        Aac , Abd = closeAngleAbs(a.angle,c.angle), closeAngleAbs(b.angle , d.angle)
        if  min(Aac,Abd) > 0.01 or max(Aac, Abd) >0.17 :
            debug( 'rectangle Failed at angles', Aac, Abd)
            return None
        notsimilarL = lambda d1,d2: abs(d1-d2)>0.20*min(d1,d2)

        pi , twopi = numpy.pi,2*numpy.pi
        angles = numpy.array( [p.angle   for p in segmentList] )
        minAngleInd = numpy.argmin( numpy.minimum( abs(angles), abs( abs(angles)-pi), abs( abs(angles)-twopi) ) )
        rotAngle = angles[minAngleInd]
        width = (segmentList[minAngleInd].length + segmentList[(minAngleInd+2)%4].length)*0.5
        height = (segmentList[(minAngleInd+1)%4].length + segmentList[(minAngleInd+3)%4].length)*0.5
        # set rectangle center as the bbox center
        x,y,w,h = computeBox( numpy.concatenate( [ p.points for p in segmentList]) )
        r = Rectangle( numpy.array( [x+w/2, y+h/2]), (width, height), rotAngle, pathGroup.listOfPaths, pathGroup.refNode)
        
        debug( ' found a rectangle !! ', a.length, b.length, c.length, d.length )
        return r


# *************************************************************
# Object manipulation functions

def toRemarkableShape( group ):
    """Test if PathGroup instance group looks like a remarkable shape (ex: Rectangle).
    if so returns a new shape instance else returns group unchanged"""
    r = Rectangle.isRectangle( group )
    if r : return r
    return group


def resetPrevNextSegment(segs):
    for i, seg in enumerate(segs[:-1]):
        s = segs[i+1]
        seg.next = s
        s.prev = seg           
    return segs

def fitSingleSegment(a):
    xmin,ymin,w,h = computeBox(a)
    inverse = w<h
    if inverse:
        a = numpy.roll(a,1,axis=1)

    seg = regLin(a)
    if inverse:
        seg.inverse()
        #a = numpy.roll(a,1,axis=0)
    return seg
        
def regLin(a , returnOnlyPars=False):
    """perform a linear regression on 2dim array a. Creates a segment object in return """
    sumX = a[:,0].sum()
    sumY = a[:,1].sum()
    sumXY = (a[:,1]*a[:,0]).sum()
    a2 = a*a
    sumX2 = a2[:,0].sum()
    sumY2 = a2[:,1].sum()
    N = a.shape[0]

    pa = (N*sumXY - sumX*sumY)/ ( N*sumX2 - sumX*sumX)
    pb = (sumY - pa*sumX) /N
    if returnOnlyPars:
        return pa,-1, pb
    return Segment(pa, -1, pb, a)


def smoothArray(a, n=2):
    count = numpy.zeros(a.shape)
    smootha = numpy.array(a)
    for i in range(n):
        count[i]=n+i+1
        count[-i-1] = n+i+1
    count[n:-n] = n+n+1
    #debug('smooth ', len(smooth[:-2]) [)
    for i in range(1,n+1):
        smootha[:-i] += a[i:]
        smootha[i:]  += a[:-i]
    #print smootha
    return smootha/count

def buildTangeants( points , averaged=True):
    tangeants = numpy.zeros( (len(points),2) )
    i=1
    tangeants[:-i] += points[i:] - points[:-i] 
    tangeants[i:]  += points[i:] - points[:-i] 

    tangeants *= 0.5
    tangeants[0] *=2
    tangeants[-1] *=2

    ## debug('points ', points)
    ## debug('buildTangeants --> ', tangeants )
    
    if averaged:
        # average over neighbours
        avTan = numpy.array(tangeants)
        avTan[:-1] += tangeants[1:]
        avTan[1:]  += tangeants[:-1]
        avTan *= 1./3
        avTan[0] *=1.5
        avTan[-1] *=1.5

    return avTan


def clusterAngles(array, dAng=0.15):
    """ array : flat array of angles
    returns [ ...,  (indi_0, indi_1),...] where each tuple are indices of cluster i
    """
    N = len(array)

    closebyAng = numpy.zeros( (N,4) , dtype=int)

    for i,a in enumerate(array):
        cb = closebyAng[i]
        cb[0] =i
        cb[2]=i
        cb[3]=i
        c=i-1
        # find number of angles within dAng in nearby positions
        while c>-1: # indices below i
            d=closeAngleAbs(a,array[c])
            if d>dAng:
                break
            cb[1]+=1                
            cb[2]=c
            c-=1
        c=i+1
        while c<N-1:# indices above i
            d=closeAngleAbs(a,array[c])
            if d>dAng:
                break
            cb[1]+=1                
            cb[3]=c
            c+=1
    closebyAng= closebyAng[numpy.argsort(closebyAng[:,1]) ]

    clusteredPos = numpy.zeros(N, dtype=int)
    clusters = []
    for cb in reversed(closebyAng):
        if clusteredPos[cb[0]]==1:
            continue
        # try to build a cluster
        minI = cb[2]
        while clusteredPos[minI]==1:
            minI+=1
        maxI = cb[3]
        while clusteredPos[maxI]==1:
            maxI-=1
        for i in range(minI, maxI+1):
            clusteredPos[i] = 1
        clusters.append( (minI, maxI) )

    return clusters
        
    
                

def adjustAllAngles(paths):
    for p in paths:
        if p.isSegment() and p.newAngle is not None:
            p.adjustToNewAngle()
    # next translate to fit end points
    tr = numpy.zeros(2)
    for p in paths[1:]:
        if p.isSegment() and p.prev.isSegment():
            tr = p.prev.pointN - p.point1
        debug(' translating ',p,' prev is', p.prev, '  ',tr, )
        p.translate(tr)

def adjustAllDistances(paths):
    for p in paths:
        if p.isSegment() and  p.newLength is not None:                
            p.adjustToNewDistance()
    # next translate to fit end points
    tr = numpy.zeros(2)
    for p in paths[1:]:
        if p.isSegment() and p.prev.isSegment():
            tr = p.prev.pointN - p.point1
        p.translate(tr)


def mergeConsecutiveParralels(segments):
    ignoreNext=False
    newList=[]
    for s in segments:
        if ignoreNext:
            ignoreNext=False
            continue
        if not s.isSegment():
            newList.append(s)
            continue
        if s.next is None:
            newList.append(s)
            continue
        if not s.next.isSegment():
            newList.append(s)
            continue
        d = closeAngleAbs(s.angle ,s.next.angle)
        if d<0.001:
            debug("merging ", s.angle ,s.next.angle )
            snew = s.mergedWithNext()
            ignoreNext=True
            newList.append(snew)
        else:
            newList.append(s)
    if len(segments)>len(newList):
        debug("merged parallel ", segments , '-->', newList)
    return newList

def parametersFromPointAngle(point, angle):
    unitv = numpy.array([ numpy.cos(angle), numpy.sin(angle) ])
    ortangle = angle+numpy.pi/2
    normal = numpy.array([ numpy.cos(ortangle), numpy.sin(ortangle) ])
    genOffset = -normal.dot(point)
    a, b = normal
    return a, b , genOffset
    
def segmentFromPointsAngle(points, angle):
    barycenter = points.sum(axis=0)/len(points)
    a, b , offset = parametersFromPointAngle(barycenter, angle)
    return Segment(a,b,offset,points)


def addPath(newList, refnode):
    """Add a node in the xml structure corresponding to the content of newList
    newList : list of Segment or Path
    refnode : xml node used as a reference, new point will be inserted a same level"""
    ele = inkex.etree.Element('{http://www.w3.org/2000/svg}path')
    ele.set('d', simplepath.formatPath(newList))
    refnode.xpath('..')[0].append(ele)
    return ele

def reformatList( refSVGPathList, paths):
    """ Returns a SVG paths list in (same format as simplepath.parsePath) from a list of our Path 'paths'.
     - Segments in paths are added in the new list
     - simple Path are retrieved from the original refSVGPathList and put in the new list (thus preserving original bezier curves)
    """
    newList = []
    first = True
    for  seg in paths:
        if seg.isSegment():
            fseg = seg.formatedSegment(first)
            newList +=fseg
        else:
            # we re-use the original path to make use of bezier curves contained in it.
            # (our Path doesn't hold this info).
            if first:
                if seg.effectiveNPoints==0: # the 1st path has been removed. ignore it.
                    continue
                pos = seg.startIndexSource
            else:
                pos = seg.startIndexSource+1
            newList += refSVGPathList[pos:pos+seg.effectiveNPoints]
        first = False
    return newList
        






# *************************************************************
# The inkscape extension
# *************************************************************

class ShapeReco(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("--title")
        self.OptionParser.add_option("-k", "--keepOrigin",
                        action="store", type="inkbool", 
                        dest="keepOrigin", default=False,
                        help="Do not replace path")
        self.defPathName = [] # debugging only



    def effect(self):

        rej='{http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd}type'
        paths = []
        for id, node in self.selected.iteritems():
            if node.tag == '{http://www.w3.org/2000/svg}path' and rej not in node.keys():                
                paths.append(node)
        if paths == []:
            paths=[self.getElementById(n) for n in self.defPathName]
        else:
            p = paths[0]

        self.extractShapes(paths)
        
        #inkex.errormsg('  options '+str(self.options.keepOrigin))

    def removeSmallEdge(self, paths, wTot,hTot):
        """Remove small Path objects which stand between 2 Segments (or at the ends of the sequence).
        Small means the bbox of the path is less then 5% of the mean of the 2 segments."""
        if len(paths)<2:
            return
        def getdiag(points):
            xmin,ymin,w,h = computeBox(p.points)
            return sqrt(w**2+h**2)
        removeSeg=[]
        def remove(p):
            removeSeg.append(p)
            if p.next : p.next.prev = p.prev
            if p.prev: p.prev.next = p.next
            p.effectiveNPoints =0
            debug('      --> remove !', p, p.length , len(p.points))
        for p in paths:
            if len(p.points)==0 :
                remove(p)
                continue
            # select only path between 2 segments
            next, prev = p.next, p.prev
            if next is None: next = prev
            if prev is None: prev = next
            if not next.isSegment() or not prev.isSegment() : continue
            diag = getdiag(p.points)

            debug(p, p.pointN, ' removing edge  diag = ', diag, p.length,  '  l=',next.length+prev.length)
            debug( '    ---> ',prev, next)
            if diag > (next.length+prev.length)*0.1 : continue
            if diag > 0.2*wTot or diag > 0.2*hTot: continue # avoid removing if significant on total
            # Avoid removing connecting path between 2 long rectangle side
            if diag > D(next.pointN, prev.point1)*0.5: continue
            remove(p)

            if next != prev:
                prev.setIntersectWithNext(next)
        debug('removed Segments ', removeSeg)
        for p in removeSeg:
            paths.remove(p)


    def distancesRatios(self, points):
        """ DEPRECATED ??
        Finds segment-like portions in a sequence of points.
        input : points is a (N,2) shaped array
        output : a list of Segment or Path objects

        The idea used here is for aligned points D(p_i,p_n) = Sum(D_i,D_i-1) = S
        So start from point 1, aggregate all following points as long as S is close to D, when this is not true or D(p_1,P_n)<D(p_1,p_n-1)
        start a new segment.
        All points in segments of length less than 4 are simply put into generic Path.
        """
        #points = self.origPoints
        Np = len(points)
        windowS= min( max( Np/10, 2) , 10 )

        presegs = []
        i=0
        def pdistance(i,j):
            return numpy.sqrt( ((points[i]-points[j])**2).sum() )

        # aggregates points into candidate segments (presegs)
        while i < Np-1:
            # start a new preseg with point i
            next = i
            #sumD = dmat[i,next]
            sumD =  0
            lastD = sumD
            curD = sumD
            r = 1
            # loop on next points as long as they're compatible with a segment.
            while r < 1.15 and curD >= lastD and next< Np-1:
                next += 1
                #curD = dmat[next-1,neaxt]
                lastD =curD
                curD = max(pdistance(i,next) , 1e-12)
                sumD += pdistance(next-1,next)
                r = sumD/curD
            # mark the start and end position of the pre-segment
            next -=1            
            presegs.append( (i,next) )
            #print (i,next) ,'  r=',r
            i = next +1
        presegs[ -1 ] = ( presegs[-1][0],presegs[-1][1]+1 ) #append the last point

        #print 'presegs ', presegs
        #  convert to Segment:
        minNpointInSeg = min(4,max(Np/20,2))
        segs = []
        for (b,e) in presegs:
            if e+1-b<4: # less than 4 points, not a Segment
                seg = Path(points[b:e+1])
            else: # comute the segment (fit a line equation)
                seg = fitSingleSegment(points[b:e+1])
                debug('  seg at ',(b,e),seg.dump())
                if seg.quality() > Segment.QUALITYCUT : # fit failed. Not a segment
                    seg = Path(points[b:e+1])
            seg.startIndexSource = b
            seg.sourcepoints = points
            segs.append( seg )
        # assign next and previous :
        for p,pnext in zip(segs[:-1] , segs[1:]):
            p.next = pnext
            pnext.prev = p
        if len(segs)>1 :segs[-1].prev = segs[-2]

            

        # merge consecutive Path objects 
        mergedpaths = []
        nseg = len(segs)
        def mergeConsecutivePath(p, pathStart=0):
            nexts = 0
            nextp = p.next
            if  p.isSegment() :
                yield p
                if p.next is None:
                    return
                nexts = p.next.startIndexSource
            elif p.next is None or p.next.isSegment() :
                newp = Path( points[pathStart:p.startIndexSource+p.effectiveNPoints] )                
                newp.startIndexSource = pathStart
                yield newp
                if p.next is None:
                    return
                nexts = p.next.startIndexSource
            else:
                nexts = pathStart
            #print ' merging cont  ... ',nextp.dump(), nexts
            for op in mergeConsecutivePath(nextp, nexts):
                yield op

        segs =[ p for p in mergeConsecutivePath(segs[0],0) ]    
        
        # merge consecutive segments if they have similar angles
        # compute delta angles
        for p in segs:
            if p.isSegment() and p.next and p.next.isSegment():                
                p.deltaAngle = abs(p.angle - p.next.angle)
            else:
                p.deltaAngle = 10
        pmin = min( segs, key = lambda x:x.deltaAngle)
        while pmin.deltaAngle < 0.13962: # degree
            spoints = points[pmin.startIndexSource:pmin.next.startIndexSource+len(pmin.next.points)]
            newSeg  = pmin.mergedWithNext(  ) 

            segs.remove(pmin)
            segs.remove(pmin.next)
            segs.append(newSeg)
            pmin = min( segs, key = lambda x:x.deltaAngle)
        
        return segs
            
            
    
    def prepareParrallelize(self,segs):
        """Group Segment by their angles (segments are grouped together if their deltAangle is within 0.02 rad)
        The angles of segments in a group are then set to the  angle of sum(unit_vector) where the sum runs on the group

        segs ; list of segments
        """
        
        closeAngleSet = {}
        for (i, seg1) in enumerate(segs):
                for seg2 in segs[i+1:]:
                    #print i, closeAngleAbs(seg1.angle , seg2.angle), abs(seg1.costheta(seg2)), '__',seg1.angle , seg2.angle
                    if abs(seg1.costheta(seg2)) > 0.99: # 0.99 ~ 8 deg
                        l1=closeAngleSet.get(seg1,set())
                        l1.add(seg1)
                        l1.add(seg2)
                        closeAngleSet[seg2] = l1
                        debug('parralelize ',seg1.dump() , seg2.dump())
        debug('parralelize : ', len(closeAngleSet))
        sign = numpy.sign
        for segs in closeAngleSet.values():            
            lsegs = list(segs)
            s0 = lsegs[0]
            signs = [sign(s0.costheta(seg )) for seg in lsegs]
            m = sum( si*s.unitv for (si,s) in zip(signs,lsegs) )
            m /= numpy.sqrt(m.dot(m))

            meanAngle0 = numpy.arccos( m[0] )*numpy.sign(m[1] )
            meanAngle1 = -numpy.arccos( -m[0] )*numpy.sign(m[1] )
            for i,s in enumerate(segs):
                #s.newAngle = meanAngle0 if signs[i]>=0 else meanAngle1
                s.newAngle =  meanAngle1 if deltaAngleAbs(meanAngle1, s.angle) < deltaAngleAbs(meanAngle0, s.angle) else meanAngle0
                debug( ' mean angles : ', s.angle , ' --> ', s.newAngle, ' XXXX ',meanAngle0 , meanAngle1)


    def prepareDistanceEqualization(self,segs, relDelta=0.1):
        """ Input segments are grouped according to their length  :
          - for each length L, find all other lengths within L*relDelta. of L.
          - Find the larger of such subgroup.
          - repeat the procedure on remaining lengths until none is left.
        Each length in a group is set to the mean length of the group

        segs : a list of segments
        relDelta : float, minimum relative distance.
        """
        lengths= sorted(  (x.tempLength() ,i) for i,x in enumerate(segs)  )
        debug( '___  lengths ',lengths)

        def findgroups(ls, startPos=0 ):
            maxgroupdelta = 0
            maxgrouppos = ()
            #print '  ccc find groups in ', ls
            for i,d in enumerate(ls):
                thisD = d[0]
                delta = thisD*relDelta
                first, last = i-1,i+1
                while first>-1:
                    if thisD-ls[first][0] > delta:
                        break
                    first-=1
                first +=1
                while last<len(ls):
                    if ls[last][0]-thisD > delta:
                        break
                    last+=1
                last -= 1
                #print i,'     --> ',last , first, delta
                if last-first > maxgroupdelta:
                    maxgroupdelta = last-first
                    maxgrouppos = (first,last)
            #print ' ====> maxgroupdelta ',maxgroupdelta, ' // ',maxgrouppos
            l = []
            if maxgroupdelta >0 :
                l += [ (maxgrouppos[0]+startPos, maxgrouppos[1]+startPos) ]
                debug( '   dist eq , ', maxgrouppos )
                l+= findgroups( ls[:maxgrouppos[0]], startPos)
                l+=findgroups( ls[maxgrouppos[1]+1:], startPos+maxgrouppos[1]+1)
            return l

        sameDistGroups = findgroups(lengths)

        if len(sameDistGroups)==1 :
            # deal with special case with low num of segments
            # --> don't let a single segment alone
            nToChange = sameDistGroups[0][1] - sameDistGroups[0][0]+1
            debug(' __ only 1 group ', len(lengths), nToChange)
            if len(lengths)-nToChange ==1:
                sameDistGroups[0] = (0,len(lengths)-1 )
                
        for f,l in sameDistGroups:
            dmean = sum( d for (d,i) in lengths[f:l+1] )/ (l-f+1)
            for d,i in lengths[f:l+1]:
                segs[i].newLength = dmean
                debug( i,' set newLength ',dmean, segs[i].length, segs[i].dumpShort())
        

    def adjustToKnownAngle(self, paths):
        """ Check current angle against remarkable angles. If close enough, change it
        paths : a list of segments"""
        for seg in paths:
            a = seg.tempAngle()
            i = (abs(knownAngle - a )).argmin()
            if abs(knownAngle[i] - a) < 0.08:
                debug( '  Known angle ', seg.tempAngle(),'  -> ', knownAngle[i]) 
                seg.newAngle = knownAngle[i]

        

    def checkForCircle(self, points, tangeants):
        """Determine if the points and their tangeants represent a circle

        Method : we consider angle of tangeant as function of lenght on path.
        For circles these are : angle = c1 x lenght + c0.
        We thus fit a line in the (angle, lenght) plane.
        If the path is made of line, angle will be constants for large section,
        thus the distance to fitted line will be bad for a large num of points.
        we use the quantiles above 80% to classify.
        
        """
        if len(points)<6:
            return False, 0
        xmin,ymin, w, h = computeBox( points)
        d2 = D2(points[0], points[-1] )
        diag2=(w*w+h*h)
        if d2 > diag2*(0.1*0.1):
            return False,0
        
        diag = sqrt(diag2)*0.5
        norms = numpy.sqrt(numpy.sum( tangeants**2, 1 ))
        #print  'norms ', norms
        #debug( 'tangeants ', tangeants)
        angles = numpy.arccos( tangeants[:,0] /norms ) *numpy.sign( tangeants[:,1] )
        #debug( 'angle = ', repr(angles))
        N = len(angles)
        angles = smoothArray(angles, n=max(N/20,2) )
        deltas =  points[1:] - points[:-1] 
        deltasD = numpy.concatenate([ [0.], numpy.sqrt(numpy.sum( deltas**2, 1 )) / diag] )


        imin = numpy.argmin(angles)
        imax = numpy.argmax(angles)
        if imin<imax:
            angles[:imin+1] +=numpy.pi*2
        else:
            angles[imin:] +=numpy.pi*2
        debug(' imin ',imin)
        angles = numpy.roll(angles, -imin)
        deltasD = numpy.roll(deltasD, -imin)
        deltasD = deltasD.cumsum()
        
        angles_vs_d = numpy.stack([ deltasD, angles] ,1)
        seg =regLin( angles_vs_d )
        dist_to_fit = numpy.abs( angles + (seg.a*angles_vs_d[:,0]+seg.c)/seg.b )
        dist_to_fit.sort()
        debug('fit q',seg.quality() , '  d90=', dist_to_fit[int(N*0.9)], ' d80=',dist_to_fit[int(N*0.8)])
        ## debug( 'angle = ', repr(angles))
        ## debug( 'deltasD=', repr(deltasD))
        ## debug( 'dist_to_fit = ', repr(dist_to_fit))        

        if dist_to_fit[int(N*0.8)] > 1.3: #empirical
            return False,0
        radius = points - numpy.array([xmin+w*0.5,ymin+h*0.5])
        radius_n = numpy.sqrt(numpy.sum( radius**2, 1 )) # normalize

        rmin = min(radius_n)
        maxi = numpy.argmax(radius_n)
        rmax = radius_n[maxi]
        anglemax = numpy.arccos( radius[maxi][0]/rmax)*numpy.sign(radius[maxi][1])
        return True, (xmin+w*0.5,ymin+h*0.5, rmin, rmax, anglemax)




    def tangeantEnvelop(self, svgCommandsList, refNode):
        a, svgCommandsList = toArray(svgCommandsList)
        tangeants = buildTangeants(a)

        newSegs = [ Segment.fromCenterAndDir( p, t ) for (p,t) in zip(a,tangeants) ]
        debug("build envelop ", newSegs[0].point1 , newSegs[0].pointN)
        clustersInd = clusterAngles( [s.angle for s in newSegs] )
        debug("build envelop cluster:  ", clustersInd)

        return TangeantEnvelop( newSegs, svgCommandsList, refNode)


    def segsFromTangeants(self,svgCommandsList, refNode):
        """Finds segments part in a list of points represented by svgCommandsList.

        The method is to build the (averaged) tangeant vectors to the curve.
        Aligned points will have tangeant with similar angle, so we cluster consecutive angles together
        to define segments.
        Then we extend segments to connected points not already part of other segments.
        Then we merge consecutive segments with similar angles.
        
        """
        sourcepoints, svgCommandsList = toArray(svgCommandsList)
        tangeants = buildTangeants(sourcepoints)

        # Check if circle -----------------------
        isCircle, res = self.checkForCircle( sourcepoints, tangeants)        
        if isCircle:
            x,y,rmin, rmax,angle = res
            if rmin/rmax>0.8:
                circ = Circle((x,y),0.5*(rmin+rmax), [], refNode )
            else:
                circ = Circle((x,y),rmin, [], refNode, rmax=rmax, angle=angle)
            circ.points = sourcepoints
            return circ
        debug("Is Circle = ", isCircle )

        # global quantities :
        x,y,wTot,hTot = computeBox(sourcepoints)
        minDim = min(wTot, hTot)
        forceClose = False#minDim*0.1 > D(a[0],a[-1])
        debug('forceClose ', forceClose)


        # cluster points by angle of their tangeants -------------
        tgSegs = [ Segment.fromCenterAndDir( p, t ) for (p,t) in zip(sourcepoints,tangeants) ]
        clustersInd = clusterAngles( [s.angle for s in tgSegs] )
        clustersInd.sort( )
        debug("build envelop cluster:  ", clustersInd)

        # build Segments from clusters 
        newSegs = []
        for imin, imax in clustersInd:
            if imin+1< imax: # consider clusters with more than 3 points
                seg = fitSingleSegment(sourcepoints[imin:imax+1])
            elif imin+1==imax: # 2 point path : we build a segment
                seg = Segment.from2Points(sourcepoints[imin], sourcepoints[imax] , sourcepoints[imin:imax+1])
            else:
                seg = Path( sourcepoints[imin:imax+1] )
            seg.startIndexSource = imin
            seg.sourcepoints = sourcepoints
            newSegs.append( seg )
        resetPrevNextSegment( newSegs )
        debug(newSegs)
        # -----------------------

        # Merge consecutive paths
        updatedSegs=[]
        def toMerge(p):
            l=[p]
            setattr(p, 'merged', True)
            if p.next and not p.next.isSegment():
                l += toMerge(p.next)
            return l
        
        for i,seg in enumerate(newSegs[:-1]):
            if seg.isSegment():
                updatedSegs.append( seg)                
                continue
            if hasattr(seg,'merged'): continue
            mergeList = toMerge(seg)
            debug('merging ', mergeList)
            for p in mergeList:
                debug(' ----> ', p.points)
            p = Path(numpy.concatenate([ p.points for p in mergeList]) )
            p.startIndexSource = mergeList[0].startIndexSource
            debug('merged == ', p.points)
            updatedSegs.append(p)

        if not hasattr(newSegs[-1],'merged'): updatedSegs.append( newSegs[-1]) 
        #debug("unmerged path", newSegs)
        debug("merged path", updatedSegs)
        newSegs = resetPrevNextSegment( updatedSegs )


        # Extend segments -----------------------------------
        class SegmentExtender:
            def nextPaths(self,seg):
                pL = []
                p = self.getNext(seg) # prev or next
                while p :
                    if p.isSegment(): break
                    pL.append(p)
                    p = self.getNext(p)
                if pL==[]:
                    return [], []
                i0=pL[0].startIndexSource
                iN = pL[-1].startIndexSource+len(pL[-1].points)+1
                return pL, seg.sourcepoints[i0:iN]
            
            def extend(self,seg):
                nextPathL, pointsToTest = self.nextPaths(seg)
                debug('extend ',self.extDir, seg , nextPathL)
                if nextPathL==[]: return seg
                distancesToLine =abs(seg.a*pointsToTest[:,0]+seg.b*pointsToTest[:,1]+seg.c)
                mergeD = seg.length*0.03
                pointsToFit, addedPoints = self.pointsToFit(seg,pointsToTest, distancesToLine, mergeD)
                if len(pointsToFit)==0:
                    return seg
                newseg = fitSingleSegment(pointsToFit)
                if newseg.quality()>0.5: # fit failed
                    return seg
                debug( '  EXTENDING ! ', len(seg.points), len(addedPoints) )
                self.removePath(seg, newseg, nextPathL, addedPoints )
                newseg.points = pointsToFit
                seg.mergedObj= newseg
                newseg.sourcepoints = seg.sourcepoints
                
                return newseg

            @staticmethod
            def extendSegments(newSegs):
                fwdExt = FwdExtender()
                bwdExt = BwdExtender()
                for seg in newSegs:
                    seg.mergedObj = seg
                for seg in sorted(newSegs, key = lambda s : s.length, reverse=True):
                    if seg.isSegment():
                        newseg=fwdExt.extend(seg)
                        seg.mergedObj = bwdExt.extend(newseg)
                updatedSegs=[seg.mergedObj for seg in newSegs if seg.mergedObj]
                return updatedSegs


        class FwdExtender(SegmentExtender):
            extDir='Fwd'
            def getNext(self, seg):
                return seg.next
            def pointsToFit(self, seg, pointsToTest, distancesToLine, mergeD):
                for i,d in reversed(list(enumerate(distancesToLine))):
                    if d<mergeD: break
                addedPoints = pointsToTest[i:]
                debug( ' ++ pointsToFit ' , mergeD, i ,len(pointsToTest), addedPoints , seg.points )
                return  numpy.concatenate([seg.points, addedPoints]), addedPoints
            def removePath(self, seg, newseg, nextPathL, pointsToFit):
                npoints = len(pointsToFit)
                acc=0
                newseg.prev = seg.prev
                for p in nextPathL:
                    if (acc+len(p.points))<=npoints:
                        p.mergedObj = None
                        acc += len(p.points)
                    else:
                        newseg.next = p
                        p.points = p.points[:(npoints-acc-len(p.points))]
                        break

        class BwdExtender(SegmentExtender):
            extDir='Bwd'
            def getNext(self, seg):
                return seg.prev
            def pointsToFit(self, seg, pointsToTest, distancesToLine, mergeD):
                for i,d in enumerate(distancesToLine):
                    if d<mergeD: break
                addedPoints = pointsToTest[:i+1]
                debug( ' ++ pointsToFit ' , mergeD, i ,len(pointsToTest), addedPoints , seg.points )
                return  numpy.concatenate([addedPoints, seg.points]), addedPoints
            def removePath(self,seg, newseg, nextPathL, pointsToFit):
                npoints = len(pointsToFit)
                acc=0
                newseg.next = seg.next                
                for p in reversed(nextPathL):
                    if (acc+len(p.points))<=npoints:
                        p.mergedObj = None
                        acc += len(p.points)
                    else:
                        newseg.prev = p        
                        p.points = p.points[(npoints-acc-len(p.points)):]                        
                        break
                

        newSegs = SegmentExtender.extendSegments( newSegs )
        debug("extended segs", newSegs)
        # ----------------------------------------

        ## # convert 1 point Path into 2 Segments to its prev and next
        ## updatedSegs =[newSegs[0]]
        ## for p in newSegs[1:-1]:
        ##     if p.isSegment() or len(p.points)>1:
        ##         updatedSegs.append(p)
        ##         continue            
        ##     newseg1 = Segment.from2Points(p.prev.points[-1], p.points[0], p.points)
        ##     newseg1.startIndexSource = p.startIndexSource
        ##     newseg2 = Segment.from2Points(p.points[0], p.next.points[0], p.points)
        ##     newseg2.startIndexSource = p.startIndexSource
        ##     updatedSegs +=[ newseg1, newseg2 ]
        ## updatedSegs.append(newSegs[-1])
        ## newSegs=resetPrevNextSegment(updatedSegs)
            

        # ---------------------------------------
        # merge consecutive segments with close angle
        updatedSegs=[]
        def toMerge(seg, mangle):
            l=[seg]
            setattr(seg, 'merged', True)
            if seg.next and seg.next.isSegment() :
                debug('merging segs ', seg.angle, ' with : ' ,seg.next.point1, seg.next.pointN, ' ang=',seg.next.angle)
                if closeAngleAbs( seg.angle, seg.next.angle) < mangle:
                    l += toMerge(seg.next,mangle)
            ## if seg.next and len(seg.next.points)==1:
            ##     l+=[seg.next] # try to merge a 1-point path
            return l
        def mergeList( segList , mangle =0.25 , q=0.5):
            updatedSegs = []
            for i,seg in enumerate(segList[:-1]):
                if not seg.isSegment() :
                    updatedSegs.append(seg)
                    continue
                if  hasattr(seg,'merged'):
                    continue
                debug(i,' inspect merge : ', seg.point1,'-',seg.pointN, seg.angle , ' q=',seg.quality())
                mList = toMerge(seg, mangle)
                debug('  --> tomerge ', len(mList))
                if len(mList)<2:
                    delattr(seg, 'merged')
                    updatedSegs.append(seg)
                    continue
                #points= seg.sourcepoints[mList[0].startIndexSource:mList[-1].startIndexSource+len(mList[-1].points)+1]
                points= numpy.concatenate( [p.points for p in mList] )
                newseg = fitSingleSegment(points)
                if newseg.quality()>q:
                    delattr(seg, 'merged')
                    updatedSegs.append(seg)
                    continue
                for p in mList:
                    setattr(seg, 'merged',True)
                newseg.sourcepoints = seg.sourcepoints
                debug('  --> post merge qual = ', newseg.quality() , seg.pointN, ' --> ', newseg.pointN, newseg.angle)
                newseg.startIndexSource = mList[0].startIndexSource
                newseg.prev = mList[0].prev
                newseg.next = mList[-1].next
                updatedSegs.append(newseg)
            if not hasattr(segList[-1], 'merged') : updatedSegs.append( segList[-1])
            return updatedSegs

        newSegs = mergeList( newSegs , mangle=0.2 )
        newSegs=resetPrevNextSegment(newSegs)
        debug(' __ 2nd angle merge')
        newSegs = mergeList( newSegs, mangle=0.35 ) # 2nd pass
        newSegs=resetPrevNextSegment(newSegs)
        ## for seg in newSegs:
        ##     debug( seg , seg.point1, seg.pointN, seg.angle)
        debug('after merge ', len(newSegs), newSegs)

        # -----------------------------------------------------
        # remove negligible Path/Segments between 2 large Segments
        self.removeSmallEdge(newSegs , wTot, hTot)
        newSegs=resetPrevNextSegment(newSegs)

        debug('after remove small ', len(newSegs),newSegs)
        # -----------------------------------------------------

        # -----------------------------------------------------
        # Extend segments to their intersections
        for p in newSegs:
            if p.isSegment() and p.next:
                p.setIntersectWithNext()
        # -----------------------------------------------------

        if forceClose:            
            newSegs[-1].setIntersectWithNext(newSegs[0])
        
        return PathGroup(newSegs, svgCommandsList, refNode)


    def simplifyPath(self, parsedpath, refNode):
        """ OBSOLETE ? Run the segment extraction technique on parsed path.
        returns a list of Segments or Paths"""
        a, parsedpath = toArray(parsedpath)
        self.currentPoints = a


        def distancesRatios():
            debug("DDD")
            self.segs = self.distancesRatios(self.currentPoints)
            self.segs.sort( key = lambda x : x.startIndexSource  )
            resetPrevNextSegment(self.segs)
            self.removeSmallEdge(self.segs)

            return self.segs

        paths = distancesRatios()
        for p in paths:
            if p.isSegment() and p.next:
                p.setIntersectWithNext()
                
        return PathGroup(paths, parsedpath, refNode )



    def extractShapes( self, nodes ):
        """The main function.
        nodes : a list of nodes"""
        analyzedNodes = []

        # convert nodes to list of segments (PathGroup) or Circle
        for n in nodes :
            parsedList = simplepath.parsePath(n.get('d'))
                      
            #g = self.simplifyPath( parsedList, n) 
            #g= self.tangeantEnvelop( parsedList, n )
            g= self.segsFromTangeants( parsedList, n )
            debug(" build group ", g, g.refNode, n)
            analyzedNodes.append( g )

        # uniformize shapes
        if 1:
            allSegs = [ p  for g in analyzedNodes for p in g.listOfPaths if p.isSegment() ]

            self.prepareParrallelize(allSegs)
            self.adjustToKnownAngle(allSegs)

            for group in analyzedNodes:
                # first pass : independently per path
                adjustAllAngles(group.listOfPaths)
                group.listOfPaths[:] = mergeConsecutiveParralels(group.listOfPaths)
                self.prepareDistanceEqualization([p for p in group.listOfPaths if p.isSegment()], 0.12)
                adjustAllDistances(group.listOfPaths)            
            ## # then 2nd global pass, with tighter criteria
            self.prepareDistanceEqualization(allSegs, 0.05)        
            for group in analyzedNodes:
                adjustAllDistances(group.listOfPaths)            

        # add new shapes in SVG document
        for group in analyzedNodes:            
            debug("final ", group.listOfPaths, group.refNode )
            # change to Rectangle if possible :
            finalshape = toRemarkableShape( group )
            ele = finalshape.addToNode( group.refNode)
            style = group.refNode.get('style')
            ele.set('style', style)
            if not self.options.keepOrigin:
                group.refNode.xpath('..')[0].remove(group.refNode)


        
if __name__ == '__main__':
    e = ShapeReco()
    e.affect()
