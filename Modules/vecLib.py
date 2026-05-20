## By: Guy Soffer (GSOF) 01/Sep/2025
__version__ = "1.0.0"
__author__ = "Guy Soffer"
__copyright__ = ""
__credits__ = [""]
__license__ = ""
__maintainer__ = ""
__email__ = "gsoffer@yahoo.com"
__status__ = "Development"

"""
Library for VECTOR (Linear-Algebra) operations without any external depandancies.
I deliberately wrote (most) of the code in functional programming style.
I wrote the code as procedures to simplify migration to C.
"""

import math

def printV(V, title="", prn=True):
    strV = str(title)
    if title != "":
        strV += ":\n"
    if (type(V) != list) and (type(V) != tuple):
        strV += "%1.3f\n"%(V)
    else:
        try:
            elm = V[0]
            strV += "["
            first = True
            for elm in V:
                if not first:
                    strV += ", "
                first = False
                strV += "%1.3f"%(elm)
            strV += "]"
        except:
            strV += str(V)
    if prn:
        print(strV+"\n")
    return strV

def dupV(V) -> list:
    """ Return a copy of vector V"""
    dim = len(V)
    O = [0]*dim
    for i, v in enumerate(V):
        O[i] = v
    return O

def intV(v) -> list:
    dim = len(v)
    if dim == 2:
        return int(v[0]), int(v[1])
    if dim == 3:
        return int(v[0]), int(v[1]), int(v[2])
    out = [0]*dim
    for i, elm in enumerate(v):
        out[i] = int(elm)
    return out

def mag2V(V) -> float:
    """ Return the magnitude of a vector """
    dim = len(V)
    if dim == 2:
        return (V[0]**2) +(V[1]**2)
    if dim == 3:
        return (V[0]**2) +(V[1]**2) +(V[2]**2)
    mag2 = 0.0
    for v in V:
        mag2 += v**2
    return mag2

def absV(V) -> float:
    """ Return the magnitude of a vector """
    return math.sqrt(mag2V(V))

def addV(V1, V2) -> list:
    """ Return the result of vector addition """
    dim = len(V1)
    if dim == 2:
        return [V1[0]+V2[0], V1[1]+V2[1]]
    if dim == 3:
        return [V1[0]+V2[0], V1[1]+V2[1], V1[2]+V2[2]]
    O = [0]*dim
    for i, (v1,v2) in enumerate(zip(V1, V2)):
        O[i] = v1+v2
    return O

def subV( V1, V2 ) -> list:
    """Add two vectors""" 
    vectorOut = [0]*len(V1)
    for i, (v1, v2) in enumerate(zip(V1, V2)):
        vectorOut[i] = v1 -v2
    return vectorOut

def dotV(V1, V2) -> float:
    product = 0.0
    for v1, v2 in zip(V1, V2):
        product += v1*v2
    return product

##def crossV3(V1, V2):
##    cross = [0]*len(V1)
##    for i,v1 in enumerate(V1):
##        for v2, j in enumerate(V2):
##            if i < j:
##                cross[i] += v1*v2
##            elif i > j:
##                cross[i] += -v1*v2
##    return cross

def crossV3( V1, V2 ) -> list:
    """Computes the cross product of two vectors"""
    vectorOut = [0]*3
    vectorOut[0]= (V1[1]*V2[2]) - (V1[2]*V2[1])
    vectorOut[1]= (V1[2]*V2[0]) - (V1[0]*V2[2])
    vectorOut[2]= (V1[0]*V2[1]) - (V1[1]*V2[0])
    return vectorOut

def proj(V1, V2):
    """Return the component of V1 on V2"""
    return scaleV(V2, dotV(V1, V2)/mag2V(V2))

def scaleV( V, scale ):
    """Multiply the vector by a scalar"""
    vectorOut = [0]*len(V)
    for i,v in enumerate(V):
        vectorOut[i] = v*scale
    return vectorOut

def normV(V) -> list:
    return scaleV(V, 1/absV(V))
    
def angleV2(V) -> float:
    """ Return the angle of 2D vector """
    return math.atan2(V[1], V[0])

def polarV2(V) -> float:
    """ Return the polar coordinates of 2D vector """
    return [absV(V), angleV2(V)]

def angleV3(v) -> list:
    r = absV3(v)
    _absV2 = absV2( v[0:2] )
    elevation = pi / 2.0                          #< For case 0 and 1
    #elevation = 0.0                                  #< For case2
    if _absV2 > 0.001:
        elevation = math.arctan( v[2] / _absV2 )  #< Case1 - https://keisan.casio.com/exec/system/1359533867 (Passed targeting test)
        #elevation = -math.arctan2( _absV2, v[2]) +pi/2  #< Case2 - https://keisan.casio.com/exec/system/1359533867
        #elevation = math.arcsin( r / v[2])               #< Case3 - https://www.mechamath.com/trigonometry/cartesian-to-spherical-coordinates-formulas-and-examples/
    azimuth = math.arctan2(v[1], v[0])
    return azimuth, elevation, r

def cartesianToPolarV3(pos):
    return angleV3(pos)

def polarToCartesianV3(azimuth, elevation, distance):
    abs_xy = distance * math.cos(elevation)
    x = abs_xy * math.cos(azimuth)
    y = abs_xy * math.sin(azimuth)
    z = distance * math.sin(elevation)
    return x, y, z

def concatenateV(*V) -> list:
    O = []
    for v in V:
        O += v    #< Concatenate all vectors
    return O

def splitV(V, sections=2) -> list:
    rows = len(V)
    N = len(V)
    secSizeList = sections
    if type(sections) == int:
        secSize = int(N/sections)
        secSizeList = [secSize]*(sections -1)
        secSizeList += [N -secSize*(sections -1)]

    O = [0]*len(secSizeList)
    stIdx = 0
    for i,size in enumerate(secSizeList):
        O[i] = V[stIdx:stIdx+size]
        stIdx += size
    return O

if __name__ == "__main__":
##    try:
##        import pysole
##        pysole.probe(runRemainingCode=True,     #< Execute the code below probe
##                     printStartupCode=True,     #< Print the command as well as it output
##                     primaryPrompt="vecLib>>> " #< 
##             )
##    except:
##        print("pip install liveConsole\nFor best interactive tool\n")

##    import pysole
##    pysole.probe(runRemainingCode=True,     #< Execute the code below probe
##                 printStartupCode=True,     #< Print the command as well as it output
##                 primaryPrompt="vecUT>>> " #< 
##         )
        
    V = [1,2,3]
    printV( V, "Test vector V" )
    a,b = splitV(V, 2)
    printV( a, "part1 Split(V, 2)" )
    printV( b, "part2 Split(V, 2)" )
    printV( concatenateV(V, V, V), "concatenateV(V,V,V)" )
    printV( scaleV(V, 1.5), "1.5V" )
    printV( dotV(V, V), "V dot V" )
    printV( dotV([1,0], [0, 1]), "0 = [1,0] dot [0,1]" )
    printV( crossV3([1,0,0], [0,1,0]), "[0,0,1] = [1,0,0] cross [0,1,0]" )
    printV( mag2V(V), "abs(V)^2" )
    printV( absV(V), "abs(V)" )

    V = [1,2,3,4,5,6,7]
    printV( V, "Test vector V" )
    printV( splitV(V, 3), "Split(V, 3)" )
    printV( splitV(V, [2,3,2]), "Split(V, [2,3,2])" )

    V = [1,2,3,4,5,6,7,8]
    printV( V, "Test vector V" )
    printV( splitV(V, 3), "Split(V, 3)" )
    printV( splitV(V, 2), "Split(V, 2)" )

    V = [1,2,3,4,5,6,7,8,9]
    printV( V, "Test vector V" )
    printV( splitV(V, 3), "Split(V, 3)" )
