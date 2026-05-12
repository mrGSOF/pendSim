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
A subset of my bigger library for MATRIX operations without any external depandancies.
I deliberately wrote (most) of the code in functional programming style.
I wrote the code as procedures to simplify migration to C.
"""

from math import cos, sin
import copy

def matrix(rows, cols, val=0) -> list:
    """ Returns the rows by cols matrix M filled with value val """
    M = [0]*rows
    for row in range(0,rows):
        M[row] = [val]*cols
    return M

def zeros(rows, cols) -> list:
    """ Returns the rows by cols zero matrix Z """
    return matrix(rows, cols, val=0)

def DCM_V2(rad) -> list:
    """ Return the 2D rotation matrix """
    cosA = cos(rad)
    sinA = sin(rad)
    return [[cosA,-sinA],[sinA,cosA]]

def rotateV2(rad, V) -> list:
    """Rotate the vector V by rad degrees"""
    return MxV(DCM_V2(rad), V)

def getCol(M, col) -> list:
    """ Returns a copy of column 'col' from the matrix 'M' """
    rows = len(M)
    V = [0]*rows
    for i, row in enumerate(M):
        V[i] = row[col]
    return V

def T(M) -> list:
    """ Returns the transposed Matrix of M """
    rows = len(M)
    try:
        cols = len(M[0])
    except:
        cols = rows
        rows = 1
        M = [M]
    O = [0]*cols
    for i in range(0,cols):
        O[i] = getCol(M,i)
    return O

def MxV(M,V) -> list:
    """
    Return the result of NxM matrix and M vector multiplication
    Matrix structure: M[row][col]
    """
    O = [0]*len(M)
    for r,row in enumerate(M):
        for m,v in zip(row,V):
            O[r] += m*v
    return O

def getM(M) -> list:
    """Return a copy of matrix M"""
    return copy.deepcopy(M)

def putCol(M, col, vCol) -> None:
    """ Write column 'col' of the matrix 'M' with vCol"""
    rows = len(M)
    for i in range(0, rows):
        M[i][col] = vCol[i]

def putRow(M, row, vRow) -> None:
    """ Write row 'row' of the matrix 'M' with vRow"""
    cols = len(M[row])
    for i in range(0, cols):
        M[row][i] = vRow[i]

def rowColRemove(A, row=-1, col=-1) -> list:
    """Return the matrix A after removing col and row"""
    B = copy.deepcopy(A)
    if row >= 0:
        del B[row]

    if col >= 0:
        for row in B:
            del row[col]
    return B 

def Mxk(M,k) -> list:
    """
    Return the result a matrix M multiplied by a scaler
    Matrix structure: M[row][col]
    """
    rows, cols = shape(M)
    Mk = zeros(rows, cols)
    for i,row in enumerate(M):
        for j,elm in enumerate(row):
            Mk[i][j] = elm*k
    return Mk

def detM2(A) -> float:
    return A[0][0]*A[1][1] - A[0][1]*A[1][0]

def detM3(A) -> float:
    A1 = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
    B1 = A[2][0]*A[1][1]*A[0][2] + A[2][1]*A[1][2]*A[0][0] + A[2][2]*A[1][0]*A[0][1]
    return A1 - B1

def inv2( A ) -> list:
    """
    Matrix structure: M[row][col]
    Return the inverted  matrix 2x2.
    """
    a = A[0][0]
    b = A[0][1]
    c = A[1][0]
    d = A[1][1]
    s = 1.0/(a*d -b*c)
    a *= s
    b *= s
    c *= s
    d *= s
    return [[d,-b],[-c,a]]

def inv3( A ) -> list:
    """
    Matrix structure: M[row][col]
    Return the inverted matrix 3x3.
    """
    invA = [[0,0,0],[0,0,0],[0,0,0]]
    adjA = 1.0/detM3(A);
    for y in range(0,3):
        for x in range(0,3):
            invA[x][y] =adjA * detM2(rowColRemove(A,y,x)) #<Also transpose the matrix
            adjA = -1 * adjA
    return invA

def inv( A ) -> list:
    rows = len(A)
    if rows == 2:
        return inv2(A)
    if rows == 3:
        return inv3(A)
    else:
        return zeros(rows, rows)
