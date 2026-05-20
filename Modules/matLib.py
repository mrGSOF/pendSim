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
Library for MATRIX (Linear-Algebra) operations without any external depandancies.
I deliberately wrote (most) of the code in functional programming style.
I wrote the code as procedures to simplify migration to C.
"""

import copy, math
from Modules.vecLib import *

def printM(M, title="") -> None:
    out = title
    if title != "":
        out += ":\n"
    try:    
        for row in M:
            out += "|"
            first = True
            for elm in row:
                if not first:
                    out += ", "
                first = False
                out += "% 1.3f"%(elm)
            out += "|\n"
    except:
        print(str(M))
    print(out)
    
def matrix(rows, cols, val=0) -> list:
    """ Returns the rows by cols matrix M filled with value val """
    M = [0]*rows
    for row in range(0,rows):
        M[row] = [val]*cols
    return M

def zeros(rows, cols) -> list:
    """ Returns the rows by cols zero matrix Z """
    return matrix(rows, cols, val=0)

def ones(rows, cols) -> list:
    """ Returns the rows by cols zero matrix Z """
    return matrix(rows, cols, val=1)

def shape(A) -> list:
    """Return the number of rows and columns of matrix A"""
    return (len(A), len(A[0]))

def I(size) -> list:
    """ Returns the size by size identity matrix I """
    I = [0]*size
    for row in range(0,size):
        I[row] = [0]*size
        I[row][row] = 1
    return I

def getM(M) -> list:
    """Return a copy of matrix M"""
    return copy.deepcopy(M)

def getCol(M, col) -> list:
    """ Returns a copy of column 'col' from the matrix 'M' """
    rows = len(M)
    V = [0]*rows
    for i, row in enumerate(M):
        V[i] = row[col]
    return V

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

def getRow(M, row) -> list:
    """Return a row from the matrix M"""
    return copy.deepcopy(M[row])

def rowColRemove(A, row=-1, col=-1) -> list:
    """Return the matrix A after removing col and row"""
    B = copy.deepcopy(A)
    if row >= 0:
        del B[row]

    if col >= 0:
        for row in B:
            del row[col]
    return B 

def swapRow(A, row1, row2) -> list:
    """Return the matrix A after swapping position of row1 and row2"""
    B = copy.deepcopy(A)
    if (row1 != row2):
        tmpRow = B[row1]
        B[row1] = B[row2]
        B[row2] = tmpRow
    return B

def swapCol(A, col1, col2) -> list:
    """Return the matrix A after swapping position of col1 and col2"""
    B = copy.deepcopy(A)
    for row in B:
        tmp = row[col1]
        row[col1] = row[col2]
        row[col2] = tmp
    return B

def concatenateM(*A) -> list:
    """Return the joined matrices A1, A2, A3..."""
    rows = len(A[0])
    M = [0]*rows
    for i in range(0,rows): #< Iterate over each row
        M[i] = []
        for a in A:
            M[i] += a[i]    #< Concatenate rows from all matrices
    return M

def hsplitM(A, sections=2) -> list:
    rows, N = shape(A)
    secSize = int(N/sections)
    secSizeList = [secSize]*(sections -1)
    secSizeList += [N -secSize*(sections -1)]
    O = [ [] for _ in range(sections) ]
    for r, row in enumerate(A):
        splitRow = splitV(row, sections=secSizeList)
        for column, section in enumerate(splitRow):
            O[column].append(section) 
    return O

def DCM_V2(rad) -> list:
    """ Return the 2D rotation matrix """
    cosA = math.cos(rad)
    sinA = math.sin(rad)
    return [[cosA,-sinA],[sinA,cosA]]

def rotateV2(rad, V) -> list:
    """Rotate the vector V by rad degrees"""
    return MxV(DCM_V2(rad), V)

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

def MxV(M,V) -> list:
    """
    Return the result of NxM matrix and M vector multiplication
    Matrix structure: M[row][col]
    """
    rows = len(M)
    O = [0]*rows
#    if rows == 2:
#        ### Only is M is 2x2
#        O[0] = M[0][0]*V[0] +(M[0][1]*V[1])
#        O[1] = M[1][0]*V[0] +(M[1][1]*V[1])    
#        return O
    for r,row in enumerate(M):
        for m,v in zip(row,V):
            O[r] += m*v
    return O

def detM2(A) -> float:
    """Returns the determinant of a 2x2 matrix"""
    return A[0][0]*A[1][1] - A[0][1]*A[1][0]

def detM3(A) -> float:
    A1 = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
    B1 = A[2][0]*A[1][1]*A[0][2] + A[2][1]*A[1][2]*A[0][0] + A[2][2]*A[1][0]*A[0][1]
    return A1 - B1

def detM(A) -> float:
    r = len(A)
    if r == 2:
        return detM2(A)
    if r == 3:
        return detM3(A)
    else:
        print("Determinant for matrices bigger than 4x4 isn't supported")
        return -1
    
def genR(a) -> list:
    """ Return the 2D rotation matrix """
    cosA = math.cos(a)
    sinA = math.sin(a)
    return [[cosA,-sinA],[sinA,cosA]]

def negR(R) -> list:
    R[0][1] *= -1
    R[1][0] *= -1
    return R

def T(M) -> list:
    """ Returns the transposed Matrix of M """
    try:
        cols = len(M[0])
    except:
        ### The input is a vector
        N = len(M) #< Number of elements
        O = [0]*N
        for i,elm in enumerate(M):
            O[i] = [elm] #< Make a column vector
        return O

    ### The input is a matrix
    O = [0]*cols
    for i in range(0,cols):
        O[i] = getCol(M,i)
    return O

def M2x2xV2( M, V ) -> list:
    """
    Multiply matrix 2x2 by vector.
    Matrix structure: M[row][col]
    """
    return [M[0][0]*V[0] +M[0][1]*V[1],
            M[1][0]*V[0] +M[1][1]*V[1]]

def M3x3xV3( M, V ) -> list:
    """
    Multiply matrix 3x3 by vector.
    Matrix structure: M[row][col]
    """
    return [M[0][0]*V[0] +M[0][1]*V[1] +M[0][2]*V[2],
            M[1][0]*V[0] +M[1][1]*V[1] +M[1][2]*V[2],
            M[2][0]*V[0] +M[2][1]*V[1] +M[2][2]*V[2]]

def MxM2x2( A, B ) -> list:
    """
    Multiply two matrices 2x2.
    Matrix structure: M[row][col]
    """
    M00 = A[0][0]*B[0][0]+ A[0][1]*B[1][0]
    M01 = A[0][0]*B[0][1]+ A[0][1]*B[1][1]
    M10 = A[1][0]*B[0][0]+ A[1][1]*B[1][0]
    M11 = A[1][0]*B[0][1]+ A[1][1]*B[1][1]
    return [[M00,M01],
            [M10,M11]]

def MxM( A, B ) -> list:
    """
    Multiply two matrices 3x3.
    Matrix structure: M[row][col]
    """
    rows = len(A)
    cols = len(B[0])
    matrixOut = zeros(rows, cols)
    for y in range(0,rows):
        for x in range(0,cols):
            for i in range(len(B)):
                matrixOut[y][x] += A[y][i]*B[i][x]
    return matrixOut

def detM2(A) -> float:
    return A[0][0]*A[1][1] - A[0][1]*A[1][0]

def detM3(A) -> float:
    A1 = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
    B1 = A[2][0]*A[1][1]*A[0][2] + A[2][1]*A[1][2]*A[0][0] + A[2][2]*A[1][0]*A[0][1]
    return A1 - B1

def detM(A) -> float:
    """Return the inverted matrix A"""
    rows, cols = shape(A)
    if (rows == 2) and (cols == rows):
        return detM2(A)
    if (rows == 3) and (cols == rows):
        return detM3(A)
    L,U,P,s = LU(A, swaps=True)
    sg = 1
    if s&1:
        sg = -1 #< Odd number of swaps
    return prodDiag(U)*sg

def addM( A, B ) -> list:
    """Add two NxM matrices""" 
    rows=len(A)
    cols=len(A[0])
    matrixOut = zeros(rows, cols)
    for row, (rowA, rowB) in enumerate(zip(A,B)):
        for col, (elmA, elmB) in enumerate(zip(rowA, rowB)):
            matrixOut[row][col] = elmA +elmB
    return matrixOut

def subM( A, B ) -> list:
    """Add two NxM matrices""" 
    rows=len(A)
    cols=len(A[0])
    matrixOut = zeros(rows, cols)
    for row, (rowA, rowB) in enumerate(zip(A,B)):
        for col, (elmA, elmB) in enumerate(zip(rowA, rowB)):
            matrixOut[row][col] = elmA -elmB
    return matrixOut

def negM( A ) -> list:
    """Negate NxM matrix""" 
    rows=len(A)
    cols=len(A[0])
    negA = zeros(rows, cols)
    for row, rowA in enumerate(A):
        for col, elmA in enumerate(rowA):
            negA[row][col] = -elmA
    return negA

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
    return Mxk( [[d,-b],[-c,a]], s )

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

def inv(A, v=False) -> list:
    """Return the inverted matrix A"""
    rows = len(A)
    cols = len(A[0])
    if (rows == 2) and (cols == rows):
        return inv2(A)
    if (rows == 3) and (cols == rows):
        return inv3(A)
    
    A = concatenateM(A, I(rows))
    H = highM(A, v)
    L = lowM(H, v)
    return hsplitM(L,2)[1]

def partialPivoting(A, stRow=0, v=False) -> list:
    """Fix ill-conditioned matrix by partial pivoting method"""
    rows = len(A)
    scale = abs(A[stRow][stRow])
    bigScaleRow = stRow
    ### Search for the row with the biggest value
    for i in range(stRow+1, rows):
        newScale = abs(A[i][stRow])
        if (newScale > scale):
            bigScaleRow = i
            scale = newScale
    return bigScaleRow

def getPermutationM(rows, row1, row2) -> list:
    """Return the Permutation matrix to swap row1 with row2"""
    P = I(rows) #< No swapping by default
    if (row1 != row2):
        ## Update the Permutation matrix
        P[row1][row2] = 1
        P[row2][row1] = 1
        P[row1][row1] = 0
        P[row2][row2] = 0
    return P

def _swapRowsM(A, row1, row2, cols, v=False) -> None:
    """swap row1 and row2 in matrix A"""
    tmpRow = A[row1][0:cols]
    A[row1][0:cols] = A[row2][0:cols]
    A[row2][0:cols] = tmpRow
    if v:
        printM(A, "After swapRow(A, %d, %d)"%(row1,row2))

##def LU(A, v=False) -> list:
##    """Return the Lower and Upper (LU) decomposition of matrix A"""
##    ### A=LU
##    rows = len(A)
##    cols = rows
##    U = copy.deepcopy(A)
##    L = I(rows)
##    for col in range(0, rows-1):
##        stRow = col
##        refRow = U[stRow]
##        for i in range(stRow +1, cols):
##            row = U[i]
##            if refRow[col] == 0.0:
##                print("Error")
##                return L,U
##            scale = row[col]/refRow[col]
##            L[i][col] = scale
##            U[i] = addV(row, scaleV(refRow, -scale))
##    return L,U

def LU(A, swaps=False, v=False) -> list:
    """Return the Lower, Upper, and Permutation (LUP) decomposition of matrix A"""
    ### PA=LU
    rows = len(A)
    cols = rows
    U = copy.deepcopy(A)
    L = I(rows)
    P = I(rows) #< Finale Permutation matrix
    s = 0
    for col in range(0, rows-1):
        M = I(rows) #< Elementry elimination matrix
        stRow = col
        swapRow = partialPivoting(U, stRow)
        if stRow != swapRow:
            s += 1
            _swapRowsM(U, stRow, swapRow, cols=cols, v=False)
            _swapRowsM(L, stRow, swapRow, cols=col,  v=False)
            _swapRowsM(P, stRow, swapRow, cols=cols, v=False)
            
        refRow = U[stRow]
        for i in range(stRow +1, cols):
            row = U[i]
            scale = row[col]/refRow[col]
            #M[i][col] = -scale
            L[i][col] = scale
            U[i] = addV(row, scaleV(refRow, -scale))
        #printM(M, "Elementry elimination matrix")
        #U = MxM(M, U) #< Eliminate rows
    if swaps:
        return L,U,P,s
    return L,U,P

def LDV(A) -> list:
    """Returns the LDV decomposition of matrix A"""
    #PA=LU=LDV
    #A=T(P)xLxDxV |where D is diagonal matrix and V is an upper triangular matrix 
    L,U,P = LU(A)
    invD = getDiag(U, lambda x:1.0/x)
    D = getDiag(U)
    V = MxM(invD, U)
    return L, D, V, P

def LDLT(A) -> list:
    """Return the LDL^T decomposition of A matrix"""
    C = Cholesky(A)
    
    if len(C) > 0:
        n = len(A)
        D = I(n)
        L = C
        for i in range(0,n):
            D[i][i] = C[i][i]**2
            #L[:][i] = C[:][i]*D[i][i]**(-1/2) #< The two lines below
            for j in range(0,n):
                L[j][i] = C[j][i]*D[i][i]**(-1/2)
        return L,D
    else:
        return [],[]        

def gramSchmidt(A) -> list:
    """Return the Gram-Schmidt transfer"""
    rows, cols = shape(A)
    Y = matrix(rows, cols, 0)
    putCol(Y, 0, normV(getCol(A,0)))
    for col in range(1,cols):
        v = getCol(A, col)
        x = dupV(v)
        for yi in range(0, col):
            v = subV(v, proj(x, getCol(Y, yi)))
        putCol(Y, col, normV(v))
    return Y

##def gramSchmidt(A) -> list:
##    """Return the Gram-Schmidt transfer"""
##    A = T(A)
##    rows, cols = shape(A)
##    Y = matrix(rows, cols, 0)
##    Y[0] = normV(A[0])
##    for row, v in enumerate(A[1:]):
##        row += 1
##        x = dupV(v)
##        for y in Y[0:row]:
##            v = subV(v, proj(x, y))
##        Y[row] = normV(v)
##    return T(Y)

def QR_gs(A) -> list:
    rows, cols = shape(A)
    Q = gramSchmidt(A)
    R = matrix(rows, cols, 0)
    for ri, row in enumerate(A):
        q = getCol(Q,ri)
        for ci in range(ri, cols):
            a = getCol(A,ci)
            R[ri][ci] = dotV(q,a)
    return Q,R

def givensrotation(a, b):
    hypot = math.sqrt(a**2 + b**2)
    cos = a / hypot
    sin = -b / hypot
    return cos, sin

def qr_givens(A):
    m, n = shape(A)
    R = dup(A)
    Q = I(m)
    for i in range(0, n - 1):
        for j in range(i + 1, m):
            cos, sin = givensrotation(R[i][i], R[j][i])
            R[i], R[j] = (R[i] * cos) + (R[j] * (-sin)), (R[i] * sin) + (R[j] * cos)
            Q[:, i], Q[:, j] = (Q[:, i] * cos) + (Q[:, j] * (-sin)), (Q[:, i] * sin) + (Q[:, j] * cos)
    return Q, R

def Cholesky(A) -> list:
    """Return the LxL^T decomposition of A matrix"""
    ## By Nick Space Cowboy
    ## A = L x T(L)
    if isPosDefinate(A):
        n = len(A)
        L = matrix(n, n, 0.0)
        for i in range(0,n):
            for j in range(0, i+1):
                sumPart = 0.0
                if i ==j:
                    for k in range(0,j):
                        sumPart += L[j][k]**2
                    L[i][j] = math.sqrt(A[j][j] -sumPart)
                else:
                    for k in range(0,j):
                        sumPart += L[i][k]*L[j][k]
                    L[i][j] = (A[i][j] -sumPart)/L[j][j]
        return L
    else:
        print("A isn't positive definate matrix")
        return []

def highM(A, v=False) -> list:
    """ Return the (U)pper triangle of M (Using Naive Gauss elimination method) """
    H = copy.deepcopy(A)
    rows = len(H)
    cols = len(H[0])
    for col in range(0, rows-1):
        stRow = col
        bigScaleRow = partialPivoting(H, stRow)
        _swapRowsM(H, stRow, bigScaleRow, cols)
        scale = H[stRow][stRow]
        if scale == 0.0:
            printM(H, "Error highM()")
            return H

        H[stRow] = scaleV(H[stRow], 1.0/scale) #< Norm row
        refRow = H[stRow]
        for i in range(stRow +1, rows):
            row = H[i]
            refRowScaled = scaleV(refRow, -row[col])
            H[i] = addV(refRowScaled, row)
    return H

def lowM(A, v=False) -> list:
    """ Return the (L)ow triangle of M (Using Naive Gauss elimination method) """
    L = copy.deepcopy(A)
    rows = len(L)
    cols = len(L[0])
    for col in range(rows-1, 0, -1):
        stRow = col
        scale = L[stRow][stRow]
        if scale == 0.0:
            printM(L, "Error lowM()")
            return L
        L[stRow] = scaleV(L[stRow], 1.0/scale) #< Norm row
        refRow = L[stRow]
        for i in range(col -1, -1, -1):
            #printM(L, "row:%d, col:%d"%(i,col))
            row = L[i]
            refRowScaled = scaleV(refRow, -row[col])
            L[i] = addV(refRowScaled, row)
    scale = L[0][0]
    if scale == 0.0:
        printM(L, "Error lowM()")
    else:
        L[0] = scaleV(L[0], 1.0/scale)
    return L

def isZero(A) -> bool:
    s = 0.0
    for row in A:
        s += sum(row)
    return round(s,10) == 0.0

def findClosestToZero(A, row) -> int:
    row = A[row]
    idx = 0
    smallest = abs(row[0])
    for i,elm in enumerate(row[1:]):
        elm = abs(elm)
        if smallest > elm:
            smallest = elm
            idx = i+1
    return idx

def isOrthogonal(A) -> bool:
    return isZero( subM(MxM(A, T(A)), I(len(A))) )

def isOrthonormal(A) -> bool:
    AxTA  = MxM(A, T(A))
    zeroA = isZero( subM(AxTA, I(len(A))) )
    detA  = round(abs(detM(A)),10)
    return (detA == 1.0) and zeroA

def isSymetric(A) -> bool:
    rows = len(A)
    cols = len(A[0])
    for r in range(0, rows):
        for c in range(0, cols):
            if r != c:
                if A[r][c] != A[c][r]:
                    return False
    return True

def isPosDefinate(A) -> bool:
    """Check if a matrix is positive definate"""
    ## By Nick Space Cowboy
    for i in range(1,len(A)):
        if (A[i][i-1]**2) >= A[i][i]:
            return False
    return True

def isSimilar(A) -> bool:
    return False

def prodDiag(A) -> float:
    p = 1.0
    for i in range(0,len(A)):
        p *= A[i][i]
    return p

def sumDiag(A) -> float:
    s = 0.0
    for i in range(0,len(A)):
        s += A[i][i]
    return s

def getDiag(A, f=lambda x:x) -> list:
    n = len(A)
    D = matrix(n, n, 0.0)
    for i in range(0,n):
        D[i][i] = f(A[i][i])
    return D
    
##def isEvenP(P) -> int:
##    return 1

def eig(A) -> list:
    return []

def evd(A) -> list:
    return []

def svd(A) -> list:
    return []

def solve(A, b) -> list:
    """Return the solution for Ax=b and the inverted matrix A"""
    invA = inv(A)
    return (MxV(invA, b), invA)

def solveLU(A, b) -> list:
    """Return the solution for Ax=b and the inverted matrix A"""
    L,U,P = LU(A)
    ## Ax=b -> PAx=Pb |PA=LU
    ## LUx=Pb -> Ly=Pb |y=Ux

    #y=L\(Pb)
    Pb = MxM(P,T(b))
    LPb = concatenateM(L,Pb)
    lastColumn = len(LPb[0])-1
    highLPb = highM(LPb)
    printM(LPb, "LPb")
    printM(highLPb, "highLPb")
    y = getCol(highM(LPb), lastColumn)

    #x=u\y
    Uy = concatenateM(U,T(y))
    lowUy = lowM(Uy)
    printM(Uy, "Uy")
    printM(lowUy, "lowUy")
    x = getCol(lowUy, lastColumn)
    return (x, L,U,P)

if __name__ == "__main__":
    import time
    def testMatrix(A, rowSwap="") -> None:
        printM(A, "\nMatrix A")
        print("det(A): %1.3f\n"%(detM(A)))
        invA = inv(A)
        L,U,P = LU(A)
        printM(L, "L")
        printM(U, "U")
        printM(P, "P {%s}"%(rowSwap))
        printM(inv(P), "Inverted P")
        printM(T(P), "Transposed P")
        printM(MxM(P,A), "P x A")
        printM(MxM(L,U), "P x A = L x U")
        printM(invA, "invA")
        printM(MxM(MxM(inv(U), inv(L)), P), "invA = invU x invL x P")
        printM(MxM(invA, A), "I = A x invA")
        print( "P Orthogonal: %s\n"%(str(isOrthogonal(P))) )
        print( "P Orthonormal: %s\n"%(str(isOrthonormal(P))) )
        print( "A Orthogonal: %s\n"%(str(isOrthogonal(A))) )

##    import pysole
##    pysole.probe(runRemainingCode=True,     #< Execute the code below probe
##                 printStartupCode=True,     #< Print the command as well as it output
##                 primaryPrompt="matUT>> "   #< Prompt
##         )

    printM( detM2(ones(2,2)), "Det2x2 of all ones is zero" )
    printM( I(3), "I(3)" )
    printM( zeros(4,3), "Zeros(4,3)" )
    printM( ones(3,4), "Ones(3,4)" )

    M = I(4)
    printM( M, "\nMatrix A" )
    M = rowColRemove(A=M, row=1, col=3)
    printM( M, "A without row 2 and col 1" )
    printM( T(M), "A transposed" )

##    A = [[1,2],
##         [5,6]]
##    B = [[5,4],
##         [1,0]]
##    printM(A, "\nMatrix A")
##    printM(B, "Matrix B")
##    printM(addM(A,B), "A+B")
##    printM(addM(A,A), "A+A")
##    printM(subM(A,A), "A+A")
##    printM(Mxk(A, -2), "-2A")
##    
##    invA = inv2(A)
##    A_x_invA = MxM(A, invA)
##    printM(invA, "invA")
##    printM(A_x_invA, "A x invA")

    A = [[1,2,3],
         [4,5,6],
         [7,8,11],
         ]
    A1 = swapRow(A,0,2)
    A2 = swapCol(A,1,2)
    
    invA = inv3(A)
    L,U,P = LU(A)
    A_x_invA = MxM(A, invA)
    printM(A, "\nMatrix A")
    printM(A1,"swapRow(A,0,2)")
    printM(A2,"swapCol(A,1,2)")
    printM(L, "L")
    printM(U, "U")
    printM(MxM(L,U), "A = L x U")
    printM(invA, "invA")
    printM(MxM(inv(U),inv(L)), "invA = invU x invL")
    printM(A_x_invA, "I = A x invA")

    A = [[1,2,3,4],
         [4,5,60,7],
         [7,8,9,11],
         [17,18,19,10],
         ]
    AA = concatenateM(A, A)
    A_A = hsplitM(AA, 2)
    A_A_A = hsplitM(AA, 3)
    lowA = lowM(A)
    highA = highM(A)
    invA = inv(A)
    L,U,P = LU(A)
    printM(A, "\nMatrix A")
    printM(L, "L")
    printM(U, "U")
    printM(MxM(L,U), "A = L x U")
    printM(invA, "invA")
    printM(MxM(inv(U), inv(L)), "invA = invU x invL")
    printM(MxM(A, invA), "I = A x invA")
    printM(AA, "concatenateM(A, A)")
    printM(A_A, "hspliteM(AA, 2)")
    printM(A_A_A, "hspliteM(AA, 3)")
    printM(lowA, "lowA")
    printM(highA, "highA")

    A = [[1,1,1],
         [4,3,-1],
         [3,5,3],
         ]
    testMatrix(A, rowSwap="3->1, 1->2, 2->3")

    A = [[0,2,1],
         [-2,-3,1],
         [3,5,3]]
    testMatrix(A, rowSwap="3->1, 1->2, 2->3")



##    A = swapRow(A,1,2)
##    printM(inv(A), "A with rows 1 and 2 swapped")

    A=[[0,1,0],
       [0,0,1],
       [1,0,0]]
    testMatrix(A, rowSwap="3->1, 1->2, 2->3")

    A=[[2,1,1],
       [4,3,3],
       [8,7,9]]
    testMatrix(A, rowSwap="3->1, 1->2, 2->2")

    A=[[2,0,-1],
       [5,1,0],
       [0,1,3]]
    testMatrix(A, rowSwap="2->1, 3->2, 1->3")

    A=[[0,0,0,2],
       [0,0,3,0],
       [0,4,0,0],
       [5,0,0,0]]
    testMatrix(A, rowSwap="4->1, 3->2, 2->3, 1->4")

    A=[[3,2,0,0],
       [4,3,0,0],
       [0,0,6,5],
       [0,0,7,6]]
    testMatrix(A, rowSwap="2->1, 1->2, 4->3, 3->4")

    A=[[0,1,0,0,0],
       [0,0,0,1,0],
       [1,0,0,0,0],
       [0,0,0,0,1],
       [0,0,1,0,0]]
    testMatrix(A, rowSwap="3->1, 1->2, 5->3, 2->4, 4->5")

    A=[[2,4,-2,-2],
       [1,2,4,-3],
       [-3,-3,8,-2],
       [-1,1,6,-3]]
    testMatrix(A, rowSwap="3->1, 1->2, 2->3, 4->4")

    A=[[2,0,2,0.6],
       [3,3,4,-2],
       [5,5,4,2],
       [-1,-2,3.4,-1]]
    testMatrix(A, rowSwap="3->1, 1->2, 4->3, 2->4")

    A=[[0,0,0,2],
       [0,-1,0,0],
       [0,0,3,0],
       [-2,0,0,0]]
    b=[2,-1,3,-2]
    dt = time.time()
    x, invA = solve(A, b)
    dt = 1000*(time.time() -dt)
    printM(A, "A (solving Ax=b)")
    printV(b, "b")
    printM(invA, "invA")
    printV(x, "x (inv(A) solution time %1.4f ms)"%(dt))

    dt = time.time()
    x, L,U,P = solveLU(A, b)
    dt = 1000*(time.time() -dt)
    printM(A, "A (solving Ax=b)")
    printV(b, "b")
    printM(L, "L")
    printM(U, "U")
    printM(P, "P")
    printV(x, "x (LU solution time %1.4f ms)"%(dt))
    printV([1,1,1,1], "Expected solution")

    A=[[2,-1,1],
       [1,1,0],
       [3,-1,-2]]
    b=[3,-1,7]
    dt = time.time()
    x, L,U,P = solveLU(A, b)
    dt = 1000*(time.time() -dt)
    printM(A, "A (solving Ax=b)")
    printM(P, "P")
    printV(x, "x (LU solution time %1.4f ms)"%(dt))
    dt = time.time()
    x, invA = solve(A, b)
    dt = 1000*(time.time() -dt)
    printV(x, "x (inv(A) solution time %1.4f ms)"%(dt))
    printV([1,-2,-1], "Expected solution")

    A = [[33,2,6],[2,73,5],[6,5,73]]
    posDef = isPosDefinate(A)
    L = Cholesky(A)
    printM(A, "A for Cholesky decomposition")
    printM(L, "A=LxT(L)")

    printM(A, "A for LDLT decomposition")
    L,D = LDLT(A)
    printM(L, "L ")
    printM(D, "D")
    printM(MxM(MxM(L,D), T(L)), "A=LxDxT(L)")

    A = [[90,81,57,62],[82,71,7,23],[72,61,93,76],[51,7,55,21]]
    L,D,V,P = LDV(A)
    ld=MxM(L,D)
    ldv=MxM(ld,V)
    ldvp=MxM(T(P),ldv)
    printM(A, "For LDV decomposition")
    printM(L, "L")
    printM(D, "D")
    printM(P, "P")
    printM(ldvp, "A = T(P)xLxDxV")

    A = [[90,2,6],[2,64,5],[6,5,28]]
    L,D = LDLT(A)
    ldlt = MxM(MxM(L,D),T(L))
    printM(A, "For LDLT decomposition")
    printM(L, "L")
    printM(D, "D")
    printM(ldlt, "A=LxDxT(L)")

    A = [[3,2,1],[3,3,3],[10,5,5]]
    printM(A, "For Gram-Schmidt transformation")
    A_gs = gramSchmidt(A)
    printM(A_gs, "Gram-Schmidt transformation")

    A = [[5,2,10],[5,2,7],[9,10,7]]
    printM(A, "For Gram-Schmidt transformation")
    A_gs = gramSchmidt(A)
    printM(A_gs, "Gram-Schmidt transformation")
    print("Orthonormal: %s"%str(isOrthonormal(A_gs)))

    A = [[6,2,3],[2,3,9],[2,3,2]]
    printM(A, "For Gram-Schmidt transformation")
    A_gs = gramSchmidt(A)
    printM(A_gs, "Gram-Schmidt transformation")
    print("Orthonormal: %s"%str(isOrthonormal(A_gs)))

    A = [[6,6,3],[1,5,1],[4,0,0]]
    printM(A, "For Gram-Schmidt transformation")
    A_gs = gramSchmidt(A)
    printM(A_gs, "Gram-Schmidt transformation")
    print("Orthonormal: %s"%str(isOrthonormal(A_gs)))

    A = [[6,6,3],[1,5,1],[4,0,0]]
    printM(A, "For QR_gs decomposition")
    Q,R = QR_gs(A)
    printM(Q, "Q")
    printM(R, "R")
    print("Orthonormal: %s"%str(isOrthonormal(Q)))

    A = [[5,2,10],[5,2,7],[9,10,7]]
    printM(A, "For QR_gs decomposition")
    Q,R = QR_gs(A)
    printM(Q, "Q")
    printM(R, "R")
    print("Orthonormal: %s"%str(isOrthonormal(Q)))
