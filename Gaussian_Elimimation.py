import math
from scipy import integrate

def rows(mat):
    "return number of rows"
    return(len(mat))

def cols(mat):
    "return number of cols"
    return(len(mat[0]))
 
def zero(m,n):
    "Create zero matrix"
    new_mat = [[0 for col in range(n)] for row in range(m)]
    return new_mat
 
def transpose(mat):
    "return transpose of mat"
    new_mat = zero(cols(mat),rows(mat))
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            new_mat[col][row] = mat[row][col]
    return(new_mat)

def dot(A,B):
    "vector dot product"
    if len(A) != len(B):
        print("dot: list lengths do not match")
        return()
    dot=0
    for i in range(len(A)):
        dot = dot + A[i]*B[i]
    return(dot)

def scalarVecMult(s,V):
    "return vector sV"
    C = []
    for i in range(len(V)):
        C.append(V[i]*s)
    return(C)  

def getCol(mat, col):
    "return column col from matrix mat"
    return([r[col] for r in mat])

def getRow(mat, row):
    "return row row from matrix mat"
    return(mat[row])

def matMult(mat1,mat2):
    "multiply two matrices"
    if cols(mat1) != rows(mat2):
        print("multiply: mismatched matrices")
        return()
    prod = zero(rows(mat1),cols(mat2))
    for row in range(rows(mat1)):
        for col in range(cols(mat2)):
            prod[row][col] = dot(mat1[row],getCol(mat2,col))
    return(prod)

def vectorQ(V):
    "mild test to see if V is a vector"
    if type(V) != type([1]):
        return(False)
    if type(V[0]) == type([1]):
        return(False)
    return(True)

def scalarMult(a,mat):
    "multiplies a scalar times a matrix"
    "If mat is a vector it returns a vector."
    if vectorQ(mat):
        return([a*m for m in mat])
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            mat[row][col] = a*mat[row][col]
    return(mat)

def addVectors(A,B):
    "add two vectors"
    if len(A) != len(B):
        print("addVectors: different lengths")
        return()
    return([A[i]+B[i] for i in range(len(A))])

def swaprows(M,i,j):
    "swap rows i and j in matrix M"
    N=copyMatrix(M)
    T = N[i]
    N[i] = N[j]
    N[j] = T
    return N

def copyMatrix(M):
    return([[M[row][col] for col in range(cols(M))]for row in
            range(rows(M))])

def addMatrices(A,B):
    "add two matrices"
    new_mat = []
    for row in range(rows(A)):
        tmp = []
        for col in range(cols(A)):
            tmp.append(A[row][col] + B[row][col])
        new_mat.append(tmp)
    return(new_mat)

def addrows(M, f, t, scale):
    "add scale times row f to row t"
    N=copyMatrix(M)
    T=addVectors(scalarMult(scale,N[f]),N[t])
    N[t]=T
    return(N)

def show(mat):
    "Print out matrix"
    for row in mat:
        print(row)

def augment(mat,vec):
    "given nxn mat and n length vector return augmented matrix"
    amat = []
    for row in range(rows(mat)):
        amat.append(mat[row]+[vec[row]])
    return(amat)

def copyVec(L):
    "return a copy of L"
    C=[k for k in L]
    return(C)


### The next two functions support checking a solution.

def getAandb(aug):
    "Returns the coef. matrix A and the vector b of Ax=b"
    m = rows(aug)
    n = cols(aug)
    A = zero(m,n-1)
    b = zero(m,1)
    for i in range(m):
        for j in range(n-1):
            A[i][j] = aug[i][j]
            
    for i in range(m):
        b[i] = aug[i][n-1]
    Aandb = [A,b]
    return(Aandb)

def checkSol_1(aug,x):
    "For aug=[A|b], returns Ax, b, and b-Ax as vectors"
    "Note: It assumes x is a vector not a matrix"
    A  = getAandb(aug)[0]
    b  = getAandb(aug)[1]
    x_col_vec = vec2colVec(x)
    Ax = matMult(A,x_col_vec)
    r  = addVectors(b,scalarMult(-1.0,colVec2vec(Ax)))
    L  = [Ax,b,r]
    return(L)

### Used for calculating the inverses
def identityMatrix(n):
    "Creates nxn identity matrix"
    I = []
    for row in range(n):
        tmp =[]
        for column in range(n):
            if row == column:
                tmp.append(1)
            else:
                tmp.append(0)
        I.append(tmp)
    return(I)

def checkRow(A, I, row):
    if A[row] != 1:
        scale = A[row]
        for i in range(len(A)):
            A[i] = (A[i]/scale)
            I[i] = (I[i]/scale)

                
##########################################
##          Naive Gaussian              ##
##########################################

def findPivotrow1(mat,col):
    "Finds index of the first row with nonzero entry on or"
    "below diagonal.  If there isn't one return(-1)."
    for row in range(col, rows(mat)):
        if mat[row][col] != 0:
            return(row)
    return(-1)


def rowReduce(M):
    "return row reduced version of M where M is the augmented matrixm [A | b]"
    N = copyMatrix(M)
    cs = cols(M)-2   # no need to consider last two cols
    rs = rows(M)
    for col in range(cs+1):
        j = findPivotrow1(N,col)
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swaprows(N,col,j)
            scale = -1.0 / N[col][col]
            for row in range(col+1,rs):                
                N=addrows(N, col, row, scale * N[row][col])
    return(N)

def inverseNaive(M):
    """ M is the nxn matrix A. """
    N = copyMatrix(M)
    I = identityMatrix(rows(N))
    cs = cols(M)-2   # no need to consider last two cols
    rs = rows(M)

    for col in range(cs+1):
        j = findPivotrow1(N,col)
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swaprows(N,col,j)
                I = swaprows(I, col,j)
            scale = -1.0 / N[col][col]
            for row in range(col+1,rs):
                tmp = scale * N[row][col]
                N=addrows(N, col, row, scale * N[row][col])
                I=addrows(I, col, row, tmp)
    for column in range(1,rs):
        checkRow(getRow(N, column), getRow(I,column), column)
        for row in range((column-1), -1, -1):
            if N[row][column] != 0:
                scale = -1.0 / N[column][column]
                tmp = scale * N[row][column]
                N=addrows(N, column, row, scale * N[row][column])
                I=addrows(I, column, row, tmp)

    checkRow(getRow(N, 0), getRow(I,0), 0)
    return(I)

##########################################
##          Partial Pivoting            ##
##########################################

def pivotRow(A, col):
    N = copyMatrix(A)
    max_val = 0
    pivrow = 0
    for i in range(col, rows(A)):
        if math.fabs(N[i][col]) > max_val:
            max_val = math.fabs(N[i][col])
            pivrow = i
    return pivrow

def rowReducePartialPivot(M):
    "return row reduced version of M"
    N = copyMatrix(M)
    cs = cols(M)-2  
    rs = rows(M)
    for col in range(cs):
        j = pivotRow(N,col)
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swaprows(N,col,j)
            for row in range(col+1,rs):
                scale = -N[row][col] / N[col][col]
                N=addrows(N, col, row, scale)
    return(N)

##########################################
##      Scaled Partial Pivoting         ##
##########################################

def maxRows(A):
    N = copyMatrix(A)
    max_rows = []

    for i in range(rows(N)):
        max_value = 0
        for j in range(cols(N)-1):
            if math.fabs(N[i][j]) > max_value:
                max_value = math.fabs(N[i][j])
        max_rows.append(max_value)
    return(max_rows)

def findPivotRowScale(column, scales, col):
    scal = 0
    pivrow = -1
    for i in range(col, len(column)):
        tmp = math.fabs(column[i])/scales[i]
        if tmp > scal:
            scal = tmp
            pivrow = i
    return(pivrow)
    
def rowReduceScaledPartialPivot(M):
    " Returns the row reduced M. M is the augmented matrix [A|b]"
    N = copyMatrix(M)
    cs = cols(M)-2  
    rs = rows(M)
    scale_arr = maxRows(N)
    for col in range(cs+1):
        j = findPivotRowScale(getCol(N,col),scale_arr, col)
        temp = scale_arr[j]
        scale_arr[j] = scale_arr[col]
        scale_arr[col] = temp
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swaprows(N,col,j)
            for row in range(col+1,rs):
                multiplier = ((-N[row][col]) / N[col][col])
                N=addrows(N, col, row, multiplier)
    return(N)

def inverseScaled(M):
    N = copyMatrix(M)
    I = identityMatrix(rows(N))
    cs = cols(M)-2  
    rs = rows(M)

    scale_arr = maxRows(N)
    for col in range(cs+1):
        j = findPivotRowScale(getCol(N,col),scale_arr, col)
        temp = scale_arr[j]
        scale_arr[j] = scale_arr[col]
        scale_arr[col] = temp
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swaprows(N,col,j)
                I = swaprows(I, col,j)
            for row in range(col+1,rs):
                multiplier = ((-N[row][col]) / N[col][col])
                N=addrows(N, col, row, multiplier)
                I=addrows(I, col, row, multiplier)
                
    for column in range(1,rs):
        checkRow(getRow(N, column), getRow(I,column), column)
        for row in range((column-1), -1, -1):
            if N[row][column] != 0:
                scale = -1.0 / N[column][column]
                tmp = scale * N[row][column]
                N=addrows(N, column, row, scale * N[row][column])
                I=addrows(I, column, row, tmp)

    checkRow(getRow(N, 0), getRow(I,0), 0)
    return(I)

##########################################
##       Gauss-Seidel Iteration         ##
##########################################

def gaussSeidel(A,B, iter_count):
    """nxn matrix A, B vector, and number of iterations iter_count.

    Prints the X vector at each iteration.
    """
    size = rows(A)
    X_old = [0. for l in range(size)] #initial guess
    X_new = [0. for l in range(size)]
    counter = 0

    while counter < iter_count:
        #print'X({0})'.format(counter)
        #print(X_old)
        for i in range(size):
            sum_upd = 0
            for j in range(i):
                sum_upd += A[i][j]*X_new[j]
            sum_old = 0
            for j in range(i+1, size):
                #print(j)
                sum_old += A[i][j]*X_old[j]
            X_i = (B[i] - sum_upd - sum_old )/A[i][i]
            X_new[i] = X_i
        X_old = copyVec(X_new)
        counter += 1
    return(X_old)

#### Calculating the solution vector

def backSub(M):
    """
    given a row reduced augmented matrix with nonzero 
    diagonal entries, returns a solution vector
    
    """
    cs = cols(M)-1 # cols not counting augmented col
    sol = [0 for i in range(cs)] # place for solution vector
    for i in range(1,cs+1):
        row = cs-i # work backwards
        sol[row] = ((M[row][cs] - sum([M[row][j]*sol[j] for
                    j in range(row+1,cs)])) / M[row][row]) 
    return(sol)


def diag_test(mat):
    """
    Returns True if no diagonal element is zero, False
    otherwise.
    
    """
    for row in range(rows(mat)):
        if mat[row][row] == 0:
            return(False)
    else:
        return(True)

def residualVector(A, x, b):
    """nxn matrix A, approx solution x, and b vector.
    Returns the residual vector of b-Ax.
    """
    Ax = []
    for i in range(rows(A)):
        tmp = dot(getRow(A, i), x)
        Ax.append(tmp)
    Ax = scalarVecMult(-1, Ax)
    return(addVectors(b, Ax))    

##### Main Functions

def ge_1(A,b):    
    """Naive Gaussian Elimination
    Given an nxn matrix A and its b vector, it returns a list. The [0]
    element is the row-reduced augmented matrix, and 
    ge_1(aug)[1] is a solution vector.  The solution
    vector is empty if there is no unique solution.
    
    """
    aug = augment(A,b)
    aug_n = rowReduce(aug)
    if diag_test(aug_n):
        sol = backSub(aug_n)
    else:
        print("\nge_1(): There is no unique solution")
        sol = []
    results = [aug_n, sol]
    return(results)

def ge_2(A,b):
    """Gaussian Elimination w/ Partial Pivoting."""
    aug = augment(A,b)
    aug_n = rowReducePartialPivot(aug)
    if diag_test(aug_n):
        sol = backSub(aug_n)
    else:
        print("\nge_2(): There is no unique solution")
        sol = []
    results = [aug_n, sol]
    return(results)

def ge_3(A,b):
    """Gaussian Elimination w/ Scaled Partial Pivoting."""
    aug = augment(A,b)
    aug_n = rowReduceScaledPartialPivot(aug)
    if diag_test(aug_n):
        sol = backSub(aug_n)
    else:
        print("\nge_3(): There is no unique solution")
        sol = []
    results = [aug_n, sol]
    return(results)

##########################################
##          Hilbert Matrices            ##
##########################################                
def hilbertMatrix(n):
    """Creates an nxn Hilbert Matrix"""
    H = []
    for i in range(1,n+1):
        tmp = []
        for j in range(1, n+1):
            tmp.append(1.0/(i+j-1))
        H.append(tmp)
    return(H)

def BMatrix(n):
    """Creates a vector of the integral of x^k*sin(pi*x) from 0 to 1 where k = 0,...,n-1"""
    B = []
    for i in range(n):
        result = integrate.quad(integ, 0.0, 1.0, args=(i))
        B.append(result[0])
    return(B)
       
def integ(x, n):
    return ((x**n)*(math.sin(math.pi*x)))

def hilbertInverse(n):
    Hinv = []
    for i in range(1,n+1):
        tmp = []
        for j in range(1, n+1):
            val = ((-1)**(i+j))*(i+j-1)*binomial(n+i-1, n-j)*binomial(n+j-1, n-i)*(binomial(i+j-2,i-1)**2)
            tmp.append(val)
        Hinv.append(tmp)
    return(Hinv)

def binomial(n,k):
    if k > n:
        n_ch_k = 0
    prod = 1.0
    for j in range(k):
        prod = prod*(1.0*(n-j)/(k-j))
    return(prod)

"""
A = [[10.0,10.0,10.0,10.0e17,10.0e17],
    [1.0, 10.0e-3, 10.0e-3,10.0e-3,1.0],
    [1.0,1.0,10.0e-3,10.0e-3,2.0],
    [1.0,1.0,1.0,10.0e-3,3.0]]

B = [ [30.00, 591400.00, 591700.00],
      [5.291000, -6.130000, 46.780000]]

C = [ [0.003000, 59.140000, 59.170000], [5.291000,-6.130000, 46.780000]]

D = [ [2.11, -4.21, 0.921, 2.01],
      [4.01, 10.2, 1.12, -3.09],
      [1.09, 0.987, 0.832, 4.21]]

E = [ [1.0,2.0,4.0,3.0,5.0,5.0],
      [3.0,5.0,3.0,1.0,2.0,6.0],
      [1.0,4.0,4.0,2.0,1.0,7.0],
      [4.0,1.0,2.0,5.0,3.0,8.0],
      [5.0,2.0,1.0,4.0,1.0,9.0]]

Ab = [[10.0,10.0,10.0,10.0e17],
    [1.0, 10.0e-3, 10.0e-3,10.0e-3],
    [1.0,1.0,10.0e-3,10.0e-3],
    [1.0,1.0,1.0,10.0e-3]]

F = [ [9,13,7],
      [2,6,8],
      [4,5,10]]
G = [ [2,4,-2],
      [4,9,-3],
      [-2,-3,7]]
H = [[1,3],
     [2,7]]


print(inverseNaive(H))
print(matMult(H,inverseScaled(H)))



print(" ")
print("Naive Gaussian")
#print(ge_1(A))
print(ge_1(A)[1])



print(" ")
print("Partial Pivoting")
#print(ge_2(A))
print(ge_2(A)[1])

print(" ")
print("Scaled Partial Pivoting")
#print(ge_3(A))
print(ge_3(A)[1])

"""
