import math
from math import *

####################
def factoringL(basis,i):
    #L[n][0] is the exponent and L[n][1] is the coefficient of the polynomial
    L=[[0.0,1.0]]
    for n in range(len(basis)):
        if(n != i):
            temp=copyList(L)
            for exponent in range(len(L)):
                L[exponent][0]=L[exponent][0]+1.0
            L.append([0.0,0.0])
            for coefficient in range(len(temp)):
                L[coefficient+1][1] = L[coefficient+1][1] + (-(basis[n])*temp[coefficient][1])
    return(L)

def copyList(M):
    return([[M[row][col] for col in range(2)]for row in
           range(len(M))])

#funtion
def f(x):
    return (1.6*math.exp(-2*x)*math.sin(3*math.pi*x))

def derF(x):
    return (1.6*((-2*math.exp(-2*x)*math.sin(3*math.pi*x))*(3*math.pi*math.exp(-2*x)*math.cos(3*math.pi*x))))

def L(basis, i, x):
    L=1
    C=1
    for n in range(len(basis)):
        if(n != i):
            L=L*(x-basis[n])
            C=C*(basis[i]-basis[n])
    return (L/C)

def derL(basis, i, x):
    C=1.0
    L=factoringL(basis,i)
    
    for n in range(len(basis)):
        if(n != i):
            C=C*1.0*(basis[i]-basis[n])
    #derivative of the polynomial
    for n in range(len(L)):
        L[n][1]=L[n][0]*L[n][1]
        L[n][0]=L[n][0]-1
    L.pop()
    solution=0
    for n in range(len(L)):
        solution= solution + (L[n][1]*(x**L[n][0]))
    return (solution/C)

######################
####Lagrange Interpolating Polynomial

def lagrangeIntPoly(basis, x):
    y=0
    for n in range(len(basis)):
        y=y+f(basis[n])* L(basis, n, x)
    return(y)


######################
####Hermite Interpolating Polynomial

def hermiteIntPoly(basis, x):
    y=0
    for n in range(len(basis)):
        l=L(basis, n, x)
        #h(x)
        y=y+ (f(basis[n])* ( (1-( (2*derL(basis,n,basis[n]) ) * (x-basis[n]))) ) * l**2 )
        #h'(x)
        y=y+ ( derF(basis[n]) * ((x-basis[n])*(l**2)) )
    return(y)


#####################
##for linear splines

def linTheta1(x):
    if(x>=0 and x < 1):
        return (1-x)
    else:
        return 0


def linTheta2(x):
    if(x>=0 and x < 1):
        return (x)
    else:
        return 0

def linearSpline(basis,i,x):
    if(i==0):
       return(linTheta1((x-basis[i]) / (basis[i+1]-basis[i])))
    elif(i==(len(basis)-1)):
        return(linTheta2((x-basis[i-1]) /(basis[i]-basis[i-1])))
    else:
        return(linTheta2((x-basis[i-1]) /(basis[i]-basis[i-1])) + linTheta1((x-basis[i]) / (basis[i+1]-basis[i])))

def linearSplineInter(basis,x):
    solution=0
    for n in range(len(basis)):
        solution=solution+(f(basis[n])* linearSpline(basis,n,x))
    return solution

#####################
##for cubic splines
def cubTheta2(x):
    if(x>=0 and x < 1):
        return (x**2 * (3-2*x))
    else:
        return 0
                
def cubTheta1(x):
    if(x>=0 and x < 1):
        return ((1-x)**2 * ((2*x)+1))
    else:
        return 0

def phi2(x):
    if(x>=0 and x < 1):
        return (x**2 * (x-1))
    else:
        return 0

def phi1(x):
    if(x>=0 and x < 1):
        return -(x*(x-1)**2)
    else:
        return 0
    
def cubicSpline(basis,i,x):
    if(i==0):
       return(cubTheta1((x-basis[i]) / (basis[i+1]-basis[i])))
    elif(i==(len(basis)-1)):
        return(cubTheta2((x-basis[i-1]) /(basis[i]-basis[i-1])))
    else:
        return(cubTheta2((x-basis[i-1]) /(basis[i]-basis[i-1])) + cubTheta1((x-basis[i]) / (basis[i+1]-basis[i])))


def derCubicSpline(basis,i,x):
    if(i==0):
       return((basis[i+1]-basis[i])* phi1((x-basis[i]) / (basis[i+1]-basis[i])))
    elif(i==(len(basis)-1)):
        return((basis[i]-basis[i-1])* phi2((x-basis[i-1]) /(basis[i]-basis[i-1])))
    else:
        return((basis[i]-basis[i-1])* phi2((x-basis[i-1]) /(basis[i]-basis[i-1])) + (basis[i+1]-basis[i])*phi1((x-basis[i]) / (basis[i+1]-basis[i])))



def cubSplineInter(basis,x):
    solution=0
    for n in range(len(basis)):
        solution=solution+(f(basis[n])* cubicSpline(basis,n,x)) + (derF(basis[n])* derCubicSpline(basis,n,x))
    return solution

###############
#test
'''
p=1.0/10
X=[0, 1.0/6, 1.0/3, 1.0/2, 7.0/12, 2.0/3, 3.0/4, 5.0/6, 11.0/12, 1]


print("Using Lagrange Interpolating Polynomials")
print(lagrangeIntPoly(X,p))
print("Using Hermite Interpolating Polynomials")
print(hermiteIntPoly(X,p))
print("Using Linear spline interpolants")
print(linearSplineInter(X,p))
print("Using cublic spline interpolants")
print(cubSplineInter(X,p))


'''


