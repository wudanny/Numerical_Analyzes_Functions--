import math
from math import *


def function(x):
    return ((math.e)**(-x**2/2))

def centeredDifference(x,h):
    i=function(x+h)-function(x-h)
    return(i/(2*h))

def realDer(x):
    return (-x*(math.e)**(-x**2/2))

def error(x):
    return realDer(1.4)-x;

'''
"Example"
H=[10**-h for h in range(21) ]
C=[centeredDifference(1.4,p) for p in H]
E=[error(b) for b in C]

for z in range(21):
    print(H[z]),
    print("         "),
    print(C[z]),
    print("         "),
    print(E[z]),
    print("         ")
    
'''

