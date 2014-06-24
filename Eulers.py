import math
from math import *
import matplotlib.pyplot as plt
from numpy import arange

def function(t,x):
    return 1.0*t-2.0*x

def eulers(h,t0,x0,end):
    solutions=[]
    t=1.0*t0
    x=1.0*x0
    solutions.append([t,x])
    while (t<1.0*end):
        x=x+h*function(t,x)
        t=t+h
        solutions.append([t,x])
    return solutions

def modifiedEulers(h,t0,x0,end):
    solutions=[]
    t=1.0*t0
    x=1.0*x0
    solutions.append([t,x])
    while (t<end):
        m1=function(t,x)
        xstar=x+h*m1
        t=t+h
        m2=function(t,xstar)
        x=x+h*((m1+m2)/2)
        solutions.append([t,x])
    return solutions

def rk4(h,t0,x0,end):
    solutions=[]
    t=1.0*t0
    x=1.0*x0
    solutions.append([t,x])
    while (t<end):
        k1=h*function(t,x)
        k2=h*function(t+h/2,x+k1/2)
        k3=h*function(t+h/2,x+k2/2)
        k4=h*function(t+h,x+k3)
        t=t+h
        x=x+((k1+2*k2+2*k3+k4)/6)
        solutions.append([t,x])
    return solutions

def exactSolution(t):
    return (5.0/4)*e**(-2*t)+(1.0/4)*(2*t-1)

'''
    "Test"
step=.05

E=eulers(step,0,1,2)
ME=modifiedEulers(step,0,1,2)
R=rk4(step,0,1,2)

print("Euler Method with h="),
print(step)
print("t"),
print("      |     "),
print("x")
for i in range(len(E)):
    print(E[i][0]),
    print("    |     "),
    print(E[i][1]),
    print("    |     "),
    print(E[i][1]-exactSolution(E[i][0]))

print
print("Modified Euler Method with h="),
print(step)
print("t"),
print("      |     "),
print("x")
for i in range(len(ME)):
    print(ME[i][0]),
    print("    |     "),
    print(ME[i][1]),
    print("    |     "),
    print(ME[i][1]-exactSolution(ME[i][0]))

print
print("RK-4 with h="),
print(step)
print("t"),
print("      |     "),
print("x")
for i in range(len(R)):
    print(R[i][0]),
    print("    |     "),
    print(R[i][1]),
    print("    |     "),
    print(R[i][1]-exactSolution(R[i][0]))




'''












    
