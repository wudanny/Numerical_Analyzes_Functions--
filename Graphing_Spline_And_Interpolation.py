
import Interpolating_And_Spline
import math
from math import *
import matplotlib.pyplot as plt
from numpy import arange

# setup domain data
delta_theta = 1.0/36
thetas = arange(0,1+delta_theta, delta_theta)


def graphLag(xCoord):
    # setup output for both functions
    interpolant = [homework2.lagrangeIntPoly(xCoord, t) for t in thetas]
    function = [Interpolating_And_Spline.f(t) for t in thetas]

    # plot the two graphs
    plt.plot(thetas, interpolant, 'red', label=r'$Lagrange Interpolating Poly$')
    plt.plot(thetas, function, color='blue', label= r'$f(x)$')



    # place a title above plot
    plt.title(r'$Lagrange Interpolating Poly$ vs $f(x)$')

    # make horizonal axis run from -pi to 2*pi while vertical axis goes from -2 to 2
    plt.axis([0, 1 , -(1.6/e), 1.6/e**(1.0/3)], 'equal')

    # put horizontal and vertical lines because matplotlib puts axes around plot
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')

    # put a legend of the two plots in lower right hand corner
    plt.legend(loc='upper right')



    # save plot to a file or show plot from IDLE shell
    plt.show()
    return

def graphHermite(xCoord):
    # setup output for both functions
    interpolant = [homework2.hermiteIntPoly(xCoord, t) for t in thetas]
    function = [Interpolating_And_Spline.f(t) for t in thetas]

    # plot the two graphs
    plt.plot(thetas, interpolant, 'red', label=r'$Hermite Interpolating Poly$')
    plt.plot(thetas, function, color='blue', label= r'$f(x)$')


    # place a title above plot
    plt.title(r'$Hermite Interpolating Poly$ vs $f(x)$')

    # make horizonal axis run from -pi to 2*pi while vertical axis goes from -2 to 2
    plt.axis([0, 1 , -(1.6/e), 1.6/e**(1.0/3)], 'equal')

    # put horizontal and vertical lines because matplotlib puts axes around plot
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')

    # put a legend of the two plots in lower right hand corner
    plt.legend(loc='upper right')



    # save plot to a file or show plot from IDLE shell
    plt.show()
    return

def graphLinearSpline(xCoord):

    # setup output for both functions
    interpolant = [homework2.linearSplineInter(xCoord, t) for t in thetas]
    function = [Interpolating_And_Spline.f(t) for t in thetas]

    # plot the two graphs
    plt.plot(thetas, interpolant, 'red', label=r'$Linear Spline$')
    plt.plot(thetas, function, color='blue', label= r'$f(x)$')


    # place a title above plot
    plt.title(r'$Linear Spline$ vs $f(x)$')

    # make horizonal axis run from -pi to 2*pi while vertical axis goes from -2 to 2
    plt.axis([0, 1 , -(1.6/e), 1.6/e**(1.0/3)], 'equal')

    # put horizontal and vertical lines because matplotlib puts axes around plot
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')

    # put a legend of the two plots in lower right hand corner
    plt.legend(loc='upper right')



    # save plot to a file or show plot from IDLE shell
    plt.show()
    return

def graphCubicSpline(): 

    # setup output for both functions
    interpolant = [homework2.cubSplineInter(X, t) for t in thetas]
    function = [homework2.f(t) for t in thetas]

    # plot the two graphs
    plt.plot(thetas, interpolant, 'red', label=r'$Cubic Spline$')
    plt.plot(thetas, function, color='blue', label= r'$f(x)$')


    # place a title above plot
    plt.title(r'$Cubic Spline$ vs $f(x)$')

    # make horizonal axis run from -pi to 2*pi while vertical axis goes from -2 to 2
    plt.axis([0, 1 , -(1.6/e), 1.6/e**(1.0/3)], 'equal')

    # put horizontal and vertical lines because matplotlib puts axes around plot
    plt.axhline(y=0, color='black')
    plt.axvline(x=0, color='black')

    # put a legend of the two plots in lower right hand corner
    plt.legend(loc='upper right')



    # save plot to a file or show plot from IDLE shell
    plt.show()
    return


    
