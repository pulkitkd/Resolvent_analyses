from numpy import pi, cos, arange, ones, zeros
import numpy as np
##############################
def clencurt(N):
    '''
    Clenshaw-Curtis quadrature
    Python version of Professor Lloyd Trefethen's matlab code
    Obtained at http://people.maths.ox.ac.uk/trefethen/clencurt.m
    author: Pratik Aghor
    '''
    theta = pi*arange(0, N+1)/N
    xc  = cos(theta) # xc = Chebyshev-Gauss-Lobatto grid

    w = zeros(N+1)
    ii = arange(1, N)
    v = ones(N-1)

    if(int(N)%2 == 0):
        w[0] = 1/((N**2)-1);
        w[N] = w[0];
        for k in range(1, int(N/2) ):
            v = v - 2*cos(2*k*theta[ii])/((4*k**2)-1)
        v = v - cos(N*theta[ii])/((N**2)-1);
    else:
        w[0] = 1/(N**2);
        w[N] = w[0];
        for k in range(1, int((N-1)/2) + 1):
            v = v - 2*cos(2*k*theta[ii])/((4*k**2)-1)

    w[ii] = 2*v/N;

    return xc, w
##############################
