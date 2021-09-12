from numpy import pi, cos, arange, ones, zeros, diag, matmul, transpose
import numpy as np
from numpy.linalg import cholesky

import dmsuite as dm
from clencurt import *
from cheb_trefethen import *
##############################
"""
Form the matrix W as an augmented form
of the energy weight matrix (M/k^{2})
(where M comes from the OS equation,
M\dot{q} = L q is the linearized system).
From W, get F by Cholesky decomp,
W = F'*F.

author: Pratik Aghor
"""
##############################
def get_F(ksq, n):
    # formulate wdiag to be the right size of M = 2n-2 x 2n-2
    
    Z = zeros((n-1, n-1))
    y, D = dm.chebdif(n, 1)         # first derivative matrix (1 x n+1 x n+1)
    D1 = D[0,:,:]                    # n+1 x n+1


    # get Chebyshev-Gauss-Lobatto grid and
    # weights for Clenshaw-Curtis quadrature
    y, w = clencurt(n)               # w is a n+1-element array
    wdiag = np.diag(w)               # wdiag is n+1 x n+1
    
    Wv = zeros((n+1, n+1))
    Weta = zeros((n+1, n+1))

    # k^2 - D^2 = (D^{T}D + k^2)
    # if first row, first col, last row, last col cropped
    Wv = ((transpose(D1) @ wdiag @ D1)/ksq + wdiag)
    Wv = Wv[1:n, 1:n] # Wv is n-1 x n-1

    Weta = (wdiag)/ksq
    Weta = Weta[1:n, 1:n] # Weta is n-1 x n-1

    W = np.block([ \
    [Wv, Z],\
    [Z, Weta] \
    ])
    print("shape(W) = ", np.shape(W))
    # print("W = ", W)

    F = cholesky(W)
    F = F.T
    # print("shape F = ", shape(F))
    return W, F
##############################
