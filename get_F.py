from numpy import pi, cos, arange, ones, zeros, diag, matmul, transpose
import numpy as np
from numpy.linalg import cholesky

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
def get_F(D1, wdiag, ksq, Ny):
    # get Chebyshev-Gauss-Lobatto grid and
    # weights for Clenshaw-Curtis quadrature
    ny = Ny + 1
    # formulate wdiag to be the right size of M = 2n x 2n
    Z = zeros((ny, ny))

    # wdiag = diag(w)

    Wv = zeros((ny, ny))
    Weta = zeros((ny, ny))

    # k^2 - D^2 = (D^{T}D + k^2)
    # if first row, first col, last row, last col cropped
    Wv = ((transpose(D1) @ wdiag @ D1)/ksq + wdiag)
    # Wv = Wv[1:ny-1, 1:ny-1] # crop first row, first col, last row, last col

    Weta = (wdiag)/ksq

    # Zcropped = zeros((ny-2, ny-2))
    Zcropped = Z

    W = np.block([ \
    [Wv, Zcropped],\
    [Zcropped, Weta] \
    ])
    # print("shape(W) = ", np.shape(W))
    # print("W = ", W)

    F = cholesky(W)
    F = F.T
    return F
##############################
