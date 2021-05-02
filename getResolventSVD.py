import numpy as np

from numpy import pi, cos, sin
from numpy.linalg.linalg import norm, inv
from cheb_trefethen import *

from clencurt import *
from get_F import *

def getResolventSVD(Re, kx, kz, omega, Ny, U):
    
    ksq = kx*kx + kz*kz
    n = Ny+1
    Reinv = 1.0/Re

    "Chebyshev grid and differentiation matrices"
    D, y = cheb(Ny)
    D2 = D@D
    D4 = D2@D2
    Z = np.zeros((n, n))
    Idn = np.eye(n)
    Id2n = np.eye(2*n)

    "Velocity field and its derivatives"
    # U = 1. - y**2

    Uy = D @ U
    Uyy = D @ Uy

    "Requirements for F"
    y, w = clencurt(Ny)
    wdiag = np.diag(w)
    
    """
    B L M matrices
    """
    M = np.block([[ksq*Idn - D2, Z], [Z, Idn]])

    Los = 1j*kx*(U)*(ksq*Idn - D2) + 1j*kx*Uyy + Reinv*(ksq**2 + D4 - 2*ksq*D2)
    Lsq = (1j*kx)*(U*Idn) + Reinv*(ksq*Idn - D2)
    L = np.block([[Los, Z],[(1j*kz)*(Uy*Idn), Lsq]])

    B = np.block([[-1j*kx*D, -ksq*Idn, -1j*kz*D],[1j*kz*Idn, Z, -1j*kx*Idn]])

    """
    Create the LHS operator in the equation
    (-i*om*M + L)[v^, eta^] = B[fu^, fv^, fw^]
    """
    LHS = -1j*omega*M + L

    """
    Apply boundary conditions on LHS and B operators
    v = dv/dy = eta = 0 at the walls
    """
    "v = 0 at the walls"
    LHS[0, :] = 0.
    LHS[0, 0] = 1.
    LHS[n-1, :] = 0.
    LHS[n-1, n-1] = 1.

    M[0, :] = 0.
    M[n-1, :] = 0.

    "dv/dy = 0 at the walls"
    LHS[1, :] = 0.
    LHS[1, 0:n] = D[0, :]
    LHS[n-2, :] = 0.
    LHS[n-2, 0:n] = D[n-1, :]

    M[1, :] = 0.
    M[n-2, :] = 0.

    "eta = 0 at the walls"
    LHS[n, :] = 0.
    LHS[n, n] = 1.
    LHS[2*n-1, :] = 0.
    LHS[2*n-1, 2*n-1] = 1.

    M[n, :] = 0
    M[2*n-1, :] = 0
    
    RHS = np.ones(2*n)
    # RHS[0] = 0
    # RHS[1] = 0
    # RHS[2*n-2] = 0
    # RHS[2*n-1] = 0
    """
    Get the Resolvent matrix defined as inverse of LHS
    LHS = -i om M + L
    Obtain H as a solution of LHS x = B y with y = ones(2*n)
       LHS x = B y
    => x = LHS^-1 B
    So obtain it by solving a linear system 
    Take the SVD of H
    """
    
    H = np.linalg.solve(LHS, M)
    F = get_F(D, wdiag, ksq, Ny)
    Finv = inv(F)

    Hs = F@H@Finv

    Psi, Sigma, PhiH = np.linalg.svd(Hs, full_matrices=False)

    return Psi, Sigma, PhiH