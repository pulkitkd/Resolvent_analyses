import numpy as np

from numpy import pi, cos, sin
from numpy.linalg.linalg import norm, inv
from cheb_trefethen import *

from clencurt import *
from noSlipBC import *
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
    M, L and B matrices
    """
    M = np.block([[ksq*Idn - D2, Z], [Z, Idn]])

    Los = 1j*kx*diag(U)@(ksq*Idn - D2) + 1j*kx*(diag(Uyy)) + Reinv*(ksq*Idn - D2)@(ksq*Idn - D2)
    Lsq = (1j*kx)*(diag(U)) + Reinv*(ksq*Idn - D2)
    L = np.block([[Los, Z],[(1j*kz)*(diag(Uy)), Lsq]])

    B = np.block([[-1j*kx*D, -ksq*Idn, -1j*kz*D],[1j*kz*Idn, Z, -1j*kx*Idn]])
    
    
    """
    Apply boundary conditions on M and L operators
    v = dv/dy = eta = 0 at the walls
    
    BCs are applied in a way to ensure invertibility of M
    """
    M, L = noSlipBC(M,L,n)
    
    """
    Create the LHS operator in the equation
    (-i*om + Minv*L)[v^, eta^] = Minv*B[fu^, fv^, fw^]
    """
    L1 = inv(M) @ L
    LHS = -1j*omega*Id2n + L1

    """
    Get the Resolvent matrix defined as inverse of LHS
    """
    
    # H = np.linalg.solve(LHS, Id2n)
    H = inv(LHS)
    
    # Scale the resolvent
    F = get_F(D, wdiag, ksq, Ny)
    Finv = inv(F)
    Hs = F@H@Finv

    Psi, Sigma, PhiH = np.linalg.svd(Hs, full_matrices=False)
    Psi = Finv @ Psi
    Phi = np.conj(np.transpose(PhiH))
    Phi = Finv @ Phi

    return Psi, Sigma, PhiH