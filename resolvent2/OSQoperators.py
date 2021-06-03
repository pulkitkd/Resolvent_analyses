import numpy as np
from numpy import diag, sin, cos

def OSQoperators(n, Re, kx, kz, U, Uy, Uyy, D1, D2, D4):
    Reinv = 1.0/Re
    ksq = kx*kx + kz*kz
    Z = np.zeros((n-1, n-1))
    Idn = np.eye(n-1)
    U = diag(U)
    Uy = diag(Uy)
    Uyy = diag(Uyy)
    
    M = np.block([[ksq*Idn - D2, Z], [Z, Idn]])

    Los = 1j*kx*U@(ksq*Idn - D2) + 1j*kx*Uyy + Reinv*(D4 - (2*ksq)*D2 + (ksq**2)*Idn)
    Lsq = (1j*kx)*U + Reinv*(ksq*Idn - D2)
    L = np.block([[Los, Z],[(1j*kz)*Uy, Lsq]])

    B = np.block([[-1j*kx*D1, -ksq*Idn, -1j*kz*D1],[1j*kz*Idn, Z, -1j*kx*Idn]])
    
    return M, L, B
