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
    deltaInv = np.linalg.inv(ksq*Idn - D2)
    
    M = np.block([[ksq*Idn - D2, Z], [Z, Idn]])

    Los = 1j*kx*U@(ksq*Idn - D2) + 1j*kx*Uyy + Reinv*(D4 - (2*ksq)*D2 + (ksq**2)*Idn)
    Lsq = (1j*kx)*U + Reinv*(ksq*Idn - D2)
    # L11 = deltaInv @ Los
    L11 = Los
    L12 = Z
    L21 = 1j*kz*Uy
    L22 = Lsq
    
    L = np.block([[L11, L12],[L21, L22]])

    B = np.block([[-1j*kx*D1, -ksq*Idn, -1j*kz*D1],[1j*kz*Idn, Z, -1j*kx*Idn]])
    
    C  = np.block([[1j*kx*D1, -1j*kz*Idn], [ksq*Idn, Z], [1j*kz*D1, 1j*kx*Idn]])
    C = (1/ksq) * C
    
    return M, L, B, C


# """
# M and L matrices
# """
# Los = 1j*kx*diag(U)@(ksq*Idn - D2) + 1j*kx*diag(Uyy) + Reinv*(D4 - 2*ksq*D2 + ksq*ksq*Idn)
# Lsq = 1j*kx*diag(U) + Reinv*(ksq*Idn - D2)
# L = np.block([[Los, np.zeros((n-1, n-1))],[1j*kz*(diag(Uy)), Lsq]])

# M = np.block([[(ksq*Idn - D2), np.zeros((n-1, n-1))], [np.zeros((n-1, n-1)), Idn ]])