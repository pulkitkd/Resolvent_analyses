import numpy as np
from scipy.linalg import eig

import matplotlib.pyplot as plt
from numpy import pi, cos, sin

from cheb_trefethen import *

"""
Inputs
"""
Re = 5000.
Reinv = 1.0/Re
Ny = 128
kx = 2.
kz = 0.
omega = 1.

ksq = kx*kx + kz*kz
n = Ny+1

"Chebyshev grid and differentiation matrices"
D, y = cheb(Ny)
D2 = D@D
D4 = D2@D2
Z = np.zeros((n, n))
Idn = np.eye(n)
Id2n = np.eye(2*n)

"Velocity field and its derivatives"
U = 1. - y*y
Uy = D @ U
Uyy = D @ Uy

"""
M and L matrices
"""
M = np.block([[ksq - D2, Z], [Z, Idn]])

Los = 1j*kx*(U)*(ksq*Idn - D2) + 1j*kx*Uyy + Reinv*(ksq**2 + D4 - 2*ksq*D2)
Lsq = (1j*kx)*(U*Idn) + Reinv*(ksq*Idn - D2)
L = np.block([[Los, Z],[(1j*kz)*(Uy*Idn), Lsq]])

"""
Boundary Conditions
v = dv/dy = eta = 0 at the walls
"""
"v = 0 at the walls"
L[0, :] = 0.
L[n-1, :] = 0.
L[0, 0] = 1.
L[n-1, n-1] = 1.


"dv/dy = 0 at the walls"
L[1, :] = 0.
L[n-2, :] = 0.
L[1, 0:n] = D[0, :]
L[n-2, 0:n] = D[n-1, :]

"eta = 0 at the walls"
L[n :] = 0.
L[2*n-1, :] = 0.
L[n, n] = 1.
L[2*n-1, 2*n-1] = 1.

"v = 0 at the walls"
M[0, :] = 0.
M[n-1, :] = 0.

"dv/dy = 0 at the walls"
M[1, :] = 0.
M[n-2, :] = 0.

"eta = 0 at the walls"
M[n, :] = 0.
M[2*n-1, :] = 0.


"""
Solve the Orr-Sommerfeld problem (a useful check on the BCs and operators)
   L x = i omega M x
=> Minv L x = i omega x
=> i omega = eigval(Minv*L) 
"""

"""
Get the eigenvalues with largest imaginary parts
"""
MinvL = np.linalg.lstsq(M, L)[0]
w, v = eig(MinvL, check_finite=True)
w = -1j*w
w = np.sort_complex(w)

"""
Plot the eigenvalues
"""
plt.scatter(np.real(w), np.imag(w))
plt.xlim((0,1))
plt.ylim((-1,0))
plt.xlabel("real")
plt.ylabel("imag")
plt.show()

"""
Get the Resolvent matrix defined as inverse of (-i omega + Minv L)
Let LHS = -i omega + Minv L
Obtain H as a psuedoinverse of LHS
Take the SVD of H
"""
Minv = np.linalg.pinv(M)
MinvL = Minv @ L

LHS = -1j*omega*Id2n + MinvL
# H = np.linalg.lstsq(LHS, np.eye(2*n))[0]
H = np.linalg.pinv(LHS)
Psi, Sigma, PhiH = np.linalg.svd(H)

"Plotting"
counts = linspace(1, 2*n, 2*n)
plt.scatter(counts,Sigma)
plt.xlabel("index")
plt.ylabel("sigma")
plt.show()

plt.plot(y, Psi[0:n, 0])
plt.xlabel("y")
plt.ylabel("Psi")
plt.show()

plt.plot(y,PhiH[0, 0:n])
plt.xlabel("y")
plt.ylabel("Phi")
plt.show()