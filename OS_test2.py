import numpy as np
from scipy.linalg import eig

import matplotlib.pyplot as plt
from numpy import pi, cos, sin, diag

from cheb_trefethen import *
from noSlipBC import *
"""
Inputs
"""
Re = 8000.
Reinv = 1.0/Re
Ny = 128
kx = 1.
kz = 0.

ksq = kx*kx + kz*kz
n = Ny+1
D, y = cheb(Ny)
D2 = D@D
D4 = D2@D2
Z = np.zeros((n, n))
Idn = np.eye(n)
U = 1. - y*y
Uy = D @ U
Uyy = D2 @ U

"""
M and L matrices
"""
Los = 1j*kx*diag(U)@(ksq*Idn - D2) + 1j*kx*diag(Uyy) + Reinv*(D4 - 2*ksq*D2 + ksq*ksq*Idn)
Lsq = 1j*kx*diag(U) + Reinv*(ksq*Idn - D2)
L = np.block([[Los, np.zeros((n, n))],[1j*kz*(diag(Uy)), Lsq]])

M = np.block([[(ksq*Idn - D2), np.zeros((n, n))], [np.zeros((n, n)), Idn ]])

# Boundary conditions
M, L = noSlipBC(M,L,n)

"""
Solve the generalized eigenvalue problem L x = i omega M x
"""
w, v = eig(L, b=1j*M, check_finite=True)

"""
Get the eigenvalues with largest imaginary parts
"""
w = w
w = np.sort_complex(w)
# print("Eigenvalues = ", w)
print("Eigenvects = ", v)


e0 = np.abs(v[0:n, 0])
e1 = np.abs(v[0:n, 1])
e2 = np.abs(v[0:n, 2])
e3 = np.abs(v[0:n, 3])
e4 = np.abs(v[0:n, 4])
e5 = np.abs(v[0:n, 5])

"""
Plot the eigenvalues
"""
plt.scatter(np.real(w), np.imag(w))
plt.xlim((0,1))
plt.ylim((-1,0.5))
plt.xlabel("real")
plt.ylabel("imag")
plt.grid()
plt.show()

# plo the eigvects
fig, axs = plt.subplots(6)
fig.suptitle('Singular response modes')
axs[0].plot(y, np.real(e0))
axs[1].plot(y, np.real(e1))
axs[2].plot(y, np.real(e2))
axs[3].plot(y, np.real(e3))
axs[4].plot(y, np.real(e4))
axs[5].plot(y, np.real(e5))
plt.show()


"---------------------------Using Matrix Inversion-----------------------------"
"""
Solve the generalized eigenvalue problem M x = i omega L x
"""
Minv = np.linalg.inv(M)
MinvL = Minv@L
w, v = eig(MinvL, check_finite=True)

"""
Get the eigenvalues with largest imaginary parts
"""

w = -1j*w
w = np.sort(w)
# print("Eigenvalues = ", w)

np.savetxt('OS_eigenvalues.dat', (np.real(w), np.imag(w)), delimiter=',')

"""
Plot the eigenvalues
"""
plt.scatter(np.real(w), np.imag(w))
plt.xlim((0,1))
plt.ylim((-1,0.5))
plt.xlabel("real")
plt.ylabel("imag")
plt.grid()
plt.show()

