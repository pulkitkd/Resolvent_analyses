import numpy as np
from scipy.linalg import eig

import matplotlib.pyplot as plt
import dmsuite as dm
from OSQoperators import *
from numpy import pi, cos, sin, diag

from cheb_trefethen import *
from noSlipBC import *
"""
Inputs
"""
n = 200
kx = 1.
kz = 0.

Re = 4000.
# Re = 1005.5
# Re = 2389.56

y, D4 =  dm.cheb4c(n)           # fourth derivative with clamped BC (n-1 X n-1)
y, D = dm.chebdif(n, 2)         # first two derivatives (2 X n+1 X n+1)
D1 = D[0,:,:]                   # first derivative (n+1 X n+1)
D2 = D[1,:,:]                   # second derivative (n+1 X n+1)

U = 1. - y*y
Uy = D1 @ U
Uyy = D2 @ U

"Apply homogeneous BCs========================================================="
y = y[1:n]
D1 = D1[1:n,1:n]
D2 = D2[1:n,1:n]
# print("size of D2 = ",np.shape(D2))

U = U[1:n]
Uy = Uy[1:n]
Uyy = Uyy[1:n]

ksq = kx**2 + kz**2
Reinv = 1 / Re
Idn = np.eye(n-1)

"""
M and L matrices
"""
Los = 1j*kx*diag(U)@(ksq*Idn - D2) + 1j*kx*diag(D2 @ U) + Reinv*(D4 - 2*ksq*D2 + ksq*ksq*Idn)
Lsq = 1j*kx*diag(U) + Reinv*(ksq*Idn - D2)
L = np.block([[Los, np.zeros((n-1, n-1))],[1j*kz*(diag(D1 @ U)), Lsq]])

M = np.block([[(ksq*Idn - D2), np.zeros((n-1, n-1))], [np.zeros((n-1, n-1)), Idn ]])

"Create operators=============================================================="
# M, L, B, C = OSQoperators(n, Re, kx, kz, U, Uy, Uyy, D1, D2, D4)

"""
Solve the generalized eigenvalue problem L x = i omega M x
"""
w, v = eig(L, b=M, check_finite=True)

"""
Get the eigenvalues with largest imaginary parts
"""
w = -1j*w
w = np.sort_complex(w)
print("Eigenvalues = ", w)

# print("Eigenvects = ", v)


e0 = np.abs(v[0:n, 0])
e1 = np.abs(v[0:n, 1])
e2 = np.abs(v[0:n, 2])
e3 = np.abs(v[0:n, 3])
e4 = np.abs(v[0:n, 4])
e5 = np.abs(v[0:n, 5])

"""
Plot the eigenvalues
"""
# plt.plot(real(w),'.')
# plt.ylim((-10,100))
plt.scatter(np.real(w), np.imag(w))
plt.xlim((0,1))
plt.ylim((-1,0.5))
plt.xlabel("real")
plt.ylabel("imag")
plt.grid()
plt.show()

# # plot the eigvects
# fig, axs = plt.subplots(6)
# fig.suptitle('Singular response modes')
# axs[0].plot(y, np.real(e0))
# axs[1].plot(y, np.real(e1))
# axs[2].plot(y, np.real(e2))
# axs[3].plot(y, np.real(e3))
# axs[4].plot(y, np.real(e4))
# axs[5].plot(y, np.real(e5))
# plt.show()


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

