import numpy as np
from numpy import pi, cos, sin
from numpy.linalg.linalg import norm, inv, cond, solve, svd
import matplotlib.pyplot as plt

import dmsuite as dm
# from cheb_trefethen import *
# from noSlipBC import *
from get_F import *
from OSQoperators import *

n = 201
kx = 1.
kz = 6.
omega = 5.
Re = 1005.5
# Re = 2389.56

y, D4 =  dm.cheb4c(n)           # fourth derivative with clamped BC (n-1 X n-1)
y96, D96 = dm.chebdif(96, 1)    # first derivative 
y, D = dm.chebdif(n, 2)         # first two derivatives (2 X n+1 X n+1)
D1 = D[0,:,:]                   # first derivative (n+1 X n+1)
D2 = D[1,:,:]                   # second derivative (n+1 X n+1)

"Velocity field and its derivatives============================================"
U96 = np.loadtxt('turbstats_Re_30000.asc', skiprows=3)
U96 = U96[:, 6]
U = np.polyfit(y96, U96, 30)
U = np.poly1d(U)
U = U(y)
Uy = D1 @ U
Uyy = D2 @ U

"Apply homogeneous BCs========================================================="
y = y[1:n]
D1 = D1[1:n,1:n]
D2 = D2[1:n,1:n]
U = U[1:n]
Uy = Uy[1:n]
Uyy = Uyy[1:n]

"Build operators==============================================================="
M, L, B = OSQoperators(n, Re, kx, kz, U, Uy, Uyy, D1, D2, D4)


"Build weight matrix to enforce energy norm of the resolvent==================="
# ksq = kx*kx + kz*kz
# y, w = clencurt(n)
# wdiag = np.diag(w)
# F = get_F(D1, wdiag, ksq, n)
# F = F[1:n,1:n]
# Finv = inv(F)

"Build a matrix to impose BCs on the RHS======================================="
Id2n = np.eye(2*n-2)
G = Id2n
G[0,0] = 0.
G[n-2,n-2] = 0.

"Get the resolvent and take its SVD============================================"
H = inv(-1j*omega*Id2n + inv(M)@L)@G

Psi, Sigma, PhiH = svd(H)

"""
Check SVD and othogonality
"""
print("check ", np.allclose(H, Psi@diag(Sigma)@PhiH))
print("check ", np.allclose(Id2n, np.transpose(np.conj(Psi)) @ Psi, rtol=1e-8))
print("check ", np.allclose(np.transpose(np.conj(Psi))@Psi, Psi@np.conj(np.transpose(Psi)), rtol=1e-8))

"Plotting======================================================================"

e0 = np.abs(Psi[0:n-1, 0])
e1 = np.abs(Psi[0:n-1, 1])
e2 = np.abs(Psi[0:n-1, 2])
e3 = np.abs(Psi[0:n-1, 3])
e4 = np.abs(Psi[0:n-1, 4])
e5 = np.abs(Psi[0:n-1, 5])

fig, axs = plt.subplots(6)
fig.suptitle('Singular response modes')
axs[0].plot(y, np.real(e0))
axs[1].plot(y, np.real(e1))
axs[2].plot(y, np.real(e2))
axs[3].plot(y, np.real(e3))
axs[4].plot(y, np.real(e4))
axs[5].plot(y, np.real(e5))
plt.show()

counts = linspace(1, 20, 20)
fig = plt.figure()
plt.semilogy(counts, Sigma[0:20],'.')
plt.title('Singular values '+str(kx)+'-'+str(kz)+'-'+str(omega))
plt.xlabel("index")
plt.ylabel("sigma")
fig.savefig('images/'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_singular_values.png')
plt.show()

# f = cos(pi*y)
# f4 = (1/pi**4)*D4@f
# fig, axs = plt.subplots(2)
# fig.suptitle('Differentiation test')
# axs[0].plot(y, f)
# axs[1].plot(y[50:n-50], f4[50:n-50])
# plt.show()

# "Chebyshev grid and differentiation matrices"
# D, y96 = cheb(96)


# D1 = D[0,:,:]                   # first derivative
# D2 = D[1,:,:]                   # second derivative



# Uy = D @ U
# Uyy = D @ Uy

