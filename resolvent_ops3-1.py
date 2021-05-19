import numpy as np
# from scipy.linalg import eig, svd, pinv2
# from scipy.sparse.linalg import inv

import matplotlib.pyplot as plt
from numpy import pi, cos, sin
from numpy.linalg.linalg import norm, inv, cond, solve

from cheb_trefethen import *
from noSlipBC import *
from get_F import *

def resolventNorm(x, y, M, L, F):
    L1 =  inv(M)@L
    omega = x + 1j*y
    LHS = -1j*omega*Id2n + L1
    H = inv(LHS)
    Finv = inv(F)
    Hs = F @ H @ Finv
    return norm(Hs)


"""
Inputs
"""
Ny = 200
Re = 1000
# Re = 1005.5
# Re = 2389.55836707658146
# lxp = 700
# lzp = 100
# c = 10

# kx = 2*pi*Re/lxp
# kz = 2*pi*Re/lzp
# omega = c/kx

kx = 1.
kz = 0.
omega = 0.05

Reinv = 1.0/Re
ksq = kx*kx + kz*kz
n = Ny+1

"Chebyshev grid and differentiation matrices"
D, y96 = cheb(96)
D, y = cheb(Ny)
D2 = D@D
D4 = D2@D2
Z = np.zeros((n, n))
Idn = np.eye(n)
Id2n = np.eye(2*n)

"Velocity field and its derivatives"
# U96 = np.loadtxt('turbstats_Re_60000.asc', skiprows=3)
# U96 = U96[:, 6]
# U = np.polyfit(y96, U96, 30)
# U = np.poly1d(U)
# U = U(y)
U = 1 - y**2

Uy = D @ U
Uyy = D @ Uy


"""
M, L and B matrices
"""
M = np.block([[ksq*Idn - D2, Z], [Z, Idn]])

Los = 1j*kx*diag(U)@(ksq*Idn - D2) + 1j*kx*(diag(Uyy)) + Reinv*(ksq*Idn - D2) @ (ksq*Idn - D2)
Lsq = (1j*kx)*(diag(U)) + Reinv*(ksq*Idn - D2)
L = np.block([[Los, Z],[(1j*kz)*(diag(Uy)), Lsq]])

B = np.block([[-1j*kx*D, -ksq*Idn, -1j*kz*D],[1j*kz*Idn, Z, -1j*kx*Idn]])



"""
Apply boundary conditions on LHS and B operators
v = dv/dy = eta = 0 at the walls
"""
M, L = noSlipBC(M,L,n)



"""
Create the LHS operator in the equation
(-i*om*Id2n + Minv*L)[v^, eta^] = Minv*B[fu^, fv^, fw^]
"""
L1 =  inv(M)@L
LHS = -1j*omega*Id2n + L1

print("cond(LHS) = ", cond(LHS))
H = inv(LHS)

"""
Scale the resolevnt operator
"""
y, w = clencurt(Ny)
wdiag = np.diag(w)
F = get_F(D, wdiag, ksq, Ny)
Finv = inv(F)
Hs = F @ H @ Finv

Psi, Sigma, PhiH = np.linalg.svd(Hs, full_matrices=False)

"""
Check SVD and othogonality
"""
print("check ", np.allclose(Hs, Psi@diag(Sigma)@PhiH))
print("check ", np.allclose(Id2n, np.transpose(np.conj(Psi)) @ Psi, rtol=1e-8))
print("check ", np.allclose(np.transpose(np.conj(Psi))@Psi, Psi@np.conj(np.transpose(Psi)), rtol=1e-8))

Phi = np.conj(np.transpose(PhiH))
Phi = Finv @ Phi
Psi = Finv @ Psi

e0 = np.abs(Psi[0:n, 0])
e1 = np.abs(Psi[0:n, 1])
e2 = np.abs(Psi[0:n, 2])
e3 = np.abs(Psi[0:n, 3])
e4 = np.abs(Psi[0:n, 4])
e5 = np.abs(Psi[0:n, 5])

f0 = np.abs(Phi[0:n, 0])
f1 = np.abs(Phi[0:n, 2])
f2 = np.abs(Phi[0:n, 4])
f3 = np.abs(Phi[0:n, 6])
f4 = np.abs(Phi[0:n, 8])
f5 = np.abs(Phi[0:n, 10])


"Plotting"
# counts = linspace(1, 20, 20)
# fig = plt.figure()
# plt.semilogy(counts, Sigma[0:20],'.')
# plt.title('Singular values '+str(kx)+'-'+str(kz)+'-'+str(omega))
# plt.xlabel("index")
# plt.ylabel("sigma")
# fig.savefig('images/'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_singular_values.png')
# plt.show()

# fig, axs = plt.subplots(6)
# fig.suptitle('Singular response modes')
# axs[0].plot(y, np.real(e0))
# axs[1].plot(y, np.real(e1))
# axs[2].plot(y, np.real(e2))
# axs[3].plot(y, np.real(e3))
# axs[4].plot(y, np.real(e4))
# axs[5].plot(y, np.real(e5))
# plt.show()

# fig, axs = plt.subplots(6)
# fig.suptitle('Singular forcing modes')
# axs[0].plot(y, f0)
# axs[1].plot(y, f1)
# axs[2].plot(y, f2)
# axs[3].plot(y, f3)
# axs[4].plot(y, f4)
# axs[5].plot(y, f5)
# plt.show()

gridpts = 50
x = np.linspace(0, 1, gridpts)
y = np.linspace(-1, 0, gridpts)

X, Y = np.meshgrid(x, y)
Z = np.zeros([gridpts, gridpts])

for i in range(0, gridpts):
    for j in range(0, gridpts):
        Z[j, i] = resolventNorm(x[i], y[j], M, L, F)
    print("progress = ",(i*100)/gridpts,"%")
        
"Save data to file"

np.savetxt('psuedospectrum.dat', Z, delimiter=',')  
np.savetxt('X.dat', X, delimiter=',')
np.savetxt('Y.dat', Y, delimiter=',')

