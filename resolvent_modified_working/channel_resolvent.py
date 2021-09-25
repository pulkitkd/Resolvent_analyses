import numpy as np
from numpy import pi, cos, sin, shape
from numpy.linalg.linalg import norm, inv, cond, solve, svd
import matplotlib.pyplot as plt
from scipy.linalg import sqrtm

import dmsuite as dm
from get_F import *
from OSQoperators import *

np.set_printoptions(suppress=True)
np.set_printoptions(precision=5)

n = 200
kx = 1.
kz = 6.
omega = 10.
Re = 2003
# Re = 1005.5
# Re = 2389.56

y, D4 = dm.cheb4c(n)           # fourth derivative with clamped BC (n-1 X n-1)
y, D = dm.chebdif(n, 2)         # first two derivatives (2 X n+1 X n+1)
D1 = D[0, :, :]                   # first derivative (n+1 X n+1)
D2 = D[1, :, :]                   # second derivative (n+1 X n+1)



"Velocity field and its derivatives============================================"
# y96, D96 = dm.chebdif(200, 1)    
# U96 = np.loadtxt('turbstats_Re_30000.asc', skiprows=3)
# U96 = np.loadtxt('velocity.txt')
# print(shape(U96), shape(y96))

# U = np.polyfit(y96, U96, 30)
# U = np.poly1d(U)
# U = U(y)
U = np.loadtxt('velocity.txt')
Uy = D1 @ U
Uyy = D2 @ U


"Apply homogeneous BCs========================================================="
y = y[1:n]
D1 = D1[1:n, 1:n]
D2 = D2[1:n, 1:n]
U = U[1:n]
Uy = Uy[1:n]
Uyy = Uyy[1:n]

"Create operators=============================================================="
M, L, B, C = OSQoperators(n, Re, kx, kz, U, Uy, Uyy, D1, D2, D4)

"Build weight matrix to enforce energy norm of the resolvent==================="
ksq = kx*kx + kz*kz
W, F = get_F(ksq, n)
F = sqrtm(W)
Finv = inv(F)

print("W = ", W)
print("diag(W) = ",diag(W))
print("shape(W) =", shape(W))

"Build a matrix to impose BCs on the RHS======================================="
Id2n = np.eye(2*n-2)


"Get the resolvent, scale it and take its SVD=================================="
H = inv(-1j*omega*Id2n - L)

print("H = ",H)
print("shape(H) = ", shape(H))

with open('H.txt','wb') as f:
    for line in H:
        np.savetxt(f, line, fmt='%.2f')

Hs = F @ H @ Finv
Psi, Sigma, PhiH = svd(Hs)


print("Singular values", Sigma[0:20])
"""
Check SVD and othogonality
"""
print("check ", np.allclose(Hs, Psi@diag(Sigma)@PhiH))
print("check ", np.allclose(Psi@np.transpose(np.conj(Psi)), Id2n))
print("check ", np.allclose(np.transpose(np.conj(Psi))@Psi, Psi@np.conj(np.transpose(Psi)), rtol=1e-8))

Psi = Finv @ Psi
Phi = np.transpose(np.conj(PhiH))
Phi = Finv @ Phi

Psi = C @ Psi
Phi = C @ Phi

"Plotting======================================================================"

# plt.plot(y,U)
# plt.title("Mean velocity profile")

Nmodes = 6

e = zeros([n-1, Nmodes])

for i in range(0, Nmodes):
    e[:, i] = Psi[2*n-2:3*n-3, i]
    print(i)
    

fig, axs = plt.subplots(Nmodes)
fig.suptitle('Singular response modes')
for i in range(0,Nmodes):
    axs[i].plot(y, np.real(e[:,i]))
    
fig.savefig('images/'+str(Re)+'-'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_response_modes.png')
plt.show()


counts = linspace(1, 20, 20)
fig = plt.figure()
plt.semilogy(counts, Sigma[0:20],'.')
plt.title('Singular values '+str(kx)+'-'+str(kz)+'-'+str(omega))
plt.xlabel("index")
plt.ylabel("sigma")
fig.savefig('images/'+str(Re)+'-'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_singular_values.png')
plt.show()
