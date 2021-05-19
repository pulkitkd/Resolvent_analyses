import numpy as np
# from scipy.linalg import eig, svd, pinv2
# from scipy.sparse.linalg import inv

import matplotlib.pyplot as plt
from numpy import pi, cos, sin
from numpy.linalg.linalg import norm, inv

from cheb_trefethen import *
from get_F import *


"""
Inputs
"""
Re = 135.35
Reinv = 1.0/Re
Ny = 96
kx = 1.
kz = 6.
omega = 0.5

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
U = np.loadtxt('Umean.asc', skiprows=1)
Uy = D @ U
Uyy = D @ Uy

y, w = clencurt(Ny)
wdiag = np.diag(w)

"""
M, L and B matrices
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

"""
Solve the Orr-Sommerfeld problem (a useful check on the BCs and operators)
   L x = i omega M x
=> Minv L x = i omega x
=> i omega = eigval(Minv*L) 
"""

# """
# Get the eigenvalues with largest imaginary parts
# """
# MinvL = np.linalg.lstsq(M, L)[0]
# w, v = eig(MinvL, check_finite=True)
# w = -1j*w
# w = np.sort_complex(w)

# """
# Plot the eigenvalues
# """
# plt.scatter(np.real(w), np.imag(w))
# plt.xlim((0,1))
# plt.ylim((-1,0))
# plt.xlabel("real")
# plt.ylabel("imag")
# plt.show()

"""
Get the Resolvent matrix defined as inverse of LHS
LHS = -i om M + L
Obtain H as a solution of LHS x = B x
!! Note that one does not yet have x on the RHS
Take the SVD of H
"""

H = np.linalg.solve(LHS, M)
F = get_F(D, wdiag, ksq, Ny)
Finv = inv(F)

Hs = F@H@Finv

Psi, Sigma, PhiH = np.linalg.svd(Hs, full_matrices=False)

"""
Check SVD and othogonality
"""
print("norm = ", norm(Hs - (Psi @ np.diag(Sigma) @ PhiH)))
print("check ", np.allclose(Hs, np.dot(Psi, np.dot(np.diag(Sigma), PhiH)),rtol=1e-10))

psi1 = Psi[:, 1]
psi2 = Psi[:, 7]

print("check ", np.allclose(Id2n, Psi @ np.transpose(Psi), rtol=1e-8))
print("norm ", norm(Id2n - Psi @ np.transpose(Psi)))

v0 = Psi[0:n, 0]
v1 = Psi[0:n, 1] 
v2 = Psi[0:n, 2]

e0 = v0*np.conj(v0)
e0 = e0/max(e0)
e1 = v1*np.conj(v1)
e1 = e1/max(e1)
e2 = v2*np.conj(v2)
e2 = e2/max(e2)


"Plotting"
counts = linspace(1, 2*n, 2*n)
fig = plt.figure()
plt.scatter(counts,Sigma)
plt.title('Singular values '+str(kx)+'-'+str(kz)+'-'+str(omega))
plt.xlabel("index")
plt.ylabel("sigma")
fig.savefig('images/'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_singular_values.png')
plt.show()

# fig1 = plt.figure()
# plt.plot(y, v0*np.conj(v0))
# plt.title('First velocity mode kinetic energy '+str(kx)+'-'+str(kz)+'-'+str(omega))
# plt.xlabel("y")
# plt.ylabel("Psi")
# fig1.savefig('images/'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_velocity'+str(0)+'.png')
# plt.show()

# fig2 = plt.figure()
# plt.plot(y, v1*np.conj(v1))
# plt.title('Second velocity mode kinetic energy'+str(kx)+'-'+str(kz)+'-'+str(omega))
# plt.xlabel("y")
# plt.ylabel("Psi")
# fig2.savefig('images/'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_velocity'+str(1)+'.png')
# plt.show()

# fig3 = plt.figure()
# plt.plot(y, v2*np.conj(v2))
# plt.title('Third velocity mode kinetic energy'+str(kx)+'-'+str(kz)+'-'+str(omega))
# plt.xlabel("y")
# plt.ylabel("Psi")
# fig3.savefig('images/'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_velocity'+str(2)+'.png')
# plt.show()

plt.show()
plt.plot(y, np.real(e0), "r", label = "first")
plt.plot(y, np.real(e1), "--g", label = "second")
plt.plot(y, np.real(e2), ":b", label = "third")
plt.xlabel("y")
plt.ylabel("Psi (response modes)")
plt.legend()
plt.show()

# plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')