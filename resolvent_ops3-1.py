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
Ny = 256
kx = 1.
kz = 1.
omega = 1

Re = 13500.35
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

"The shear velocity u_tau"
ustar = 0.045116673574063015
Rec = 3000.
# Rec = 2367.56
 
"Velocity field and its derivatives"
U96 = np.loadtxt('Umean.asc', skiprows=1)
print("U_mean_max = ",np.max(U96))
print("Re_tau = ", Rec*ustar/1.)
U = np.polyfit(y96, U96, 30)
U = np.poly1d(U)
# U = U(y)
U = 1 - y**2

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

H = np.linalg.solve(LHS, M)
F = get_F(D, wdiag, ksq, Ny)
Finv = inv(F)

Hs = F @ H @ Finv

Psi, Sigma, PhiH = np.linalg.svd(Hs, full_matrices=False)

Psi = Finv @ Psi
Phi = np.conj(np.transpose(PhiH))
Phi = Finv @ Phi
print("Sigma = ", Sigma)
print("np.dot(Sigma,Sigma) =", np.dot(Sigma,Sigma))
print("Sigma @ Sigma =", Sigma @ Sigma)


"""
Check SVD and othogonality
"""
print("check ", np.allclose(Hs, np.dot(Psi, np.dot(np.diag(Sigma), PhiH)), rtol=1e-10))

psi1 = Psi[:, 1]
psi2 = Psi[:, 7]

print("check ", np.allclose(Id2n, np.transpose(np.conj(Psi)) @ Psi, rtol=1e-8))
print("check ", np.allclose(np.transpose(np.conj(Psi))@Psi, Psi@np.conj(np.transpose(Psi)), rtol=1e-8))
print("Check perpendicular ", norm(psi1 @ np.conj(psi2)))

v0 = Psi[0:n, 0]
v1 = Psi[0:n, 1] 
v2 = Psi[0:n, 2]

e0 = v0*np.conj(v0)
e0 = e0/max(e0)
e1 = v1*np.conj(v1)
e1 = e1/max(e1)
e2 = v2*np.conj(v2)
e2 = e2/max(e2)


# print("Product of vectors = ",np.linalg.norm(np.dot(v0,v2)))
# print("Norm of the vector = ",np.linalg.norm(np.dot(v2,v2)))


"Plotting"
counts = linspace(1, 2*n, 2*n)
fig = plt.figure()
plt.plot(counts, Sigma,'.')
plt.title('Singular values '+str(kx)+'-'+str(kz)+'-'+str(omega))
plt.xlabel("index")
plt.ylabel("sigma")
fig.savefig('images/'+str(kx)+'-'+str(kz)+'-'+str(omega)+'_singular_values.png')
plt.show()

plt.show()
plt.plot(y, np.real(e0), "r", label = "first")
plt.plot(y, np.real(e1), "--g", label = "second")
plt.plot(y, np.real(e2), ":b", label = "third")
plt.xlabel("y")
plt.ylabel("Psi (response modes)")
plt.legend()
plt.show()

