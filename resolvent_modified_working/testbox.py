import numpy as np
from numpy import pi, cos, sin
from numpy.linalg.linalg import norm, inv, cond, solve
import matplotlib.pyplot as plt

import dmsuite as dm
from cheb_trefethen import *
from noSlipBC import *
from get_F import *


# Ny = 9
# kx = 1
# kz = 1
# ksq = kx**2 + kz**2
# D, y = cheb(Ny)
# y, w = clencurt(Ny)
# wdiag = np.diag(w)

# # show(wdiag)

# print(f"{D=} {np.shape(D)=}")

# W, F = get_F(D, wdiag, ksq, Ny)

# print(f"{F=} {np.shape(F)=}")

# import numpy as np
# import matplotlib.pyplot as plt
"""
ncheb = 32; mder = 2; pi = np.pi
x, D = dm.chebdif(ncheb, mder)        # first two derivatives
D1 = D[0,:,:]                         # first derivative
D2 = D[1,:,:]                         # second derivative
y = np.sin(2 * pi * x)                      # function at Chebyshev nodes
yd = 2 * pi * np.cos(2 * pi * x)        # theoretical first derivative
ydd = - 4 * pi ** 2 * np.sin(2 * pi * x)  # theoretical second derivative
fig, axe = plt.subplots(3, 1, sharex=True)
axe[0].plot(x, y)
axe[0].set_ylabel(r'$y$')
axe[1].plot(x, yd, '-')
axe[1].plot(x, np.dot(D1, y), 'o')
axe[1].set_ylabel(r'$y^{\prime}$')
axe[2].plot(x, ydd, '-')
axe[2].plot(x, np.dot(D2, y), 'o')
axe[2].set_xlabel(r'$x$')
axe[2].set_ylabel(r'$y^{\prime\prime}$')
plt.show()
"""
n = 200
kx = 1.
kz = 6.
omega = 0.5
Re = 1000


y, D4 =  dm.cheb4c(n)           # fourth derivative with clamped BC (n-1 X n-1)
y96, D96 = dm.chebdif(200, 1)    # first derivative 
y, D = dm.chebdif(n, 2)         # first two derivatives (2 X n+1 X n+1)
D1 = D[0,:,:]                   # first derivative (n+1 X n+1)
D2 = D[1,:,:]                   # second derivative (n+1 X n+1)

print("shape of D = ", shape(D1))
D1 = D1[1:n, 1:n]
print("shape of D = ", shape(D1))

U96 = np.loadtxt('velocity.txt')
U = np.polyfit(y96, U96, 30)
U = np.poly1d(U)
U = U(y)
print(shape(U),shape(U96))
print(norm(U - U96))
plt.plot(U,y,U96,y96)
plt.show()

# f = cos(pi*y)
# f4 = (1/pi**4)*D4@f
# fig, axs = plt.subplots(2)
# fig.suptitle('Differentiation test')
# axs[0].plot(y, f)
# axs[1].plot(y[50:n-50], f4[50:n-50])
# plt.show()


###############################################################################
"""
clencurt(Ny) 
- outputs a (2, Ny+1) sized arrays
- gives the clencurt grid points and their weights
"""