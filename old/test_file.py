import numpy as np
import matplotlib.pyplot as plt
from cheb_trefethen import *

Ny = 96
kx = 1.
kz = 1.
omega = 1

Re = 8445.6169437088125
Reinv = 1.0/Re
ksq = kx*kx + kz*kz
n = Ny+1

"Chebyshev grid and differentiation matrices"
D, y96 = cheb(96)
D, y = cheb(Ny)

U96 = np.loadtxt('turbstats_Re_30000.asc', skiprows=3)
U96uu = U96[:, 0]
U96vv = U96[:, 4]
U96ww = U96[:, 5]


plt.plot(y,U96uu)
# plt.plot(y,U96vv)
plt.plot(y,U96ww)
plt.xlabel("y")
plt.ylabel("U(y)")
# plt.legend()
plt.show()

#%%
import numpy as np
from numpy import diag
from cheb_trefethen import *
n = 4
D, y = cheb(n-1)
D2 = D@D
Idn = np.eye(n)
U = 1.0 - y*y
u = 1. - y*y
Uy = D@U
# print(y)
# print(D)
print(U[0:3])
# print(np.diag(U)*(Idn-D2))

# print(np.matmul(D, D) - D@D)
# %%
