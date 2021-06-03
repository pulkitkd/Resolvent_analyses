import numpy as np
# from scipy.linalg import eig, svd, pinv2
# from scipy.sparse.linalg import inv

import matplotlib.pyplot as plt
from numpy import pi, cos, sin
from numpy.linalg.linalg import norm, inv

from cheb_trefethen import *
from getResolventSVD import *

"""
Inputs
"""
# Re = 1005.5
# Re = 2389.558367
Re = 10000
Ny = 200
kx = 1.
kz = 0.
omega = 0.1
steps = 25
U96 = np.loadtxt('turbstats_Re_60000.asc', skiprows=3)
U96 = U96[:, 6]

singVals = zeros([steps,1])
variable = zeros(steps)
n = Ny + 1;
D, y96 = cheb(96)
D, y = cheb(Ny)
# U = np.polyfit(y96, U96, 30)
# U = np.poly1d(U)
# U = U(y)
U = 1 - y**2
Umax = np.max(U)

# stepsize = (1.1*kx*Umax - omega)/steps
omega2 = 1.5
omega1 = 0.
stepsize = (omega2-omega1)/steps
omega = omega1

print("steps = ", steps)
print("stepsize = ", stepsize)
print("Umax = ", Umax)


# for i in range(0, steps):
#    c = omega/kx
#    # print("omega = ", omega,"\n")
#    print("c/Umax = ", c/Umax, "\n")
#    variable[i] = c/Umax
#    Psi, Sigma, PhiH = getResolventSVD(Re, kx, kz, omega, Ny, U)
#    singVals[i,0] = Sigma[0];
#    singVals[i,1] = Sigma[2];
#    singVals[i,2] = Sigma[4];
#    singVals[i,3] = Sigma[6];
#    singVals[i,4] = Sigma[8];
#    omega = omega + stepsize;

for i in range(0, steps):
   # c = omega/kx
   # print("omega = ", omega,"\n")
   print("omega = ", omega, "\n")
   variable[i] = omega
   Psi, Sigma, PhiH = getResolventSVD(Re, kx, kz, omega, Ny, U)
   singVals[i,0] = Sigma[0];
   # singVals[i,1] = Sigma[2];
   # singVals[i,2] = Sigma[4];
   # singVals[i,3] = Sigma[6];
   # singVals[i,4] = Sigma[8];
   omega = omega + stepsize;


"Plotting"

# plt.loglog(variable, singVals)
plt.plot(variable,singVals)
plt.xlabel("c/Umax")
plt.ylabel("sigma")
# plt.ylim(0, 100)
plt.title('Singular values kx='+str(kx)+' kz='+str(kz))
plt.show()


