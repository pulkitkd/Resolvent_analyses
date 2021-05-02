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
Re = 135.35
Ny = 99
kx = 1.
kz = 1.
omega = 0.05
steps = 50
U96 = np.loadtxt('Umean.asc', skiprows=1)


singVals = zeros([steps,5])
variable = zeros(steps)
n = Ny + 1;
D, y96 = cheb(96)
D, y = cheb(Ny)
U = np.polyfit(y96, U96, 30)
U = np.poly1d(U)
U = U(y)
Umax = np.max(U)

stepsize = (1.1*kx*Umax - omega)/steps

print("steps = ", steps)
print("stepsize = ", stepsize)
print("Umax = ", Umax)


for i in range(0, steps):
   c = omega/kx
   print("omega = ", omega,"\n")
   print("c/Umax = ", c/Umax, "\n")
   variable[i] = c/Umax
   Psi, Sigma, PhiH = getResolventSVD(Re, kx, kz, omega, Ny, U)
   singVals[i,0] = Sigma[0];
   singVals[i,1] = Sigma[1];
   singVals[i,2] = Sigma[2];
   singVals[i,3] = Sigma[3];
   singVals[i,4] = Sigma[4];
   omega = omega + stepsize;


"Plotting"

plt.loglog(variable, singVals)
plt.xlabel("c/Umax")
plt.ylabel("sigma")
plt.title('Singular values kx='+str(kx)+' kz='+str(kz))
plt.show()


