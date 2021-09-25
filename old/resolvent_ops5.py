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
# Re = 1005.5 #Re_tau
Re = 2389.558367 #Re_tau
Ny = 128

lambdax1p = 10.
lambdax2p = 1.0e+5
lambdaz1p = 10.
lambdaz2p = 1.0e+5

yp = 100.
y1 = yp/Re
U96 = np.loadtxt('turbstats_Re_60000.asc', skiprows=3)
U96 = U96[:, 6]
D, y96 = cheb(96)
D, y = cheb(Ny)
U = np.polyfit(y96, U96, 30)
U = np.poly1d(U)
c = U(y1)
U = U(y)

steps = 30

E = zeros([steps, steps])
X = zeros(steps)
Z = zeros(steps)
n = Ny + 1;
Umax = np.max(U)

dlx = (lambdax2p - lambdax1p) / steps
dlz = (lambdaz2p - lambdaz1p) / steps

print("dlx, dlz = ", dlx," and ", dlz)
print("y1, yp = ", y1," and ", yp)
print("c = ", c)
print("c/Umax = ", c/Umax)

lambdaxp = lambdax1p
lambdazp = lambdaz1p

for i in range(0, steps):
   for j in range(0, steps):
      kx = 2*pi*Re/lambdaxp
      kz = 2*pi*Re/lambdazp
      omega = 2*pi*c*Re/lambdaxp
      Psi, Sigma, PhiH = getResolventSVD(Re, kx, kz, omega, Ny, U)
      sigma1 = Sigma[0]**2 + Sigma[1]**2
      totalE = Sigma @ Sigma
      E[i][j] = sigma1/totalE
      X[j] = lambdaxp;
      Z[i] = lambdazp;
      lambdaxp = lambdaxp + dlx
   lambdaxp = lambdax1p
   lambdazp = lambdazp + dlz
   print("lambdazp = ", lambdazp)

"Save data to file"

np.savetxt('energy.dat', E, delimiter=',')
np.savetxt('X.dat', X, delimiter=',')   
np.savetxt('Z.dat', Z, delimiter=',')   

