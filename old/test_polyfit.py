import numpy as np
from cheb_trefethen import *
import matplotlib.pyplot as plt
"""
# Specify the grid on which interpolated U is desired
Ny = 256

# Load data(U) to be interpolated and the grid(y96) on which it is defined
U96 = np.loadtxt('Umean.asc', skiprows=1)
D, y96 = cheb(96)

# Specify the grid(y) on which interpolation is desired
D, y = cheb(Ny)

# polyfit will give the best fit polynomial coefficients of specified degree
U = np.polyfit(y96, U96, 4)

# poly1d will give a function to which we may input an x and get function value
U = np.poly1d(U)

# U_intp is the interpolated velocity
U_intp = U(y)

print("Norm of the difference = ", np.linalg.norm(U(y96)-U96))

# plot and check the 'eye' norm
plt.scatter(y96, U96)
plt.plot(y, U_intp,'r')
plt.xlabel("y")
plt.ylabel("U")
plt.show()
"""
nu = 0.00033333333333333332
Re = 135.3499
y = np.loadtxt('y.asc', skiprows=1)
yp = np.loadtxt('yp.asc', skiprows=1)
ustar = 0.045116673574063015
# print("Re_tau = yp/y = ", yp/y)
print("ustar/nu = ", ustar/nu)

