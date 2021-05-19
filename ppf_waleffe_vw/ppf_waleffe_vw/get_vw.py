"""
Ref: On a self-sustaining process in shear flows - Waleffe, Phy.of Fluids, 1997
1. For a given gamma, get p by solving p tan(p) + gamma tanh (gamma) = 0
2. With gamma and p, construct vHat(y) ~ cos (py)/cos p - cosh (gamma y)/cosh gamma
3. find vHatmax
4. vHat = vHat/vHatmax so max(vHat) = 1
5. construct V = V0 vHat(y) cos(gamma z)
6. Obtain W from incompressibility dV/dy + dW/dz = 0, W periodic in z, W = 0 at y = 1, -1
7. Plot VW contour plots
8. Convert V and W into channelflow lingo

Author: Pratik Aghor
"""
from scipy import optimize
import numpy as np
from numpy import sin, cos, tan, arccos, cosh, tanh, arccosh, pi, arange
from numpy.fft import fft, ifft
import cheb
from cheb import *

import matplotlib
import matplotlib.pyplot as plt

#####################
"""
functions:
"""
def get_p(gamma_, guess):
    """
    solves ptan(p) + gamma_ tanh(gamma_) =  0
    this ensures dV/dy = 0 at y = +1, -1
    """
    f = lambda p : p*tan(p) + gamma_*tanh(gamma_)

    root = optimize.newton(f, guess)

    return root

#####################
"""
parameters
"""
gamma = 5.0/3.0
guess = 1.3
p = get_p(gamma, guess)
print ("p = ", p, "gamma = ", gamma)
R = 400 # Reynolds number
V0 = 1 # (p**2 + gamma**2)/R # amplitude, see Waleffe, 1997
#####################
"""
discretization and grid

y is [-1, 1] Chebyshev-Gauss-Lobatto grid (length = Ny+1)

z is the Physical domain [-pi/2, pi/2)
zc is the computational domain [0, 2\pi)

NOTE: We do not include the last point here, because periodic BC's are assumed.
Also, fft assumes periodic BC's and it would be erroneous to do the calculation
inclding the last point.

While saving and plotting, I pad up the solution as well as get the z_full, in
order to get the solution on the whole domain [-pi/2, pi/2]

xc = ax + b
"""
Ny = 16; Nz = 16;
D, y = cheb(Ny)

zc_min = 0;
zc_max = 2*pi;
a = 2.0;
b = pi;

dzc = (zc_max-zc_min)/(Nz-1)
zc = arange(zc_min, zc_max, dzc)
z = (zc - b*np.ones(Nz-1))/a

zc_full = arange(zc_min, zc_max+dzc, dzc)
z_full = (zc_full - b*np.ones(Nz))/a

kz = np.zeros(Nz);
kz[0:int(Nz/2)] = arange(0, int(Nz/2) ); kz[int(Nz/2+1):] = arange(int(-Nz/2+1),0,1);
print("kz = ", kz)

kzc = np.zeros(Nz-1);
for i in range(0, Nz-1):
	if i < int(Nz/2):
		kzc[i] = kz[i];
	else:
		kzc[i] = kz[i+1];
print("kzc = ", kzc)
#####################
"""
construct V
index j for y, k for z
"""
#####################
vHat = (cos(p*y)/cos(p)) - (cosh(gamma*y)/cosh(gamma))
idx = np.argmax(abs(vHat)) # find vHatmax
vHat = vHat/vHat[idx] # normalize

# print("vHat = \n", vHat)

fig = plt.figure(0)  # Create a figure instance
ax = fig.gca()
ax.plot(y, vHat)
ax.set_xlabel(r'$y$')  # Set x label
ax.set_ylabel(r'$\hat{V}$')  # Set y label
plt.savefig('vHat.png')
#####################
Z, Y = np.meshgrid(z_full, y)


V = np.zeros((Ny+1, Nz))
for j in range(0, Ny+1):
    for k in range(0, Nz):
        V[j, k] = vHat[j]*cos(gamma*z_full[k])
V = V0*V

# print("V = \n", V)
#####################
"""
construct W from continuity: dV/dy + dW/dz = 0,
W periodic in z, W = 0 at y = 1, -1

We take fft of V_y, divide by ik and take ifft to get W
fft(W_z) = i*kzc*fft(W)
dW/dz = -dV/dy
i*kzc*fft(W) = -fft_V_y
fft(W) = i*fft_V_y/kzc

W_periodic is used for [-\pi/2, pi/2)
W is actully W_full [-\pi/2, pi/2]

index j for y, k for z
"""
V_y = np.zeros((Ny+1, Nz))
for k in range(0, Nz):
    V_y[:, k] = np.matmul(D, V[:, k])
#####################
W_periodic = np.zeros((Ny+1, Nz-1))
fft_W = np.zeros(Nz-1, dtype = complex)
W = np.zeros((Ny+1, Nz))
for j in range(0, Ny+1):
    fft_V_y = fft(V_y[j, :])
    for k in range(1, Nz-1):
        fft_W[k] = 1j*fft_V_y[k]/kzc[k]
    W_periodic[j, :] = np.real(ifft(fft_W))

W[:, 0:Nz-1] = W_periodic[:, 0:Nz-1]
W[:, Nz-1] = W[:, 0] # periodic bc
#####################

#####################
"""
plotting
"""
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()
ax.quiver(Z, Y, W, V)
ax.set_xlabel(r'$z$')  # Set x label
ax.set_ylabel(r'$y$')  # Set y label
plt.savefig('vw.png')
#####################

fig = plt.figure(2)  # Create a figure instance
ax = fig.gca()
ax.contour(Z, Y, V)
ax.set_xlabel(r'$z$')  # Set x label
ax.set_ylabel(r'$y$')  # Set y label
plt.savefig('v_contours.png')
#####################
