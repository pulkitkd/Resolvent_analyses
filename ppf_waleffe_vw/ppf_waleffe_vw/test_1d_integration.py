import numpy as np
from numpy import sin, cos, tan, arccos, cosh, tanh, arccosh, pi, arange
from numpy.fft import fft, ifft
#####################
"""
give dW/dz = -sin(z) and obtain W = cos(z) + c
"""
"""
discretization and grid

z is the Physical domain [-pi/2, pi/2)
zc is the computational domain [0, 2\pi)

NOTE: We do not include the last point here, because periodic BC's are assumed.
Also, fft assumes periodic BC's and it would be erroneous to do the calculation
inclding the last point.

While saving and plotting, I pad up the solution as well as get the z_full, in
order to get the solution on the whole domain [-pi/2, pi/2]

zc = az + b
"""

Nz = 16;

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
fft(W_z) = i*kzc*fft(W)
i*kzc*fft(W) = fft_W_1z
fft(W) = -i*fft_W_1z/kzc

W_periodic is used for [-\pi/2, pi/2)
W is actully W_full [-\pi/2, pi/2]
"""
W_1z = -sin(zc)

fft_W_1z = fft(W_1z)
fft_W = np.zeros(Nz-1, dtype = complex)
W_periodic = np.zeros(Nz-1)
W = np.zeros(Nz)

for k in range(1, Nz-1):
	fft_W[k] = -1j*fft_W_1z[k]/kzc[k]

W_periodic = np.real(ifft(fft_W))
W[0:Nz-1] = W_periodic[0:Nz-1]
W[Nz-1] = W[0]
#####################
"""
plotting
"""
import matplotlib
import matplotlib.pyplot as plt

fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()
ax.plot(z_full, W, '-o')
ax.plot(z_full, cos(zc_full))

ax.set_xlabel(r'$W$')  # Set x label
ax.set_ylabel(r'$z$')  # Set y label
plt.savefig('test_W.png')
#####################
