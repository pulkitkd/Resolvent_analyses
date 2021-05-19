import numpy as np
import matplotlib.pyplot as plt

kx = np.loadtxt('X.dat')
kz = np.loadtxt('Z.dat')
[X, Z] = np.meshgrid(kx, kz)
Y = np.loadtxt('energy.dat', delimiter=',')
# Ymax = np.max(Y)
# Y = Y / Ymax
# print("New max Y ",np.max(Y))

fig, ax = plt.subplots(1, 1)
cf = ax.contourf(X, Z, Y, interpolation='spline')
fig.colorbar(cf, ax=ax)
plt.yscale('log')
plt.xscale('log')
# ax.colorbar()
ax.set_title('Contour Plot')
# plt.contourf(Y, cmap=plt.cm.RdBu, vmin=abs(E).min(), vmax=abs(E).max())
plt.show()