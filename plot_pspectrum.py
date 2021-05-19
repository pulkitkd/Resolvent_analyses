import numpy as np

import matplotlib.pyplot as plt

Z = np.loadtxt('psuedospectrum.dat', delimiter=',')
Y = np.loadtxt('Y.dat', delimiter=',')
X = np.loadtxt('X.dat', delimiter=',')

w = np.loadtxt('OS_eigenvalues.dat', delimiter=',')

logLimit1 = -10
logLimit2 = 10
noOfContours = 2*(logLimit2 - logLimit1 + 1)
print(np.logspace(logLimit1, logLimit2, noOfContours))

plt.scatter(w[0,:], w[1,:])
plt.xlim(0,1)
plt.ylim(-1,0)
plt.contour(X, Y, Z, np.logspace(logLimit1, logLimit2, noOfContours), cmap='RdGy');
plt.show()

