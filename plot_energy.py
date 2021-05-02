import numpy as np
import matplotlib.pyplot as plt

E = np.loadtxt('energy.dat', delimiter=',')
print("E = ", E)
plt.imshow(E, cmap=plt.cm.RdBu, vmin=abs(E).min(), vmax=abs(E).max())
plt.xrange(1,1e+5)
plt.yscale('log')
plt.xscale('log')
plt.show()