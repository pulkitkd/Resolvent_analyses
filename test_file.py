import numpy as np

A = np.matrix([[1 ,2, 0],[0, 4, 3],[1, 2, 3]])
b = np.matrix([[8,2,3],[18,1,8],[9,8,9]])
  
sol = np.linalg.lstsq(A, A)[0]

print("sol = ", sol)
#%%
import numpy as np
a = np.array([1,2,3])
print(a@a)
# %%
# A = np.random.rand(3,3)
import numpy as np
A = [[10000000, 1],[0, 0.1]]
psi, sigma, phih = np.linalg.svd(A)
print("Psi = ",psi[1]@psi[1])
# %%
