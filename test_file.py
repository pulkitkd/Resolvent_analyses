import numpy as np

A = np.matrix([[1 ,2, 0],[0, 4, 3],[1, 2, 3]])
b = np.matrix([[8,2,3],[18,1,8],[9,8,9]])
  
sol = np.linalg.lstsq(A, A)[0]

print("sol = ", sol)