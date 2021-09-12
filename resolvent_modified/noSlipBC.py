from cheb_trefethen import  *

def noSlipBC(D, M, L, n, a=50):
    
    # v^ = 0
    L[0, :] = 0.
    L[0, 0] = 1.
    L[n, :] = 0.
    L[n, n] = 1.

    M[0, :] = 0.
    M[0, 0] = a
    M[n, :] = 0.
    M[n, n] = a

    # eta^ = 0
    L[n+1, :] = 0.
    L[n+1, n+1] = 1.
    L[2*n+1, :] = 0.
    L[2*n+1, 2*n+1] = 1.

    M[n+1, :] = 0.
    M[n+1, n+1] = a
    M[2*n+1, :] = 0.
    M[2*n+1, 2*n+1] = a

    # dv^ / dy = 0
    L[1, :] = 0.
    L[n-1, :] = 0.
    L[1, 0:n+1] = a*D[0, :]
    L[n-1, 0:n+1] = a*D[n, :]

    M[1, :] = 0.
    M[n-1, :] = 0.
    M[1, 0:n+1] = D[0, :]
    M[n-2, 0:n+1] = D[n, :]
    
    return M, L
