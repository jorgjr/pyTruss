
import numpy as np
import scipy.sparse.linalg as spl

def solver_static(a):

    #a: object that contains the parameters for the system of equation

    f = a.bF[a.idx_rebuild[0]]
    aux = a.K[a.idx_rebuild[0],:]
    k = aux[:,a.idx_rebuild[0]]
    u = spl.minres(k,f)
    a.U = np.copy(a.bU)
    for c,v in enumerate(a.idx_rebuild[0]):
        a.U[v] = u[0][c]
    a.F = np.asmatrix(a.K.toarray())*np.asmatrix(a.U)
    return a.F, a.U
