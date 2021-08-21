
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

def delete_from_csr(mat, row_indices=[], col_indices=[]):
    if not isinstance(mat, sp.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    rows = []
    cols = []
    if row_indices:
        rows = list(row_indices)
    if col_indices:
        cols = list(col_indices)
    if len(rows) > 0 and len(cols) > 0:
        row_mask = np.ones(mat.shape[0], dtype=bool)
        row_mask[rows] = False
        col_mask = np.ones(mat.shape[1], dtype=bool)
        col_mask[cols] = False
        return mat[row_mask][:,col_mask]
    elif len(rows) > 0:
        mask = np.ones(mat.shape[0], dtype=bool)
        mask[rows] = False
        return mat[mask]
    elif len(cols) > 0:
        mask = np.ones(mat.shape[1], dtype=bool)
        mask[cols] = False
        return mat[:,mask]
    else:
        return mat

def solver_static(a):

    #a: object that contains the parameters for the system of equation

    f = np.delete(a.bF,a.idx_reduction[0],0)
    k = delete_from_csr(a.K, a.idx_reduction[0].tolist(), a.idx_reduction[0].tolist())
    u = spl.minres(k,f)
    a.U = np.copy(a.bU)
    for c,v in enumerate(a.idx_rebuild[0]):
        a.U[v] = u[0][c]
    a.F = np.asmatrix(a.K.toarray())*np.asmatrix(a.U)
    return a.F, a.U
