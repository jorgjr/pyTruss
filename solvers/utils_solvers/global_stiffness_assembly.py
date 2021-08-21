
import numpy as np
import scipy.sparse as sp

def global_stiffness_assembly(a,b,c,d,e):

    #a: node tags of the element
    #b: elemental stiffness matrix
    #c: total number of elements
    #d: dof
    #e: number of node by element

    ndof = d*e
    auxm = np.reshape(a,(-1,1)).astype(int)
    bb = len(auxm)
    row = np.zeros((bb,bb)).astype(int)
    col = np.zeros((bb,bb)).astype(int)
    for i in range(d):
        auxmm = (d*auxm)+i
        arow = np.repeat(auxmm,ndof,axis=1)
        row = np.hstack((row,arow))
        col = np.hstack((col,auxmm))
    row = np.reshape(row[:,bb::],(1,-1))
    col = np.tile(np.reshape(col[:,bb::],(c,ndof)),ndof)
    K = sp.csr_matrix((np.ravel(b), (np.ravel(row), np.ravel(col))))
    return K, row, col
