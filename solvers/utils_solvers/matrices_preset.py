
import numpy as np

class Matrices:
  def __init__(self):
    self.fe = []
    self.ke = []
    self.ue = []
    self.fn = []
    self.stresses = []
    self.F = []
    self.K = []
    self.U = []
    self.bF = []
    self.bK = []
    self.bU = []
    self.rol = []
    self.col = []
    self.idx_reduction = []
    self.idx_rebuild = []

def matrices_preset(a,b,c,d):

    # a: total number of nodes
    # b: total number of elements
    # c: number of DOFs by nodes
    # d: number of nodes by each element

    matrices = Matrices()
    cd = c*d
    matrices.ke = np.zeros([b,cd,cd])
    matrices.stresses = np.zeros([b,1])
    matrices.fe = np.zeros([b,1])
    return matrices
