
import numpy as np

def truss(a,b,c,d,e):

    #a: object for the elements properties
    #b: object for the boundary conditions
    #c: number of nodes
    #d: dof
    #e: gmsh object

    F = np.zeros((c*d,1))
    U = np.zeros((c*d,1))
    U.fill(np.nan)
    for cc,v in enumerate(b.tag):
        b_tags = e.getNodesForPhysicalGroup(b.dim[cc],v)[0]-1
        if b.type[cc] == 0:
            for j in b_tags.tolist():
                U[3*j] = b.xyz[cc][0]
                U[3*j+1] = b.xyz[cc][1]
                U[3*j+2] = b.xyz[cc][2]
        elif b.type[cc] == 1:
            lb_tags = len(b_tags)
            xn = b.xyz[cc][0]/lb_tags
            yn = b.xyz[cc][1]/lb_tags
            zn = b.xyz[cc][2]/lb_tags
            for j in b_tags.tolist():
                F[3*j] = xn
                F[3*j+1] = yn
                F[3*j+2] = zn
    return F, U
