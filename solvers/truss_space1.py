
import numpy as np
import solvers.utils_solvers.matrices_preset as m_preset
import solvers.utils_solvers.global_stiffness_assembly as gs_assembly
import solvers.utils_solvers.eqn_solver as eqn_s
import solvers.utils_solvers.set_boundary_conditions as set_bc

class ElementsTruss:
  def __init__(self):
    self.tag = []
    self.ntags = []
    self.L = []
    self.thetax = []
    self.thetay = []
    self.thetaz = []
    self.A = []
    self.EAL = []
    self.Cxyz = []

def truss_solver(d_truss,b_conditions,model):

    # the elemental stiffness matrix of the linear truss
    # will be 6x6, because each element of the truss has 2 nodes
    # and 3 dofs for each node, so (2*3)x(2*3)

    # the global stiffnes matrix of the linear truss 
    # will be (3*nodes)x(3*nodes), so it needs to be evaluated
    # for each analysis, because the number of the nodes
    # in the domain can change

    def truss_preset(d_truss,model):
        
        e_truss = ElementsTruss()
        dof = 3
        num_nodel = 2
        nodeTags, _, _ = model.mesh.getNodes(-1,-1)
        nodeTags = np.uint(nodeTags-1)
        num_nod = len(nodeTags)
        coord = np.zeros((num_nod,3))
        for i in nodeTags:
            coord[i] = model.mesh.getNode(np.uint(i+1))[0]
        p_ent = model.getEntities(0)
        nodeTags = nodeTags.tolist()
        for i in p_ent:
            idx = nodeTags.index(np.uint(model.mesh.getNodes(0,i[1])[0][0]-1))
            model.setEntityName(0,i[1],'#' + str(idx+1))
        e_truss.tag, e_truss.ntags = model.mesh.getElementsByType(1,-1)
        num_el = len(e_truss.tag)
        e_truss.tag = np.array(e_truss.tag)
        e_truss.ntags = np.uint(np.array(e_truss.ntags).reshape((num_el,2))-1)
        e_truss.A = np.zeros([num_el])
        e_truss.EAL = np.zeros([num_el])
        e_truss.Cxyz = np.zeros([num_el,6])
        m_truss = m_preset.matrices_preset(num_nod,num_el,dof,num_nodel)
        for c,v in enumerate(e_truss.ntags):
            dx = coord[v][1][0]-coord[v][0][0]
            dy = coord[v][1][1]-coord[v][0][1]
            dz = coord[v][1][2]-coord[v][0][2]
            e_truss.L.append(np.sqrt(dx**2+dy**2+dz**2))
            e_truss.thetax.append(np.arccos((dx/e_truss.L[c]))*(180/np.pi))
            e_truss.thetay.append(np.arccos((dy/e_truss.L[c]))*(180/np.pi))
            e_truss.thetaz.append(np.arccos((dz/e_truss.L[c]))*(180/np.pi))
        return num_nod, coord, e_truss, num_el, m_truss

    def truss_elemental_stiffness(d_truss,e_truss,m_truss,model):
        aux = e_truss.tag.tolist()
        for c,i in enumerate(d_truss.tag):
            entTags = model.getEntitiesForPhysicalGroup(1,i)
            EA = d_truss.E[c]*d_truss.A[c]
            for v in entTags:
                d_tags = model.mesh.getElementsByType(1,v)[0]
                for j in d_tags:
                    idx = aux.index(j)
                    x = e_truss.thetax[idx]*np.pi/180
                    w = e_truss.thetay[idx]*np.pi/180
                    v = e_truss.thetaz[idx]*np.pi/180
                    Cx = np.cos(x)
                    Cy = np.cos(w)
                    Cz = np.cos(v)
                    y = np.matrix([[Cx*Cx, Cx*Cy, Cx*Cz], [Cy*Cx, Cy*Cy, Cy*Cz], [Cz*Cx, Cz*Cy, Cz*Cz]])
                    e_truss.A[idx] = d_truss.A[c]
                    e_truss.EAL[idx] = EA/e_truss.L[idx]
                    e_truss.Cxyz[idx] = [np.negative(Cx), np.negative(Cy), np.negative(Cz), Cx, Cy, Cz]
                    m_truss.ke[idx] = e_truss.EAL[idx]*np.vstack((np.hstack((y,np.negative(y))),np.hstack((np.negative(y),y))))
        return e_truss, m_truss

    def truss_elemental_values(e_truss,m_truss,num_el,num_nod,a):

        #a: number of dof

        m_truss.ue = np.zeros((num_el,6,1))
        #m_truss.fn = np.zeros((num_el,6,1))
        for c,v in enumerate(e_truss.ntags):
            k = 0
            for i in v:
                for j in range(a):
                    m_truss.ue[c][np.int(a*k+j)] = m_truss.U[np.int(a*i+j)]
                    #m_truss.fn[c][np.int(a*k+j)] = m_truss.F[np.int(a*i+j)]
                k += 1
        d_nodal = np.zeros((1,num_nod,3))[0]
        f_nodal = np.zeros((1,num_nod,3))[0]
        ncoord = np.zeros((1,num_nod,3))[0]
        for i in range(num_nod):
            for j in range(a):
                auxd = m_truss.U[np.int(3*i+j)]
                auxf = m_truss.F[np.int(3*i+j)]
                d_nodal[i][j] = auxd
                f_nodal[i][j] = auxf
                ncoord[i][j] = coord[i][j]+auxd
        for i in range(num_el):
            m_truss.fe[i] = np.dot((e_truss.EAL[i]*e_truss.Cxyz[i]),m_truss.ue[i])[0]
            m_truss.stresses[i] = m_truss.fe[i]/e_truss.A[i]
        return m_truss.fe, m_truss.ue, m_truss.stresses, d_nodal, f_nodal, ncoord

    num_nod, coord, e_truss, num_el, m_truss = truss_preset(d_truss,model)
    e_truss, m_truss = truss_elemental_stiffness(d_truss,e_truss,m_truss,model)
    m_truss.K, m_truss.rol, m_truss.col = gs_assembly.global_stiffness_assembly(e_truss.ntags,m_truss.ke,num_el,3,2)
    m_truss.bF, m_truss.bU = set_bc.truss(e_truss,b_conditions,num_nod,3,model.mesh)
    m_truss.idx_rebuild = np.where(np.isnan(m_truss.bU))
    m_truss.F, m_truss.U = eqn_s.solver_static(m_truss)
    m_truss.fe, m_truss.ue, m_truss.stresses, d_nodal, f_nodal, ncoord = truss_elemental_values(e_truss,m_truss,num_el,num_nod,3)
    return m_truss.fe, d_nodal, f_nodal, m_truss.stresses, e_truss.tag, e_truss.ntags, num_el, num_nod, coord, ncoord, e_truss.L

def truss_post(f_elements,d_nodal,s_elements,e_tags,n_el,n_nod,m_Name,n_steps,views):
    
    a = d_nodal/n_steps
    views.addModelData(0,0,m_Name,"NodeData",np.arange(1,n_nod+1),a*0)
    for i in range(1,np.uint(n_steps+1)):
        views.addModelData(0,i,m_Name,"NodeData",np.arange(1,n_nod+1),a*(i))
    views.addModelData(1,0,m_Name,"ElementData",e_tags,f_elements.tolist())
    views.addModelData(2,0,m_Name,"ElementData",e_tags,s_elements.tolist())
    views.write(0, m_Name + "_displacementData.pos")
    views.write(1, m_Name + "_forceData.pos")
    views.write(2, m_Name + "_stressesData.pos")

def plot_report(e_ntags,e_L,n_nod,d_nodal,f_nodal,s_elements,f_elements,m_Name,model,views):

    import os
    import string
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import pdfkit
    from decimal import Decimal
    from PyPDF2 import PdfFileMerger

    phyg = model.getPhysicalGroups(1)
    p_ent = model.getEntities(0)
    ent = []
    for i in phyg:
        aux = []
        for j in model.getEntitiesForPhysicalGroup(i[0],i[1]):
            aux.append((1,j))
        ent.append(aux)
    alph = list(string.ascii_uppercase)
    nalph = alph.copy()
    n_ent = []
    pdfs = []
    aux = e_ntags.tolist()
    j = 0
    l = 0
    for vi in ent:
        k = 0
        for v in vi:
            fig = plt.figure(constrained_layout=True)
            gs = gridspec.GridSpec(2, 2, figure=fig)
            ax0 = fig.add_subplot(gs[0, :])
            ax1 = fig.add_subplot(gs[1, 0])
            ax2 = fig.add_subplot(gs[1, 1])
            n_ent.append(nalph[j]+str(k))
            model.setEntityName(1,v[1],n_ent[-1])
            k += 1
            nodeTags = np.uint(model.mesh.getNodes(1,v[1],includeBoundary = True)[0]-1).tolist()
            idx = aux.index(nodeTags)
            x2 = e_L[idx]
           # Displacement
            y1 = np.sqrt(d_nodal[nodeTags[0]][0]**2+d_nodal[nodeTags[0]][1]**2+d_nodal[nodeTags[0]][2]**2)
            y2 = np.sqrt(d_nodal[nodeTags[1]][0]**2+d_nodal[nodeTags[1]][1]**2+d_nodal[nodeTags[1]][2]**2)
            ax0.set_title('Displacement', fontsize=10)
            ax0.set_xlabel('Truss Length', fontsize=10)
            ax0.set_ylabel(' ', fontsize=10)
            ax0.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0), useOffset = True)
            ax0.text(0, 1.2, '2D Plot of the entity: ' + n_ent[-1],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax0.transAxes)
            ax0.text(1, 1.2, 'At node ' + '#' + str(nodeTags[0]+1) + ': ' + '{:.2e}'.format(Decimal(str(y1))),
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax0.transAxes)
            ax0.text(1, 1.15, 'At node ' + '#' + str(nodeTags[1]+1) + ': '+ '{:.2e}'.format(Decimal(str(y2))),
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax0.transAxes)
            ax0.plot([0,x2],[y1,y2])
            # Stress
            y1 = s_elements[idx][0]
            y2 = y1
            ax1.plot([0,x2],[y1,y2])
            ax1.set_title('Stress: ' + '{:.2e}'.format(Decimal(str(y1))), fontsize=10)
            ax1.set_ylabel(' ', fontsize=10)
            ax1.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0), useOffset = True)
            # Forces
            y1 = f_elements[idx][0]
            y2 = y1
            ax2.plot([0,x2],[y1,y2])
            ax2.set_title('Force: ' + '{:.2e}'.format(Decimal(str(y1))), fontsize=10)
            ax2.set_ylabel(' ', fontsize=10)
            ax2.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0), useOffset = True)
            pdf_file = n_ent[-1] + '_2d_plot.pdf'
            fig.set_size_inches(11.69,8.27)
            fig.tight_layout()
            plt.savefig(pdf_file, dpi=150)
            pdfs.append(pdf_file)
            #plt.show()
            ax0.cla()
            ax1.cla()
            ax2.cla()
        if j == 25: nalph = np.char.add(nalph,alph[l]); l += 1
        if l == 25: l = 0
        j += 1 if j <= 24 else -25
    auxn = np.uint(np.ceil(n_nod/40)) # number of pages
    auxd = np.resize(d_nodal,(40*auxn,1,3))
    plt.figure(1)
    plt.clf()
    for i in range(auxn):
        plt.suptitle('Displacement Node Data')
        plt.text(.1, .95, 'Tag of the Node | [X Y Z]'+'\n')
        plt.axis('off')
        k = np.uint(i*40)
        for c,v in enumerate(auxd[k:np.uint(k+40),0,:]):
            plt.text(.1, .9-(0.025*c), '#'+(str((c+1+k))+' | '+str(v)+'\n'))
            print(c+1+k)
            if c+1+k == n_nod:
                break
        pdf_file = str(i) + '_d_nodal.pdf'
        plt.figure(1).set_size_inches(11.69,8.27)
        plt.tight_layout()
        plt.savefig(pdf_file, dpi=150)
        plt.clf()
        pdfs.append(pdf_file)
    auxd = np.resize(f_nodal,(40*auxn,1,3))
    plt.figure(1)
    plt.clf()
    for i in range(auxn):
        plt.suptitle('Force Node Data')
        plt.text(.1, .95, 'Tag of the Node | [X Y Z]'+'\n')
        plt.axis('off')
        k = np.uint(i*40)
        for c,v in enumerate(auxd[k:np.uint(k+40),0,:]):
            plt.text(.1, .9-(0.025*c), '#'+(str((c+1+k))+' | '+str(v)+'\n'))
            print(c+1+k)
            if c+1+k == n_nod:
                break
        pdf_file = str(i) + '_f_nodal.pdf'
        plt.figure(1).set_size_inches(11.69,8.27)
        plt.tight_layout()
        plt.savefig(pdf_file, dpi=150)
        plt.clf()
        pdfs.append(pdf_file)
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    merger.write(m_Name + "_plot_report.pdf")
    merger.close()
    for i in pdfs:
        os.remove(i)