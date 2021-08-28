
import gmsh
import sys
import time
import numpy as np

import solvers.truss_space1
import utils.tictoc as tictoc
import utils.gui as gui

class DomainTruss:
  def __init__(self):
    self.tag = []
    self.E = []
    self.A = []
    self.name = []

class BoundaryConditions:
  def __init__(self):
    self.dim = []
    self.tag = []
    self.type = []
    self.xyz = []
    self.name = []

class PostOptions:
  def __init__(self):
    self.type = []
    self.model = []
    self.check = False

class Change:
  def __init__(self):
    self.x = 0
    self.y = 0
    self.z = 0
    self.bound = 0
    self.truss = 0

class Check:
  def __init__(self):
    self.dim = -1
    self.bound_name = []
    self.bound_value = []
    self.truss_name = []
    self.truss_value = []

class PhysicalGroups:
  def __init__(self):
    self.name = []
    self.tag = []

def set_timer():
    # This function creates the TicTocs variables, enabling the measurement
    # of the times elapsed between the calling of tictoc.tic() and tictoc.toc().
    global TicToc, TicToc2
    TicToc = tictoc.TicTocGenerator()
    TicToc2 = tictoc.TicTocGenerator()

def checkForEvent():
    # check if an action is requested

    global domain_truss, boundary_conditions, post_options, change, check, physical_groups, \
        f_elements, d_nodal, s_elements, e_tags, e_ntags, n_el, n_nod, coord, ncoord, \
        n_steps, e_L, f_nodal

    while 1:
        root_path = "0Modules/0Selection Mode"
        dim = np.uint(gmsh.onelab.getNumber(root_path)[0])
        if gmsh.model.getCurrent() and dim != check.dim:
            aux = gmsh.model.getPhysicalGroups(dim)
            auxs = len(aux)
            check.dim = dim
            physical_groups.name = []
            physical_groups.tag = []
            for i in aux:
                physical_groups.name.append(gmsh.model.getPhysicalName(dim,i[1]))
                physical_groups.tag.append(i[1])
            gmsh.fltk.lock()
            gui.update_domains(physical_groups.name,np.arange(auxs).tolist())
            gmsh.fltk.unlock()
            gmsh.fltk.update()
            break
        else:
            break

    while 1:
        root_path = "0Modules/Boundary Conditions/"
        if gmsh.onelab.getNumber(root_path + "0X/1Free X")[0] == 0:
            if change.x == 1:
                change.x = 0
                gmsh.fltk.lock()
                gui.x_val()
                gmsh.fltk.unlock()
                gmsh.fltk.update()
            break
        elif gmsh.onelab.getNumber(root_path + "0X/1Free X")[0] == 1:
            if change.x == 0:
                change.x = 1
                gmsh.fltk.lock()
                gui.x_free()
                gmsh.fltk.unlock()
                gmsh.fltk.update()
            break

    while 1:
        root_path = "0Modules/Boundary Conditions/"
        if gmsh.onelab.getNumber(root_path + "1Y/1Free Y")[0] == 0:
            if change.y == 1:
                change.y = 0
                gmsh.fltk.lock()
                gui.y_val()
                gmsh.fltk.unlock()
                gmsh.fltk.update()
            break
        elif gmsh.onelab.getNumber(root_path + "1Y/1Free Y")[0] == 1:
            if change.y == 0:
                change.y = 1
                gmsh.fltk.lock()
                gui.y_free()
                gmsh.fltk.unlock()
                gmsh.fltk.update()
            break

    while 1:
        root_path = "0Modules/Boundary Conditions/"
        if gmsh.onelab.getNumber(root_path + "2Z/1Free Z")[0] == 0:
            if change.z == 1:
                change.z = 0
                gmsh.fltk.lock()
                gui.z_val()
                gmsh.fltk.unlock()
                gmsh.fltk.update()
            break
        elif gmsh.onelab.getNumber(root_path + "2Z/1Free Z")[0] == 1:
            if change.z == 0:
                change.z = 1
                gmsh.fltk.lock()
                gui.z_free()
                gmsh.fltk.unlock()
                gmsh.fltk.update()
            break

    while 1:
        root_path = "0Modules/Boundary Conditions/3Check Boundaries/"
        aux_check = np.uint(gmsh.onelab.getNumber(root_path + "4Physical ID")[0])
        if check.bound_name:
            if aux_check != change.bound:
                b_id = check.bound_name[aux_check]
                idx = boundary_conditions.name.index(b_id)
                x = str(boundary_conditions.xyz[idx][0])
                y = str(boundary_conditions.xyz[idx][1])
                z = str(boundary_conditions.xyz[idx][2])
                t = boundary_conditions.type[idx]
                gmsh.onelab.setString(root_path + "0Type",["Displacement"]) if t == 0 else gmsh.onelab.setString(root_path + "0Type",["Force"])
                gmsh.onelab.setString(root_path + "1X Value",[x])
                gmsh.onelab.setString(root_path + "2Y Value",[y])
                gmsh.onelab.setString(root_path + "3Z Value",[z])
                change.bound = aux_check
                gmsh.fltk.update()
                break
            else:
                break
        else:
            break

    while 1:
        root_path = "0Modules/Solver/0Check Properties/"
        aux_check = np.uint(gmsh.onelab.getNumber(root_path + "2Physical ID")[0])
        if check.truss_name:
            if aux_check != change.truss:
                b_id = check.truss_name[aux_check]
                idx = domain_truss.name.index(b_id)
                A = str(domain_truss.A[idx])
                E = str(domain_truss.E[idx])
                gmsh.onelab.setString(root_path + "0Cross-sectional Area",[A])
                gmsh.onelab.setString(root_path + "1Young's Modulus",[E])
                gmsh.onelab.setString(root_path + "2Physical ID",[b_id])
                change.truss = aux_check
                gmsh.fltk.update()
                break
            else:
                break
        else:
            break

    action = gmsh.onelab.getString("ONELAB/Action")

    if len(action) < 1:
        # no action requested
        pass

    elif action[0] == "check_boundary_conditions":
        gmsh.onelab.setString("ONELAB/Action", [""])
        root_path = "0Modules/Boundary Conditions/3Check Boundaries/4Physical ID"
        dim = np.uint(gmsh.onelab.getNumber("0Modules/0Selection Mode"))
        gmsh.fltk.setStatusMessage("Please select the Entity to check the boundaries (press 'q' to quit)", True)
        _, ent = gmsh.fltk.selectEntities(dim)
        check.bound_name = []
        check.bound_value = []
        if ent:
            aux = 0
            b_tags = gmsh.model.getPhysicalGroupsForEntity(dim,ent[0][1])
            for c,v in enumerate(b_tags):
                if gmsh.model.getPhysicalName(dim,v) in boundary_conditions.name:
                    check.bound_name.append(gmsh.model.getPhysicalName(dim,v))
                    check.bound_value.append(c-aux)
                else:
                    aux += 1
            if check.bound_name:
                gmsh.fltk.lock()
                gui.update_check(check.bound_name,check.bound_value,root_path)
                gmsh.fltk.unlock()
                gmsh.fltk.update()
                change.bound = 1
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("", True)
            else:
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("No Boundary Condition applied to the selected Entity yet", True)
        else:
            if gmsh.fltk.isAvailable() == 0: return 0
            gmsh.fltk.setStatusMessage("", True)

    elif action[0] == "set_boundary_conditions":
        gmsh.onelab.setString("ONELAB/Action", [""])
        root_path = "0Modules/Boundary Conditions/"
        dim = np.uint(gmsh.onelab.getNumber("0Modules/0Selection Mode"))
        n = np.nan
        auxx = gmsh.onelab.getNumber(root_path + "0X/1Free X")[0]
        auxy = gmsh.onelab.getNumber(root_path + "1Y/1Free Y")[0]
        auxz = gmsh.onelab.getNumber(root_path + "2Z/1Free Z")[0]
        if physical_groups.name:
            b_id = physical_groups.name[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])]
            if b_id in boundary_conditions.name:
                idx = boundary_conditions.name.index(b_id)
                boundary_conditions.tag.pop(idx)
                boundary_conditions.dim.pop(idx)
                boundary_conditions.type.pop(idx)
                boundary_conditions.xyz.pop(idx)
                boundary_conditions.name.pop(idx)
            boundary_conditions.tag.append(physical_groups.tag[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])])
            boundary_conditions.dim.append(dim)
            boundary_conditions.type.append(gmsh.onelab.getNumber(root_path + "4Type")[0])
            boundary_conditions.xyz.append([gmsh.onelab.getNumber(root_path + "0X/0X Value")[0] if auxx == 0 else n,
                                            gmsh.onelab.getNumber(root_path + "1Y/0Y Value")[0] if auxy == 0 else n,
                                            gmsh.onelab.getNumber(root_path + "2Z/0Z Value")[0] if auxz == 0 else n])
            boundary_conditions.name.append(b_id)
            if gmsh.fltk.isAvailable() == 0: return 0
            gmsh.fltk.setStatusMessage("Boundary Condition added at Physical Group: " + b_id, True)
        else:
            if gmsh.fltk.isAvailable() == 0: return 0
            gmsh.fltk.setStatusMessage("", True)
        
    elif action[0] == "remove_boundary_conditions":
        gmsh.onelab.setString("ONELAB/Action", [""])
        dim = np.uint(gmsh.onelab.getNumber("0Modules/0Selection Mode"))
        if physical_groups.name:
            b_id = physical_groups.name[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])]
            if b_id in boundary_conditions.name:
                idx = boundary_conditions.name.index(b_id)
                boundary_conditions.tag.pop(idx)
                boundary_conditions.dim.pop(idx)
                boundary_conditions.type.pop(idx)
                boundary_conditions.xyz.pop(idx)
                boundary_conditions.name.pop(idx)
                gmsh.model.removePhysicalGroups([(dim[0],physical_groups.tag[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])])])
                gmsh.fltk.update()
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("Boundary " + b_id + " removed", True)
            else:
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("No Boundary Condition applied to the Physical Group yet", True)
        else:
            if gmsh.fltk.isAvailable() == 0: return 0
            gmsh.fltk.setStatusMessage("", True)

    elif action[0] == "check_material_properties_truss":
        gmsh.onelab.setString("ONELAB/Action", [""])
        root_path = "0Modules/Solver/0Check Properties/2Physical ID"
        dim = np.uint(1)
        gmsh.fltk.setStatusMessage("Please select the Truss to check the properties (press 'q' to quit)", True)
        _, ent = gmsh.fltk.selectEntities(dim)
        check.truss_name = []
        check.truss_value = []
        if ent:
            aux = 0
            b_tags = gmsh.model.getPhysicalGroupsForEntity(dim,ent[0][1])
            for c,v in enumerate(b_tags):
                if gmsh.model.getPhysicalName(dim,v) in domain_truss.name:
                    check.truss_name.append(gmsh.model.getPhysicalName(dim,v))
                    check.truss_value.append(c-aux)
                else:
                    aux += 1
            if check.truss_name:
                gmsh.fltk.lock()
                gui.update_check(check.truss_name,check.truss_value,root_path)
                gmsh.fltk.unlock()
                gmsh.fltk.update()
                change.truss = 1
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("", True)
            else:
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("No Material Properties applied to the selected Truss yet", True)
        else:
            if gmsh.fltk.isAvailable() == 0: return 0
            gmsh.fltk.setStatusMessage("", True)

    elif action[0] == "set_material_properties_truss":
        gmsh.onelab.setString("ONELAB/Action", [""])
        root_path = "0Modules/Solver/"
        dim = np.uint(gmsh.onelab.getNumber("0Modules/0Selection Mode"))
        if dim == 1:
            if physical_groups.name:
                b_id = physical_groups.name[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])]
                if b_id in domain_truss.name:
                    idx =  domain_truss.name.index(b_id)
                    domain_truss.tag.pop(idx)
                    domain_truss.E.pop(idx)
                    domain_truss.A.pop(idx)
                    domain_truss.name.pop(idx)
                domain_truss.tag.append(physical_groups.tag[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])])
                domain_truss.E.append(gmsh.onelab.getNumber(root_path + "2Young's Modulus")[0])
                domain_truss.A.append(gmsh.onelab.getNumber(root_path + "1Cross-sectional Area")[0])
                domain_truss.name.append(b_id)
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("Material Property added at Physical Group: " + b_id, True)
            else:
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("", True)
        else:
            if gmsh.fltk.isAvailable() == 0: return 0
            gmsh.fltk.setStatusMessage("The Physical Group doesn't belong to a Truss", True)

    elif action[0] == "remove_material_properties_truss":
        gmsh.onelab.setString("ONELAB/Action", [""])
        dim = np.uint(gmsh.onelab.getNumber("0Modules/0Selection Mode"))
        if physical_groups.name:
            b_id = physical_groups.name[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])]
            if b_id in domain_truss.name:
                idx = domain_truss.name.index(b_id)
                domain_truss.tag.pop(idx)
                domain_truss.E.pop(idx)
                domain_truss.A.pop(idx)
                domain_truss.name.pop(idx)
                gmsh.model.removePhysicalGroups([(dim[0],physical_groups.tag[np.uint(gmsh.onelab.getNumber("0Modules/1Physical ID")[0])])])
                gmsh.fltk.update()
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("Material Property " + b_id + " removed", True)
            else:
                if gmsh.fltk.isAvailable() == 0: return 0
                gmsh.fltk.setStatusMessage("No Material Property applied to the Physical Group yet", True)
        else:
            if gmsh.fltk.isAvailable() == 0: return 0
            gmsh.fltk.setStatusMessage("", True)

    elif action[0] == "solver_truss":
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.fltk.setStatusMessage("Solving...", True)
        gmsh.fltk.lock()
        gmsh.model.mesh.generate(1)
        f_elements,d_nodal,f_nodal,s_elements,e_tags,e_ntags,n_el,n_nod,coord,ncoord,e_L = solvers.truss_space1.truss_solver(domain_truss,boundary_conditions,gmsh.model)
        post_options.type = "Truss"
        post_options.model = gmsh.model.getCurrent()
        post_options.check = True
        n_steps = gmsh.onelab.getNumber("0Modules/Post-processing/0Steps Data/0Number of Steps")[0]
        solvers.truss_space1.truss_post(f_elements,d_nodal,s_elements,e_tags,n_el,n_nod,post_options.model,n_steps,gmsh.view)
        solvers.truss_space1.plot_report(e_ntags,e_L,n_nod,d_nodal,f_nodal,s_elements,f_elements,post_options.model,gmsh.model,gmsh.view)
        gmsh.fltk.unlock()
        gmsh.fltk.update()
        if gmsh.fltk.isAvailable() == 0: return 0
        gmsh.fltk.setStatusMessage("Analysis finished", True)

    elif action[0] == "update_steps":
        gmsh.onelab.setString("ONELAB/Action", [""])
        if post_options.check:
            gmsh.fltk.lock()
            n_steps = gmsh.onelab.getNumber("0Modules/Post-processing/0Steps Data/0Number of Steps")[0]
            solvers.truss_space1.truss_post(f_elements,d_nodal,s_elements,e_tags,n_el,n_nod,post_options.model,n_steps,gmsh.view)
            gmsh.fltk.unlock()
            gmsh.fltk.update()
        if gmsh.fltk.isAvailable() == 0: return 0
        gmsh.fltk.setStatusMessage("Steps updated", True)

    elif action[0] == "clear":
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.clear()
        domain_truss = DomainTruss()
        boundary_conditions = BoundaryConditions()
        post_options = PostOptions()
        change = Change()
        check = Check()
        physical_groups = PhysicalGroups()
        if gmsh.fltk.isAvailable() == 0: return 0
        gmsh.fltk.setStatusMessage("Data cleared", True)

    return 1

set_timer()
tictoc.tic(TicToc)
gmsh.initialize()
domain_truss = DomainTruss()
boundary_conditions = BoundaryConditions()
post_options = PostOptions()
change = Change()
check = Check()
physical_groups = PhysicalGroups()
if len(sys.argv) > 1:
    gmsh.open(sys.argv[1])
gmsh.fltk.lock()
gui.set_gui()
gmsh.fltk.unlock()
gmsh.fltk.update()
tictoc.toc(TicToc)

if "-nopopup" not in sys.argv:
    gmsh.fltk.initialize()
    gmsh.fltk.closeTreeItem("0Modules/Boundary Conditions")
    gmsh.fltk.closeTreeItem("0Modules/Boundary Conditions/0X")
    gmsh.fltk.closeTreeItem("0Modules/Boundary Conditions/1Y")
    gmsh.fltk.closeTreeItem("0Modules/Boundary Conditions/2Z")
    gmsh.fltk.closeTreeItem("0Modules/Boundary Conditions/3Check Boundaries")
    gmsh.fltk.closeTreeItem("0Modules/Solver")
    gmsh.fltk.closeTreeItem("0Modules/Solver/0Check Properties")
    gmsh.fltk.closeTreeItem("0Modules/Post-processing/0Steps Data")
    gmsh.fltk.closeTreeItem("0Modules/Post-processing/1[2D] Plot Data")
    gmsh.fltk.update()
    while gmsh.fltk.isAvailable() and checkForEvent():
        gmsh.fltk.wait()

gmsh.finalize()
