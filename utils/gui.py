
import gmsh
import json

def set_gui():
    # This function set the configuration for the gmsh GUI.

    gmsh.option.setNumber("Mesh.MeshSizeMin", 
                          gmsh.option.getNumber("Mesh.MeshSizeMax"))

    gmsh.option.setNumber("General.ShowModuleMenu", 1)
    gmsh.option.setNumber("General.MenuWidth", 340)
    gmsh.option.setNumber("General.ColorScheme", 0)
    gmsh.option.setNumber("General.FltkColorScheme", 0)
    gmsh.option.setNumber('Geometry.Tolerance', 1e-4)
    gmsh.option.setNumber('Geometry.Points', 1)
    gmsh.option.setNumber('Geometry.Lines', 1)
    gmsh.option.setNumber('Geometry.Surfaces', 1)
    gmsh.option.setNumber('Geometry.Volumes', 1)
    gmsh.option.setNumber('Mesh.Points', 0)
    gmsh.option.setNumber('Mesh.PointType',1)
    gmsh.option.setNumber('Mesh.SurfaceEdges', 1)
    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.VolumeEdges', 1)
    gmsh.option.setNumber('Mesh.VolumeFaces', 1)
    gmsh.option.setNumber('Mesh.ColorCarousel', 0)
    gmsh.option.setColor('Mesh.Points', 0,0,210,255)
    gmsh.option.setColor('Mesh.Triangles', 210,210,0,255)
    gmsh.option.setColor('Mesh.Quadrangles', 210,210,0,255)
    gmsh.option.setColor('Mesh.Tetrahedra', 210,210,0,255)
    gmsh.option.setColor('Mesh.Hexahedra', 210,210,0,255)
    gmsh.option.setColor('Mesh.Prisms', 210,210,0,255)
    gmsh.option.setColor('Mesh.Pyramids', 210,210,0,255)
    gmsh.option.setColor('Mesh.Trihedra', 210,210,0,255)
    gmsh.view.add("Displacement",0)
    gmsh.option.setNumber('View[0].Visible', 0)
    gmsh.view.add("Force",1)
    gmsh.option.setNumber('View[1].Visible', 0)
    gmsh.view.add("Stresses",2)
    gmsh.option.setNumber('View[2].Visible', 0)

    param_gui = """
    [
      {
        "type":"number",
        "name":"0Modules/0Selection Mode",
        "values":[0],
        "choices":[0, 1, 2, 3],
        "valueLabels":{"Point": 0, "Curve": 1, "Surface": 2, "Body": 3}
      },
      {
        "type":"number",
        "name":"0Modules/1Physical ID",
        "visible":false
      },
      {
        "type":"string",
        "name":"ONELAB/Button",
        "values":["Clear Data", "clear"],
        "visible":false
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/0X/0X Value",
        "values":[0],
        "min":-1e9,
        "max":1e9
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/0X/1Free X",
        "values":[0],
        "choices":[0, 1]
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/1Y/0Y Value",
        "values":[0],
        "min":-1e9,
        "max":1e9
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/1Y/1Free Y",
        "values":[0],
        "choices":[0, 1]
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/2Z/0Z Value",
        "values":[0],
        "min":-1e9,
        "max":1e9
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/2Z/1Free Z",
        "values":[0],
        "choices":[0, 1]
      },
      {
        "type":"string",
        "name":"0Modules/Boundary Conditions/3Check Boundaries/0Type",
        "values":[" "],
        "readOnly":true,
        "attributes":{"Highlight":"AliceBlue"}
      },
      {
        "type":"string",
        "name":"0Modules/Boundary Conditions/3Check Boundaries/1X Value",
        "values":[" "],
        "readOnly":true,
        "attributes":{"Highlight":"AliceBlue"}
      },
      {
        "type":"string",
        "name":"0Modules/Boundary Conditions/3Check Boundaries/2Y Value",
        "values":[" "],
        "readOnly":true,
        "attributes":{"Highlight":"AliceBlue"}
      },
      {
        "type":"string",
        "name":"0Modules/Boundary Conditions/3Check Boundaries/3Z Value",
        "values":[" "],
        "readOnly":true,
        "attributes":{"Highlight":"AliceBlue"}
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/3Check Boundaries/4Physical ID",
        "values":[0],
        "choices":[0],
        "valueLabels":{"Empty": 0}
      },
      {
        "type":"string",
        "name":"0Modules/Boundary Conditions/3Check Boundaries/5Select Entity",
        "values":["check_boundary_conditions"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      },
      {
        "type":"number",
        "name":"0Modules/Boundary Conditions/4Type",
        "values":[0],
        "choices":[0, 1],
        "valueLabels":{"Displacement": 0, "Force": 1}
      },
      {
        "type":"string",
        "name":"0Modules/Boundary Conditions/6Set Boundary Conditions",
        "values":["set_boundary_conditions"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      },
      {
        "type":"string",
        "name":"0Modules/Boundary Conditions/7Remove Boundary",
        "values":["remove_boundary_conditions"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      },
      {
        "type":"string",
        "name":"0Modules/Solver/0Check Properties/0Cross-sectional Area",
        "values":[" "],
        "readOnly":true,
        "attributes":{"Highlight":"AliceBlue"}
      },
      {
        "type":"string",
        "name":"0Modules/Solver/0Check Properties/1Young's Modulus",
        "values":[" "],
        "readOnly":true,
        "attributes":{"Highlight":"AliceBlue"}
      },
      {
        "type":"number",
        "name":"0Modules/Solver/0Check Properties/2Physical ID",
        "values":[0],
        "choices":[0],
        "valueLabels":{"Empty": 0}
      },
      {
        "type":"string",
        "name":"0Modules/Solver/0Check Properties/3Select Entity",
        "values":["check_material_properties_truss"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      },
      {
        "type":"number",
        "name":"0Modules/Solver/1Cross-sectional Area",
        "values":[0],
        "min":-1e9,
        "max":1e9
      },
      {
        "type":"number",
        "name":"0Modules/Solver/2Young's Modulus",
        "values":[0],
        "min":-1e9,
        "max":1e9
      },
      {
        "type":"string",
        "name":"0Modules/Solver/3Set Material Properties",
        "values":["set_material_properties_truss"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      },
      {
        "type":"string",
        "name":"0Modules/Solver/4Remove Properties",
        "values":["remove_material_properties_truss"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      },
      {
        "type":"string",
        "name":"0Modules/Solver/5Run Solver",
        "values":["solver_truss"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      },
      {
        "type":"number",
        "name":"0Modules/Post-processing/0Steps Data/0Number of Steps",
        "values":[10],
        "min":0,
        "max":1e9,
        "step":1
      },
      {
        "type":"string",
        "name":"0Modules/Post-processing/0Steps Data/1Update Steps",
        "values":["update_steps"],
        "attributes":{"Macro":"Action", "Aspect":"Button"}
      }
    ]"""

    gmsh.onelab.set(param_gui)

def x_val():
    gmsh.onelab.set("""{
        "type":"number",
        "name":"0Modules/Boundary Conditions/0X/0X Value",
        "values":[0],
        "min":-1e9,
        "max":1e9,
        "visible":true
        }""")

def x_free():
    gmsh.onelab.set("""{
        "type":"number",
        "name":"0Modules/Boundary Conditions/0X/0X Value",
        "values":[0],
        "min":-1e9,
        "max":1e9,
        "visible":false
        }""")

def y_val():
    gmsh.onelab.set("""{
        "type":"number",
        "name":"0Modules/Boundary Conditions/1Y/0Y Value",
        "values":[0],
        "min":-1e9,
        "max":1e9,
        "visible":true
        }""")

def y_free():
    gmsh.onelab.set("""{
        "type":"number",
        "name":"0Modules/Boundary Conditions/1Y/0Y Value",
        "values":[0],
        "min":-1e9,
        "max":1e9,
        "visible":false
        }""")

def z_val():
    gmsh.onelab.set("""{
        "type":"number",
        "name":"0Modules/Boundary Conditions/2Z/0Z Value",
        "values":[0],
        "min":-1e9,
        "max":1e9,
        "visible":true
        }""")

def z_free():
    gmsh.onelab.set("""{
        "type":"number",
        "name":"0Modules/Boundary Conditions/2Z/0Z Value",
        "values":[0],
        "min":-1e9,
        "max":1e9,
        "visible":false
        }""")

def update_domains(a,b):

    vLabels = dict(zip(a,b))
    param = json.loads(gmsh.onelab.get("0Modules/1Physical ID"))
    param["values"] = [0]
    param["choices"] = b
    param["valueLabels"] = vLabels
    param["visible"] = (True if a else False)
    gmsh.onelab.set(json.dumps(param))

def update_check(a,b,c):
    
    vLabels = dict(zip(a,b))
    param = json.loads(gmsh.onelab.get(c))
    param["values"] = [0]
    param["choices"] = b
    param["valueLabels"] = vLabels
    gmsh.onelab.set(json.dumps(param))
