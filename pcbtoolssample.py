import os
import sys

import gmsh

sys.path.append(os.environ["HOME"]+"/pcb-tools")

from gerber import PCB
from gerber.render import theme
from gerber.render.cairo_backend import GerberCairoContext

path = "./sample1"

gmsh.initialize()
ctx = GerberCairoContext()
pcb = PCB.from_directory(path)
ctx.render_layers(pcb.top_layers, path + '/../pcb_top.png', theme.THEMES['OSH Park'], max_width=800, max_height=600)
surfaces_index = []
fused = ctx.gmsh_surf.pop(0)
#print(ctx.gmsh_surf)
#print(ctx.gmsh_hole)
"""
for surface in ctx.gmsh_surf:
	sf = gmsh.model.occ.fuse([(2,fused)], [(2, surface)], removeObject=True, removeTool=True)
	print(sf)
	gmsh.model.occ.synchronize()
	fused = sf[0][0][1]
for surface in ctx.gmsh_surf:
	surfaces_index.append((2, surface))
#surfaces_index = ([(2, surface)] for surface in ctx.gmsh_surf)
print(surfaces_index)
surfaces = gmsh.model.occ.fragment(surfaces_index, [], removeObject=True, removeTool=True)
#while len(ctx.gmsh_surf) > 1:
#	ctx.gmsh_surf.append(gmsh.model.occ.fuse([(2, ctx.gmsh_surf.pop(0))],[(2, ctx.gmsh_surf.pop(0))]))
holes_index = []
for hole in ctx.gmsh_hole:
	holes_index.append((2,hole))
holes = gmsh.model.occ.fragment(holes_index, [], removeObject=True, removeTool=True)
print("print fragments")
print(surfaces)
print(surfaces[0][-1][1])
print(holes)
print(holes[0][-1][1])
#holes = gmsh.model.occ.fuse([holes_index.pop(0)], holes_index, removeObject=True, removeTool=True)
#holes = gmsh.model.occ.fuse(([(2, surface)] for surface in ctx.gmsh_hole))
gmsh.model.occ.cut([(2,surfaces[0][-1][1])], [(2,holes[0][-1][1])], removeObject=True, removeTool=True)
"""

surfaces_index = []
for s_iter in ctx.gmsh_surf:
	surface = s_iter[0]
	holes = s_iter[1]
	for hole in holes:
		#print([[(2, surface)],[(2,hole)]])
		s = gmsh.model.occ.cut([(2, surface)],[(2,hole)], removeObject=True, removeTool=False)
		#print(s)
		surface = s[0][0][1]
		gmsh.model.occ.remove([(2,hole)], recursive=True)
	for hole in ctx.gmsh_hole:
		#print([[(2, surface)],[(2,hole[0])]])
		s = gmsh.model.occ.cut([(2, surface)],[(2,hole[0])], removeObject=True, removeTool=False)
		#print(s)
		surface = s[0][0][1]
	surfaces_index.append((2, surface))
for hole in ctx.gmsh_hole:
	gmsh.model.occ.remove([(2,hole[0])], recursive=True)
surfaces = gmsh.model.occ.fragment(surfaces_index, [], removeObject=True, removeTool=True)

gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(2,[x[1] for x in surfaces[0]],0)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin",ctx.lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax",ctx.lc)
gmsh.model.mesh.generate(2)
gmsh.fltk.run()
gmsh.write("t1_.msh")


from dolfinx.io import gmshio
from mpi4py import MPI

gmsh_model_rank = 0
gdim = 2
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)
from dolfinx import fem
V = fem.FunctionSpace(domain, ("CG", 1))
import ufl
from petsc4py.PETSc import ScalarType
x = ufl.SpatialCoordinate(domain)
"""
set boundary conditions
"""

import re
import os
filepath = ""
for file in os.listdir(path):
    if(re.match(".*(-PTH.drl)$",file)):
        filepath = file

drills = []
holes = []
radius = 0
print(path + "/" + filepath)

with open(path + "/" + filepath) as f:
    for s_line in f:
        #print(s_line)
        result = re.match("T\d+C\d[.]\d+",s_line)
        if(result):
            print(s_line)
            print(s_line.split("C"))
            drills.append(s_line.split("C"))
        for drill in drills:
            if(re.match(drill[0],s_line)):
                print(s_line)
                radius = drill[1].split("\n")[0]
        if(re.match("X\d+[.]\d+Y[-]\d+[.]\d+",s_line)):
            print(s_line.split("X")[1].split("\n")[0].split("Y"))
            holes.append([radius, s_line.split("X")[1].split("\n")[0].split("Y")])

print(holes)


import numpy as np

bcs = []

for hole in holes:
	def on_boundary(x):
		return np.isclose(np.sqrt((x[0] - float(hole[1][0]) * 25.4)**2 + (x[1] - float(hole[1][1]) * 25.4)**2), float(hole[0]) * 25.4 / 2)
	boundary_dofs = fem.locate_dofs_geometrical(V, on_boundary)
	bcs.append(fem.dirichletbc(ScalarType(1), boundary_dofs, V))

"""
def on_boundary1(x):
    return np.isclose(np.sqrt((x[0] - 17.78)**2 + (x[1] + 20.32)**2), 0.2)
	#return np.isclose(np.sqrt((x[0] - 0.7*25.4)**2 + (x[1] - 0.8*25.4)**2), 0.157 * 25.4)
def on_boundary2(x):
    return np.isclose(np.sqrt((x[0] - 27.94)**2 + (x[1] + 20.32)**2), 0.2)
    #return np.isclose(np.sqrt((x[0] - 1.1*25.4)**2 + (x[1] - 0.8*25.4)**2), 0.157 * 25.4)

boundary_dofs1 = fem.locate_dofs_geometrical(V, on_boundary1)
bc1 = fem.dirichletbc(ScalarType(1), boundary_dofs1, V)

boundary_dofs2 = fem.locate_dofs_geometrical(V, on_boundary2)
bc2 = fem.dirichletbc(ScalarType(0), boundary_dofs2, V)
bcs = [bc1,bc2]
"""

f = fem.Constant(domain, ScalarType(0))
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx
problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

import pyvista
import dolfinx


u_topology, u_cell_types, u_geometry = dolfinx.plot.create_vtk_mesh(V)
u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
u_grid.point_data["u"] = uh.x.array.real
u_grid.set_active_scalars("u")
u_plotter = pyvista.Plotter()
u_plotter.add_mesh(u_grid, show_edges=True)
u_plotter.view_xy()
u_plotter.show()

warped = u_grid.warp_by_scalar()
plotter2 = pyvista.Plotter()
plotter2.add_mesh(warped, show_edges=True, show_scalar_bar=True)
if not pyvista.OFF_SCREEN:
    plotter2.show()



gmsh.finalize()
ctx.gmsh_surf = []
ctx.gmsh_hole = []


gmsh.initialize()

"""
gmsh.model.occ.synchronize()
gdim = 2
gmsh.model.addPhysicalGroup(gdim, [membrane], 1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0.05)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax",0.05)
gmsh.model.mesh.generate(gdim)
"""
ctx.render_layers(pcb.bottom_layers, path + '/../pcb_bottom.png', theme.THEMES['OSH Park'], max_width=800, max_height=600)

surfaces_index = []
for s_iter in ctx.gmsh_surf:
	surface = s_iter[0]
	holes = s_iter[1]
	for hole in holes:
		#print([[(2, surface)],[(2,hole)]])
		s = gmsh.model.occ.cut([(2, surface)],[(2,hole)], removeObject=True, removeTool=False)
		#print(s)
		surface = s[0][0][1]
		gmsh.model.occ.remove([(2,hole)], recursive=True)
	for hole in ctx.gmsh_hole:
		#print([[(2, surface)],[(2,hole[0])]])
		s = gmsh.model.occ.cut([(2, surface)],[(2,hole[0])], removeObject=True, removeTool=False)
		#print(s)
		surface = s[0][0][1]
	surfaces_index.append((2, surface))
for hole in ctx.gmsh_hole:
	gmsh.model.occ.remove([(2,hole[0])], recursive=True)
surfaces = gmsh.model.occ.fragment(surfaces_index, [], removeObject=True, removeTool=True)

gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(2,[x[1] for x in surfaces[0]],0)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin",ctx.lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax",ctx.lc)
gmsh.model.mesh.generate(2)
gmsh.fltk.run()
gmsh.write("t2_.msh")


from dolfinx.io import gmshio
from mpi4py import MPI

gmsh_model_rank = 0
gdim = 2
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)
from dolfinx import fem
V = fem.FunctionSpace(domain, ("CG", 1))
import ufl
from petsc4py.PETSc import ScalarType
x = ufl.SpatialCoordinate(domain)
"""
set boundary conditions
"""

import re
import os
filepath = ""
for file in os.listdir(path):
    if(re.match(".*(-PTH.drl)$",file)):
        filepath = file

drills = []
holes = []
radius = 0
print(path + "/" + filepath)

with open(path + "/" + filepath) as f:
    for s_line in f:
        #print(s_line)
        result = re.match("T\d+C\d[.]\d+",s_line)
        if(result):
            print(s_line)
            print(s_line.split("C"))
            drills.append(s_line.split("C"))
        for drill in drills:
            if(re.match(drill[0],s_line)):
                print(s_line)
                radius = drill[1].split("\n")[0]
        if(re.match("X\d+[.]\d+Y[-]\d+[.]\d+",s_line)):
            print(s_line.split("X")[1].split("\n")[0].split("Y"))
            holes.append([radius, s_line.split("X")[1].split("\n")[0].split("Y")])

print(holes)


import numpy as np

bcs = []

for hole in holes:
	def on_boundary(x):
		return np.isclose(np.sqrt((x[0] - float(hole[1][0]) * 25.4)**2 + (x[1] - float(hole[1][1]) * 25.4)**2), float(hole[0]) * 25.4 / 2)
	boundary_dofs = fem.locate_dofs_geometrical(V, on_boundary)
	bcs.append(fem.dirichletbc(ScalarType(1), boundary_dofs, V))

"""
def on_boundary1(x):
    return np.isclose(np.sqrt((x[0] - 17.78)**2 + (x[1] + 20.32)**2), 0.2)
	#return np.isclose(np.sqrt((x[0] - 0.7*25.4)**2 + (x[1] - 0.8*25.4)**2), 0.157 * 25.4)
def on_boundary2(x):
    return np.isclose(np.sqrt((x[0] - 27.94)**2 + (x[1] + 20.32)**2), 0.2)
    #return np.isclose(np.sqrt((x[0] - 1.1*25.4)**2 + (x[1] - 0.8*25.4)**2), 0.157 * 25.4)

boundary_dofs1 = fem.locate_dofs_geometrical(V, on_boundary1)
bc1 = fem.dirichletbc(ScalarType(1), boundary_dofs1, V)

boundary_dofs2 = fem.locate_dofs_geometrical(V, on_boundary2)
bc2 = fem.dirichletbc(ScalarType(0), boundary_dofs2, V)
bcs = [bc1,bc2]
"""

f = fem.Constant(domain, ScalarType(0))
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx
problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

import pyvista
import dolfinx


u_topology, u_cell_types, u_geometry = dolfinx.plot.create_vtk_mesh(V)
u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
u_grid.point_data["u"] = uh.x.array.real
u_grid.set_active_scalars("u")
u_plotter = pyvista.Plotter()
u_plotter.add_mesh(u_grid, show_edges=True)
u_plotter.view_xy()
u_plotter.show()

warped = u_grid.warp_by_scalar()
plotter2 = pyvista.Plotter()
plotter2.add_mesh(warped, show_edges=True, show_scalar_bar=True)
if not pyvista.OFF_SCREEN:
    plotter2.show()



gmsh.finalize()
ctx.gmsh_surf = []
ctx.gmsh_hole = []

gmsh.initialize()

ctx.render_layers(pcb.copper_layers + pcb.drill_layers, path + '/../pcb_transparent_copper.png', theme.THEMES['Transparent Copper'], max_width=800, max_height=600)

gmsh.finalize()

