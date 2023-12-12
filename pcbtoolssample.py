import os
import sys

import gmsh
import numpy as np

sys.path.append(os.environ["HOME"]+"/pcb-tools")

from gerber import PCB
from gerber.render import theme
from gmsh_backend import GerberGmshContext
import dolfinx.fem.petsc
from dolfinx.io import gmshio
from mpi4py import MPI
import pyvista


path = "./sample1"

gmsh.initialize()
ctx = GerberGmshContext()
pcb = PCB.from_directory(path)

surfaces_index = []

ctx.render_layers(pcb.bottom_layers, path + '/../pcb_bottom.png', theme.THEMES['OSH Park'], max_width=800, max_height=600)
surfaces_index = []

print(ctx.gmsh_surf)
print(ctx.gmsh_hole)

for s_iter in ctx.gmsh_surf:
	surface = s_iter[0]
	holes = s_iter[1]
	for hole in holes:
		print([[(2, surface)],[(2,hole)]])
		s = gmsh.model.occ.cut([(2, surface)],[(2,hole)], removeObject=True, removeTool=False)
		print(s)
		surface = s[0][0][1]
		gmsh.model.occ.remove([(2,hole)], recursive=True)
	for hole in ctx.gmsh_hole:
		print([[(2, surface)],[(2,hole[0])]])
		s = gmsh.model.occ.cut([(2, surface)],[(2,hole[0])], removeObject=True, removeTool=False)
		print(s)
		surface = s[0][0][1]
	surfaces_index.append((2, surface))
for hole in ctx.gmsh_hole:
	gmsh.model.occ.remove([(2,hole[0])], recursive=True)
surfaces = gmsh.model.occ.fragment(surfaces_index, [], removeObject=True, removeTool=True)
gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(2,[x[1] for x in surfaces[0]],0)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin",ctx.lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax",ctx.lc)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)
gmsh.fltk.run()
gmsh.write("t2_.msh")



gmsh_model_rank = 0
gdim = 2
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)
from dolfinx import fem
import dolfinx

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
    if(re.match(".*(.DRD)$",file)):
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



potential = np.zeros(len(holes))
potential[1] = 1
boundaries = []
bcs = []

for i, hole in enumerate(holes):
	def on_boundary(x):
		return np.isclose(np.sqrt((x[0] - float(hole[1][0]) * 25.4)**2 + (x[1] - float(hole[1][1]) * 25.4)**2), float(hole[0]) * 25.4 / 2)
	boundary_dofs = fem.locate_dofs_geometrical(V, on_boundary)
	boundaries.append((i + 1, lambda x: np.isclose(np.sqrt((x[0] - float(hole[1][0]) * 25.4)**2 + (x[1] - float(hole[1][1]) * 25.4)**2), float(hole[0]) * 25.4 / 2)))
	bcs.append(fem.dirichletbc(ScalarType(potential[i]), boundary_dofs, V))

facet_indices, facet_markers = [], []
fdim = domain.topology.dim - 1

for (marker, locator) in boundaries:
	facets = dolfinx.mesh.locate_entities(domain, fdim, locator)
	facet_indices.append(facets)
	facet_markers.append(np.full_like(facets, marker))
facet_indices = np.hstack(facet_indices).astype(np.int32)
facet_markers = np.hstack(facet_markers).astype(np.int32)
sorted_facets = np.argsort(facet_indices)
facet_tag = dolfinx.mesh.meshtags(domain, fdim, facet_indices[sorted_facets], facet_markers[sorted_facets])


f = fem.Constant(domain, ScalarType(0))
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx
problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

sigma_Cu = 2.0115 * 10**3

n = ufl.FacetNormal(domain)
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tag)
current = dolfinx.fem.assemble_scalar(dolfinx.fem.form(ufl.dot(ufl.grad(u),n)*ds)) * sigma_Cu
print(current)




u_topology, u_cell_types, u_geometry = dolfinx.plot.vtk_mesh(V)
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
