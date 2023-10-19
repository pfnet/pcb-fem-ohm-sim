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
print(ctx.gmsh_surf)
print(ctx.gmsh_hole)
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
for surface in ctx.gmsh_surf:
	for hole in ctx.gmsh_hole:
		s = gmsh.model.occ.cut([(2, surface)],[(2,hole)], removeObject=True, removeTool=False)
		surface = s[0][0][1]
	surfaces_index.append((2, surface))

surfaces = gmsh.model.occ.fragment(surfaces_index, [], removeObject=True, removeTool=True)

gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(2,[x[1] for x in surfaces[0]],0)
gmsh.model.mesh.generate(2)
gmsh.fltk.run()
gmsh.write("t1_.msh")

ctx.gmsh_surf = []
ctx.gmsh_hole = []


"""
gmsh.model.occ.synchronize()
gdim = 2
gmsh.model.addPhysicalGroup(gdim, [membrane], 1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0.05)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax",0.05)
gmsh.model.mesh.generate(gdim)
"""
ctx.render_layers(pcb.bottom_layers, path + '/../pcb_bottom.png', theme.THEMES['OSH Park'], max_width=800, max_height=600)
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)
gmsh.fltk.run()
gmsh.write("t2_.msh")
ctx.render_layers(pcb.copper_layers + pcb.drill_layers, path + '/../pcb_transparent_copper.png', theme.THEMES['Transparent Copper'], max_width=800, max_height=600)


