import openmc
from openmc.model import RightCircularCylinder as RCC
from openmc.model import RectangularParallelepiped as RPP


r0 = 7.5
r1 = 15
r2 = 30
r3 = 45
r4 = 48
r5 = 50
r6 = 60
r7 = 100

z0 = 0
z1 = 10
z2 = 110
z3 = 320
z4 = 645
z5 = 660
z6 = 690
z7 = 700

x0 = 150
y0 = 150

source_reg = -RCC((0, 0, z0), z1 - z0, r7)
frame_reg = +RCC((0, 0, z2), z5 - z2, r5) & -RCC((0, 0, z2), z5 - z2, r7)
plug_reg = +RCC((0, 0, z2), z3 - z2, r0) & -RCC((0, 0, z2), z3 - z2, r4)
back_plate_reg = -RCC((0, 0, z4), z5 - z4, r4)
void_reg0 = -RCC((0, 0, z1), z2 - z1, r7)
void_reg1 = +RCC((0, 0, z2), z5 - z2, r4) & -RCC((0, 0, z2), z5 - z2, r5)
void_reg2 = -RCC((0, 0, z2), z3 - z2, r0)
void_reg3 = -RCC((0, 0, z3), z4 - z3, r4)
void_reg4 = -RCC((0, 0, z5), z6 - z5, r7)
tally_reg0 = -RCC((0, 0, z6), z7 - z6, r1) 
tally_reg1 = +RCC((0, 0, z6), z7 - z6, r1) & -RCC((0, 0, z6), z7 - z6, r2)
tally_reg2 = +RCC((0, 0, z6), z7 - z6, r2) & -RCC((0, 0, z6), z7 - z6, r3)
tally_reg3 = +RCC((0, 0, z6), z7 - z6, r3) & -RCC((0, 0, z6), z7 - z6, r6)
tally_reg4 = +RCC((0, 0, z6), z7 - z6, r6) & -RCC((0, 0, z6), z7 - z6, r7)
domain_reg = -RPP(-x0, x0, -y0, y0, z0, z7, boundary_type='vacuum') & \
    +openmc.ZCylinder(r=r7)

source_cell = openmc.Cell(region=source_reg)
frame_cell = openmc.Cell(region=frame_reg)
plug_cell = openmc.Cell(region=plug_reg)
back_plate_cell = openmc.Cell(region=back_plate_reg)

void_cell0 = openmc.Cell(region=void_reg0)
void_cell1 = openmc.Cell(region=void_reg1)
void_cell2 = openmc.Cell(region=void_reg2)
void_cell3 = openmc.Cell(region=void_reg3)
void_cell4 = openmc.Cell(region=void_reg4)
void_cells = [void_cell0, void_cell1, void_cell2, void_cell3, void_cell4]

tally_cell0 = openmc.Cell(region=tally_reg0)
tally_cell1 = openmc.Cell(region=tally_reg1)
tally_cell2 = openmc.Cell(region=tally_reg2)
tally_cell3 = openmc.Cell(region=tally_reg3)
tally_cell4 = openmc.Cell(region=tally_reg4)
tally_cells = [tally_cell0, tally_cell1, tally_cell2, tally_cell3, tally_cell4]

domain_cell = openmc.Cell(region=domain_reg)

cells = [source_cell, frame_cell, plug_cell, back_plate_cell]
cells += void_cells + tally_cells + [domain_cell]

univ = openmc.Universe(cells=cells)
geom = openmc.Geometry(root=univ)
geom.export_to_xml(path='geometry2.xml', remove_surfs=True)
