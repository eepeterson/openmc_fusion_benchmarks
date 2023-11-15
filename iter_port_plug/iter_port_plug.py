import numpy as np
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

# Build regions
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

# water, formula: H2O, density: 1.0 g/cm3
water = openmc.Material(name='water')
water.add_elements_from_formula('H2O')
water.set_density('g/cm3', 1.0)

# Stainless Steel 316L
steel = openmc.Material(name="steel_316L")
steel.add_element("Fe", 62.045, "wo")
steel.add_element("Cr", 18., "wo")
steel.add_element("Ni", 14., "wo")
steel.add_element("Mo", 3., "wo")
steel.add_element("Mn", 2., "wo")
steel.add_element("Si", 0.75, "wo")
steel.add_element("N", 0.1, "wo")
steel.add_element("C", 0.03, "wo")
steel.add_element("P", 0.045, "wo")
steel.add_element("S", 0.03, "wo")
steel.set_density("g/cm3", 7.93) 

steel_water_mix = openmc.Material.mix_materials([steel, water],
                                                fracs=[0.7988456, 0.2011544],
                                                percent_type='vo')

# Build cells
source_cell = openmc.Cell(region=source_reg)
frame_cell = openmc.Cell(region=frame_reg, fill=steel)
plug_cell = openmc.Cell(region=plug_reg, fill=steel_water_mix)
back_plate_cell = openmc.Cell(region=back_plate_reg, fill=steel)

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

# Make sure ZPlanes at z=0 and z=700 are vacuum boundaries
for cell in cells:
    for surf in cell.region.get_surfaces().values():
        if isinstance(surf, openmc.ZPlane):
            if np.isclose(surf.z0, 0) or np.isclose(surf.z0, z7):
                surf.boundary_type = 'vacuum'

univ = openmc.Universe(cells=cells)
geometry = openmc.Geometry(root=univ)
geometry.merge_surfaces = True

# Settings
settings = openmc.Settings(run_mode='fixed source',
                           particles=int(1e6),
                           batches=1000,
                           weight_windows_on=True)
energy = openmc.stats.Discrete([14.08e6], [1])
zspace = openmc.stats.Uniform(z0, z1)
phispace = openmc.stats.Uniform(0, 2*np.pi)
rspace = openmc.stats.PowerLaw(0, r7, n=1)
space = openmc.stats.CylindricalIndependent(rspace, phispace, zspace)
source = openmc.IndependentSource(space=space, energy=energy)
settings.source = source
settings.weight_windows = openmc.wwinp_to_wws('iter_port_plug_cadis.wwinp')

# Build filters and tallies
mesh = openmc.CylindricalMesh(r_grid=np.linspace(0, r7, 201),
                              z_grid=np.linspace(z0, z7, 176))
energy_filter = openmc.EnergyFilter.from_group_structure('VITAMIN-J-175')
mesh_filter = openmc.MeshFilter(mesh)
particle_filter = openmc.ParticleFilter(['neutron'])

t1 = openmc.Tally()
t1.filters = [particle_filter, mesh_filter, energy_filter]
t1.scores = ['flux']
tallies = openmc.Tallies([t1])

model = openmc.model.Model(geometry=geometry,
                           settings=settings,
                           tallies=tallies)
model.run()
