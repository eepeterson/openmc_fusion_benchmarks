import openmc
import openmc.lib


s19 = openmc.YCylinder(0, 0, 17.7)
s48 = openmc.Plane(0.7071067811865476, 6.123233995736766e-17, 0.7071067811865476, 11.45)
s58 = openmc.Plane(0.7071067811865476, 6.123233995736766e-17, 0.7071067811865476, 14.35)
s62 = openmc.Plane(6.123233995736766e-17, 1.0, 6.123233995736766e-17, 1.5999999999999999)
s64 = openmc.Plane(6.123233995736766e-17, 1.0, 6.123233995736766e-17, -1.3)
s65 = openmc.Plane(-0.7071067811865475, 6.123233995736766e-17, 0.7071067811865476, 1.45)
s68 = openmc.Plane(-0.7071067811865475, 6.123233995736766e-17, 0.7071067811865476, -1.45)
s101 = openmc.YPlane(5.6, boundary_type='vacuum')

region = +s64 & -s101 & -s19 & +s48 & (+s58 | +s62 | -s64 | +s65 | -s68)
cell = openmc.Cell(region=region)

univ = openmc.Universe(cells=[cell])
geometry = openmc.Geometry(root=univ)
settings = openmc.Settings(particles=100, batches=100)
model = openmc.model.Model(geometry=geometry, settings=settings)
model.export_to_model_xml()

openmc.lib.init()
print(cell.bounding_box)
print(openmc.lib.cells[1].bounding_box)
openmc.lib.finalize()
