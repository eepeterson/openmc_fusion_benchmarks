import json
import openmc

model = openmc.Model.from_model_xml('converted_model.xml')
model.geometry.merge_surfaces = True

settings = openmc.Settings(run_mode='fixed source',
                           particles=1_000_000,
                           batches=50,
                           output={'tallies':False})

source = openmc.IndependentSource()
source.energy = openmc.stats.Discrete([14.1e6], [1])
source.space = openmc.stats.Point((0, 0.001, 0))
settings.source = source

wws = openmc.wwinp_to_wws('fng_dose.fwcadis.wwinp')
settings.weight_windows = wws
settings.weight_windows_on = True

model.settings = settings

def apply_volumes(model, filename='../cell_volumes.json', material=False):
    # Read volumes from JSON
    with open(filename, 'r') as fh:
        volumes = json.load(fh)

    # Add volumes to cells/materials
    cells = model.geometry.get_all_cells()
    for uid, volume in volumes.items():
        cell = cells[int(uid)]
        if material:
            cell.fill = cell.fill.clone()
            cell.fill.depletable = True
            cell.fill.volume = volume
        else:
            cell.volume = volume

    # Make sure new materials in geometry get exported
    if material:
        model.materials.clear()

apply_volumes(model)

bc_surfs = [21, 92, 93, 94, 95, 124]
cells = model.geometry.get_all_cells()
for surf in cells[654].region.get_surfaces().values():
    surf.boundary_type = 'vacuum'

cells.pop(229)
cells.pop(241)
cells.pop(654)
cell_ids = list(cells.keys())

cell_filter = openmc.CellFilter(cell_ids)
energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')

flux_tally = openmc.Tally(name='flux_tally')
flux_tally.scores = ['flux']
flux_tally.filters = [cell_filter, energy_filter]
model.tallies = openmc.Tallies([flux_tally])

model.export_to_model_xml('model.xml')
