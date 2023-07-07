import openmc

from dose_cells import dose_cell_ids
import data_config


model = openmc.Model.from_xml()
model.materials.clear()
cells = model.geometry.get_all_cells()
dose_cells = [cells[uid] for uid in dose_cell_ids]

dose_materials = []
for cell in dose_cells:
    cell.fill = cell.fill.clone()
    cell.fill.depletable = True
    dose_materials.append(cell.fill)

# Determine bounding box for volume calcuation
combined_region = openmc.Union([cell.region for cell in dose_cells])
lower_left, upper_right = combined_region.bounding_box

# Run a volume calculation
vol_calc = openmc.VolumeCalculation(dose_materials, 10_000_000, lower_left, upper_right)
model.settings.volume_calculations = [vol_calc]
model.calculate_volumes()

# Export new model
model.materials.clear()  # Make sure new materials in geometry get exported
model.export_to_xml()
