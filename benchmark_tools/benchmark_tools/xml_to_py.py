import openmc
import re


def print_header(filename="openmc_model.py"):
    """creates the openmc Python API input file and prints the header"""
    with open(str(filename), "w", encoding="utf-8") as fh:
        fh.write("#%%\n")
        fh.write("import openmc\n")


def materials_to_py(
    filename="openmc_model.py", xml_file="materials.xml", export_to_xml=True
):
    """reads materials.xml file and gives back the corresponding Python API commands"""
    # read the xml file
    materials = openmc.Materials.from_xml(xml_file)

    # open the Python API file to write the model
    with open(str(filename), "a", encoding="utf-8") as fh:
        # insert comments in the file
        fh.write("\n#%%\n")
        fh.write("\n# MATERIALS\n\n")

        # write materials in the Python file
        for m in materials:
            # store material id, density etc.
            m_id = str(m.id)
            m_d = str(m.get_mass_density())
            m_du = str(m.density_units)
            m_nuc = m.get_nuclide_densities()
            matname = "mat_" + m_id
            # instantiate openmc material in the Python file
            fh.write("# " + m.name + "\n")
            fh.write(matname + " = openmc.Material(material_id=" + m_id)
            if m.name:
                fh.write(", name='" + str(m.name) + "'")
            fh.write(")\n")
            fh.write(matname + ".set_density('" + m_du + "', " + m_d + ")\n")
            # add nuclides, density and others to each material
            for k in m_nuc.keys():
                p = str(m_nuc[k].percent)
                pt = str(m_nuc[k].percent_type)
                fh.write(
                    matname + ".add_nuclide('" + k + "', " + p + ", '" + pt + "')\n"
                )
            if m.volume:
                m_v = str(m.volume)
                fh.write(matname + ".volume" + " = " + m_v + "\n")
            if m.temperature:
                m_t = str(m.temperature)
                fh.write(matname + ".temperature" + " = " + m_t + "\n")

        # collect all materials in openmc.Materials
        fh.write("\n# create materials instance\n")
        line = "materials = openmc.Materials(["
        for m in materials:
            m_id = str(m.id)
            line += "mat_" + m_id + ", "
        line = line[:-2] + "])\n"

        fh.write(line)

        if export_to_xml:
            fh.write("\nmaterials.export_to_xml()\n")


def geometry_to_py(
    filename="openmc_model.py", xml_file="geometry.xml", export_to_xml=True
):
    """reads geometry.xml file and gives back the corresponding Python API commands"""

    surfs_dict = {
        "x-plane": "XPlane",
        "y-plane": "YPlane",
        "z-plane": "ZPlane",
        "plane": "Plane",
        "quadric": "Quadric",
        "x-cylinder": "XCylinder",
        "y-cylinder": "YCylinder",
        "z-cylinder": "ZCylinder",
        "sphere": "Sphere",
        "x-cone": "XCone",
        "y-cone": "YCone",
        "z-cone": "ZCone",
        "x-torus": "XTorus",
        "y-torus": "lYTorus",
        "z-torus": "ZTorus",
    }

    # read xml file
    geometry = openmc.Geometry.from_xml(xml_file)

    # get surfaces, cells and universes
    # currently does not supports lattices and macrobodies
    surfs = geometry.get_all_surfaces()
    cells = geometry.get_all_cells()
    universes = geometry.get_all_universes()

    # open the Python API file to write the model
    with open(str(filename), "a", encoding="utf-8") as fh:
        # insert comments in the file
        fh.write("\n#%%\n")
        fh.write("\n# GEOMETRY\n")

        # write surfaces in the Python file
        fh.write("\n# surfaces\n")
        for surf in surfs:
            # collect surfaces id, type etc.
            s = surfs[surf]
            s_id = str(s.id)
            s_type = s.type
            s_bt = s.boundary_type
            s_coeffs = s.coefficients

            # instantiate openmc surface in the Python file
            surf_class = surfs_dict[s_type]
            line = "surf_" + s_id + " = "
            line += "openmc." + surf_class + "("
            line += "surface_id=" + s_id + ", "

            # add surface coefficients
            for k in s_coeffs.keys():
                line += k + "=" + str(s_coeffs[k])
                line += ", "

            # ad boundary type
            if s_bt != "transmission":
                line += "boundary_type='"
                line += s_bt + "'" + ", "

            # add name if applicable
            if s.name:
                line += "name="
                line += "'" + s.name + "'" + ", "

            line = line[:-2] + ")\n"

            fh.write(line)

        # write regions in the Python file
        fh.write("\n# regions\n")
        for i, c in cells.items():
            # collect cells/region info
            c_id = str(c.id)
            c_region = str(c.region)

            # reconstruct csg commands
            reg_surf_ids = re.findall(r"\b\d+\b", c_region)
            c_region = re.sub(r"\b\d+\b", "+surf_{}", c_region).format(*reg_surf_ids)
            c_region = (
                c_region.replace("-+", "-").replace(" ", " & ").replace("& | &", "|")
            )

            line = "region_" + c_id
            line += " = " + c_region + "\n"

            fh.write(line)

        # write cells in the Python file
        fh.write("\n# cells\n")
        for i, c in cells.items():
            # collect cells id, material etc.
            c_id = str(c.id)

            # instantiate openmc cell in the Python file
            line = "cell_" + c_id + " = " + "openmc.Cell("
            line += "cell_id=" + c_id + ", "
            line += "region=region_" + c_id + ", "

            # fill the cell with the material
            line += "fill="
            if c.fill:
                c_mat_id = str(c.fill.id)
                line += "mat_" + c_mat_id
            else:
                line += "None"
            line += ", "

            # add optional elements
            if c.name:
                line += "name=" + "'" + c.name + "'" + ", "
            if c.temperature:
                line += "temperature=" + c.temperature + ", "
            if c.volume:
                line += "volume=" + c.volume + ", "

            line = line[:-2] + ")\n"

            fh.write(line)

        # collect all cells in the root universe
        line = "universe = openmc.Universe(cells=["

        for i, c in cells.items():
            c_id = str(c.id)
            line += "cell_" + c_id + ", "
        line = line[:-2] + "])\n"

        fh.write("\n# create root universe\n")
        fh.write(line)

        fh.write("\n# create geometry instance\n")
        fh.write("geometry = openmc.Geometry(universe)")

        if export_to_xml:
            fh.write("\ngeometry.export_to_xml()\n")


if __name__ == "__main__":
    print_header()
    materials_to_py()
    geometry_to_py()
