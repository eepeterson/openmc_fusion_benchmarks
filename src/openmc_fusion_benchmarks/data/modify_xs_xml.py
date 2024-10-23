import openmc
import os


_this_script_dir = os.path.dirname(os.path.abspath(__file__))
_model_xml_path = os.path.join(_this_script_dir, 'model_xs.xml')


def get_env_variable(var_name: str) -> str:
    """Get the value of an environment variable.

    Parameters
    ----------
    var_name : str
        The name of the environment variable.

    Returns
    -------
    str
        The value of the environment variable.
    """
    value = os.getenv(var_name)
    if value is None:
        return f"Environment variable '{var_name}' is not set."
    return value


def rewrite_xs_xml(new_xs_file: str = 'cross_sections_mod.xml'):
    """Rewrite the cross_sections.xml file locally, in order not to modify the
    original file.

    Parameters
    ----------
    new_xs_file : str, optional
        name of the new xs file to be created,
        by default 'cross_sections_mod.xml'
    """
    # read xs file
    model_xs_file = get_env_variable('OPENMC_CROSS_SECTIONS')
    myxs = openmc.data.DataLibrary.from_xml(model_xs_file)
    myxs.export_to_xml(new_xs_file)


def create_xs_dict(xs_h5_file: str, nuclide: str, type: str = 'neutron') -> dict:
    """Create a dictionary with the correct format for an cross_sections.xml file
    for a given target nuclide and the path to its corresponding h5 xs file.


    Parameters
    ----------
    xs_h5_file : str
        path to the target nuclide h5 xs file
    nuclide : str
        nuclide name in gnds format
    type : str, optional
        type of incident particle (i.e. neutron, photon), by default 'neutron'

    Returns
    -------
    dict
        dictionary with the path to the h5 xs file, the type of incident
        particle and the target nuclide
    """
    return {'path': xs_h5_file, 'type': type, 'materials': [nuclide]}


def perturb_xs_xml(xs_file: str, xs_h5_file: str, nuclide: str):
    """Perturb the cross_sections.xml file by replacing the original path of a
    nuclide xs with a new path to a new (perturbed) h5 xs file.

    Parameters
    ----------
    xs_file : str
        name of the xs xml file to perturb
    xs_h5_file : str
        path to the new (perturbed) h5 xs file
    nuclide : str
        nuclide name in gnds format
    """
    # read xs file
    myxs = openmc.data.DataLibrary.from_xml(xs_file)
    try:
        myxs.libraries.get_by_material(nuclide)['path'] = xs_h5_file
    except TypeError:
        myxs.append(create_xs_dict(xs_h5_file, nuclide))
    # export to modified cross_sections xml file
    myxs.export_to_xml(xs_file)
