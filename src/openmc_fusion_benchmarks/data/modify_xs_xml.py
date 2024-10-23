import openmc
import os


_this_script_dir = os.path.dirname(os.path.abspath(__file__))
_model_xml_path = os.path.join(_this_script_dir, 'model_xs.xml')

def get_env_variable(var_name):
    # Retrieve the value of the environment variable
    value = os.getenv(var_name)
    if value is None:
        return f"Environment variable '{var_name}' is not set."
    return value


def rewrite_xs_xml(model_xs_file:str=_model_xml_path, new_xs_file:str='cross_sections_mod.xml'):
    # read xs file
    model_xs_file = get_env_variable('OPENMC_CROSS_SECTIONS')
    myxs = openmc.data.DataLibrary.from_xml(model_xs_file)
    myxs.export_to_xml(new_xs_file)


def create_xs_dict(xs_h5_file:str, nuclide: str, type:str='neutron'):
    return {'path': xs_h5_file, 'type': type, 'materials': [nuclide]}

def perturb_xs_xml(xs_file, xs_h5_file, nuclide: str):
    # read xs file
    myxs = openmc.data.DataLibrary.from_xml(xs_file)
    try:
        myxs.libraries.get_by_material(nuclide)['path'] = xs_h5_file
    except TypeError:
        myxs.append(create_xs_dict(xs_h5_file, nuclide))
    # export to modified cross_sections xml file
    myxs.export_to_xml(xs_file)