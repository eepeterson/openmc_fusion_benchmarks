from .data_conventions import get_nuclide_zaid, get_nuclide_gnds
import os
import glob
import argparse
import sandy
import openmc
import openmc.data


def remove_ace_files(directory, lib_name):
    # List of extensions to remove
    extensions = ['03c', '03c.xsd', lib_name]

    # Loop through each extension and remove matching files
    for ext in extensions:
        files_to_remove = glob.glob(os.path.join(directory, f'*.{ext}'))
        for file_path in files_to_remove:
            try:
                os.remove(file_path)
                print(f"Removed {file_path}")
            except Exception as e:
                print(f"Error removing {file_path}: {e}")


def get_ace_files(nsamples, lib_name, nuclide, reaction, nprocesses, error):
    # convert nuclide to ZAID
    nuclide_zaid = get_nuclide_zaid(nuclide)
    tape = sandy.get_endf6_file(lib_name, "xs", nuclide_zaid*10, to_file=True)

    # Generate perturbations
    perturbations = tape.get_perturbations(
        nsamples,
        njoy_kws=dict(
            err=error,
            chi=False,
            mubar=False,
            xs=True,
            nubar=False,
            verbose=True,
            errorr33_kws=dict(mt=reaction)
        ),
    )

    # Apply perturbations and generate ACE files
    ace_files = tape.apply_perturbations(
        perturbations,
        processes=nprocesses,
        njoy_kws=dict(err=error),
        to_ace=True,
        to_file=True,
        ace_kws=dict(
            err=error,
            temperature=294,
            verbose=True,
            purr=True,
            heatr=True,
            thermr=True,
            gaspr=True,
            groupr=True,
            errorr=True,
            heatr_kws={'local': True}
        ),
        verbose=True,
    )

    # Optionally return or process 'outs' if needed
    return ace_files


def ace_to_hdf5(nsamples, lib_name, nuclide, remove_ace=True):
    # convert nuclide to gnds and ZAID
    nuclide_gnds = get_nuclide_gnds(nuclide)
    nuclide_zaid = get_nuclide_zaid(nuclide)

    directory = f"{nuclide_gnds}_{lib_name}"
    if not os.path.exists(directory):
        os.makedirs(directory)

    for n in range(0, nsamples):
        acefile = f"{nuclide_zaid}_{n}.03c"
        h5file = f"{directory}/{nuclide_gnds}_{n}_{lib_name}.h5"

        # Convert ACE to HDF5 using OpenMC
        try:
            nuc_data = openmc.data.IncidentNeutron.from_ace(acefile)
            print(f'Writing to {os.getcwd()}/{h5file}')
            nuc_data.export_to_hdf5(h5file)
        except FileNotFoundError:
            print(f'Error: ACE file {acefile} not found.')
        except Exception as e:
            print(f'Error processing {acefile}: {e}')

    if remove_ace:
        remove_ace_files('', lib_name)


def perturb_to_hdf5(nsamples, lib_name, nuclide, reaction, nprocesses=1, error=.001):
    get_ace_files(nsamples, lib_name, nuclide, reaction, nprocesses, error)
    ace_to_hdf5(nsamples,  lib_name, nuclide, remove_ace=True)


def main():
    parser = argparse.ArgumentParser(
        description="Generate ACE files with perturbed nuclear data.")
    parser.add_argument("-n", "--nuclide", required=True,
                        help="Nuclide identifier (ZA number)")
    parser.add_argument("-xs", "--lib_name", type=str, default=1,
                        help="Cross section choice ['JEFF_32', 'JEFF_33', 'ENDFB_71', 'ENDFB_80', 'JENDL_40U', 'IRDFF_2']")
    parser.add_argument("-r", "--reaction", type=int,
                        nargs='+', help="MT value for specific reaction")
    parser.add_argument("-e", "--error", type=float,
                        default=0.001, help="Error tolerance for processing")
    parser.add_argument("-p", "--nprocesses", type=int, default=1,
                        help="Number of processes for parallel processing")
    parser.add_argument("-ns", "--nsamples", type=int,
                        default=1, help="Number of perturbations")

    args = parser.parse_args()

    # Call function with arguments from command line
    perturb_to_hdf5(args.nsamples, args.lib_name, args.nuclide, args.reaction,
                    args.nprocesses)


if __name__ == "__main__":
    main()
