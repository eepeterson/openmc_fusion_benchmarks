import sys
import argparse
from pathlib import Path


def main():

    parser = argparse.ArgumentParser(
        description='Input for Cubit conversion script.')
    parser.add_argument("abaqus_file", type=Path,
                        help='Path to the Abaqus input file')
    parser.add_argument("-p", "--cubit_path", type=Path, required=False,
                        help="Path to Cubit executable (required if PATH does not include Cubit)")
    parser.add_argument("-o", "--output", type=Path, required=False, default="exodus_out.exo",
                        help="Name/Location of the output file (Default is current/directory/exodus_out.exo)")
    parser.add_argument("-f", "--force", action='store_true', default=False,
                        help="If present, overwrite a pre-existing output files.")

    # Should either add cubit to path beforehand or require input for the line below
    # IMP Could add try/catch
    # IMP could try to avoid importing sys
    args = parser.parse_args()

    # Makes sure cubit is importable
    if args.cubit_path is not None:
        sys.path.append(str(args.cubit_path))

    convert_groups_to_blocks(args.abaqus_file, args.output, args.force)

# if __name__ == "__main__":
#     main()
