import json
import gdown
import importlib


LIB_PATH = importlib.resources.files(
    "openmc_fusion_benchmarks.lib")


def download_geometry(benchmark_name: str, file_format: str, run_option: str = None):

    filepath = LIB_PATH / "dagmc_geometries.json"
    with open(filepath, "r") as f:
        data = json.load(f)

        # Your Google Drive file link
    if run_option is not None:
        url = data[benchmark_name][run_option][file_format]
    else:
        url = data[benchmark_name][file_format]

    # Extract the file ID from the URL
    file_id = url.split("/d/")[1].split("/")[0]
    download_url = f"https://drive.google.com/uc?export=download&id={file_id}"

    # Download the file
    gdown.download(download_url, quiet=False, use_cookies=False)
