import openmc
from pathlib import Path


class ResultsFromOpenmc:

    def __init__(self, filename: str, path: str):

        self.filename = filename
        source_folder = Path(path)
        self.myfile = source_folder / filename

    def get_tally_dataframe(self, tally_name):
        pass
