import h5py
import pandas as pd
from pathlib import Path


class ResultsFromDatabase:

    def __init__(self, filename: str, path: str = 'results_database'):

        self.filename = filename
        source_folder = Path(path)
        self.myfile = source_folder / filename

    def list_tallies(self):
        with h5py.File(self.myfile) as f:
            print(f.keys())

    def get_tally_dataframe(self, tally_name: str):
        with h5py.File(self.myfile) as f:
            return f[tally_name+'/table'][()]

    def get_tally_xaxis(self, tally_name: str):
        with h5py.File(self.myfile) as f:
            return f[tally_name+'/table'].attrs['x_axis']

    @property
    def literature_info(self):
        with h5py.File(self.myfile) as f:
            try:
                return f.attrs['literature']
            except KeyError:
                return 'n/a'

    @property
    def when(self):
        with h5py.File(self.myfile) as f:
            try:
                return f.attrs['when']
            except KeyError:
                return 'n/a'

    @property
    def where(self):
        with h5py.File(self.myfile) as f:
            try:
                return f.attrs['where']
            except KeyError:
                return 'n/a'

    @property
    def code_version(self):
        with h5py.File(self.myfile) as f:
            try:
                return f.attrs['code_version']
            except KeyError:
                return 'n/a'

    @property
    def xs_library(self):
        with h5py.File(self.myfile) as f:
            try:
                return f.attrs['xs_library']
            except KeyError:
                return 'n/a'

    @staticmethod
    def print_code_info(self):
        print(f'Code version:{self.code_version} \n',
              f'XS library: {self.xs_library} \n')

    def print_all_info(self):
        print(f'When: {self.when} \n',
              f'Where: {self.where} \n',
              f'Code version:{self.code_version} \n',
              f'XS library: {self.xs_library} \n'
              f'Literature: {self.literature_info} \n')
