import h5py
import openmc
from pathlib import Path
from typing import Iterable


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


class ResultsFromOpenmc:

    def __init__(self, filename: str, path: str):

        self.filename = filename
        source_folder = Path(path)
        self.myfile = source_folder / filename
        self.statepoint = openmc.StatePoint(self.myfile)

    def get_tally_dataframe(self, tally_name: str, normalize_over: Iterable = None):

        tally_dataframe = self.statepoint.get_tally(
            tally_name).get_pandas_dataframe()

        if normalize_over:
            tally_dataframe['mean'] = tally_dataframe['mean'] / normalize_over
            tally_dataframe['std. dev.'] = tally_dataframe['std. dev.'] / \
                normalize_over

        return tally_dataframe

    @property
    def get_openmc_version(self):
        return self.statepoint.version

    @property
    def get_particles_per_batch(self):
        return format(self.statepoint.n_particles, '.2e')

    @property
    def get_batches(self):
        return self.statepoint.n_batches

    def tally_to_hdf(self, filename: str, tally_name: str, normalize_over: Iterable, xs_library: str, x_axis: str = None, path_to_database: str = '../results_database'):

        path = Path(path_to_database)
        path = path / filename

        tally_df = self.get_tally_dataframe(
            tally_name, normalize_over=normalize_over)

        tally_df.to_hdf(filename, tally_name, mode='a',
                        format='table', data_columns=True, index=False)

        code_version = 'openmc-' + '.'.join(map(str, self.get_openmc_version))

        with h5py.File(filename, 'a') as f:
            f[tally_name + '/table'].attrs['x_axis'] = x_axis
            f.attrs['code_version'] = code_version
            f.attrs['xs_library'] = xs_library
            f.attrs['batches'] = self.get_batches
            f.attrs['particles_per_batch'] = self.get_particles_per_batches
