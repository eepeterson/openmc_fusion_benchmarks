import h5py
import pandas as pd


class TallyFromDatabase:
    def __init__(self, filename: str, tally_name: str):
        self.filename = filename
        self.tally_name = tally_name

        with h5py.File(self.filename) as f:
            self.tally = f[tally_name+'/table']

    def get_pandas_dataframe(self):
        return pd.DataFrame(self.tally[()])

    @property
    def get_xaxis_label(self):
        return self.tally.attrs['x_axis']


class ResultFromDatabase:

    def __init__(self, filename: str):

        self.filename = filename

    def list_tallies(self):
        with h5py.File(self.filename) as f:
            print(f.keys())

    def get_tally(self, tally_name: str):
        return TallyFromDatabase(self.filename, tally_name)

    @property
    def literature_info(self):
        with h5py.File(self.filename) as f:
            return f.attrs['literature']

    @property
    def when(self):
        with h5py.File(self.filename) as f:
            return f.attrs['when']

    @property
    def where(self):
        with h5py.File(self.filename) as f:
            return f.attrs['where']

    @property
    def code_version(self):
        with h5py.File(self.filename) as f:
            return f.attrs['code_version']

    @property
    def xs_library(self):
        with h5py.File(self.filename) as f:
            return f.attrs['xs_library']

    @property
    def print_code_info(self):
        print(f'Code version:{self.code_version} \n',
              f'XS library: {self.xs_library} \n')

    def print_all_info(self):
        print(f'When: {self.when} \n',
              f'Where: {self.where} \n',
              f'Code version:{self.code_version} \n',
              f'XS library: {self.xs_library} \n'
              f'Literature: {self.literature_info} \n')
