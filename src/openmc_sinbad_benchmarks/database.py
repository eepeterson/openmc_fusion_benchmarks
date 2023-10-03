import h5py


class ResultFromDatabase:

    def __init__(self, filename: str):

        self.filename = filename

    def get_tally(self, tally_name: str):
        with h5py.File(self.filename) as f:
            return f[tally_name+'/table'][()]

    def get_columns(self, tally_name: str):
        pass

    def get_xaxis_label(self, tally_name: str):
        pass
