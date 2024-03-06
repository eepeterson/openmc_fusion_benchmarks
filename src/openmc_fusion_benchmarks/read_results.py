import h5py
import openmc
from pathlib import Path
from typing import Iterable
import pandas as pd

_del_columns = ['cell', 'particle', 'nuclide', 'score', 'energyfunction']


def build_hdf_filename(software_name: str, software_version: Iterable, xs_library: str) -> str:
    """Builds the name for the hdf file to be stored in a results_database folder

    Parameters
    ----------
    software_name : str
        name of the software used for the simulation
    software_version : Iterable
        list or tuple with the software version
    xs_library : str
        name of the nuclear data library used in the simulation

    Returns
    -------
    str
        name of the hdf file
    """

    filename = software_name.strip().replace(
        ' ', '').replace('.', '').replace('-', '').lower()
    filename += '-' + '-'.join(map(str, software_version)) + '_'
    filename += xs_library.strip().replace(' ', '').replace('.',
                                                            '').replace('-', '').lower()
    filename += '.h5'

    return filename


class ResultsFromDatabase:
    """This class takes in a hdf file and its path and generates a generic
    object by reading it. It is specifically designed for hdf files present
    in the "results_database" folders of the benchmark models in order to be
    able to read, postprocess and plot the results from previous simulations
    that have been stored there
    """

    def __init__(self, filename: str, path: str = ''):
        """ResultsFromDatabase class constructor

        Parameters
        ----------
        filename : str
            Name of the hdf file present in the results_database folder
        path : str, optional
            path to the hdf file, by default 'results_database'
        """

        self.filename = filename
        source_folder = Path(path)

        # merge path to the file
        self._myfile = source_folder / filename

    def list_tallies(self):
        """Prints the names of all the tallies available in the hdf file
        """
        with h5py.File(self._myfile) as f:
            print(f.keys())

    def get_tally_dataframe(self, tally_name: str):
        """Retrieves the results of a given tally in a Pandas DataFrame format.
        It relies on the openmc.Statepoint().get_tally().get_pandas_dataframe()
        method

        Parameters
        ----------
        tally_name : str
            Exact name of the tally in the hdf file 

        Returns
        -------
        pandas.DataFrame
            DataFrame with tally results
        """
        with h5py.File(self._myfile) as f:
            df = pd.DataFrame(f[tally_name+'/table'][()]).drop(columns='index')
            # decode hdf5 strings to strings if necessary
            try:
                df[self.get_tally_xaxis(tally_name)] = [el.decode()
                                                        for el in df[self.get_tally_xaxis(tally_name)]]
            except:
                pass

            return df

    def get_tally_xaxis(self, tally_name: str):
        """Retrieves the string with the exact name associated to the tally
        x-axis. It is necessary because each result of each benchmark has a
        different data classification that can vary from depth in the shield,
        energy bin, distance from source, detector position etc. Retrieving this
        string allows to directly use it as string to get the column with that lists
        the data classification (x-axis for plotting purposes)

        Parameters
        ----------
        tally_name : str
            Exact name of the tally in the hdf file 

        Returns
        -------
        str
            Name used for the dataframe column with the x-axis info
        """
        with h5py.File(self._myfile) as f:
            return f[tally_name+'/table'].attrs['x_axis']

    @property
    def literature_info(self):
        """Retrieves the peer reviewed publication(s) associated with the data
        provided in the hdf file

        Returns
        -------
        str
            DOI or link to the publication
        """
        with h5py.File(self._myfile) as f:
            try:
                return f.attrs['literature']
            except KeyError:
                return 'n/a'

    @property
    def when(self):
        """Retrieves the experiment or simulation place if provided

        Returns
        -------
        str
            Experiment or simulation execution year
        """
        with h5py.File(self._myfile) as f:
            try:
                return f.attrs['when']
            except KeyError:
                return 'n/a'

    @property
    def where(self):
        """Retrieves the experiment or simulation place if provided

        Returns
        -------
        str
            Experiment or simulation execution place
        """
        with h5py.File(self._myfile) as f:
            try:
                return f.attrs['where']
            except KeyError:
                return 'n/a'

    @property
    def code_version(self):
        """Retrieves the code version if the hdf file refers to a simulation's
        results and the info was provided

        Returns
        -------
        str
            Code version
        """
        with h5py.File(self._myfile) as f:
            try:
                return f.attrs['code_version']
            except KeyError:
                return 'n/a'

    @property
    def xs_library(self):
        """Retrieves the nuclear data library used if the hdf file refers
        to a simulation's results and the info was provided

        Returns
        -------
        str
            Nuclear data library name
        """
        with h5py.File(self._myfile) as f:
            try:
                return f.attrs['xs_library']
            except KeyError:
                return 'n/a'

    def print_code_info(self):
        """Prints all the code info available
        """
        print(f'Code version:{self.code_version}\n',
              f'XS library: {self.xs_library}\n')

    def print_all_info(self):
        """Prints all info available
        """
        print(f'Info:\n',
              f'When: {self.when}\n',
              f'Where: {self.where}\n',
              f'Code version: {self.code_version}\n',
              f'XS library: {self.xs_library}\n',
              f'Literature: {self.literature_info}\n')


class ResultsFromOpenmc:
    """Similarly to ResultsFromDatabase, this class creates an object
    containing the results of a fresh new openmc simulation. It extracts
    from an openmc statepoint.h5 file.
    If openmc results have already been stored in an hdf file in the
    results_database folder it is necessary to use ResultsFromDatabase class
    """

    def __init__(self, statepoint_file: str = 'statepoint.100.h5', path: str = ''):
        """ResultsFromOpenmc class constructor

        Parameters
        ----------
        statepoint_file : str, optional
            name of the statepoint.h5 file, by default 'statepoint.100.h5'
        path : str, optional
            path to the statepoint file, by default 'results'
        """
        self.statepoint_file = statepoint_file
        source_folder = Path(path)
        # merge path to statepoint file
        self._myfile = source_folder / statepoint_file
        # open statepoint file with openmc
        self.statepoint = openmc.StatePoint(self._myfile)

    def get_tally_dataframe(self, tally_name: str, normalize_over: Iterable = None):
        """Retrieves the results of a given tally in a Pandas DataFrame format.
        It relies on the openmc.Statepoint().get_tally().get_pandas_dataframe()
        method

        Parameters
        ----------
        tally_name : str
            Exact name of the tally as defined in the openmc model
        normalize_over : Iterable, optional
            Some openmc tallies (e.g. cell tally, surface tally) need to be normalized 
            by their filter dimension (e.g cell volume, surface area), by default None

        Returns
        -------
        pandas.DataFrame
            DataFrame with tally results
        """
        # extract tally in dataframe format from statepoint file
        tally_dataframe = self.statepoint.get_tally(
            name=tally_name).get_pandas_dataframe()

        # normalize tally over tally filter dimension (e.g. cell volume, surface area etc.)
        if normalize_over is not None:
            tally_dataframe['mean'] = tally_dataframe['mean'] / normalize_over
            tally_dataframe['std. dev.'] = tally_dataframe['std. dev.'] / \
                normalize_over

        return tally_dataframe

    @property
    def get_openmc_version(self):
        """Retrieves openmc's version used in the simulation

        Returns
        -------
        tuple
            openmc version
        """
        return self.statepoint.version

    @property
    def get_particles_per_batch(self):
        """Retrieves the number of particles per batch defined in the simulation

        Returns
        -------
        float
            Number of particle per batch
        """
        return format(self.statepoint.n_particles, '.2e')

    @property
    def get_batches(self):
        """Retrieves the number of batches defined in the simulation

        Returns
        -------
        int
            Number of batches
        """
        return self.statepoint.n_batches

    def tally_to_hdf(self, tally_name: str, normalize_over: Iterable, xs_library: str, xaxis_name: str,
                     xaxis_list: Iterable = None, path_to_database: str = '../results_database', when: str = 'n/a',
                     where: str = 'n/a'):
        """Stores the openmc tally in a hdf file for the results_database folder

        Parameters
        ----------
        tally_name : str
            Exact name of the tally to store in the hdf file. It should be the
            same as the openmc tally
        normalize_over : Iterable
            Some openmc tallies (e.g. cell tally, surface tally) need to be normalized 
            by their filter dimension (e.g cell volume, surface area), by default None
        xs_library : str
            Name of the nuclear data library used for the simulation
        xaxis_name : str, optional
            name of the x_axis to store in order to be retrieved with the
            ResultsFromDatabase.get_tally_xaxis method, by default None
        xaxis_list : Iterable
            list of elements to apply to the df columns identified by the xaxis_name
        path_to_database : str, optional
            path to the results_database folder for storing the new hdf file,
            by default '../results_database'
        when : str, optional
            Can be the year(s) (YYYY-YYYY) or the month and year (Month, YYYY) of the model run
        where : str, optional
            Name of the institution that run the simulation
        """

        hdf_file = build_hdf_filename(
            'openmc', self.get_openmc_version, xs_library)
        path = Path(path_to_database)
        # merge path to hdf file
        path_to_file = path / hdf_file

        # extract tally in dataframe format from statepoint file
        tally_df = self.get_tally_dataframe(
            tally_name, normalize_over=normalize_over)

        # rework the dataframe dropping useless columns
        for c in _del_columns:
            if c in tally_df.columns:
                tally_df = tally_df.drop(columns=c)

        # add xaxis columns if required
        if xaxis_list is not None:
            tally_df.insert(loc=0, column=xaxis_name, value=xaxis_list)

        # write the tally in the hdf file
        tally_df.to_hdf(path_to_file, tally_name, mode='a',
                        format='table', data_columns=True, index=False)
        code_version = 'openmc-' + '.'.join(map(str, self.get_openmc_version))

        # write attributes to the hdf file
        with h5py.File(path_to_file, 'a') as f:
            f[tally_name + '/table'].attrs['x_axis'] = xaxis_name
            f.attrs['code_version'] = code_version
            f.attrs['xs_library'] = xs_library
            f.attrs['batches'] = self.get_batches
            f.attrs['particles_per_batch'] = self.get_particles_per_batch
            f.attrs['when'] = str(when)
            f.attrs['where'] = where
