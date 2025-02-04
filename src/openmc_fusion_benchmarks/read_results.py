import h5py
import openmc
from pathlib import Path
from typing import Iterable
import pandas as pd

_del_columns = ['cell', 'particle', 'nuclide', 'score', 'energyfunction']


def to_hdf(df: pd.DataFrame, file: str, tally_name: str, xs_library: str = None,
           xaxis_name: str = None,
           when: str = 'n/a', where: str = 'n/a', code_version: str = None,
           batches: int = None, particles_per_batch: int = None, literature: int = 'n/a'):
    """Stores a DataFrame to a given hdf5 file. Useful function to generate new 
    hdf5 files results.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of results to store in the hdf5 file
    file : str
        name of the hdf5 file of results. Can include the path to the file
    tally_name : str
        name of the tally to store
    xs_library : str, optional
        name of the nuclear data library used in the simulation, 
        by default None
    xaxis_name : str, optional
        name of the x_axis to store in order to be retrieved with the
        ResultsFromDatabase.get_tally_xaxis method, by default None
    when : str, optional
        Can be the year(s) (YYYY-YYYY) or the month and year (Month, YYYY)
        of the model run, by default 'n/a'
    where : str, optional
        ame of the institution that run the simulation/experiment,
        by default 'n/a'
    code_version : str, optional
        version of the code used (if simulation result), by default None
    batches : int, optional
        number of batches simulated (if simulation result, 
        assuming openmc particles/batches structure), by default None
    particles_per_batch : int, optional
        number of particles per batch simulated (if simulation result, 
        assuming openmc particles/batches structure), by default None
    literature : int, optional
        title/DOI/link if results associated to a publication,
        by default None
    """

    filepath = Path(file)
    # write the tally in the hdf file
    df.to_hdf(filepath, tally_name, mode='a',
              format='table', data_columns=True, index=False)

    # write attributes to the hdf file
    with h5py.File(filepath, 'a') as f:
        f[tally_name + '/table'].attrs['x_axis'] = xaxis_name
        f.attrs['when'] = str(when)
        f.attrs['where'] = where
        if code_version is not None:
            f.attrs['code_version'] = code_version
        if xs_library is not None:
            f.attrs['xs_library'] = xs_library
        if batches is not None:
            f.attrs['batches'] = batches
        if particles_per_batch is not None:
            f.attrs['particles_per_batch'] = particles_per_batch
        if literature is not None:
            f.attrs['literature_info'] = literature


def build_hdf_filename(code_name: str, code_version: Iterable, xs_library: str) -> str:
    """Builds the name for the hdf file to be stored in a results_database folder.

    Parameters
    ----------
    code_name : str
        name of the code used for the simulation
    code_version : Iterable
        list or tuple with the code version
    xs_library : str
        name of the nuclear data library used in the simulation

    Returns
    -------
    str
        name of the hdf file
    """

    filename = code_name.strip().replace(
        ' ', '').replace('.', '').replace('-', '').lower()
    filename += '-' + '-'.join(map(str, code_version)) + '_'
    filename += xs_library.strip().replace(' ', '').replace('.',
                                                            '').replace('-', '').lower()
    filename += '.h5'

    return filename


class ResultsFromDatabase:
    """This class takes in a hdf file and its path and generates a generic
    object by reading it. It is specifically designed for hdf files present
    in the "results_database" folders of the benchmark models in order to be
    able to read, postprocess and plot the results from previous simulations
    that have been stored there.
    """

    def __init__(self, file: str):
        """ResultsFromDatabase class constructor

        Parameters
        ----------
        file : str
            Name of the hdf file present in the results_database folder.
            Can include the path to the file
        """

        self.filename = file.strip().split('/')
        self.filepath = Path(file)

    def list_tallies(self):
        """Prints the names of all the tallies available in the hdf file
        """
        with h5py.File(self.filepath) as f:
            print(f.keys())

    def get_tally_dataframe(self, tally_name: str) -> pd.DataFrame:
        """Retrieves the results of a given tally in a Pandas DataFrame format.
        It relies on the openmc.Statepoint().get_tally().get_pandas_dataframe()
        method

        Parameters
        ----------
        tally_name : str
            Exact name of the tally in the hdf file 

        Returns
        -------
        pd.DataFrame
            DataFrame with tally results
        """
        with h5py.File(self.filepath) as f:
            df = pd.DataFrame(f[tally_name+'/table'][()]).drop(columns='index')
            # decode hdf5 strings to strings if necessary
            try:
                df[self.get_tally_xaxis(tally_name)] = [el.decode()
                                                        for el in df[self.get_tally_xaxis(tally_name)]]
            except:
                pass

            return df

    def get_tally_xaxis(self, tally_name: str) -> str:
        """Retrieves the string with the exact name associated to the tally
        x-axis. It is necessary because each result of each benchmark has a
        different data classification that can vary from depth in the shield,
        energy bin, distance from source, detector position etc. Retrieving this
        string allows to directly use it as string to get the column with that lists
        the data classification (x-axis for plotting purposes).

        Parameters
        ----------
        tally_name : str
            Exact name of the tally in the hdf file 

        Returns
        -------
        str
            Name used for the dataframe column with the x-axis info
        """
        with h5py.File(self.filepath) as f:
            return f[tally_name+'/table'].attrs['x_axis']

    @property
    def literature_info(self) -> str:
        """Retrieves the peer reviewed publication(s) associated with the data
        provided in the hdf file.

        Returns
        -------
        str
            DOI or link to the publication
        """
        with h5py.File(self.filepath) as f:
            try:
                return f.attrs['literature']
            except KeyError:
                return 'n/a'

    @property
    def when(self) -> str:
        """Retrieves the experiment or simulation place if provided.

        Returns
        -------
        str
            Experiment or simulation execution year
        """
        with h5py.File(self.filepath) as f:
            try:
                return f.attrs['when']
            except KeyError:
                return 'n/a'

    @property
    def where(self) -> str:
        """Retrieves the experiment or simulation place if provided.

        Returns
        -------
        str
            Experiment or simulation execution place
        """
        with h5py.File(self.filepath) as f:
            try:
                return f.attrs['where']
            except KeyError:
                return 'n/a'

    @property
    def code_version(self) -> str:
        """Retrieves the code version if the hdf file refers to a simulation's
        results and the info was provided.

        Returns
        -------
        str
            Code version
        """
        with h5py.File(self.filepath) as f:
            try:
                return f.attrs['code_version']
            except KeyError:
                return 'n/a'

    @property
    def xs_library(self) -> str:
        """Retrieves the nuclear data library used if the hdf file refers
        to a simulation's results and the info was provided.

        Returns
        -------
        str
            Nuclear data library name
        """
        with h5py.File(self.filepath) as f:
            try:
                return f.attrs['xs_library']
            except KeyError:
                return 'n/a'

    def print_code_info(self):
        """Prints all the code info available
        """
        print(f'Code version: {self.code_version}\n',
              f'XS library: {self.xs_library}\n')

    def print_all_info(self):
        """Prints all info available.
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
    results_database folder it is necessary to use ResultsFromDatabase class.
    """

    def __init__(self, file: str = 'statepoint.100.h5'):
        """ResultsFromOpenmc class constructor.

        Parameters
        ----------
        file : str, optional
            name of the statepoint.h5 file. Can include the path to the file,
            by default 'statepoint.100.h5'
        """
        self.filename = file.strip().split('/')
        self.filepath = Path(file)
        # open statepoint file with openmc
        self.statepoint = openmc.StatePoint(self.filepath)

    def list_tallies(self):
        """Prints the names of all the tallies available in the statepoint.h5
        """
        for k in self.statepoint.tallies.keys():
            print(self.statepoint.tallies[k].name)

    def get_tally_dataframe(self, tally_name: str, normalize_over: Iterable = None) -> pd.DataFrame:
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
        pd.DataFrame
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
    def get_openmc_version(self) -> tuple:
        """Retrieves openmc's version used in the simulation.

        Returns
        -------
        tuple
            openmc version
        """
        return self.statepoint.version

    @property
    def get_particles_per_batch(self) -> float:
        """Retrieves the number of particles per batch defined in the simulation.

        Returns
        -------
        float
            Number of particle per batch
        """
        return format(self.statepoint.n_particles, '.2e')

    @property
    def get_batches(self) -> int:
        """Retrieves the number of batches defined in the simulation.

        Returns
        -------
        int
            Number of batches
        """
        return self.statepoint.n_batches

    def tally_to_hdf(self, tally_name: str, normalize_over: Iterable, xs_library: str, xaxis_name: str,
                     xaxis_list: Iterable = None, path_to_database: str = '../results_database', when: str = 'n/a',
                     where: str = 'n/a', literature: int = None):
        """Stores the openmc tally in a hdf file for the results_database folder.

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

        filename = build_hdf_filename(
            'openmc', self.get_openmc_version, xs_library)
        file = path_to_database + '/' + filename

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

            if tally_df.columns.nlevels > 1:
                #Case where the mesh voxel indices wind up in a multi-level
                # column. Drop a level and then remove these new identically-named
                # columsn.
                tally_df.columns = tally_df.columns.droplevel(1)            
                duplicate_cols = tally_df.columns[tally_df.columns.duplicated()]
                tally_df.drop(columns=duplicate_cols, inplace=True)
            
            
        code_version = 'openmc-' + '.'.join(map(str, self.get_openmc_version))

        to_hdf(tally_df, file, tally_name, xs_library, xaxis_name, when, where,
               code_version, self.get_batches, self.get_particles_per_batch, literature)
