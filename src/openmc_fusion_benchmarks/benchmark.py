"""Module for defining and managing benchmarks"""
import openmc
from .cloud_interface import download_geometry
# from openmc_fusion_benchmarks import StatePoint
# from openmc_fusion_benchmarks import get_statepoint_path
from functools import wraps


class Benchmark:
    def __init__(self, name: str):
        self.name = name

    def get_model(self, geometry_type: str) -> openmc.Model:
        """Dynamically import and return the model object from benchmarks/{benchmark_name}/model.py"""

        if geometry_type not in ['csg', 'cad']:
            raise ValueError(
                'Invalid geometry type can be either "csg" or "cad"')

        try:
            module_path = f"openmc_fusion_benchmarks.benchmarks.{self.name}.benchmark_module"
            benchmark_module = __import__(module_path, fromlist=['model'])
            # Retrieve the model function or class from the module
            benchmark_func = benchmark_module.model

            model = (
                benchmark_func(geometry_type=geometry_type,
                               run_option=self.run_option)
                if hasattr(self, "run_option")
                else benchmark_func(geometry_type=geometry_type)
            )

            # Wrap `run()` only if geometry_type == 'cad'
            if geometry_type == "cad" and hasattr(model, "run") and callable(model.run):
                model.run = _wrap_run(self.download_h5m_file, model.run)

            return model

        except ModuleNotFoundError:
            raise ValueError(
                f"Model {self.model_name} not found in myrepo.models")

    # def statepoint(self) -> StatePoint:
    #     sp_path = get_statepoint_path(self.geometry_type)

    def download_step_file(self, cwd: str = None):
        download_geometry(self.name, 'step', self.run_option, cwd)

    def download_rtt_file(self, cwd: str = None):
        download_geometry(self.name, 'rtt', self.run_option, cwd)

    def download_h5m_file(self, cwd: str = None):
        download_geometry(self.name, 'h5m', self.run_option, cwd)

    def download_weight_windows(self, cwd: str = None):
        # file needs to go on drive with the rest

        # download ww file:
        # return openmc.wwinp_to_wws(path/to/ww_file)
        pass

    def _run_and_store(self):
        pass


def _wrap_run(download_files_func, original_run):
    """Standalone function to wrap `run()` and ensure files are downloaded first."""
    @wraps(original_run)
    def wrapped_run(*args, **kwargs):
        cwd = kwargs.get("cwd", ".")  # Extract cwd argument (default: ".")
        download_files_func(cwd)  # Ensure files are downloaded
        return original_run(*args, **kwargs)  # Call the original method
    return wrapped_run


class FngStr(Benchmark):
    def __init__(self, run_option: str = 'onaxis'):
        super().__init__("fng_str")

        self.run_option = run_option


class FngW(Benchmark):
    def __init__(self, run_option: str = 'reaction_rates'):
        super().__init__("fng_w")

        self.run_option = run_option


class Oktavian(Benchmark):
    def __init__(self, run_option: str = 'Al'):
        super().__init__("oktavian")


class FnsDuct(Benchmark):
    def __init__(self):
        super().__init__("fns_duct")


class FnsCleanW(Benchmark):
    def __init__(self):
        super().__init__("fns_clean_w")


class BenchmarkDatabase:
    @staticmethod
    def get_benchmark(name: str, **kwargs):
        if name == "fng_str":
            return FngStr(**kwargs)
        elif name == "fng_w":
            return FngW(**kwargs)
        elif name == "oktavian":
            return Oktavian(**kwargs)
        elif name == "fns_duct":
            return FnsDuct(**kwargs)
        elif name == "fns_clean_w":
            return FnsCleanW(**kwargs)
        else:
            return Benchmark(name)
