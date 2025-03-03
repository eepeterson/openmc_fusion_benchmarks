"""Module for defining and managing benchmarks"""
import openmc
from .cloud_interface import download_geometry
# from openmc_fusion_benchmarks import StatePoint
# from openmc_fusion_benchmarks import get_statepoint_path
import importlib


class Benchmark:
    def __init__(self, name: str, geometry_type: str):
        self.name = name
        self.geometry_type = geometry_type

        if geometry_type not in ['csg', 'dagmc']:
            raise ValueError(
                'Invalid geometry type can be either "csg" or "dagmc"')

    def model(self) -> openmc.Model:
        """Dynamically import and return the model object from benchmarks/{benchmark_name}/model.py"""
        try:
            module_path = f"openmc_fusion_benchmarks.benchmarks.{self.name}.benchmark_module"
            benchmark_module = importlib.import_module(module_path)
            # Retrieve the model function or class from the module
            benchmark_func = benchmark_module.model

            # Check if 'run_option' exists in the instance, and pass it if available
            if hasattr(self, "run_option"):
                return benchmark_func(run_option=self.run_option)
            else:
                return benchmark_func()

        except ModuleNotFoundError:
            raise ValueError(
                f"Model {self.model_name} not found in myrepo.models")

    # def statepoint(self) -> StatePoint:
    #     sp_path = get_statepoint_path(self.geometry_type)

    def get_step_file(self):
        download_geometry(self.name, 'step', self.run_option)

    def get_rtt_file(self):
        download_geometry(self.name, 'rtt', self.run_option)

    def get_h5m_file(self):
        download_geometry(self.name, 'h5m', self.run_option)

    # def get_cad_file(self, file_format: str = 'step'):

    #     if file_format not in ['step', 'rtt', 'h5m']:
    #         raise ValueError(
    #             'Invalid file format, can be "step", "rtt" or "h5m"')

    #     download_geometry(self.name, file_format, self.run_option)

    def get_weight_windows(self):
        # file needs to go on drive with the rest

        # download ww file:
        # return openmc.wwinp_to_wws(path/to/ww_file)
        pass

    def _run_and_store(self):
        pass


class FngStr(Benchmark):
    def __init__(self, geometry_type: str, run_option: str = 'onaxis'):
        super().__init__("fng_str", geometry_type)

        self.run_option = run_option


class FngW(Benchmark):
    def __init__(self, geometry_type: str, run_option: str = 'reaction_rates'):
        super().__init__("fng_w", geometry_type)

        self.run_option = run_option


class Oktavian(Benchmark):
    def __init__(self, geometry_type: str, run_option: str = 'Al'):
        super().__init__("oktavian", geometry_type)


class FnsDuct(Benchmark):
    def __init__(self, geometry_type: str):
        super().__init__("fns_duct", geometry_type)


class FnsCleanW(Benchmark):
    def __init__(self, geometry_type: str):
        super().__init__("fns_clean_w", geometry_type)


class BenchmarkDatabase:
    @staticmethod
    def get_benchmark(name: str, geometry_type: str, **kwargs):
        if name == "fng_str":
            return FngStr(geometry_type, **kwargs)
        elif name == "fng_w":
            return FngW(geometry_type, **kwargs)
        elif name == "oktavian":
            return Oktavian(geometry_type, **kwargs)
        elif name == "fns_duct":
            return FnsDuct(geometry_type, **kwargs)
        elif name == "fns_clean_w":
            return FnsCleanW(geometry_type, **kwargs)
        else:
            return Benchmark(name, geometry_type)
