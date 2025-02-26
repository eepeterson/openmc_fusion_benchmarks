import openmc
from openmc_fusion_benchmarks import StatePoint
from openmc_fusion_benchmarks import get_statepoint_path

class Benchmark:
    def __init__(self, name, geometry_type:str):
        self.name = name
        self.geometry_type = geometry_type

        if geometry_type not in ['csg', 'cad']:
            raise ValueError('Invalid geometry type can be either "csg" or "cad"')

    def model(self) -> openmc.Model:
        pass
    
    def statepoint(self) -> StatePoint:
        sp_path = get_statepoint_path(self.geometry_type)
    
    def _run_and_store(self):
        pass

class FngStr(Benchmark):
    def __init__(self, geometry_type:str, run_option:str='onaxis'):
        super().__init__("fng_str", geometry_type)

class FngW(Benchmark):
    def __init__(self, geometry_type:str, run_option:str='reaction_rates'):
        super().__init__("fng_w", geometry_type)

class Oktavian(Benchmark):
    def __init__(self, geometry_type:str, run_option:str='Al'):
        super().__init__("oktavian", geometry_type)

class FnsDuct(Benchmark):
    def __init__(self, geometry_type:str):
        super().__init__("fns_duct", geometry_type)

class FnsCleanW(Benchmark):
    def __init__(self, geometry_type:str):
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