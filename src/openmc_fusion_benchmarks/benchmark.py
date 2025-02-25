import openmc

class Benchmark:
    def __init__(self, geometry_type:str):
        self.geometry_type = geometry_type

        if geometry_type not in ['csg', 'cad']:
            raise ValueError('Invalid geometry type can be either "csg" or "cad"')

    def model(self) -> openmc.Model:
        pass
    
    def statepoint(self) -> openmc.StatePoint:
        pass
    
    def _run_and_store(self):
        pass