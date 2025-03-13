import pytest
import openmc
import openmc_fusion_benchmarks as ofb


def test_benchmark_init():
    benchmark = ofb.Benchmark('test')
    assert benchmark.name == 'test'

# def test_benchmark_get_model():
#     benchmark = ofb.Benchmark('test')
#     model = benchmark.get_model('csg')
#     assert model is not None
#     assert type(model) == openmc.Model
#     assert hasattr(model, 'settings')
#     assert hasattr(model, 'geometry')
#     assert hasattr(model, 'materials')
#     assert hasattr(model, 'tallies')
#     assert hasattr(model, 'run')
#     assert callable(model.run)
