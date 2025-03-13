import openmc_fusion_benchmarks


def test_version():
    """Simple check that to tst the  .__version__ attribute exists"""
    version = openmc_fusion_benchmarks.__version__
    assert isinstance(version, str)
