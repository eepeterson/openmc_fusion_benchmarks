import pytest
from openmc_fusion_benchmarks import build_hdf_filename

def test_dummy():
    assert 1 == 1


def test_build_hdf_filename():

    sample1 = build_hdf_filename('test', (0,0,0), 'test')
    sample2 = build_hdf_filename('openmc', (0,15,0), 'fendl-3.2b')
    sample3 = build_hdf_filename('mcnp', (4,6,0), 'endf-b-8.1')


    assert sample1 == 'test-0-0-0_test.h5'
    assert sample2 == 'openmc-0-15-0_fendl-3.2b.h5'
    assert sample3 == 'mcnp-4-6-0_endf-b-8.1.h5'