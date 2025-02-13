import pytest

import openmc_fusion_benchmarks as ofb
import matplotlib.pyplot as plt

def test_add_floor_ceiling():

    checklist_lin = [1,2,3,4,5]
    checklist_log = [1,2,3,4,5,6,7,8,9,10]

    fig, ax = plt.subplots()

    min_value,max_value = ofb.add_floor_ceiling(ax,checklist_lin, scale='lin', gap=0)
    assert min_value == 1
    assert max_value == 5

    min_value,max_value = ofb.add_floor_ceiling(ax,checklist_lin, scale='lin', gap=5)
    assert min_value == -4
    assert max_value == 10

    min_value,max_value = ofb.add_floor_ceiling(ax,checklist_log, scale='log', gap=0)
    assert min_value == 0
    assert max_value == 1

    min_value,max_value = ofb.add_floor_ceiling(ax,checklist_log, scale='log', gap=1)
    assert min == -1
    assert max == 2

    plt.close(fig)