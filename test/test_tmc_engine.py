import importlib.util
import pytest
from unittest.mock import Mock, patch, MagicMock
import openmc
import numpy as np
from pathlib import Path
from openmc_fusion_benchmarks.uq.tmc_engine import tmc_engine

try:
    OPENMC_AVAILABLE = importlib.util.find_spec("openmc") is not None
except ValueError:
    OPENMC_AVAILABLE = False


@pytest.fixture
def mock_openmc_model():
    """Create a mock OpenMC model."""
    model = Mock(spec=openmc.Model)
    model.run = Mock(return_value="statepoint.001.h5")
    return model


@pytest.fixture
def temp_results_dir(tmp_path, monkeypatch):
    """Set up a temporary directory for test results."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_tmc_engine_basic_structure(mock_openmc_model, temp_results_dir):
    """Test that tmc_engine accepts correct parameters."""
    with patch('openmc_fusion_benchmarks.uq.tmc_engine.get_nuclide_gnds') as mock_gnds, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_to_hdf5') as mock_perturb, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_xs_xml') as mock_perturb_xml, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.StatePoint') as mock_sp, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine._save_result') as mock_save, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.config', {'cross_sections': ''}, create=True):
        
        # Setup mocks
        mock_gnds.return_value = 'H1'
        
        # Create mock tally
        mock_tally = Mock()
        mock_tally.name = 'test_tally'
        mock_df = Mock()
        mock_df.drop = Mock(return_value=mock_df)
        mock_df.columns = ['mean', 'std. dev.']
        mock_df.values = np.array([[1.0, 0.1]])
        mock_df.shape = (1, 2)
        mock_tally.get_pandas_dataframe = Mock(return_value=mock_df)
        
        # Setup statepoint mock
        mock_statepoint = Mock()
        mock_statepoint.tallies = [1]
        mock_statepoint.get_tally = Mock(return_value=mock_tally)
        mock_sp.return_value = mock_statepoint
        
        # Create directory that perturb_xs_xml expects
        xs_dir = temp_results_dir / "H1_ENDF/B-VIII.0"
        xs_dir.mkdir(parents=True)
        xs_file = xs_dir / "H1_0_ENDF" / "B-VIII.0.h5"
        xs_file.parent.mkdir(parents=True, exist_ok=True)
        xs_file.touch()
        
        # Run tmc_engine with minimal realizations
        tmc_engine(
            model=mock_openmc_model,
            realizations=1,
            lib_name='ENDF/B-VIII.0',
            nuclide='H1',
            perturb_xs=True
        )
        
        # Verify key functions were called
        assert mock_gnds.called
        assert mock_perturb.called
        assert mock_openmc_model.run.called


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_tmc_engine_without_perturbation(mock_openmc_model, temp_results_dir):
    """Test tmc_engine with perturb_xs=False."""
    with patch('openmc_fusion_benchmarks.uq.tmc_engine.get_nuclide_gnds') as mock_gnds, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_to_hdf5') as mock_perturb, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_xs_xml') as mock_perturb_xml, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.StatePoint') as mock_sp, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine._save_result') as mock_save, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.config', {'cross_sections': ''}, create=True):
        
        mock_gnds.return_value = 'U238'
        
        # Create mock tally
        mock_tally = Mock()
        mock_tally.name = 'test_tally'
        mock_df = Mock()
        mock_df.drop = Mock(return_value=mock_df)
        mock_df.columns = ['mean', 'std. dev.']
        mock_df.values = np.array([[2.0, 0.2]])
        mock_df.shape = (1, 2)
        mock_tally.get_pandas_dataframe = Mock(return_value=mock_df)
        
        mock_statepoint = Mock()
        mock_statepoint.tallies = [1]
        mock_statepoint.get_tally = Mock(return_value=mock_tally)
        mock_sp.return_value = mock_statepoint
        
        # Create expected directory
        xs_dir = temp_results_dir / "U238_ENDF/B-VIII.0"
        xs_dir.mkdir(parents=True)
        xs_file = xs_dir / "U238_0_ENDF" / "B-VIII.0.h5"
        xs_file.parent.mkdir(parents=True, exist_ok=True)
        xs_file.touch()
        
        # Run without perturbation
        tmc_engine(
            model=mock_openmc_model,
            realizations=1,
            lib_name='ENDF/B-VIII.0',
            nuclide='U238',
            perturb_xs=False
        )
        
        # Verify perturbation was NOT called
        assert not mock_perturb.called
        # But model.run should still be called
        assert mock_openmc_model.run.called


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_tmc_engine_multiple_realizations(mock_openmc_model, temp_results_dir):
    """Test that tmc_engine runs multiple realizations."""
    with patch('openmc_fusion_benchmarks.uq.tmc_engine.get_nuclide_gnds') as mock_gnds, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_to_hdf5') as mock_perturb, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_xs_xml') as mock_perturb_xml, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.StatePoint') as mock_sp, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine._save_result') as mock_save, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.config', {'cross_sections': ''}, create=True):
        
        mock_gnds.return_value = 'H1'
        
        # Create mock tally
        mock_tally = Mock()
        mock_tally.name = 'test_tally'
        mock_df = Mock()
        mock_df.drop = Mock(return_value=mock_df)
        mock_df.columns = ['mean', 'std. dev.']
        mock_df.values = np.array([[1.5, 0.15]])
        mock_df.shape = (1, 2)
        mock_tally.get_pandas_dataframe = Mock(return_value=mock_df)
        
        mock_statepoint = Mock()
        mock_statepoint.tallies = [1]
        mock_statepoint.get_tally = Mock(return_value=mock_tally)
        mock_sp.return_value = mock_statepoint
        
        # Create directories for multiple realizations
        base_dir = temp_results_dir / "H1_ENDF/B-VIII.0"
        base_dir.mkdir(parents=True)
        for i in range(3):
            xs_file = base_dir / f"H1_{i}_ENDF" / "B-VIII.0.h5"
            xs_file.parent.mkdir(parents=True, exist_ok=True)
            xs_file.touch()
        
        # Run with 3 realizations
        tmc_engine(
            model=mock_openmc_model,
            realizations=3,
            lib_name='ENDF/B-VIII.0',
            nuclide='H1',
            perturb_xs=True
        )
        
        # Verify model.run was called 3 times
        assert mock_openmc_model.run.call_count == 3
        # Verify _save_result was called 3 times (once per realization)
        assert mock_save.call_count == 3


@pytest.mark.skipif(not OPENMC_AVAILABLE, reason="OpenMC not installed")
def test_tmc_engine_benchmark_mode(mock_openmc_model, temp_results_dir):
    """Test tmc_engine in benchmark mode (_is_benchmark=True)."""
    with patch('openmc_fusion_benchmarks.uq.tmc_engine.get_nuclide_gnds') as mock_gnds, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_to_hdf5') as mock_perturb, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.perturb_xs_xml') as mock_perturb_xml, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.StatePoint') as mock_sp, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine._openmc_to_ofb') as mock_ofb, \
         patch('openmc_fusion_benchmarks.uq.tmc_engine.openmc.config', {'cross_sections': ''}, create=True):
        
        mock_gnds.return_value = 'H1'
        
        mock_statepoint = Mock()
        mock_sp.return_value = mock_statepoint
        
        # Create directory
        xs_dir = temp_results_dir / "H1_ENDF/B-VIII.0"
        xs_dir.mkdir(parents=True)
        xs_file = xs_dir / "H1_0_ENDF" / "B-VIII.0.h5"
        xs_file.parent.mkdir(parents=True, exist_ok=True)
        xs_file.touch()
        
        # Run in benchmark mode
        tmc_engine(
            model=mock_openmc_model,
            realizations=1,
            lib_name='ENDF/B-VIII.0',
            nuclide='H1',
            perturb_xs=False,
            _is_benchmark=True,
            _mesh='dummy_mesh',
            _spec_tallies='dummy_tallies'
        )
        
        # Verify _openmc_to_ofb was called instead of regular save
        assert mock_ofb.called
        mock_ofb.assert_called_once()
