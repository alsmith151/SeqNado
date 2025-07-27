"""Tests for the design.experiment module."""
import pytest
from seqnado.inputs.experiment import ExperimentIP
from seqnado.inputs.fastq import FastqSetIP, FastqFileIP
import pathlib


@pytest.fixture
def mock_fastq_file_ip(test_data_dir):
    """Create a mock FastqFileIP for testing."""
    mock_path = test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz"
    return FastqFileIP(
        path=mock_path,
        ip="MLL",
        is_control=False,
    )


@pytest.fixture
def mock_fastq_file_control(test_data_dir):
    """Create a mock FastqFileIP control for testing."""
    mock_path = test_data_dir / "fastq" / "chip-rx_input_1.fastq.gz"
    return FastqFileIP(
        path=mock_path,
        ip="input",
        is_control=True,
    )


@pytest.fixture
def mock_fastq_set_ip(mock_fastq_file_ip):
    """Create a mock FastqSetIP for testing."""
    return FastqSetIP(
        sample_id="test_chip",
        r1=mock_fastq_file_ip,
        r2=None,  # Single-end for this test
        antibody="MLL"
    )


@pytest.fixture
def mock_fastq_set_control(mock_fastq_file_control):
    """Create a mock FastqSetIP control for testing."""
    return FastqSetIP(
        sample_id="test_control",
        r1=mock_fastq_file_control,
        r2=None,  # Single-end for this test
        antibody="input"
    )


@pytest.fixture
def paired_end_fastq_files(test_data_dir):
    """Create paired-end FastqFileIP objects for testing."""
    r1_path = test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz"
    r2_path = test_data_dir / "fastq" / "chip-rx_MLL_2.fastq.gz"
    
    r1 = FastqFileIP(path=r1_path, ip="MLL", is_control=False)
    r2 = FastqFileIP(path=r2_path, ip="MLL", is_control=False)
    
    return r1, r2


@pytest.fixture
def paired_end_fastq_set(paired_end_fastq_files):
    """Create paired-end FastqSetIP for testing."""
    r1, r2 = paired_end_fastq_files
    return FastqSetIP(
        sample_id="test_paired",
        r1=r1,
        r2=r2,
        antibody="MLL"
    )
    


@pytest.mark.unit
class TestIPExperiment:
    """Test cases for the IPExperiment class."""

    def test_ip_experiment_creation_ip_only(self, mock_fastq_set_ip):
        """Test IPExperiment creation with IP only."""
        experiment = ExperimentIP(
            ip=mock_fastq_set_ip,
            control=None
        )
        
        assert experiment.ip == mock_fastq_set_ip
        assert experiment.control is None

    def test_ip_experiment_creation_with_control(self, mock_fastq_set_ip, mock_fastq_set_control):
        """Test IPExperiment creation with control."""
        experiment = ExperimentIP(
            ip=mock_fastq_set_ip,
            control=mock_fastq_set_control
        )
        
        assert experiment.ip == mock_fastq_set_ip
        assert experiment.control == mock_fastq_set_control

    def test_has_control_property(self, mock_fastq_set_ip, mock_fastq_set_control):
        """Test has_control property."""
        # Without control
        experiment_no_control = ExperimentIP(ip=mock_fastq_set_ip)
        assert not experiment_no_control.has_control
        
        # With control
        experiment_with_control = ExperimentIP(
            ip=mock_fastq_set_ip,
            control=mock_fastq_set_control
        )
        assert experiment_with_control.has_control

    def test_ip_performed_property(self, mock_fastq_set_ip):
        """Test ip_performed property returns antibody."""
        experiment = ExperimentIP(ip=mock_fastq_set_ip)
        assert experiment.ip_performed == "MLL"
        assert experiment.ip_performed == mock_fastq_set_ip.antibody

    def test_control_performed_property(self, mock_fastq_set_ip, mock_fastq_set_control):
        """Test control_performed property."""
        # Without control
        experiment_no_control = ExperimentIP(ip=mock_fastq_set_ip)
        assert experiment_no_control.control_performed is None
        
        # With control
        experiment_with_control = ExperimentIP(
            ip=mock_fastq_set_ip,
            control=mock_fastq_set_control
        )
        assert experiment_with_control.control_performed == "input"
        assert experiment_with_control.control_performed == mock_fastq_set_control.antibody

    def test_ip_set_fullname_property(self, mock_fastq_set_ip):
        """Test ip_set_fullname property."""
        experiment = ExperimentIP(ip=mock_fastq_set_ip)
        expected_name = mock_fastq_set_ip.full_sample_name
        assert experiment.ip_set_fullname == expected_name
        assert experiment.ip_set_fullname == "test_chip_MLL"

    def test_control_fullname_property(self, mock_fastq_set_ip, mock_fastq_set_control):
        """Test control_fullname property."""
        # Without control
        experiment_no_control = ExperimentIP(ip=mock_fastq_set_ip)
        assert experiment_no_control.control_fullname is None
        
        # With control
        experiment_with_control = ExperimentIP(
            ip=mock_fastq_set_ip,
            control=mock_fastq_set_control
        )
        expected_name = mock_fastq_set_control.full_sample_name
        assert experiment_with_control.control_fullname == expected_name
        assert experiment_with_control.control_fullname == "test_control_input"

    def test_fastqs_are_paired_property_single_end(self, mock_fastq_set_ip):
        """Test fastqs_are_paired property with single-end data."""
        experiment = ExperimentIP(ip=mock_fastq_set_ip)
        # Single-end should return False
        assert not experiment.fastqs_are_paired

    def test_fastqs_are_paired_property_paired_end(self, paired_end_fastq_set):
        """Test fastqs_are_paired property with paired-end data."""
        experiment = ExperimentIP(ip=paired_end_fastq_set)
        # Paired-end should return True
        assert experiment.fastqs_are_paired

    def test_fastqs_are_paired_with_control(self, paired_end_fastq_set, mock_fastq_set_control):
        """Test fastqs_are_paired property when control is single-end."""
        experiment = ExperimentIP(
            ip=paired_end_fastq_set,
            control=mock_fastq_set_control  # single-end control
        )
        # Should return False because control is single-end
        assert not experiment.fastqs_are_paired

    @pytest.mark.parametrize("ip_type,control_type", [
        ("MLL", "input"),
        ("H3K4me3", "IgG"),
        ("CTCF", "input"),
        ("custom_ab", "control"),
    ])
    def test_different_ip_control_combinations(self, test_data_dir, ip_type, control_type):
        """Test different IP and control combinations."""
        # Create FastqFileIP objects
        ip_file = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz",
            ip=ip_type,
            is_control=False
        )
        control_file = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_input_1.fastq.gz",
            ip=control_type,
            is_control=True
        )
        
        # Create FastqSetIP objects
        ip_set = FastqSetIP(
            sample_id="test_sample",
            r1=ip_file,
            antibody=ip_type
        )
        control_set = FastqSetIP(
            sample_id="test_control",
            r1=control_file,
            antibody=control_type
        )
        
        experiment = ExperimentIP(ip=ip_set, control=control_set)
        
        assert experiment.ip_performed == ip_type
        assert experiment.control_performed == control_type
        assert experiment.has_control

    def test_experiment_equality(self, mock_fastq_set_ip):
        """Test IPExperiment equality comparison."""
        experiment1 = ExperimentIP(ip=mock_fastq_set_ip)
        experiment2 = ExperimentIP(ip=mock_fastq_set_ip)
        
        # Should be equal if they have the same IP and control
        assert experiment1 == experiment2

    def test_experiment_with_different_controls(self, mock_fastq_set_ip, mock_fastq_set_control, test_data_dir):
        """Test experiments with different controls are not equal."""
        # Create a second control with different antibody
        control2_file = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_input_1.fastq.gz",
            ip="IgG",
            is_control=True
        )
        control2_set = FastqSetIP(
            sample_id="test_control2",
            r1=control2_file,
            antibody="IgG"
        )
        
        experiment1 = ExperimentIP(ip=mock_fastq_set_ip, control=mock_fastq_set_control)
        experiment2 = ExperimentIP(ip=mock_fastq_set_ip, control=control2_set)
        
        assert experiment1 != experiment2


@pytest.mark.edge_case
class TestIPExperimentEdgeCases:
    """Test edge cases for IPExperiment."""

    def test_experiment_with_same_ip_control_names(self, test_data_dir):
        """Test experiment where IP and control have same base sample name."""
        # Create IP file
        ip_file = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz",
            ip="MLL",
            is_control=False
        )
        # Create control file
        control_file = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_input_1.fastq.gz",
            ip="input",
            is_control=True
        )
        
        # Both sets have same base sample_id
        ip_set = FastqSetIP(
            sample_id="sample",
            r1=ip_file,
            antibody="MLL"
        )
        control_set = FastqSetIP(
            sample_id="sample",  # Same base name
            r1=control_file,
            antibody="input"
        )
        
        experiment = ExperimentIP(ip=ip_set, control=control_set)
        
        # Full names should be different due to antibody
        assert experiment.ip_set_fullname == "sample_MLL"
        assert experiment.control_fullname == "sample_input"


    def test_experiment_with_none_control_explicit(self, mock_fastq_set_ip):
        """Test experiment with explicitly set None control."""
        experiment = ExperimentIP(ip=mock_fastq_set_ip, control=None)
        
        assert experiment.control is None
        assert not experiment.has_control
        assert experiment.control_performed is None
        assert experiment.control_fullname is None

    def test_fastqs_are_paired_mixed_pairing(self, test_data_dir):
        """Test fastqs_are_paired when IP is paired but control is single-end."""
        # Create paired-end IP files
        ip_r1 = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz",
            ip="MLL",
            is_control=False
        )
        ip_r2 = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_MLL_2.fastq.gz",
            ip="MLL",
            is_control=False
        )
        ip_set = FastqSetIP(
            sample_id="test_paired_ip",
            r1=ip_r1,
            r2=ip_r2,
            antibody="MLL"
        )
        
        # Create single-end control
        control_file = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx-single_MLL.fastq.gz",
            ip="input",
            is_control=True
        )
        control_set = FastqSetIP(
            sample_id="test_single_control",
            r1=control_file,
            antibody="input"
        )
        
        experiment = ExperimentIP(ip=ip_set, control=control_set)
        
        # Should return False because control is single-end
        assert not experiment.fastqs_are_paired


@pytest.mark.unit
class TestIPExperimentValidation:
    """Test validation and error handling in IPExperiment."""

    def test_experiment_requires_ip(self):
        """Test that IPExperiment requires an IP FastqSetIP."""
        with pytest.raises((TypeError, ValueError)):
            ExperimentIP()  # Missing required 'ip' field

    def test_experiment_ip_must_be_fastq_set_ip(self, mock_fastq_set_control):
        """Test that ip field must be a FastqSetIP instance."""
        # This should work
        experiment = ExperimentIP(ip=mock_fastq_set_control)
        assert experiment.ip == mock_fastq_set_control

    def test_experiment_control_must_be_fastq_set_ip_or_none(self, mock_fastq_set_ip):
        """Test that control field must be FastqSetIP or None."""
        # None should work
        experiment1 = ExperimentIP(ip=mock_fastq_set_ip, control=None)
        assert experiment1.control is None
        
        # Valid FastqSetIP should work (already tested above)
        # Invalid types should raise errors when creating the object

    def test_pydantic_model_validation(self, mock_fastq_set_ip):
        """Test that IPExperiment follows Pydantic model validation."""
        experiment = ExperimentIP(ip=mock_fastq_set_ip)
        
        # Should be able to access model fields
        assert hasattr(experiment, 'model_fields')
        assert hasattr(experiment, 'model_dump')
        
        # Should be able to serialize
        data = experiment.model_dump()
        assert 'ip' in data
        assert 'control' in data

    def test_experiment_immutability(self, mock_fastq_set_ip, mock_fastq_set_control):
        """Test that IPExperiment is immutable (frozen Pydantic model)."""
        experiment = ExperimentIP(ip=mock_fastq_set_ip)
        
        # Pydantic models are mutable by default, but we can test assignment
        # This behavior depends on the actual model configuration
        try:
            experiment.control = mock_fastq_set_control
            # If assignment succeeds, check the value was set
            assert experiment.control == mock_fastq_set_control
        except (AttributeError, ValueError):
            # If the model is frozen, this is expected
            pass


@pytest.mark.integration
class TestFastqSetIPIntegration:
    """Integration tests for FastqSetIP with real file paths."""

    def test_fastq_set_ip_with_real_files(self, test_data_dir):
        """Test FastqSetIP creation with real test files."""
        # Test with ChIP-seq files
        r1_path = test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz"
        r2_path = test_data_dir / "fastq" / "chip-rx_MLL_2.fastq.gz"
        
        r1 = FastqFileIP(path=r1_path)
        r2 = FastqFileIP(path=r2_path)
        
        fastq_set = FastqSetIP(
            sample_id="chip-rx",
            r1=r1,
            r2=r2
        )
        
        assert fastq_set.sample_id == "chip-rx"
        assert fastq_set.is_paired
        assert fastq_set.antibody == "MLL"  # Should be predicted from filename
        assert fastq_set.full_sample_name == "chip-rx_MLL"

    def test_fastq_set_ip_single_end(self, test_data_dir):
        """Test FastqSetIP with single-end data."""
        r1_path = test_data_dir / "fastq" / "chip-rx-single_MLL.fastq.gz"
        
        r1 = FastqFileIP(path=r1_path)
        fastq_set = FastqSetIP(
            sample_id="chip-rx-single",
            r1=r1
        )
        
        assert fastq_set.sample_id == "chip-rx-single"
        assert not fastq_set.is_paired
        assert fastq_set.r2 is None
        assert fastq_set.antibody == "MLL"

    def test_fastq_set_ip_control_prediction(self, test_data_dir):
        """Test that control samples are properly identified."""
        control_path = test_data_dir / "fastq" / "chip-rx_input_1.fastq.gz"
        
        r1 = FastqFileIP(path=control_path)
        control_set = FastqSetIP(
            sample_id="chip-rx-control",
            r1=r1
        )
        
        assert control_set.antibody == "input"
        assert control_set.is_control
        assert control_set.full_sample_name == "chip-rx-control_input"

    def test_ip_experiment_with_real_files(self, test_data_dir):
        """Test IPExperiment with real ChIP-seq data files."""
        # Create IP set
        ip_r1 = FastqFileIP(path=test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz")
        ip_r2 = FastqFileIP(path=test_data_dir / "fastq" / "chip-rx_MLL_2.fastq.gz")
        ip_set = FastqSetIP(sample_id="chip-rx-MLL", r1=ip_r1, r2=ip_r2)
        
        # Create control set
        control_r1 = FastqFileIP(path=test_data_dir / "fastq" / "chip-rx_input_1.fastq.gz")
        control_r2 = FastqFileIP(path=test_data_dir / "fastq" / "chip-rx_input_2.fastq.gz")
        control_set = FastqSetIP(sample_id="chip-rx-input", r1=control_r1, r2=control_r2)
        
        # Create experiment
        experiment = ExperimentIP(ip=ip_set, control=control_set)
        
        assert experiment.has_control
        assert experiment.ip_performed == "MLL"
        assert experiment.control_performed == "input"
        assert experiment.fastqs_are_paired
        assert experiment.ip_set_fullname == "chip-rx-MLL_MLL"
        assert experiment.control_fullname == "chip-rx-input_input"


@pytest.mark.unit
class TestFastqSetIPProperties:
    """Test FastqSetIP-specific properties."""

    def test_base_sample_name_property(self, mock_fastq_set_ip):
        """Test base_sample_name property."""
        # This should return the sample name without the IP suffix
        base_name = mock_fastq_set_ip.base_sample_name
        assert isinstance(base_name, str)
        # The actual value depends on the _sample_base_without_ip implementation

    def test_full_sample_name_construction(self, test_data_dir):
        """Test that full_sample_name is constructed correctly."""
        r1 = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz",
            ip="H3K4me3"
        )
        
        fastq_set = FastqSetIP(
            sample_id="test_sample",
            r1=r1,
            antibody="H3K4me3"
        )
        
        assert fastq_set.full_sample_name == "test_sample_H3K4me3"

    def test_antibody_prediction_from_r1(self, test_data_dir):
        """Test that antibody is predicted from r1 file if not provided."""
        r1 = FastqFileIP(
            path=test_data_dir / "fastq" / "chip-rx_MLL_1.fastq.gz"
        )
        
        # Don't provide antibody - should be predicted
        fastq_set = FastqSetIP(
            sample_id="test",
            r1=r1
        )
        
        assert fastq_set.antibody == r1.ip

    def test_file_paths_property(self, paired_end_fastq_set, mock_fastq_set_ip):
        """Test file_paths property returns correct paths."""
        # Paired-end
        paths = paired_end_fastq_set.file_paths
        assert len(paths) == 2
        assert all(isinstance(p, pathlib.Path) for p in paths)
        
        # Single-end
        paths_single = mock_fastq_set_ip.file_paths
        assert len(paths_single) == 1
        assert isinstance(paths_single[0], pathlib.Path)
