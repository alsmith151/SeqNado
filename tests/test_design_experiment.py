"""Tests for the design.experiment module."""
import pytest
from seqnado.design.experiment import IPExperiment
import pathlib


@pytest.fixture
def mock_fastq_file_ip():
    """Create a mock FastqFileIP for testing."""
    # Using a mock path since we don't need real files for these tests
    mock_path = pathlib.Path("/test/path/sample_chip_1.fastq.gz")
    
    # Create a minimal mock that has the required attributes
    class MockFastqFileIP:
        def __init__(self, path):
            self.path = path
            self.sample_base_without_ip = "sample"
            self.is_control = False
            self.read_number = 1
    
    return MockFastqFileIP(mock_path)


@pytest.fixture
def mock_fastq_set_ip(mock_fastq_file_ip):
    """Create a mock FastqSetIP for testing."""
    class MockFastqSetIP:
        def __init__(self, name, r1, r2=None, ip_performed="chip"):
            self.name = name
            self.r1 = r1
            self.r2 = r2
            self.ip_performed = ip_performed
    
    return MockFastqSetIP("sample", mock_fastq_file_ip)


@pytest.mark.unit
class TestIPExperiment:
    """Test cases for the IPExperiment class."""

    def test_ip_experiment_creation_ip_only(self, mock_fastq_set_ip):
        """Test IPExperiment creation with IP only."""
        experiment = IPExperiment(
            ip=mock_fastq_set_ip,
            control=None
        )
        
        assert experiment.ip == mock_fastq_set_ip
        assert experiment.control is None

    def test_ip_experiment_creation_with_control(self, mock_fastq_set_ip):
        """Test IPExperiment creation with control."""
        control_set = type(mock_fastq_set_ip)(
            "sample", 
            mock_fastq_set_ip.r1, 
            ip_performed="input"
        )
        
        experiment = IPExperiment(
            ip=mock_fastq_set_ip,
            control=control_set
        )
        
        assert experiment.ip == mock_fastq_set_ip
        assert experiment.control == control_set

    def test_has_control_property(self, mock_fastq_set_ip):
        """Test has_control property."""
        # Without control
        experiment_no_control = IPExperiment(ip=mock_fastq_set_ip)
        assert not experiment_no_control.has_control
        
        # With control
        control_set = type(mock_fastq_set_ip)(
            "sample", 
            mock_fastq_set_ip.r1, 
            ip_performed="input"
        )
        experiment_with_control = IPExperiment(
            ip=mock_fastq_set_ip,
            control=control_set
        )
        assert experiment_with_control.has_control

    def test_ip_performed_property(self, mock_fastq_set_ip):
        """Test ip_performed property."""
        experiment = IPExperiment(ip=mock_fastq_set_ip)
        assert experiment.ip_performed == "chip"

    def test_control_performed_property(self, mock_fastq_set_ip):
        """Test control_performed property."""
        # Without control
        experiment_no_control = IPExperiment(ip=mock_fastq_set_ip)
        assert experiment_no_control.control_performed is None
        
        # With control
        control_set = type(mock_fastq_set_ip)(
            "sample", 
            mock_fastq_set_ip.r1, 
            ip_performed="input"
        )
        experiment_with_control = IPExperiment(
            ip=mock_fastq_set_ip,
            control=control_set
        )
        assert experiment_with_control.control_performed == "input"

    def test_ip_set_fullname_property(self, mock_fastq_set_ip):
        """Test ip_set_fullname property."""
        experiment = IPExperiment(ip=mock_fastq_set_ip)
        expected_name = f"{mock_fastq_set_ip.name}_{mock_fastq_set_ip.ip_performed}"
        assert experiment.ip_set_fullname == expected_name

    def test_control_fullname_property(self, mock_fastq_set_ip):
        """Test control_fullname property."""
        # Without control
        experiment_no_control = IPExperiment(ip=mock_fastq_set_ip)
        assert experiment_no_control.control_fullname is None
        
        # With control
        control_set = type(mock_fastq_set_ip)(
            "sample", 
            mock_fastq_set_ip.r1, 
            ip_performed="input"
        )
        experiment_with_control = IPExperiment(
            ip=mock_fastq_set_ip,
            control=control_set
        )
        expected_name = f"{control_set.name}_{control_set.ip_performed}"
        assert experiment_with_control.control_fullname == expected_name

    @pytest.mark.parametrize("ip_type,control_type", [
        ("chip", "input"),
        ("cat", "igg"),
        ("h3k4me3", "input"),
        ("custom_ab", "control"),
    ])
    def test_different_ip_control_combinations(self, mock_fastq_file_ip, ip_type, control_type):
        """Test different IP and control combinations."""
        # Create mock sets with different IP types
        ip_set = type('MockFastqSetIP', (), {
            'name': 'sample',
            'r1': mock_fastq_file_ip,
            'r2': None,
            'ip_performed': ip_type
        })()
        
        control_set = type('MockFastqSetIP', (), {
            'name': 'sample',
            'r1': mock_fastq_file_ip,
            'r2': None,
            'ip_performed': control_type
        })()
        
        experiment = IPExperiment(ip=ip_set, control=control_set)
        
        assert experiment.ip_performed == ip_type
        assert experiment.control_performed == control_type
        assert experiment.has_control

    def test_experiment_equality(self, mock_fastq_set_ip):
        """Test IPExperiment equality comparison."""
        experiment1 = IPExperiment(ip=mock_fastq_set_ip)
        experiment2 = IPExperiment(ip=mock_fastq_set_ip)
        
        # Should be equal if they have the same IP and control
        assert experiment1 == experiment2

    def test_experiment_with_different_controls(self, mock_fastq_set_ip):
        """Test experiments with different controls are not equal."""
        control1 = type(mock_fastq_set_ip)("sample", mock_fastq_set_ip.r1, ip_performed="input")
        control2 = type(mock_fastq_set_ip)("sample", mock_fastq_set_ip.r1, ip_performed="igg")
        
        experiment1 = IPExperiment(ip=mock_fastq_set_ip, control=control1)
        experiment2 = IPExperiment(ip=mock_fastq_set_ip, control=control2)
        
        assert experiment1 != experiment2


@pytest.mark.edge_case
class TestIPExperimentEdgeCases:
    """Test edge cases for IPExperiment."""

    def test_experiment_with_same_ip_control_names(self, mock_fastq_file_ip):
        """Test experiment where IP and control have same base name."""
        ip_set = type('MockFastqSetIP', (), {
            'name': 'sample',
            'r1': mock_fastq_file_ip,
            'ip_performed': 'chip'
        })()
        
        control_set = type('MockFastqSetIP', (), {
            'name': 'sample',  # Same name
            'r1': mock_fastq_file_ip,
            'ip_performed': 'input'
        })()
        
        experiment = IPExperiment(ip=ip_set, control=control_set)
        
        # Full names should be different due to IP type
        assert experiment.ip_set_fullname == "sample_chip"
        assert experiment.control_fullname == "sample_input"

    def test_experiment_string_representation(self, mock_fastq_set_ip):
        """Test string representation of IPExperiment."""
        experiment = IPExperiment(ip=mock_fastq_set_ip)
        str_repr = str(experiment)
        
        # Should contain key information
        assert "IPExperiment" in str_repr
        assert mock_fastq_set_ip.name in str_repr
        assert mock_fastq_set_ip.ip_performed in str_repr
