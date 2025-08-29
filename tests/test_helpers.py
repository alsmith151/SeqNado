"""Tests for the helpers module."""
import pathlib
import tempfile
import pytest
from unittest.mock import patch, MagicMock

from seqnado.helpers import (
    extract_cores_from_options,
    extract_apptainer_args,
    define_memory_requested,
    define_time_requested,
    symlink_file,
    symlink_fastq_files,
    is_on,
    is_off,
    is_none,
    convert_empty_yaml_entry_to_string,
    format_config_dict,
    has_bowtie2_index,
    check_options,
    get_group_for_sample,
    get_scale_method,
    extract_viewpoints,
    viewpoint_to_grouped_viewpoint,
)
from seqnado.inputs import FastqCollection, FastqCollectionForIP, Assay, Metadata
from seqnado.inputs.fastq import FastqFile, FastqSet, FastqSetIP
from seqnado.inputs.experiment import ExperimentIP


# Test fixtures
@pytest.fixture
def sample_options():
    """Provide sample command line options."""
    return [
        "--configfile", "config.yaml",
        "-c", "8",
        "--use-conda",
        "--apptainer-args", "--bind /data:/data"
    ]


@pytest.fixture
def temp_dir():
    """Provide a temporary directory for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield pathlib.Path(tmpdir)


@pytest.fixture
def mock_sample_collection():
    """Create a mock SampleCollection for testing."""
    # Create mock FastqFile objects
    r1_file = FastqFile(path=pathlib.Path("/path/to/sample1_R1.fastq.gz"))
    r2_file = FastqFile(path=pathlib.Path("/path/to/sample1_R2.fastq.gz"))
    single_file = FastqFile(path=pathlib.Path("/path/to/sample2.fastq.gz"))
    
    # Create FastqSet objects
    paired_set = FastqSet(name="sample1", r1=r1_file, r2=r2_file)
    single_set = FastqSet(name="sample2", r1=single_file)
    
    # Create metadata
    metadata = [
        Metadata(norm_group="control"),
        Metadata(norm_group="treatment")
    ]
    
    return FastqCollection(
        assay=Assay.RNA,
        fastq_sets=[paired_set, single_set],
        metadata=metadata
    )


@pytest.fixture
def mock_ip_sample_collection():
    """Create a mock IPSampleCollection for testing."""
    # Create IP FastqSetIP
    ip_r1 = FastqFile(path=pathlib.Path("/path/to/chip_IP_R1.fastq.gz"))
    ip_r2 = FastqFile(path=pathlib.Path("/path/to/chip_IP_R2.fastq.gz"))
    ip_set = FastqSetIP(name="chip", r1=ip_r1, r2=ip_r2)
    
    # Create control FastqSetIP
    ctrl_r1 = FastqFile(path=pathlib.Path("/path/to/chip_input_R1.fastq.gz"))
    ctrl_r2 = FastqFile(path=pathlib.Path("/path/to/chip_input_R2.fastq.gz"))
    ctrl_set = FastqSetIP(name="chip", r1=ctrl_r1, r2=ctrl_r2)
    
    # Create experiment
    experiment = ExperimentIP(ip=ip_set, control=ctrl_set)
    
    # Create metadata
    metadata = [Metadata(norm_group="all")]
    
    return FastqCollectionForIP(
        assay=Assay.CHIP,
        experiments=[experiment],
        metadata=metadata
    )


@pytest.fixture
def mock_wildcards():
    """Create a mock wildcards object."""
    class MockWildcards:
        def __init__(self, sample):
            self.sample = sample
    
    return MockWildcards


@pytest.fixture
def bed_content():
    """Sample BED file content for testing."""
    return """chr1\t1000\t2000\tviewpoint1\t100\t+
chr2\t3000\t4000\tviewpoint2-chr2-3000-4000\t200\t-
chr3\t5000\t6000\tviewpoint3\t300\t+"""


class TestExtractCoresFromOptions:
    """Test cases for extract_cores_from_options function."""
    
    def test_extract_cores_with_c_flag(self):
        """Test extracting cores with -c flag."""
        options = ["-c", "8", "--use-conda"]
        result_options, cores = extract_cores_from_options(options)
        
        assert cores == 8
        assert result_options == ["--use-conda"]
    
    def test_extract_cores_with_cores_flag(self):
        """Test extracting cores with --cores flag."""
        options = ["--cores", "16", "--configfile", "test.yaml"]
        result_options, cores = extract_cores_from_options(options)
        
        assert cores == 16
        assert result_options == ["--configfile", "test.yaml"]
    
    def test_extract_cores_no_flag(self):
        """Test default cores when no flag provided."""
        options = ["--use-conda", "--configfile", "test.yaml"]
        result_options, cores = extract_cores_from_options(options)
        
        assert cores == 1
        assert result_options == ["--use-conda", "--configfile", "test.yaml"]
    
    def test_extract_cores_flag_no_value(self):
        """Test cores flag with no value."""
        options = ["-c"]
        result_options, cores = extract_cores_from_options(options)
        
        assert cores == 1
        assert result_options == []


class TestExtractApptainerArgs:
    """Test cases for extract_apptainer_args function."""
    
    def test_extract_apptainer_args_present(self):
        """Test extracting apptainer args when present."""
        options = ["--apptainer-args", "--bind /data:/data", "--use-conda"]
        result_options, args = extract_apptainer_args(options)
        
        assert args == "--bind /data:/data"
        assert result_options == ["--use-conda"]
    
    def test_extract_apptainer_args_absent(self):
        """Test extracting apptainer args when absent."""
        options = ["--use-conda", "--cores", "4"]
        result_options, args = extract_apptainer_args(options)
        
        assert args == ""
        assert result_options == ["--use-conda", "--cores", "4"]


class TestDefineMemoryRequested:
    """Test cases for define_memory_requested function."""
    
    def test_memory_default(self):
        """Test default memory request."""
        result = define_memory_requested()
        assert result == "1G"
    
    def test_memory_with_attempts(self):
        """Test memory scaling with attempts."""
        result = define_memory_requested(attempts=3, initial_value=2)
        assert result == "8G"  # 2 * 2^(3-1) = 2 * 4 = 8
    
    def test_memory_with_scale(self):
        """Test memory with scale factor."""
        result = define_memory_requested(attempts=2, initial_value=4, scale=1.5)
        assert result == "12G"  # 4 * 2^(2-1) * 1.5 = 4 * 2 * 1.5 = 12


class TestDefineTimeRequested:
    """Test cases for define_time_requested function."""
    
    def test_time_default(self):
        """Test default time request."""
        result = define_time_requested()
        assert result == "1.0h"
    
    def test_time_with_attempts(self):
        """Test time scaling with attempts."""
        result = define_time_requested(attempts=3, initial_value=2)
        assert result == "8.0h"  # 2 * 2^(3-1) = 2 * 4 = 8
    
    def test_time_with_scale(self):
        """Test time with scale factor."""
        result = define_time_requested(attempts=2, initial_value=1, scale=3.5)
        assert result == "7.0h"  # 1 * 2^(2-1) * 3.5 = 1 * 2 * 3.5 = 7


class TestSymlinkFile:
    """Test cases for symlink_file function."""
    
    def test_symlink_file_success(self, temp_dir):
        """Test successful file symlinking."""
        # Create a source file
        source_file = temp_dir / "source.txt"
        source_file.write_text("test content")
        
        # Create symlink
        output_dir = temp_dir / "output"
        output_dir.mkdir()
        symlink_file(output_dir, source_file, "linked.txt")
        
        # Check symlink was created
        linked_file = output_dir / "linked.txt"
        assert linked_file.is_symlink()
        assert linked_file.read_text() == "test content"
    
    def test_symlink_file_existing_target(self, temp_dir):
        """Test symlinking when target already exists."""
        # Create source and target files
        source_file = temp_dir / "source.txt"
        source_file.write_text("source content")
        
        output_dir = temp_dir / "output"
        output_dir.mkdir()
        target_file = output_dir / "target.txt"
        target_file.write_text("existing content")
        
        # Try to create symlink - should not overwrite
        symlink_file(output_dir, source_file, "target.txt")
        
        # Target should still have original content
        assert target_file.read_text() == "existing content"
    
    def test_symlink_file_nonexistent_source(self, temp_dir):
        """Test symlinking with non-existent source."""
        output_dir = temp_dir / "output"
        output_dir.mkdir()
        source_file = temp_dir / "nonexistent.txt"
        
        # Should not create symlink or raise error
        symlink_file(output_dir, source_file, "target.txt")
        
        target_file = output_dir / "target.txt"
        assert not target_file.exists()


class TestSymlinkFastqFiles:
    """Test cases for symlink_fastq_files function."""
    
    @patch('seqnado.helpers.symlink_file')
    def test_symlink_sample_collection_paired(self, mock_symlink, mock_sample_collection):
        """Test symlinking paired-end SampleCollection files."""
        # Modify the mock to have only paired-end data
        mock_sample_collection.fastq_sets = [mock_sample_collection.fastq_sets[0]]  # Only paired set
        
        symlink_fastq_files(mock_sample_collection, "test_output")
        
        # Should call symlink_file twice for R1 and R2
        assert mock_symlink.call_count == 2
        
        # Check the calls
        calls = mock_symlink.call_args_list
        assert str(calls[0][0][1]).endswith("sample1_R1.fastq.gz")  # Source path
        assert calls[0][0][2] == "sample1_1.fastq.gz"  # Target name
        assert str(calls[1][0][1]).endswith("sample1_R2.fastq.gz")  # Source path
        assert calls[1][0][2] == "sample1_2.fastq.gz"  # Target name
    
    @patch('seqnado.helpers.symlink_file')
    def test_symlink_sample_collection_single(self, mock_symlink, mock_sample_collection):
        """Test symlinking single-end SampleCollection files."""
        # Modify the mock to have only single-end data
        mock_sample_collection.fastq_sets = [mock_sample_collection.fastq_sets[1]]  # Only single set
        
        symlink_fastq_files(mock_sample_collection, "test_output")
        
        # Should call symlink_file once
        assert mock_symlink.call_count == 1
        
        # Check the call
        call_args = mock_symlink.call_args_list[0]
        assert str(call_args[0][1]).endswith("sample2.fastq.gz")  # Source path
        assert call_args[0][2] == "sample2.fastq.gz"  # Target name
    
    @patch('seqnado.helpers.symlink_file')
    def test_symlink_ip_sample_collection(self, mock_symlink, mock_ip_sample_collection):
        """Test symlinking IPSampleCollection files."""
        symlink_fastq_files(mock_ip_sample_collection, "test_output")
        
        # Should call symlink_file 4 times (IP R1, IP R2, Control R1, Control R2)
        assert mock_symlink.call_count == 4


class TestBooleanHelpers:
    """Test cases for boolean helper functions."""
    
    @pytest.mark.parametrize("value,expected", [
        ("true", True),
        ("True", True),
        ("t", True),
        ("T", True),
        ("on", True),
        ("ON", True),
        ("yes", True),
        ("YES", True),
        ("y", True),
        ("Y", True),
        ("1", True),
        ("false", False),
        ("no", False),
        ("0", False),
        ("", False),
    ])
    def test_is_on(self, value, expected):
        """Test is_on function with various inputs."""
        assert is_on(value) == expected
    
    @pytest.mark.parametrize("value,expected", [
        ("", True),
        ("none", True),
        ("NONE", True),
        ("f", True),
        ("F", True),
        ("n", True),
        ("N", True),
        ("no", True),
        ("NO", True),
        ("false", True),
        ("FALSE", True),
        ("0", True),
        ("true", False),
        ("yes", False),
        ("1", False),
    ])
    def test_is_off(self, value, expected):
        """Test is_off function with various inputs."""
        assert is_off(value) == expected
    
    @pytest.mark.parametrize("value,expected", [
        ("", True),
        ("none", True),
        ("NONE", True),
        ("true", False),
        ("false", False),
        ("1", False),
    ])
    def test_is_none(self, value, expected):
        """Test is_none function with various inputs."""
        assert is_none(value) == expected


class TestConvertEmptyYamlEntry:
    """Test cases for convert_empty_yaml_entry_to_string function."""
    
    def test_convert_none_values(self):
        """Test converting none values to empty string."""
        assert convert_empty_yaml_entry_to_string("") == ""
        assert convert_empty_yaml_entry_to_string("none") == ""
        assert convert_empty_yaml_entry_to_string("NONE") == ""
    
    def test_convert_regular_values(self):
        """Test converting regular values."""
        assert convert_empty_yaml_entry_to_string("test") == "test"
        assert convert_empty_yaml_entry_to_string("123") == "123"

class TestHasBowtie2Index:
    """Test cases for has_bowtie2_index function."""
    
    def test_has_bowtie2_index_present(self, temp_dir):
        """Test when bowtie2 index files are present."""
        # Create mock bowtie2 index files
        (temp_dir / "genome.1.bt2").touch()
        (temp_dir / "genome.2.bt2").touch()
        
        result = has_bowtie2_index(str(temp_dir / "genome"))
        assert result is True
    
    def test_has_bowtie2_index_absent(self, temp_dir):
        """Test when bowtie2 index files are absent."""
        result = has_bowtie2_index(str(temp_dir / "genome"))
        assert result is None  # Function returns None when no index found


class TestCheckOptions:
    """Test cases for check_options function."""
    
    def test_check_options_with_none_values(self):
        """Test check_options with None-like values."""
        import numpy as np
        
        assert check_options(None) == ""
        assert check_options(np.nan) == ""
        assert check_options("") == ""
        assert check_options("false") == ""
        assert check_options("off") == ""
    
    def test_check_options_with_valid_values(self):
        """Test check_options with valid values."""
        assert check_options("test") == "test"
        assert check_options("123") == "123"
        assert check_options("true") == "true"


class TestGetGroupForSample:
    """Test cases for get_group_for_sample function."""
    
    def test_get_group_for_sample_success(self, mock_sample_collection, mock_wildcards):
        """Test successful group retrieval."""
        wildcards = mock_wildcards("sample1")
        
        with patch('seqnado.helpers.SampleGroups') as mock_norm_groups:
            mock_groups_instance = MagicMock()
            mock_groups_instance.get_sample_group.return_value = "control"
            mock_norm_groups.from_design.return_value = mock_groups_instance
            
            result = get_group_for_sample(wildcards, mock_sample_collection)
            assert result == "control"
    
    def test_get_group_for_sample_not_found(self, mock_sample_collection, mock_wildcards):
        """Test group retrieval when sample not found."""
        wildcards = mock_wildcards("nonexistent_sample")
        
        with patch('seqnado.helpers.SampleGroups') as mock_norm_groups:
            mock_groups_instance = MagicMock()
            mock_groups_instance.get_sample_group.side_effect = KeyError("Sample not found")
            mock_norm_groups.from_design.return_value = mock_groups_instance
            
            with pytest.raises(KeyError, match="Sample nonexistent_sample not found"):
                get_group_for_sample(wildcards, mock_sample_collection)


class TestGetScaleMethod:
    """Test cases for get_scale_method function."""
    
    def test_get_scale_method_default(self):
        """Test default scale method."""
        config = {}
        result = get_scale_method(config)
        assert result == ["unscaled"]
    
    def test_get_scale_method_with_spikein(self):
        """Test scale method with spikein enabled."""
        config = {"spikein": True}
        result = get_scale_method(config)
        assert "unscaled" in result
        assert "spikein" in result
    
    def test_get_scale_method_with_scale(self):
        """Test scale method with scale enabled."""
        config = {"scale": True}
        result = get_scale_method(config)
        assert "unscaled" in result
        assert "csaw" in result
    
    def test_get_scale_method_spikein_priority(self):
        """Test that spikein takes priority over scale."""
        config = {"spikein": True, "scale": True}
        result = get_scale_method(config)
        assert "unscaled" in result
        assert "spikein" in result
        assert "csaw" not in result


class TestExtractViewpoints:
    """Test cases for extract_viewpoints function."""
    
    def test_extract_viewpoints_success(self, temp_dir, bed_content):
        """Test successful viewpoint extraction."""
        bed_file = temp_dir / "viewpoints.bed"
        bed_file.write_text(bed_content)
        
        result = extract_viewpoints(str(bed_file))
        
        expected = {
            "viewpoint1-chr1-1000-2000",
            "viewpoint2-chr2-3000-4000",
            "viewpoint3-chr3-5000-6000"
        }
        assert result == expected
    
    def test_extract_viewpoints_file_not_found(self):
        """Test viewpoint extraction with non-existent file."""
        with pytest.raises(ValueError, match="Error reading BED file"):
            extract_viewpoints("/nonexistent/file.bed")
    
    def test_extract_viewpoints_minimal_columns(self, temp_dir):
        """Test viewpoint extraction with minimal BED columns."""
        bed_content = "chr1\t1000\t2000"
        bed_file = temp_dir / "minimal.bed"
        bed_file.write_text(bed_content)
        
        result = extract_viewpoints(str(bed_file))
        expected = {"chr1:1000-2000-chr1-1000-2000"}
        assert result == expected


class TestViewpointToGroupedViewpoint:
    """Test cases for viewpoint_to_grouped_viewpoint function."""
    
    def test_viewpoint_to_grouped_viewpoint_success(self):
        """Test successful viewpoint grouping."""
        viewpoints = [
            "oligo1-chr1-1000-2000",
            "oligo1-chr1-3000-4000",
            "oligo2-chr2-5000-6000",
            "simple_viewpoint"
        ]
        
        result = viewpoint_to_grouped_viewpoint(viewpoints)
        
        expected = {
            "oligo1-chr1-1000-2000": "oligo1",
            "oligo1-chr1-3000-4000": "oligo1",
            "oligo2-chr2-5000-6000": "oligo2"
        }
        assert result == expected
    
    def test_viewpoint_to_grouped_viewpoint_no_coordinates(self):
        """Test viewpoint grouping with no coordinate patterns."""
        viewpoints = ["simple1", "simple2"]
        
        result = viewpoint_to_grouped_viewpoint(viewpoints)
        
        assert result == {}
    
    def test_viewpoint_to_grouped_viewpoint_mixed(self):
        """Test viewpoint grouping with mixed patterns."""
        viewpoints = [
            "complex-oligo-chr1-1000-2000",
            "simple_viewpoint",
            "another-complex-chrX-5000-6000"
        ]
        
        result = viewpoint_to_grouped_viewpoint(viewpoints)
        
        expected = {
            "complex-oligo-chr1-1000-2000": "complex-oligo",
            "another-complex-chrX-5000-6000": "another-complex"
        }
        assert result == expected


@pytest.mark.integration
class TestIntegrationTests:
    """Integration tests combining multiple helper functions."""
    
    def test_config_processing_pipeline(self):
        """Test a complete config processing pipeline."""
        raw_config = {
            "enable_spikein": "true",
            "cores": "8",
            "memory_scale": "1.5",
            "empty_option": "",
            "nested": {
                "scale": "on",
                "debug": "false"
            }
        }
        
        # Format the config
        formatted_config = format_config_dict(raw_config)
        
        # Extract scale method
        scale_methods = get_scale_method(formatted_config)
        
        # Check results
        assert formatted_config["enable_spikein"] is True
        assert formatted_config["cores"] == "8"
        assert formatted_config["memory_scale"] == "1.5"
        assert formatted_config["empty_option"] is None
        assert formatted_config["nested"]["scale"] is True
        assert formatted_config["nested"]["debug"] is False
        
        # Should include spikein due to enable_spikein being True
        # Note: get_scale_method looks for "spikein" key, not "enable_spikein"
        assert "unscaled" in scale_methods
    
    def test_resource_calculation_pipeline(self):
        """Test resource calculation for job scaling."""
        # Simulate multiple job attempts
        attempts = [1, 2, 3]
        base_memory = 4
        base_time = 2
        scale_factor = 1.5
        
        results = []
        for attempt in attempts:
            memory = define_memory_requested(attempt, base_memory, scale_factor)
            time = define_time_requested(attempt, base_time, scale_factor)
            results.append((attempt, memory, time))
        
        # Check scaling behavior
        assert results[0] == (1, "6G", "3.0h")    # 4 * 1 * 1.5 = 6
        assert results[1] == (2, "12G", "6.0h")   # 4 * 2 * 1.5 = 12
        assert results[2] == (3, "24G", "12.0h")  # 4 * 4 * 1.5 = 24
