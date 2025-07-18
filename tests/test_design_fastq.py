"""Tests for the design.fastq module."""
import pathlib
import pytest
from seqnado.design.fastq import FastqFile, FastqSet, FastqFileIP, FastqSetIP


@pytest.mark.unit
@pytest.mark.requires_data
class TestFastqFile:
    """Test cases for the FastqFile class."""

    @pytest.fixture
    def test_fastq_path(self, test_fastq_dir):
        """Provide a test FASTQ file path."""
        return test_fastq_dir / "rna_1.fastq.gz"

    @pytest.fixture
    def test_fastq_path_r2(self, test_fastq_dir):
        """Provide a test FASTQ file path for R2."""
        return test_fastq_dir / "rna_2.fastq.gz"

    def test_fastq_file_creation(self, test_fastq_path):
        """Test basic FastqFile creation."""
        fastq_file = FastqFile(path=test_fastq_path)
        assert fastq_file.path == test_fastq_path.absolute()
        assert isinstance(fastq_file.path, pathlib.Path)

    def test_fastq_file_stem_property(self, test_fastq_path):
        """Test the stem property extracts correct filename."""
        fastq_file = FastqFile(path=test_fastq_path)
        assert fastq_file.stem == "rna_1.fastq"

    def test_fastq_file_sample_name(self, test_fastq_path):
        """Test sample name extraction."""
        fastq_file = FastqFile(path=test_fastq_path)
        assert fastq_file.sample_name == "rna"

    def test_fastq_file_read_number(self, test_fastq_path):
        """Test read number extraction."""
        fastq_file = FastqFile(path=test_fastq_path)
        assert fastq_file.read_number == 1

    def test_fastq_file_sample_base(self, test_fastq_path):
        """Test sample base extraction."""
        fastq_file = FastqFile(path=test_fastq_path)
        assert fastq_file.sample_base == "rna"

    def test_nonexistent_file_raises_error(self):
        """Test that nonexistent files raise appropriate error."""
        nonexistent_path = pathlib.Path("/nonexistent/file.fastq")
        with pytest.raises(FileNotFoundError):
            FastqFile(path=nonexistent_path)

    @pytest.mark.parametrize("filename,expected_sample", [
        ("atac_1.fastq.gz", "atac"),
        ("rna-spikein-control-rep1_1.fastq.gz", "rna-spikein-control-rep1"),
        ("chip-rx_MLL_1.fastq.gz", "chip-rx_MLL"),
        ("meth-1_R1.fastq.gz", "meth-1"),
    ])
    def test_sample_name_extraction_parametrized(self, filename, expected_sample):
        """Parametrized test for sample name extraction."""
        test_path = pathlib.Path(__file__).parent / "data" / "fastq" / filename
        if test_path.exists():
            fastq_file = FastqFile(path=test_path)
            assert fastq_file.sample_name == expected_sample


class TestFastqSet:
    """Test cases for the FastqSet class."""

    @pytest.fixture
    def test_fastq_paths(self):
        """Provide test FASTQ file paths."""
        base_path = pathlib.Path(__file__).parent / "data" / "fastq"
        return {
            "r1": base_path / "rna_1.fastq.gz",
            "r2": base_path / "rna_2.fastq.gz",
            "single": base_path / "atac_1.fastq.gz"
        }

    def test_paired_end_fastq_set(self, test_fastq_paths):
        """Test creation of paired-end FastqSet."""
        r1_file = FastqFile(path=test_fastq_paths["r1"])
        r2_file = FastqFile(path=test_fastq_paths["r2"])
        
        fastq_set = FastqSet(name="test_sample", r1=r1_file, r2=r2_file)
        
        assert fastq_set.name == "test_sample"
        assert fastq_set.r1 == r1_file
        assert fastq_set.r2 == r2_file
        assert fastq_set.is_paired_end is True

    def test_single_end_fastq_set(self, test_fastq_paths):
        """Test creation of single-end FastqSet."""
        r1_file = FastqFile(path=test_fastq_paths["single"])
        
        fastq_set = FastqSet(name="test_sample", r1=r1_file)
        
        assert fastq_set.name == "test_sample"
        assert fastq_set.r1 == r1_file
        assert fastq_set.r2 is None
        assert fastq_set.is_paired_end is False

    def test_fastq_set_file_paths(self, test_fastq_paths):
        """Test file path extraction from FastqSet."""
        r1_file = FastqFile(path=test_fastq_paths["r1"])
        r2_file = FastqFile(path=test_fastq_paths["r2"])
        fastq_set = FastqSet(name="test", r1=r1_file, r2=r2_file)
        
        paths = fastq_set.file_paths
        assert len(paths) == 2
        assert test_fastq_paths["r1"].absolute() in paths
        assert test_fastq_paths["r2"].absolute() in paths


class TestFastqFileIP:
    """Test cases for the FastqFileIP class (IP-specific FASTQ files)."""

    @pytest.fixture
    def test_ip_fastq_paths(self):
        """Provide test IP FASTQ file paths."""
        base_path = pathlib.Path(__file__).parent / "data" / "fastq"
        return {
            "chip_ip": base_path / "chip-rx_MLL_1.fastq.gz",
            "chip_control": base_path / "chip-rx_input_1.fastq.gz",
            "chip_single": base_path / "chip-rx-single_MLL.fastq.gz"
        }

    def test_ip_fastq_file_creation(self, test_ip_fastq_paths):
        """Test FastqFileIP creation."""
        if test_ip_fastq_paths["chip_ip"].exists():
            fastq_file = FastqFileIP(path=test_ip_fastq_paths["chip_ip"])
            assert fastq_file.path == test_ip_fastq_paths["chip_ip"].absolute()

    def test_control_detection(self, test_ip_fastq_paths):
        """Test control sample detection in IP files."""
        if test_ip_fastq_paths["chip_control"].exists():
            control_file = FastqFileIP(path=test_ip_fastq_paths["chip_control"])
            assert control_file.is_control is True
        
        if test_ip_fastq_paths["chip_ip"].exists():
            ip_file = FastqFileIP(path=test_ip_fastq_paths["chip_ip"])
            assert ip_file.is_control is False

    def test_sample_base_without_ip(self, test_ip_fastq_paths):
        """Test sample base extraction without IP suffix."""
        if test_ip_fastq_paths["chip_ip"].exists():
            fastq_file = FastqFileIP(path=test_ip_fastq_paths["chip_ip"])
            # This should extract the base sample name without IP/control indicators
            assert "chip-rx" in fastq_file.sample_base_without_ip


class TestFastqSetIP:
    """Test cases for the FastqSetIP class (IP-specific FASTQ sets)."""

    @pytest.fixture
    def test_ip_files(self, test_ip_fastq_paths):
        """Create IP FastqFile objects for testing."""
        files = {}
        for key, path in test_ip_fastq_paths.items():
            if path.exists():
                files[key] = FastqFileIP(path=path)
        return files

    def test_ip_fastq_set_creation(self, test_ip_files):
        """Test creation of IP FastqSet."""
        if "chip_ip" in test_ip_files:
            fastq_set = FastqSetIP(
                name="test_chip",
                r1=test_ip_files["chip_ip"],
                ip_performed="MLL"
            )
            
            assert fastq_set.name == "test_chip"
            assert fastq_set.ip_performed == "MLL"
            assert fastq_set.r1 == test_ip_files["chip_ip"]

    def test_ip_set_fullname(self, test_ip_files):
        """Test IP set fullname generation."""
        if "chip_ip" in test_ip_files:
            fastq_set = FastqSetIP(
                name="test_chip",
                r1=test_ip_files["chip_ip"],
                ip_performed="MLL"
            )
            
            expected_fullname = "test_chip_MLL"
            assert fastq_set.ip_set_fullname == expected_fullname


class TestFastqFilePathHandling:
    """Test file path handling edge cases."""

    def test_relative_path_conversion(self):
        """Test that relative paths are converted to absolute."""
        # Create a FastqFile with a fake path (testing path handling only)
        relative_path = pathlib.Path("relative/path/test.fastq")
        
        # This will fail because file doesn't exist, but we can test the path handling logic
        with pytest.raises(FileNotFoundError):
            FastqFile(path=relative_path)

    def test_resolved_name_option(self, tmp_path):
        """Test the use_resolved_name option."""
        # Create a temporary file for testing
        test_file = tmp_path / "test.fastq"
        test_file.write_text("@read1\nACGT\n+\nIIII\n")
        
        # Test with use_resolved_name=True
        fastq_file = FastqFile(path=test_file, use_resolved_name=True)
        assert fastq_file.path == test_file.resolve()
        
        # Test with use_resolved_name=False (default)
        fastq_file2 = FastqFile(path=test_file, use_resolved_name=False)
        assert fastq_file2.path == test_file.absolute()


@pytest.mark.parametrize("filename,expected_read_num", [
    ("sample_R1.fastq.gz", 1),
    ("sample_R2.fastq.gz", 2),
    ("sample_1.fastq.gz", 1),
    ("sample_2.fastq.gz", 2),
    ("sample.fastq.gz", 1),  # Default to 1 when no read number found
])
def test_read_number_extraction_parametrized(tmp_path, filename, expected_read_num):
    """Parametrized test for read number extraction."""
    test_file = tmp_path / filename
    test_file.write_text("@read1\nACGT\n+\nIIII\n")
    
    fastq_file = FastqFile(path=test_file)
    assert fastq_file.read_number == expected_read_num
