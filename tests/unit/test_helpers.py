"""Tests for seqnado.helpers module."""

import glob
import os
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch, mock_open

import pandas as pd
import pytest

from seqnado import DataScalingTechnique
from seqnado.helpers import (
    define_memory_requested,
    define_time_requested,
    extract_apptainer_args,
    extract_cores_from_options,
    extract_viewpoints,
    get_scale_method,
    pepe_silvia,
    remove_unwanted_run_files,
    run_batch_job_on_error,
    viewpoint_to_grouped_viewpoint,
)


class TestExtractCoresFromOptions:
    """Tests for extract_cores_from_options function."""

    def test_extract_cores_with_c_flag(self):
        """Test extraction with -c flag."""
        options = ["-c", "8", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 8
        assert result_options == ["--other", "value"]

    def test_extract_cores_with_cores_flag(self):
        """Test extraction with --cores flag."""
        options = ["--cores", "16", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 16
        assert result_options == ["--other", "value"]

    def test_extract_cores_no_flag(self):
        """Test default behavior when no core flag is present."""
        options = ["--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 1
        assert result_options == ["--other", "value"]

    def test_extract_cores_flag_no_value(self):
        """Test behavior when core flag is present but no value given."""
        options = ["-c", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 1
        # When int() conversion fails, falls back to ValueError handler which doesn't remove -c
        assert result_options == ["-c", "--other", "value"]

    def test_extract_cores_multiple_flags(self):
        """Test that first occurrence of core flag is used."""
        options = ["-c", "4", "--other", "value"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 4
        assert result_options == ["--other", "value"]

    def test_extract_cores_flag_at_end(self):
        """Test behavior when core flag is at the end of options (IndexError)."""
        options = ["--other", "value", "-c"]
        result_options, cores = extract_cores_from_options(options)
        assert cores == 1
        # The -c flag should be removed when IndexError occurs
        assert "-c" not in result_options
        assert result_options == ["--other", "value"]


class TestExtractApptainerArgs:
    """Tests for extract_apptainer_args function."""

    def test_extract_apptainer_args_present(self):
        """Test extraction when apptainer args are present."""
        options = ["--apptainer-args", "-B /data:/data", "--other", "value"]
        result_options, apptainer_args = extract_apptainer_args(options)
        assert apptainer_args == "-B /data:/data"
        assert result_options == ["--other", "value"]

    def test_extract_apptainer_args_absent(self):
        """Test extraction when apptainer args are absent."""
        options = ["--other", "value"]
        result_options, apptainer_args = extract_apptainer_args(options)
        assert apptainer_args == ""
        assert result_options == ["--other", "value"]

    def test_extract_apptainer_args_empty_string(self):
        """Test extraction with empty string argument."""
        options = ["--apptainer-args", "", "--other", "value"]
        result_options, apptainer_args = extract_apptainer_args(options)
        assert apptainer_args == ""
        assert result_options == ["--other", "value"]


class TestDefineMemoryRequested:
    """Tests for define_memory_requested function."""

    def test_default_memory(self):
        """Test default memory allocation (1 attempt, 1G initial)."""
        result = define_memory_requested()
        assert result == "1G"

    def test_memory_scaling_with_attempts(self):
        """Test memory doubles with each attempt."""
        assert define_memory_requested(attempts=1, initial_value=2) == "2G"
        assert define_memory_requested(attempts=2, initial_value=2) == "4G"
        assert define_memory_requested(attempts=3, initial_value=2) == "8G"
        assert define_memory_requested(attempts=4, initial_value=2) == "16G"

    def test_memory_with_scale_factor(self):
        """Test memory with scale factor."""
        result = define_memory_requested(attempts=1, initial_value=4, scale=2.5)
        assert result == "10G"

    def test_memory_avoids_decimals(self):
        """Test that memory values are always integers."""
        result = define_memory_requested(attempts=1, initial_value=3, scale=1.5)
        assert result == "4G"  # int(3 * 1.5) = 4
        assert "." not in result


class TestDefineTimeRequested:
    """Tests for define_time_requested function."""

    def test_default_time(self):
        """Test default time allocation (1 attempt, 1 hour)."""
        result = define_time_requested()
        assert result == "1.0h"

    def test_time_scaling_with_attempts(self):
        """Test time doubles with each attempt."""
        assert define_time_requested(attempts=1, initial_value=2) == "2.0h"
        assert define_time_requested(attempts=2, initial_value=2) == "4.0h"
        assert define_time_requested(attempts=3, initial_value=2) == "8.0h"

    def test_time_with_scale_factor(self):
        """Test time with scale factor."""
        result = define_time_requested(attempts=1, initial_value=4, scale=2.5)
        assert result == "10.0h"

    def test_time_allows_decimals(self):
        """Test that time can have decimal values."""
        result = define_time_requested(attempts=1, initial_value=3, scale=1.5)
        assert result == "4.5h"


class TestPepeSilvia:
    """Tests for pepe_silvia function."""

    def test_pepe_silvia_returns_url(self):
        """Test that pepe_silvia returns a URL."""
        with patch("builtins.print") as mock_print:
            result = pepe_silvia()
            assert "http" in result
            assert ".jpg" in result
            mock_print.assert_called_once_with("PEPE SILVIA")


class TestGetScaleMethod:
    """Tests for get_scale_method function."""

    def test_scale_method_default(self):
        """Test default scale method (unscaled only)."""
        config = {}
        result = get_scale_method(config)
        assert result == ["unscaled"]  # Returned as lowercase string from .value

    def test_scale_method_with_spikein(self):
        """Test scale method with spikein."""
        config = {"spikein": True}
        result = get_scale_method(config)
        assert "unscaled" in result
        assert "spikein" in result
        assert len(result) == 2

    def test_scale_method_with_scale(self):
        """Test scale method with csaw scaling."""
        config = {"scale": True}
        result = get_scale_method(config)
        assert "unscaled" in result
        assert "csaw" in result
        assert len(result) == 2

    def test_scale_method_spikein_takes_precedence(self):
        """Test that spikein takes precedence over scale."""
        config = {"spikein": True, "scale": True}
        result = get_scale_method(config)
        assert "spikein" in result
        assert "csaw" not in result


class TestRemoveUnwantedRunFiles:
    """Tests for remove_unwanted_run_files function."""

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_slurm_files(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of slurm output files."""
        mock_glob.side_effect = [
            ["slurm-123.out", "slurm-456.out"],  # slurm files
            [],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = False

        remove_unwanted_run_files()

        assert mock_remove.call_count == 2

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_sps_files(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of sps files."""
        mock_glob.side_effect = [
            [],  # slurm files
            ["sps-temp1", "sps-temp2"],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = False

        remove_unwanted_run_files()

        assert mock_remove.call_count == 2

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_simg_files(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of simg files."""
        mock_glob.side_effect = [
            [],  # slurm files
            [],  # sps files
            ["container.simg"],  # simg files
        ]
        mock_isdir.return_value = False

        remove_unwanted_run_files()

        mock_remove.assert_called_once_with("container.simg")

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_remove_directory(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test removal of directories."""
        mock_glob.side_effect = [
            [],  # slurm files
            ["sps-temp-dir"],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = True

        remove_unwanted_run_files()

        mock_rmtree.assert_called_once_with("sps-temp-dir")
        mock_remove.assert_not_called()

    @patch("glob.glob")
    @patch("os.path.isdir")
    @patch("os.remove")
    @patch("shutil.rmtree")
    def test_handles_removal_errors(self, mock_rmtree, mock_remove, mock_isdir, mock_glob):
        """Test that removal errors are caught and don't crash."""
        mock_glob.side_effect = [
            ["slurm-123.out"],  # slurm files
            [],  # sps files
            [],  # simg files
        ]
        mock_isdir.return_value = False
        mock_remove.side_effect = PermissionError("Cannot remove")

        # Should not raise exception
        with patch("builtins.print") as mock_print:
            remove_unwanted_run_files()
            mock_print.assert_called()


class TestRunBatchJobOnError:
    """Tests for run_batch_job_on_error function."""

    @patch("subprocess.run")
    @patch("builtins.open", new_callable=mock_open)
    def test_creates_slurm_script(self, mock_file, mock_subprocess):
        """Test that slurm script is created with correct content."""
        email = "test@example.com"
        mock_subprocess.return_value = Mock()

        run_batch_job_on_error(email)

        mock_file.assert_called_once_with("error_notification.sh", "w")
        written_content = "".join(
            call.args[0] for call in mock_file().write.call_args_list
        )
        assert email in written_content
        assert "SBATCH" in written_content
        assert "seqnado_error_notification" in written_content

    @patch("subprocess.run")
    @patch("builtins.open", new_callable=mock_open)
    def test_submits_sbatch_job(self, mock_file, mock_subprocess):
        """Test that sbatch command is executed."""
        email = "test@example.com"
        mock_subprocess.return_value = Mock()

        run_batch_job_on_error(email)

        mock_subprocess.assert_called_once_with(
            ["sbatch", "error_notification.sh"], check=True
        )

    @patch("subprocess.run")
    @patch("builtins.open", new_callable=mock_open)
    def test_handles_sbatch_failure(self, mock_file, mock_subprocess):
        """Test handling of sbatch submission failure."""
        email = "test@example.com"
        mock_subprocess.side_effect = subprocess.CalledProcessError(1, "sbatch")

        # Should not raise exception
        run_batch_job_on_error(email)


class TestExtractViewpoints:
    """Tests for extract_viewpoints function."""

    def test_extract_viewpoints_basic_bed(self, tmp_path):
        """Test extraction from basic BED file."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "chr1\t1000\t2000\tViewpoint1\t100\t+\n"
        bed_content += "chr2\t3000\t4000\tViewpoint2\t200\t-\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 2
        assert "Viewpoint1-chr1-1000-2000" in result
        assert "Viewpoint2-chr2-3000-4000" in result

    def test_extract_viewpoints_with_coordinates_in_name(self, tmp_path):
        """Test extraction when name already contains coordinates."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "chr1\t1000\t2000\tVP1-chr1-1000-2000\t100\t+\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 1
        # Should use the existing coordinate format
        assert "VP1-chr1-1000-2000" in result

    def test_extract_viewpoints_three_column_bed(self, tmp_path):
        """Test extraction from minimal 3-column BED file."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "chr1\t1000\t2000\n"
        bed_content += "chr2\t3000\t4000\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 2
        # Name should be auto-generated from coordinates
        assert any("chr1:1000-2000" in vp for vp in result)
        assert any("chr2:3000-4000" in vp for vp in result)

    def test_extract_viewpoints_with_comments(self, tmp_path):
        """Test extraction with comment lines in BED file."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "# This is a comment\n"
        bed_content += "chr1\t1000\t2000\tViewpoint1\t100\t+\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 1
        assert "Viewpoint1-chr1-1000-2000" in result

    def test_extract_viewpoints_invalid_file(self):
        """Test error handling for invalid file."""
        with pytest.raises(ValueError, match="Error reading BED file"):
            extract_viewpoints("/nonexistent/file.bed")


class TestViewpointToGroupedViewpoint:
    """Tests for viewpoint_to_grouped_viewpoint function."""

    def test_group_single_viewpoint(self):
        """Test grouping of single viewpoint with coordinates."""
        viewpoints = ["VP1-chr1-1000-2000"]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result == {"VP1-chr1-1000-2000": "VP1"}

    def test_group_multiple_viewpoints(self):
        """Test grouping of multiple viewpoints."""
        viewpoints = [
            "VP1-chr1-1000-2000",
            "VP1-chr2-3000-4000",
            "VP2-chr3-5000-6000",
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result["VP1-chr1-1000-2000"] == "VP1"
        assert result["VP1-chr2-3000-4000"] == "VP1"
        assert result["VP2-chr3-5000-6000"] == "VP2"

    def test_group_viewpoints_different_chromosomes(self):
        """Test grouping with X, Y, M, MT chromosomes."""
        viewpoints = [
            "VP1-chrX-1000-2000",
            "VP2-chrY-3000-4000",
            "VP3-chrM-5000-6000",
            "VP4-chrMT-7000-8000",
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result["VP1-chrX-1000-2000"] == "VP1"
        assert result["VP2-chrY-3000-4000"] == "VP2"
        assert result["VP3-chrM-5000-6000"] == "VP3"
        assert result["VP4-chrMT-7000-8000"] == "VP4"

    def test_group_viewpoints_without_coordinates(self):
        """Test that viewpoints without coordinate pattern are not included."""
        viewpoints = [
            "VP1-chr1-1000-2000",
            "SimpleViewpoint",  # No coordinate pattern
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert "VP1-chr1-1000-2000" in result
        assert "SimpleViewpoint" not in result

    def test_group_empty_list(self):
        """Test grouping with empty viewpoint list."""
        viewpoints = []
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result == {}

    def test_group_viewpoints_complex_names(self):
        """Test grouping with complex viewpoint names containing hyphens."""
        viewpoints = [
            "Complex-VP-Name-chr1-1000-2000",
            "Another-Complex-VP-chr2-3000-4000",
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result["Complex-VP-Name-chr1-1000-2000"] == "Complex-VP-Name"
        assert result["Another-Complex-VP-chr2-3000-4000"] == "Another-Complex-VP"
