"""Helper functions for pipeline tests.

This package contains reusable helper functions for pipeline tests:
- config_utils.py: TestContext and test path management
- data.py: Data download and preparation functions
- genome.py: Genome resource management
- multi_assay.py: Multi-assay pytest fixtures
- project.py: Project initialization and configuration
- utils.py: Common utilities (download, extract, patterns, config setup)
"""

from .config_utils import TestContext, make_test_paths
from .data import ensure_all_genome_data, ensure_fastqs_present
from .genome import ensure_genome_resources
from .multi_assay import multi_assay_configs, multi_assay_run_directory
from .project import (
    copy_fastqs,
    create_config_yaml,
    create_design_file,
    init_seqnado_project,
    patch_config_yaml,
)
from .utils import (
    download_with_retry,
    extract_tar,
    get_fastq_pattern,
    setup_genome_config,
)

__all__ = [
    # config_utils.py
    "TestContext",
    "make_test_paths",
    # data.py
    "ensure_all_genome_data",
    "ensure_fastqs_present",
    # genome.py
    "ensure_genome_resources",
    # multi_assay.py
    "multi_assay_configs",
    "multi_assay_run_directory",
    # project.py
    "copy_fastqs",
    "create_config_yaml",
    "create_design_file",
    "init_seqnado_project",
    "patch_config_yaml",
    # utils.py
    "download_with_retry",
    "extract_tar",
    "get_fastq_pattern",
    "setup_genome_config",
]
