from __future__ import annotations

import shutil
from pathlib import Path

import pytest

# All helpers are now in the helpers package
from helpers.config_utils import TestContext, make_test_paths
from helpers.data import ensure_fastqs_present
from helpers.genome import ensure_genome_resources
from helpers.project import (
    create_config_yaml,
    create_design_file,
    init_seqnado_project,
)
from helpers.utils import get_fastq_pattern


# Fixture for all genome resources
@pytest.fixture(scope="session")
def genome_resources(test_context: TestContext) -> dict:
    """
    Provides a function to get genome resources for an assay.
    Returns a function that takes an assay and returns a dict of resource paths.
    """

    def _get_resources(assay: str) -> dict:
        return ensure_genome_resources(test_context.test_paths.genome, assay)

    return _get_resources


@pytest.fixture(scope="function")
def seqnado_run_dir(config_yaml_for_testing: Path) -> Path:
    """Return the seqnado run directory."""
    return Path(config_yaml_for_testing).parent


@pytest.fixture(scope="session")
def test_profile_path(test_context: TestContext) -> Path:
    """Return the path to the test profile for Snakemake."""
    return test_context.test_paths.test_profile


@pytest.fixture(scope="session")
def cores(test_context: TestContext) -> int:
    """Return the number of cores to use for pipeline tests."""
    return test_context.cores


@pytest.fixture(scope="session")
def test_context(pytestconfig, tmp_path_factory) -> TestContext:
    """Session-scoped test context with paths and settings."""

    selected_assays = pytestconfig.getoption("--assays")
    if isinstance(selected_assays, str):
        selected_assays = [a.strip() for a in selected_assays.split(",") if a.strip()]
    elif not selected_assays:
        selected_assays = []
    make_test_paths(Path(__file__).resolve())

    return TestContext(pytestconfig, tmp_path_factory)


def ensure_seqnado_init(test_context, genome_resources, assay, monkeypatch):
    """
    Initialize seqnado project with required genome resources.
    Returns: (run_directory, resources) tuple
    """
    run_directory = test_context.run_directory(assay)
    resources = genome_resources(assay)  # Use genome_resources as a parameter

    # Ensure genome configuration is initialized
    init_seqnado_project(
        run_directory=run_directory,
        assay=assay,
        resources=resources,
        test_data_path=test_context.test_paths.test_data,
        monkeypatch=monkeypatch,
    )

    # Verify genome configuration file exists
    genome_config_file = run_directory / ".config" / "seqnado" / "genome_config.json"
    if not genome_config_file.exists():
        raise FileNotFoundError(f"Genome configuration file not found at {genome_config_file}")

    return run_directory, resources


@pytest.fixture(scope="function")
def config_yaml_for_testing(
    test_context: TestContext,
    assay: str,
    monkeypatch: pytest.MonkeyPatch,
    genome_resources,
) -> Path:
    """Create and patch config YAML for testing."""
    # Initialize seqnado project first and get the run directory
    run_directory, resources = ensure_seqnado_init(test_context, genome_resources, assay, monkeypatch)

    # Now create config YAML using the same run directory
    config_path = create_config_yaml(run_directory, assay, monkeypatch, resources)
    return config_path


@pytest.fixture(scope="function")
def design(test_context: TestContext, assay: str, seqnado_run_dir: Path) -> Path:
    """Generate design file for the assay in the seqnado run directory."""
    assay_type = test_context.assay_type(assay)

    # Download FASTQ files
    ensure_fastqs_present(test_context.test_paths.fastq, [assay])
    fastq_source_dir = test_context.test_paths.fastq

    # Use the correct pattern from get_fastq_pattern
    pattern = get_fastq_pattern(assay_type)
    fastqs_to_copy = list(fastq_source_dir.glob(pattern))

    if not fastqs_to_copy:
        raise FileNotFoundError(
            f"No FASTQ files found for assay '{assay_type}' in {fastq_source_dir}"
        )

    # Move FASTQs to the run directory
    fastq_dest_dir = seqnado_run_dir / "fastqs"
    fastq_dest_dir.mkdir(parents=True, exist_ok=True)
    for fq in fastqs_to_copy:
        shutil.copy2(fq, fastq_dest_dir / fq.name)

    # Generate design file
    design_file = create_design_file(
        run_directory=seqnado_run_dir,
        assay=assay_type,
    )

    return design_file


@pytest.fixture(scope="function")
def multi_assay_configs(
    tmp_path_factory, test_context, monkeypatch, request, genome_resources
):
    """
    For a list of assays, set up config and metadata files for each assay.
    Returns: {assay: {"config": config_path, "metadata": metadata_path}}
    """
    # Get the list of assays from the test function's parameters
    # When a test is parametrized with multi_assays, we need to extract it from the test node
    if hasattr(request, "node") and hasattr(request.node, "callspec"):
        multi_assays = request.node.callspec.params.get(
            "multi_assays", ["atac", "chip", "meth", "rna", "snp"]
        )
    else:
        # Default to all assays if not parametrized
        multi_assays = ["atac", "chip", "meth", "rna", "snp"]
    # Set up a run directory for the multi-assay test
    run_dir = tmp_path_factory.mktemp("multi_assay_run")
    configs = {}

    # Initialize seqnado once in the shared run_dir for all assays
    # We need to collect all resources needed for all assays first
    all_resources = {}
    for assay in multi_assays:
        all_resources[assay] = genome_resources(assay)

    # Merge all resources to ensure genome_config.json has everything needed
    # (e.g., STAR index for RNA, Bowtie2 index for ATAC, etc.)
    # Priority: RNA-specific resources (star_index, RNA GTF) take precedence
    merged_resources = {}
    rna_assay = None
    for assay in multi_assays:
        if "rna" in assay.lower():
            rna_assay = assay
            break

    # If RNA is present, use its resources as the base (to get RNA-specific GTF and STAR index)
    if rna_assay:
        merged_resources.update(all_resources[rna_assay])

    # Then merge in resources from other assays (won't overwrite RNA-specific ones)
    for assay in multi_assays:
        for key, value in all_resources[assay].items():
            if key not in merged_resources or merged_resources[key] is None:
                merged_resources[key] = value

    # Initialize with merged resources from all assays
    # Check if MCC is in the list of assays and use it for init if present
    # (init_seqnado_project needs to know about MCC to set up viewpoints)
    init_assay = next((a for a in multi_assays if "mcc" in a.lower()), multi_assays[0])
    init_seqnado_project(
        run_directory=run_dir,
        assay=init_assay,
        resources=merged_resources,
        test_data_path=test_context.test_paths.test_data,
        monkeypatch=monkeypatch,
    )

    # Ensure FASTQ files are present for all requested assays
    # Do this AFTER init_seqnado_project to avoid the directory being cleaned up
    fastq_path = test_context.test_paths.fastq
    ensure_fastqs_present(fastq_path, multi_assays)

    # Add assay-specific genome configs for each assay (except the one used for init)
    # Each assay gets its own genome config entry with assay-specific resources
    from helpers.utils import setup_genome_config
    genome_config_file = run_dir / ".config" / "seqnado" / "genome_config.json"
    genes_bed = test_context.test_paths.test_data / "hg38_genes.bed"

    for assay in multi_assays:
        if assay != init_assay:
            # Add this assay's specific configuration
            setup_genome_config(
                genome_config_file,
                star_index=all_resources[assay].get("star_index"),
                bt2_index=all_resources[assay]["bt2_index"],
                chromsizes=all_resources[assay]["chromsizes"],
                gtf=all_resources[assay]["gtf"],
                blacklist=all_resources[assay]["blacklist"],
                genes_bed=genes_bed,
                fasta=all_resources[assay].get("fasta"),
                assay=assay,
            )

    for assay in multi_assays:
        assay_type = test_context.assay_type(assay)

        # Create config YAML in the run directory
        # The genome config was already set up by init_seqnado_project with merged resources
        config_yaml = create_config_yaml(run_dir, assay, monkeypatch, all_resources[assay])

        # Now copy FASTQs to the directory where the config was created
        fastq_dest_dir = config_yaml.parent / "fastqs"
        fastq_dest_dir.mkdir(parents=True, exist_ok=True)

        # Copy FASTQ files for this assay
        fastq_source_dir = test_context.test_paths.fastq
        pattern = get_fastq_pattern(assay_type)
        fastqs_to_copy = list(fastq_source_dir.glob(pattern))

        if not fastqs_to_copy:
            raise FileNotFoundError(
                f"No FASTQ files found for assay '{assay_type}' in {fastq_source_dir}"
            )

        for fq in fastqs_to_copy:
            shutil.copy2(fq, fastq_dest_dir / fq.name)

        # Generate design file in the same directory as the config file
        design_file = create_design_file(
            run_directory=config_yaml.parent,
            assay=assay_type,
        )

        configs[assay] = {"config": config_yaml, "metadata": design_file}
    return configs


@pytest.fixture(scope="function")
def multi_assay_run_directory(multi_assay_configs: dict[str, dict[str, Path]]) -> Path:
    """
    Return the run directory for multi-assay tests.
    This extracts the run directory from the config paths created by multi_assay_configs.
    """
    # Get any config path and extract the run directory from it
    # All configs share the same run directory
    first_config = next(iter(multi_assay_configs.values()))
    config_path = first_config["config"]
    # The run directory is the parent of the config file
    # config is at: run_dir/<project_name>/config_<assay>.yaml
    return config_path.parent
