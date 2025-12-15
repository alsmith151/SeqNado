from __future__ import annotations

import shutil
import site
import subprocess
import sys
from pathlib import Path

import pytest
import yaml
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
        raise FileNotFoundError(
            f"Genome configuration file not found at {genome_config_file}"
        )

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
    run_directory, resources = ensure_seqnado_init(
        test_context, genome_resources, assay, monkeypatch
    )

    # Copy hg38_genes.bed from tests/data to test_output/data if it doesn't exist
    genes_bed_source = (
        test_context.test_paths.repo / "tests" / "data" / "hg38_genes.bed"
    )
    genes_bed = test_context.test_paths.test_data / "hg38_genes.bed"
    if not genes_bed.exists() and genes_bed_source.exists():
        shutil.copy2(genes_bed_source, genes_bed)

    # Copy plotting_coordinates.bed from tests/data to test_output/data if it doesn't exist
    plot_coords_source = (
        test_context.test_paths.repo / "tests" / "data" / "plotting_coordinates.bed"
    )
    plot_coords = test_context.test_paths.test_data / "plotting_coordinates.bed"
    if not plot_coords.exists() and plot_coords_source.exists():
        shutil.copy2(plot_coords_source, plot_coords)

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
def multiomics_configs(
    tmp_path_factory, test_context, monkeypatch, request, genome_resources
):
    """
    For a list of assays, set up config and metadata files for each assay.
    Returns: {assay: {"config": config_path, "metadata": metadata_path}}
    """
    # Get the list of assays from the test function's parameters
    # When a test is parametrized with multiomics, we need to extract it from the test node
    if hasattr(request, "node") and hasattr(request.node, "callspec"):
        multiomics = request.node.callspec.params.get(
            "multiomics", ["atac", "chip", "meth", "rna", "snp"]
        )
    else:
        # Default to all assays if not parametrized
        multiomics = ["atac", "chip", "meth", "rna", "snp"]
    # Set up a run directory for the Multiomic test
    run_dir = tmp_path_factory.mktemp("multiomics_run")
    configs = {}

    # Initialize seqnado once in the shared run_dir for all assays
    # We need to collect all resources needed for all assays first
    all_resources = {}
    for assay in multiomics:
        all_resources[assay] = genome_resources(assay)

    # Merge all resources to ensure genome_config.json has everything needed
    # (e.g., STAR index for RNA, Bowtie2 index for ATAC, etc.)
    # Priority: RNA-specific resources (star_index, RNA GTF) take precedence
    merged_resources = {}
    rna_assay = None
    for assay in multiomics:
        if "rna" in assay.lower():
            rna_assay = assay
            break

    # If RNA is present, use its resources as the base (to get RNA-specific GTF and STAR index)
    if rna_assay:
        merged_resources.update(all_resources[rna_assay])

    # Then merge in resources from other assays (won't overwrite RNA-specific ones)
    for assay in multiomics:
        for key, value in all_resources[assay].items():
            if key not in merged_resources or merged_resources[key] is None:
                merged_resources[key] = value

    # Initialize with merged resources from all assays
    # Check if MCC is in the list of assays and use it for init if present
    # (init_seqnado_project needs to know about MCC to set up viewpoints)
    init_assay = next((a for a in multiomics if "mcc" in a.lower()), multiomics[0])
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
    ensure_fastqs_present(fastq_path, multiomics)

    # Add assay-specific genome configs for each assay (except the one used for init)
    # Each assay gets its own genome config entry with assay-specific resources
    from helpers.utils import setup_genome_config

    genome_config_file = run_dir / ".config" / "seqnado" / "genome_config.json"

    # Copy hg38_genes.bed from tests/data to test_output/data if it doesn't exist
    genes_bed_source = (
        test_context.test_paths.repo / "tests" / "data" / "hg38_genes.bed"
    )
    genes_bed = test_context.test_paths.test_data / "hg38_genes.bed"
    if not genes_bed.exists() and genes_bed_source.exists():
        shutil.copy2(genes_bed_source, genes_bed)

    # Copy plotting_coordinates.bed from tests/data to test_output/data if it doesn't exist
    plot_coords_source = (
        test_context.test_paths.repo / "tests" / "data" / "plotting_coordinates.bed"
    )
    plot_coords = test_context.test_paths.test_data / "plotting_coordinates.bed"
    if not plot_coords.exists() and plot_coords_source.exists():
        shutil.copy2(plot_coords_source, plot_coords)

    for assay in multiomics:
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

    # Create assay-specific fastq subdirectories and copy FASTQ files
    # This creates the fastqs/{assay}/ structure expected by multiomics mode
    for assay in multiomics:
        assay_type = test_context.assay_type(assay)

        # Create assay-specific fastq subdirectory for multiomics structure
        # fastqs/{assay}/
        fastq_dest_dir = run_dir / "fastqs" / assay
        fastq_dest_dir.mkdir(parents=True, exist_ok=True)

        # Copy FASTQ files for this assay into the assay-specific directory
        fastq_source_dir = test_context.test_paths.fastq
        pattern = get_fastq_pattern(assay_type)
        fastqs_to_copy = list(fastq_source_dir.glob(pattern))

        if not fastqs_to_copy:
            raise FileNotFoundError(
                f"No FASTQ files found for assay '{assay_type}' in {fastq_source_dir}"
            )

        for fq in fastqs_to_copy:
            shutil.copy2(fq, fastq_dest_dir / fq.name)

    # Set up apptainer bind mounts for genome data
    # This is needed so the container can access genome files (bt2 indexes, etc.)

    # Get the genome directory that needs to be mounted
    first_assay = multiomics[0]
    genome_dir = Path(all_resources[first_assay]["bt2_index"]).parent.parent.resolve()
    test_data_dir = genome_dir.parent.resolve()

    # Find and update the test profile configuration
    seqnado_paths = [p for p in sys.path if "seqnado" in p and "site-packages" in p]
    if not seqnado_paths:
        # Fall back to searching site-packages
        for site_pkg in site.getsitepackages():
            test_profile_config = (
                Path(site_pkg)
                / "seqnado"
                / "workflow"
                / "envs"
                / "profiles"
                / "profile_test"
                / "config.v8+.yaml"
            )
            if test_profile_config.exists():
                break
    else:
        test_profile_config = (
            Path(seqnado_paths[0])
            / "seqnado"
            / "workflow"
            / "envs"
            / "profiles"
            / "profile_test"
            / "config.v8+.yaml"
        )

    # Update the profile config with apptainer-args
    if test_profile_config.exists():
        with open(test_profile_config) as f:
            profile_config = yaml.safe_load(f)

        # Bind mount test_data_dir so genome files are accessible in the container
        bind_arg = f"--bind {test_data_dir}:{test_data_dir}"
        if "apptainer-args" in profile_config:
            # Check if this bind mount is already present
            if str(test_data_dir) not in profile_config["apptainer-args"]:
                profile_config["apptainer-args"] += f" {bind_arg}"
        else:
            profile_config["apptainer-args"] = bind_arg

        with open(test_profile_config, "w") as f:
            yaml.dump(profile_config, f, sort_keys=False)

    # Set up fastq_screen.conf before running seqnado config
    genome_path = genome_dir
    fastq_screen_source = genome_path / "fastq_screen.conf"

    # Copy to .config/seqnado/ where seqnado expects it
    seqnado_config_dir = run_dir / ".config" / "seqnado"
    fastq_screen_dest = seqnado_config_dir / "fastq_screen.conf"

    # Ensure the fastq_screen.conf exists at the source location
    if fastq_screen_source.exists():
        shutil.copy2(fastq_screen_source, fastq_screen_dest)
    else:
        # If it doesn't exist at genome level, create a minimal config
        fastq_screen_dest.parent.mkdir(parents=True, exist_ok=True)
        with open(fastq_screen_dest, "w") as f:
            # Add entries for all assays' bt2 indexes
            seen_indexes = set()
            for assay in multiomics:
                bt2_index = all_resources[assay]["bt2_index"]
                # Avoid duplicate entries for the same index
                if bt2_index not in seen_indexes:
                    seen_indexes.add(bt2_index)
                    # Use hg38 as the database name (or assay name for special cases)
                    db_name = "hg38"
                    f.write(f"DATABASE\t{db_name}\t{bt2_index}\n")

    # Now use seqnado config in multiomics mode (--no-interactive)
    # This will create config_*.yaml for each assay and config_multiomics.yaml
    result = subprocess.run(
        ["seqnado", "config", "--no-interactive", "--no-make-dirs"],
        cwd=run_dir,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"seqnado config (multiomics) failed:\nSTDERR: {result.stderr}\nSTDOUT: {result.stdout}"
    )

    # Update generated configs to enable fastq_screen with correct paths
    import yaml

    for assay in multiomics:
        config_file = run_dir / f"config_{assay}.yaml"
        with open(config_file) as f:
            config = yaml.safe_load(f)

        # Enable fastq_screen and set the config path
        if "qc" in config:
            config["qc"]["run_fastq_screen"] = True
        if "genome" in config:
            config["genome"]["fastq_screen_config"] = str(fastq_screen_dest)
        if (
            "third_party_tools" in config
            and "fastq_screen" in config["third_party_tools"]
        ):
            config["third_party_tools"]["fastq_screen"]["config"] = str(
                fastq_screen_dest
            )

        # Write updated config
        with open(config_file, "w") as f:
            yaml.dump(config, f, sort_keys=False)

    # Use seqnado design in multiomics mode (--no-interactive)
    # This will create metadata_*.csv for each assay
    result = subprocess.run(
        ["seqnado", "design", "--no-interactive", "--accept-all-defaults"],
        cwd=run_dir,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"seqnado design (multiomics) failed:\nSTDERR: {result.stderr}\nSTDOUT: {result.stdout}"
    )

    # Collect the generated config and metadata files
    for assay in multiomics:
        config_yaml = run_dir / f"config_{assay}.yaml"
        design_file = run_dir / f"metadata_{assay}.csv"

        assert config_yaml.exists(), (
            f"config_{assay}.yaml not created by seqnado config"
        )
        assert design_file.exists(), (
            f"metadata_{assay}.csv not created by seqnado design"
        )

        configs[assay] = {"config": config_yaml, "metadata": design_file}

    # Verify multiomics config was created
    multiomics_config = run_dir / "config_multiomics.yaml"
    assert multiomics_config.exists(), (
        "config_multiomics.yaml not created by seqnado config"
    )

    return configs


@pytest.fixture(scope="function")
def multiomics_run_directory(multiomics_configs: dict[str, dict[str, Path]]) -> Path:
    """
    Return the run directory for Multiomic tests.
    This extracts the run directory from the config paths created by multiomics_configs.
    """
    first_config = next(iter(multiomics_configs.values()))
    config_path = first_config["config"]
    return config_path.parent
