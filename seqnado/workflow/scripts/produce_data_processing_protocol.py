"""Generate data processing protocol text for GEO submission.

This script builds a protocol description by analyzing the pipeline configuration
and detecting which tools and processing steps were actually used.
"""

import logging
import re
import subprocess
from pathlib import Path
from typing import Any, Optional

# Try to import version from installed package, fallback to reading _version.py
try:
    from seqnado._version import __version__ as seqnado_version
except (ImportError, ModuleNotFoundError):
    # When running in Singularity container, seqnado may not be installed
    # Try to read version from _version.py file
    try:
        version_file = Path(__file__).parent.parent.parent / "_version.py"
        if version_file.exists():
            version_locals = {}
            exec(version_file.read_text(), version_locals)
            seqnado_version = version_locals.get("__version__", "unknown")
        else:
            seqnado_version = "unknown"
    except Exception:
        seqnado_version = "unknown"

# ============================================================================
# Configuration and Logger Setup
# ============================================================================

assay = snakemake.params.assay
config = snakemake.config
third_party_tools = config.get("third_party_tools") or {}
assay_config = config.get("assay_config") or {}


def configure_logger() -> logging.Logger:
    """Set up logging with optional file output."""
    logger = logging.getLogger("geo_protocol")
    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    handlers = [logging.StreamHandler()]

    log_target = None
    try:
        log_target = snakemake.log[0]
    except (AttributeError, TypeError, IndexError):
        log_target = None

    if log_target:
        log_path = Path(log_target)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path))

    for handler in handlers:
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    logger.propagate = False
    return logger


logger = configure_logger()


# ============================================================================
# Utility Functions
# ============================================================================


def safe_get(config_dict: dict, *keys, default=None) -> Any:
    """Safely navigate nested dictionary, handling None values.

    Args:
        config_dict: The dictionary to navigate
        *keys: Sequence of keys to navigate
        default: Default value if key doesn't exist or is None

    Returns:
        The value at the nested location or default
    """
    current = config_dict
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
        if current is None:
            return default
    return current if current is not None else default


def listify(value: Any) -> list:
    """Convert value to list format."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def tool_configured(tool_name: str) -> bool:
    """Check if a tool is configured in third_party_tools."""
    tool_block = third_party_tools.get(tool_name)
    if tool_block is None:
        return False
    if isinstance(tool_block, dict):
        return True
    return bool(tool_block)


def tool_arguments(*keys, default=None) -> str:
    """Get command-line arguments for a tool from config.

    Returns the command-line arguments, or 'default parameters' if empty/None.
    """
    value = safe_get(third_party_tools, *keys, default=default)
    # Handle empty strings or None - return "default parameters"
    if value is None or (isinstance(value, str) and not value.strip()):
        return "default parameters"
    return value


# ============================================================================
# Tool Version Detection
# ============================================================================

# Tool version detection configuration: (command, pattern, fallback_version)
TOOL_VERSIONS = {
    "fastqc": ("fastqc --version", r"v(\d+\.\d+\.\d+)", "0.12.1"),
    "trim_galore": ("trim_galore --version", r"version (\d+\.\d+\.\d+)", "0.6.10"),
    "star": ("STAR --version", r"(\d+\.\d+\.\d+[a-z]*)", "2.7.11a"),
    "bowtie2": ("bowtie2 --version", r"version (\d+\.\d+\.\d+)", "2.5.4"),
    "featureCounts": ("featureCounts -v", r"v(\S+)", "2.0.3"),
    "bamCoverage": ("bamCoverage --version", r"bamCoverage (\S+)", "3.5.6"),
    "picard": (
        "LC_ALL=C.UTF-8 picard MarkDuplicates --version",
        r"Version:?\s*(\d+\.\d+\.\d+)",
        "3.4.0",
    ),
    "macs2": ("macs2 --version", r"(\d[\d\.]+)$", "2.2.9"),
    "methyldackel": ("MethylDackel --version", r"(\d+\.\d+\.\d+)", "0.6.1"),
    "bcftools": ("bcftools --version", r"bcftools (\d+\.\d+)", "1.17"),
    "bamnado": ("bamnado --version", r"(\d+\.\d+\.\d+)", "0.1.0"),
    "homer": (
        # Try multiple locations where HOMER config might be found
        "{ cat $(dirname $(which homer))/../share/homer/config.txt 2>/dev/null || "
        "cat $CONDA_PREFIX/share/homer/config.txt 2>/dev/null || "
        "cat $HOME/.homer/config.txt 2>/dev/null; } | grep '^homer' | awk '{print $2}' | head -1",
        r"v(\d+\.\d+)",
        "5.1",
    ),
    "seacr": (None, None, "1.3"),  # No version command available
    "lanceotron": (None, None, "1.2.7"),  # No version command available
    "minimap2": ("minimap2 --version", r"(\d+\.\d+)", "2.30"),
    "cooler": ("cooler --version", r"(\d+\.\d+\.\d+)", "0.10.4"),
    "samtools": ("samtools --version", r"samtools (\d+\.\d+)", "1.22"),
    "salmon": ("salmon --version", r"salmon (\d+\.\d+\.\d+)", "1.10.0"),
}


def capture_version(
    cmd: str,
    pattern: Optional[str] = None,
    group: int = 0,
    fallback: Optional[str] = None,
    description: Optional[str] = None,
) -> str:
    """Execute command to capture tool version.

    Args:
        cmd: Shell command to execute
        pattern: Regex pattern to extract version
        group: Regex group number to extract
        fallback: Fallback version if command fails
        description: Tool description for logging

    Returns:
        Version string
    """
    label = description or cmd
    logger.debug("Running version command: %s", cmd)
    try:
        output = (
            subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            .decode()
            .strip()
        )
    except subprocess.CalledProcessError as exc:
        exc_output = exc.output.decode().strip() if exc.output else ""
        if exc_output and pattern:
            match = re.search(pattern, exc_output)
            if match:
                logger.warning(
                    "Command '%s' exited with status %s but produced parsable output.",
                    cmd,
                    exc.returncode,
                )
                return match.group(group)
        if fallback is not None:
            logger.warning(
                "Command '%s' failed (%s). Falling back to '%s'",
                cmd,
                exc_output or exc,
                fallback,
            )
            return fallback
        logger.error(
            "Unable to determine %s version from command '%s'. Output: %s",
            label,
            cmd,
            exc_output,
        )
        raise

    if pattern:
        match = re.search(pattern, output)
        if match:
            return match.group(group)
        if fallback is not None:
            logger.warning(
                "Pattern '%s' not found in output of '%s'. Using fallback version: %s",
                pattern,
                cmd,
                fallback,
            )
            return fallback
        logger.warning(
            "Pattern '%s' not found in output of '%s'. Using raw output.", pattern, cmd
        )
    return output


def get_tool_version(tool: str) -> str:
    """Get version for a tool using configuration."""
    if tool not in TOOL_VERSIONS:
        logger.warning("Unknown tool '%s', returning 'unknown'", tool)
        return "unknown"
    cmd, pattern, fallback = TOOL_VERSIONS[tool]
    if cmd is None:  # Tool doesn't have version command (like SEACR, lanceotron)
        logger.info(
            "Tool '%s' has no version command, using fallback version: %s",
            tool,
            fallback,
        )
        return fallback
    return capture_version(cmd, pattern, group=1, fallback=fallback, description=tool)


# ============================================================================
# Protocol Builder Classes
# ============================================================================


class ProtocolSection:
    """Base class for protocol sections."""

    def __init__(
        self, config: dict, assay_config: dict, tools: dict, logger: logging.Logger
    ):
        self.config = config
        self.assay_config = assay_config
        self.tools = tools
        self.logger = logger
        self.versions = {}

    def is_applicable(self) -> bool:
        """Check if this section should be included in the protocol."""
        return True

    def collect_versions(self) -> dict:
        """Collect tool versions needed for this section."""
        return {}

    def generate_text(self) -> list[str]:
        """Generate protocol text for this section."""
        return []


class QCSection(ProtocolSection):
    """Quality control steps."""

    def collect_versions(self) -> dict:
        return {
            "fastqc": get_tool_version("fastqc"),
            "trim_galore": get_tool_version("trim_galore"),
        }

    def generate_text(self) -> list[str]:
        trim_options = tool_arguments("trim_galore", "trim", "command_line_arguments")
        if not trim_options:
            self.logger.warning(
                "No trim_galore options found; using 'default parameters'."
            )
            trim_options = "default parameters"

        return [
            f"FASTQ files were quality checked using FastQC v{self.versions['fastqc']}.",
            f"Adapter sequences were removed and low quality reads were trimmed using trim_galore v{self.versions['trim_galore']} with the following parameters: {trim_options}.",
        ]


class AlignmentSection(ProtocolSection):
    """Read alignment steps."""

    ALIGNER_LABELS = {"star": "STAR", "bowtie2": "bowtie2"}

    def _detect_aligner(self) -> str:
        """Determine which aligner is configured."""
        index_type = safe_get(self.config, "genome", "index", "type")
        if isinstance(index_type, str):
            normalized = index_type.strip().lower()
            if normalized in self.ALIGNER_LABELS:
                return normalized

        # Fallback to checking configured tools
        for candidate in ("star", "bowtie2"):
            if tool_configured(candidate):
                return candidate

        return "bowtie2"

    def collect_versions(self) -> dict:
        aligner = self._detect_aligner()
        return {
            aligner: get_tool_version(aligner),
            "samtools": get_tool_version("samtools"),
        }

    def generate_text(self) -> list[str]:
        aligner_key = self._detect_aligner()
        aligner_label = self.ALIGNER_LABELS.get(aligner_key, aligner_key)
        aligner_version = self.versions.get(aligner_key, "unknown")
        samtools_version = self.versions.get("samtools", "unknown")

        aligner_options = tool_arguments(aligner_key, "align", "command_line_arguments")
        if not aligner_options:
            self.logger.warning(
                "No %s options found in config; using 'default parameters'.",
                aligner_label,
            )
            aligner_options = "default parameters"

        genome_name = self.config.get("genome", {}).get("name", "reference genome")

        lines = [
            f"Reads were aligned to the reference genome {genome_name} using {aligner_label} v{aligner_version} with the following parameters: {aligner_options}.",
            f"Aligned reads were sorted by coordinate using samtools sort v{samtools_version}.",
        ]

        # Add blacklist removal if enabled
        if self.config.get("remove_blacklist", False):
            blacklist_file = self.config.get("genome", {}).get("blacklist", "")
            if blacklist_file:
                lines.append(
                    "Reads mapping to blacklisted regions were removed using bedtools intersect."
                )

        return lines


class PCRDuplicatesSection(ProtocolSection):
    """PCR duplicate handling."""

    def is_applicable(self) -> bool:
        return self.config.get("pcr_duplicates", {}).get("strategy") == "remove"

    def collect_versions(self) -> dict:
        if self.is_applicable():
            return {"picard": get_tool_version("picard")}
        return {}

    def generate_text(self) -> list[str]:
        picard_options = tool_arguments(
            "picard",
            "mark_duplicates",
            "command_line_arguments",
            default="default parameters",
        )
        return [
            f"PCR duplicate reads were removed using Picard MarkDuplicates v{self.versions['picard']} with the following parameters: {picard_options}."
        ]


class Tn5ShiftSection(ProtocolSection):
    """Tn5 transposase shifting for ATAC-seq."""

    def is_applicable(self) -> bool:
        return self.assay_config.get("tn5_shift", False)

    def collect_versions(self) -> dict:
        if self.is_applicable():
            return {"bamnado": get_tool_version("bamnado")}
        return {}

    def generate_text(self) -> list[str]:
        return [
            f"Alignments were offset by +4bp (positive strand) and -5bp (negative strand) using bamnado v{self.versions['bamnado']} to account for the Tn5 transposase adaptor insertions."
        ]


class FilteringSection(ProtocolSection):
    """Alignment filtering."""

    def is_applicable(self) -> bool:
        samtools_filter_options = tool_arguments(
            "samtools", "view", "command_line_arguments", default=""
        )
        if (
            not samtools_filter_options
            or samtools_filter_options == "default parameters"
        ):
            return False

        # Check if options contain actual filtering flags
        filtering_flags = [
            "-q",
            "-f ",
            "-F ",
            "-G",
            "-L",
            "-M",
            "-r",
            "-R",
            "-d",
            "-D",
            "-s",
        ]
        return any(flag in samtools_filter_options for flag in filtering_flags)

    def generate_text(self) -> list[str]:
        samtools_filter_options = tool_arguments(
            "samtools", "view", "command_line_arguments", default=""
        )

        # Extract only filtering flags (-f, -F, -q and their values)
        filtering_parts = []

        # Match -f, -F, -q followed by their values
        for pattern in [r"-f\s+\S+", r"-F\s+\S+", r"-q\s+\S+"]:
            matches = re.findall(pattern, samtools_filter_options)
            filtering_parts.extend(matches)

        if filtering_parts:
            filter_params = " ".join(filtering_parts)
        else:
            filter_params = "default parameters"

        # Check if format conversion flags are present to adjust phrasing
        has_format_flags = any(
            flag in samtools_filter_options for flag in ["-b", "-S", "-h", "-H"]
        )

        if has_format_flags:
            return [
                f"Reads were converted to BAM format and filtered using samtools view with the following parameters: {filter_params}."
            ]
        else:
            return [
                f"Alignments were filtered using samtools view with the following parameters: {filter_params}."
            ]


class PeakCallingSection(ProtocolSection):
    """Peak calling for ChIP-seq, ATAC-seq, etc."""

    def is_applicable(self) -> bool:
        call_peaks = self.assay_config.get("call_peaks", False)
        peak_methods = listify(self.assay_config.get("peak_calling", {}).get("method"))
        return call_peaks and bool(peak_methods)

    # Configuration mapping for peak callers
    PEAK_CALLER_CONFIG = {
        "lanceotron": {
            "version_key": "lanceotron",
            "display_name": "lanceotron",
            "config_path": ("lanceotron", "call_peaks", "command_line_arguments"),
            "format_params": lambda opts: f"with the following parameters: {opts}",
        },
        "macs2": {
            "version_key": "macs2",
            "display_name": "MACS2",
            "config_path": ("macs", "call_peaks", "command_line_arguments"),
            "format_params": lambda opts: f"with the following parameters: {opts}",
            "aliases": ["macs"],
        },
        "homer": {
            "version_key": "homer",
            "display_name": "HOMER findPeaks",
            "config_path": ("homer", "find_peaks", "command_line_arguments"),
            "format_params": lambda opts: f"with the following parameters: {opts}",
        },
        "seacr": {
            "version_key": "seacr",
            "display_name": "SEACR",
            "config_path": None,  # Special case: uses threshold instead
            "format_params": lambda threshold: f"with threshold: {threshold}",
        },
    }

    def collect_versions(self) -> dict:
        versions = {}
        if not self.is_applicable():
            return versions

        peak_methods = listify(self.assay_config.get("peak_calling", {}).get("method"))
        peak_methods_normalized = {str(method).lower() for method in peak_methods}

        for method_key, config in self.PEAK_CALLER_CONFIG.items():
            # Check if this method or any of its aliases is in the config
            if method_key in peak_methods_normalized:
                versions[config["version_key"]] = get_tool_version(
                    config["version_key"]
                )
            elif "aliases" in config:
                for alias in config["aliases"]:
                    if alias in peak_methods_normalized:
                        versions[config["version_key"]] = get_tool_version(
                            config["version_key"]
                        )
                        break

        return versions

    def generate_text(self) -> list[str]:
        peak_methods = listify(self.assay_config.get("peak_calling", {}).get("method"))
        peak_methods_normalized = [str(method).lower() for method in peak_methods]
        lines = []

        for method_name in peak_methods_normalized:
            # Find the config for this method (check aliases too)
            config = None
            for key, cfg in self.PEAK_CALLER_CONFIG.items():
                if key == method_name or method_name in cfg.get("aliases", []):
                    config = cfg
                    break

            if not config:
                continue

            version = self.versions.get(config["version_key"], "unknown")
            display_name = config["display_name"]

            # Get parameters
            if config["config_path"]:
                params = tool_arguments(
                    *config["config_path"], default="default parameters"
                )
            else:
                # Special case for SEACR which uses threshold
                params = safe_get(third_party_tools, "seacr", "threshold", default=0.01)

            param_text = config["format_params"](params)
            lines.append(
                f"Peak calling was performed using {display_name} v{version} {param_text}."
            )

        return lines


class BigWigSection(ProtocolSection):
    """BigWig track generation."""

    # Configuration mapping for BigWig generators
    BIGWIG_TOOL_CONFIG = {
        "deeptools": {
            "version_key": "bamCoverage",
            "display_name": "deepTools bamCoverage",
            "config_path": ("deeptools", "bam_coverage", "command_line_arguments"),
        },
        "homer": {
            "version_key": "homer",
            "display_name": "HOMER makeBigWig",
            "config_path": ("homer", "make_bigwig", "command_line_arguments"),
        },
        "bamnado": {
            "version_key": "bamnado",
            "display_name": "bamnado",
            "config_path": ("bamnado", "bam_coverage", "command_line_arguments"),
        },
    }

    def is_applicable(self) -> bool:
        create_bigwigs = self.assay_config.get("create_bigwigs", True)
        bigwig_methods = listify(
            self.assay_config.get("bigwigs", {}).get("pileup_method")
        )
        return create_bigwigs and bool(bigwig_methods)

    def collect_versions(self) -> dict:
        versions = {}
        if not self.is_applicable():
            return versions

        bigwig_methods = listify(
            self.assay_config.get("bigwigs", {}).get("pileup_method")
        )
        bigwig_methods_normalized = {str(method).lower() for method in bigwig_methods}

        for method_key, config in self.BIGWIG_TOOL_CONFIG.items():
            if method_key in bigwig_methods_normalized:
                versions[config["version_key"]] = get_tool_version(
                    config["version_key"]
                )

        return versions

    def generate_text(self) -> list[str]:
        bigwig_config = self.assay_config.get("bigwigs", {})
        bigwig_methods = listify(bigwig_config.get("pileup_method"))
        bigwig_methods_normalized = [str(method).lower() for method in bigwig_methods]
        binsize = bigwig_config.get("bin_size", 10)
        lines = []

        for method_name in bigwig_methods_normalized:
            config = self.BIGWIG_TOOL_CONFIG.get(method_name)

            if config:
                version = self.versions.get(config["version_key"], "unknown")
                display_name = config["display_name"]
                params = tool_arguments(
                    *config["config_path"], default="default parameters"
                )
                lines.append(
                    f"BigWig files were generated using {display_name} v{version} with the following parameters: {params}."
                )
            else:
                # Unknown tool - generic message
                lines.append(
                    f"Coverage tracks were generated using {method_name} (bin size {binsize} bp)."
                )

        return lines


class RNAQuantificationSection(ProtocolSection):
    """RNA-seq quantification."""

    # Configuration mapping for RNA quantification tools
    RNA_QUANT_CONFIG = {
        "feature_counts": {
            "version_key": "featureCounts",
            "tool_check": "subread",
            "preamble": "Strands were separated using the --filterRNAstrand option.",
            "display_name": "featureCounts",
            "verb": "quantified",
            "config_path": ("subread", "feature_counts", "command_line_arguments"),
        },
        "salmon": {
            "version_key": "salmon",
            "tool_check": "salmon",
            "preamble": None,
            "display_name": "Salmon",
            "verb": "quantified",
            "config_path": ("salmon", "quant", "command_line_arguments"),
        },
    }

    def _get_method(self) -> str:
        """Get the configured RNA quantification method."""
        rna_quant = self.assay_config.get("rna_quantification") or {}
        return (rna_quant.get("method") or "").lower()

    def is_applicable(self) -> bool:
        method = self._get_method()
        return method in self.RNA_QUANT_CONFIG

    def collect_versions(self) -> dict:
        versions = {}
        method = self._get_method()
        config = self.RNA_QUANT_CONFIG.get(method)

        if config and tool_configured(config["tool_check"]):
            versions[config["version_key"]] = get_tool_version(config["version_key"])

        return versions

    def generate_text(self) -> list[str]:
        lines = []
        method = self._get_method()
        config = self.RNA_QUANT_CONFIG.get(method)

        if not config or config["version_key"] not in self.versions:
            return lines

        # Add preamble if specified
        if config["preamble"]:
            lines.append(config["preamble"])

        # Add main quantification line
        params = tool_arguments(*config["config_path"], default="default parameters")
        version = self.versions[config["version_key"]]
        display_name = config["display_name"]

        if method == "salmon":
            lines.append(
                f"Transcript abundances were {config['verb']} using {display_name} v{version} with the following parameters: {params}."
            )
        else:
            lines.append(
                f"Alignments were {config['verb']} using {display_name} v{version} with the following parameters: {params}."
            )

        return lines


class MethylationSection(ProtocolSection):
    """Methylation calling."""

    def is_applicable(self) -> bool:
        return tool_configured("methyldackel")

    def collect_versions(self) -> dict:
        if self.is_applicable():
            return {"methyldackel": get_tool_version("methyldackel")}
        return {}

    def generate_text(self) -> list[str]:
        methyldackel_options = tool_arguments(
            "methyldackel", "command_line_arguments", default="default parameters"
        )
        return [
            f"Methylation was quantified using MethylDackel v{self.versions['methyldackel']} with the following parameters: {methyldackel_options}."
        ]


class VariantCallingSection(ProtocolSection):
    """SNP/variant calling."""

    def is_applicable(self) -> bool:
        return tool_configured("bcftools")

    def collect_versions(self) -> dict:
        if self.is_applicable():
            return {"bcftools": get_tool_version("bcftools")}
        return {}

    def generate_text(self) -> list[str]:
        return [
            f"Variants were called using bcftools mpileup and bcftools call v{self.versions['bcftools']}.",
            "Multi-allelic variants were split using bcftools norm.",
        ]


class MCCSection(ProtocolSection):
    """Micro-C/Capture-C processing."""

    def is_applicable(self) -> bool:
        return self.assay_config.get("mcc") is not None

    def collect_versions(self) -> dict:
        versions = {}
        if self.is_applicable():
            if tool_configured("minimap2"):
                versions["minimap2"] = get_tool_version("minimap2")
            if tool_configured("cooler"):
                versions["cooler"] = get_tool_version("cooler")
            versions["samtools"] = get_tool_version("samtools")
            versions["bamnado"] = get_tool_version("bamnado")
            if tool_configured("lanceotron_mcc") or tool_configured("lanceotron"):
                versions["lanceotron"] = get_tool_version("lanceotron")
        return versions

    def generate_text(self) -> list[str]:
        mcc_config = self.assay_config.get("mcc", {})
        viewpoints_file = mcc_config.get("viewpoints", "viewpoints BED file")
        resolutions = mcc_config.get("resolutions", [100])
        exclusion_zone = mcc_config.get("exclusion_zone", 500)

        lines = [
            f"Viewpoint sequences were extracted and indexed using samtools faidx v{self.versions.get('samtools', 'unknown')}.",
        ]

        minimap2_options = tool_arguments(
            "minimap2", "align", "command_line_arguments", default="-x sr -a -k 8 -w 1"
        )
        lines.append(
            f"Reads were aligned to viewpoint sequences using minimap2 v{self.versions.get('minimap2', 'unknown')} with parameters: {minimap2_options}."
        )

        lines.extend(
            [
                "Viewpoint-aligned reads were extracted and their unmapped mates were re-aligned to the reference genome.",
                "Genome-mapped reads were combined to create viewpoint-to-genome contact maps.",
                "Ligation junctions were identified and sorted by genomic position.",
                f"Exclusion regions within {exclusion_zone} bp of each viewpoint were defined to avoid self-ligation artifacts.",
                f"Contact matrices were generated at resolutions of {', '.join(map(str, resolutions))} bp using cooler v{self.versions.get('cooler', 'unknown')}.",
                f"Multi-resolution contact matrices (.mcool files) were created using cooler zoomify v{self.versions.get('cooler', 'unknown')} for interactive visualization.",
            ]
        )

        # Add lanceotron peak calling if configured
        if "lanceotron" in self.versions:
            lines.append(
                f"Chromatin interactions were identified using lanceotron v{self.versions.get('lanceotron', 'unknown')} peak calling on viewpoint contact profiles."
            )

        lines.append(
            f"BigWig coverage tracks were generated using bamnado v{self.versions['bamnado']} for each viewpoint to visualize interaction landscapes."
        )

        return lines


# ============================================================================
# Main Protocol Builder
# ============================================================================


class ProtocolBuilder:
    """Main protocol builder that coordinates all sections."""

    # Define the order of sections for the protocol
    SECTION_CLASSES = [
        QCSection,
        AlignmentSection,
        PCRDuplicatesSection,
        Tn5ShiftSection,
        FilteringSection,
        BigWigSection,
        RNAQuantificationSection,
        MethylationSection,
        VariantCallingSection,
        MCCSection,  # MCC must come before PeakCalling to avoid conflicts
        PeakCallingSection,
    ]

    def __init__(self, config: dict, assay_config: dict, logger: logging.Logger):
        self.config = config
        self.assay_config = assay_config
        self.logger = logger
        self.sections = []
        self.all_versions = {}

    def build(self) -> str:
        """Build the complete protocol text."""
        self.logger.info("Generating GEO protocol content for assay '%s'", assay)

        # Initialize all sections
        for section_class in self.SECTION_CLASSES:
            section = section_class(
                self.config, self.assay_config, third_party_tools, self.logger
            )
            self.sections.append(section)

        # Collect versions from applicable sections
        for section in self.sections:
            if section.is_applicable():
                versions = section.collect_versions()
                self.all_versions.update(versions)
                section.versions = self.all_versions  # Share versions across sections

        # Log collected versions
        for tool, version in sorted(self.all_versions.items()):
            self.logger.info("%s version: %s", tool, version)

        # Generate text from all applicable sections
        content_lines = []

        # Add SeqNado version header
        content_lines.append(f"Data was processed using SeqNado v{seqnado_version}.")

        for section in self.sections:
            if section.is_applicable():
                section_text = section.generate_text()
                content_lines.extend(section_text)

        # Join and clean up the content
        content = "\n".join(content_lines)
        content = self._clean_content(content)

        return content

    def _clean_content(self, content: str) -> str:
        """Clean up the generated protocol text."""
        content = (
            content.strip()
            .replace("the following parameters: False", "with default parameters")
            .replace("\n\n", "\n")
        )
        # Remove empty lines
        content = "\n".join(
            [line.strip() for line in content.splitlines() if line.strip()]
        )
        return content


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    try:
        builder = ProtocolBuilder(config, assay_config, logger)
        protocol_text = builder.build()

        with open(snakemake.output[0], "w") as f:
            f.write(protocol_text)

        logger.info("Protocol successfully written to %s", snakemake.output[0])
    except Exception as e:
        logger.error("Failed to generate protocol: %s", e, exc_info=True)
        raise
