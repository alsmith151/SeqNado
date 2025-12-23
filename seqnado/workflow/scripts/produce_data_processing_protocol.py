import logging
import re
import subprocess
from pathlib import Path


assay = snakemake.params.assay
config = snakemake.config


def configure_logger() -> logging.Logger:
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


def capture_version(cmd: str, pattern: str | None = None, group: int = 0,
                    fallback: str | None = None, description: str | None = None) -> str:
    label = description or cmd
    logger.debug("Running version command: %s", cmd)
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode().strip()
    except subprocess.CalledProcessError as exc:
        exc_output = exc.output.decode().strip() if exc.output else ""
        if fallback is not None:
            logger.warning("Command '%s' failed (%s). Falling back to '%s'", cmd, exc_output or exc, fallback)
            return fallback
        logger.error("Unable to determine %s version from command '%s'. Output: %s", label, cmd, exc_output)
        raise

    if pattern:
        match = re.search(pattern, output)
        if match:
            return match.group(group)
        logger.warning("Pattern '%s' not found in output of '%s'. Using raw output.", pattern, cmd)
    return output


def fastqc_version() -> str:
    return capture_version("fastqc --version", r'v\d+\.\d+\.\d+')


def trim_galore_version() -> str:
    return capture_version("trim_galore --version", r'version \d+\.\d+\.\d+')


def star_version() -> str:
    return capture_version("STAR --version")


def bowtie2_version() -> str:
    return capture_version("bowtie2 --version", r'bowtie2-align-s version \d+\.\d+\.\d+')


def featureCounts_version() -> str:
    return capture_version("featureCounts -v", r'featureCounts v(\S+)', group=1, description="featureCounts")


def bamCoverage_version() -> str:
    return capture_version("bamCoverage --version")


def picard_version() -> str:
    version = capture_version("picard MarkDuplicates --version", fallback="3.1.1", description="Picard MarkDuplicates")
    return version.removeprefix('Version:')


def lanceotron_version() -> str:
    return capture_version("lanceotron --version")


def macs2_version() -> str:
    return capture_version("macs2 --version", r'(\d[\d\.]+)$', group=1)


logger = configure_logger()
logger.info("Generating GEO protocol content for assay '%s'", assay)

versions = {
    "fastqc": fastqc_version(),
    "trim_galore": trim_galore_version(),
    "bamCoverage": bamCoverage_version(),
}

trim_options = config.get('trim_galore', {}).get('command_line_arguments')
if trim_options is None:
    logger.warning("No trim_galore options found in config; using 'default parameters'.")
    trim_options = "default parameters"

if assay == 'RNA':
    versions["star"] = star_version()
else:
    versions["bowtie2"] = bowtie2_version()

if config.get('remove_pcr_duplicates_method'):
    versions["picard"] = picard_version()

if assay == 'RNA':
    versions["featureCounts"] = featureCounts_version()

peak_methods = config.get('peak_calling_method', [])
if assay in ['ChIP', 'ATAC', 'CUT&TAG'] and config.get('call_peaks'):
    if "lanceotron" in peak_methods:
        versions["lanceotron"] = lanceotron_version()
    if 'macs2' in peak_methods:
        versions["macs2"] = macs2_version()

for tool, version in versions.items():
    logger.info("%s version: %s", tool, version)

aligner = 'STAR' if assay == 'RNA' else 'bowtie2'
aligner_version = versions['star'] if assay == 'RNA' else versions['bowtie2']

aligner_key = 'star' if assay == 'RNA' else 'bowtie2'
aligner_options = config.get(aligner_key, {}).get('command_line_arguments')
if aligner_options is None:
    logger.warning("No %s options found in config; using 'default parameters'.", aligner_key)
    aligner_options = "default parameters"

content = f"""
FASTQ files were quality checked using FastQC version {versions['fastqc']}.
Adapter sequences were removed and low quality reads were trimmed using trim_galore {versions['trim_galore']} using the following parameters: {trim_options}.
Reads were aligned to the reference genome {config['genome']['name']} using {aligner} v{aligner_version} with the following parameters: {aligner_options}.
"""

if config.get('remove_pcr_duplicates_method'):
    picard_options = config.get('picard', {}).get('command_line_arguments', 'default parameters')
    content += f"""Duplicate reads were removed using Picard MarkDuplicates v{versions['picard']} with the following parameters: {picard_options}"""


deeptools_options = config.get('deeptools', {}).get('command_line_arguments', 'default parameters')
content += f"""
BigWig files were generated using deepTools {versions['bamCoverage']} with the following parameters: {deeptools_options}.
"""

if assay == 'RNA':
    content += """Strands were separated using the --filterRNAstrand option. """
    featurecounts_options = config.get('featurecounts', {}).get('command_line_arguments', 'default parameters')
    content += f"""Alignments were quantified using featureCounts v{versions['featureCounts']} with the following parameters: {featurecounts_options}"""


if assay in ['ChIP', 'ATAC', 'CUT&TAG'] and config.get('call_peaks'):
    if "lanceotron" in peak_methods:
        lanceotron_ver = versions.get('lanceotron', '1.2.6')
        lanceotron_options = config.get('lanceotron', {}).get('command_line_arguments', 'default parameters')
        content += f"""Peak calling was performed using lanceotron v{lanceotron_ver} with the following parameters: {lanceotron_options}"""

    if 'macs2' in peak_methods:
        macs_options = config.get('macs2', {}).get('command_line_arguments', 'default parameters')
        content += f"""Peak calling was performed using MACS2 v{versions.get('macs2', 'unknown')} with the following parameters: {macs_options}"""


content = content.strip().replace('the following parameters: False', 'with default parameters').replace('/n/n', '/n')
content = "\n".join([line.strip('\n') for line in content.splitlines() if line.strip()])


with open(snakemake.output[0], 'w') as f:
    f.write(content)


logger.info("Protocol successfully written to %s", snakemake.output[0])

