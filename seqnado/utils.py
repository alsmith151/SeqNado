import hashlib
import os
import pathlib
import re
from collections import defaultdict
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import snakemake
from loguru import logger
from pydantic import BaseModel, Field, computed_field, validator
from pydantic.dataclasses import dataclass
from snakemake.io import expand

FILETYPE_TO_DIR_MAPPING = {
    "tag": "tag_dirs",
    "bigwig": "bigwigs/deeptools",
    "bam": "aligned",
}

FILETYPE_TO_EXTENSION_MAPPING = {"tag": "/", "bigwig": ".bigWig", "bam": ".bam"}


def is_on(param: str) -> bool:
    """
    Returns True if parameter in "on" values
    On values:
        - true
        - t
        - on
        - yes
        - y
        - 1
    """
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_off(param: str):
    """Returns True if parameter in "off" values"""
    values = ["", "None", "none", "f", "n", "no", "false", "0"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_none(param: str) -> bool:
    """Returns True if parameter is none"""
    values = ["", "none"]
    if str(param).lower() in values:
        return True
    else:
        return False


def convert_empty_yaml_entry_to_string(param: str) -> str:
    """
    Converts empty yaml entries to string
    """
    if is_none(param):
        return ""
    else:
        return param


def format_config_dict(config: Dict) -> Dict:
    """
    Formats the config dictionary to ensure that all entries are strings.

    """
    for key, value in config.items():
        if isinstance(value, dict):
            config[key] = format_config_dict(value)
        else:
            entry = convert_empty_yaml_entry_to_string(value)

            if is_on(entry):
                config[key] = True
            elif is_off(entry):
                config[key] = False
            elif is_none(entry):
                config[key] = False
            else:
                config[key] = entry

    return config


def get_fastq_files(path: str, recursive=False) -> pd.DataFrame:
    files = (
        pathlib.Path(path).glob("**/*.fastq.gz")
        if recursive
        else pathlib.Path(path).glob("*.fastq.gz")
    )
    return files


def has_bowtie2_index(prefix: str) -> bool:
    """
    Checks if bowtie2 index is present.
    """

    path_prefix = pathlib.Path(prefix).resolve()
    path_dir = path_prefix.parent
    path_prefix_stem = path_prefix.stem

    bowtie2_indicies = list(path_dir.glob(f"{path_prefix_stem}*.bt2"))

    if len(bowtie2_indicies) > 0:
        return True


def check_options(value: object):
    if value in [None, np.nan, ""]:
        return ""
    elif is_off(value):
        return ""
    else:
        return value


def define_output_files(
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"],
    chip_spikein_normalisation: bool = False,
    sample_names: list = None,
    pileup_method: list = None,
    peak_calling_method: list = None,
    make_bigwigs: bool = False,
    call_peaks: bool = False,
    make_heatmaps: bool = False,
    make_ucsc_hub: bool = False,
    call_snps: bool = False,
    annotate_snps: bool = False,
    run_deseq2: bool = False,
    **kwargs,
) -> list:
    """Define output files for the pipeline"""

    analysis_output = [
        "seqnado_output/qc/fastq_raw_qc.html",
        "seqnado_output/qc/fastq_trimmed_qc.html",
        "seqnado_output/qc/alignment_raw_qc.html",
        "seqnado_output/qc/alignment_filtered_qc.html",
        "seqnado_output/design.csv",
    ]
    assay_output = []

    if kwargs["remove_pcr_duplicates_method"] == "picard":
        analysis_output.append("seqnado_output/qc/library_complexity_qc.html")

    if make_heatmaps:
            assay_output.extend(
                [
                    "seqnado_output/heatmap/heatmap.pdf",
                    "seqnado_output/heatmap/metaplot.pdf",
                ]
            )

    if make_ucsc_hub:
        hub_dir = kwargs["ucsc_hub_details"].get("directory")
        hub_name = kwargs["ucsc_hub_details"].get("name")
        hub_file = os.path.join(hub_dir, f"{hub_name}.hub.txt")
        analysis_output.append(hub_file)

    if assay in ["ChIP", "ATAC"]:
        if chip_spikein_normalisation:
            if assay == "ChIP":
                assay_output.extend(
                    [
                        "seqnado_output/qc/full_fastqscreen_report.html",
                        "seqnado_output/normalisation_factors.tsv",
                    ]
                )

        if make_bigwigs and pileup_method:
            assay_output.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}.bigWig",
                    sample=sample_names,
                    method=pileup_method,
                )
            )

        if call_peaks and peak_calling_method:
            if assay == "ChIP":
                assay_output.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{ip}.bed",
                        ip=kwargs["sample_names_ip"],
                        method=peak_calling_method,
                    )
                )
            else:
                assay_output.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{sample}.bed",
                        sample=sample_names,
                        method=peak_calling_method,
                    )
                )

    elif assay == "RNA":
        if make_bigwigs and pileup_method:
            assay_output.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}_{strand}.bigWig",
                    sample=sample_names,
                    method=pileup_method,
                    strand=["plus", "minus"],
                )
            )

        if run_deseq2:
            project_id = kwargs["deseq2"].get("project_id")
            assay_output.append(f"DESeq2_{project_id}.html")

        assay_output.extend(
            [
                "seqnado_output/feature_counts/read_counts.tsv",
                *expand(
                    "seqnado_output/aligned/{sample}.bam",
                    sample=sample_names,
                ),
            ]
        )

    elif assay == "SNP":
        if call_snps:
            assay_output.expand(
                "seqnado_output/variant/{sample}_filtered.anno.vcf.gz",
                sample=sample_names,
            )

        if annotate_snps:
            assay_output.append(
                expand(
                    "seqnado_output/variant/{sample}_filtered.stats.txt",
                    sample=sample_names,
                ),
            )

    analysis_output.extend(assay_output)

    return analysis_output



class FastqFile(BaseModel):
    path: pathlib.Path

    def model_post_init(self, *args):
        self.path = pathlib.Path(self.path).resolve()

        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} does not exist.")

    @computed_field
    @property
    def sample_name(self) -> str:
        return pathlib.Path(str(self.path).removesuffix(".gz")).stem

    @computed_field
    @property
    def sample_base(self) -> str:
        to_sub = {
            r"_\S\d+_": "_",
            r"_L00\d_": "_",
            r"_R?[12](_001)?$": "",
        }

        base = self.sample_name
        for pattern, rep in to_sub.items():
            base = re.sub(pattern, rep, base)

        return base

    @computed_field
    @property
    def read_number(self) -> Optional[int]:
        """
        Return the read number of the fastq file.

        Checks multiple regex patterns to find the read number.

        """

        regex_std_illumina_paired = re.compile(r".*_R?([12])(_001)?")
        regexes = [regex_std_illumina_paired]

        for regex in regexes:
            match = regex.match(self.sample_name)
            if match:
                return int(match.group(1))

        logger.warning(
            f"Could not find read number for {self.sample_name} assuming unpaired."
        )

    @computed_field
    @property
    def is_paired(self) -> bool:
        """
        Return True if the fastq file is paired.

        """
        return True if self.read_number else False

    @computed_field
    @property
    def is_lane(self) -> bool:
        """
        Return True if the fastq file is lane split.

        """
        return "_L00" in self.sample_name

    def __lt__(self, other):
        return self.path < other.path

    def __gt__(self, other):
        return self.path > other.path

    def __eq__(self, other):
        return self.path == other.path


class FastqFileIP(FastqFile):
    ip: str = Field(default=None, description="IP performed on the sample")
    is_control: bool = Field(default=None, description="Is the sample a control")

    def model_post_init(self, *args):
        if self.ip is None:
            self.ip = self.predict_ip()

        if self.is_control is None:
            self.is_control = self.predict_is_control()

    def predict_ip(self) -> Optional[str]:
        """
        Predict the IP performed on the sample.

        Uses the sample base to predict the IP performed on the sample.

        """
        try:
            return self.sample_base.split("_")[-1]
        except IndexError:
            logger.warning(f"Could not predict IP for {self.sample_base}")
            return None

    def sample_name_without_antibody(self) -> str:
        """
        Return the sample name without the antibody name.

        """
        return re.sub(f"_{self.antibody}_", "_", self.sample_name)

    def predict_is_control(self) -> bool:
        """
        Return True if the fastq file is an input.

        """

        input_substrings = ["input", "mock", "igg", "control"]
        return any([substring in self.ip.lower() for substring in input_substrings])

    @computed_field
    @property
    def sample_base_without_ip(self) -> str:
        """
        Return the sample base without the antibody name.

        """
        return re.sub(
            f"(_{self.ip})?(_S\\d+)?(_L00\\d)?(_R?[12])?(_001)?",
            "",
            self.sample_name,
        )


class AssayNonIP(BaseModel):
    name: str = Field(default=None, description="Name of the assay")
    r1: FastqFile
    r2: Optional[FastqFile] = None
    metadata: Optional[dict] = None

    @property
    def fastq_paths(self):
        return [self.r1.path, self.r2.path] if self.is_paired else [self.r1.path]

    @property
    def is_paired(self):
        return self.r2 is not None

    @classmethod
    def from_fastq_files(cls, fq: List[FastqFile], **kwargs):
        """
        Create a SampleInfo object from a list of FastqFiles.

        """

        sample_name = fq[0].sample_base

        if len(fq) == 1:
            return cls(name=sample_name, r1=fq[0], **kwargs)
        elif len(fq) == 2:
            return cls(name=sample_name, r1=fq[0], r2=fq[1], **kwargs)
        else:
            raise ValueError(f"Invalid number of fastq files for {sample_name}.")


class AssayIP(AssayNonIP):
    name: str = Field(default=None, description="Name of the assay")
    r1: FastqFileIP
    r2: Optional[FastqFileIP] = None
    metadata: Optional[dict] = None


class ExperimentIP(BaseModel):
    ip_files: AssayIP
    control_files: AssayIP = None
    name: str = Field(default=None, description="Name of the assay")
    ip: str = Field(default=None, description="IP performed on the sample")
    control: Optional[str] = Field(
        default=None, description="IP performed on the control"
    )
    metadata: Optional[dict] = None

    def model_post_init(self, *args) -> None:
        if self.name is None:
            self.name = (
                f"{self.ip_files.r1.sample_base_without_ip}_{self.ip_files.r1.ip}"
            )

        if self.ip is None:
            self.ip = self.ip_files.r1.ip

        if self.control is None and self.control_files is not None:
            self.control = self.control_files.r1.ip

    @classmethod
    def from_fastq_files(cls, fq: List[FastqFileIP], **kwargs):
        """
        Create a SampleInfo object from a list of FastqFiles.

        """

        fq_ip = [f for f in fq if not f.is_control]
        fq_control = [f for f in fq if f.is_control]

        assert len(fq_ip) > 0, "No IP files found"
        assert len(fq_ip) < 3, "Too many IP files found"
        assert len(fq_control) < 3, "Too many control files found"

        if len(fq_control) == 0:
            return cls(ip_files=AssayIP.from_fastq_files(fq_ip, **kwargs), **kwargs)
        else:
            return cls(
                ip_files=AssayIP.from_fastq_files(fq_ip, **kwargs),
                control_files=AssayIP.from_fastq_files(fq_control, **kwargs),
                **kwargs,
            )


class Design(BaseModel):
    assays: Dict[str, AssayNonIP] = Field(
        default_factory=dict,
        description="Dictionary of assay classes keyed by sample name",
    )

    @computed_field
    @property
    def sample_names(self) -> List[str]:
        return list(self.assays.keys())

    @computed_field
    @property
    def fastq_paths(self) -> List[pathlib.Path]:
        paths = []
        for assay in self.assays.values():
            paths.append(assay.r1.path)
            if assay.r2 is not None:
                paths.append(assay.r2.path)

        return paths

    def query(self, sample_name: str) -> AssayNonIP:
        return self.assays[sample_name]

    @classmethod
    def from_fastq_files(cls, fq: List[FastqFile], **kwargs):
        """
        Create a SampleInfo object from a list of FastqFiles.

        """

        sample_names = set([f.sample_base for f in fq])
        assays = {}
        for sample_name in sample_names:
            assays[sample_name] = AssayNonIP.from_fastq_files(
                [f for f in fq if f.sample_base == sample_name], **kwargs
            )

        return cls(assays=assays, **kwargs)

    @classmethod
    def from_directory(
        cls, path: pathlib.Path, metadata: Dict[str, Any] = None, **kwargs
    ):
        """
        Create a SampleInfo object from a directory of fastq files.

        """

        fq = list(pathlib.Path(path).glob("*.fastq.gz"))
        fq = sorted([FastqFile(path=f) for f in fq])
        return cls.from_fastq_files(fq, metadata=metadata, **kwargs)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, simplified: bool = True, **kwargs):
        assays = {}
        for assay_name, row in df.iterrows():
            if simplified:
                metadata = {}
                for k, v in row.items():
                    if k not in ["r1", "r2"]:
                        metadata[k] = v
                assays[assay_name] = AssayNonIP(
                    name=assay_name,
                    r1=FastqFile(path=row["r1"]),
                    r2=FastqFile(path=row["r2"]),
                    metadata=metadata,
                )
            else:
                raise NotImplementedError("Not implemented")
        return cls(assays=assays, **kwargs)

    def to_dataframe(self, simplify: bool = True) -> pd.DataFrame:
        """
        Return the SampleInfo object as a pandas DataFrame.

        """

        if not simplify:
            df = pd.DataFrame.from_dict(
                {k: v.model_dump() for k, v in self.assays.items()}, orient="index"
            )
        else:
            data = {}
            for assay_name, assay in self.assays.items():
                data[assay_name] = {}
                data[assay_name]["r1"] = assay.r1.path
                if assay.r2 is not None:
                    data[assay_name]["r2"] = assay.r2.path

                if assay.metadata is not None:
                    for k, v in assay.metadata.items():
                        data[assay_name][k] = v

            df = pd.DataFrame.from_dict(data, orient="index")

        return df


class DesignIP(BaseModel):
    assays: Dict[str, ExperimentIP] = Field(
        default_factory=dict,
        description="Dictionary of experiment classes keyed by sample name",
    )

    @computed_field
    @property
    def sample_names_ip(self) -> List[str]:
        sample_names = set()
        for experiment in self.assays.values():
            sample_names.add(experiment.ip_files.name)

        return list(sample_names)

    @property
    def sample_names_control(self) -> List[str]:
        sample_names = set()
        for experiment in self.assays.values():
            if experiment.control is not None:
                sample_names.add(experiment.control_files.name)

        return list(sample_names)

    @property
    def sample_names(self) -> List[str]:
        return self.sample_names_ip + self.sample_names_control

    @property
    def ip_names(self) -> List[str]:
        return list(set([experiment.ip for experiment in self.assays.values()]))

    @property
    def control_names(self) -> List[str]:
        return list(set([experiment.control for experiment in self.assays.values()]))

    @computed_field
    @property
    def fastq_paths(self) -> List[pathlib.Path]:
        paths = []
        for experiment in self.assays.values():
            paths.extend(experiment.ip_files.fastq_paths)

            if experiment.control_files is not None:
                paths.extend(experiment.control_files.fastq_paths)

        return paths

    def query(self, sample_name: str, ip: str, control: str = None) -> ExperimentIP:
        """
        Extract an experiment from the design.
        """

        for experiment in self.assays.values():
            if experiment.name == f"{sample_name}_{ip}":
                if control is not None:
                    if experiment.control == control:
                        return experiment
                else:
                    return experiment

        raise ValueError(
            f"Could not find experiment with sample name {sample_name} and ip {ip} and control {control}"
        )

    @classmethod
    def from_fastq_files(cls, fq: List[FastqFileIP], **kwargs):
        """
        Create a SampleInfo object from a list of FastqFiles.

        """

        samples = defaultdict(list)
        for f in fq:
            samples[f.sample_base_without_ip].append(f)

        assays = {}
        for sample_name, sample in samples.items():
            assays[sample_name] = ExperimentIP.from_fastq_files(sample, **kwargs)

        return cls(assays=assays, **kwargs)

    @classmethod
    def from_directory(
        cls, path: pathlib.Path, metadata: Dict[str, Any] = None, **kwargs
    ):
        """
        Create a SampleInfo object from a directory of fastq files.

        """

        fq = list(pathlib.Path(path).glob("*.fastq.gz"))
        fq = sorted([FastqFileIP(path=f) for f in fq])
        return cls.from_fastq_files(fq, metadata=metadata, **kwargs)

    def to_dataframe(self, simplify: bool = True):
        """
        Return the SampleInfo object as a pandas DataFrame.

        """

        data = {}
        for experiment_name, experiment in self.assays.items():
            data[experiment_name] = {}
            data[experiment_name]["ip_r1"] = experiment.ip_files.r1.path
            if experiment.ip_files.r2 is not None:
                data[experiment_name]["ip_r2"] = experiment.ip_files.r2.path

            if experiment.control_files is not None:
                data[experiment_name]["control_r1"] = experiment.control_files.r1.path
                if experiment.control_files.r2 is not None:
                    data[experiment_name][
                        "control_r2"
                    ] = experiment.control_files.r2.path

            if experiment.metadata is not None:
                for k, v in experiment.metadata.items():
                    data[experiment_name][k] = v

            data[experiment_name]["ip"] = experiment.ip
            data[experiment_name]["control"] = experiment.control

        df = pd.DataFrame.from_dict(data, orient="index")

        return df

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, simplified: bool = True, **kwargs):
        experiments = {}
        for experiment_name, row in df.iterrows():
            if simplified:
                metadata = {}
                for k, v in row.items():
                    if k not in [
                        "ip_r1",
                        "ip_r2",
                        "control_r1",
                        "control_r2",
                        "ip",
                        "control",
                    ]:
                        metadata[k] = v

                ip = row["ip"]
                control = row["control"]

                if "control_r1" not in row:
                    experiments[experiment_name] = ExperimentIP(
                        ip_files=AssayIP(
                            r1=FastqFileIP(path=row["ip_r1"]),
                            r2=FastqFileIP(path=row["ip_r2"]),
                        ),
                        ip=ip,
                        control=None,
                        metadata=metadata,
                    )
                else:
                    experiments[experiment_name] = ExperimentIP(
                        ip_files=AssayIP(
                            r1=FastqFileIP(path=row["ip_r1"]),
                            r2=FastqFileIP(path=row["ip_r2"]),
                        ),
                        control_files=AssayIP(
                            r1=FastqFileIP(path=row["control_r1"]),
                            r2=FastqFileIP(path=row["control_r2"]),
                        ),
                        ip=ip,
                        control=control,
                        metadata=metadata,
                    )
            else:
                raise NotImplementedError("Not implemented")

        return cls(assays=experiments, **kwargs)


def symlink_files(
    output_dir: pathlib.Path, assay: Union[AssayNonIP, AssayIP], assay_name: str
):
    r1_path_new = pathlib.Path(f"{output_dir}/{assay_name}_1.fastq.gz")
    r2_path_new = pathlib.Path(f"{output_dir}/{assay_name}_2.fastq.gz")

    if not r1_path_new.exists():
        try:
            r1_path_new.symlink_to(assay.r1.path.resolve())
        except FileExistsError:
            logger.warning(f"Symlink for {r1_path_new} already exists.")

    if assay.r2 and not r2_path_new.exists():
        try:
            r2_path_new.symlink_to(assay.r2.path.resolve())
        except FileExistsError:
            logger.warning(f"Symlink for {r2_path_new} already exists.")


def symlink_fastq_files(
    design: Union[Design, DesignIP], output_dir: str = "seqnado_output/fastqs/"
) -> None:
    """
    Symlink the fastq files to the output directory.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(design, Design):
        for assay_name, assay in design.assays.items():
            symlink_files(output_dir, assay, assay_name)

    elif isinstance(design, DesignIP):
        for experiment_name, experiment in design.assays.items():
            assay = experiment.ip_files
            assay_name = assay.name
            symlink_files(output_dir, assay, assay_name)

            if experiment.control_files:
                assay = experiment.control_files
                assay_name = assay.name
                symlink_files(output_dir, assay, assay_name)


def define_output_files(
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"],
    sample_names: list = None,
    pileup_method: list = None,
    peak_calling_method: list = None,
    make_bigwigs: bool = False,
    call_peaks: bool = False,
    make_heatmaps: bool = False,
    make_ucsc_hub: bool = False,
    call_snps: bool = False,
    annotate_snps: bool = False,
    **kwargs,
) -> list:
    """Define output files for the pipeline"""

    analysis_output = [
        "seqnado_output/qc/fastq_raw_qc.html",
        "seqnado_output/qc/fastq_trimmed_qc.html",
        "seqnado_output/qc/alignment_raw_qc.html",
        "seqnado_output/qc/alignment_filtered_qc.html",
        "seqnado_output/design.csv",
    ]
    assay_output = []

    if kwargs["remove_pcr_duplicates_method"] == "picard":
        analysis_output.append("seqnado_output/qc/library_complexity_qc.html")

    if make_ucsc_hub:
        hub_dir = kwargs["ucsc_hub_details"].get("directory")
        hub_name = kwargs["ucsc_hub_details"].get("name")
        hub_file = os.path.join(hub_dir, f"{hub_name}.hub.txt")
        analysis_output.append(hub_file)

    if assay in ["ChIP", "ATAC"]:
        if make_bigwigs and pileup_method:
            assay_output.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}.bigWig",
                    sample=sample_names,
                    method=pileup_method,
                )
            )

        if call_peaks and peak_calling_method:
            if assay == "ChIP":
                assay_output.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{ip}.bed",
                        ip=kwargs["ip"],
                        method=peak_calling_method,
                    )
                )
            else:
                assay_output.extend(
                    expand(
                        "seqnado_output/peaks/{method}/{sample}.bed",
                        sample=sample_names,
                        method=peak_calling_method,
                    )
                )

        if make_heatmaps:
            assay_output.extend(
                expand(
                    "seqnado_output/heatmap/{method}/{sample}.png",
                    sample=sample_names,
                    method=pileup_method,
                )
            )

    elif assay == "RNA":
        if make_bigwigs and pileup_method:
            assay_output.extend(
                expand(
                    "seqnado_output/bigwigs/{method}/{sample}_{strand}.bigWig",
                    sample=sample_names,
                    method=pileup_method,
                    strand=["plus", "minus"],
                )
            )

        if kwargs["run_deseq2"]:
            project_id = kwargs["deseq2"].get("project_id")
            assay_output.append(f"DESeq2_{project_id}.html")

        assay_output.extend(
            [
                "seqnado_output/feature_counts/read_counts.tsv",
                *expand(
                    "seqnado_output/aligned/{sample}.bam",
                    sample=sample_names,
                ),
            ]
        )

    elif assay == "SNP":
        if call_snps:
            assay_output.expand(
                "seqnado_output/variant/{sample}_filtered.anno.vcf.gz",
                sample=sample_names,
            )

        if annotate_snps:
            assay_output.append(
                expand(
                    "seqnado_output/variant/{sample}_filtered.stats.txt",
                    sample=sample_names,
                ),
            )

    analysis_output.extend(assay_output)

    return analysis_output
