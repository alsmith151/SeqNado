import os
import pathlib
import re
import sys
from typing import Any, Dict, List, Literal, LiteralString, Optional, Union

import numpy as np
import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, computed_field, field_validator
from snakemake.io import expand
import pandera
from pandera.typing import Index, DataFrame, Series


def is_path(path: Optional[Union[str, pathlib.Path]]) -> Optional[pathlib.Path]:
    if isinstance(path, str):
        p = pathlib.Path(path)
    elif isinstance(path, pathlib.Path):
        p = path
    else:
        p = None

    if p is not None:
        return True
    else:
        return False


class FastqFile(BaseModel):
    path: pathlib.Path

    def model_post_init(self, *args):
        self.path = pathlib.Path(self.path).resolve()

        if not self.path.exists() or str(self.path) in ["-", ".", "", None]:
            raise FileNotFoundError(f"{self.path} does not exist.")

    @property
    def stem(self):
        return pathlib.Path(str(self.path).removesuffix(".gz")).stem

    @computed_field
    @property
    def sample_name(self) -> str:
        name = pathlib.Path(str(self.path).removesuffix(".gz")).stem
        if name.endswith("_001"):
            name = name.removesuffix("_001")
        return name

    @computed_field
    @property
    def sample_base(self) -> str:
        to_sub = {
            r"_S\d+_": "_",
            r"_L00\d_": "_",
            r"_R?[12](_001)?$": "_",
            r"__": "_",
            r"_$": "",
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
        pattern = f"(_{self.ip})?(_S\\d+)?(_L00\\d)?(_R?[12])?(_001)?"
        base = re.sub(
            pattern,
            "",
            self.sample_name,
        )

        return base


class Metadata(BaseModel):
    deseq2: Optional[str] = None
    merge: Optional[str] = None
    scale_group: Union[str, int] = "all"

    @field_validator("deseq2", "merge")
    @classmethod
    def prevent_none(cls, v):
        none_vals = [
            None,
            "None",
            "none",
            "null",
            "Null",
            "NULL",
            ".",
            "",
            "NA",
            np.nan,
        ]
        if any([v == n for n in none_vals]):
            assert v is not None, "None is not allowed when setting metadata"
        return v


class FastqSet(BaseModel):
    name: str = Field(default=None, description="Name of the assay")
    r1: FastqFile
    r2: Optional[FastqFile] = None

    @property
    def fastq_paths(self):
        return [self.r1.path, self.r2.path] if self.is_paired else [self.r1.path]

    @property
    def is_paired(self):
        if self.r2 is None:
            return False
        elif self.r2.path.is_file():
            return True
        else:
            return False

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


class FastqSetIP(FastqSet):
    name: str = Field(default=None, description="Name of the sample set")
    r1: FastqFileIP
    r2: Optional[FastqFileIP] = None

    @property
    def ip_or_control_name(self) -> str:
        return self.r1.ip

    @property
    def sample_name(self) -> str:
        return f"{self.name}_{self.ip_or_control_name}"

    @property
    def is_control(self) -> bool:
        return self.r1.is_control


class IPExperiment(BaseModel):
    ip: FastqSetIP
    control: Optional[FastqSetIP] = None

    @property
    def has_control(self) -> bool:
        return self.control is not None

    @property
    def ip_set_fullname(self) -> str:
        return self.ip.sample_name

    @property
    def control_set_fullname(self) -> str:
        return self.control.sample_name

    @property
    def ip_set_ip_performed(self) -> str:
        return self.ip.ip_or_control_name

    @property
    def control_set_ip_performed(self) -> str:
        return self.control.ip_or_control_name

    @computed_field
    @property
    def fastqs_are_paired(self) -> bool:

        ip = self.ip.is_paired
        control = self.control.is_paired if self.control else True
        return ip and control


class DataFrameDesign(pandera.DataFrameModel):
    sample_name: Series[str]
    r1: Series[str] = pandera.Field(coerce=True)
    r2: Series[str] = pandera.Field(coerce=True)
    scale_group: Series[str]
    deseq2: Optional[Series[str]] = pandera.Field()
    merge: Optional[Series[str]] = pandera.Field()


class DataFrameDesignIP(pandera.DataFrameModel):
    sample_name: Series[str]
    ip: Series[str] = pandera.Field(coerce=True)
    control: Optional[Series[str]] = pandera.Field(coerce=True, nullable=True)
    ip_r1: Series[str] = pandera.Field(coerce=True)
    ip_r2: Series[str] = pandera.Field(coerce=True)
    control_r1: Optional[Series[str]] = pandera.Field(coerce=True, nullable=True)
    control_r2: Optional[Series[str]] = pandera.Field(coerce=True, nullable=True)
    scale_group: Series[str]


class Design(BaseModel):
    fastq_sets: List[FastqSet]
    metadata: List[Metadata]

    @property
    def sample_names(self) -> List[str]:
        return [f.name for f in self.fastq_sets]

    @property
    def fastq_paths(self) -> List[pathlib.Path]:
        paths = []
        for fastq_set in self.fastq_sets:
            paths.append(fastq_set.r1.path)
            if fastq_set.r2 is not None:
                paths.append(fastq_set.r2.path)
        return paths

    def query(self, sample_name: str) -> FastqSet:
        """
        Extract a sample pair of fastq files from the design.
        """

        for fastq_set in self.fastq_sets:
            if fastq_set.name == sample_name:
                return fastq_set

        raise ValueError(f"Could not find sample with name {sample_name}")

    @classmethod
    def from_fastq_files(cls, fq: List[Union[str, pathlib.Path]], **kwargs):
        """
        Generate a Design object from a list of FastqFiles.
        """
        import pandas as pd

        fq = sorted([pathlib.Path(f) for f in fq])

        # Collate the fastq files by sample name
        df = pd.DataFrame(fq, columns=["path"])
        df = df.assign(
            fastq_files=lambda x: x["path"].apply(lambda y: FastqFile(path=y)),
            read_number=lambda x: x["fastq_files"].apply(lambda y: y.read_number),
            sample_stem=lambda x: x["fastq_files"].apply(lambda y: y.stem),
            sample_base=lambda x: x["fastq_files"].apply(lambda y: y.sample_base),
        )

        # Create the fastq sets
        fastq_sets = []
        for sample_name, group in df.groupby("sample_base"):
            if group.shape[0] == 1:
                fq_set = FastqSet(
                    name=sample_name, r1=group["fastq_files"].iloc[0], **kwargs
                )
            elif group.shape[0] == 2:
                fq_set = FastqSet(
                    name=sample_name,
                    r1=group["fastq_files"].iloc[0],
                    r2=group["fastq_files"].iloc[1],
                    **kwargs,
                )
            else:
                raise ValueError(
                    f"Invalid number of fastq files ({group.shape[0]}) for {sample_name}"
                )

            fastq_sets.append(fq_set)

        return cls(
            fastq_sets=fastq_sets,
            metadata=[Metadata(scale_group="all") for _ in fastq_sets],
        )

    def to_dataframe(self) -> pd.DataFrame:
        """
        Return the Design object as a pandas DataFrame.
        """

        data = []
        for fastq_set, metadata in zip(self.fastq_sets, self.metadata):
            row = {
                "sample_name": fastq_set.name,
                "r1": fastq_set.r1.path,
                "r2": fastq_set.r2.path if fastq_set.r2 is not None else None,
            }

            for k, v in metadata.model_dump(exclude_none=True).items():
                row[k] = v

            data.append(row)

        df = pd.DataFrame(data).sort_values("sample_name")

        return DataFrameDesign.validate(df)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, **kwargs):
        """
        Create a Design object from a pandas DataFrame.
        """
        df = DataFrameDesign.validate(df)

        fastq_sets = []
        metadata = []

        non_metadata_keys = ["sample_name", "r1", "r2"]

        for _, row in df.iterrows():
            if row["r2"] is not None:
                fastq_sets.append(
                    FastqSet(
                        name=row["sample_name"],
                        r1=FastqFile(path=row["r1"]),
                        r2=FastqFile(path=row["r2"]),
                        **kwargs,
                    )
                )
            else:
                fastq_sets.append(
                    FastqSet(
                        name=row["sample_name"], r1=FastqFile(path=row["r1"]), **kwargs
                    )
                )

            metadata.append(
                Metadata(**{k: v for k, v in row.items() if k not in non_metadata_keys})
            )

        return cls(fastq_sets=fastq_sets, metadata=metadata, **kwargs)


class DesignIP(BaseModel):
    experiments: List[IPExperiment]
    metadata: List[Metadata]

    @property
    def sample_names_ip(self) -> List[str]:
        return [f.ip_set_fullname for f in self.experiments]

    @property
    def sample_names_control(self) -> List[str]:
        names = set()
        for f in self.experiments:
            if f.has_control:
                names.add(f.control_set_fullname)
        return list(names)

    @property
    def sample_names(self) -> List[str]:
        return sorted([*self.sample_names_ip, *self.sample_names_control])

    @property
    def ips_performed(self) -> List[str]:
        ip = set()
        for f in self.experiments:
            ip.add(f.ip_set_ip_performed)
        return list(ip)

    @property
    def controls_performed(self) -> List[str]:
        control = set()
        for f in self.experiments:
            if f.has_control:
                control.add(f.control_set_ip_performed)
        return list(control)

    def query(self, name: str) -> FastqSetIP:
        """
        Extracts a pair of fastq files from the design.
        """
        ip_names = set(f.ip_set_fullname for f in self.experiments)
        control_names = set(
            f.control_set_fullname for f in self.experiments if f.has_control
        )

        if name in ip_names or name in control_names:
            for experiment in self.experiments:
                if experiment.ip_set_fullname == name:
                    return experiment.ip
                elif experiment.has_control and experiment.control_set_fullname == name:
                    return experiment.control
        else:
            raise ValueError(f"Could not find sample with name {name}")

    @classmethod
    def from_fastq_files(cls, fq: List[Union[str, pathlib.Path]], **kwargs):
        """
        Generate a Design object from a list of FastqFiles.
        """
        import pandas as pd

        fq = sorted([pathlib.Path(f) for f in fq])

        # Collate the fastq files by sample name
        df = pd.DataFrame(fq, columns=["path"])
        df = df.assign(
            fastq_files=lambda x: x["path"].apply(lambda y: FastqFileIP(path=y)),
            read_number=lambda x: x["fastq_files"].apply(lambda y: y.read_number),
            sample_stem=lambda x: x["fastq_files"].apply(lambda y: y.stem),
            sample_base=lambda x: x["fastq_files"].apply(lambda y: y.sample_base),
            sample_base_without_ip=lambda x: x["fastq_files"].apply(
                lambda y: y.sample_base_without_ip
            ),
            is_control=lambda x: x["fastq_files"].apply(lambda y: y.is_control),
        )

        # Break the dataframe into IP and control files
        ip_files = df.query("is_control == False")
        control_files = df.query("is_control == True")[["path", "sample_base_without_ip"]]

        # Merge the IP and control files using the sample base without the IP
        df = ip_files.merge(
            control_files,
            on=["sample_base_without_ip"],
            suffixes=("_ip", "_control"),
            how="left",
        ).assign(
            has_control=lambda x: x["path_control"].notnull(),
        ).drop_duplicates("read_number")

        # Group the files by the sample base
        experiments = []
        for base, group in df.groupby("sample_base"):
            if group.shape[0] == 1:
                # Single end experiment
                ip = FastqSetIP(name=base, r1=FastqFileIP(path=group["path_ip"].iloc[0]))

                if group["path_control"].iloc[0]:
                    control = FastqSetIP(name=base, r1=FastqFileIP(path=group["path_control"].iloc[0]))
                else:
                    control = None

                experiments.append(IPExperiment(ip=ip, control=control, **kwargs))

            elif group.shape[0] == 2:
                # Paired end experiment
                ip = FastqSetIP(
                    name=base,
                    r1=FastqFileIP(path=group["path_ip"].iloc[0]),
                    r2=FastqFileIP(path=group["path_ip"].iloc[1]),
                )

                if group["has_control"].iloc[0]:
                    control = FastqSetIP(name=base, r1=FastqFileIP(path=group["path_control"].iloc[0]), r2=FastqFileIP(path=group["path_control"].iloc[1]))
                else:
                    control = None

                experiments.append(IPExperiment(ip=ip, control=control, **kwargs))

            else:
                raise ValueError(
                    f"Invalid number of fastq files ({group.shape[0]}) for {base}"
                )

        return cls(
            experiments=experiments,
            metadata=[Metadata(scale_group="all") for _ in experiments],
        )

    def to_dataframe(self) -> pd.DataFrame:
        """
        Return the Design object as a pandas DataFrame.
        """

        data = []
        for experiment, metadata in zip(self.experiments, self.metadata):
            row = {
                "sample_name": experiment.ip.sample_name,
                "ip": experiment.ip.ip_or_control_name,
                "control": experiment.control.ip_or_control_name if experiment.control else None,
                "ip_r1": experiment.ip.r1.path,
                "ip_r2": experiment.ip.r2.path if experiment.ip.r2 else None,
                "control_r1": experiment.control.r1.path if experiment.control else None,
                "control_r2": experiment.control.r2.path if experiment.control else None,
            }

            for k, v in metadata.model_dump(exclude_none=True).items():
                row[k] = v

            data.append(row)

        df = pd.DataFrame(data).sort_values("sample_name")

        return DataFrameDesignIP.validate(df)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, **kwargs):
        """
        Create a Design object from a pandas DataFrame.
        """
        df = DataFrameDesignIP.validate(df)

        experiments = []
        metadata = []

        non_metadata_keys = [
            "sample_name",
            "ip",
            "control",
            "ip_r1",
            "ip_r2",
            "control_r1",
            "control_r2",
        ]

        for _, row in df.iterrows():
            ip = FastqSetIP(
                r1=FastqFileIP(path=row["ip_r1"]),
                r2=FastqFileIP(path=row["ip_r2"]) if row["ip_r2"] else None,
            )
            control = (
                FastqSetIP(
                    r1=FastqFileIP(path=row["control_r1"]),
                    r2=FastqFileIP(path=row["control_r2"]) if row["control_r2"] else None,
                )
                if row["control_r1"]
                else None
            )

            experiments.append(IPExperiment(ip=ip, control=control, **kwargs))
            metadata.append(
                Metadata(**{k: v for k, v in row.items() if k not in non_metadata_keys})
            )

        return cls(experiments=experiments, metadata=metadata, **kwargs)


class NormGroup(BaseModel):
    """
    Class to handle normalisation groups.
    """

    group: Optional[Union[str, int]] = "all"
    samples: List[str]
    reference_sample: Optional[str] = None

    @classmethod
    def from_design(
        cls,
        design: Union[Design, DesignIP],
        reference_sample: Optional[str] = None,
        subset_column: Optional[str] = "scale_group",
        subset_value: Optional[List[str]] = None,
    ):
        df = design.to_dataframe()

        if subset_value:
            df = df.query(f"{subset_column} in {subset_value}")

        samples = df.index.tolist()
        reference_sample = reference_sample or df.index[0]
        return cls(
            samples=samples,
            reference_sample=reference_sample,
            group=subset_value[0] or "all",
        )


class NormGroups(BaseModel):
    """
    Class to handle normalisation groups.
    """

    groups: List[NormGroup]

    @classmethod
    def from_design(
        cls,
        design: Union[Design, DesignIP],
        reference_sample: Optional[str] = None,
        subset_column: Optional[str] = "scale_group",
    ):
        df = design.to_dataframe()

        # If the subset column is in the design
        # make the groups based on the subset column
        if subset_column in df.columns:
            subset_values = df[subset_column].drop_duplicates()

            return cls(
                groups=[
                    NormGroup.from_design(
                        design,
                        reference_sample,
                        subset_column,
                        [
                            subset_value,
                        ],
                    )
                    for subset_value in subset_values
                ]
            )

        else:  # If not then just make one group with all the samples
            return cls(
                groups=[
                    NormGroup(
                        group="all",
                        samples=design.sample_names,
                    )
                ]
            )

    @property
    def sample_groups(self) -> Dict[str, List[str]]:
        return {group.group: group.samples for group in self.groups}

    @property
    def group_samples(self) -> Dict[str, List[str]]:
        return {
            sample: group.group for group in self.groups for sample in group.samples
        }

    def get_sample_group(self, sample: str) -> str:
        return self.group_samples[sample]

    def get_grouped_samples(self, group: str) -> List[str]:
        return self.sample_groups[group]


class QCFiles(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]
    fastq_screen: bool = False
    library_complexity: bool = False

    @property
    def default_files(self) -> List[str]:
        return [
            "seqnado_output/qc/fastq_raw_qc.html",
            "seqnado_output/qc/fastq_trimmed_qc.html",
            "seqnado_output/qc/alignment_raw_qc.html",
            "seqnado_output/qc/alignment_filtered_qc.html",
            "seqnado_output/qc/full_qc_report.html",
        ]

    @property
    def fastq_screen_files(self) -> List[str]:
        return ["seqnado_output/qc/full_fastqscreen_report.html"]

    @property
    def library_complexity_files(self) -> List[str]:
        return ["seqnado_output/qc/library_complexity_qc.html"]

    @computed_field
    @property
    def files(self) -> List[str]:
        files = self.default_files
        if self.fastq_screen:
            files.extend(self.fastq_screen_files)
        if self.library_complexity:
            files.extend(self.library_complexity_files)
        return files


class BigWigFiles(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]
    names: List[str]
    pileup_method: Union[
        Literal["deeptools", "homer"], List[Literal["deeptools", "homer"]]
    ] = None
    make_bigwigs: bool = False
    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw", "merged"]] = None
    prefix: Optional[str] = "seqnado_output/bigwigs/"

    def model_post_init(self, __context: Any) -> None:
        if isinstance(self.pileup_method, str):
            self.pileup_method = [self.pileup_method]

    @property
    def bigwigs_non_rna(self):
        scale_methods = (
            ["unscaled", self.scale_method] if self.scale_method else ["unscaled"]
        )

        return expand(
            self.prefix + "{method}/{scale}/{sample}.bigWig",
            sample=self.names,
            scale=scale_methods,
            method=self.pileup_method,
        )

    @property
    def bigwigs_rna(self):

        scale_methods = (
            ["unscaled", self.scale_method] if self.scale_method else ["unscaled"]
        )

        return expand(
            self.prefix + "{method}/{scale}/{sample}_{strand}.bigWig",
            sample=self.names,
            scale=scale_methods,
            method=self.pileup_method,
            strand=["plus", "minus"],
        )

    @computed_field
    @property
    def files(self) -> List[str]:
        if self.make_bigwigs:
            if self.assay == "RNA":
                return self.bigwigs_rna
            else:
                return self.bigwigs_non_rna
        else:
            return []


class PeakCallingFiles(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]
    names: List[str]
    peak_calling_method: Union[
        Literal["macs", "homer", "lanceotron", "seacr", False],
        List[Literal["macs", "homer", "lanceotron", "seacr"]],
    ] = None
    call_peaks: bool = False
    prefix: Optional[str] = "seqnado_output/peaks/"

    @property
    def peak_files(self) -> List[str]:
        return expand(
            self.prefix + "{method}/{sample}.bed",
            sample=self.names,
            method=self.peak_calling_method,
        )

    @computed_field
    @property
    def files(self) -> List[str]:
        if self.call_peaks:
            return self.peak_files
        else:
            return []


class HeatmapFiles(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]

    @property
    def heatmap_files(self) -> List[str]:
        return [
            "seqnado_output/heatmap/heatmap.pdf",
            "seqnado_output/heatmap/metaplot.pdf",
        ]

    @computed_field
    @property
    def files(self) -> List[str]:
        return self.heatmap_files


class HubFiles(BaseModel):
    hub_dir: pathlib.Path
    hub_name: str
    make_ucsc_hub: bool = False

    @computed_field
    @property
    def hub_txt(self) -> pathlib.Path:
        return self.hub_dir / f"{self.hub_name}.hub.txt"

    @computed_field
    @property
    def files(self) -> List[str]:
        if self.make_ucsc_hub:
            return [str(self.hub_txt)]
        else:
            return []

    @classmethod
    def from_dict(cls, d: dict, create_ucsc_hub: bool = False):
        return cls(
            hub_dir=d["directory"],
            hub_name=d["name"],
            make_ucsc_hub=create_ucsc_hub,
        )


class SpikeInFiles(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]
    sample_names: List[str]
    chip_spikein_normalisation: bool = False

    @property
    def norm_factors(self):
        return "seqnado_output/resources/normalisation_factors.tsv"

    @computed_field
    @property
    def files(self) -> List[str]:
        if self.chip_spikein_normalisation:
            return [self.norm_factors]
        else:
            return []


class Output(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]
    run_design: Union[Design, DesignIP]
    sample_names: List[str]

    make_bigwigs: bool = False
    pileup_method: Optional[
        Union[Literal["deeptools", "homer"], List[Literal["deeptools", "homer"]]]
    ] = None
    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw"]] = None

    make_heatmaps: bool = False
    make_ucsc_hub: bool = False

    ucsc_hub_details: Optional[Dict[str, Any]] = None

    fastq_screen: bool = False
    library_complexity: bool = False

    @property
    def merge_bigwigs(self):
        return "merge" in self.run_design.to_dataframe().columns

    @property
    def design_dataframe(self):
        return self.run_design.to_dataframe()

    @property
    def design(self):
        return ["seqnado_output/design.csv"]

    @property
    def bigwigs(self):
        bwf_samples = BigWigFiles(
            assay=self.assay,
            names=self.sample_names,
            make_bigwigs=self.make_bigwigs,
            pileup_method=self.pileup_method,
            scale_method=self.scale_method,
        )
        if self.merge_bigwigs:
            bwf_merged = BigWigFiles(
                assay=self.assay,
                names=self.design_dataframe["merge"].unique().tolist(),
                make_bigwigs=self.make_bigwigs,
                pileup_method="deeptools",
                scale_method="merged",
            )

            files = bwf_samples.files + bwf_merged.files
        else:
            files = bwf_samples.files

        return files or []

    @property
    def heatmaps(self):
        hmf = HeatmapFiles(assay=self.assay, make_heatmaps=self.make_heatmaps)
        return hmf.files

    @property
    def ucsc_hub(self):
        hbf = HubFiles.from_dict(
            self.ucsc_hub_details,
            create_ucsc_hub=self.make_ucsc_hub,
        )
        return hbf

    @property
    def bigbed(self) -> List[str]:
        bb = []
        for peak_file in self.peaks:
            bed = pathlib.Path(peak_file)
            bigbed = bed.with_suffix(".bigBed")
            bb.append(bigbed)
        return bb


class RNAOutput(Output):
    assay: Literal["RNA"]
    project_name: str
    run_deseq2: bool = False
    rna_quantification: Optional[Literal["feature_counts", "salmon"]] = None

    @property
    def counts(self):
        if self.rna_quantification == "feature_counts":
            return ["seqnado_output/readcounts/feature_counts/read_counts.tsv"]
        elif self.rna_quantification == "salmon":
            return ["seqnado_output/readcounts/salmon/salmon_counts.csv"]

    @property
    def deseq2(self):
        if self.run_deseq2:
            return [f"deseq2_{self.project_name}.html"]

    @property
    def peaks(self):
        return []

    @computed_field
    @property
    def files(self) -> List[str]:

        files = []
        files.extend(
            QCFiles(
                assay=self.assay,
                fastq_screen=self.fastq_screen,
                library_complexity=self.library_complexity,
            ).files
        )

        for file_list in (
            self.bigwigs,
            self.heatmaps,
            self.ucsc_hub.files,
            self.counts,
            self.design,
        ):
            if file_list:
                files.extend(file_list)

        if self.run_deseq2:
            files.append(self.deseq2)

        return files


class NonRNAOutput(Output):
    assay: Union[Literal["ChIP"], Literal["ATAC"]]
    call_peaks: bool = False
    peak_calling_method: Optional[
        Union[
            Literal["macs", "homer", "lanceotron", False],
            List[Literal["macs", "homer", "lanceotron"]],
        ]
    ] = None

    @property
    def merge_peaks(self):
        return "merge" in self.design_dataframe.columns

    @property
    def merged_peaks(self):
        return PeakCallingFiles(
            assay=self.assay,
            names=self.design_dataframe["merge"].unique().tolist(),
            call_peaks=self.call_peaks,
            peak_calling_method="lanceotron",
            prefix="seqnado_output/peaks/merged/",
        )

    @computed_field
    @property
    def peaks(self) -> List[str]:
        pcf_samples = PeakCallingFiles(
            assay=self.assay,
            names=self.sample_names,
            call_peaks=self.call_peaks,
            peak_calling_method=self.peak_calling_method,
        )

        if self.merge_peaks:
            pcf_merged = self.merged_peaks

            files = pcf_samples.files + pcf_merged.files
        else:
            files = pcf_samples.files

        return files or []

    @computed_field
    @property
    def files(self) -> List[str]:
        files = []
        files.extend(
            QCFiles(
                assay=self.assay,
                fastq_screen=self.fastq_screen,
                library_complexity=self.library_complexity,
            ).files
        )

        for file_list in (
            self.bigwigs,
            self.heatmaps,
            self.ucsc_hub.files,
            self.peaks,
            self.design,
        ):
            if file_list:
                files.extend(file_list)

        return files


class ATACOutput(NonRNAOutput):
    assay: Literal["ATAC"]


class ChIPOutput(NonRNAOutput):
    assay: Literal["ChIP"]
    ip_names: List[str]
    control_names: List[str]
    call_peaks: bool = False
    peak_calling_method: Optional[
        Union[
            Literal["macs", "homer", "lanceotron", "seacr", False],
            List[Literal["macs", "homer", "lanceotron", "seacr"]],
        ]
    ] = None
    chip_spikein_normalisation: bool = False
    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw"]] = None

    @property
    def peaks(self):
        ip_sample_names = [
            s
            for s in self.sample_names
            if any([c not in s for c in self.control_names])
        ]
        pcf_samples = PeakCallingFiles(
            assay=self.assay,
            names=ip_sample_names,
            call_peaks=self.call_peaks,
            peak_calling_method=self.peak_calling_method,
        )

        if self.merge_peaks:
            pcf_merged = self.merged_peaks
            return pcf_samples.files + pcf_merged.files

        else:
            return pcf_samples.files

    @property
    def spikeins(self):
        sif = SpikeInFiles(
            assay=self.assay,
            sample_names=self.ip_names,
            chip_spikein_normalisation=self.chip_spikein_normalisation,
        )
        return sif.files

    @computed_field
    @property
    def files(self) -> List[str]:
        files = []
        files.extend(
            QCFiles(
                assay=self.assay,
                fastq_screen=self.fastq_screen,
                library_complexity=self.library_complexity,
            ).files
        )

        for file_list in (
            self.bigwigs,
            self.heatmaps,
            self.ucsc_hub.files,
            self.peaks,
            self.spikeins,
            self.design,
        ):
            if file_list:
                files.extend(file_list)

        return files
