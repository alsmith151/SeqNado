import os
import pathlib
import re
import sys
from enum import Enum
from typing import Any, Dict, List, Literal, LiteralString, Optional, Union

import numpy as np
import pandas as pd
import pandera
from loguru import logger
from pandera.typing import DataFrame, Index, Series
from pydantic import BaseModel, Field, computed_field, field_validator, validator
from snakemake.io import expand


def predict_organism(genome: str) -> str:
    if "hg" in genome:
        return "Homo sapiens"
    elif "mm" in genome:
        return "Mus musculus"


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
    use_resolved_name: bool = False

    def model_post_init(self, *args):
        if self.use_resolved_name:
            self.path = pathlib.Path(self.path).resolve()
        else:
            self.path = pathlib.Path(self.path).absolute()

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

    @validator("ip")
    def allow_na_or_nan(cls, v):
        if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
            return v
        if not isinstance(v, str):
            raise ValueError("ip must be a string, None, pd.NA, or np.nan")
        return v


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
    ip: str = Field(default=None, description="IP performed on the sample")

    @property
    def ip_or_control_name(self) -> str:
        return self.ip if self.ip else self.r1.ip

    @property
    def sample_name(self) -> str:
        return f"{self.name}_{self.ip_or_control_name}"

    @property
    def sample_name_base(self) -> str:
        return self.r1.sample_base_without_ip

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
    def control_fullname(self) -> str:
        return self.control.sample_name

    @property
    def ip_performed(self) -> str:
        return self.ip.ip_or_control_name

    @property
    def control_performed(self) -> str:
        return self.control.ip_or_control_name

    @property
    def fastqs_are_paired(self) -> bool:
        ip = self.ip.is_paired
        control = self.control.is_paired if self.control else True
        return ip and control


class DataFrameDesign(pandera.DataFrameModel):
    sample_name: Series[str]
    r1: Series[str] = pandera.Field(coerce=True)
    r2: Series[str] = pandera.Field(coerce=True, nullable=True)
    scale_group: Series[str]
    deseq2: Optional[Series[str]] = pandera.Field()
    merge: Optional[Series[str]] = pandera.Field()


class DataFrameDesignIP(pandera.DataFrameModel):
    sample_name: Series[str]
    ip: Series[str] = pandera.Field(coerce=True)
    control: Optional[Series[str]] = pandera.Field(coerce=True, nullable=True)
    ip_r1: Series[str] = pandera.Field(coerce=True)
    ip_r2: Series[str] = pandera.Field(coerce=True, nullable=True)
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

    @classmethod
    def from_directory(cls, directory: Union[pathlib.Path, str], **kwargs):
        """
        Create a Design object from a directory of fastq files.
        """
        directory = pathlib.Path(directory)
        file_patterns = ["*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"]
        fastq_files = sorted(
            [f for pattern in file_patterns for f in directory.glob(pattern)]
        )

        if len(fastq_files) == 0:
            raise FileNotFoundError(f"No fastq files found in {directory}")

        return cls.from_fastq_files(fastq_files, **kwargs)

    def to_geo_dataframe(
        self, assay: Literal["ATAC", "RNA", "SNP"], pipeline_config: dict
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame with the GEO metadata.
        """
        geo_samples = []
        for sample_row in self.to_dataframe().itertuples():
            if assay != "RNA":
                processed_data_file = [f"{sample_row.sample_name}.bw"]

            else:
                processed_data_file = [
                    "read_counts.tsv",
                    f"{sample_row.sample_name}_plus.bw",
                    f"{sample_row.sample_name}_minus.bw",
                ]

            sample = GEOSample(
                assay=assay,
                library_name=f"{sample_row.sample_name}",
                title=f"{sample_row.sample_name}",
                organism=predict_organism(pipeline_config["genome"]["name"]),
                cell_line=None,
                cell_type=None,
                antibody=None,
                genotype=None,
                treatment=sample_row.treatment
                if hasattr(sample_row, "treatment")
                else None,
                time=sample_row.time if hasattr(sample_row, "time") else None,
                single_or_paired="paired-end" if sample_row.r2 else "single",
                instrument_model=pipeline_config.get(
                    "instrument_model", "Illumina NovaSeq X"
                ),
                description=None,
                processed_data_file=processed_data_file,
                raw_file=[
                    pathlib.Path(sample_row.r1).name,
                    pathlib.Path(sample_row.r2).name,
                ]
                if sample_row.r2
                else [pathlib.Path(sample_row.r1).name],
            )

            geo_samples.append(sample)

        df_samples = GEOSamples(samples=geo_samples).to_dataframe()
        return df_samples


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
                names.add(f.control_fullname)
        return list(names)

    @property
    def sample_names(self) -> List[str]:
        return sorted([*self.sample_names_ip, *self.sample_names_control])

    @property
    def ips_performed(self) -> List[str]:
        ip = set()
        for f in self.experiments:
            ip.add(f.ip_performed)
        return list(ip)

    @property
    def controls_performed(self) -> List[str]:
        control = set()
        for f in self.experiments:
            if f.has_control:
                control.add(f.control_performed)
        return list(control)

    def query(
        self, sample_name: str, full_experiment: bool = False
    ) -> Union[FastqSetIP, Dict[str, FastqSetIP]]:
        """
        Extracts a pair of fastq files from the design.
        """
        ip_names = set(f.ip_set_fullname for f in self.experiments)
        control_names = set(
            f.control_fullname for f in self.experiments if f.has_control
        )
        is_control = False

        experiment_files = dict()
        is_control = False

        experiment_files = dict()

        if sample_name in ip_names or sample_name in control_names:
            for experiment in self.experiments:
                if experiment.ip_set_fullname == sample_name:
                    experiment_files["ip"] = experiment.ip
                    experiment_files["control"] = experiment.control

                    experiment_files["ip"] = experiment.ip
                    experiment_files["control"] = experiment.control

                elif (
                    experiment.has_control
                    and experiment.control_fullname == sample_name
                ):
                    is_control = True
                    experiment_files["ip"] = experiment.ip
                    experiment_files["control"] = experiment.control
                    is_control = True
                    experiment_files["ip"] = experiment.ip
                    experiment_files["control"] = experiment.control
        else:
            raise ValueError(f"Could not find sample with name {sample_name}")

        if full_experiment:
            return experiment_files
        else:
            return (
                experiment_files["ip"]
                if not is_control
                else experiment_files["control"]
            )

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
        control_files = df.query("is_control == True")[
            ["path", "sample_base_without_ip"]
        ]

        # Merge the IP and control files using the sample base without the IP
        df = ip_files.merge(
            control_files,
            on=["sample_base_without_ip"],
            suffixes=("_ip", "_control"),
            how="outer",
        ).assign(
            has_control=lambda x: x["path_control"].notnull(),
        )

        # Group the files by the sample base
        experiments = []
        for base, group in df.groupby("sample_base"):
            name_without_ip = group["sample_base_without_ip"].iloc[0]

            match group.shape[0]:
                case 1:
                    # Single end experiment no control
                    ip = FastqSetIP(
                        name=name_without_ip,
                        r1=FastqFileIP(path=group["path_ip"].iloc[0]),
                    )
                    control = None

                case 2:
                    # Paired end experiment no control
                    ip = FastqSetIP(
                        name=name_without_ip,
                        r1=FastqFileIP(path=group["path_ip"].iloc[0]),
                        r2=FastqFileIP(path=group["path_ip"].iloc[1]),
                    )
                    control = None

                case 4:
                    # Paired end experiment with control

                    # | path_ip | path_control |
                    # |---------|--------------|
                    # | r1      | r1           |
                    # | r1      | r2           |
                    # | r2      | r1           |
                    # | r2      | r2           |

                    ip = FastqSetIP(
                        name=name_without_ip,
                        r1=FastqFileIP(path=group["path_ip"].iloc[0]),
                        r2=FastqFileIP(path=group["path_ip"].iloc[2]),
                    )
                    control = FastqSetIP(
                        name=name_without_ip,
                        r1=FastqFileIP(path=group["path_control"].iloc[0]),
                        r2=FastqFileIP(path=group["path_control"].iloc[1]),
                    )

                case _:
                    raise ValueError(
                        f"Invalid number of fastq files ({group.shape[0]}) for {name_without_ip}"
                    )

            experiments.append(IPExperiment(ip=ip, control=control, **kwargs))

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
                "sample_name": experiment.ip.name,
                "ip": experiment.ip.ip_or_control_name,
                "control": (
                    experiment.control.ip_or_control_name
                    if experiment.control
                    else None
                ),
                "ip_r1": experiment.ip.r1.path,
                "ip_r2": experiment.ip.r2.path if experiment.ip.r2 else None,
                "control_r1": (
                    experiment.control.r1.path if experiment.control else None
                ),
                "control_r2": (
                    experiment.control.r2.path
                    if experiment.control and experiment.control.r2
                    else None
                ),
            }
            for k, v in metadata.model_dump(exclude_none=True).items():
                row[k] = v

            data.append(row)

        df = pd.DataFrame(data).sort_values("sample_name")

        return DataFrameDesignIP.validate(df)

    def to_geo_dataframe(
        self, assay: Literal["ChIP", "CAT"], pipeline_config: dict
    ) -> pd.DataFrame:
        """
        Create a pandas DataFrame with the GEO metadata.
        """
        geo_samples = []
        for sample_row in self.to_dataframe().itertuples():
            for ip_type in ["ip", "control"]:
                processed_data_file = [
                    f"{sample_row.sample_name}_{getattr(sample_row, ip_type)}.bigWig"
                ]

                raw_files = [
                    pathlib.Path(getattr(sample_row, f"{ip_type}_r1")).name,
                    pathlib.Path(getattr(sample_row, f"{ip_type}_r2")).name,
                ]

                sample = GEOSample(
                    assay=assay,
                    library_name=f"{sample_row.sample_name}",
                    title=f"{sample_row.sample_name}",
                    organism=predict_organism(pipeline_config["genome"]["name"]),
                    cell_line=None,
                    cell_type=None,
                    antibody=getattr(sample_row, ip_type),
                    genotype=None,
                    treatment=sample_row.treatment
                    if hasattr(sample_row, "treatment")
                    else None,
                    time=sample_row.time if hasattr(sample_row, "time") else None,
                    single_or_paired="paired-end" if sample_row.r2 else "single",
                    instrument_model=pipeline_config.get(
                        "instrument_model", "Illumina NovaSeq X"
                    ),
                    description=None,
                    processed_data_file=processed_data_file,
                    raw_file=raw_files,
                )

                geo_samples.append(sample)

        df_samples = GEOSamples(samples=geo_samples).to_dataframe()
        return df_samples

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, **kwargs):
        """
        Create a Design object from a pandas DataFrame.
        """

        df = df.fillna("")
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
                name=row["sample_name"],
                ip=row["ip"],
                r1=FastqFileIP(path=row["ip_r1"]),
                r2=FastqFileIP(path=row["ip_r2"]) if row["ip_r2"] else None,
            )
            control = (
                FastqSetIP(
                    name=row["sample_name"],
                    ip=row["control"],
                    r1=FastqFileIP(path=row["control_r1"])
                    if row["control_r1"]
                    else None,
                    r2=(
                        FastqFileIP(path=row["control_r2"])
                        if row["control_r2"]
                        else None
                    ),
                )
                if row["control_r1"]
                else None
            )

            experiments.append(IPExperiment(ip=ip, control=control, **kwargs))
            metadata.append(
                Metadata(**{k: v for k, v in row.items() if k not in non_metadata_keys})
            )
        return cls(experiments=experiments, metadata=metadata, **kwargs)

    @classmethod
    def from_directory(cls, directory: Union[pathlib.Path, str], **kwargs):
        """
        Create a Design object from a directory of fastq files.
        """
        directory = pathlib.Path(directory)
        file_patterns = ["*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"]
        fastq_files = sorted(
            [f for pattern in file_patterns for f in directory.glob(pattern)]
        )

        if len(fastq_files) == 0:
            raise FileNotFoundError(f"No fastq files found in {directory}")

        return cls.from_fastq_files(fastq_files, **kwargs)

    @property
    def fastq_paths(self) -> List[pathlib.Path]:
        paths = []
        for experiment in self.experiments:
            paths.extend(experiment.ip.fastq_paths)
            if experiment.control:
                paths.extend(experiment.control.fastq_paths)
        return paths


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
        include_controls: bool = False,
    ):
        if isinstance(design, Design):
            df = (
                design.to_dataframe()
                .assign(sample_fullname=lambda df: df.sample_name)
                .set_index("sample_fullname")
            )
        elif isinstance(design, DesignIP) and not include_controls:
            df = (
                design.to_dataframe()
                .assign(sample_fullname=lambda df: df.sample_name + "_" + df.ip)
                .set_index("sample_fullname")
            )
        elif isinstance(design, DesignIP) and include_controls:
            df_ip = (
                design.to_dataframe()
                .assign(sample_fullname=lambda df: df.sample_name + "_" + df.ip)
                .set_index("sample_fullname")
            )
            df_control = (
                design.to_dataframe()
                .query("control.notnull()")
                .assign(sample_fullname=lambda df: df.sample_name + "_" + df.control)
                .set_index("sample_fullname")
            )
            df = pd.concat([df_ip, df_control])

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
        include_controls: bool = False,
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
                        include_controls,
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


def generate_fastq_raw_names(
    sample_name: str, is_paired: bool = True
) -> Dict[str, List[str]]:
    """
    Get the fastq files for a sample.
    """

    if is_paired:
        exts = ["_1.fastq.gz", "_2.fastq.gz"]
    else:
        exts = [".fastq.gz"]
    fq = [f"{sample_name}{ext}" for ext in exts]
    return {sample_name: fq}


class GEOFiles(BaseModel):
    make_geo_submission_files: bool

    assay: Literal["ChIP", "ATAC", "RNA", "SNP", "CUT&TAG", "METH"]
    sample_names: List[str]
    config: dict
    design: pd.DataFrame
    extensions_allowed: List[str] = [".txt", ".bigWig", ".bed", ".tsv", ".vcf.gz"]
    processed_files: Optional[List[Union[str, pathlib.Path]]] = None

    class Config:
        arbitrary_types_allowed = True

    @property
    def md5sums(self):
        return [
            "seqnado_output/geo_submission/md5sums.txt",
            "seqnado_output/geo_submission/raw_data_checksums.txt",
            "seqnado_output/geo_submission/processed_data_checksums.txt",
            "seqnado_output/geo_submission/samples_table.txt",
            "seqnado_output/geo_submission/protocol.txt",
        ]

    @property
    def upload_directory(self):
        return pathlib.Path("seqnado_output/geo_submission") / self.assay.replace(
            "&", "_and_"
        )  # Fix for CUT&TAG

    @property
    def upload_instructions(self):
        return pathlib.Path("seqnado_output/geo_submission") / "upload_instructions.txt"

    @property
    def processed_data_files(self) -> pd.DataFrame:
        wanted_exts = self.extensions_allowed
        unwanted_files = [*self.md5sums]

        # Create a DataFrame with the processed files
        df = pd.Series([pathlib.Path(p) for p in self.processed_files]).to_frame("path")
        df = df.assign(
            ext=lambda x: x["path"].apply(lambda x: x.suffix),
        )

        # Filter the DataFrame to only include the wanted files
        df = df.query("ext in @wanted_exts")
        df = df.query("~path.astype('str').isin(@unwanted_files)")
        df = df.query('~path.astype("str").str.contains("hub.txt")')

        # Add the sample name, method, and normalisation columns
        df = df.assign(name=lambda x: x["path"].apply(lambda x: x.stem))
        df = df.assign(
            normalisation=lambda x: np.where(
                x["ext"] != ".bed", x["path"].apply(lambda x: x.parts[-2]), ""
            )
        )
        df = df.assign(
            method=lambda x: np.where(
                x["ext"] == ".bigWig",
                x["path"].apply(lambda x: x.parts[-3]),
                x["path"].apply(lambda x: x.parts[-2]),
            ),
        )
        df = df.assign(
            normalisation=lambda df: df.normalisation.str.replace(
                "unscaled", ""
            ).str.replace("spikein", "reference-normalised")
        )
        df = df.sort_values(by=["name", "ext", "method", "normalisation"])

        # Add the output file name and file type columns
        df = df.assign(
            output_file_name=lambda x: (
                x["name"] + "_" + x["method"] + "_" + x["normalisation"] + x["ext"]
            ).str.replace("_.", "."),
            file_type=lambda df: np.select(
                [
                    df["ext"] == ".bigWig",
                    df["ext"] == ".bed",
                    df["ext"] == ".tsv",
                    df["ext"] == ".vcf.gz",
                ],
                ["signal", "peaks", "counts", "variants"],
                default="other",
            ),
        )

        df = df[df.file_type != "other"]

        return df

    @property
    def processed_data_per_sample(self) -> Dict[str, List[str]]:
        return (
            self.processed_data_files.groupby("name")["output_file_name"]
            .apply(list)
            .to_dict()
        )

    @property
    def raw_files(self) -> Dict[str, List[str]]:
        fastq = dict()
        for row in self.design.itertuples():
            if not self.assay in ["ChIP", "CUT&TAG"]:
                sample_name = (
                    row.sample_name
                    # if self.assay not in ["ChIP", "CUT&TAG"]
                    # else f"{row.sample_name}_{row.ip}"
                )

                is_paired = False
                if hasattr(row, "r2") and row.r2:
                    is_paired = True
                # elif hasattr(row, "ip_r2") and row.ip_r2:
                #     is_paired = True

                fqs = generate_fastq_raw_names(sample_name, is_paired)
                fastq.update(fqs)

            else:
                has_control = hasattr(row, "control") and row.control
                sample_name = f"{row.sample_name}_{row.ip}"
                is_paired = hasattr(row, "ip_r2") and row.ip_r2
                fqs = generate_fastq_raw_names(sample_name, is_paired)
                fastq.update(fqs)

                if has_control:
                    control_name = f"{row.sample_name}_{row.control}"
                    is_paired = hasattr(row, "control_r2") and row.control_r2
                    fqs = generate_fastq_raw_names(control_name, is_paired)
                    fastq.update(fqs)
        return fastq

    @property
    def metadata(self):
        """
        Create a pandas DataFrame with the GEO metadata.
        """
        geo_samples = []

        organism = predict_organism(self.config["genome"]["name"])
        instrument_model = self.config.get("instrument_model", "Illumina NovaSeq X")

        for sample_row in self.design.itertuples():
            if self.assay not in ["ChIP", "CUT&TAG"]:
                sample_name = sample_row.sample_name
                library_name = sample_name
                title = sample_name
                antibody = None
                is_paired_end = hasattr(sample_row, "r2")

                geo_sample = GEOSample(
                    assay=self.assay,
                    library_name=library_name,
                    title=title,
                    organism=organism,
                    cell_line=None,
                    cell_type=None,
                    antibody=antibody,
                    genotype=None,
                    treatment=sample_row.treatment
                    if hasattr(sample_row, "treatment")
                    else None,
                    time=sample_row.time if hasattr(sample_row, "time") else None,
                    single_or_paired="paired-end" if is_paired_end else "single",
                    instrument_model=instrument_model,
                    description=None,
                    raw_file=[
                        pathlib.Path(p).name for p in self.raw_files[library_name]
                    ],
                    processed_data_file=[
                        str(p)
                        for p in self.processed_data_per_sample.get(library_name, [])
                    ],
                )

                geo_samples.append(geo_sample)

            elif self.assay in ["ChIP", "CUT&TAG"] and not sample_row.control:
                sample_name = f"{sample_row.sample_name}_{sample_row.ip}"
                library_name = sample_name
                title = sample_name
                is_paired_end = hasattr(sample_row, "ip_r2") and sample_row.ip_r2
                antibody = sample_row.ip

                geo_sample = GEOSample(
                    assay=self.assay,
                    library_name=library_name,
                    title=title,
                    organism=organism,
                    cell_line=None,
                    cell_type=None,
                    antibody=antibody,
                    genotype=None,
                    treatment=sample_row.treatment
                    if hasattr(sample_row, "treatment")
                    else None,
                    time=sample_row.time if hasattr(sample_row, "time") else None,
                    single_or_paired="paired-end" if is_paired_end else "single",
                    instrument_model=instrument_model,
                    description=None,
                    raw_file=[
                        pathlib.Path(p).name for p in self.raw_files[library_name]
                    ],
                    processed_data_file=[
                        str(p)
                        for p in self.processed_data_per_sample.get(library_name, [])
                    ],
                )

                geo_samples.append(geo_sample)

            elif self.assay in ["ChIP", "CUT&TAG"] and sample_row.control:
                for ip_type in ["ip", "control"]:
                    sample_name = (
                        f"{sample_row.sample_name}_{getattr(sample_row, ip_type)}"
                    )
                    library_name = sample_name
                    title = sample_name
                    antibody = getattr(sample_row, ip_type)
                    is_paired_end = hasattr(sample_row, f"{ip_type}_r2") and getattr(
                        sample_row, f"{ip_type}_r2"
                    )

                    geo_sample = GEOSample(
                        assay=self.assay,
                        library_name=library_name,
                        title=title,
                        organism=organism,
                        cell_line=None,
                        cell_type=None,
                        antibody=antibody,
                        genotype=None,
                        treatment=sample_row.treatment
                        if hasattr(sample_row, "treatment")
                        else None,
                        time=sample_row.time if hasattr(sample_row, "time") else None,
                        single_or_paired="paired-end" if is_paired_end else "single",
                        instrument_model=instrument_model,
                        description=None,
                        raw_file=[
                            pathlib.Path(p).name for p in self.raw_files[library_name]
                        ],
                        processed_data_file=[
                            str(p)
                            for p in self.processed_data_per_sample.get(
                                library_name, []
                            )
                        ],
                    )

                    geo_samples.append(geo_sample)

        df_samples = GEOSamples(samples=geo_samples).to_dataframe()
        return df_samples

    @property
    def files(self) -> List[str]:
        if self.make_geo_submission_files:
            return [
                *self.md5sums,
                self.upload_directory,
                self.upload_instructions,
                "seqnado_output/geo_submission/.validated",
            ]
        else:
            return []


class QCFiles(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP", "CUT&TAG", "METH"]
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
    assay: Literal["ChIP", "ATAC", "RNA", "SNP", "CUT&TAG", "METH"]
    names: List[str]
    pileup_method: Union[
        Literal["deeptools", "homer", False],
        List[
            Literal[
                "deeptools",
                "homer",
            ]
        ],
        Literal["deeptools", "homer", False],
        List[
            Literal[
                "deeptools",
                "homer",
            ]
        ],
    ] = None
    make_bigwigs: bool = False
    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw", "merged"]] = None
    include_unscaled: bool = True
    prefix: Optional[str] = "seqnado_output/bigwigs/"

    def model_post_init(self, __context: Any) -> None:
        if isinstance(self.pileup_method, str):
            self.pileup_method = [self.pileup_method]

        if self.include_unscaled and not self.scale_method:
            self.scale_method = [
                "unscaled",
            ]
        elif self.include_unscaled and self.scale_method:
            self.scale_method = ["unscaled", self.scale_method]
        else:
            self.scale_method = [self.scale_method]

    @property
    def bigwigs_non_rna(self):
        return expand(
            self.prefix + "{method}/{scale}/{sample}.bigWig",
            sample=self.names,
            scale=self.scale_method,
            method=self.pileup_method,
        )

    @property
    def bigwigs_rna(self):
        return expand(
            self.prefix + "{method}/{scale}/{sample}_{strand}.bigWig",
            sample=self.names,
            scale=self.scale_method,
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
    assay: Literal["ChIP", "ATAC", "CUT&TAG"]
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
    assay: Literal["ChIP", "ATAC", "RNA", "CUT&TAG"]
    make_heatmaps: bool = False
    make_heatmaps: bool = False

    @property
    def heatmap_files(self) -> List[str]:
        return [
            "seqnado_output/heatmap/heatmap.pdf",
            "seqnado_output/heatmap/metaplot.pdf",
        ]

    @computed_field
    @property
    def files(self) -> List[str]:
        if self.make_heatmaps:
            return self.heatmap_files
        else:
            return []


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
    assay: Literal["ChIP", "ATAC", "RNA", "CUT&TAG"]
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


class PlotFiles(BaseModel):
    perform_plotting: bool = False
    plotting_coordinates: Optional[Union[str, pathlib.Path]] = None
    plotting_format: Literal["svg", "png", "pdf"] = "svg"

    def get_plot_names(self):
        import pyranges as pr

        plots = []

        try:
            coords = pr.read_bed(str(self.plotting_coordinates))
            outdir = pathlib.Path("seqnado_output/genome_browser_plots/")
            for region in coords.df.itertuples():
                fig_name = (
                    f"{region.Chromosome}-{region.Start}-{region.End}"
                    if not hasattr(region, "Name") and not region.Name
                    else region.Name
                )
                plots.append(outdir / f"{fig_name}.{self.plotting_format}")
        
        except FileNotFoundError:
            logger.warning(
                f"Could not find plotting coordinates file: {self.plotting_coordinates}"
            )

        return plots

    @property
    def files(self) -> List[str]:
        if self.plotting_coordinates:
            return self.get_plot_names()
        else:
            return []


class Output(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP", "CUT&TAG"]
    config: dict
    run_design: Union[Design, DesignIP]
    sample_names: List[str]

    make_bigwigs: bool = False
    pileup_method: Union[
        Literal["deeptools", "homer", False],
        List[Literal["deeptools", "homer"]],
    ] = None

    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw"]] = None

    make_heatmaps: bool = False
    make_ucsc_hub: bool = False

    ucsc_hub_details: Optional[Dict[str, Any]] = None

    fastq_screen: bool = False
    library_complexity: bool = False

    geo_submission_files: bool = False

    perform_plotting: bool = False
    plotting_format: Literal["svg", "png", "pdf"] = "svg"
    plotting_coordinates: Optional[Union[str, pathlib.Path]] = None

    # Correct plotting_coordinates type as it may be False
    @validator("plotting_coordinates", pre=True)
    def validate_plotting_coordinates(cls, v):
        if v is False:
            return None
        return v

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
                include_unscaled=False,
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

    @property
    def plots(self):
        if self.perform_plotting:
            pf = PlotFiles(
                plotting_coordinates=self.plotting_coordinates, plotting_format=self.plotting_format
            )
            return pf.files
        else:
            return []

    @property
    def geo_files(self):
        return GEOFiles(
            make_geo_submission_files=self.geo_submission_files,
            assay=self.assay,
            sample_names=self.sample_names,
            design=self.design_dataframe,
            config=self.config,
        )


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

        files.extend(self.geo_files.files)

        for file_list in (
            self.bigwigs,
            self.heatmaps,
            self.ucsc_hub.files,
            self.counts,
            self.design,
            self.plots,
        ):
            if file_list:
                files.extend(file_list)

        if self.run_deseq2:
            files.append(self.deseq2)

        return files


class NonRNAOutput(Output):
    assay: Literal["ChIP", "ATAC", "CUT&TAG"]
    consensus_counts: bool = False
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

    @property
    def merged_counts(self):
        # Get the merged counts file if peaks are merged and consensus counts are requested
        if self.merge_peaks and self.consensus_counts:
            groups = self.merged_peaks.names
            count_files = [
                f"seqnado_output/readcounts/featurecounts/{group}_counts.tsv"
                for group in groups
            ]
            return count_files
        else:
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

        files.extend(self.geo_files.files)

        for file_list in (
            self.bigwigs,
            self.heatmaps,
            self.ucsc_hub.files,
            self.peaks,
            self.design,
            self.plots,
            self.merged_counts,
        ):
            if file_list:
                files.extend(file_list)

        return files


class ATACOutput(NonRNAOutput):
    assay: Literal["ATAC"]


class IPOutput(NonRNAOutput):
    assay: Literal["ChIP", "CUT&TAG"]
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
            if not any([c in s for c in self.control_names])
            if not any([c in s for c in self.control_names])
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

    @property
    def files(self) -> List[str]:
        files = super().files

        for file_list in (
            self.bigwigs,
            self.heatmaps,
            self.ucsc_hub.files,
            self.peaks,
            self.spikeins,
            self.design,
            self.plots,
            self.merged_counts,
        ):
            if file_list:
                files.extend(file_list)

        return list(set(files))


class SNPOutput(Output):
    assay: Literal["SNP"]
    call_snps: bool = False
    sample_names: List[str]
    make_ucsc_hub: bool = False
    snp_calling_method: Optional[
        Union[
            Literal["bcftools", "deepvariant", False],
            List[Literal["bcftools", "deepvariant"]],
        ]
    ] = None

    @property
    def design(self):
        return ["seqnado_output/design.csv"]

    @property
    def snp_files(self) -> List[str]:
        if self.call_snps:
            return expand(
                "seqnado_output/variant/{method}/{sample}.vcf.gz",
                sample=self.sample_names,
                method=self.snp_calling_method,
            )
        else:
            return []

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
            self.snp_files,
            self.design,
        ):
            if file_list:
                files.extend(file_list)

        if self.call_snps:
            files.append(self.snp_files)

        return files


class METHOutput(Output):
    assay: Literal["METH"]
    call_methylation: bool = False
    sample_names: List[str]
    config: dict
    genomes: List[str]
    make_ucsc_hub: bool = False
    methylation_assay: Optional[
        Union[
            Literal["bisulfite", "taps"],
            List[Literal["bisulfite", "taps"]],
        ]
    ] = None

    @property
    def design(self):
        return ["seqnado_output/design.csv"]

    @property
    def meth_split_bams(self) -> List[str]:
        if self.call_methylation:
            return expand(
                "seqnado_output/aligned/spikein/{sample}_{genome}.bam",
                sample=self.sample_names,
                genome=self.genomes,
            )
        return []
    
    @property
    def meth_files(self) -> List[str]:
        if self.call_methylation and "taps" not in self.methylation_assay:
            return expand(
                "seqnado_output/methylation/methyldackel/{sample}_{genome}_CpG.bedGraph",
                sample=self.sample_names,
                genome=self.genomes,
            )
        return []

    @property
    def taps_files(self) -> List[str]:
        if self.call_methylation and "taps" in self.methylation_assay:
            return expand(
                "seqnado_output/methylation/methyldackel/{sample}_{genome}_CpG_TAPS.bedGraph",
                sample=self.sample_names,
                genome=self.genomes,
            )
        return []

    @property
    def methylation_bias(self) -> List[str]:
        """
        Get the methylation bias files. and seqnado_output/methylation/methylation_conversion.tsv"""
        if self.call_methylation:
            return expand(
                "seqnado_output/methylation/methyldackel/bias/{sample}_{genome}.txt",
                sample=self.sample_names,
                genome=self.genomes,
            ), "seqnado_output/methylation/methylation_conversion.tsv"

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
            self.meth_files,
            self.taps_files,
            self.methylation_bias,
            self.design,
        ):
            if file_list:
                files.extend(file_list)

        return files


class MCCOutput(Output):
    assay: Literal["MCC"]
    sample_names: List[str]
    
    viewpoint_oligos: List[str]
    viewpoints_grouped: List[str]

    config: dict
    make_ucsc_hub: bool = False

    resolutions: List[int] = [100]

    @property
    def design(self):
        return ["seqnado_output/design.csv"]

    @property
    def cooler_files(self) -> List[str]:
        return expand(
            "seqnado_output/mcc/{sample}/{viewpoint}.mcool",
            sample=self.sample_names,
            viewpoint=self.viewpoints_grouped,
        )


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
            self.cooler_files,
            self.design,
        ):
            if file_list:
                files.extend(file_list)

        return files



class Molecule(Enum):
    rna_total = "total RNA"
    rna_polya = "polyA RNA"
    rna_cytoplasmic = "cytoplasmic RNA"
    rna_nuclear = "nuclear RNA"
    dna_genomic = "genomic DNA"
    protein = "protein"
    other = "other"


class GEOSample(BaseModel):
    assay: Literal["ATAC", "ChIP", "CUT&TAG", "RNA", "SNP"]
    library_name: str
    title: str
    organism: Literal["Homo sapiens", "Mus musculus"]
    cell_line: Optional[str] = None
    cell_type: Optional[str] = None
    antibody: Optional[str] = None
    genotype: Optional[str] = None
    treatment: Optional[str] = None
    time: Optional[str] = None
    single_or_paired: Literal["single", "paired-end"]
    instrument_model: str
    description: Optional[str] = None
    processed_data_file: List[str]
    raw_file: List[str]

    @computed_field
    @property
    def molecule(self) -> Molecule:
        if self.assay == "ATAC":
            return Molecule.dna_genomic
        elif self.assay == "RNA" and any(
            n in self.title.lower() for n in ["tt-seq", "nasc", "point"]
        ):
            return Molecule.rna_nuclear
        elif self.assay == "RNA":
            return Molecule.rna_polya
        elif self.assay == "SNP":
            return Molecule.dna_genomic
        else:
            return Molecule.dna_genomic

    @property
    def to_series(self):
        data = {
            "library name": self.library_name,
            "title": self.title,
            "organism": self.organism,
            "cell line": self.cell_line,
            "cell type": self.cell_type,
            "ChIP antibody": self.antibody,
            "molecule": self.molecule.value,
            "single or paired-end": self.single_or_paired,
            "instrument model": self.instrument_model,
            "description": self.description,
        }

        processed_data = {
            f"processed data file {i}": f
            for i, f in enumerate(self.processed_data_file)
        }

        raw_data = {f"raw file {i}": f for i, f in enumerate(self.raw_file)}

        data.update(processed_data)
        data.update(raw_data)

        if self.assay in ["ATAC", "RNA", "SNP"]:
            del data["ChIP antibody"]

        return pd.Series(data)


class GEOSamples(BaseModel):
    samples: List[GEOSample]

    def to_dataframe(self):
        df = pd.concat([s.to_series for s in self.samples], axis=1).T
        df.columns = df.columns.str.replace(r"\s\d+$", "", regex=True).str.strip()
        return df
