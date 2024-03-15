import pathlib
import re
from typing import Any, Dict, List, Optional, Union, Literal

import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, computed_field
from snakemake.io import expand


class FastqFile(BaseModel):
    path: pathlib.Path

    def model_post_init(self, *args):
        self.path = pathlib.Path(self.path).resolve()

        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} does not exist.")

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

    @property
    def is_control(self) -> bool:
        return self.r1.is_control


class ExperimentIP(BaseModel):
    ip_files: AssayIP
    control_files: AssayIP = None
    name: str = Field(default=None, description="Name of the assay")
    ip: str = Field(default=None, description="IP performed on the sample")
    control: Optional[str] = Field(
        default=None, description="IP performed on the control"
    )
    metadata: Optional[dict] = None

    def model_post_init(self, *args):
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
        Create a Experiment object from a list of FastqFiles.

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

        if all([s is None for s in sample_names]):
            return []
        else:
            return list(sample_names)

    @property
    def sample_names(self) -> List[str]:
        return self.sample_names_ip + self.sample_names_control

    @property
    def ip_names(self) -> List[str]:
        return list(set([experiment.ip for experiment in self.assays.values()]))

    @property
    def control_names(self) -> List[str]:
        names = list(set([experiment.control for experiment in self.assays.values()]))
        if all([s is None for s in names]):
            return []
        else:
            return names

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

        ## Run through the list and pair upt the read1 and read2 files
        ## If there is only one file, then it is the read1 file
        import itertools

        # Collate the fastq files by sample name
        fq = sorted(fq)
        fastq_collated = dict()
        for f in fq:
            if f.sample_base not in fastq_collated:
                fastq_collated[f.sample_base] = dict()
                fastq_collated[f.sample_base][f.read_number or 1] = f
            else:
                fastq_collated[f.sample_base][f.read_number] = f

        # Create the assays
        assays = {}
        for sample_name, fastq_files in fastq_collated.items():
            assays[sample_name] = AssayIP(
                name=sample_name, r1=fastq_files[1], r2=fastq_files.get(2), **kwargs
            )

        # Create the experiments
        experiments = {}

        for base, assay in itertools.groupby(
            assays.values(), lambda x: x.r1.sample_base_without_ip
        ):
            assay = list(assay)

            if len(assay) == 1:
                experiments[assay[0].name] = ExperimentIP(ip_files=assay[0], **kwargs)
            elif len(assay) == 2 and any([a.is_control for a in assay]):
                ip = [a for a in assay if not a.is_control][0]
                control = [a for a in assay if a.is_control][0]
                experiments[ip.name] = ExperimentIP(
                    ip_files=ip, control_files=control, **kwargs
                )
            elif len(assay) >= 2 and not any([a.is_control for a in assay]):
                for a in assay:
                    experiments[a.name] = ExperimentIP(ip_files=a, **kwargs)

            elif len(assay) >= 2 and any([a.is_control for a in assay]):
                logger.warning(f"Multiple controls for {assay[0].name}")
                logger.warning("Will generate all possible combinations")
                ip = [a for a in assay if not a.is_control]
                control = [a for a in assay if a.is_control]

                for combination in itertools.product(ip, control):
                    experiments[combination[0].name] = ExperimentIP(
                        ip_files=combination[0], control_files=combination[1], **kwargs
                    )

        return cls(assays=experiments, **kwargs)

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
                # Add the metadata
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

                # Add the experiment
                ip = row["ip"]
                control = row["control"]

                if "control_r1" not in row:
                    experiments[experiment_name] = ExperimentIP(
                        ip_files=AssayIP(
                            name=experiment_name,
                            r1=FastqFileIP(path=row["ip_r1"]),
                            r2=(
                                FastqFileIP(path=row["ip_r2"])
                                if "ip_r2" in row
                                else None
                            ),
                        ),
                        ip=ip,
                        control=None,
                        metadata=metadata,
                    )
                else:
                    experiments[experiment_name] = ExperimentIP(
                        ip_files=AssayIP(
                            name=experiment_name,
                            r1=FastqFileIP(path=row["ip_r1"]),
                            r2=(
                                FastqFileIP(path=row["ip_r2"])
                                if "ip_r2" in row
                                else None
                            ),
                        ),
                        control_files=AssayIP(
                            name=experiment_name,
                            r1=FastqFileIP(path=row["control_r1"]),
                            r2=(
                                FastqFileIP(path=row["control_r2"])
                                if "control_r2" in row
                                else None
                            ),
                        ),
                        ip=ip,
                        control=control,
                        metadata=metadata,
                    )
            else:
                raise NotImplementedError("Not implemented")

        return cls(assays=experiments, **kwargs)


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
        if not self.files:
            files = self.default_files
        else:
            files = self.files

        if self.fastq_screen:
            files.extend(self.fastq_screen_files)
        if self.library_complexity:
            files.extend(self.library_complexity_files)
        return files


class BigWigFiles(BaseModel):
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]
    names: List[str]
    pileup_method: List[Literal["deeptools", "homer"]] = None
    make_bigwigs: bool = False
    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw"]] = None

    @property
    def bigwigs_non_rna(self):
        scale_methods = (
            ["unscaled", self.scale_method] if self.scale_method else ["unscaled"]
        )

        return expand(
            "seqnado_output/bigwigs/{method}/{scale}/{sample}.bigWig",
            sample=self.names,
            scale=scale_methods,
            method=self.pileup_method,
        )

    @property
    def bigwigs_rna(self):
        return expand(
            "seqnado_output/bigwigs/{method}/{scale}/{sample}_{strand}.bigWig",
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
    assay: Literal["ChIP", "ATAC", "RNA", "SNP"]
    names: List[str]
    peak_calling_method: List[Literal["macs", "homer", "lanceotron", "seacr"]] = None
    call_peaks: bool = False

    @property
    def peak_files(self) -> List[str]:
        return expand(
            "seqnado_output/peaks/{method}/{sample}.bed",
            sample=self.names,
            method=self.peak_calling_method,
        )

    @computed_field
    @property
    def files(self) -> List[str]:
        if self.call_peaks:
            return self.peak_files


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
    pileup_method: List[Literal["deeptools", "homer"]] = None
    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw"]] = None

    make_heatmaps: bool = False
    make_ucsc_hub: bool = False

    ucsc_hub_details: Optional[Dict[str, Any]] = None

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
                pileup_method=self.pileup_method,
                scale_method="rpkm",
            )

            return bwf_samples.files + bwf_merged.files
        else:
            return bwf_samples.files

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
        return hbf.files


class RNAOutput(Output):
    assay: Literal["RNA"]
    project_name: str
    run_deseq2: bool = False

    @property
    def counts(self):
        return "seqnado_output/feature_counts/read_counts.tsv"

    @property
    def deseq2(self):
        if self.run_deseq2:
            return [f"deseq2_{self.project_name}.html"]

    @computed_field
    @property
    def files(self) -> List[str]:
        files = self.bigwigs + self.heatmaps + self.ucsc_hub + self.counts + self.design
        if self.run_deseq2:
            files.append(self.deseq2)

        return files


class NonRNAOutput(Output):
    assay: Union[Literal["ChIP"], Literal["ATAC"]]
    call_peaks: bool = False
    peak_calling_method: List[Literal["macs", "homer", "lanceotron"]] = None

    @property
    def merge_peaks(self):
        return "merge" in self.design_dataframe.columns

    @property
    def merged_peaks(self):
        return PeakCallingFiles(
            assay=self.assay,
            names=self.design_dataframe["merge"].unique().tolist(),
            call_peaks=self.call_peaks,
            peak_calling_method=self.peak_calling_method,
        )

    @property
    def peaks(self):
        pcf_samples = PeakCallingFiles(
            assay=self.assay,
            names=self.sample_names,
            call_peaks=self.call_peaks,
            peak_calling_method=self.peak_calling_method,
        )

        if self.merge_peaks:
            pcf_merged = self.merged_peaks

            return pcf_samples.files + pcf_merged.files
        else:
            return pcf_samples.files

    @computed_field
    @property
    def files(self) -> List[str]:
        files = (
            self.bigwigs + self.heatmaps + self.ucsc_hub + self.peaks + self.run_design
        )
        return files


class ATACOutput(NonRNAOutput):
    assay: Literal["ATAC"]


class ChIPOutput(NonRNAOutput):
    assay: Literal["ChIP"]
    ip_names: List[str]
    control_names: List[str]
    call_peaks: bool = False
    peak_calling_method: List[Literal["macs", "homer", "lanceotron", "seacr"]] = None
    chip_spikein_normalisation: bool = False
    scale_method: Optional[Literal["cpm", "rpkm", "spikein", "csaw"]] = None

    @property
    def peaks(self):
        pcf_samples = PeakCallingFiles(
            assay=self.assay,
            names=self.ip_names,
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
        files = (
            self.bigwigs
            + self.heatmaps
            + self.ucsc_hub
            + self.peaks
            + self.spikeins
            + self.design
        )
        return files
