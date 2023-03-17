import os
import pathlib
import pandas as pd
import re
from typing import List
from collections import defaultdict

pd.set_option("mode.chained_assignment", None)


def sample_names_follow_convention(
    df: pd.DataFrame, name_column: str = "basename"
) -> bool:
    naming_pattern_paired = r"(.*)_(.*)_R?[12].fastq(.gz)?"
    naming_pattern_single = r"(.*)_(.*).fastq(.gz)?"

    return (
        df[name_column].str.match(naming_pattern_paired)
        | df[name_column].str.match(naming_pattern_single)
    ).all()


class ChipseqFastqSamples:
    def __init__(self, design):

        # Expected columns: sample, antibody, fq1, fq2, control
        self.design = design
        self.design = self.design.assign(
            paired=(~self.design[["fq1", "fq2"]].isna().any(axis=1))
        )

    @classmethod
    def from_files(cls, files: List) -> "ChipseqFastqSamples":

        df = pd.DataFrame(files, columns=["fn"])

        df[["sample", "read"]] = (
            df["fn"].apply(str).str.extract("(?!.*/)?(.*)_.*_R?([12]).fastq.gz")
        )

        df["sample"] = df["sample"].apply(lambda p: pathlib.Path(p).name)
        df["read"] = "fq" + df["read"]

        df["antibody"] = df["fn"].astype(str).str.split("_").str[-2]

        df = (
            df.pivot(columns="read", index=["sample", "antibody"])
            .droplevel(level=0, axis=1)
            .reset_index()
        )

        df_input = df.loc[df["antibody"].str.lower().str.contains("input")]
        df_input = df_input.assign(
            control=df_input["sample"] + "_" + df_input["antibody"]
        )
        df_ip = df.loc[~df["antibody"].str.lower().str.contains("input")]
        df = df_ip.merge(df_input[["sample", "control"]], on="sample", how="left")

        return cls(design=df)

    @property
    def fastq_ip_files(self):

        fastq_files = []

        for sample in self.design.itertuples():
            if sample.paired:
                for ii, fq in enumerate([sample.fq1, sample.fq2]):
                    fastq_files.append(fq)
            else:
                for fq in getattr(sample, "fq1"):
                    fastq_files.append(fq)
        return sorted(fastq_files)

    @property
    def fastq_control_files(self):

        fastq_files = []

        for sample in self.design.itertuples():
            path = pathlib.Path(sample.fq1).parent
            fqs = [fq for fq in path.glob(f"{sample.control}*.fastq.gz")]
            for fq in fqs:
                fastq_files.append(str(fq))

        return sorted(list(set(fastq_files)))

    @property
    def fastq_files(self):
        return sorted([*self.fastq_ip_files, *self.fastq_control_files])

    @property
    def sample_names_all(self):
        samples_ip = (
            self.design["sample"].apply(pathlib.Path).apply(lambda p: p.name)
            + "_"
            + self.design["antibody"]
        )
        samples_control = pd.Series(
            self.design["control"]
            .dropna()
            .apply(lambda p: pathlib.Path(p).name)
            .unique()
        )
        return pd.concat([samples_ip, samples_control]).to_list()

    @property
    def sample_names_ip(self):
        samples_ip = (
            self.design["sample"].apply(pathlib.Path).apply(lambda p: p.name)
            + "_"
            + self.design["antibody"]
        )
        return samples_ip

    @property
    def sample_names_control(self):
        samples_control = pd.Series(
            self.design["control"]
            .dropna()
            .apply(lambda p: pathlib.Path(p).name)
            .unique()
        )
        return samples_control

    @property
    def antibodies(self):
        return self.design["antibody"].unique()

    @property
    def paired_ip_and_control(self):
        _design = self.design.assign(
            treatment=lambda df: df["sample"] + "_" + df["antibody"]
        )
        return _design.set_index("treatment")["control"].to_dict()


    def _translate_control_samples(self):
        fq_translation = {}
        for fq in self.fastq_control_files:
            for control in self.design["control"]:
                if str(control) in fq:
                    
                    read = re.match(r".*/?.*_R?([12])(?:_001)?.fastq.gz", fq).group(1)
                    fq_translation[f"{control}_{read}.fastq.gz"] = os.path.abspath(fq)

        return fq_translation
    
    def _translate_ip_samples(self):

        fq_translation = {}
        for sample in self.design.itertuples():
            for read, fq in enumerate([sample.fq1, sample.fq2]):

                if os.path.exists(fq):
                    fq_translation[f"{sample.sample}_{sample.antibody}_{read + 1}.fastq.gz"] = os.path.realpath(fq)

        return fq_translation
    
    @property
    def translation(self):
        """Create a dictionary with the fastq files and their new names"""
        fq_translation = {}
        fq_translation.update(self._translate_ip_samples())
        fq_translation.update(self._translate_control_samples())
        return fq_translation

    
