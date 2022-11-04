import os
from typing import Dict, List
import pathlib
import pybedtools
import pandas as pd
import re
import glob
import numpy as np
from typing import List, Literal
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


class ChipseqFastqSamples():
    def __init__(self, design):
        
        # Expected columns: sample, antibody, fq1, fq2, control
        self.design = design
        self.design = self.design.assign(paired=(~self.design[["fq1", "fq2"]].isna().any(axis=1)))
 
    
    @classmethod
    def from_files(cls, files: List) -> (pd.DataFrame):

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
        df_input = df_input.assign(control=df_input["sample"] + "_" + df_input["antibody"])
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
        samples_ip = (self.design["sample"].apply(pathlib.Path).apply(lambda p: p.name) + "_" + self.design["antibody"])
        samples_control = pd.Series(self.design["control"].dropna().apply(lambda p: pathlib.Path(p).name).unique())
        return pd.concat([samples_ip, samples_control]).to_list()
    
    @property
    def sample_names_ip(self):
        samples_ip = (self.design["sample"].apply(pathlib.Path).apply(lambda p: p.name) + "_" + self.design["antibody"])
        return samples_ip
    
    @property
    def sample_names_control(self):
        samples_control = pd.Series(self.design["control"].dropna().apply(lambda p: pathlib.Path(p).name).unique())
        return samples_control
    
    @property 
    def antibodies(self):
        return self.design["antibody"].unique()
    
    @property
    def paired_ip_and_control(self):
        _design = self.design.assign(treatment=lambda df: df["sample"] + "_" + df["antibody"])
        return _design.set_index("treatment")["control"].to_dict()
        
    
    def symlink_fastq_files(self, outdir="fastq"):
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        # Control samples
        for fq in self.fastq_control_files:
            for control in self.design["control"]:
                if str(control) in fq:
                    src = os.path.abspath(fq)
                    read = re.match(r".*/?.*_R?([12])(?:_001)?.fastq.gz", fq).group(1)
                    dest = os.path.join(outdir, f"{control}_{read}.fastq.gz")
                    
                    try:
                        os.symlink(src, dest)
                    except FileExistsError:
                        pass

        # IP samples
        fq_links = defaultdict(list) 
        for sample in self.design.itertuples():
            for read, fq in enumerate([sample.fq1, sample.fq2]):
                
                if os.path.exists(fq):
                    src = os.path.abspath(fq)
                    dest = os.path.join(outdir, f"{sample.sample}_{sample.antibody}_{read + 1}.fastq.gz")
                    fq_links[f"fq{read + 1}"].append(dest)
                    
                    try:
                        os.symlink(src, dest)
                    except FileExistsError:
                        pass
            
        self.design["fq1"] = fq_links["fq1"]
        self.design["fq2"] = fq_links["fq2"]
    

