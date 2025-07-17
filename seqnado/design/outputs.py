



    def to_geo_dataframe(
        self,
        assay: Literal["ATAC", "RNA", "SNP"],
        pipeline_config: dict[str, Any],
    ) -> pd.DataFrame:
        """
        Generate GEO-compliant metadata table for submission.

        Args:
            assay: One of "ATAC", "RNA", "SNP".
            pipeline_config: Must contain genome.name and optionally instrument_model.

        Returns:
            pandas DataFrame ready for GEO.
        """
        # Build GEOSample entries
        geo_samples: list[GEOSample] = []
        df = self.to_dataframe()
        for row in df.itertuples(index=False):
            # Determine processed files
            if assay == "RNA":
                proc_files = ["read_counts.tsv", f"{row.sample_name}_plus.bw", f"{row.sample_name}_minus.bw"]
            else:
                proc_files = [f"{row.sample_name}.bw"]

            raw_files = [pathlib.Path(row.r1).name]
            if row.r2:
                raw_files.append(pathlib.Path(row.r2).name)

            sample = GEOSample(
                assay=assay,
                library_name=row.sample_name,
                title=row.sample_name,
                organism=predict_organism(pipeline_config["genome"]["name"]),
                cell_line=None,
                cell_type=None,
                antibody=None,
                genotype=None,
                treatment=getattr(row, "treatment", None),
                time=getattr(row, "time", None),
                single_or_paired=("paired-end" if row.r2 else "single"),
                instrument_model=pipeline_config.get("instrument_model", "Illumina NovaSeq X"),
                description=None,
                processed_data_file=proc_files,
                raw_file=raw_files,
            )
            geo_samples.append(sample)

        return GEOSamples(samples=geo_samples).to_dataframe()




    def to_geo_dataframe(
        self,
        assay: Literal["ChIP", "CAT"],
        pipeline_config: dict[str, Any],
    ) -> pd.DataFrame:
        """
        Generate GEO metadata for IP/control experiments.

        Each sample (IP and control) becomes a row.
        """
        geo_samples: list[GEOSample] = []
        df = self.to_dataframe()
        for rec in df.itertuples(index=False):
            for side in ("ip", "control"):  # produce both
                label = getattr(rec, side)
                if label is None:
                    continue

                r1 = getattr(rec, f"{side}_r1")
                r2 = getattr(rec, f"{side}_r2")
                raw_files = [pathlib.Path(r1).name]
                if r2:
                    raw_files.append(pathlib.Path(r2).name)

                proc_file = f"{rec.sample_name}_{label}.bigWig"
                sample = GEOSample(
                    assay=assay,
                    library_name=rec.sample_name,
                    title=rec.sample_name,
                    organism=predict_organism(pipeline_config["genome"]["name"]),
                    cell_line=None,
                    cell_type=None,
                    antibody=label,
                    genotype=None,
                    treatment=getattr(rec, "treatment", None),
                    time=getattr(rec, "time", None),
                    single_or_paired=("paired-end" if r2 else "single"),
                    instrument_model=pipeline_config.get("instrument_model", "Illumina NovaSeq X"),
                    description=None,
                    processed_data_file=[proc_file],
                    raw_file=raw_files,
                )
                geo_samples.append(sample)

        return GEOSamples(samples=geo_samples).to_dataframe()
