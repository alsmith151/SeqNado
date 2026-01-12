[Back to Index](index.md)

# Pipeline

After successful Design of a SeqNado run [Design](design.md)

The SeqNado pipeline processes multiomics data through the following steps:

1. **Preprocessing**: Quality control and trimming
2. **Alignment**: Mapping reads to the genome
3. **Analysis**: Generating peaks, quantifications, and visualizations

## Input/Output
- Input: FASTQ files
- Output: Processed data in `seqnado_output/`

## Helper Functions
SeqNado includes helper functions for:
- Extracting cores from options (`extract_cores_from_options` in `helpers.py`)
- File handling and logging

These utilities streamline the pipeline execution.



## Next Steps

SeqNado output will be stored in `seqnado_output/`