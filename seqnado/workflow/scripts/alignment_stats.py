import pandas as pd
from loguru import logger

logger.add(snakemake.log[0], level="INFO")

wide_data = {}
for f in snakemake.input:
    sample = f.split('/')[-1].replace('alignment_stats_', '').replace('.tsv', '')
    df = pd.read_csv(f, sep='\t')
    wide_data[sample] = df.set_index('Step')['Read Count']

combined_df = pd.DataFrame(wide_data).T.reset_index()
combined_df.rename(columns={'index': 'sample'}, inplace=True)
combined_df.to_csv(snakemake.output[0], sep='\t', index=False)
logger.info("Alignment stats table written to {}".format(snakemake.output[0]))
