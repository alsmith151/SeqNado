import os
import sys
import argparse
import pandas as pd

def main(tsv, header=False, col_names=None, output=sys.stdout, query=None, output_header=False):
   
    df = pd.read_csv(sys.stdin if tsv == '-' else tsv,
                     sep='\t',
                     header=0 if header else None,
                     names=col_names if col_names else None)
    if query:
        df = df.query(query)

    df.to_csv(output, sep='\t', header=output_header, index=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tsv', help='tsv to process')
    parser.add_argument('--header', help='File header', default=None, action='store_true')
    parser.add_argument('--col_names', help='Column names', default=None, nargs='+')
    parser.add_argument('-o', '--output', help='Filename to write. Default stdout', default=sys.stdout)
    parser.add_argument('-q', '--query', help='Query to apply to dataframe')
    parser.add_argument('--output_header', help='Output the header or not', action='store_true', default=False)
    args = parser.parse_args()

    main(**vars(args))
