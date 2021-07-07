import os
import sys
import pandas as pd
import numpy as np
import argparse

p = argparse.ArgumentParser()
p.add_argument('-d', '--destination_directory', help='Directory to place renamed symlinks', default='.')
args = p.parse_args()




def main():
    fnames = os.listdir()
    df = (pd.DataFrame(fnames, columns=['fn'])
            .loc[lambda df: df['fn'].str.contains(r'.*.fastq.gz')]
            .assign(base=lambda df: df['fn'].str.split('_').str[1:].str.join('_'),
                    dest=lambda df: args.destination_directory.rstrip('/') + '/' + df['base'],
                    fp=lambda df: [os.path.abspath(fn) for fn in df['fn']])
            .replace('', np.nan)
            .dropna(how='any')
         )
    
    for i, (src, dest) in df[['fp', 'dest']].iterrows():
        try:
            os.symlink(src, dest)
        except Exception as e:
            print(e)   


if __name__ == '__main__':
    main()
