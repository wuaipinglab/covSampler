'''
compress metadata for webserver
'''

import os
import pandas as pd
import subprocess

DATE = '2022-05-05'

DIRPATH = os.path.join('/datadirectory/', DATE)

rawdata_path = os.path.join(DIRPATH, 'rawdata/')
for root, dirs, files in os.walk(rawdata_path):
    for file in files:
        if file.startswith('metadata') and file.endswith('.tsv'):
            meta_path = os.path.join(root, file)
infos_path = os.path.join(DIRPATH, 'infos.tsv')
outtsv_path = os.path.join(DIRPATH, 'split_metadata.tsv')


def main():
    seqs = []
    with open(infos_path) as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith('ID'):
                seq_id = line.split('\t')[0]
                seqs.append(seq_id)

    meta = pd.read_csv(meta_path, delimiter='\t', index_col=2)
    meta_split = meta[meta.index.isin(seqs)][['date', 'region_exposure', 'pango_lineage']]
    meta_split.to_csv(outtsv_path, sep='\t')
    subprocess.run('gzip '+outtsv_path, shell=True)


if __name__ == '__main__':
    main()
