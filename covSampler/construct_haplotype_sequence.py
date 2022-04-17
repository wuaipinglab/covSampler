'''
construct haplotype sequence for each sequence
haplotype sequence: a short pseudo sequence composed of key sites of genome
'''

import pandas as pd
from covSampler import snps_path, key_site_path, haplotype_sequence_path
from get_nonsynonymous import read_seq


def main():
    seqs = read_seq(snps_path)
    key_sites = pd.read_csv(key_site_path)['Site'].astype(str).values.tolist()
    haplotype_sequence = {}
    for i in seqs:
        haplotype_sequence[i] = []
        for snp in seqs[i]:
            if snp.split('_')[0] in key_sites and snp.split('_')[1] != '-':
                haplotype_sequence[i].append(snp)

    with open(haplotype_sequence_path, 'w') as f:
        for i in haplotype_sequence:
            f.write(i + ':' + ','.join(haplotype_sequence[i]) + '\n')


if __name__ == '__main__':
    main()
