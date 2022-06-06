'''
construct haplotype sequence for each sequence
haplotype sequence: a short pseudo sequence composed of key sites of genome
'''

import argparse
import pandas as pd


def read_seq(file):
    '''
    read sequence snps from snps.txt

    :param file: str, snps file path
    
    :return: dict, key: accession id; value: snps
    '''
    seqs = {}
    with open(file) as f:
        for line in f.readlines():
            seq_id = line.split(':')[0]
            snps = line.strip().split(':')[1].split(',')
            seqs[seq_id] = snps

    return seqs


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Construct haplotype sequences')
    parser.add_argument('--sites', required=True, help='Key sites file')
    parser.add_argument('--snps', required=True, help='SNPs file')
    parser.add_argument('--output', required=True, help='Haplotype sequences file')
    args = parser.parse_args()

    seqs = read_seq(args.snps)
    key_sites = pd.read_csv(args.sites)['Site'].astype(str).values.tolist()
    haplotype_sequence = {}
    for i in seqs:
        haplotype_sequence[i] = []
        for snp in seqs[i]:
            if snp.split('_')[0] in key_sites and snp.split('_')[1] != '-':
                haplotype_sequence[i].append(snp)

    with open(args.output, 'w') as f:
        for i in haplotype_sequence:
            f.write(i + ':' + ','.join(haplotype_sequence[i]) + '\n')


if __name__ == '__main__':
    main()
