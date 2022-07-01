'''
construct haplotype sequence for each strain
haplotype sequence: a short pseudo sequence composed of key genome sites
'''

import argparse
import pandas as pd


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Construct haplotype sequences')
    parser.add_argument('--strains', required=True, help='Strains file')
    parser.add_argument('--nextclade-tsv', required=True, help='Nextclade tsv file')
    parser.add_argument('--sites', required=True, help='Key sites file')
    parser.add_argument('--output', required=True, help='Haplotype sequences file')
    args = parser.parse_args()

    # get strains
    strains = []
    with open(args.strains) as f:
        for n, line in enumerate(f.readlines()):
            if n != 0:
                strains.append(line.strip())

    # get substitutions
    nextclade = pd.read_csv(args.nextclade_tsv, delimiter='\t', usecols=['seqName', 'substitutions'], index_col='seqName')
    nextclade = nextclade[nextclade.index.isin(strains)]
    nextclade.fillna('', inplace=True)

    seq_substitutions = dict(
        [(strain, substitutions.split(',')) for strain, substitutions in nextclade['substitutions'].to_dict().items()]
        )
    
    # get key sites
    key_sites = []
    with open(args.sites) as f:
        for n, line in enumerate(f.readlines()):
            if n != 0:
                key_sites.append(line.strip())

    # get haplotype sequences
    haplotype_sequences = {}
    for i in seq_substitutions:
        haplotype_sequences[i] = []
        for substitution in seq_substitutions[i]:
            if substitution[1:-1] in key_sites:
                haplotype_sequences[i].append(substitution[1:-1]+'_'+substitution[-1])

    with open(args.output, 'w') as f:
        for i in haplotype_sequences:
            f.write(i + ':' + ','.join(haplotype_sequences[i]) + '\n')


if __name__ == '__main__':
    main()
