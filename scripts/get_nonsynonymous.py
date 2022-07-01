'''
1. get nonsynonymous mutations in SARS-CoV-2 whole genome
2. get strains in each group (group -> nonsynonymous mutation + continent)
'''

import argparse
import warnings
from Bio.Seq import Seq
import pandas as pd


def get_new_codon(codon, change_base, change_base_pos):
    '''
    get codon after single site mutation

    :param codon: str, codon before mutation
    :param change_base: str, mutated base
    :param change_base_pos: int, position of mutated base in codon [1-3]

    :return: str, codon after single site mutation
    '''
    if change_base_pos == 1:
        return change_base + codon[1:]
    elif change_base_pos == 2:
        return codon[0] + change_base + codon[2]
    elif change_base_pos == 3:
        return codon[:2] + change_base


def get_nonsynonymous_mutation(genome):
    '''
    get nonsynonymous mutations in SARS-CoV-2 whole genome

    :param genome: df, SARS-CoV-2 genome information (gene, product, aaPos, peptidePos, aa, genomePos, nucleotide)

    :return: list, nonsynonymous mutations in SARS-CoV-2 whole genome
    '''
    nonsynonymous = []
    peptidepos = genome['peptidePos'].drop_duplicates()
    for i in peptidepos:
        codon = ''.join(list(genome[genome['peptidePos'] == i]['nucleotide']))
        codon_site = list(genome[genome['peptidePos'] == i]['genomePos'])
        for n, m in enumerate(codon, 1):
            base_list = ['A', 'T', 'C', 'G']
            base_list.remove(m)
            for base in base_list:
                codon_new = get_new_codon(codon, base, n)
                if Seq(codon_new).translate() != Seq(codon).translate():
                    nonsynonymous.append(m+str(codon_site[n-1])+base)
    return nonsynonymous


def main():
    warnings.filterwarnings('ignore')

    # command line interface
    parser = argparse.ArgumentParser(description='Get nonsynonymous')
    parser.add_argument('--genome', required=True, help='SARS-CoV-2 genome file')
    parser.add_argument('--strains', required=True, help='Strains file')
    parser.add_argument('--nextclade-tsv', required=True, help='Nextclade tsv file')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--output', required=True, help='Nonsynonymous file')
    args = parser.parse_args()

    # get nonsynonymous mut in sars-cov-2 genome
    genome = pd.read_csv(args.genome, index_col=0)
    nonsynonymous = get_nonsynonymous_mutation(genome)

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
    
    # get continent
    meta = pd.read_csv(args.metadata, delimiter='\t', usecols=['strain', 'region_exposure'], index_col='strain')
    meta = meta[meta.index.isin(strains)]
    seq_continent = meta['region_exposure'].to_dict()
    
    # get groups (key: nonsynonymous mut + continent, value: strains)
    groups = {}
    for i in seq_substitutions:
        continent = seq_continent[i]
        seq_nonsynonymous = list(set(seq_substitutions[i]) & set(nonsynonymous))
        for n in seq_nonsynonymous:
            groups.setdefault(n+'_'+continent, []).append(i)

    with open(args.output, 'w') as f:
        for group in groups:
            f.write(group+':'+','.join(groups[group])+'\n')


if __name__ == '__main__':
    main()
