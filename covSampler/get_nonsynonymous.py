'''
1. get nonsynonymous mutations in SARS-CoV-2 whole genome
2. get accession id of sequences in each group (group -> nonsynonymous mutation + continent)
'''

import warnings
from Bio.Seq import Seq
import pandas as pd
from covSampler import genome_path, meta_path, snps_path, nonsynonymous_path


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
        num = 0
        for m in codon:
            num += 1
            base_list = ['A', 'T', 'C', 'G']
            base_list.remove(m)
            for base in base_list:
                codon_new = get_new_codon(codon, base, num)
                if Seq(codon_new).translate() != Seq(codon).translate():
                    nonsynonymous.append(str(codon_site[num-1]) + '_' + base)
    return nonsynonymous


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


def write_new_file(file, groups):
    with open(file, 'w') as f:
        for group in groups:
            f.write(group+':'+','.join(groups[group])+'\n')


def main():
    warnings.filterwarnings('ignore')
    genome = pd.read_csv(genome_path, index_col=0)
    nonsynonymous = get_nonsynonymous_mutation(genome)
    
    seqs = read_seq(snps_path)
    seq_continent = pd.read_csv(meta_path, delimiter='\t', index_col=2)['region_exposure'].to_dict()
    
    groups = {}
    for i in seqs:
        # get groups (nonsynonymous mutation + continent) that the sequence is in
        continent = seq_continent[i]
        seq_nonsynonymous = list(set(seqs[i]) & set(nonsynonymous))
        seq_nonsynonymous.sort()
        for n in seq_nonsynonymous:
            groups.setdefault(n+'_'+continent, []).append(i)

    write_new_file(nonsynonymous_path, groups)


if __name__ == '__main__':
    main()
