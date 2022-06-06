'''
construct the divergent pathways of global SARS-CoV-2 sequences
divergent pathways <--> networks,
-> nodes: SARS-CoV-2 sequences
-> edges: a pair of sequences with
          1. same division
          2. date interval <= 14
          3. hamming distance between their haplotype sequences <= 1
          (haplotype sequence is a short pseudo sequence composed of key sites of genome)
'''

import argparse
import warnings
from multiprocessing import Pool
from datetime import datetime
import pandas as pd


def read_haplotype_sequence(file):
    '''
    read haplotype sequences (constructed by combining pre-calculated
    SARS-CoV-2 key sites) of sequences from haplotype_sequence.txt

    :param file: str, haplotype sequence file path

    :return: dict, key: accession id; value: snps at key sites
    '''
    seq_vs = {}
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if line.strip().split(':')[1] == '':
                seq_vs[line.split(':')[0]] = []
            else:
                seq_vs[line.split(':')[0]] = line.strip().split(':')[1].split(',')
    return seq_vs


def group_by_division(seq_vs, meta):
    '''
    group sequences by division (province / state)

    :param seq_vs: dict, key: accession id; value: snps at key sites
    :param meta: df, metadata

    :return: list, accession id of sequences in each group
    '''
    seq_groups = {}
    for i in seq_vs:
        division = meta.loc[i, 'region_exposure'] \
                   + '_' + meta.loc[i, 'country_exposure'] \
                   + '_' + meta.loc[i, 'division_exposure']
        seq_groups.setdefault(division, []).append(i)
    seq_groups = list(sorted(seq_groups.values(), key=lambda x: len(x), reverse=True))
    return seq_groups


def calculate_day(date_1, date_2):
    '''
    calculate number of days between two dates

    :param date_1: str, date one (YYYY-MM-DD)
    :param date_2: str, date two (YYYY-MM-DD)

    :return: int, number of days between two dates
    '''
    year_1 = int(date_1.split('-')[0])
    month_1 = int(date_1.split('-')[1])
    day_1 = int(date_1.split('-')[2])
    year_2 = int(date_2.split('-')[0])
    month_2 = int(date_2.split('-')[1])
    day_2 = int(date_2.split('-')[2])
    return (datetime(year_2, month_2, day_2) - datetime(year_1, month_1, day_1)).days


def is_linked(seq_one, seq_two):
    '''
    determine if the hamming distance between the haplotype sequences of two sequences <= 1

    :param seq_one: snps at key sites of sequence one
    :param seq_two: snps at key sites of sequence two

    :return: boolean, the hamming distance between the haplotype sequences of two sequences <= 1?
    '''
    a = set(seq_one) & set(seq_two)
    b = set(seq_one) | set(seq_two)
    c = b - a
    d = list(set([m.split('_')[0] for m in c]))
    if len(d) <= 1:
        return True


def construct_pathway(seq_groups_thread, submeta, seq_vs):
    '''
    construct divergent pathways
    Note: the divergent pathways are similar with networks, where
          nodes: SARS-CoV-2 sequences
          edges: a pair of sequences with
                1. same division
                2. date interval <= 14
                3. hamming distance between their haplotype sequences <= 1
          * not all links (edges) are calculated, we just cluster the sequences by the theory of networks
            the clustering results of simplified calculation and full calculation are the same

    :param seq_groups_thread: list, group(s) of sequences (group <-> division)
    :param submeta: dict, key: accession id; value: date
    :param seq_vs: dict, key: accession; value: snps at key sites
    
    :return: list, divergent pathways
    '''
    subpathways_thread = []
    # construct divergent pathways of each group (division)
    for seq_group in seq_groups_thread:
        # sort sequence by date
        seq_group = sorted(seq_group, key=lambda x: submeta[x])
        subpathways = []

        # add the first sequence as the first divergent pathway to the empty divergent pathways list
        subpathways.append([seq_group[0]])

        # determine the relationship between each subsequent sequence and sequences in existing divergent pathways
        # relationship -> 1. date interval 2. hamming distance between haplotype sequences
        for i in seq_group[1:]:
            pathway_linked = []
            for pathway in subpathways:
                for j in pathway[::-1]:
                    if calculate_day(submeta[j], submeta[i]) > 14:
                        break
                    if is_linked(seq_vs[i], seq_vs[j]):
                        pathway_linked.append(pathway)
                        break

            # there is no link between the sequence and sequences in existing divergent pathways
            if len(pathway_linked) == 0:
                subpathways.append([i])
            
            # there is (are) link(s) between the sequence and sequence(s) in one existing divergent pathway
            elif len(pathway_linked) == 1:
                pathway_linked[0].append(i)

            # there are links between the sequence and sequences in more than one existing divergent pathway
            else:
                pathway_merged = [i]
                for pathway in pathway_linked:
                    pathway_merged.extend(pathway)
                    subpathways.remove(pathway)
                pathway_merged = sorted(pathway_merged, key=lambda x: submeta[x])
                subpathways.append(pathway_merged)

        subpathways_thread.extend(subpathways)

    return subpathways_thread


def write_new_file(pathway_file, pathways):
    with open(pathway_file, 'w') as f:
        f.write('ID'+ ',' + 'Path' + '\n')
        pw = 0
        for pathway in pathways:
            pw += 1
            for i in pathway:
                f.write(i + ',' + str(pw) + '\n')


def main():
    warnings.filterwarnings('ignore')

    # command line interface
    parser = argparse.ArgumentParser(description='Construct divergent pathways')
    parser.add_argument('--threads', type=int, required=True, help='Number of threads')
    parser.add_argument('--haplotypes', required=True, help='Haplotype sequences file')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--output', required=True, help='Divergent pathways file')
    args = parser.parse_args()

    meta = pd.read_csv(args.metadata, delimiter='\t', index_col=2)
    seq_vs = read_haplotype_sequence(args.haplotypes)
    seq_groups = group_by_division(seq_vs, meta)
    
    # calculate with multiprocessing
    # each of the top x groups (with more sequences) is calculated with one seperate thread
    # x = threads - 1
    # other groups are calculated with the remaining one thread
    seq_groups_threads = []
    seq_groups_merged = []
    for n, i in enumerate(seq_groups, start=1):
        if n < args.threads:
            seq_groups_threads.append([i])
        else:
            seq_groups_merged.append(i)
    seq_groups_threads.append(seq_groups_merged)

    pathways = []
    p = Pool(args.threads)
    result_list = []
    for i in seq_groups_threads:
        # get accession id, date and snps at key sites of sequences calculated with the thread
        ids = [k for j in i for k in j]
        submeta = meta[meta.index.isin(ids)]['date'].to_dict()
        sub_seq_vs = dict([(key, seq_vs.get(key, None)) for key in ids])
        result = p.apply_async(construct_pathway, args=(i, submeta, sub_seq_vs))
        result_list.append(result)
    p.close()
    p.join()
    for result in result_list:
        pathways.extend(result.get())

    pathways = list(sorted(pathways, key=lambda x: len(x), reverse=True))

    write_new_file(args.output, pathways)


if __name__ == '__main__':
    main()
