'''
construct the divergent pathways of all strains
divergent pathways <--> networks,
-> nodes: SARS-CoV-2 strains
-> edges: a pair of strains with
          1. same division
          2. date interval <= 14
          3. Hamming distance between their haplotype sequences <= 1
'''

import argparse
import warnings
from multiprocessing import Pool
from datetime import datetime
import pandas as pd


def calculate_day(date_1, date_2):
    '''
    calculate number of days between two dates

    :param date_1: str, date one (YYYY-MM-DD)
    :param date_2: str, date two (YYYY-MM-DD)

    :return: int, number of days between two dates
    '''
    year_1, month_1, day_1 = map(int, date_1.split('-'))
    year_2, month_2, day_2 = map(int, date_2.split('-'))
    return (datetime(year_2, month_2, day_2) - datetime(year_1, month_1, day_1)).days


def is_linked(seq_one, seq_two):
    '''
    determine if the Hamming distance between the haplotype sequences of two strains <= 1

    :param seq_one: snps at key sites of sequence one
    :param seq_two: snps at key sites of sequence two

    :return: boolean, the Hamming distance between the haplotype sequences of two strains <= 1?
    '''
    a = set(seq_one) & set(seq_two)
    b = set(seq_one) | set(seq_two)
    c = b - a
    d = list(set([m.split('_')[0] for m in c]))
    if len(d) <= 1:
        return True


def construct_pathway(seq_divisions_thread, dates, seq_hap):
    '''
    construct divergent pathways
    Note: the divergent pathways are similar with networks, where
          nodes: SARS-CoV-2 strains
          edges: a pair of strains with
                1. same division
                2. date interval <= 14
                3. Hamming distance between their haplotype sequences <= 1
          * Not all links (edges) are calculated, we just cluster the sequences by the theory of networks.
            Simplified and fully computed clustering results are the same.

    :param seq_divisions_thread: list, group(s) of sequences (group <-> division)
    :param dates: dict, key: strain; value: date
    :param seq_hap: dict, key: strain; value: snps at key sites
    
    :return: list, divergent pathways
    '''
    subpathways_thread = []
    
    # construct divergent pathways for each division
    for seq_division in seq_divisions_thread:
        seq_division = sorted(seq_division, key=lambda x: dates[x])
        subpathways = []

        # add the first strain as the first divergent pathway to the empty divergent pathways list
        subpathways.append([seq_division[0]])

        # determine the relationship between each subsequent strain and strains in existing divergent pathways
        # relationship -> 1. date interval 2. Hamming distance between haplotype sequences
        for i in seq_division[1:]:
            pathway_linked = []
            for pathway in subpathways:
                for j in pathway[::-1]:
                    if calculate_day(dates[j], dates[i]) > 14:
                        break
                    if is_linked(seq_hap[i], seq_hap[j]):
                        pathway_linked.append(pathway)
                        break

            # there is no link between this strain and strains in existing divergent pathways
            if len(pathway_linked) == 0:
                subpathways.append([i])
            
            # there is (are) link(s) between this strain and strain(s) in one existing divergent pathway
            elif len(pathway_linked) == 1:
                pathway_linked[0].append(i)

            # there are links between this strain and strains in more than one existing divergent pathways
            else:
                combined_pathway = [i]
                for pathway in pathway_linked:
                    combined_pathway.extend(pathway)
                    subpathways.remove(pathway)
                combined_pathway = sorted(combined_pathway, key=lambda x: dates[x])
                subpathways.append(combined_pathway)

        subpathways_thread.extend(subpathways)

    return subpathways_thread


def main():
    warnings.filterwarnings('ignore')

    # command line interface
    parser = argparse.ArgumentParser(description='Construct divergent pathways')
    parser.add_argument('--threads', type=int, required=True, help='Number of threads')
    parser.add_argument('--haplotypes', required=True, help='Haplotype sequences file')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--output', required=True, help='Divergent pathways file')
    args = parser.parse_args()

    # get haplotype seq
    seq_hap = {}
    with open(args.haplotypes) as f:
        lines = f.readlines()
        for line in lines:
            if line.strip().split(':')[1] == '':
                seq_hap[line.split(':')[0]] = []
            else:
                seq_hap[line.split(':')[0]] = line.strip().split(':')[1].split(',')

    # get meta
    meta_cols = ['strain', 'date', 'region_exposure', 'country_exposure', 'division_exposure']
    meta = pd.read_csv(args.metadata, delimiter='\t', usecols=meta_cols, index_col='strain')
    meta = meta[meta.index.isin(seq_hap)]
    
    # divide strains by division (province / state)
    location_dict = meta[['region_exposure', 'country_exposure', 'division_exposure']].to_dict()
    seq_divisions = {}
    for i in seq_hap:
        division = location_dict['region_exposure'][i]+'_'+location_dict['country_exposure'][i] +'_'+location_dict['division_exposure'][i]
        seq_divisions.setdefault(division, []).append(i)
    
    seq_divisions = sorted(seq_divisions.values(), key=lambda x: len(x), reverse=True)
    
    # multiprocessing
    #
    # each of the top x divisions (with more sequences) is calculated with one seperate thread
    # x = threads number - 1
    # other divisions are calculated with the remaining one thread
    #
    seq_divisions_threads = []
    seq_divisions_remaining = []
    for n, i in enumerate(seq_divisions, start=1):
        if n < args.threads:
            seq_divisions_threads.append([i])
        else:
            seq_divisions_remaining.append(i)
    seq_divisions_threads.append(seq_divisions_remaining)

    # construct divergent pathways
    pathways = []
    p = Pool(args.threads)
    result_list = []
    for i in seq_divisions_threads:
        ids = [k for j in i for k in j]
        dates = meta[meta.index.isin(ids)]['date'].to_dict()
        sub_seq_hap = dict([(key, seq_hap.get(key, None)) for key in ids])
        result = p.apply_async(construct_pathway, args=(i, dates, sub_seq_hap))
        result_list.append(result)
    p.close()
    p.join()
    for result in result_list:
        pathways.extend(result.get())

    pathways = list(sorted(pathways, key=lambda x: len(x), reverse=True))

    with open(args.output, 'w') as f:
        f.write('ID'+ ',' + 'Path' + '\n')
        for n, pathway in enumerate(pathways, 1):
            for i in pathway:
                f.write(i + ',' + str(n) + '\n')


if __name__ == '__main__':
    main()
