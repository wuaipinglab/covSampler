'''
get key sites of SARS-CoV-2 genome
key sites <- genome sites of nonsynonymous mutations that increased in number
             of weekly emerging genomes for four consecutive weeks on at least
             one continent were defined as key sites
'''

import argparse
import warnings
from datetime import datetime
import pandas as pd


def read_mutations(file):
    '''
    read sequences in each group (nonsynonymous mutation + continent) from nonsynonymous.txt

    :param file: str, nonsynonymous file path

    :return: dict, key: group (nonsynonymous mutation + continent); value: accession id
    '''
    seqs_with_mutation = {}
    with open(file) as f:
        for line in f.readlines():
            mutation = line.split(':')[0]
            ids = line.strip().split(':')[1].split(',')
            seqs_with_mutation[mutation] = ids
    return seqs_with_mutation


def get_genome_counts(file):
    with open(file) as f:
        genome_counts = len(f.readlines())
    return genome_counts


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


def is_continuous_increase(subseqs, meta, count_threshold):
    '''
    determine if the number of the nonsynonymous mutation increases for four consecutive weeks in one continent

    :param subseqs: list, accession id of sequences in the group (nonsynonymous mutation + continent)
    :param meta: df, metadata
    :param count_threshold: int, the minimum threshold for the number of nonsynonymous mutation
                                 that can be designated as key nonsynomymous mutation

    :return: boolean, the number of the nonsynonymous mutation increases in one continent for four consecutive weeks?
    '''
    # get number of sequences per date
    dates = {}
    for i in subseqs:
        date = meta.loc[i, 'date']
        dates[date] = dates.get(date, 0) + 1

    # get the earliest date of sequences in the group
    first_date = sorted(dates.items(), key=lambda x: x[0])[0][0]

    # get number of sequences per day
    days = {}
    for i in dates:
        day = calculate_day(first_date, i)+1
        days[day] = dates[i]
    for j in range(1, max(days.keys())):
        days.setdefault(j, 0)
    days = dict(sorted(days.items(), key=lambda x: x[0]))

    # get number of sequences per period (seven days / one week)
    counts_period = []
    count_period = 0
    for i in days:
        count_period += days[i]
        if i % 7 == 0 or i == max(days.keys()):
            counts_period.append(count_period)
            count_period = 0

    # determine if the number of the nonsynonymous mutation increases for four consecutive weeks in one continent
    if len(counts_period) >= 4:
        for i in range(0, len(counts_period)-3):
            if count_threshold < counts_period[i] < counts_period[i+1] < counts_period[i+2] < counts_period[i+3]:
                return True
    return False


def write_new_file(file, key_sites):
    with open(file, 'w') as f:
        f.write('Site'+'\n')
        for i in key_sites:
            f.write(i+'\n')


def main():
    warnings.filterwarnings('ignore')
    
    # command line interface
    parser = argparse.ArgumentParser(description='Get key sites')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--nonsynonymous', required=True, help='Nonsynonymous file')
    parser.add_argument('--snps', required=True, help='SNPs file')
    parser.add_argument('--output', required=True, help='Key sites file')
    args = parser.parse_args()

    seqs_with_mutation = read_mutations(args.nonsynonymous)
    meta = pd.read_csv(args.metadata, delimiter='\t', index_col=2)
    # calculate the minimum threshold for the number of key nonsynonymous
    # mutations in the first week of four consecutive increasing weeks

    # threshold = total number of global sequences / 50,000
    # | total number | 100,000 | 1,000,000 | 10,000,000 |
    # |  threshold   |    2    |     20    |     200    |
    
    count_threshold = get_genome_counts(args.snps) / 50000
    
    key_sites = []
    for mutation in seqs_with_mutation:
        if mutation.split('_')[0] not in key_sites:
            if is_continuous_increase(seqs_with_mutation[mutation], meta, count_threshold):
                key_sites.append(mutation.split('_')[0])
    
    key_sites = sorted(key_sites, key=lambda x:int(x))

    write_new_file(args.output, key_sites)

if __name__ == '__main__':
    main()
