'''
get key sites in SARS-CoV-2 genome
key sites <- genome sites of nonsynonymous muts that increased in proportion
             for four consecutive weeks on at least one continent
'''

import argparse
import warnings
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


def is_key_site(seqs, seq_date, first_date, day_count, count_threshold):
    '''
    determine if the proportion of the nonsynonymous mut increases for four consecutive weeks in the continent

    :param seqs: list, strains in the group (group -> nonsynonymous mut + continent)
    :param seq_date: dict, key: strain; value: date
    :param first_date: str, date of the earliest strain in whole date set
    :param day_count: dict, key: day; value: number of strains in the continent
    :param count_threshold: int, the minimum threshold for the absolute number of nonsynonymous muts per week

    :return: boolean, the proportion of the nonsynonymous mut increases for four consecutive weeks in the continent?
    '''
    # get number of the nonsynonymous mut per day
    days = {}
    for i in seqs:
        day = calculate_day(first_date, seq_date[i])
        days[day] = days.get(day, 0) + 1

    for j in range(min(days.keys()), max(days.keys())):
        days.setdefault(j, 0)

    days = dict(sorted(days.items(), key=lambda x: x[0]))

    # get proportion of the nonsynonymous mut per period (seven days / one week)
    proportions = []
    count_7d = 0
    for n, i in enumerate(days, 1):
        count_7d += days[i]
        if n % 7 == 0 or i == max(days.keys()):
            if count_7d != 0:
                proportions.append((count_7d, count_7d / sum([day_count[d] for d in day_count if (i-6<=d<=i)])))
            else:
                proportions.append((0, 0))
            count_7d = 0

    # determine if the proportion of the nonsynonymous mut grows for three consecutive weeks in one continent
    if len(proportions) >= 4:
        for i in range(0, len(proportions)-3):
            if count_threshold < min(proportions[i][0], proportions[i+1][0], proportions[i+2][0], proportions[i+3][0]):
                if proportions[i][1] < proportions[i+1][1] < proportions[i+2][1] < proportions[i+3][1]:
                    return True
    return False


def main():
    warnings.filterwarnings('ignore')
    
    # command line interface
    parser = argparse.ArgumentParser(description='Get key sites')
    parser.add_argument('--strains', required=True, help='Strains file')
    parser.add_argument('--nonsynonymous', required=True, help='Nonsynonymous file')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--output', required=True, help='Key sites file')
    args = parser.parse_args()

    # get strains
    strains = []
    with open(args.strains) as f:
        for n, line in enumerate(f.readlines()):
            if n != 0:
                strains.append(line.strip())

    # get strains within each group (group -> nonsynonymous mut + continent)
    seq_groups = {}
    with open(args.nonsynonymous) as f:
        for line in f.readlines():
            group = line.split(':')[0]
            ids = line.strip().split(':')[1].split(',')
            seq_groups[group] = ids

    # get meta
    meta = pd.read_csv(args.metadata, delimiter='\t', usecols=['strain', 'date', 'region_exposure'], index_col='strain')
    meta = meta[meta.index.isin(strains)]
    
    first_date = meta['date'].drop_duplicates().sort_values().values[0]

    meta_continents = {}
    for c in meta['region_exposure'].drop_duplicates():
        submeta = meta[meta['region_exposure'] == c]
        meta_continents[c] = {}
        meta_continents[c]['seq_date'] = submeta['date'].to_dict()
        meta_continents[c]['day_count'] = dict(
            [(calculate_day(first_date, date), count) for date, count in submeta['date'].value_counts().to_dict().items()]
            )

    # get the minimum threshold for the absolute number of nonsynonymous muts per week
    #
    # threshold = the total number of sequences in the original data set / 50,000
    # | total number | 100,000 | 1,000,000 | 10,000,000 |
    # |  threshold   |    2    |     20    |     200    |
    #
    count_threshold = len(strains) / 50000
    
    # get key sites
    key_sites = []
    for group in seq_groups:
        if group.split('_')[0][1:-1] not in key_sites:
            subseqs = seq_groups[group]
            meta_continent = meta_continents[group.split('_')[-1]]
            if is_key_site(
                subseqs, meta_continent['seq_date'], first_date, meta_continent['day_count'], count_threshold
                ):
                key_sites.append(group.split('_')[0][1:-1])

    key_sites = sorted(key_sites, key=lambda x: int(x))

    with open(args.output, 'w') as f:
        f.write('Site'+'\n')
        for i in key_sites:
            f.write(i+'\n')


if __name__ == '__main__':
    main()
