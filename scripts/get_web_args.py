'''
prepare the content of each parameter in the web interface
'''

import os
import argparse
from datetime import datetime
import pandas as pd


def read_infos(file):
    '''
    read information of each sequence from infos.tsv

    :param file: str, info file path

    :return: dict, key: accession id; value: infos of sequences
    '''
    infos = {}
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith('ID'):
                lineList = line.strip().split('\t')
                seqID = lineList[0]
                seqRegion = lineList[1]
                seqCountry = lineList[2]
                seqDivision = lineList[3]
                seqDate = lineList[4]
                seqPL = lineList[5]
                seqNC = lineList[6]
                seqGC = lineList[7]
                seqNT = lineList[8]
                seqAA = lineList[9]
                infos[seqID] = {
                    'region': seqRegion,
                    'country': seqCountry,
                    'division': seqDivision,
                    'date': seqDate,
                    'pangoLineage': seqPL,
                    'nextstrainClade': seqNC,
                    'gisaidClade': seqGC,
                    'nt': seqNT,
                    'aa': seqAA
                }
    
    return infos


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Get web args')
    parser.add_argument('--infos', required=True, help='Infos file')
    parser.add_argument('--genome', required=True, help='SARS-CoV-2 genome file')
    parser.add_argument('--output', required=True, help='Web args directory')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    infos = read_infos(args.infos)
    locations = []
    dates = []
    pango_lineages = []
    nextstrain_clades = []
    gisaid_clades = []
    for i in infos:
        for location in [infos[i]['region'], infos[i]['region']+'/'+infos[i]['country'], infos[i]['region']+'/'+infos[i]['country']+'/'+infos[i]['division']]:
            if location not in locations:
                locations.append(location)
        dates.append(infos[i]['date'])
        pango_lineages.append(infos[i]['pangoLineage'])
        nextstrain_clades.append(infos[i]['nextstrainClade'])
        gisaid_clades.append(infos[i]['gisaidClade'])

    # location
    locations.sort()
    with open(os.path.join(args.output, 'locations.txt'), 'w') as f:
        f.write('Locations'+'\n')
        f.write('Global'+'\n')
        for l in locations:
            f.write(l+'\n')

    # dates
    dates.sort()
    dates = [datetime.strftime(x,'%Y-%m-%d') for x in list(pd.date_range(start=dates[0], end=dates[-1]))]
    with open(os.path.join(args.output, 'dates.txt'), 'w') as f:
        f.write('Dates'+'\n')
        for d in dates:
            f.write(d+'\n')
    
    # pango lineage
    pango_lineages = list(set(pango_lineages))
    pango_lineages.sort()
    with open(os.path.join(args.output, 'pango_lineages.txt'), 'w') as f:
        f.write('Pango_lineages'+'\n')
        for l in pango_lineages:
            f.write(l+'\n')

    # nextstrain clade
    nextstrain_clades = list(set(nextstrain_clades))
    if 'nan' in nextstrain_clades:
        nextstrain_clades.remove('nan')
    nextstrain_clades.sort()
    with open(os.path.join(args.output, 'nextstrain_clades.txt'), 'w') as f:
        f.write('Nextstrain_clades'+'\n')
        for nc in nextstrain_clades:
            f.write(nc+'\n')

    # gisaid clade
    gisaid_clades = list(set(gisaid_clades))
    if 'nan' in gisaid_clades:
        gisaid_clades.remove('nan')
    gisaid_clades.sort()
    with open(os.path.join(args.output, 'gisaid_clades.txt'), 'w') as f:
        f.write('Gisaid_clades'+'\n')
        for gc in gisaid_clades:
            f.write(gc+'\n')

    # variants of concern and variants of interest (WHO)
    who_variants = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Omicron', 'Lambda', 'Mu']
    with open(os.path.join(args.output, 'who_variants.txt'), 'w') as f:
        f.write('WHO_variants'+'\n')
        for v in who_variants:
            f.write(v+'\n')

    # amino acid substitutions
    amino_acid = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S', 'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'stop', 'del']
    df = pd.read_csv(args.genome, index_col=0).drop_duplicates(subset=['product', 'aaPos'])
    sarscov2_AA = []
    for i in df.index:
        for aa in amino_acid:
            sarscov2_aa = df.loc[i, 'product'] + '_' + str(df.loc[i, 'aaPos']) + aa
            sarscov2_AA.append(sarscov2_aa)
    with open(os.path.join(args.output, 'amino_acid.txt'), 'w') as f:
        f.write('AA'+'\n')
        for a in sarscov2_AA:
            f.write(a+'\n')

    # nucleotide substitutions
    GENOME_LENTH = 29891
    nucleotide = ['A', 'T', 'G', 'C']
    sarscov2_nt = []
    for gsite in range(1, GENOME_LENTH+1):
        for nt in nucleotide:
            sarscov2_nt.append(str(gsite)+nt)
    with open(os.path.join(args.output, 'nucleotide.txt'), 'w') as f:
        f.write('Nucleotide'+'\n')
        for n in sarscov2_nt:
            f.write(n+'\n')


if __name__ == '__main__':
    main()
