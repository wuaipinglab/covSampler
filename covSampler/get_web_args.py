'''
prepare the content of each parameter in the web interface
'''

import os
from datetime import datetime
import pandas as pd
from covSampler import genome_path, infos_path, args_path, arg_location_path, arg_date_path, arg_pangoLineage_path, arg_nextstrainClade_path, arg_gisaidClade_path, arg_WHO_path, arg_aa_path, arg_nt_path


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
    if not os.path.exists(args_path):
        os.makedirs(args_path)

    infos = read_infos(infos_path)
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
    with open(arg_location_path, 'w') as f:
        f.write('Locations'+'\n')
        f.write('Global'+'\n')
        for l in locations:
            f.write(l+'\n')

    # dates
    dates.sort()
    dates = [datetime.strftime(x,'%Y-%m-%d') for x in list(pd.date_range(start=dates[0], end=dates[-1]))]
    with open(arg_date_path, 'w') as f:
        f.write('Dates'+'\n')
        for d in dates:
            f.write(d+'\n')
    
    # pango lineage
    pango_lineages = list(set(pango_lineages))
    pango_lineages.sort()
    with open(arg_pangoLineage_path, 'w') as f:
        f.write('Pango_lineages'+'\n')
        for l in pango_lineages:
            f.write(l+'\n')

    # nextstrain clade
    nextstrain_clades = list(set(nextstrain_clades))
    if 'nan' in nextstrain_clades:
        nextstrain_clades.remove('nan')
    nextstrain_clades.sort()
    with open(arg_nextstrainClade_path, 'w') as f:
        f.write('Nextstrain_clades'+'\n')
        for nc in nextstrain_clades:
            f.write(nc+'\n')

    # gisaid clade
    gisaid_clades = list(set(gisaid_clades))
    if 'nan' in gisaid_clades:
        gisaid_clades.remove('nan')
    gisaid_clades.sort()
    with open(arg_gisaidClade_path, 'w') as f:
        f.write('Gisaid_clades'+'\n')
        for gc in gisaid_clades:
            f.write(gc+'\n')

    # variants of concern and variants of interest (WHO)
    who_variants = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Omicron', 'Lambda', 'Mu']
    with open(arg_WHO_path, 'w') as f:
        f.write('WHO_variants'+'\n')
        for v in who_variants:
            f.write(v+'\n')

    # amino acid substitutions
    amino_acid = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S', 'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'stop', 'del']
    df = pd.read_csv(genome_path, index_col=0).drop_duplicates(subset=['product', 'aaPos'])
    sarscov2_AA = []
    for i in df.index:
        for aa in amino_acid:
            sarscov2_aa = df.loc[i, 'product'] + '_' + str(df.loc[i, 'aaPos']) + aa
            sarscov2_AA.append(sarscov2_aa)
    with open(arg_aa_path, 'w') as f:
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
    with open(arg_nt_path, 'w') as f:
        f.write('Nucleotide'+'\n')
        for n in sarscov2_nt:
            f.write(n+'\n')


if __name__ == '__main__':
    main()
