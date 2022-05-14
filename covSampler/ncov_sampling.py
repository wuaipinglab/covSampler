'''
subsampling
'''

import os
import re
import time
import argparse
import pandas as pd


def get_target_ids(file, args, genome_path):
    '''
    get accession id of sequences in selected range (location, date and variants)

    :param file: str, infos file path
    :param args: argparse, args
    :param genome_path: str, SARS-CoV-2 genome file path

    :return: dict, key: accession id; value: infos of sequences
             list, accession id of sequences in selected range
    '''
    # read infos from infos.tsv
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
    
    target_ids = []
    # add accession id of sequences in selected location
    if args.location == 'Global':
        target_ids = list(infos.keys())
    elif args.location.count('/') == 0:
        for i in infos:
            if infos[i]['region'] == args.location:
                target_ids.append(i)
    elif args.location.count('/') == 1:
        for i in infos:
            if infos[i]['region'] == args.location.split('/')[0] and infos[i]['country'] == args.location.split('/')[1]:
                target_ids.append(i)
    elif args.location.count('/') == 2:
        for i in infos:
            if infos[i]['region'] == args.location.split('/')[0] and infos[i]['country'] == args.location.split('/')[1] and infos[i]['division'] == args.location.split('/')[2]:
                target_ids.append(i)
    
    # remove accession id of sequences not in selected date range
    dropDate = []
    for i in target_ids:
        if not args.dateStart <= infos[i]['date'] <= args.dateEnd:
            dropDate.append(i)
    
    target_ids = list(set(target_ids)-set(dropDate))
    
    # remove accession id of sequences not belong to selected variant(s)
    if args.variants is not None:
        for variant in args.variants:
            dropVariant = []
            # lineage
            if variant.startswith('Lineage'):
                # VOC or VOI (example: Lineage/WHO/Alpha)
                if variant.split('/')[1] == 'WHO':
                    PANGO_WHO = {
                        'B.1.1.7': 'Alpha',
                        'Q': 'Alpha',
                        'B.1.351': 'Beta',
                        'P.1': 'Gamma',
                        'B.1.617.2': 'Delta',
                        'AY': 'Delta',
                        'B.1.1.529': 'Omicron',
                        'BA': 'Omicron',
                        'C.37': 'Lambda',
                        'B.1.621': 'Mu'
                    }
                    for i in target_ids:
                        isVar = False
                        for l in PANGO_WHO:
                            if infos[i]['pangoLineage'] == l or infos[i]['pangoLineage'].startswith(l+'.'):
                                lWHO = PANGO_WHO[l]
                                if lWHO == variant.split('/')[2]:
                                    isVar = True
                        if not isVar:
                            dropVariant.append(i)
                # pango lineage (example: Lineage/Pango_lineage/B.1.1.7)
                elif variant.split('/')[1] == 'Pango_lineage':
                    for i in target_ids:
                        if infos[i]['pangoLineage'] != variant.split('/')[2]:
                            dropVariant.append(i)
                # nextstrain clade (example: Lineage/Nextstrain_clade/20I (Alpha, V1))
                elif variant.split('/')[1] == 'Nextstrain_clade':
                    for i in target_ids:
                        if infos[i]['nextstrainClade'] != variant.split('/')[2]:
                            dropVariant.append(i)
                # gisaid clade (example: Lineage/Gisaid_clade/G)
                elif variant.split('/')[1] == 'Gisaid_clade':
                    for i in target_ids:
                        if infos[i]['gisaidClade'] != variant.split('/')[2]:
                            dropVariant.append(i)
            # site
            elif variant.startswith('Site'):
                genome = pd.read_csv(genome_path, index_col=0)
                # nucleotide site
                if variant.split('/')[1] == 'Nucleotide':
                    genomePos = int(variant.split('/')[2][:-1])
                    ref = genome[genome['genomePos'] == genomePos]['nucleotide'].values[0]
                    # ref (example: Site/Nucleotide/23403A)
                    if variant.split('/')[2][-1] == ref:
                        for i in target_ids:
                            if re.search('[^0-9]'+str(genomePos)+'[^0-9]', infos[i]['nt']):
                                dropVariant.append(i)
                    # substitution (example: Site/Nucleotide/23403G)
                    else:
                        for i in target_ids:
                            if not re.search('[^0-9]'+variant.split('/')[2], infos[i]['nt']):
                                dropVariant.append(i)
                # amino acid site
                elif variant.split('/')[1] == 'Amino_acid':
                    product = variant.split('/')[2].split('_')[0]
                    if not (variant.split('/')[2].endswith('del') or variant.split('/')[2].endswith('stop')):
                        aaPos = int(variant.split('/')[2].split('_')[1][:-1])
                        ref = genome[(genome['product'] == product) & (genome['aaPos'] == aaPos)]['aa'].values[0]
                        # ref (example: Site/Amino_acid/Spike_614D)
                        if variant.split('/')[2][-1] == ref:
                            for i in target_ids:
                                if re.search(product+'_'+ref+str(aaPos)+'[^0-9]', infos[i]['aa']):
                                    dropVariant.append(i)
                        # substitution (example: Site/Amino_acid/Spike_614G)
                        else:
                            for i in target_ids:
                                if product+'_'+ref+variant.split('/')[2].split('_')[1] not in infos[i]['aa']:
                                    dropVariant.append(i)
                    else:
                        # deletion (example: Site/Amino_acid/Spike_69del)
                        if variant.split('/')[2].endswith('del'):
                            aaPos = int(variant.split('/')[2].split('_')[1][:-3])
                        # stop (example: Site/Amino_acid/NS8_Q27stop)
                        elif variant.split('/')[2].endswith('stop'):
                            aaPos = int(variant.split('/')[2].split('_')[1][:-4])
                        ref = genome[(genome['product'] == product) & (genome['aaPos'] == aaPos)]['aa'].values[0]
                        for i in target_ids:
                            if product+'_'+ref+variant.split('/')[2].split('_')[1] not in infos[i]['aa']:
                                dropVariant.append(i)

            target_ids = list(set(target_ids)-set(dropVariant))
    
    target_ids = sorted(target_ids, key=lambda x:int(x.split('_')[2]))
    
    return infos, target_ids
            

def read_haplotype_sequence(file, target_ids):
    '''
    read haplotype sequences (constructed by combining pre-calculated SARS-CoV-2
    key sites) of sequences in selected range from haplotype_sequence.txt

    :param file: str, haplotype sequence file path
    :param target_ids: list, accession id of sequences in selected range

    :return: dict, key: accession id; value: snps at key sites
    '''
    all_seqs = {}
    seqs = {}
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if line.strip().split(':')[1] == '':
                all_seqs[line.split(':')[0]] = ' '
            else:
                all_seqs[line.split(':')[0]] = line.strip().split(':')[1]
    for i in target_ids:
        seqs[i] = all_seqs[i]
    return seqs


def read_path(file, target_ids, infos):
    '''
    read divergent pathways of sequences in selected range from divergent_pathway.csv
    divide divergent pathways by continent

    :param file: str, divergent pathways file path
    :param target_ids: list, accession id of sequences in selected range
    :param infos: dict, key: accession id; value: infos of sequences
    
    :return: dict, key: continent; value: divergent pathways
    '''
    paths_in_continent = {}
    seq_path = {}
    target_seq_path = {}
    target_path = {}
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('EPI_ISL'):
                seq_path[line.split(',')[0]] = line.strip().split(',')[1]
        for i in target_ids:
            target_seq_path[i] = seq_path[i]
        for i in target_seq_path:
            target_path.setdefault(target_seq_path[i], []).append(i)
        for p in list(target_path.values()):
            continent = infos[p[0]]['region']
            paths_in_continent.setdefault(continent, []).append(p)
    for c in paths_in_continent:
        paths_in_continent[c] = sorted(paths_in_continent[c], key=lambda x: len(x), reverse=True)

    return paths_in_continent


def calculate_continent_sample_number(required_sample_num, seqs, infos):
    '''
    distribute the number of subsamples to each continent

    :param required_sample_num: int, number of subsamples
    :param seqs: dict, key: accession id; value: snps at key sites
    :param infos: dict, key: accession id; value: infos of sequences
    
    :return: dict, key: continent; value: number of subsamples in the continent
    '''
    s = required_sample_num
    continent_genome_num = {}
    for i in seqs:
        continent = infos[i]['region']
        continent_genome_num[continent] = continent_genome_num.get(continent, 0) + 1

    continent_genome_num = dict(sorted(continent_genome_num.items(), key=lambda x: x[1], reverse=True))
    
    continent_sample_number = {}
    for i in continent_genome_num:
        continent_sample_number[i] = 0
    while s > 0:
        for i in continent_sample_number:
            if s == 0:
                break
            if continent_sample_number[i] < continent_genome_num[i]:
                continent_sample_number[i] += 1
                s -= 1
    return continent_sample_number


def com_sampling(required_sample_num, seqs, infos, paths_in_continent):
    '''
    comprehensive subsampling

    :param required_sample_num: int, number of subsamples
    :param seqs: dict, key: accession id; value: snps at key sites
    :param infos: dict, key: accession id; value: infos of sequences
    :param paths_in_continent: dict, key: continent; value: divergent pathways

    :return: list, accession id of subsamples
    '''
    continent_sample_number = calculate_continent_sample_number(required_sample_num, seqs, infos)
    com_samples = []
    # sampling in each continent
    for continent in paths_in_continent:
        sample_number_in_continent = continent_sample_number[continent]
        # distribute the number of subsamples (in the continent) to each divergent pathways
        sample_number_in_paths = []
        for i in range(0, len(paths_in_continent[continent])):
            sample_number_in_paths.append(0)
        while sample_number_in_continent > 0:
            for i in range(0, len(sample_number_in_paths)):
                if sample_number_in_continent == 0:
                    break
                if sample_number_in_paths[i] < len(paths_in_continent[continent][i]):
                    sample_number_in_paths[i] += 1
                    sample_number_in_continent -= 1
        # sampling in each divergent pathway
        for i in range(0, len(paths_in_continent[continent])):
            path = paths_in_continent[continent][i]
            # distribute the number of subsamples (in the divergent pathway) to each month
            seq_in_months = {}
            for seq_id in path:
                month = infos[seq_id]['date'][:7]
                seq_in_months.setdefault(month, []).append(seq_id)
            seq_in_months = dict(sorted(seq_in_months.items(), key=lambda x: x[0], reverse=True))
            
            sample_number_in_path = sample_number_in_paths[i]
            sample_number_in_months = {}

            for m in seq_in_months:
                sample_number_in_months[m] = 0

            while sample_number_in_path > 0:
                for m in seq_in_months:
                    if sample_number_in_path == 0:
                        break
                    if sample_number_in_months[m] < len(seq_in_months[m]):
                        sample_number_in_months[m] += 1
                        sample_number_in_path -= 1
            # sampling in each month
            for m in seq_in_months:
                seq_in_month = seq_in_months[m]
                # distribute the number of subsamples (in the month) to each haplotype
                seq_in_haplotype_sequences = {}
                for seq_id in seq_in_month:
                    haplotype_sequence = seqs[seq_id]
                    seq_in_haplotype_sequences.setdefault(haplotype_sequence, []).append(seq_id)

                seq_in_haplotype_sequences = dict(sorted(seq_in_haplotype_sequences.items(),
                                                       key=lambda x: len(x[0].split(',')), reverse=True))

                sample_number_in_month = sample_number_in_months[m]
                sample_number_in_haplotype_sequences = {}

                for v in seq_in_haplotype_sequences:
                    sample_number_in_haplotype_sequences[v] = 0

                while sample_number_in_month > 0:
                    for v in seq_in_haplotype_sequences:
                        if sample_number_in_month == 0:
                            break
                        if sample_number_in_haplotype_sequences[v] < len(seq_in_haplotype_sequences[v]):
                            sample_number_in_haplotype_sequences[v] += 1
                            sample_number_in_month -= 1
                # sampling in each haplotype
                for v in seq_in_haplotype_sequences:
                    seq_in_haplotype_sequence = seq_in_haplotype_sequences[v]
                    sample_number_in_haplotype_sequence = sample_number_in_haplotype_sequences[v]
                    seq_in_haplotype_sequence = list(sorted(seq_in_haplotype_sequence,
                                                          key=lambda x: int(x.split('_')[2]), reverse=True))
                    com_samples.extend(seq_in_haplotype_sequence[:sample_number_in_haplotype_sequence])
    
    com_samples = sorted(com_samples, key=lambda x: int(x.split('_')[2]))
    
    return com_samples


def calculate_continent_genome_num_proportion(seqs, infos):
    '''
    calculate the proportion and number of sequences in each continent
    
    :param seqs: dict, key: accession id; value: snps at key sites
    :param infos: dict, key: accession id; value: infos of sequences

    :return: dict, key: continent; value: proportion and number of sequences in the continent
    '''
    continent_genome_num = {}
    for i in seqs:
        continent = infos[i]['region']
        continent_genome_num[continent] = continent_genome_num.get(continent, 0) + 1

    continent_genome_num_proportion = {}
    for i in continent_genome_num:
        continent_genome_num_proportion[i] = {'proportion': continent_genome_num[i] / len(seqs), 'count': continent_genome_num[i]}
    return continent_genome_num_proportion


def calculate_continent_threshold(required_sample_num, continent_genome_num_proportion, paths_in_continent):
    '''
    1. distribute the number of subsamples to each continent
       -> make sure the proportion of subsamples in each continent is
          consistent with the proportion of sequences in selected range
          in each continent
    2. calculate the threshold of each continent
       in each continent,
       -> the number of subsamples in each divergent pathway = the number of
          all sequences in seleted range in the divergent pathway / threshold
       -> the threshold starts with 0 and increments by 1,
          until the number of subsamples in all divergent pathways <= the
          number of subsamples in the continent (distributed in step 1)
    
    :param required_sample_num: int, number of subsamples
    :param continent_genome_num_proportion:
           dict, key: continent;
                 value: proportion and number of sequences in the continent
    :param paths_in_continent: dict, key: continent; value: divergent pathways

    :return: dict, key: continent;
                   value: threshold and corresponding number of subsamples in
                          all divergent pathways in the continent
    '''
    continent_threshold = {}
    for continent in paths_in_continent:
        proportion = continent_genome_num_proportion[continent]['proportion']
        required_sample_num_in_continent = required_sample_num * proportion

        threshold = 1
        while True:
            actual_sample_num_in_continent = 0
            for path in paths_in_continent[continent]:
                actual_sample_num_in_continent += int(len(path)/threshold)
                if actual_sample_num_in_continent > required_sample_num_in_continent:
                    break
            if actual_sample_num_in_continent <= required_sample_num_in_continent:
                break
            else:
                threshold += 1
        continent_threshold[continent] = {'threshold': threshold, 'actual_num': actual_sample_num_in_continent}
    
    return continent_threshold


def calculate_continent_extra_sample_num(required_sample_num, continent_genome_num_proportion, continent_threshold):
    '''
    calculate the number of extra subsamples of each continent after threshold calculation

    :param required_sample_num: int, number of subsamples
    :param continent_genome_num_proportion: 
           dict, key: continent;
                 value: proportion and number of sequences in the continent
    :param continent_threshold: 
           dict, key: continent;
                 value: threshold and corresponding number of subsamples in
                        all divergent pathways in the continent
    
    :return: dict, key: continent;
                   value: the number of extra subsamples of the continent
                          after threshold calculation
    '''
    continent_extra_sample_num = {}
    sample_number_with_threshold = 0
    for i in continent_threshold:
        sample_number_with_threshold += continent_threshold[i]['actual_num']

    extra_sample_num = required_sample_num - sample_number_with_threshold

    continent_list_sorted = [i[0] for i in sorted(continent_threshold.items(), key=lambda x: x[1]['actual_num'])]

    for i in continent_list_sorted:
        continent_extra_sample_num[i] = 0

    while extra_sample_num > 0:
        for i in continent_list_sorted:
            if extra_sample_num == 0:
                break
            if continent_threshold[i]['actual_num'] + continent_extra_sample_num[i] < continent_genome_num_proportion[i]['count']:
                continent_extra_sample_num[i] += 1
                extra_sample_num -= 1

    return continent_extra_sample_num


def rep_sampling(required_sample_num, seqs, infos, paths_in_continent):
    '''
    representative subsampling

    :param required_sample_num: int, number of subsamples
    :param seqs: dict, key: accession id; value: snps at key sites
    :param infos: dict, key: accession id; value: infos of sequences
    :param paths_in_continent: dict, key: continent; value: divergent pathways

    :return: list, accession id of subsamples
    '''
    continent_genome_num_proportion = calculate_continent_genome_num_proportion(seqs, infos)
    continent_threshold = calculate_continent_threshold(required_sample_num, continent_genome_num_proportion, paths_in_continent)
    continent_extra_sample_num = calculate_continent_extra_sample_num(required_sample_num, continent_genome_num_proportion, continent_threshold)

    rep_samples = []
    # sampling in each continent
    for continent in paths_in_continent:
        # distribute the number of subsamples (in the continent) to each divergent pathways
        sample_number_in_paths = []
        for path in paths_in_continent[continent]:
            sample_number_in_paths.append(int(len(path)/continent_threshold[continent]['threshold']))
        extra_sample_number_in_continent = continent_extra_sample_num[continent]
        while extra_sample_number_in_continent > 0:
            for i in range(0, len(sample_number_in_paths)):
                if extra_sample_number_in_continent == 0:
                    break
                if sample_number_in_paths[i] < len(paths_in_continent[continent][i]):
                    sample_number_in_paths[i] += 1
                    extra_sample_number_in_continent -= 1
        # sampling in each divergent pathway
        for i in range(0, len(paths_in_continent[continent])):
            path = paths_in_continent[continent][i]
            # distribute the number of subsamples (in the divergent pathway) to each month
            seq_in_months = {}
            for seq_id in path:
                month = infos[seq_id]['date'][:7]
                seq_in_months.setdefault(month, []).append(seq_id)
            seq_in_months = dict(sorted(seq_in_months.items(), key=lambda x: x[0], reverse=True))

            sample_number_in_path = sample_number_in_paths[i]
            sample_number_in_months = {}

            for m in seq_in_months:
                sample_number_in_months[m] = 0

            while sample_number_in_path > 0:
                for m in seq_in_months:
                    if sample_number_in_path == 0:
                        break
                    if sample_number_in_months[m] < len(seq_in_months[m]):
                        sample_number_in_months[m] += 1
                        sample_number_in_path -= 1
            # sampling in each month
            for m in seq_in_months:
                seq_in_month = seq_in_months[m]
                # distribute the number of subsamples (in the month) to each haplotype
                seq_in_haplotype_sequences = {}
                for seq_id in seq_in_month:
                    haplotype_sequence = seqs[seq_id]
                    seq_in_haplotype_sequences.setdefault(haplotype_sequence, []).append(seq_id)

                seq_in_haplotype_sequences = dict(sorted(seq_in_haplotype_sequences.items(),
                                                       key=lambda x: len(x[1]), reverse=True))

                sample_number_in_month = sample_number_in_months[m]
                sample_number_in_haplotype_sequences = {}

                for v in seq_in_haplotype_sequences:
                    sample_number_in_haplotype_sequences[v] = 0

                while sample_number_in_month > 0:
                    for v in seq_in_haplotype_sequences:
                        if sample_number_in_month == 0:
                            break
                        if sample_number_in_haplotype_sequences[v] < len(seq_in_haplotype_sequences[v]):
                            sample_number_in_haplotype_sequences[v] += 1
                            sample_number_in_month -= 1
                # sampling in each haplotype
                for v in seq_in_haplotype_sequences:
                    seq_in_haplotype_sequence = seq_in_haplotype_sequences[v]
                    sample_number_in_haplotype_sequence = sample_number_in_haplotype_sequences[v]
                    seq_in_haplotype_sequence = list(sorted(seq_in_haplotype_sequence,
                                                          key=lambda x: int(x.split('_')[2]), reverse=True))
                    rep_samples.extend(seq_in_haplotype_sequence[:sample_number_in_haplotype_sequence])
    
    rep_samples = sorted(rep_samples, key=lambda x: int(x.split('_')[2]))
    
    return rep_samples


def write_new_file(file, samples, DATE, args, target_ids):
    with open(file, 'w') as f:
        f.write('## File time: '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+' (Beijing Time)'+'\n')
        f.write('## Data version: '+DATE+'\n')
        f.write('## Location of samples: '+args.location+'\n')
        f.write('## Start date of samples: '+args.dateStart+'\n')
        f.write('## End date of samples: '+args.dateEnd+'\n')
        if args.variants is not None:
            for variant in args.variants:
                f.write('## Variant of samples: '+variant.replace('//','')+'\n') # .replace('//','') -> front end issue
        f.write('## Number of all samples in range: '+str(len(target_ids))+'\n')
        f.write('## Sampling characteristic: '+args.characteristic+'\n')
        if args.size > len(target_ids):
            f.write('## Required sample size: '+str(args.size)+'\n')
            f.write('## WARNING!!! Required sample size > number of all samples in range'+'\n')
            f.write('## Actual sample size: '+str(len(target_ids))+'\n')
        else:
            f.write('## Sample size: '+str(args.size)+'\n')
        f.write('# ID'+'\n')
        for i in samples:
            f.write(i+'\n')


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='A subsampling method for large-scale SARS-CoV-2 genomes')
    parser.add_argument('--dirpath', required=True, help='Data directory path')
    parser.add_argument('--location', required=True, help='Location of subsamples')
    parser.add_argument('--dateStart', required=True, help='Start date of subsamples')
    parser.add_argument('--dateEnd', required=True, help='End date of subsamples')
    parser.add_argument('--variants', action='append', help='Variants of subsamples')
    parser.add_argument('--size', type=int, required=True, help='Number of subsamples')
    parser.add_argument('--characteristic', required=True, help='Characteristic of subsampling')
    parser.add_argument('--output', required=True, help='Output directory path')
    args = parser.parse_args()

    # data version
    DATE = '2022-05-05'

    # files required for subsampling
    DIRPATH = args.dirpath
    genome_path = os.path.join(DIRPATH, 'rawdata', 'SARS_CoV_2.csv')
    seq_info_path = os.path.join(DIRPATH, 'infos.tsv')
    haplotype_sequence_path = os.path.join(DIRPATH, 'haplotype_sequence.txt')
    divergent_pathway_path = os.path.join(DIRPATH, 'divergent_pathway.csv')

    infos, target_ids = get_target_ids(seq_info_path, args, genome_path)
    
    # compare number of all sequences in selected range and number of required subsamples
    if args.size > len(target_ids):
        required_sample_num = len(target_ids)
    else:
        required_sample_num = args.size

    seqs = read_haplotype_sequence(haplotype_sequence_path, target_ids)

    paths_in_continent = read_path(divergent_pathway_path, target_ids, infos)

    # characteristic of subsampling
    if args.characteristic == 'Comprehensive':
        samples = com_sampling(required_sample_num, seqs, infos, paths_in_continent)
    elif args.characteristic == 'Representative':
        samples = rep_sampling(required_sample_num, seqs, infos, paths_in_continent)

    write_new_file(os.path.join(args.output, 'samples.txt'), samples, DATE, args, target_ids)


if __name__ == '__main__':
    main()    
