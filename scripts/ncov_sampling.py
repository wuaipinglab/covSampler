'''
subsampling
'''

import os
import time
import argparse


def get_target_ids(args):
    '''
    get strains in selected range (location, date and variants)

    :param args: argparse, args

    :return: dict, key: strain; value: infos
             list, strains in selected range
    '''
    # get infos
    infos = {}
    with open(args.infos) as f:
        lines = f.readlines()
        for n, line in enumerate(lines):
            if n != 0:
                seq_items = line.strip().split('\t')
                seq_id = seq_items[0]
                seq_region = seq_items[1]
                seq_country = seq_items[2]
                seq_division = seq_items[3]
                seq_date = seq_items[4]
                seq_pangolineage = seq_items[5]
                seq_nextstrainclade = seq_items[6]
                seq_substitutions = seq_items[7]
                seq_aasubstitutions = seq_items[8]
                seq_aadeletions = seq_items[9]
                infos[seq_id] = {
                    'region': seq_region,
                    'country': seq_country,
                    'division': seq_division,
                    'date': seq_date,
                    'pangoLineage': seq_pangolineage,
                    'nextstrainClade': seq_nextstrainclade,
                    'substitutions': seq_substitutions,
                    'aaSubstitutions': seq_aasubstitutions,
                    'aaDeletions': seq_aadeletions
                }
                if args.return_genbank_accession:
                    seq_accession = seq_items[10]
                    infos[seq_id]['genbankAccession'] = seq_accession

    target_ids = []

    # add strains in selected location
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
    
    # remove strains not in selected date range
    to_rm_date = []
    for i in target_ids:
        if not args.date_start <= infos[i]['date'] <= args.date_end:
            to_rm_date.append(i)
    
    target_ids = list(set(target_ids)-set(to_rm_date))
    
    # remove strains != selected variant(s)
    if args.variants is not None:
        for variant in args.variants:
            to_rm_variant = []
            # lineage
            if variant.startswith('Lineage'):
                # variants of concern or variants of interest (example: Lineage/WHO/Alpha)
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
                        'BC': 'Omicron',
                        'BD': 'Omicron',
                        'BE': 'Omicron',
                        'BF': 'Omicron',
                        'BG': 'Omicron',
                        'BH': 'Omicron',
                        'BJ': 'Omicron',
                        'BK': 'Omicron',
                        'BL': 'Omicron',
                        'BM': 'Omicron',
                        'BN': 'Omicron',
                        'BP': 'Omicron',
                        'BQ': 'Omicron',
                        'BR': 'Omicron',
                        'BS': 'Omicron',
                        'BT': 'Omicron',
                        'BU': 'Omicron',
                        'BV': 'Omicron',
                        'BW': 'Omicron',
                        'BY': 'Omicron',
                        'BZ': 'Omicron',
                        'CA': 'Omicron',
                        'CB': 'Omicron',
                        'CC': 'Omicron',
                        'CD': 'Omicron',
                        'CE': 'Omicron',
                        'CF': 'Omicron',
                        'CG': 'Omicron',
                        'CH': 'Omicron',
                        'CJ': 'Omicron',
                        'CK': 'Omicron',
                        'CL': 'Omicron',
                        'CM': 'Omicron',
                        'CN': 'Omicron',
                        'CP': 'Omicron',
                        'CQ': 'Omicron',
                        'CR': 'Omicron',
                        'CS': 'Omicron',
                        'CT': 'Omicron',
                        'CU': 'Omicron',
                        'CV': 'Omicron',
                        'CW': 'Omicron',
                        'CY': 'Omicron',
                        'CZ': 'Omicron',
                        'DA': 'Omicron',
                        'DB': 'Omicron',
                        'DC': 'Omicron',
                        'DD': 'Omicron',
                        'DE': 'Omicron',
                        'DF': 'Omicron',
                        'DG': 'Omicron',
                        'DH': 'Omicron',
                        'DJ': 'Omicron',
                        'DK': 'Omicron',
                        'DL': 'Omicron',
                        'DM': 'Omicron',
                        'DN': 'Omicron',
                        'DP': 'Omicron',
                        'B.1.427': 'Epsilon',
                        'B.1.429': 'Epsilon',
                        'P.2': 'Zeta',
                        'B.1.525': 'Eta',
                        'P.3': 'Theta',
                        'B.1.526': 'Iota',
                        'B.1.617.1': 'Kappa',
                        'C.37': 'Lambda',
                        'B.1.621': 'Mu',
                        'BB': 'Mu'
                    }
                    for i in target_ids:
                        is_who_var = False
                        for l in PANGO_WHO:
                            if infos[i]['pangoLineage'] == l or infos[i]['pangoLineage'].startswith(l+'.'):
                                l_who = PANGO_WHO[l]
                                if l_who == variant.split('/')[2]:
                                    is_who_var = True
                        if not is_who_var:
                            to_rm_variant.append(i)
                # pango lineage (example: Lineage/Pango_lineage/B.1.1.7)
                elif variant.split('/')[1] == 'Pango_lineage':
                    for i in target_ids:
                        if infos[i]['pangoLineage'] != variant.split('/')[2]:
                            to_rm_variant.append(i)
                # nextstrain clade (example: Lineage/Nextstrain_clade/20I (Alpha, V1))
                elif variant.split('/')[1] == 'Nextstrain_clade':
                    for i in target_ids:
                        if infos[i]['nextstrainClade'] != variant.split('/')[2]:
                            to_rm_variant.append(i)
            # site
            elif variant.startswith('Site'):
                # nucleotide substitution (example: Site/Nucleotide/A23403G)
                if variant.split('/')[1] == 'Nucleotide':
                    for i in target_ids:
                        if variant.split('/')[2] not in infos[i]['substitutions']:
                            to_rm_variant.append(i)
                # amino acid
                elif variant.split('/')[1] == 'Amino_acid':
                    # amino acid substitution (example: Site/Amino_acid/S:D614G)
                    if variant[-1] != '-':
                        for i in target_ids:
                            if variant.split('/')[2] not in infos[i]['aaSubstitutions']:
                                to_rm_variant.append(i)
                    # amino acid deletion (example: Site/Amino_acid/S:H69-)
                    else:
                        for i in target_ids:
                            if variant.split('/')[2] not in infos[i]['aaDeletions']:
                                to_rm_variant.append(i)

            target_ids = list(set(target_ids)-set(to_rm_variant))
    
    target_ids.sort()
    
    return infos, target_ids


def get_temporally_even_distribution(required_sample_num, target_ids, infos):
    '''
    divide strains in selected range by month and calculate required subsample number of each month (temporally even)
    
    :param required_sample_num: int, number of all subsamples
    :param target_ids: list, all strains in selected range
    :param infos: dict, key: strain; value: infos

    :return: dict, key: month; value: required subsample number and target ids
    '''
    temporally_even_distribution = {}
    # get target ids (each month)
    for i in target_ids:
        month = infos[i]['date'][:7]
        if month not in temporally_even_distribution:
            temporally_even_distribution[month] = {'required_sample_num': 0, 'target_ids': []}
        temporally_even_distribution[month]['target_ids'].append(i)
    
    temporally_even_distribution = dict(sorted(temporally_even_distribution.items(), key=lambda x: x[0], reverse=True))
    
    # get required sample number (each month)
    while required_sample_num > 0:
        for m in temporally_even_distribution:
            if required_sample_num == 0:
                break
            if temporally_even_distribution[m]['required_sample_num'] < len(temporally_even_distribution[m]['target_ids']):
                temporally_even_distribution[m]['required_sample_num'] += 1
                required_sample_num -= 1
    
    return temporally_even_distribution


def read_haplotype_sequence(file, target_ids):
    '''
    read haplotype sequences of strains in selected range

    :param file: str, haplotype sequence file path
    :param target_ids: list, strains in selected range

    :return: dict, key: strain; value: snps at key sites
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
    read divergent pathways of strains in selected range and divide these divergent pathways by continent

    :param file: str, divergent pathways file path
    :param target_ids: list, strains in selected range
    :param infos: dict, key: strain; value: infos
    
    :return: dict, key: continent; value: divergent pathways
    '''
    paths_in_continent = {}
    seq_path = {}
    target_seq_path = {}
    target_path = {}
    with open(file) as f:
        lines = f.readlines()
        for n, line in enumerate(lines):
            if n != 0:
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
    assign the number of subsamples to each continent

    :param required_sample_num: int, number of subsamples
    :param seqs: dict, key: strain; value: snps at key sites
    :param infos: dict, key: strain; value: infos
    
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
    :param seqs: dict, key: strain; value: snps at key sites
    :param infos: dict, key: strain; value: infos
    :param paths_in_continent: dict, key: continent; value: divergent pathways

    :return: list, subsamples
    '''
    continent_sample_number = calculate_continent_sample_number(required_sample_num, seqs, infos)
    com_samples = []
    # sampling in each continent
    for continent in paths_in_continent:
        sample_number_in_continent = continent_sample_number[continent]
        # assign the number of subsamples (in the continent) to each divergent pathway
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
            # assign the number of subsamples (in the divergent pathway) to each month
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
                # assign the number of subsamples (in the month) to each haplotype
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
                    sample_number_in_haplotype_sequence = sample_number_in_haplotype_sequences[v]
                    if sample_number_in_haplotype_sequence != 0:
                        seq_in_haplotype_sequence = seq_in_haplotype_sequences[v]
                        seq_in_haplotype_sequence.sort()
                        samples_in_haplotype_sequence = []
                        sampling_step = int(len(seq_in_haplotype_sequence)/sample_number_in_haplotype_sequence)
                        for s in range(0, len(seq_in_haplotype_sequence), sampling_step):
                            samples_in_haplotype_sequence.append(seq_in_haplotype_sequence[s])
                            if len(samples_in_haplotype_sequence) == sample_number_in_haplotype_sequence:
                                break
                        com_samples.extend(samples_in_haplotype_sequence)
    
    return com_samples


def calculate_continent_genome_num_proportion(seqs, infos):
    '''
    calculate the proportion and number of selected strains in each continent
    
    :param seqs: dict, key: strain; value: snps at key sites
    :param infos: dict, key: strain; value: infos

    :return: dict, key: continent; value: proportion and number of selected strains
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
    1. assign the number of subsamples to each continent
       -> make sure the proportion of subsamples in each continent is consistent
          with the proportion of strains in selected range in the continent
    2. calculate the threshold of each continent
       in each continent,
       -> the number of subsamples in each divergent pathway = the number of
          all strains in selected range in the divergent pathway / threshold
       -> the threshold starts with 0 and increments by 1,
          until the number of subsamples in all divergent pathways <= the
          number of subsamples in the continent (assigned in step 1)
    
    :param required_sample_num: int, number of subsamples
    :param continent_genome_num_proportion:
           dict, key: continent;
                 value: proportion and number of selected strains
    :param paths_in_continent: dict, key: continent; value: divergent pathways

    :return: dict, key: continent;
                   value: threshold and corresponding number of subsamples in all divergent pathways
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
                 value: proportion and number of selected strains
    :param continent_threshold: 
           dict, key: continent;
                 value: threshold and corresponding number of subsamples in all divergent pathways
    
    :return: dict, key: continent;
                   value: number of extra subsamples after threshold calculation
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
    :param seqs: dict, key: strain; value: snps at key sites
    :param infos: dict, key: strain; value: infos
    :param paths_in_continent: dict, key: continent; value: divergent pathways

    :return: list, subsamples
    '''
    continent_genome_num_proportion = calculate_continent_genome_num_proportion(seqs, infos)
    continent_threshold = calculate_continent_threshold(required_sample_num, continent_genome_num_proportion, paths_in_continent)
    continent_extra_sample_num = calculate_continent_extra_sample_num(required_sample_num, continent_genome_num_proportion, continent_threshold)

    rep_samples = []
    # sampling in each continent
    for continent in paths_in_continent:
        # assign the number of subsamples (in the continent) to each divergent pathway
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
            # assign the number of subsamples (in the divergent pathway) to each month
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
                # assign the number of subsamples (in the month) to each haplotype
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
                    sample_number_in_haplotype_sequence = sample_number_in_haplotype_sequences[v]
                    if sample_number_in_haplotype_sequence != 0:
                        seq_in_haplotype_sequence = seq_in_haplotype_sequences[v]
                        seq_in_haplotype_sequence.sort()
                        samples_in_haplotype_sequence = []
                        sampling_step = int(len(seq_in_haplotype_sequence)/sample_number_in_haplotype_sequence)
                        for s in range(0, len(seq_in_haplotype_sequence), sampling_step):
                            samples_in_haplotype_sequence.append(seq_in_haplotype_sequence[s])
                            if len(samples_in_haplotype_sequence) == sample_number_in_haplotype_sequence:
                                break
                        rep_samples.extend(samples_in_haplotype_sequence)
        
    return rep_samples


def get_version():
    about = {}
    thisdir = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(thisdir, '__version__.py')) as f:
        exec(f.read(), about)
    return 'v'+about['__version__']


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Subsampling')
    parser.add_argument('--infos', required=True, help='Infos file')
    parser.add_argument('--haplotypes', required=True, help='Haplotype sequences file')
    parser.add_argument('--divergent-pathways', required=True, help='Divergent pathways file')
    parser.add_argument('--description', help='Description recorded in the output file')
    parser.add_argument('--location', required=True, help='Location of subsamples')
    parser.add_argument('--date-start', required=True, help='Start date of subsamples')
    parser.add_argument('--date-end', required=True, help='End date of subsamples')
    parser.add_argument('--variants', action='append', help='Variants of subsamples')
    parser.add_argument('--size', type=int, required=True, help='Number of subsamples')
    parser.add_argument('--characteristic', required=True, help='Characteristic of subsampling')
    parser.add_argument('--temporally-even', action='store_true', help='Temporally even subsampling')
    parser.add_argument('--return-genbank-accession', action='store_true', help='Return genbank accession')
    parser.add_argument('--output', required=True, help='Subsamples file')
    args = parser.parse_args()

    # get infos and strains in selected range
    infos, target_ids = get_target_ids(args)

    # determine the number of subsamples
    required_sample_num = min(args.size, len(target_ids))

    if args.temporally_even:
        # get temporally even subsamples (monthly)
        temporally_even_distribution = get_temporally_even_distribution(required_sample_num, target_ids, infos)
        samples = []
        for m in temporally_even_distribution:
            seqs = read_haplotype_sequence(args.haplotypes, temporally_even_distribution[m]['target_ids'])
            paths_in_continent = read_path(args.divergent_pathways, temporally_even_distribution[m]['target_ids'], infos)
            # subsampling
            if args.characteristic == 'Comprehensive':
                samples_month = com_sampling(temporally_even_distribution[m]['required_sample_num'], seqs, infos, paths_in_continent)
            elif args.characteristic == 'Representative':
                samples_month = rep_sampling(temporally_even_distribution[m]['required_sample_num'], seqs, infos, paths_in_continent)
            samples.extend(samples_month)
    else:
        seqs = read_haplotype_sequence(args.haplotypes, target_ids)
        paths_in_continent = read_path(args.divergent_pathways, target_ids, infos)
        # subsampling
        if args.characteristic == 'Comprehensive':
            samples = com_sampling(required_sample_num, seqs, infos, paths_in_continent)
        elif args.characteristic == 'Representative':
            samples = rep_sampling(required_sample_num, seqs, infos, paths_in_continent)

    samples.sort()

    with open(args.output, 'w') as f:
        if args.description is not None:
            f.write('## Description: '+args.description+'\n')
        f.write('## covSampler version: '+get_version()+'\n')
        f.write('## File time: '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n')
        f.write('## Location of subsamples: '+args.location+'\n')
        f.write('## Start date of subsamples: '+args.date_start+'\n')
        f.write('## End date of subsamples: '+args.date_end+'\n')
        if args.variants is not None:
            for variant in args.variants:
                f.write('## Variant of subsamples: '+variant.replace('//','')+'\n') # .replace('//','') -> front end issue
        f.write('## Number of all sequences in range: '+str(len(target_ids))+'\n')
        if args.size > len(target_ids):
            f.write('## Required subsample size: '+str(args.size)+'\n')
            f.write('## WARNING!!! Required subsample size > number of all sequences in range'+'\n')
            f.write('## Actual subsample size: '+str(len(target_ids))+'\n')
        else:
            f.write('## Subsample size: '+str(args.size)+'\n')
        f.write('## Subsampling characteristic: '+args.characteristic+'\n')
        if args.temporally_even:
            f.write('## Temporally even: True'+'\n')
        else:
            f.write('## Temporally even: False'+'\n')
        f.write('# ID'+'\n')
        if args.return_genbank_accession:
            for i in samples:
                f.write(infos[i]['genbankAccession']+'\n')
        else:
            for i in samples:
                f.write(i+'\n')


if __name__ == '__main__':
    main()    
