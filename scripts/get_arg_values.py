'''
get all available arg values
'''

import os
import argparse
from datetime import datetime
import pandas as pd
from Bio import SeqIO


def main():
    # command line interface
    parser = argparse.ArgumentParser(description='Get all available arg values')
    parser.add_argument('--infos', required=True, help='Infos file')
    parser.add_argument('--reference-dir', required=True, help='Reference sequence directory')
    parser.add_argument('--output', required=True, help='Args directory')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

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
                seq_aasubsitutions = seq_items[8]
                seq_aadeletions = seq_items[9]
                infos[seq_id] = {
                    'region': seq_region,
                    'country': seq_country,
                    'division': seq_division,
                    'date': seq_date,
                    'pangoLineage': seq_pangolineage,
                    'nextstrainClade': seq_nextstrainclade,
                    'substitutions': seq_substitutions,
                    'aaSubsitutions': seq_aasubsitutions,
                    'aaDeletions': seq_aadeletions
                }

    locations = []
    dates = []
    pango_lineages = []
    nextstrain_clades = []
    for i in infos:
        region = infos[i]['region']
        country = infos[i]['region']+'/'+infos[i]['country']
        division = infos[i]['region']+'/'+infos[i]['country']+'/'+infos[i]['division']
        for location in [region, country, division]:
            if location not in locations:
                locations.append(location)
        dates.append(infos[i]['date'])
        pango_lineages.append(infos[i]['pangoLineage'])
        nextstrain_clades.append(infos[i]['nextstrainClade'])

    # locations
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
    
    # pango lineages
    pango_lineages = list(set(pango_lineages))
    pango_lineages.sort()
    with open(os.path.join(args.output, 'pango_lineages.txt'), 'w') as f:
        f.write('Pango_lineages'+'\n')
        for l in pango_lineages:
            f.write(l+'\n')

    # nextstrain clades
    nextstrain_clades = list(set(nextstrain_clades))
    nextstrain_clades.sort()
    with open(os.path.join(args.output, 'nextstrain_clades.txt'), 'w') as f:
        f.write('Nextstrain_clades'+'\n')
        for n in nextstrain_clades:
            f.write(n+'\n')

    # variants of concern and variants of interest (WHO)
    who_variants = [
        'Alpha', 'Beta', 'Gamma', 'Delta', 'Omicron', 'Epsilon', 'Zeta', 'Eta', 'Theta', 'Iota', 'Kappa', 'Lambda', 'Mu'
        ]
    with open(os.path.join(args.output, 'who_variants.txt'), 'w') as f:
        f.write('WHO_variants'+'\n')
        for v in who_variants:
            f.write(v+'\n')
    
    # nucleotide substitutions
    nucleotide = ['A', 'T', 'G', 'C']
    substitutions = []
    with open(os.path.join(args.reference_dir, 'sequence.fasta')) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            sequence = str(seq_record.seq)
            for site, ref_nuc in enumerate(sequence, 1):
                for nuc in nucleotide:
                    if nuc != ref_nuc:
                        substitutions.append(ref_nuc+str(site)+nuc)
                        
    with open(os.path.join(args.output, 'nucleotide.txt'), 'w') as f:
        f.write('Nucleotide'+'\n')
        for s in substitutions:
            f.write(s+'\n')

    # amino acid substitutions
    amino_acid = [
        'G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S', 'T', 'C', 'M', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', '*', '-'
        ]
    sc2_genes = ['ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'N']
    aa_substitutions = []
    for gene in sc2_genes:
        gene_fasta_file = 'sequence.gene.'+gene+'.fasta'
        with open(os.path.join(args.reference_dir, gene_fasta_file)) as f:
            for seq_record in SeqIO.parse(f, 'fasta'):
                gene_sequence = str(seq_record.seq)
                for site, ref_aa in enumerate(gene_sequence, 1):
                    for aa in amino_acid:
                        if aa != ref_aa:
                            aa_substitutions.append(gene+':'+ref_aa+str(site)+aa)

    with open(os.path.join(args.output, 'amino_acid.txt'), 'w') as f:
        f.write('AA'+'\n')
        for a in aa_substitutions:
            f.write(a+'\n')


if __name__ == '__main__':
    main()
