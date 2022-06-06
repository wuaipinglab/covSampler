'''
prepare region, country, division, date, pango lineage, nextstrain clade,
gisaid clade, nucleotide substitutions and amino acid substitutions of
each sequence for subsampling parameter selection
'''

import argparse
import warnings
import pandas as pd


def filter_seq(file):
    '''
    filter sequences in divergent pathways that contain the number of sequences <= 2
    
    :param file: str, divergent pathways file path

    :return: list, filtered accession id
    '''
    seqs_filtered = []
    pathways = {}
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith('ID'):
                seq_id = line.split(',')[0]
                pathway = line.strip().split(',')[1]
                pathways.setdefault(pathway, []).append(seq_id)
    for p in pathways:
        if len(pathways[p]) > 2:
            seqs_filtered.extend(pathways[p])
    seqs_filtered = sorted(seqs_filtered, key=lambda x: int(x.split('_')[2]))

    return seqs_filtered


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


def write_new_file(file, infos):
    with open(file, 'w') as f:
        f.write('ID'+'\t'+'region'+'\t'+'country'+'\t'+'division'+'\t'+'date'+'\t' +
                'pangoLineage'+'\t'+'nextstrainClade'+'\t'+'gisaidClade'+'\t'+'nt'+'\t'+'aa'+'\n')
        for i in infos:
            f.write(
                i+'\t'
                + infos[i]['region']+'\t'
                + infos[i]['country']+'\t'
                + infos[i]['division']+'\t'
                + infos[i]['date']+'\t'
                + infos[i]['pangoLineage']+'\t'
                + infos[i]['nextstrainClade']+'\t'
                + infos[i]['gisaidClade']+'\t'
                + infos[i]['nt']+'\t'
                + infos[i]['aa']+'\n'
            )


def main():
    warnings.filterwarnings('ignore')

    # command line interface
    parser = argparse.ArgumentParser(description='Get infos')
    parser.add_argument('--divergent-pathways', required=True, help='Divergent pathways file')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--snps', required=True, help='SNPs file')
    parser.add_argument('--surveillance', required=True, help='Variant surveillance file')
    parser.add_argument('--output', required=True, help='Infos file')
    args = parser.parse_args()

    seqs_filtered = filter_seq(args.divergent_pathways)    
    meta = pd.read_csv(args.metadata, delimiter='\t', index_col=2)
    
    infos = {}
    for i in seqs_filtered:
        infos[i] = {}
        infos[i]['region'] = meta.loc[i, 'region_exposure']
        infos[i]['country'] = meta.loc[i, 'country_exposure']
        infos[i]['division'] = meta.loc[i, 'division_exposure']
        infos[i]['date'] = meta.loc[i, 'date']
        infos[i]['pangoLineage'] = meta.loc[i, 'pango_lineage']
        infos[i]['nextstrainClade'] = str(meta.loc[i, 'Nextstrain_clade'])
        infos[i]['gisaidClade'] = str(meta.loc[i, 'GISAID_clade'])

    # get nucleotide substitutions from snps.txt
    snps = read_seq(args.snps)
    for i in seqs_filtered:
        snps_formatted = []
        for snp in snps[i]:
            snps_formatted.append(snp.replace('_', ''))
        infos[i]['nt'] = '('+','.join(snps_formatted)+')'

    # get amino acid substitutions from variant surveillance
    vs = pd.read_csv(args.surveillance, delimiter='\t', index_col=0)
    for i in seqs_filtered:
        infos[i]['aa'] = vs.loc[i, 'AA Substitutions']

    write_new_file(args.output, infos)


if __name__ == '__main__':
    main()
