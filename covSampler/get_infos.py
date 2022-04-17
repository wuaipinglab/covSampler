'''
prepare region, country, division, date, pango lineage, nextstrain clade,
gisaid clade, nucleotide substitutions and amino acid substitutions of
each sequence for subsampling parameter selection
'''

import warnings
import pandas as pd
from covSampler import divergent_pathway_path, meta_path, snps_path, variant_surveillance_path, infos_path
from get_nonsynonymous import read_seq


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
    seqs_filtered = filter_seq(divergent_pathway_path)    
    meta = pd.read_csv(meta_path, delimiter='\t', index_col=2)
    
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
    snps = read_seq(snps_path)
    for i in seqs_filtered:
        snps_formatted = []
        for snp in snps[i]:
            snps_formatted.append(snp.replace('_', ''))
        infos[i]['nt'] = '('+','.join(snps_formatted)+')'

    # get amino acid substitutions from variant surveillance
    vs = pd.read_csv(variant_surveillance_path, delimiter='\t', index_col=0)
    for i in seqs_filtered:
        infos[i]['aa'] = vs.loc[i, 'AA Substitutions']

    write_new_file(infos_path, infos)


if __name__ == '__main__':
    main()
