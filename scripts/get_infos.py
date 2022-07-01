'''
get region, country, division, date, pango lineage, nextstrain clade and mutations of each strain
'''

import argparse
import warnings
import pandas as pd


def main():
    warnings.filterwarnings('ignore')

    # command line interface
    parser = argparse.ArgumentParser(description='Get infos')
    parser.add_argument('--divergent-pathways', required=True, help='Divergent pathways file')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--nextclade-tsv', required=True, help='Nextclade tsv file')
    parser.add_argument('--include-genbank-accession', action='store_true', help='Include genbank accession')
    parser.add_argument('--output', required=True, help='Infos file')
    args = parser.parse_args()

    # read and filter strains from divergent pathway file
    seqs_filtered = []
    pathways = {}
    with open(args.divergent_pathways) as f:
        lines = f.readlines()
        for n, line in enumerate(lines):
            if n != 0:
                seq_id = line.split(',')[0]
                pathway = line.strip().split(',')[1]
                pathways.setdefault(pathway, []).append(seq_id)
    for p in pathways:
        if len(pathways[p]) > 2:
            seqs_filtered.extend(pathways[p])

    # get meta
    meta_cols = ['strain', 'date', 'region_exposure', 'country_exposure', 'division_exposure', 'pango_lineage']
    if args.include_genbank_accession:
        meta_cols.append('genbank_accession')
    meta = pd.read_csv(args.metadata, delimiter='\t', usecols=meta_cols, index_col='strain')
    meta = meta[meta.index.isin(seqs_filtered)]

    # get nextstrain clade and mutations
    nextclade_cols = ['seqName', 'clade', 'substitutions', 'aaSubstitutions', 'aaDeletions']
    nextclade = pd.read_csv(args.nextclade_tsv, delimiter='\t', usecols=nextclade_cols, index_col='seqName')
    nextclade = nextclade[nextclade.index.isin(seqs_filtered)]
    nextclade.fillna('NA', inplace=True)

    # get infos
    infos = {}
    for i in seqs_filtered:
        infos[i] = {}
        infos[i]['region'] = meta.loc[i, 'region_exposure']
        infos[i]['country'] = meta.loc[i, 'country_exposure']
        infos[i]['division'] = meta.loc[i, 'division_exposure']
        infos[i]['date'] = meta.loc[i, 'date']
        infos[i]['pangoLineage'] = meta.loc[i, 'pango_lineage']
        infos[i]['nextstrainClade'] = nextclade.loc[i, 'clade']
        infos[i]['substitutions'] = nextclade.loc[i, 'substitutions']
        infos[i]['aaSubstitutions'] = nextclade.loc[i, 'aaSubstitutions']
        infos[i]['aaDeletions'] = nextclade.loc[i, 'aaDeletions']
        if args.include_genbank_accession:
            infos[i]['genbankAccession'] = meta.loc[i, 'genbank_accession']

    with open(args.output, 'w') as f:
        if args.include_genbank_accession:
            f.write('ID'+'\t'+'region'+'\t'+'country'+'\t'+'division'+'\t'+'date'+'\t'+'pangoLineage'+'\t'
                +'nextstrainClade'+'\t'+'substitutions'+'\t'+'aaSubstitutions'+'\t'+'aaDeletions'+'\t'
                +'genbankAccession'+'\n')
            for i in infos:
                f.write(
                    i + '\t'
                    + infos[i]['region'] + '\t'
                    + infos[i]['country'] + '\t'
                    + infos[i]['division'] + '\t'
                    + infos[i]['date'] + '\t'
                    + infos[i]['pangoLineage'] + '\t'
                    + infos[i]['nextstrainClade'] + '\t'
                    + infos[i]['substitutions'] + '\t'
                    + infos[i]['aaSubstitutions'] + '\t'
                    + infos[i]['aaDeletions'] + '\t'
                    + infos[i]['genbankAccession'] + '\n'
                )
        else:
            f.write('ID'+'\t'+'region'+'\t'+'country'+'\t'+'division'+'\t'+'date'+'\t'+'pangoLineage'+'\t'
                +'nextstrainClade'+'\t'+'substitutions'+'\t'+'aaSubstitutions'+'\t'+'aaDeletions'+'\n')
            for i in infos:
                f.write(
                    i + '\t'
                    + infos[i]['region'] + '\t'
                    + infos[i]['country'] + '\t'
                    + infos[i]['division'] + '\t'
                    + infos[i]['date'] + '\t'
                    + infos[i]['pangoLineage'] + '\t'
                    + infos[i]['nextstrainClade'] + '\t'
                    + infos[i]['substitutions'] + '\t'
                    + infos[i]['aaSubstitutions'] + '\t'
                    + infos[i]['aaDeletions'] + '\n'
                )


if __name__ == '__main__':
    main()
