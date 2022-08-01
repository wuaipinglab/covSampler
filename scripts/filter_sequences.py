import argparse
import warnings
import pandas as pd
from Bio import Phylo


def main():
    warnings.filterwarnings('ignore')

    # command line interface
    parser = argparse.ArgumentParser(description='Filter sequences')
    parser.add_argument('--metadata', required=True, help='Metadata file')
    parser.add_argument('--covizu-tree', help='CoVizu time-scaled tree file')
    parser.add_argument('--nextclade-tsv', required=True, help='Nextclade tsv file')
    parser.add_argument('--include-genbank-accession', action='store_true', help='Include genbank accession')
    parser.add_argument('--output', required=True, help='Remaining strains file')
    args = parser.parse_args()

    # --- filter meta ---
    meta = pd.read_csv(args.metadata, delimiter='\t')

    meta_cols = ['strain', 'date', 'region_exposure', 'country_exposure', 'division_exposure', 'pango_lineage']
    if args.include_genbank_accession:
        meta_cols.append('genbank_accession')
        meta = meta[meta['host'] == 'Homo sapiens']

    # filter meta: na
    meta = meta.dropna(axis=0, subset=meta_cols)

    # filter meta: "?"
    for col in meta_cols:
        meta = meta[meta[col] != '?']

    # filter meta: unclear date
    meta['date'] = meta['date'].astype(str)
    meta = meta[meta['date'].str.contains('20\d\d-\d\d-\d\d', regex=True)]

    # filter meta: pango lineage not in covizu time-scaled tree (optional)
    if args.covizu_tree is not None:
        tree = Phylo.read(args.covizu_tree, 'newick')
        lineages_in_tree = [str(nodes) for nodes in tree.get_terminals()]
        meta = meta[meta['pango_lineage'].isin(lineages_in_tree)]

    # --- filter nextclade ---
    nextclade = pd.read_csv(args.nextclade_tsv, delimiter='\t')
    nextclade = nextclade[
        (nextclade['errors'].isna()) & (nextclade['warnings'].isna()) & (nextclade['failedGenes'].isna())
        ]

    # --- get remaining strains ---
    remaining_strains = nextclade[nextclade['seqName'].isin(meta['strain'])]['seqName']

    with open(args.output, 'w') as f:
        f.write('Strain'+'\n')
        for i in remaining_strains:
            f.write(i+'\n')


if __name__ == '__main__':
    main()
