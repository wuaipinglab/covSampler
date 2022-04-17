'''
1. filter sequences with incomplete metadata
2. filter sequences with low quality genomes
3. perform mutation calling
'''

import re
import warnings
from multiprocessing import Pool
import pandas as pd
from Bio import SeqIO, Phylo
from covSampler import REF, CORES, tree_path, meta_path, variant_surveillance_path, sequence_path, snps_path


def read_seq(file, meta, info_seqs, REF, times):
    '''
    read sequences from fasta file
    filter sequences with incomplete metadata
    filter sequences with low quality genomes
    
    :param file: str, msa file path
    :param meta: df, filtered metadata
    :param vs: df, filtered variant surveillance
    :param REF: str, reference sequence accession id (EPI_ISL_402125)
    :param times: int, loop execution times

    :return: boolean, calculate all sequences?
             dict, key: accession id; value: sequence
    '''
    # calculate the index range of sequences in this loop
    nstart = 500000 * (times - 1)
    nend = 500000 * times

    # read sequences from msa file
    msa_seqs = []
    seqs = {}
    with open(file) as f:
        for n, seq_record in enumerate(SeqIO.parse(f, 'fasta')):
            if seq_record.description != 'reference_sequence':
                seq_id = seq_record.description.split('|')[1]
                if nstart <= n < nend or seq_id == REF:
                    msa_seqs.append(seq_id)
                    seqs[seq_id] = str(seq_record.seq)

    # check if calculate all sequences
    if nend > n:
        finish_all_seq = True
    else:
        finish_all_seq = False

    # remove sequences without complete information
    miss_info_seqs = list(set(msa_seqs) - set(info_seqs))
    for miss_info_seq in miss_info_seqs:
        seqs.pop(miss_info_seq)

    # remove sequences with unclear date (clear date format: YYYY-MM-DD)
    ambiguous_time_seqs = []
    for seq_id in seqs:
        if not re.search('20\d\d-\d\d-\d\d', meta[seq_id]):
            ambiguous_time_seqs.append(seq_id)
    for ambiguous_time_seq in ambiguous_time_seqs:
        seqs.pop(ambiguous_time_seq)

    return finish_all_seq, seqs


def get_genome_sites(sequence):
    '''
    get the index of SARS-CoV-2 whole genome sites from reference sequence
    
    :param sequence: str, reference sequence (not trimmed)
    
    :return: list, index of SARS-CoV-2 whole genome sites
    '''
    sites = []
    n = 0
    for i in sequence:
        if i != '-':
            sites.append(n)
        n += 1
    return sites


def trim_seq(sequence, sites):
    '''
    trim sequence

    :param sequence: str, sequence to be trimmed
    :param sites: list, index of SARS-CoV-2 whole genome sites

    :return: str, trimed sequence
    '''
    seq_trimed_list = []
    for i in sites:
        seq_trimed_list.append(sequence[i])
    seq_trimed = ''.join(seq_trimed_list)
    return seq_trimed


def is_qualified(sequence):
    if sequence.count('A') + sequence.count('T') + sequence.count('G') + sequence.count('C') >= 27000:
        return True
    else:
        return False


def clean_seq(seq_id, sequence, sites, ref_genome):
    '''
    1. trim sequence
    2. remove sequence with low quality genome ('A'+'T'+'G'+'C' < 27,000)
    3. perform mutation calling

    :param seq_id: str, sequence accession id
    :param sequence: str, sequence (not trimmed)
    :param sites: list, index of SARS-CoV-2 whole genome sites
    :param ref_genome: str, reference sequence (trimmed)

    :return: dict, key: accession id; value: snps
    '''
    try:
        seq_trimed = trim_seq(sequence, sites)
        if is_qualified(seq_trimed):
            snps = get_snps(seq_trimed, ref_genome)
            return {seq_id: snps}
        return False
    except:
        return False


def get_snps(sequence, ref_genome):
    '''
    perform mutation calling

    :param sequence: str, query sequence
    :param ref_genome: str, reference sequence

    :return: list, snps of query sequence
    '''
    snps = []
    for gsite in range(0, len(ref_genome)):
        if sequence[gsite] != ref_genome[gsite] and sequence[gsite] in ['A', 'T', 'G', 'C', '-']:
            snps.append(str(gsite+1)+'_'+sequence[gsite])
    return snps


def write_new_file(file, seq_snps):
    with open(file, 'w') as f:
        for i in seq_snps:
            f.write(i+':'+','.join(seq_snps[i])+'\n')


def main():
    warnings.filterwarnings('ignore')
    # filter sequences in metadata file based on the following criteria:
    # 1. without date, region exposure, country exposure, division exposure and pango lineage
    # 2. with non-human hosts
    # 3. with the pango lineage not in the latest covizu phylogenetic tree (https://filogeneti.ca/CoVizu/)
    # filter sequences in variant surveillance file based on the following criteria:
    # 1. without aa substitutions information
    tree = Phylo.read(tree_path, "newick")
    lineages_in_tree = [str(nodes) for nodes in tree.get_terminals()]
    meta = pd.read_csv(meta_path, delimiter='\t', index_col=2).dropna(axis=0, subset=['date', 'region_exposure', 'country_exposure', 'division_exposure', 'pango_lineage'])
    meta = meta[(meta['host'] == 'Human') & (meta['pango_lineage'].isin(lineages_in_tree))]
    variant_surveillance = pd.read_csv(variant_surveillance_path, delimiter='\t').dropna(axis=0, subset=['AA Substitutions'])
    
    # get accession id of sequences with complete information (in filtered metadata and in filtered variant surveillance)
    info_seqs = list(set(meta.index.tolist()) & set(variant_surveillance['Accession ID'].values.tolist()))

    # perform mutation calling
    # calculate 500,000 sequences each time (because of the limited computing resources)
    finish_all_seq = False
    times = 0
    seq_snps = {}
    while not finish_all_seq:
        times += 1
        finish_all_seq, seqs = read_seq(sequence_path, meta['date'].to_dict(), info_seqs, REF, times)
        sites = get_genome_sites(seqs[REF])
        ref_genome = trim_seq(seqs[REF], sites)

        result_list = []
        p = Pool(CORES)
        for i in seqs:
            result = p.apply_async(clean_seq, args=(i, seqs[i], sites, ref_genome))
            result_list.append(result)
        p.close()
        p.join()
        for result in result_list:
            if result.get():
                seq_snps.update(result.get())

    seq_snps = dict(sorted(seq_snps.items(), key=lambda x: int(x[0].split('_')[2])))

    write_new_file(snps_path, seq_snps)

if __name__ == '__main__':
    main()
