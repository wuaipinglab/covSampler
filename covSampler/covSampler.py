import os
import math
import time
import subprocess

DATE = '2022-03-30' # data version
REF = 'EPI_ISL_402125' # GISAID accession id of SARS-CoV-2 reference genome 
CORES = 12

# --- data ---
DIRPATH = os.path.join('/datadirectory/', DATE)
rawdata_path = os.path.join(DIRPATH, 'rawdata/')
for root, dirs, files in os.walk(rawdata_path):
    for file in files:
        if file == 'timetree.nwk':
            tree_path = os.path.join(root, file)
        if file.startswith('metadata') and file.endswith('.tsv'):
            meta_path = os.path.join(root, file)
        if file == 'variant_surveillance.tsv':
            variant_surveillance_path = os.path.join(root, file)
        if file.startswith('msa') and file.endswith('.fasta'):
            sequence_path = os.path.join(root, file)
        if file == 'SARS_CoV_2.csv':
            genome_path = os.path.join(root, file)
snps_path = os.path.join(DIRPATH, 'snps.txt')
nonsynonymous_path = os.path.join(DIRPATH, 'nonsynonymous.txt')
key_site_path = os.path.join(DIRPATH, 'key_site.csv')
haplotype_sequence_path = os.path.join(DIRPATH, 'haplotype_sequence.txt')
divergent_pathway_path = os.path.join(DIRPATH, 'divergent_pathway.csv')
infos_path = os.path.join(DIRPATH, 'infos.tsv')
args_path = os.path.join(DIRPATH, 'args')
arg_location_path = os.path.join(args_path, 'locations.txt')
arg_date_path = os.path.join(args_path, 'dates.txt')
arg_pangoLineage_path = os.path.join(args_path, 'pango_lineages.txt')
arg_nextstrainClade_path = os.path.join(args_path, 'nextstrain_clades.txt')
arg_gisaidClade_path = os.path.join(args_path, 'gisaid_clades.txt')
arg_WHO_path = os.path.join(args_path, 'who_variants.txt')
arg_aa_path = os.path.join(args_path, 'amino_acid.txt')
arg_nt_path = os.path.join(args_path, 'nucleotide.txt')

# --- scripts ---
scripts_path = '/scriptsdirecrory/'
clean_genome_script = os.path.join(scripts_path, 'clean_genome.py')
get_nonsynonymous_script = os.path.join(scripts_path, 'get_nonsynonymous.py')
get_key_sites_script = os.path.join(scripts_path, 'get_key_sites.py')
construct_haplotype_sequence_script = os.path.join(scripts_path, 'construct_haplotype_sequence.py')
construct_divergent_pathway_script = os.path.join(scripts_path, 'construct_divergent_pathway.py')
get_infos_script = os.path.join(scripts_path, 'get_infos.py')
get_web_args_script = os.path.join(scripts_path, 'get_web_args.py')


def main():
    scripts = [
        clean_genome_script, 
        get_nonsynonymous_script, 
        get_key_sites_script, 
        construct_haplotype_sequence_script, 
        construct_divergent_pathway_script,
        get_infos_script,
        get_web_args_script
        ]
    for script in scripts:
        step = script.split('/')[-1].split('.')[0].replace('_', ' ').upper()
        print(' '+'-'*50)
        print('|'+' '*int((50-len(step))/2)+step+' '*math.ceil((50-len(step))/2)+'|')
        print(' '+'-'*50)
        print('Start Time: ' + time.asctime(time.localtime(time.time())))
        subprocess.run('python3 '+script, shell=True)
        print('Finish Time: ' + time.asctime(time.localtime(time.time())))


if __name__ == '__main__':
    main()
