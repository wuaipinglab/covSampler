import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser(description='Compress metadata')
    parser.add_argument('--infos', required=True, help='Infos file')
    parser.add_argument('--output', required=True, help='Compressed metadata file')
    args = parser.parse_args()

    seqs = {}
    with open(args.infos) as f:
        lines = f.readlines()
        for n, line in enumerate(lines):
            if n != 0:
                seq_items = line.strip().split('\t')
                seq_id = seq_items[10]
                seq_date = seq_items[4]
                seq_region = seq_items[1]                
                seq_pangolineage = seq_items[5]
                seqs[seq_id] = {
                    'date': seq_date,
                    'region_exposure': seq_region,
                    'pango_lineage': seq_pangolineage
                }

    with open(args.output, 'w') as f:
        f.write('strain'+'\t'+'date'+'\t'+'region_exposure'+'\t'+'pango_lineage'+'\n')
        for i in seqs:
            f.write(i+'\t'+seqs[i]['date']+'\t'+seqs[i]['region_exposure']+'\t'+seqs[i]['pango_lineage']+'\n')

    subprocess.run('gzip '+args.output, shell=True)


if __name__ == '__main__':
    main()
