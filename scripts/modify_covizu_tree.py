import argparse
import warnings
import requests
from ete3 import Tree


def main():
    warnings.filterwarnings('ignore')

    # command line interface
    parser = argparse.ArgumentParser(description='Add missing lineages to covizu tree')
    parser.add_argument('--covizu-tree', required=True, help='CoVizu time-scaled tree file')
    parser.add_argument('--output', required=True, help='Modified tree file')
    args = parser.parse_args()

    # add lineages to tree
    t = Tree(args.covizu_tree, format=1)
    leaves = [leaf.name for leaf in t.get_leaves()]
    url = 'https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt'
    resp = requests.get(url)
    content = resp.content.decode()
    for i in content.split('\n')[1:]:
        l = i.split('\t')[0]
        if (not l.startswith('*')) and (l != ''):
            if l not in leaves:
                t.add_child(name=l, dist=1)

    # write modified tree
    t.write(outfile=args.output)


if __name__ == '__main__':
    main()
