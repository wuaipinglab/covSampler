#### Data directory structure
```
├── divergent_pathways.csv
├── haplotype_sequences.txt
├── infos.tsv # region, country, administrative division, date, pango lineage, nextstrain clade, gisaid clade, nucleotide substitutions and amino acid substitutions of each sequence
├── key_sites.csv
├── nonsynonymous.txt # accession id of sequences in each group (group -> nonsynonymous mutation + continent)
│
├── rawdata
│   ├── metadata.tsv
│   ├── sequences.fasta
│   ├── SARS_CoV_2.csv
│   ├── timetree.nwk
│   └── variant_surveillance.tsv
│
└── snps.txt
```