#### Data directory structure
```
├── args # args in the webserver
│   ├── amino_acid.txt
│   ├── dates.txt
│   ├── gisaid_clades.txt
│   ├── locations.txt
│   ├── nextstrain_clades.txt
│   ├── nucleotide.txt
│   ├── pango_lineages.txt
│   └── who_variants.txt
│
├── divergent_pathway.csv
├── haplotype_sequence.txt
├── infos.tsv # region, country, administrative division, date, pango lineage, nextstrain clade, gisaid clade, nucleotide substitutions and amino acid substitutions of each sequence
├── key_site.csv
├── nonsynonymous.txt # accession id of sequences in each group (group -> nonsynonymous mutation + continent)
│
├── rawdata
│   ├── metadata.tsv
│   ├── msa.fasta
│   ├── SARS_CoV_2.csv
│   ├── timetree.nwk
│   └── variant_surveillance.tsv
│
└── snps.txt
```