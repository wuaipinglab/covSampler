# Data preparation guide

To use covSampler to analyze your own data, you’ll need to prepare two files:

  1. A FASTA file with viral genomic sequences.

  2. A corresponding TSV file with metadata describing each sequence.

## Format your sequence data

Prepare your nucleotide sequences in a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format file named `sequences.fasta`.

You can see a formatted example sequence file [here](../data/example_project/rawdata/sequences.fasta).

## Format your metadata

Prepare your metadata in a [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) format file named `metadata.tsv`.

A metadata file must include the following fields:

|  fields  |  Description  |  Format  |
|  -  |  -  |  -  |
|  strain  |  Sequence name |  The strain values in the metadata file must match them in the fasta file  |
|  date  |  Collection date  |  YYYY-MM-DD (Ambiguous value is unacceptable)  |
|  region_exposure  |  Continent |  Africa / Asia / Europe / North America / Oceania / South America  |
|  country_exposure  |  Country |  Country  |
|  division_exposure  |  Administrative division  |  Division  |
|  pango_lineage*  |  Viral lineage under the [Pango nomenclature](https://cov-lineages.org/index.html)  |  See the [lastest Pango lineage list](https://cov-lineages.org/lineage_list.html)  |

\* Currently covSampler workflow does not include Pango lineage assignment. You can perform the Pango lineage assignment using [pangolin](https://cov-lineages.org/resources/pangolin.html) or [nextclade](https://clades.nextstrain.org/).

You can see a formatted example metadata file [here](../data/example_project/rawdata/metadata.tsv).

## Create your project data directory

All data are in the `data/` directory. The raw data and intermediate data of each project will be stored in its corresponding directory.

For a new project (here named `tutorial_project`):

1. Create your project data folder in `data/`.
   
2. Create `rawdata/` folder in `data/tutorial_project`.
   
3. Move your sequence data and metadata into `data/turotial_project/rawdata/` folder.

  Now, the `data/` directory structure should look like this:

  ```
  data
  ├── README.md
  ├── example_project
  │   └── rawdata
  │       ├── metadata.tsv
  │       └── sequences.fasta
  └── tutorial_project
      └── rawdata
          ├── metadata.tsv
          └── sequences.fasta
  ```

## What's next?

[Run covSampler with your data](../my_profiles/README.md)
