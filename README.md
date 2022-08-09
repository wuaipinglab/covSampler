# covSampler

Phylogenetic analysis has been widely used to describe, display and infer the evolutionary patterns of viruses. The unprecedented accumulation of SARS-CoV-2 genomes has provided valuable materials for the real-time study of SARS-CoV-2 evolution. However, the large number of SARS-CoV-2 genome sequences also poses great challenges for phylogenetic analysis.

Here, we developed a subsampling method named **covSampler** based on the spatiotemporal distribution and genetic variation of SARS-CoV-2 genomes to subsampling the large-scale SARS-CoV-2 genome data sets.

## Overview

<p align="center">
<img src="img/workflow.svg" width="80%" height="80%" />
</p>

### Workflow (Figure A)

#### Determination of sites of spreading mutations of SARS-CoV-2 (Figure B)

- Genome sites of nonsynonymous mutations that increased in frequency per week for four consecutive weeks on at least one continent are defined as sites of spreading mutations.
- The haplotype sequence of each genome is constructed according to these sites of spreading mutations.

#### Construction of divergent pathways (Figure C)

- Divergent pathways in covSampler are network-like results of viral clustering, constructed by connecting each pair of viral sequences that are geographically close `(from the same administrative division)`, temporally close `(collected <= 14 days apart)`, and genetically similar `(Hamming distance <= 1 between their haplotype sequences)`.
- Network with only internal links (i.e. individual cluster) --> divergent pathway.
- Each divergent pathway reflects the local dynamic transmission and evolution of viruses over a period of time.

#### Representative and comprehensive subsampling based on the divergent pathways (Figure D)

- During subsampling, genomes in the original data set are hierarchically divided into subsampling units, which are groups containing genomes with(in) the same continent, divergent pathway, month and haplotype.
- Then, each subsampling unit is assigned an desired number of subsamples. Representative subsampling and comprehensive subsampling are performed using different strategies to assign the desired number of subsamples.

## Usage

### Download the covSampler workflow

```
git clone https://github.com/wuaipinglab/covSampler.git
```

### Create and activate the covSampler conda environment

```
cd covSampler
conda env create -f environment.yaml
conda activate covsampler
```

### Prepare data

To use covSampler to analyze your own data, you’ll need to prepare two files:

  1. A FASTA file with viral genomic sequences.

  2. A corresponding TSV file with metadata describing each sequence.

See the [data preparation guide](docs/tutorial/prep_data.md) for detailed instructions.

### Run covSampler

After data preparation, you can run the covSampler workflow to subsampling. The workflow consists of two main parts:

  1. Data processing for subsampling.

  2. Subsampling.

For the same data set, once you've processed the data, you can perform subsampling multiple times by adjusting  parameters (location, date range, subsampling characteristic et al.) to get the subsamples.

See this [tutorial](docs/tutorial/run_workflow.md) for detailed instructions.

## Web application

We also provide a [web application](https://www.covsampler.net/) of covSampler. It allows users to perform global SARS-CoV-2 subsampling. The sequences used in the web application are sourced from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/).

## Acknowledgements

We gratefully acknowledge [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/) and all the authors, originating and submitting laboratories of the SARS-CoV-2 sequences for sharing their work in open databases.
