# How to run covSampler?

You will learn how to run the workflow with your data.

## Prerequisites

1. Create and activate covSampler conda environment.

2. Prepare your sequence data and metadata.

## Be familiar with workflow configuration profile

Before we start running the workflow, let's get familiar with the workflow configuration profile.

You can configure the workflow by specifying values in the workflow configuration profile.

All profiles are in `my_profiles/` directory. Each project has its corresponding profile.

For the example project, its corresponding profile is in `my_profiles/example_profile/`. You can use it to better understand the workflow configuration profile.

### What's the config.yaml?

The contents of `my_profiles/example_project/config.yaml`:

```
configfile: my_profiles/example_project/parameters.yaml

cores: 2

printshellcmds: True
```

#### configfile
* type: str
* description: Configuration file of this project workflow
* format: my_profiles/\<project_name\>/parameter.yaml

#### cores
* type: int
* description: Maximum number of cores you want to use in the workflow

#### printshellcmds
* type: bool
* description: Print the commands

### What's the parameters.yaml?

The contents of `my_profiles/example_project/parameters.yaml`:

```
data_directory:
    data/example_project

sequence_min_length:
    27000

covizu_tree:
    False

genbank_accession:
    False

start_subsampling:
    False

subsampling:
    output_path:
        results/example_project/subsamples.txt
    description:
        Example subsampling
    location:
        Global
    date_start:
        2020-01-01
    date_end:
        2022-01-01
    variants:
        - Lineage/WHO/Alpha
        - Lineage/Nextstrain_clade/20I (Alpha, V1)
        - Lineage/Pango_lineage/B.1.1.7
        - Site/Nucleotide/A23403G
        - Site/Amino_acid/S:N501Y
        - Site/Amino_acid/S:H69-
        - Site/Amino_acid/ORF7a:Q62*
    size:
        500
    characteristic:
        Representative
    temporally_even:
        False
```

#### data_directory
* type: str
* description: Project data directory
* format: data/\<project_name\>

#### sequence_min_length
* type: int
* description: Sequences with length < min length will be removed

#### covizu_tree
* type: bool
* description: Remove sequences whose corresponding pango lineages are not in the provided [CoVizu](https://filogeneti.ca/CoVizu/) time-scaled tree. This parameter is mainly designed for the workflow for covSampler web application. In general, keep it = `False`.

#### genbank_accession
* type: bool
* description: Return genbank accession. This parameter is designed for the workflow for covSampler web application. Please keep it = `False`.

#### start_subsampling
* type: bool
* description: Start subsampling
        
        This pipeline consists of two parts: 1) data processing and 2) subsampling.

        1) Before subsampling, you should process the data for subsampling.
            In this stage, please set "start_subsampling" parameter = "False".

        2) After data processing, you can perform subsampling.
            In this stage, please set "start_subsampling" parameter = "True".

#### subsampling - output_path
* type: str
* description: Output path of subsamples
* format: results/\<project_name\>/\<subsamples_file_name\>.txt

#### subsampling - description
* type: False or str
* description: Description recorded in the output file
* examples:
    * `False`
    * `Global subsampling 2022-01-01`
    * `My subsampling project`

#### subsampling - location
* type: str
* description: Location of subsamples
* format: \<Continent\>/(\<Country\>)/(\<Division\>)
* examples:
    * `Global`
    * `Europe`
    * `Europe/United Kingdom`
    * `Europe/United Kingdom/England`
* note: For each project, after performing the data processing part, all available values will be recorded in the`data/<project_name>/args/locations.txt`

#### subsampling - date_start / date_end
* type: str
* description: Date range of subsamples
* format: YYYY-MM-DD
* examples:
    * `2020-01-01`
    * `2022-01-01`
* note: For each project, after performing the data processing part, all available values will be recorded in the`data/<project_name>/args/dates.txt`

#### subsampling - variants
* type: False or str
* description: Variant of subsamples. There are three submodules in the `variants` module: `Nonspecific`, `Lineage`, and `Site`. `Nonspecific` means no restrictions on lineage or mutation of subsamples. In `Lineage` submodule, `WHO (VOC or VOI)`, `Pango lineage` and `Nextstrain clade` are available. In `Site` submodule, `Nucleotide` and `Amino Acid` are available.
* format:
    1. Nonspecific: False
    2. Lineage:
       1. WHO: - Lineage/WHO/\<VOC or VOI\>
       2. Pango lineage: - Lineage/Pango_lineage/\<Pango lineage\>
       3. Nextstrain clade: - Lineage/Nextstrain_clade/\<Nextstrain clade\>
    3. Site:
       1. Nucleotide: - Site/Nucleotide/\<ref\>\<site\>\<mut\>
       2. Amino acid: - Site/Amino_acid/\<gene\>:\<ref\>\<site\>\<mut\> ("-" = deletion, "*" = stop)
* examples:
    * False
    * \- Lineage/WHO/Alpha
    * \- Lineage/Pango_lineage/B.1.1.7
    * \- Lineage/Nextstrain_clade/20I (Alpha, V1)
    * \- Site/Nucleotide/A23403G
    * \- Site/Amino_acid/S:N501Y
    * \- Site/Amino_acid/S:H69-
    * \- Site/Amino_acid/ORF7a:Q62*
* note:
    1. For each project, after performing the data processing part, all available values will be recorded in the `data/<project_name>/args/who_variants.txt`, `data/<project_name>/args/pango_lineages.txt`, `data/<project_name>/args/nextstrain_clades.txt`, `data/<project_name>/args/nucleotide.txt` and `data/<project_name>/args/amino_acid.txt`
    2. You can use just one query or put multiple queries together:
            
            # --- one query example ---
            variants:
            - Lineage/WHO/Alpha

            # --- multiple queries example ---
            variants:
            - Lineage/WHO/Alpha
            - Lineage/Nextstrain_clade/20I (Alpha, V1)
            - Lineage/Pango_lineage/B.1.1.7
            - Site/Nucleotide/A23403G
            - Site/Amino_acid/S:N501Y
            - Site/Amino_acid/S:H69-
            - Site/Amino_acid/ORF7a:Q62*
            
    3. [The SARS-CoV-2 genome map](../defaults/genome_map/README.md) may be helpful when specifying the amino acid mutation.

#### subsampling - size
* type: int
* description: Number of subsamples

#### subsampling - characteristic
* type: str
* description: Two subsampling characteristic, `Comprehensive` and `Representative`, are available. You can choose one of them according to the application scenario.
    
      - Comprehensive:
      1) Geographic distribution: Evenly distributed across continents.
      2) Temporal distribution: Similar with the original data set.
      3) Genetic variation: Higher genetic diversity (compared to representative subsampling).
      4) Application scenario: 
          a. When the geographic bias (at the continent level) is unacceptable.
          b. When exploring the phylogeny of specific viruses in a rich and diverse genetic background.
          c. When exploring the phylogeny of specific viruses with low similarity to widely prevalent variants in the original data set.

      - Representative:
      1) Geographic distribution: Consistent with the original data set.
      2) Temporal distribution: Similar with the original data set.
      3) Genetic variation: Similar with the original data set.
      4) Application scenario:
          When the subsamples are expected to approximately reflect the spatiotemporal distribution and key phylogenetic relationships of genomes in the original data set.

      In addition, you can also perform subsampling multiple times with different parameters (characteristic and/or range) to get the final subsamples. 

* format: Representative or Comprehensive

#### subsampling - temporally_even
* type: bool
* description: The number of subsamples is the same for each month
* note: 
  1. The subsamples within each month still fit the selected Comprehensive or Representative distribution characteristic.
  2. The temporally even option will take a few more minutes.

## Configure your project workflow for data processing

Before subsampling, you should process the data for subsampling.

1. Create your project (here named `tutorial_project`) profile folder in `my_profiles/` directory.
   
2. Copy `config.yaml` and `parameters.yaml` files from `my_profiles/example_project` to `my_profiles/tutorial_project`.

    Now, the `my_profiles/` directory structure should look like this:

    ```
    my_profiles
    ├── README.md
    ├── example_project
    │   ├── config.yaml
    │   └── parameters.yaml
    └── tutorial_project
        ├── config.yaml
        └── parameters.yaml
    ```

3. Change parameters in `my_profiles/tutorial_project/config.yaml` and `my_profiles/tutorial_project/parameters.yaml`.

       my_profiles/tutorial_project/config.yaml:
       
           configfile: change path to your project parameter.yaml path (my_profiles/tutorial_project/parameter.yaml)

       my_profiles/tutorial_project/parameter.yaml:
     
           data_directory: change path to your project data path (data/tutorial_project)
       
           start_subsampling: False

           subsampling: you can change these parameters after data processing
        
       Other parameters not mentioned can be adjusted as required.
        
## Data processing

Change directory to the `covSampler` directory if you are not there.

```
cd covSampler
```

Run the command for data processing.

```
snakemake --profile my_profiles/tutorial_project
```

The workflow can take a while to run. When the sequence numbers are in the millions, the workflow may run for several days.

## Configure your project workflow for subsampling

After data processing, you can change the parameters for subsampling.

    my_profiles/tutorial_project/parameter.yaml:

        start_subsampling: True

        subsampling: you can change these parameters as required

## Subsampling

Change directory to the `covSampler` directory if you are not there.

```
cd covSampler
```

Run the command for subsampling. (This is the second time running the command, the first run is for data processing.)

```
snakemake --profile my_profiles/tutorial_project
```

After calculation, you can see the result subsamples in your specified output file path.

## Subsampling the same data set multiple times

If you want to subsample the same data set multiple times, you don't need to re-process the data.

You can just change the subsampling parameters in `my_profiles/tutorial_project/parameter.yaml`. Don't forget to change the output file path.

Then, change directory to the `covSampler` directory and run the command again to get your new subsamples.

```
snakemake --profile my_profiles/tutorial_project
```
