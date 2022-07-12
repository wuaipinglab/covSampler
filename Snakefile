import os

all_input = []
all_input.append(os.path.join(config["data_directory"], "args"))
if config.get("start_subsampling"):
    all_input.append(config["subsampling"]["output_path"])

rule all:
    input:
        all_input        

rule prepare_nextclade_dataset:
    output:
        nextclade_dataset = directory(os.path.join(config["data_directory"], "nextclade_sc2_dataset"))
    params:
        name = "sars-cov-2"
    shell:
        """
        nextclade dataset get --name {params.name} --output-dir {output.nextclade_dataset}
        """

rule run_nextclade:
    input:
        sequences = os.path.join(config["data_directory"], "rawdata/sequences.fasta"),
        nextclade_dataset = os.path.join(config["data_directory"], "nextclade_sc2_dataset")
    output:
        nextclade_tsv = os.path.join(config["data_directory"], "nextclade_output/nextclade.tsv")
    params:
        min_length = config["sequence_min_length"],
        nextclade_dir = os.path.join(config["data_directory"], "nextclade_output")
    threads: 8
    shell:
        """
        nextclade \
            --jobs {threads} \
            --input-fasta {input.sequences} \
            --input-dataset {input.nextclade_dataset} \
            --min-length {params.min_length} \
            --output-dir {params.nextclade_dir} \
            --output-tsv {output.nextclade_tsv}  
        """

rule filter_sequences:
    input:
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv"),
        nextclade_tsv = os.path.join(config["data_directory"], "nextclade_output/nextclade.tsv")
    output:
        strains = os.path.join(config["data_directory"], "strains.txt")
    params:
        genbank_accession = "--include-genbank-accession" if config.get("genbank_accession") else "",
        covizu_tree = "--covizu-tree "+os.path.join(config["data_directory"], "rawdata/timetree.nwk") if config.get("covizu_tree") else ""
    shell:
        """
        python3 scripts/filter_sequences.py \
            --metadata {input.metadata} \
            --nextclade-tsv {input.nextclade_tsv} \
            {params.genbank_accession} \
            {params.covizu_tree} \
            --output {output.strains}
        """

rule get_nonsynonymous:
    input:
        genome = "defaults/genome.csv",
        strains = os.path.join(config["data_directory"], "strains.txt"),
        nextclade_tsv = os.path.join(config["data_directory"], "nextclade_output/nextclade.tsv"),
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv")
    output:
        nonsynonymous = os.path.join(config["data_directory"], "nonsynonymous.txt")
    shell:
        """
        python3 scripts/get_nonsynonymous.py \
            --genome {input.genome} \
            --strains {input.strains} \
            --nextclade-tsv {input.nextclade_tsv} \
            --metadata  {input.metadata} \
            --output {output.nonsynonymous}
        """

rule get_key_sites:
    input:
        strains = os.path.join(config["data_directory"], "strains.txt"),
        nonsynonymous = os.path.join(config["data_directory"], "nonsynonymous.txt"),
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv")
    output:
        key_sites = os.path.join(config["data_directory"], "key_sites.txt")
    shell:
        """
        python3 scripts/get_key_sites.py \
            --strains {input.strains} \
            --nonsynonymous {input.nonsynonymous} \
            --metadata  {input.metadata} \
            --output {output.key_sites}
        """

rule construct_haplotype_sequences:
    input:
        strains = os.path.join(config["data_directory"], "strains.txt"),
        nextclade_tsv = os.path.join(config["data_directory"], "nextclade_output/nextclade.tsv"),
        key_sites = os.path.join(config["data_directory"], "key_sites.txt")
    output:
        haplotype_sequences = os.path.join(config["data_directory"], "haplotype_sequences.txt")
    shell:
        """
        python3 scripts/construct_haplotype_sequences.py \
            --strains {input.strains} \
            --nextclade-tsv {input.nextclade_tsv} \
            --sites {input.key_sites} \
            --output {output.haplotype_sequences}
        """
        
rule construct_divergent_pathways:
    input:
        haplotype_sequences = os.path.join(config["data_directory"], "haplotype_sequences.txt"),
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv")
    output:
        divergent_pathways = os.path.join(config["data_directory"], "divergent_pathways.csv")
    threads: 8
    shell:
        """
        python3 scripts/construct_divergent_pathways.py \
            --threads {threads} \
            --haplotypes {input.haplotype_sequences} \
            --metadata  {input.metadata} \
            --output {output.divergent_pathways}
        """

rule get_infos:
    input:
        divergent_pathways = os.path.join(config["data_directory"], "divergent_pathways.csv"),
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv"),
        nextclade_tsv = os.path.join(config["data_directory"], "nextclade_output/nextclade.tsv")
    output:
        infos = os.path.join(config["data_directory"], "infos.tsv")
    params:
        genbank_accession = "--include-genbank-accession" if config.get("genbank_accession") else ""
    shell:
        """
        python3 scripts/get_infos.py \
            --divergent-pathways {input.divergent_pathways} \
            --metadata  {input.metadata} \
            --nextclade-tsv {input.nextclade_tsv} \
            {params.genbank_accession} \
            --output {output.infos}
        """

rule get_arg_values:
    input:
        infos = os.path.join(config["data_directory"], "infos.tsv"),
        reference_dir = "defaults/reference"
    output:
        args = directory(os.path.join(config["data_directory"], "args"))
    shell:
        """
        python3 scripts/get_arg_values.py \
            --infos {input.infos} \
            --reference-dir {input.reference_dir} \
            --output {output.args}
        """

rule subsampling:
    input:
        infos = os.path.join(config["data_directory"], "infos.tsv"),
        haplotype_sequences = os.path.join(config["data_directory"], "haplotype_sequences.txt"),
        divergent_pathways = os.path.join(config["data_directory"], "divergent_pathways.csv")
    output:
        subsamples = config["subsampling"]["output_path"]
    params:
        description = "--description "+"\""+str(config["subsampling"]["description"])+"\"" if config["subsampling"].get("description") else "",
        location = "\""+config["subsampling"]["location"]+"\"",
        date_start = config["subsampling"]["date_start"],
        date_end = config["subsampling"]["date_end"],
        variants = "".join([" --variants "+"\""+v+"\"" for v in config["subsampling"]["variants"]]) if config["subsampling"].get("variants") else "",
        size = config["subsampling"]["size"],
        characteristic = config["subsampling"]["characteristic"],
        temporally_even = "--temporally-even" if config["subsampling"].get("temporally_even") else "",
        genbank_accession = "--return-genbank-accession" if config.get("genbank_accession") else ""
    shell:
        """
        python3 scripts/ncov_sampling.py \
            --infos {input.infos} \
            --haplotypes {input.haplotype_sequences} \
            --divergent-pathways {input.divergent_pathways} \
            {params.description} \
            --location {params.location} \
            --date-start {params.date_start} \
            --date-end {params.date_end} \
            {params.variants} \
            --size {params.size} \
            --characteristic {params.characteristic} \
            {params.temporally_even} \
            {params.genbank_accession} \
            --output {output.subsamples}
        """
