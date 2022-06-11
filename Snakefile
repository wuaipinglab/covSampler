rule all:
    input:
        [os.path.join(config["data_directory"], "args"), config["subsampling"]["subsamples_path"]] if config.get("webserver_args") else config["subsampling"]["subsamples_path"]
        
rule clean_genome:
    input:
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv"),
        msa = os.path.join(config["data_directory"], "rawdata/msa.fasta"),
        surveillance = os.path.join(config["data_directory"], "rawdata/variant_surveillance.tsv"),
        tree = os.path.join(config["data_directory"], "rawdata/timetree.nwk")
    output:
        snps = os.path.join(config["data_directory"], "snps.txt")
    params:
        reference = config["reference_sequence"]
    threads: 8
    shell:
        """
        python3 scripts/clean_genomes.py \
        --threads {threads} \
        --metadata  {input.metadata} \
        --msa {input.msa} \
        --surveillance  {input.surveillance} \
        --tree  {input.tree} \
        --reference {params.reference} \
        --output {output.snps}
        """

rule get_nonsynonymous:
    input:
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv"),
        sc2 = os.path.join(config["data_directory"], "rawdata/SARS_CoV_2.csv"),
        snps = os.path.join(config["data_directory"], "snps.txt")
    output:
        nonsynonymous = os.path.join(config["data_directory"], "nonsynonymous.txt")
    shell:
        """
        python3 scripts/get_nonsynonymous.py \
        --metadata  {input.metadata} \
        --genome {input.sc2} \
        --snps {input.snps} \
        --output {output.nonsynonymous}
        """

rule get_key_sites:
    input:
        metadata = os.path.join(config["data_directory"], "rawdata/metadata.tsv"),
        nonsynonymous = os.path.join(config["data_directory"], "nonsynonymous.txt"),
        snps = os.path.join(config["data_directory"], "snps.txt")
    output:
        key_sites = os.path.join(config["data_directory"], "key_sites.csv")
    shell:
        """
        python3 scripts/get_key_sites.py \
        --metadata  {input.metadata} \
        --nonsynonymous {input.nonsynonymous} \
        --snps {input.snps} \
        --output {output.key_sites}
        """

rule construct_haplotype_sequence:
    input:
        key_sites = os.path.join(config["data_directory"], "key_sites.csv"),
        snps = os.path.join(config["data_directory"], "snps.txt")
    output:
        haplotype_sequences = os.path.join(config["data_directory"], "haplotype_sequences.txt")
    shell:
        """
        python3 scripts/construct_haplotype_sequences.py \
        --sites {input.key_sites} \
        --snps {input.snps} \
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
        snps = os.path.join(config["data_directory"], "snps.txt"),
        surveillance = os.path.join(config["data_directory"], "rawdata/variant_surveillance.tsv")
    output:
        infos = os.path.join(config["data_directory"], "infos.tsv")
    shell:
        """
        python3 scripts/get_infos.py \
        --divergent-pathways {input.divergent_pathways} \
        --metadata  {input.metadata} \
        --snps {input.snps} \
        --surveillance  {input.surveillance} \
        --output {output.infos}
        """

rule get_web_args:
    input:
        infos = os.path.join(config["data_directory"], "infos.tsv"),
        sc2 = os.path.join(config["data_directory"], "rawdata/SARS_CoV_2.csv")
    output:
        args = directory(os.path.join(config["data_directory"], "args"))
    shell:
        """
        python3 scripts/get_web_args.py \
        --infos {input.infos} \
        --genome {input.sc2} \
        --output {output.args}
        """

rule subsampling:
    input:
        divergent_pathways = os.path.join(config["data_directory"], "divergent_pathways.csv"),
        haplotype_sequences = os.path.join(config["data_directory"], "haplotype_sequences.txt"),
        infos = os.path.join(config["data_directory"], "infos.tsv"),
        sc2 = os.path.join(config["data_directory"], "rawdata/SARS_CoV_2.csv")
    output:
        subsamples = config["subsampling"]["subsamples_path"]
    params:
        description = "--description " + "\""+str(config["description"])+"\"" if config.get("description") else "",
        location = "\""+config["subsampling"]["location"]+"\"",
        date_start = config["subsampling"]["date_start"],
        date_end = config["subsampling"]["date_end"],
        variants = "".join([" --variants "+"\""+v+"\"" for v in config["subsampling"]["variants"]]) if config["subsampling"].get("variants") else "",
        size = config["subsampling"]["size"],
        characteristic = config["subsampling"]["characteristic"],
        temporally_even = "--temporally-even" if config["subsampling"].get("temporally_even") else ""
    shell:
        """
        python3 scripts/ncov_sampling.py \
        --divergent-pathways {input.divergent_pathways} \
        --genome {input.sc2} \
        --haplotypes {input.haplotype_sequences} \
        --infos {input.infos} \
        {params.description} \
        --location {params.location} \
        --date-start {params.date_start} \
        --date-end {params.date_end} \
        {params.variants} \
        --size {params.size} \
        --characteristic {params.characteristic} \
        {params.temporally_even} \
        --output {output.subsamples}
        """
