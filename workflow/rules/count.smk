rule create_fasta:
    input:
        csv=csv,
    output:
        fasta=fasta,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/create_fasta.log",
    script:
        "../scripts/csv_to_fasta.py"


rule bowtie2_index:
    input:
        ref=fasta,
    output:
        multiext(
            "resources/bowtie2_index/barcodes",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        extra="",
    log:
        "logs/bowtie2/index.log",
    threads: 6
    resources:
        runtime=20,
    wrapper:
        f"{wrapper_version}/bio/bowtie2/build"


rule cutadapt:
    input:
        f"reads/{{sample}}.fastq{EXT}",
    output:
        fastq=temp("results/trimmed/{sample}.fastq.gz"),
        qc="results/trimmed/{sample}.qc.txt",
    params:
        extra=f"-q 20 {cut_adapt_arg()}",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 4  # Set desired number of threads here
    resources:
        runtime=25,
    wrapper:
        f"{wrapper_version}/bio/cutadapt/se"


rule count_barcodes:
    input:
        fq="results/trimmed/{sample}.fastq.gz",
        idx=multiext(
            "resources/bowtie2_index/barcodes",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        temp("results/count/{sample}.barcode.counts.txt"),
    params:
        mm=config["bowtie2"]["mismatch"],
        idx=lambda wc, input: input["idx"][0].replace(".1.bt2", ""),
        extra=config["bowtie2"]["extra"],
    threads: 8
    resources:
        runtime=45,
    log:
        "logs/count/{sample}.log",
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/count_barcodes.sh"


rule create_count_table:
    input:
        files=expand("results/count/{sample}.barcode.counts.txt", sample=SAMPLES),
        csv=csv,
    output:
        "results/count/counts-aggregated.tsv",
    threads: 1
    resources:
        runtime=10,
    conda:
        "../envs/stats.yaml"
    log:
        "logs/count/aggregate_counts.log",
    script:
        "../scripts/create_count_table.py"
