rule bowtie2_index:
    input:
        fasta = fasta,
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
        "logs/bowtie2/index.log"
    threads: 6
    resources:
        runtime=20
    wrapper:
        f"{wrapper_version}/bio/bowtie2/build"


rule cutadapt:
    input:
        "reads/{sample}.fastq.gz",
    output:
        fastq=temp("results/trimmed/{sample}.fastq.gz"),
        qc="results/trimmed/{sample}.qc.txt",
    params:
        extra=f"-q 20 {cut_adapt_arg()}",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 4  # set desired number of threads here
    resources:
        runtime=25
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
        "results/count/{sample}.barcode.counts.txt"
    params:
        mm=config["mismatch"],
        idx=lambda wc, input: input["idx"][0].replace(".1.bt2", ""),
    threads: 8
    resources:
        runtime=45
    log:
        "logs/count/{sample}.log"
    conda:
        "../envs/count.yaml"
    shell:
        "bowtie2 --no-hd "
        "-p {threads} "
        "-t -N {params.mm} "
        "-x {params.idx} "
        "-U {input.fq} "
        "2> {log} | "
        "sed '/XS:/d' | "
        "cut -f 3 | "
        "sort | "
        "uniq -c | "
        "sed '1d;s/^ *//' "
        "> {output}"


rule create_count_table:
    input:
        files=expand("results/count/{sample}.barcode.counts.txt", sample=SAMPLES)
    output:
        "results/count/counts-aggregated.tsv"
    params:
        fa=fasta,
        csv=csv,
    conda:
        "../envs/count.yaml"
    log:
        "logs/count/aggregate_counts.log"
    script:
        "../scripts/join.py"
