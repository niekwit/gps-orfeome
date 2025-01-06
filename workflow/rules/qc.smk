rule fastqc:
    input:
        "results/trimmed/{sample}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip",
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 2
    resources:
        runtime=20,
        mem_mb = 2048,
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        report("results/qc/multiqc.html", caption="../report/multiqc.rst", category="MultiQC"),
    params:
        extra="",  # Optional: extra parameters for multiqc
    threads: 4
    resources:
        runtime=20,
        mem_mb = 2048,
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        f"{wrapper_version}/bio/multiqc"


rule plot_alignment_rate:
    input:
        expand("logs/count/{sample}.log", sample=SAMPLES)
    output:
        report("results/qc/alignment-rates.pdf", caption="../report/alignment-rates.rst", category="Alignment rates")
    log:
        "logs/plot-alignment-rate.log"
    threads: 1
    resources:
        runtime=5,
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_alignment_rate.R"


rule plot_coverage:
    input:
        tsv="results/count/counts-aggregated.tsv",
        fasta=fasta,
    output:
        report("results/qc/sequence-coverage.pdf", caption="../report/plot-coverage.rst", category="Sequence coverage")
    params:
        bin_number=config["bin_number"]
    threads: 1
    resources:
        runtime=5,
    log:
        "logs/plot-coverage.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_coverage.R"


rule plot_missed_barcodes:
    input:
        "results/count/counts-aggregated.tsv"
    output:
        report("results/qc/missed-barcodes.pdf", caption="../report/missed-barcodes.rst", category="Missed barcodes")
    params:
        bin_number=config["bin_number"]
    threads: 1
    resources:
        runtime=5,
    log:
        "logs/missed-rgrnas.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot_missed_barcodes.R"