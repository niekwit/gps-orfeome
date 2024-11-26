if config["bin_number"] == 1:
    rule mageck:
        input: 
            counts="results/count/counts-aggregated.tsv",
        output:
            rnw="results/mageck/temp/{comparison}/{comparison}_summary.Rnw",
            gs="results/mageck/temp/{comparison}/{comparison}.gene_summary.txt",
            ss="results/mageck/temp/{comparison}/{comparison}.sgrna_summary.txt",
            norm="results/mageck/temp/{comparison}/{comparison}.normalized.txt",
        params:
            control_genes=mageck_control(),
            dir_name=lambda wc, output: os.path.dirname(output["rnw"]),
            test_sample=lambda wc: wc.comparison.split("_vs_")[0],
            control_sample=lambda wc: wc.comparison.split("_vs_")[1],
            extra=config["mageck"]["extra_mageck_arguments"],
        threads: 2
        resources:
            runtime=20
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck/{comparison}.log"
        shell:
            "mageck test "
            "--normcounts-to-file "
            "--count-table {input.counts} "
            "--treatment-id {params.test_sample} "
            "--control-id {params.control_sample} "
            "--output-prefix {params.dir_name}/{wildcards.comparison} "
            "{params.control_genes} "
            "{params.extra} "
            "2> {log}"


    rule barcode_rank_plot:
        input:
            "results/temp/mageck/{comparison}/{comparison}.sgrna_summary.txt",
        output:
            report("results/mageck_plots/{comparison}/barcode_rank.pdf", caption="../report/barcoderank.rst", category="MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "Barcode rank plot"})
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/barcode_rank_{comparison}.log"
        script:
            "../scripts/plot_barcoderank.R"


    # Change all references of sgRNA to barcode in all relevant files
    rule rename_to_barcode:
        input:
            gs="results/mageck/temp/{comparison}/{comparison}.gene_summary.txt",
            ss="results/mageck/temp/{comparison}/{comparison}.sgrna_summary.txt",
            norm="results/mageck/temp/{comparison}/{comparison}.normalized.txt",
        output:
            gs="results/mageck/{comparison}/{comparison}.gene_summary.txt",
            ss="results/mageck/{comparison}/{comparison}.barcode_summary.txt",
            norm="results/mageck/{comparison}/{comparison}.normalized.txt",
        threads: 1
        resources:
            runtime=5
        log:
            "logs/rename_to_barcode/{comparison}.log"
        conda:
            "../envs/count.yaml"
        script:
            "../scripts/rename_to_barcode.py"
    

    rule rename_to_barcode_count_file:
        input:
            "results/count/counts-aggregated.tsv",
        output:
            "results/count/barcode-counts-aggregated.tsv",
        threads: 1
        resources:
            runtime=5
        log:
            "logs/rename_to_barcode/count_file.log"
        conda:
            "../envs/count.yaml"
        shell:
            "sed 's/sgRNA/barcode/' {input} > {output} 2> {log}"
else:
    rule psi:
        input:
            counts="results/count/counts-aggregated.tsv",
        output:
            expand("results/psi/{comparison}.csv", comparison=COMPARISONS),
        threads: 1
        resources:
            runtime=10
        conda:
            "../envs/stats.yaml"
        log:
            "logs/calculate_psi_{comparison}.log"
        script:
            "../scripts/calculate_psi.py"


    rule plot_barcode_profiles:
        input:
            csv="results/psi/{comparison}.csv",
        output:
            outdir=dir("results/psi_plots/{comparison}/"),
        threads: 1
        resources:
            runtime=30
        conda:
            "../envs/stats.yaml"
        log:
            "logs/plot_psi_{comparison}.log"
        script:
            "../scripts/plot_barcode_profiles.R"

