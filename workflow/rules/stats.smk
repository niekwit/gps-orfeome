if config["bin_number"] == 1:
    rule mageck:
        input: 
            counts="results/count/counts-aggregated.tsv",
        output:
            rnw="results/mageck/{comparison}/{comparison}_summary.Rnw",
            gs=report("results/mageck/{comparison}/{comparison}.gene_summary.txt", caption="../report/mageck.rst", category="MAGeCK"),
            ss="results/mageck/{comparison}/{comparison}.sgrna_summary.txt",
            norm="results/mageck/{comparison}/{comparison}.normalized.txt",
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
            "results/mageck/{comparison}/{comparison}.sgrna_summary.txt",
        output:
            report("results/mageck_plots/{comparison}/{comparison}.barcoderank.pdf", caption="../report/sgrank.rst", category="MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "sgrank plot"})
        params:
            fdr=config["mageck"]["fdr"],
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        log:
            "logs/mageck_plots/sgrank_{comparison}.log"
        script:
            "../scripts/sgrank_plot.R"


    # Change all references of sgRNA to barcode in all relevant files
    rule rename_to_barcode:
        input:
            #"results/mageck/{comparison}/{comparison}.sgrna_summary.txt",
        output:
            #"results/mageck/{comparison}/{comparison}.sgrna_summary.barcode.txt"
        threads: 1
        resources:
            runtime=5
        conda:
            "../envs/stats.yaml"
        script:
            "../scripts/rename_to_barcode.py"

else:
    rule calculate_psi:
        input:
            counts="results/count/counts-aggregated.tsv",
        output:
            "results/psi/psi_all.csv"
        params:
            csv=csv,
        threads: 1
        resources:
            runtime=10
        conda:
            "../envs/stats.yaml"
        log:
            "logs/calculate_psi.log"
        script:
            "../scripts/calculate_psi.R"


    rule find_hits:
        input:
            psi="results/psi/psi_all.csv",
        output:
            hits="results/psi/hits_{comparison}.csv",
        params:
