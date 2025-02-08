if config["bin_number"] == 1:
    if config["mageck"]["run"]:
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


        rule lfc_plots:
            input:
                "results/mageck/temp/{comparison}/{comparison}.gene_summary.txt",
            output:
                pos=report("results/mageck_plots/{comparison}/{comparison}.lfc_pos.pdf", caption="../report/lfc_pos.rst", category="MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "lfc plot enriched genes"}),
                neg=report("results/mageck_plots/{comparison}/{comparison}.lfc_neg.pdf", caption="../report/lfc_neg.rst", category="MAGeCK plots", subcategory="{comparison}", labels={"Comparison":"{comparison}","Figure": "lfc plot depleted genes"}),
            threads: 1
            resources:
                runtime=5
            conda:
                "../envs/stats.yaml"
            log:
                "logs/mageck_plots/lfc_{comparison}.log"
            script:
                "../scripts/lfc_plots.R"


        rule barcode_rank_plot:
            input:
                "results/mageck/temp/{comparison}/{comparison}.sgrna_summary.txt",
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
                "../envs/stats.yaml"
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
                "../envs/stats.yaml"
            shell:
                "sed 's/sgRNA/barcode/' {input} > {output} 2> {log}"

    if config["drugz"]["run"]:
        rule install_drugz:
            output:
                directory("resources/drugz"),
            log:
                "logs/drugz/install.log"
            threads: 1
            resources:
                runtime=5
            conda:
                "../envs/stats.yaml"
            shell:
                "git clone https://github.com/hart-lab/drugz.git {output} 2> {log}"

        
        rule drugz:
            input:
                counts="results/count/counts-aggregated.tsv",
                drugz="resources/drugz",
            output:
                report("results/drugz/{comparison}.txt", caption="../report/drugz.rst", category="DrugZ", subcategory="{comparison}", labels={"Comparison":"{comparison}", "Figure":"DrugZ output"})
            params:
                test=lambda wc, output: wc.comparison.split("_vs_")[0].replace("-", ","),
                control=lambda wc, output: wc.comparison.split("_vs_")[1].replace("-", ","),
                extra=config["drugz"]["extra"]
            threads: 2
            resources:
                runtime=15
            conda:
                "../envs/stats.yaml"
            log:
                "logs/drugz/{comparison}.log"
            shell:
                "python {input.drugz}/drugz.py "
                "-i {input.counts} "
                "-c {params.control} "
                "-x {params.test} "
                "{params.extra} "
                "-o {output} 2> {log} "

else:
    rule calculate_psi:
        input:
            counts="results/count/counts-aggregated.tsv",
        output:
            csv="results/psi/hit-th{ht}_sd-th{st}_prop_th{pt}_penalty{p}/{comparison}.csv",
            ranked="results/psi/hit-th{ht}_sd-th{st}_prop_th{pt}_penalty{p}/{comparison}_ranked.csv",
        threads: 1
        resources:
            runtime=10
        conda:
            "../envs/stats.yaml"
        log:
            "logs/calculate_psi/{comparison}/hit-th{ht}_sd-th{st}_prop_th{pt}_penalty{p}.log"
        script:
            "../scripts/calculate_psi.py"


    rule plot_barcode_profiles:
        input:
            csv="results/psi/hit-th{ht}_sd-th{st}_prop_th{pt}_penalty{p}/{comparison}.csv",
        output:
            flag=temp("results/psi_plots/hit-th{ht}_sd-th{st}_prop_th{pt}_penalty{p}/{comparison}/plotting_done.txt"),
        params:
            outdir=lambda wc,output: os.path.dirname(output["flag"]),
        threads: 18,
        resources:
            runtime=30,
        conda:
            "../envs/stats.yaml",
        log:
            "logs/plot_psi/hit-th{ht}_sd-th{st}_prop_th{pt}_penalty{p}_{comparison}.log",
        script:
            "../scripts/plot_barcode_profiles.R"

