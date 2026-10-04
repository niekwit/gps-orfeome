if config["bin_number"] == 1:
    if config["mageck"]["run"]:

        # category: Analysis
        rule mageck:
            input:
                counts="results/count/counts-aggregated.tsv",
            output:
                rnw="results/mageck/temp/{comparison}/{comparison}_summary.Rnw",
                gs="results/mageck/temp/{comparison}/{comparison}.gene_summary.txt",
                ss="results/mageck/temp/{comparison}/{comparison}.sgrna_summary.txt",
                norm="results/mageck/temp/{comparison}/{comparison}.normalized.txt",
            log:
                "logs/mageck/{comparison}.log",
            conda:
                "../envs/stats.yaml"
            threads: 2
            resources:
                runtime=20,
            params:
                control_genes=mageck_control(),
                dir_name=lambda wc, output: os.path.dirname(output["rnw"]),
                test_sample=lambda wc: wc.comparison.split("_vs_")[0],
                control_sample=lambda wc: wc.comparison.split("_vs_")[1],
                extra=config["mageck"]["extra_mageck_arguments"],
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

        # category: Analysis
        rule lfc_plots:
            input:
                "results/mageck/temp/{comparison}/{comparison}.gene_summary.txt",
            output:
                pos=report(
                    "results/mageck_plots/{comparison}/{comparison}.lfc_pos.pdf",
                    caption="../report/lfc_pos.rst",
                    category="MAGeCK plots",
                    subcategory="{comparison}",
                    labels={
                        "Comparison": "{comparison}",
                        "Figure": "lfc plot enriched genes",
                    },
                ),
                neg=report(
                    "results/mageck_plots/{comparison}/{comparison}.lfc_neg.pdf",
                    caption="../report/lfc_neg.rst",
                    category="MAGeCK plots",
                    subcategory="{comparison}",
                    labels={
                        "Comparison": "{comparison}",
                        "Figure": "lfc plot depleted genes",
                    },
                ),
            log:
                "logs/mageck_plots/lfc_{comparison}.log",
            conda:
                "../envs/stats.yaml"
            threads: 1
            resources:
                runtime=5,
            script:
                "../scripts/plot_lfc.R"

        # category: Analysis
        rule barcode_rank_plot:
            input:
                "results/mageck/temp/{comparison}/{comparison}.sgrna_summary.txt",
            output:
                report(
                    "results/mageck_plots/{comparison}/barcode_rank.pdf",
                    caption="../report/barcoderank.rst",
                    category="MAGeCK plots",
                    subcategory="{comparison}",
                    labels={
                        "Comparison": "{comparison}",
                        "Figure": "Barcode rank plot",
                    },
                ),
            log:
                "logs/mageck_plots/barcode_rank_{comparison}.log",
            conda:
                "../envs/stats.yaml"
            threads: 1
            resources:
                runtime=5,
            script:
                "../scripts/plot_barcoderank.R"

        # category: Analysis
        rule rename_to_barcode:
            input:
                gs="results/mageck/temp/{comparison}/{comparison}.gene_summary.txt",
                ss="results/mageck/temp/{comparison}/{comparison}.sgrna_summary.txt",
                norm="results/mageck/temp/{comparison}/{comparison}.normalized.txt",
            output:
                gs="results/mageck/{comparison}/{comparison}.gene_summary.txt",
                ss="results/mageck/{comparison}/{comparison}.barcode_summary.txt",
                norm="results/mageck/{comparison}/{comparison}.normalized.txt",
            log:
                "logs/rename_to_barcode/{comparison}.log",
            conda:
                "../envs/stats.yaml"
            threads: 1
            resources:
                runtime=5,
            script:
                "../scripts/rename_to_barcode.py"

        # category: Analysis
        rule rename_to_barcode_count_file:
            input:
                "results/count/counts-aggregated.tsv",
            output:
                "results/count/barcode-counts-aggregated.tsv",
            log:
                "logs/rename_to_barcode/count_file.log",
            conda:
                "../envs/stats.yaml"
            threads: 1
            resources:
                runtime=5,
            shell:
                "sed 's/sgRNA/barcode/' {input} > {output} 2> {log}"

    if config["drugz"]["run"]:

        # category: Analysis
        rule install_drugz:
            output:
                directory("resources/drugz"),
            log:
                "logs/drugz/install.log",
            conda:
                "../envs/stats.yaml"
            threads: 1
            resources:
                runtime=5,
            shell:
                "git clone https://github.com/hart-lab/drugz.git {output} 2> {log}"

        # category: Analysis
        rule drugz:
            input:
                counts="results/count/counts-aggregated.tsv",
                drugz="resources/drugz",
            output:
                report(
                    "results/drugz/{comparison}.txt",
                    caption="../report/drugz.rst",
                    category="DrugZ",
                    subcategory="{comparison}",
                    labels={"Comparison": "{comparison}", "Figure": "DrugZ output"},
                ),
            log:
                "logs/drugz/{comparison}.log",
            conda:
                "../envs/stats.yaml"
            threads: 2
            resources:
                runtime=15,
            params:
                test=lambda wc, output: wc.comparison.split("_vs_")[0].replace("-", ","),
                control=lambda wc, output: wc.comparison.split("_vs_")[1].replace(
                    "-", ","
                ),
                extra=config["drugz"]["extra"],
            shell:
                "python {input.drugz}/drugz.py "
                "-i {input.counts} "
                "-c {params.control} "
                "-x {params.test} "
                "{params.extra} "
                "-o {output} 2> {log} "

else:

    # category: Analysis
    rule calculate_psi:
        input:
            counts="results/count/counts-aggregated.tsv",
        output:
            csv="results/psi/hit-th{ht}_prop_th{pt}/{comparison}_barcode.summary.csv",
            ranked="results/psi/hit-th{ht}_prop_th{pt}/{comparison}_gene.summary.csv",
            sums=temp("results/psi/hit-th{ht}_prop_th{pt}/{comparison}_sums.csv"),
        log:
            "logs/calculate_psi/{comparison}/hit-th{ht}_prop_th{pt}.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        resources:
            runtime=10,
        script:
            "../scripts/calculate_psi.py"

    # category: Analysis
    rule plot_barcode_profiles:
        input:
            ranked="results/psi/hit-th{ht}_prop_th{pt}/{comparison}_gene.summary.csv",
            proportions="results/psi/hit-th{ht}_prop_th{pt}/{comparison}_barcode.summary.csv",
        output:
            d=report(
                directory(
                    "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}/destabilised/"
                ),
                patterns=["{name}.pdf"],
                caption="../report/profiles.rst",
                category="Barcode profiles {comparison}",
                subcategory="Destabilised",
            ),
            s=report(
                directory(
                    "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}/stabilised/"
                ),
                patterns=["{name}.pdf"],
                caption="../report/profiles.rst",
                category="Barcode profiles {comparison}",
                subcategory="Stabilised",
            ),
            flag=temp(
                touch(
                    "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}/plotting_done.txt"
                )
            ),
        log:
            "logs/plot_psi/hit-th{ht}_prop_th{pt}_{comparison}.log",
        conda:
            "../envs/stats.yaml"
        threads: 18
        resources:
            runtime=60,
        params:
            outdir=lambda wc, output: os.path.dirname(output["flag"]),
            bin_number=config["bin_number"],
        script:
            "../scripts/plot_barcode_profiles.R"

    # category: Analysis
    rule plot_dpsi_rank:
        input:
            ranked="results/psi/hit-th{ht}_prop_th{pt}/{comparison}_gene.summary.csv",
        output:
            pdf=report(
                "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}_dpsi_rank.pdf",
                caption="../report/dpsi_rank.rst",
                category="PSI rank plots",
                subcategory="{comparison}",
                labels={
                    "Comparison": "{comparison}",
                    "Figure": "Ranked dPSI values",
                },
            ),
        log:
            "logs/plot_psi/dpsi_rank_hit-th{ht}_prop_th{pt}_{comparison}.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        resources:
            runtime=5,
        script:
            "../scripts/plot_dpsi_rank.R"

    # category: Analysis
    rule plot_histograms:
        input:
            csv="results/psi/hit-th{ht}_prop_th{pt}/{comparison}_barcode.summary.csv",
            sums="results/psi/hit-th{ht}_prop_th{pt}/{comparison}_sums.csv",
        output:
            psi=report(
                "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}_psi_histogram.pdf",
                caption="../report/histograms.rst",
                category="Histograms",
                subcategory="{comparison}",
                labels={
                    "Comparison": "{comparison}",
                    "Figure": "Histogram of PSI values",
                },
            ),
            dpsi=report(
                "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}_dpsi_histogram.pdf",
                caption="../report/histograms.rst",
                category="Histograms",
                subcategory="{comparison}",
                labels={
                    "Comparison": "{comparison}",
                    "Figure": "Histogram of dPSI values",
                },
            ),
            dpsi_sd=report(
                "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}_dpsi_sd_histogram.pdf",
                caption="../report/histograms.rst",
                category="Histograms",
                subcategory="{comparison}",
                labels={
                    "Comparison": "{comparison}",
                    "Figure": "Histogram of dPSI SD values",
                },
            ),
            sob=report(
                "results/psi_plots/hit-th{ht}_prop_th{pt}/{comparison}_sob_histogram.pdf",
                caption="../report/histograms.rst",
                category="Histograms",
                subcategory="{comparison}",
                labels={
                    "Comparison": "{comparison}",
                    "Figure": "Histogram of SOB values",
                },
            ),
        log:
            "logs/plot_histograms/hit-th{ht}_prop_th{pt}_{comparison}.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        resources:
            runtime=10,
        params:
            sob_threshold=config["psi"]["sob_threshold"],
        script:
            "../scripts/plot_histograms.R"

    # category: Analysis
    rule merge_gene_summary_data:
        input:
            ranks=expand(
                "results/psi/hit-th{{ht}}_prop_th{{pt}}/{comparison}_gene.summary.csv",
                comparison=COMPARISONS,
            ),
        output:
            "results/psi/hit-th{ht}_prop_th{pt}/gene.summary_all_conditions.csv",
        log:
            "logs/merge_rank_data_all_conditions/hit-th{ht}_prop_th{pt}.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        resources:
            runtime=5,
        script:
            "../scripts/merge_gene_summary_data.py"

    # category: Analysis
    rule plot_heatmap:
        input:
            csv="results/psi/hit-th{ht}_prop_th{pt}/gene.summary_all_conditions.csv",
        output:
            pdf=report(
                "results/psi_plots/hit-th{ht}_prop_th{pt}/heatmap.pdf",
                caption="../report/heatmap.rst",
                category="Heatmap multi conditions",
                subcategory="{ht}_{pt}",
                labels={
                    "Hit threshold": "{ht}",
                    "Proportion threshold": "{pt}",
                    "Figure": "Heatmap of dPSI values",
                },
            ),
            csv="results/psi_plots/hit-th{ht}_prop_th{pt}/heatmap_data.csv",
        log:
            "logs/plot_heatmap/hit-th{ht}_prop_th{pt}.log",
        conda:
            "../envs/stats.yaml"
        threads: 1
        resources:
            runtime=10,
        params:
            bin_number=config["bin_number"],
            clusters=config["psi"]["heatmap"]["clusters"],
            rownames=config["psi"]["heatmap"]["rownames"],
            fontsize=config["psi"]["heatmap"]["row_font_size"],
            width=config["psi"]["heatmap"]["width"],
            height=config["psi"]["heatmap"]["height"],
        script:
            "../scripts/plot_heatmap.R"
