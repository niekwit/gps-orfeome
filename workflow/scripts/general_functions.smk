import datetime
import os, glob, sys, re
import socket, platform
from importlib.metadata import version, PackageNotFoundError

import pandas as pd
from gpsw import __container_image_version__

from snakemake.logging import logger
from snakemake.shell import shell


def get_package_version():
    try:
        return version("gpsw")
    except PackageNotFoundError:
        # Package is not installed (e.g., running locally without installation)
        # For development purposes
        return "0.0.0-dev"


# Workflow version
VERSION = get_package_version()
DOCKER_VERSION = __container_image_version__
wrapper_version = "v5.2.1"

## This header can help with debugging, extend later with more info
# https://github.com/moiexpositoalonsolab/grenepipe/blob/53912b0749f90c7c29ecde2855657318e1269020/workflow/rules/initialize.smk
indent = 21

# Get some info on the platform and OS
pltfrm = platform.platform() + "\n" + (" " * indent) + platform.version()
try:
    # Not available in all versions, so we need to catch this
    ld = platform.linux_distribution()
    if len(ld):
        pltfrm += "\n" + (" " * indent) + ld
    del ld
except:
    pass
try:
    # Mac OS version comes back as a nested tuple?!
    # Need to merge the tuples...
    def merge_tuple(x, bases=(tuple, list)):
        for e in x:
            if type(e) in bases:
                for e in merge_tuple(e, bases):
                    yield e
            else:
                yield e

    mv = " ".join(merge_tuple(platform.mac_ver()))
    if not mv.isspace():
        pltfrm += "\n" + (" " * indent) + mv
    del mv, merge_tuple
except:
    pass

# Get nicely wrapped command line
cmdline = sys.argv[0]
for i in range(1, len(sys.argv)):
    if sys.argv[i].startswith("-"):
        cmdline += "\n" + (" " * indent) + sys.argv[i]
    else:
        cmdline += " " + sys.argv[i]

logger.info(
    "============================================================================="
)
logger.info("")
# https://texteditor.com/multiline-text-art/
logger.info(r"                     ██████╗ ██████╗ ███████╗██╗    ██╗")
logger.info(r"                    ██╔════╝ ██╔══██╗██╔════╝██║    ██║")
logger.info(r"                    ██║  ███╗██████╔╝███████╗██║ █╗ ██║")
logger.info(r"                    ██║   ██║██╔═══╝ ╚════██║██║███╗██║")
logger.info(r"                    ╚██████╔╝██║     ███████║╚███╔███╔╝")
logger.info(r"                     ╚═════╝ ╚═╝     ╚══════╝ ╚══╝╚══╝ ")
logger.info("")
logger.info(
    "  Date:              " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
)
logger.info("  Platform:          " + pltfrm)
logger.info("  Python:            " + str(sys.version.split(" ")[0]))
logger.info("  Snakemake:         " + snakemake.__version__)
logger.info("  GPSW:              " + VERSION)
logger.info("  Wrapper:           " + wrapper_version)
# if apptainer argument is used, print the version of the docker image

logger.info("  Docker image:      " + DOCKER_VERSION)
logger.info("  Command:           " + cmdline)
logger.info("")

# Print contents of profile if parsed from command line
# i.e. if --profile flag used
if "--profile" in sys.argv:
    profile = sys.argv[sys.argv.index("--profile") + 1]
    file_ = os.path.join(profile, "config.yaml")
    logger.info("  Using profile:     " + file_)
    counter = 0
    with open(file_, "r") as file:
        for line in file:
            if counter == 0:
                logger.info("  Profile contents:  " + line.strip())
                counter += 1
            else:
                logger.info(" " * indent + line.strip())
                counter += 1

logger.info(
    "============================================================================="
)
logger.info("")


def targets():
    TARGETS = [
        "results/qc/multiqc.html",
        "results/qc/alignment-rates.pdf",
        "results/qc/sequence-coverage.pdf",
        "results/qc/missed-barcodes.pdf",
    ]
    if config["bin_number"] == 1:
        if config["mageck"]["run"]:
            TARGETS.extend(
                [
                    expand(
                        "results/mageck/temp/{comparison}/{comparison}.gene_summary.txt",
                        comparison=COMPARISONS,
                    ),
                    "results/count/barcode-counts-aggregated.tsv",
                    expand(
                        "results/mageck/{comparison}/{comparison}.gene_summary.txt",
                        comparison=COMPARISONS,
                    ),
                    expand(
                        "results/mageck/{comparison}/{comparison}.barcode_summary.txt",
                        comparison=COMPARISONS,
                    ),
                    expand(
                        "results/mageck/{comparison}/{comparison}.normalized.txt",
                        comparison=COMPARISONS,
                    ),
                    expand(
                        "results/mageck_plots/{comparison}/barcode_rank.pdf",
                        comparison=COMPARISONS,
                    ),
                    expand(
                        "results/mageck_plots/{comparison}/{comparison}.lfc_{val}.pdf",
                        comparison=COMPARISONS,
                        val=["pos", "neg"],
                    ),
                ]
            )
        if config["drugz"]["run"]:
            TARGETS.extend(
                [
                    expand("results/drugz/{comparison}.txt", comparison=COMPARISONS),
                ]
            )
    else:
        TARGETS.extend(
            [
                expand(
                    "results/psi_plots/hit-th{ht}_sd-th{st}_prop_th{pt}_pen_th{pnth}/{comparison}/plotting_done.txt",
                    zip,
                    comparison=COMPARISONS,
                    ht=HIT_TH,
                    st=SD_TH,
                    pt=PROP_TH,
                    pnth=PEN_TH,
                ),
                expand(
                    "results/psi_plots/hit-th{ht}_sd-th{st}_prop_th{pt}_pen_th{pnth}/{comparison}_dotplot.pdf",
                    zip,
                    comparison=COMPARISONS,
                    ht=HIT_TH,
                    st=SD_TH,
                    pt=PROP_TH,
                    pnth=PEN_TH,
                ),
            ]
        )
        if multiple_conditions(COMPARISONS):
            TARGETS.extend(
                [
                    "results/qc/pca_plot.pdf",
                    expand(
                        "results/psi_plots_multi_conditions/hit-th{ht}_sd-th{st}_prop_th{pt}_pen_th{pnth}/plotting_done.txt",
                        zip,
                        ht=HIT_TH,
                        st=SD_TH,
                        pt=PROP_TH,
                        pnth=PEN_TH,
                    ),
                    expand(
                        "results/psi/hit-th{ht}_sd-th{st}_prop_th{pt}_pen_th{pnth}/gene.summary_all_conditions.csv",
                        zip,
                        ht=HIT_TH,
                        st=SD_TH,
                        pt=PROP_TH,
                        pnth=PEN_TH,
                    ),
                ]
            )

    return TARGETS


def csv():
    """
    Locate CSV file
    """
    csv = glob.glob("resources/*.csv")

    # Check if csv file is present
    assert len(csv) == 1, "Only one CSV file should be present in resources directory"
    return csv[0], csv[0].replace(".csv", ".fasta")


def cut_adapt_arg():
    """
    Generates Cutadapt argument for removing vector sequences
    and extra arguments specified in config.yml
    """
    five_prime = f"-g {config['cutadapt']['five_prime_adapter']}"
    three_prime = f"-a {config['cutadapt']['three_prime_adapter']}"
    length = config["cutadapt"]["barcode_length"]
    extra = config["cutadapt"]["extra"]

    if length == 0:
        return f"{five_prime} {three_prime} {extra}"
    else:
        length = f"--length {length}"
        return f"{five_prime} {length} {extra}"


def sample_names():
    """
    Get sample names from fastq files and check for invalid characters
    """
    bin_number = config["bin_number"]

    test_samples = config["conditions"]["test"]
    control_samples = config["conditions"]["control"]
    if len(test_samples) == 0 or len(control_samples) == 0:
        raise ValueError(
            "No test or control samples found in config file. Please check your config file."
        )
    assert len(test_samples) == len(
        control_samples
    ), "Number of test and control samples should be equal"

    # Check for invalid characters in sample names (non-alphanumeric characters or _)
    sample_names = list(set(test_samples + control_samples))
    invalid_sample_names = []
    for name in sample_names:
        if not re.search("^[a-zA-Z0-9_]*$", name):
            invalid_sample_names.append(name)

    if len(invalid_sample_names) > 0:
        invalid_sample_names = "\n".join(invalid_sample_names)
        print("ERROR: Only alphanumeric characters and _ are allowed in sample names")
        print(f"Invalid character(s) found in sample name(s):\n{invalid_sample_names}")
        sys.exit(1)

    # Check if fastq files are named properly and exist
    if bin_number == 1:
        fastq = [f"reads/{x}.fastq{EXT}" for x in sample_names]
        fastq_not_found = [x for x in fastq if not os.path.exists(x)]
        if len(fastq_not_found) > 0:
            fastq_not_found = "\n".join(fastq_not_found)
            (f"ERROR: Fastq file(s) not found:\n{fastq_not_found}")
            sys.exit(1)
        else:
            return sample_names
    else:
        fastq = [
            f"reads/{x}_{y}.fastq{EXT}"
            for x in sample_names
            for y in range(1, bin_number + 1)
        ]
        fastq_not_found = [x for x in fastq if not os.path.exists(x)]
        if len(fastq_not_found) > 0:
            fastq_not_found = "\n".join(fastq_not_found)
            print(f"ERROR: Fastq file(s) not found:\n{fastq_not_found}")
            sys.exit(1)

        else:
            sample_names = [
                f"{x}_{y}" for x in sample_names for y in range(1, bin_number + 1)
            ]
            return sample_names


def wildcard_values():
    """
    Load comparisons for MAGeCK/PSI
    """
    test_samples = config["conditions"]["test"]
    control_samples = config["conditions"]["control"]
    assert len(test_samples) == len(
        control_samples
    ), "Number of test and control samples should be equal"

    COMPARISONS = []
    for t, c in zip(test_samples, control_samples):
        COMPARISONS.append(f"{t}_vs_{c}")

    hit_th = config["psi"]["hit_threshold"]
    sd_th = config["psi"]["sd_threshold"]
    prop_th = config["psi"]["proportion_threshold"]
    pen_th = config["psi"]["penalty_factor"]
    threshold_list = [hit_th, sd_th, prop_th, pen_th]

    # All threshold lists should be the same length
    assert (
        len(hit_th) == len(sd_th) == len(prop_th) == len(pen_th)
    ), "Threshold lists are not the same length"

    # As not all permutations are used (zip argument is used with Snakemake expand),
    # the COMPARISONS list is zipped with the threshold lists and should be the same length
    # as the threshold lists.
    # Repeat all th values for each comparison
    extended_comparisons = []
    for i in COMPARISONS:
        extended_comparisons.extend([i] * len(hit_th))

    extended_thresholds = []
    for i in threshold_list:
        extended_thresholds.append(i * len(COMPARISONS))

    return (
        extended_comparisons,
        extended_thresholds[0],
        extended_thresholds[1],
        extended_thresholds[2],
        extended_thresholds[3],
    )


def mageck_control():
    """
    Load control genes for MAGeCK
    """
    if (
        config["mageck"]["mageck_control_barcodes"] == "all"
    ):  # Use all genes as controls
        control = ""
    else:  # Use genes from file set in config
        file = config["mageck"]["mageck_control_barcodes"]

        # Check if file exists
        assert os.path.exists(file), f"Control gene file ({file}) does not exist"
        control = f"--control-gene {file}"

    return control


def multiple_conditions(comparisons):
    """
    Check if multiple comparisons are ent
    """
    # Count unique comparisons
    unique_comparisons = set(comparisons)
    if len(unique_comparisons) > 1:
        return True
    else:
        return False


def get_extension():
    """
    Get extension for fastq files
    Can be .gz or bz2
    """
    files = glob.glob("reads/*fastq*")

    # Check if fastq files are present
    assert len(files) > 0, "No fastq files found in reads directory"

    # Check if all fastq files have the same extension
    extensions = set([os.path.splitext(x)[1] for x in files])
    assert len(extensions) == 1, "All fastq files should have the same extension"

    return extensions.pop()
