import argparse
import os, sys, re
import tarfile
import tempfile
import shutil
import subprocess
from urllib.request import urlretrieve

from . import utils


def fetch_code(args):
    """_summary_"""
    if args.tag == "latest":

        tag = utils.get_latest_release_tag()
        if tag is None:
            print("Error fetching the latest release tag.")
            return
    else:
        tag = args.tag

    print(f"Fetching code of release tag {tag}")

    # Check if the specified directory exists, if not create it
    if not os.path.exists(args.directory):
        print(f"Directory {args.directory} does not exist. Creating it.")
        os.makedirs(args.directory)

    # Download the tar.gz file in temp dir and untar to specified directory
    with tempfile.TemporaryDirectory(dir=args.directory) as tmpdirname:
        # Validate tag format to prevent URL manipulation
        pattern = r"^[vV](0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)$"
        if not re.match(pattern, tag):
            print(
                f"Error: Invalid tag format '{tag}'. Tags should only contain alphanumeric characters, dots, hyphens, and underscores."
            )
            sys.exit(1)

        url = "https://github.com/niekwit/gps-orfeome"
        target_url = f"{url}/archive/refs/tags/{tag}.tar.gz"
        download_file = os.path.join(tmpdirname, "target_repo.tar.gz")
        print(f"Downloading {target_url} to {download_file}")
        urlretrieve(target_url, download_file)

        # Unpack the whole tar.gz file
        print("Unpacking the downloaded tar.gz file")
        with tarfile.open(download_file) as tar:
            tar.extractall(path=tmpdirname)

        # Specify the directories to copy
        if args.test_data:
            dirs = [".test/config", ".test/reads", ".test/resources", "workflow"]
            end_message = "Workflow is ready to run with test data"
        else:
            dirs = ["config", "workflow"]

            # Create empty resources directory if it doesn't exist
            resources_dir = os.path.join(args.directory, "resources")
            os.makedirs(resources_dir, exist_ok=True)

            end_message = "Please configure the workflow in config/config.yaml and provide a CSV file with barcode information in resources/"

        # Get the unpacked directory name
        tag_no_v = re.sub("^v", "", tag)  # Remove the leading 'v' from the tag
        dirs = [os.path.join(tmpdirname, f"gps-orfeome-{tag_no_v}", d) for d in dirs]

        # Copy tar directories of unpacked archive to the specified directory
        for d in dirs:
            shutil.copytree(
                d, os.path.join(args.directory, os.path.basename(d)), dirs_exist_ok=True
            )

    print("Done!")
    print(end_message)


def run_workflow(args):
    """

    Args:
        args (_type_): _description_
    """
    # --- Dry-run workflow ---
    if args.dry_run:
        utils.dry_run(args)

    # --- Create rule graph ---
    # this functions as a dry-run to catch any errors, so no need
    # to run snakemake -np again (still make option to run a dry-run only)
    utils.create_rule_graph()

    # --- Run workflow ---
    # Create Snakemake command
    command = ["snakemake"]
    command.extend(utils.profile_arg(args))
    if args.snakemake_args is not None:
        # Remove leading and trailing quotes
        snakemake_args = args.snakemake_args.strip('"').strip("'")
        command.extend(snakemake_args.split())

    if not args.quiet:
        command.append("-p")
    else:
        command.extend(["--quiet", "all"])

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Snakemake: {e}")
        sys.exit(1)

    # --- Create report ---
    try:
        print("Creating report...")
        subprocess.run(
            ["snakemake", "--report", "report.html", "--quiet", "all"], check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error creating report: {e}")
        sys.exit(1)
    print("Report (report.html) created successfully!")


def main():
    # get version for the --version flag
    package_version = utils.get_package_version()

    # set up the argument parser
    parser = argparse.ArgumentParser(
        description="GPSW: A tool for analysing and processing Global Protein Stability Profiling data.",
        prog="gpsw",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {package_version}",
    )
    # define subparsers *before* parsing args
    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")

    # --- Fetch ---
    parser_fetch = subparsers.add_parser(
        "fetch",
        description="Fetch GPSW code from a specific release from https://github.com/niekwit/gps-orfeome.",
    )
    parser_fetch.add_argument(
        "-t",
        "--tag",
        type=str,
        default="latest",
        help="Release tag of the code to download",
    )
    parser_fetch.add_argument(
        "--test-data",
        action="store_true",
        help="Include test data in the download (and relating config and resources directories)",
    )
    parser_fetch.add_argument(
        "-d",
        "--directory",
        type=str,
        default=".",
        help="Directory to download the code to target diectory (default: current directory)",
    )
    parser_fetch.set_defaults(func=fetch_code)

    # --- Run ---
    parser_run = subparsers.add_parser(
        "run", description="Run the GPSW pipeline and create report."
    )

    parser_run.add_argument(
        "-p",
        "--profile",
        type=str,
        default=None,
        help="Path to Snakemake profile YAML file (only has to be provided at first run) (OPTIONAL, use value None for no profile)",
    )
    parser_run.add_argument(
        "--snakemake-args",
        type=str,
        help="""Extra Snakemake arguments (should be '"double-quoted"')""",
    )
    parser_run.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        default=False,
        help="Dry-run workflow only",
    )
    parser_run.add_argument(
        "-q", "--quiet", action="store_true", default=False, help="Run GPSW quietly"
    )

    parser_run.set_defaults(func=run_workflow)

    # parse arguments *after* defining all subparsers
    args = parser.parse_args()

    # --- Execute the corresponding function ---
    # This requires adding set_defaults(func=...) to each subparser
    # and defining those functions (fetch_code, dry_run_pipeline, run_pipeline)
    if hasattr(args, "func"):
        args.func(args)
    elif args.command is None and not any(
        arg for arg in sys.argv[1:] if arg != "--version"
    ):
        # If no subcommand was given (and it wasn't just --version) print help
        parser.print_help()
    # add logic here based on args.command if not using set_defaults
    # elif args.command == 'fetch':
    #    fetch_code(args)
    # elif args.command == 'dry-run':
    #     ... etc ...


if __name__ == "__main__":
    main()
