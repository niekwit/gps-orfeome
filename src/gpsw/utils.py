import os, sys
import subprocess
import configparser
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

import pydot
from github import Github


def get_package_version():
    try:
        return version("gpsw")
    except PackageNotFoundError:
        # Package is not installed (e.g., running locally without installation)
        # For development purposes
        return "v0.0.0-dev"


def get_latest_release_tag():
    try:
        g = Github()

        repository = g.get_user("niekwit").get_repo("gps-orfeome")
        latest_release = repository.get_latest_release()
        return latest_release.tag_name
    except Exception as e:
        print(f"Error fetching release using PyGithub: {e}")
        return None


def dry_run(args):
    """
    Run a dry-run of the workflow using Snakemake.
    """
    print("Performing dry-run of the workflow")
    command = ["snakemake", "-n"]
    if not args.quiet:
        command.append("-p")
    else:
        command.extend(["--quiet", "all"])

    # Run command and check for errors
    try:
        subprocess.run(command, check=True)
        print("Dry-run completed successfully!")
        sys.exit(0)
    except subprocess.CalledProcessError as e:
        # If it was run quietly, we need to check the error message
        # so run dry-run again without --quiet
        if args.quiet:
            subprocess.run(["snakemake", "-np"])
        sys.exit(1)


def create_rule_graph():
    """
    Generates a rule graph for a Snakemake workflow and saves it as a PDF.
    This function creates an 'images' directory if it does not exist, runs Snakemake to generate a rule graph in DOT format,
    filters out specific lines from the output, and then uses pydot to convert the DOT file into a PDF.
    Requirements:
        - Snakemake must be installed and accessible from the command line.
        - pydot must be installed.
        - The function should be run from the root directory of the Snakemake workflow.
    Outputs:
        - images/rulegraph.dot: The DOT file representing the rule graph.
        - images/rulegraph.pdf: The PDF visualization of the rule graph.
    """
    # Create images directory if it doesn't exist
    os.makedirs("images", exist_ok=True)

    # Create the rule graph
    dot_file = "images/rulegraph.dot"
    result = subprocess.run(
        ["snakemake", "--quiet", "all", "--forceall", "--rulegraph"],
        capture_output=True,
        text=True,
        check=True,
    )
    # Remove "rule all" connections:
    # Filter out lines containing "-> 0" or "0[label = "all"]"
    filtered_lines = [
        line
        for line in result.stdout.split("\n")
        if not ("-> 0" in line or '0[label = "all"' in line)
    ]
    with open(dot_file, "w") as f:
        f.write("\n".join(filtered_lines))

    # Use pydot to create a pdf
    graphs = pydot.graph_from_dot_file(dot_file)
    graph = graphs[0]
    graph.write_pdf("images/rulegraph.pdf")


def profile_arg(args):
    """
    Determine the Snakemake profile argument based on configuration.

    Args:
        args: Command line arguments containing profile information

    Returns:
        list: Profile arguments for Snakemake or empty list if no profile
    """

    if args.profile == "None":
        # No profile to be used with Snakemake
        return [""]

    # Check if GPSW configuration file exists
    config_dir = Path.home() / ".gpsw"
    config_file_path = config_dir / "config.ini"

    if os.path.exists(config_file_path):
        config = configparser.ConfigParser()
        config.read(config_file_path)

        # Check if the profile is set in the config file
        if "profile" in config["DEFAULT"]:
            profile = config["DEFAULT"]["profile"]

        else:
            # Fetch profile and store it in the config file
            # This should not be needed, but just in case
            if args.profile is not None:
                profile = args.profile
                config["DEFAULT"] = {
                    "profile": profile,
                }
                with open(config_file_path, "w") as configfile:
                    config.write(configfile)
        return ["--profile", profile]

    else:
        if args.profile is not None:
            os.makedirs(config_dir, exist_ok=True)
            profile = args.profile

            # Create the config file if it doesn't exist
            config = configparser.ConfigParser()
            config["DEFAULT"] = {
                "profile": profile,
            }
            with open(config_file_path, "w") as configfile:
                config.write(configfile)
            return ["--profile", profile]
        else:
            raise ValueError(
                "No profile provided at first run of GPSW. Please provide a Snakemake profile."
            )
