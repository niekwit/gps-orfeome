import os, sys, re
import subprocess
import configparser
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError
from collections import defaultdict

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


def extract_rule_categories(rules_dir="workflow/rules"):
    """
    Extract rule-to-category mappings from all .smk files in a rules directory.
    Looks for '# category: ...' comments immediately above 'rule <name>:' lines.

    Returns:
        dict: {rule_name: category}
    """
    rule_categories = {}
    current_category = "Other"

    for root, _, files in os.walk(rules_dir):
        for filename in files:
            if filename.endswith(".smk"):
                filepath = os.path.join(root, filename)
                with open(filepath, "r") as f:
                    for line in f:
                        line_stripped = line.strip()
                        if line_stripped.startswith("# category:"):
                            current_category = line_stripped.split(":", 1)[1].strip()
                        elif line_stripped.startswith("rule "):
                            m = re.match(r"rule\s+(\w+)\s*:", line_stripped)
                            if m:
                                rule_name = m.group(1)
                                rule_categories[rule_name] = current_category
                                current_category = "Other"  # Reset
    return rule_categories


def create_rule_graph():
    """
    Generate a visual graph of Snakemake workflow rules and their dependencies.

    This function runs Snakemake with the `--rulegraph` option to obtain a
    directed graph of all rules in the workflow. It captures and processes
    the raw dot-format output, filtering out edges and nodes related to the
    special "all" target node to simplify the graph.

    The filtered graph definition is saved as a DOT file (`images/rulegraph.dot`),
    which is then parsed to generate a PDF visualization (`images/rulegraph.pdf`).

    The resulting graph helps visualize the relationships and dependencies
    between workflow rules, aiding in debugging and workflow understanding.

    The function ensures the output directory exists and handles subprocess
    execution, file writing, and graph rendering.

    Raises:
        subprocess.CalledProcessError: If the Snakemake command fails.
        OSError: If file or directory operations fail.
    """

    os.makedirs("images", exist_ok=True)

    # Get category mapping from Snakefile
    rule_categories = extract_rule_categories("workflow/rules")

    # Run Snakemake and get DOT
    result = subprocess.run(
        ["snakemake", "--quiet", "all", "--forceall", "--rulegraph"],
        capture_output=True,
        text=True,
        check=True,
    )

    # Filter out rule all
    lines = [
        line
        for line in result.stdout.splitlines()
        if not ("-> 0" in line or '0[label = "all"' in line)
    ]

    # Parse nodes and edges
    node_lines = []
    edge_lines = []
    for line in lines:
        if "->" in line:
            edge_lines.append(line)
        elif re.match(r"\s*\d+\[label", line):
            node_lines.append(line)

    # Extract rule name and ID from node lines
    rule_to_id = {}
    id_to_rule = {}
    for line in node_lines:
        m = re.match(r'\s*(\d+)\[label\s*=\s*"([^"]+)"', line)
        if m:
            rule_id, rule_name = m.groups()
            rule_to_id[rule_name] = rule_id
            id_to_rule[rule_id] = rule_name

    # Group rules into categories
    clusters = defaultdict(list)
    for rule_id, rule_name in id_to_rule.items():
        category = rule_categories.get(rule_name, "Other")
        clusters[category].append((rule_id, rule_name))

    # Assign background colors to clusters
    palette = [
        "#F0FAF0",
        "#FAF0F0",
        "#F0F0FA",
        "#FFF8E1",
        "#E6F0FF",
        "#FDEDEC",
        "#F3F3F3",
    ]
    category_colors = {
        category: palette[i % len(palette)]
        for i, category in enumerate(sorted(clusters))
    }

    # Write new DOT
    dot_path = "images/rulegraph.dot"
    with open(dot_path, "w") as f:
        f.write("digraph snakemake_dag {\n")
        f.write("    rankdir=TB;\n")
        f.write("    bgcolor=white;\n")
        f.write("    margin=0.2;\n")
        f.write(
            '    node [shape=box, style="rounded,filled", fontname=Helvetica, fontsize=10, penwidth=1.5, color=black, fillcolor=white];\n'
        )
        f.write('    edge [color="#888888", penwidth=1.2];\n\n')

        for category, nodes in clusters.items():
            fill = category_colors.get(category, "#F5F5F5")
            safe_id = re.sub(r"\W+", "_", category.lower())
            f.write(f"    subgraph cluster_{safe_id} {{\n")
            f.write(f'        label = "{category}";\n')
            f.write('        style = "filled,dashed";\n')
            f.write("        color = lightgrey;\n")
            f.write(f'        fillcolor = "{fill}";\n')
            for rule_id, rule_name in nodes:
                f.write(f'        {rule_id} [label="{rule_name}"];\n')
            f.write("    }\n\n")

        for line in edge_lines:
            f.write("    " + line + "\n")

        f.write("}\n")

    # Optional: create PDF
    graphs = pydot.graph_from_dot_file(dot_path)
    graphs[0].write_pdf("images/rulegraph.pdf")


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
