``GPSW`` commands
================================================================================

.. code-block:: shell

   $ gpsw --version
   gpsw 0.6.3

   $ gpsw --help
   usage: gpsw [-h] [--version] {fetch,run} ...

   GPSW: A tool for analysing and processing Global Protein Stability Profiling
   data.

   positional arguments:
     {fetch,run}  Sub-command help

   options:
     -h, --help   show this help message and exit
     --version    show programs version number and exit

   $ gpsw fetch --help
   usage: gpsw fetch [-h] [-t TAG] [--test-data] [-d DIRECTORY]

   Fetch GPSW code from a specific release from https://github.com/niekwit/gps-
   orfeome.

   options:
     -h, --help            show this help message and exit
     -t TAG, --tag TAG     Release tag of the code to download
     --test-data           Include test data in the download (and relating config
                           and resources directories)
     -d DIRECTORY, --directory DIRECTORY
                           Directory to download the code to target directory
                           (default: current directory)

   $ gpsw run --help
   usage: gpsw run [-h] [-p PROFILE] [--snakemake-args SNAKEMAKE_ARGS] [-d] [-q]

   Run the GPSW pipeline and create report.

   options:
     -h, --help            show this help message and exit
     -p PROFILE, --profile PROFILE
                           Path to Snakemake profile YAML file (only has to be
                           provided at first run) (OPTIONAL, use value None for
                           no profile)
     --snakemake-args SNAKEMAKE_ARGS
                           Extra Snakemake arguments (should be '"double-
                           quoted"')
     -d, --dry-run         Dry-run workflow only
     -q, --quiet           Run GPSW quietly