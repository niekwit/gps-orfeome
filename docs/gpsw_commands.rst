``GPSW`` commands
================================================================================

Version
--------

To display the `GPSW` version:

.. parsed-literal::

   $ gpsw --version
   gpsw |version|

Help
--------

To display the `GPSW` help message, use the ``--help`` option. This will show you the available sub-commands and options for `GPSW`. The help message provides a brief overview of the tool and its functionalities:

.. code-block:: shell

   $ gpsw --help
   usage: gpsw [-h] [--version] {fetch,run} ...

   GPSW: A tool for analysing and processing Global Protein Stability Profiling
   data.

   positional arguments:
     {fetch,run}  Sub-command help

   options:
     -h, --help   show this help message and exit
     --version    show programs version number and exit


`GPSW` fetch
--------------------

The sub-command ``fetch`` is used to download the GPSW code and analysis-specific configuration files from a specific release on GitHub:

.. code-block:: shell

   $ gpsw fetch --help
   usage: gpsw fetch [-h] [-t TAG] [--test-data] [-d DIRECTORY]

   Fetch GPSW code from a specific release from https://github.com/niekwit/gps-
   orfeome.

   options:
     -h, --help            show this help message and exit
     -t TAG, --tag TAG     Release tag of the code to download (default: latest)
     --test-data           Include test data in the download (and relating config
                           and resources directories)
     -d DIRECTORY, --directory DIRECTORY
                           Directory to download the code to target directory
                           (default: current directory)

.. note::
   If you specify the ``TAG`` as anything other than `latest` and `0.6.3` or earlier, make sure that the `gpsw` conda environment was made with the same version of `GPSW`. For example, if you want to use the `0.6.3` release, you should create the conda environment with the command:
   
   .. code-block:: shell

      $ conda create -n gpsw bioconda::gpsw=0.6.3 pandas=2.2.3 pygments=2.19.1


`GPSW` run
--------------------

The sub-command ``run`` is used to run the GPSW pipeline and create a report:

.. code-block:: shell

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

