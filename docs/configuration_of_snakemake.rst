Configuration of Snakemake
================================================================================

To set up a profile for custom `Snakemake` command line arguments, create a new profile (`config.yaml`) in ``$HOME/.config/snakemake/standard/``:

.. code-block:: yaml

   cores: 40 # Maximum number of cores to use
   latency-wait: 10
   use-conda: True
   rerun-incomplete: True
   printshellcmds: True
   show-failed-logs: True
   use-apptainer: True

It is recommend to set both ``use-conda`` and ``use-apptainer`` to True, as this will ensure that the workflow runs in a reproducible environment.

When running `GPSW` on a cluster it is recommended to also add the following settings:

.. code-block:: yaml

   executor: slurm # for SLURM cluster
   jobs: 100 # Maximum number of jobs to run in parallel
   apptainer-args: "--bind '/rds/gps/'" # if analysis in not in /home/$USER
   local-cores: 2 # Limit core usage for local rules
   default-resources:
         slurm_partition: <PARTITION>
         slurm_account: <ACCOUNT>

.. note::
   If you are running the analysis using `Apptainer` outside of your home directory, you need to bind the directory where the analysis is located. For example, if your analysis is in `/rds/gps/`, you can use `apptainer-args: "--bind '/rds/gps/'"`.