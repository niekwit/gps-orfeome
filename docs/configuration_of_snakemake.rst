Configuration of Snakemake
================================================================================

To set up a profile for custom `Snakemake` command line arguments, create a new profile (`config.yaml`) in ``$HOME/.config/snakemake/standard/``:

.. code-block:: yaml

   cores: 40
   latency-wait: 10
   use-conda: True
   rerun-incomplete: True
   printshellcmds: True
   show-failed-logs: True
   use-apptainer: True