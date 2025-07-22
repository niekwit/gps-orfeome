Installation of required software
================================================================================

It is recommended to run `GPSW` on a Linux-based system (e.g. Ubuntu).

Make sure you have `Conda <https://docs.conda.io/projects/conda/en/latest/index.html>`_ installed.

Installation of the latest version via Conda
--------------------------------------------------------------------------------

The easiest wasy to install `GPSW` is via Conda, which will handle all dependencies automatically.:

.. code-block:: shell

   $ conda install -c bioconda gpsw


.. dropdown:: Installation of older versions
   :icon: info
   :color: primary
   
   To install versions older than "0.7.0" use the following command:

   .. code-block:: shell

      $ conda create -n gpsw=0.6.3 bioconda::gpsw pandas=2.2.3  pygments=2.19.1 conda=24.7.1

   Otherwise, run:

   .. code-block:: shell

      $ conda create -n gpsw=VERSION

   All available releases can be found on the `GitHub releases page <https://github.com/niekwit/gps-orfeome/releases>`_.


Installation of development version
--------------------------------------------------------------------------------

First, create a yaml file with the environment configuration, e.g. `gpsw.yaml`:

.. code-block:: yaml
   :substitutions:

   name: gpsw
   channels:
      - conda-forge
      - bioconda
      - defaults
   dependencies: 
      - python=3.12
      - pygithub=2.6.1
      - pydot=3.0.4
      - apptainer=1.4.0
      - snakemake-minimal=|snakemake_version|
      - numpy=2.2.6
      - pandas=2.2.3
      - pygments
      - conda=24.7.1

Then, create the environment using the yaml:

.. code-block:: shell
   :substitutions:

   $ conda env create -f gpsw.yaml

.. important::
   If you want to use Apptainer, it is essential to install Snakemake |snakemake_version|, as other versions might not work with the pre-build image.

To install `GPSW`:

.. code-block:: shell

   $ cd /path/to/clone/gpsw
   $ git clone https://github.com/niekwit/gps-orfeome.git
   $ cd gps-orfeome
   $ pip install -e .

At a later point, `GPSW` can be updated:

.. code-block:: shell

   $ git pull

Installation of Snakemake plugins
--------------------------------------------------------------------------------

Snakemake plugins are required for the workflow to run correctly in certain environments. Please follow the instructions `here <https://snakemake.github.io/snakemake-plugin-catalog/index.html>`_ to install the required plugins.