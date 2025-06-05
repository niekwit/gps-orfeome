Installation of required software
================================================================================

It is recommended to run `GPSW` on a Linux-based system (e.g. Ubuntu).

Make sure you have `Conda <https://docs.conda.io/projects/conda/en/latest/index.html>`_ installed.

Installation of stable version via Conda
--------------------------------------------------------------------------------
.. code-block:: shell

   $ conda create -n gpsw bioconda::gpsw pandas=2.2.3  pygments=2.19.1

This will install the stable version of `GPSW` and all of its dependencies. This is the recommended way to install `GPSW`.

Installation of development version
--------------------------------------------------------------------------------
First, create a Conda env with the dependencies:

.. code-block:: shell

   $ conda create -n gpsw snakemake=8.25.5 apptainer=1.4.0 pandas=2.2.3  pygments=2.19.1

.. note::
   If you want you use Apptainer, it is essential to install Snakemake v8.25.5, as later versions might not work with the pre-build image.

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
For running `GPSW` on