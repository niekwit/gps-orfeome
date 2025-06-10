Run ``GPSW`` with your own data
================================================================================

Fetching workflow code
--------------------------------------------------------------------------------

Download workflow code as follows:

.. code-block:: shell

   $ cd /path/to/analysis/dir
   $ gpsw fetch


Configuring GPSW
--------------------------------------------------------------------------------

To configure GPSW, edit the `config/config.yaml` file. This file contains various parameters that control the workflow's behavior, as explained in the :doc:`workflow settings <workflow_settings>` section


Running the workflow
--------------------------------------------------------------------------------
To initiate the workflow, run the following command:

.. code-block:: shell

   $ gpsw run
