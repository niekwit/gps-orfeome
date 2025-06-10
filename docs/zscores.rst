z-score calculation
--------------------------------------------------------------------------------

With `bin_number` greater than 1, `GPSW` will perform a protein stability analysis using Protein Stability Index (PSI, :math:`\Psi`) as a metric. This is calculated as follows:

.. math::

   \Psi=\sum_{i=1}^nR_i \times i

where:

- :math:`R_i` is the proportion of the Illumina reads present for an ORF in that given subpopulation :math:`i`.
- :math:`n` is the number of bins.
- :math:`i` is the bin number.

Between two conditions (test and control), the :math:`\Delta\Psi` is calculated as:

.. math::

   \Delta\Psi = \Psi_{test} - \Psi_{control}

Next, :math:`\Delta\Psi` is converted to a robust z-score:

.. math::
   z = \frac{\Delta\Psi_i - \text{median}(\Delta\Psi_j)}{k \times \text{median}(|\Delta\Psi_j - \text{median}(\Delta\Psi_j)|)}


where:

- :math:`\Delta\Psi_i` represents a single :math:`\Delta\Psi` value for a given ORF :math:`i`.
- :math:`\Delta\Psi_j` represents :math:`\Delta\Psi` values for all ORFs.
- :math:`k` is a standard scaling constant (:math:`1.4826`). It is approximately :math:`1/(\Phi^{-1}(3/4))`, where :math:`\Phi^{-1}` is the inverse of the cumulative distribution function for a standard normal distribution. 

The z-score of ORFs with a low number of `good barcodes` (see :ref:`here <good_barcodes>` for definition) is corrected, as follows:

.. math::

   z_{c} =
   \begin{cases}
     \frac{z}{\sqrt{ \left( 1 + \frac{m - n}{p} \right) }} & \text{if } n < m \\
     z & \text{if } n \ge m
   \end{cases}

Where:

- :math:`z_{c}` is the corrected :math:`z`.
- :math:`z` is the z-score.
- :math:`n` is the number of `good barcodes`.
- :math:`m` is the median of `good barcodes` of all ORFs.
- :math:`p` is a user-defined penalty factor (`penalty_factor` in `config.yaml`).

A final z-score correction is applied to correct for intra-ORF variability:

.. math::

   z_{c}' = \begin{cases}
   \frac{z_{c}}{\sigma_i} & \text{if } \sigma_i > 0 \\
   \frac{z_{c}}{\epsilon} & \text{if } \sigma_i = 0
   \end{cases} \times \frac{|\Delta\Psi_i|}{h}

Where:

- :math:`z_{c}'` is the final corrected z-score.
- :math:`\sigma_{i}` is the standard deviation of :math:`\Delta\Psi` values of an individual ORF.
- :math:`h` is a user-defined, absolute, :math:`\Delta\Psi` threshold for calling a hit.
- :math:`|\Delta\Psi_i|` is the absolute value of :math:`\Delta\Psi` for the individual ORF.
- :math:`\epsilon` is the lowest :math:`\sigma_i` of all ORFs (to avoid division by zero).

z-score scaling
--------------------------------------------------------------------------------
We next scaled the z-scores to a range of -128 to -2 for negative z-scores and 2 to 128 for positive z-scores, followed by log2 transformation. As the z-score's direction is important, we are scaling the positive and negative z-scores separately.

Positive z-score scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The scaled positive values, :math:`z_{spos}`, are calculated and log2 transformed as follows:

.. math::

   z_{spos} = \log2(L_{pos} + \frac{z_{corr}' - \min(D_{pos})}{\max(D_{pos}) - \min(D_{pos})} \times (U_{pos} - L_{pos}))

Where:

- :math:`z_{spos}` is the scaled positive z-score.
- :math:`z_{corr}'` is the corrected z-score.
- :math:`D_{pos}` represents all the positive values among all :math:`z_{corr}'` values.
- :math:`L_{pos}` is the desired lower bound for the scaled positive values (2).
- :math:`U_{pos}` is the desired upper bound for the scaled positive values (128).

Negative z-score scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The scaled negative values, :math:`z_{sneg}`, are calculated and log2 transformed as follows:

.. math::

   z_{sneg} = -\log2(L_{neg} + \frac{(z_{corr}' - \min(D_{neg}))}{(\max(D_{neg}) - \min(D_{neg}))} \times (U_{neg} - L_{neg}))

Where:

- :math:`z_{sneg}` is the scaled negative z-score.
- :math:`z_{corr}'` is the corrected z-score.
- :math:`D_{neg}` represents all the negative values among all :math:`z_{corr}'` values.
- :math:`L_{neg}` is the desired lower bound for the scaled negative values (-128).
- :math:`U_{neg}` is the desired upper bound for the scaled negative values (-2).


