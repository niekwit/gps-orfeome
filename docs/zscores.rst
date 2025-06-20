.. _zscore:

Robust z-score derivation, corrections and scaling
--------------------------------------------------------------------------------

Initial z-score derivation
===========================

``GPSW`` converts the :math:`\Delta\Psi_i` for individual ORFs to robust z-scores. This is done for two reasons:

1. Standardization and Comparability: converting :math:`\Delta\Psi_i` to a z-score standardizes the data, expressing each value in terms of its distance from the central tendency relative to the spread. This allows for meaningful comparisons of a specific metric value across different contexts or datasets, even if their original scales or distributions vary. It also helps in identifying how "unusual" an individual data point is.
2. Robustness to Outliers: while a standard z-score uses the mean and standard deviation, which are sensitive to extreme values, our robust z-score employs the median and median absolute deviation (MAD). The median is resistant to outliers as a measure of central tendency, and the MAD provides a robust estimate of data variability. This ensures that the resulting z-score is a more reliable and stable indicator of deviation, particularly crucial in datasets where outliers are present or the data is not perfectly normally distributed.

The raw robust z-score is calculated as follows:

.. math::
   z_{raw} = \frac{\Delta\Psi_i - \text{median}(\Delta\Psi_j)}{k \times \text{median}(|\Delta\Psi_j - \text{median}(\Delta\Psi_j)|)}


Where:

- :math:`\Delta\Psi_i` represents a single :math:`\Delta\Psi` value for a given ORF :math:`i`.
- :math:`\Delta\Psi_j` represents :math:`\Delta\Psi` values for all ORFs.
- :math:`k` is a standard scaling constant (:math:`1.4826`). It is approximately :math:`1/(\Phi^{-1}(3/4))`, where :math:`\Phi^{-1}` is the inverse of the cumulative distribution function for a standard normal distribution. 

z-score corrections
=====================

We next apply several corrections to the z-scores to account for specific characteristics of the data:

1. Correction for low number of `good barcodes`: ORFs with a low number of `good barcodes` (see :ref:`here <good_barcodes>` for definition) are corrected to avoid skewing the z-score.
2. Intra-ORF variability.
3. Low :math:`\Delta\Psi_i` values

Good barcodes correction
^^^^^^^^^^^^^^^^^^^^^^^^^^

The z-score of ORFs with a low number of `good barcodes` (see :ref:`here <good_barcodes>` for definition) is corrected, as follows:

.. math::

   z_{barcode\_corrected} =
   \begin{cases}
     \frac{z_{raw}}{\sqrt{ \left( 1 + \frac{m - n}{p} \right) }} & \text{if } n < m \\
     z_{raw} & \text{if } n \ge m
   \end{cases}

Where:

- :math:`z_{raw}` is the uncorrected z-score.
- :math:`n` is the number of `good barcodes`.
- :math:`m` is the median of `good barcodes` of all ORFs.
- :math:`p` is a user-defined penalty factor (`penalty_factor` in `config.yaml`).


Intra-ORF variability correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each ORF, the :math:`\Delta\Psi` values are calculated for each individual good barcode. The intra-ORF variability is then assessed by calculating the standard deviation (:math:`\sigma`) of these :math:`\Delta\Psi` values. This standard deviation is then used to correct the z-score for intra-ORF variability:

.. math::

   z_{variability\_corrected} = \begin{cases}
   \frac{z_{barcode\_corrected}}{\sigma_i} & \text{if } \sigma_i > 0 \\
   \frac{z_{barcode\_corrected}}{\epsilon} & \text{if } \sigma_i = 0
   \end{cases}

Where:

- :math:`\sigma_{i}` is the standard deviation of :math:`\Delta\Psi` values of an individual ORF.
- :math:`\epsilon` is the lowest :math:`\sigma_i` of all ORFs.

In rare cases, an ORF may have no variability in its :math:`\Delta\Psi` values (i.e. :math:`\sigma = 0`). In such cases, we apply a small :math:`\epsilon` value to avoid division by zero. :math:`\epsilon` is the lowest :math:`\sigma` of all ORFs, ensuring that the z-score remains defined.


Low :math:`\Delta\Psi_i` values correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some ORFs may have a high z-score but a low :math:`\Delta\Psi_i` value. This can occur if there is very low intra ORF variability, leading to a high z-score despite a low :math:`\Delta\Psi_i`, which can skew the results. To address this, we apply a correction based on the absolute value of :math:`\Delta\Psi_i` and a user-defined threshold for calling a hit:

.. math::

   z_{final} = z_{variability\_corrected} \times \frac{|\Delta\Psi_i|}{h}  

Where:

- :math:`|\Delta\Psi_i|` is the absolute value of :math:`\Delta\Psi` for the individual ORF.
- :math:`h` is a user-defined, absolute, :math:`\Delta\Psi` threshold for calling a hit (defined in `config.yaml`).

Applying this correction ensures that ORFs with very low :math:`\Delta\Psi_i` values are penalised, preventing them from having a disproportionately high z-score.

z-score scaling
=====================

We next scale the z-scores to a range of -128 to -2 for negative z-scores and 2 to 128 for positive z-scores, followed by log2 transformation. As the z-score's direction is important, we are scaling the positive and negative z-scores separately. Scaling the z-scores allows for more consistent plotting of the results, and it ensures that the z-scores are within a defined range, making them easier to compare across different datasets.


Positive z-score scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The scaled positive values, :math:`z_{scaled,pos}`, are calculated and log2 transformed as follows:

.. math::

   z_{scaled,pos} = \log2(L_{pos} + \frac{z_{final,pos} - \min(D_{pos})}{\max(D_{pos}) - \min(D_{pos})} \times (U_{pos} - L_{pos}))

Where:

- :math:`z_{final,pos}` is a positive corrected z-score.
- :math:`D_{pos}` represents all the positive values among all :math:`z_{final}` values.
- :math:`L_{pos}` is the desired lower bound for the scaled positive values (2).
- :math:`U_{pos}` is the desired upper bound for the scaled positive values (128).


Negative z-score scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The scaled negative values, :math:`z_{scaled,neg}`, are calculated and log2 transformed as follows:

.. math::

   z_{scaled,neg} = -\log2(L_{neg} + \frac{(z_{final,neg} - \min(D_{neg}))}{(\max(D_{neg}) - \min(D_{neg}))} \times (U_{neg} - L_{neg}))

Where:

- :math:`z_{final,neg}` is a negative corrected z-score.
- :math:`D_{neg}` represents all the negative values among all :math:`z_{final}` values.
- :math:`L_{neg}` is the desired lower bound for the scaled negative values (-128).
- :math:`U_{neg}` is the desired upper bound for the scaled negative values (-2).

