Global Protein Stability Profiling 
==================================================

Global Protein Stability (GPS) Profiling is a powerful genetic method to study the stability of proteins in a single cell.

It relies on a retroviral reporter construct, which contains a single promoter that drives an expression cassette containing an internal ribosome entry site (IRES). This allows the expression of two fluorescent proteins: DsRed and EGPF fused to a protein or peptide of interest (POI). DsRed fluorescence serves as an internal control, while EGFP fluorescence is dependent on the stability of the POI. Importantly, DsRed and EGFP-POI proteins should be produced at a constant ratio because they are translated from the same mRNA.

The EGFP/DsRed ratio, quantifiable by flowcytometry, serves as a direct readout for POI stability. Perturbations that specifically influence EGFP-POI's stability will predictably alter the EGFP-POI abundance without affecting DsRed levels. This selective change in EGFP-X, but not DsRed, will manifest as a measurable shift in the EGFP/DsRed ratio.

During the GPS profiling experiment, cells are sorted using FACS into multiple bins based on their EGFP/DsRed ratio (see Figure 1). Each bin represents a distinct subpopulation of cells with varying EGFP/DsRed ratios. Cell from these bins are then sequenced to determine the abundance of each ORF in each subpopulation. 


.. figure:: images/psi-flow-plot.png
   :alt: Sorting for PSI analysis

   Figure 1: Typical cell sort strategy for PSI analysis (provided by Hudson Coates)


To quantify the stability of each individual ORF in the experiment, the protein stability index (:math:`\Psi`) is calculated for each ORF across the sorted bins. The :math:`\Psi` is a measure of how the abundance of the EGFP-POI changes relative to DsRed across the different bins, reflecting the stability of the POI. It is calculated as follows:

.. math::

   \Psi=\sum_{i=1}^nR_i \times i

where:

- :math:`R_i` is the proportion of the Illumina reads present for an ORF in that given subpopulation :math:`i`.
- :math:`n` is the number of bins.
- :math:`i` is the bin number.


Between two experimental conditions (e.g., a test condition and a control condition), the difference in protein stability index is computed for each individual ORF:

.. math::

   \Delta\Psi = \Psi_{test} - \Psi_{control}

Negative :math:`\Delta\Psi` values indicate that the ORF is less stable in the test condition compared to the control, while positive values indicate greater stability in the test condition.

:math:`\Delta\Psi` values are generated for each barcode of an individual ORF, after which the mean is calculated, :math:`\Delta\Psi_i`.

In the :ref:`next section <zscore>`, we will discuss how :math:`\Delta\Psi_i` values are converted to robust z-scores, which standardize the data and allow for meaningful comparisons across different datasets.

