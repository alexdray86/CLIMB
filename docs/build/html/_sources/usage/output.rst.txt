Output 
___________________

*pyTEnrich* writes two tables per input bed file:
    - One table for TE **subfamily** enrichment
    - One table for TE **family** enrichment

They both have similar format :

.. list-table:: pyTEnrich output (bedfile_name_SUBFAM.tsv)
   :widths: 8 8 8 8 8 8 8 8 8
   :header-rows: 1

   * -
     - n_i
     - n_T.te
     - n_T.peak
     - expected
     - fc.expected
     - comparison
     - padj.final
     - significance
   * - SVA_D
     - 934
     - 1433
     - 4395
     - 2.74
     - 250.2
     - peak_in_te
     - 0
     - \*\*\*\*
   * - LTR83
     - 4
     - 724
     - 4395
     - 0.79
     - 2.8
     - te_in_peak
     - 0.773
     - n.s
   * - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...

*Columns interpretations*


      First columns represents TE family/subfamily name (Note that this column is not named in output)

**n_i**
      The overlap found between input bed file and this TE family/subfamily. The size of the overlap rely on Bedtools intersect and the options used for the intersection.

**n_T.te**
      Total number of loci in the TE subfamily/family. 

**n_T.peak**
      Total number of intervals in input bed files. Typically, the total number of peaks in ChIP-seq.

**expected**
      Expected overlap between TE group and input bed file, according to a probability derived from TE genome occupancy.
   
**fc.expected**
      Fold change between the observed overlap **n_i** and the expected one **expected**. 

**comparison**
    Type of comparison that was done for enrichment analysis. Can be either *te_in_peak* or *peak_in_te*. For more explanation consult Details methods section.

**padj.final**
    Adjusted p-value. 

**significance**
    Symbol representing p-value as either

    +--------------+--------------------------------------+
    | **n.s**      |   No significance                    |
    +--------------+--------------------------------------+
    | **\***       | :math:`0.01 < p-val \leq 0.05`       |
    +--------------+--------------------------------------+
    | **\*\***     | :math:`10^{-3} < p-val \leq 10^{-2}` |
    +--------------+--------------------------------------+
    | **\*\*\***   | :math:`10^{-4} < p-val \leq 10^{-3}` |
    +--------------+--------------------------------------+
    | **\*\*\*\*** | :math:`p-val \leq 10^{-4}`           |
    +--------------+--------------------------------------+


