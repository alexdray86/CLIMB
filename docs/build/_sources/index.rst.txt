.. pyTEnrich documentation master file, created by
   sphinx-quickstart on Wed Dec 30 11:29:02 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: images/CLIMB_logo_design.png

.. toctree::
   :maxdepth: 1
   :caption: Contents: 
   
   usage/installation.rst
   usage/execution.rst
   usage/detailmethods.rst
   usage/results.rst
   
Overview of CLIMB method
=======================

**CLIMB** employs an innovative computational approach to unravel the cellular composition and expression patterns of bulk samples, utilizing a reference scRNA-seq dataset as a guide for reconstructing bulk expression. The resultant single-cell to bulk sample mapping is used to deconvolute bulk samples into cell-subtype abundance and cell-subtype expression, providing valuable insight into the underlying cellular heterogeneity of the samples.

.. image:: images/CLIMB_sketch1.png

**CLIMB** differentiates itself from other methods by not relying on a signature matrix for bulk deconvolution. Instead, it fits coefficients at the single-cell level and then aggregates them at the cell-type level. For further information on the CLIMB method, see `CLIMB method` section.
