Execution
_________

**Load CLIMB package** 

Once CLIMB is installed, you should load CLIMB into your R environment with : 

.. code-block:: R
    
    library(climb)

**Generating input objects for CLIMB**

To be able to use CLIMB, you need to create two **ExpressionSet** object (requires **Biobase** package), with the same shape as requested by MUSIC or BisqueRNA method [ref, ref]. You will need one **ExpressionSet** for your scRNA-seq reference dataset, containing genes as row and cells as columns, with the cell type labels named **cellType**. Similarly, your bulk data should be given as an **ExpressionSet** object, with genes as rows and mixtures as columns. Hereafter, an example to build the **ExpressionSet** object from R matrices: 

.. code-block:: R

    library(climb) ; library(Biobase)
    scRNA.es = ExpressionSet(scRNAseq_count_matrix)
    scRNA.es$cellType = celltype_labels_vector
    bulk.es = ExpressionSet(bulk_count_matrix)

**Minimal run of CLIMB method**

To run CLIMB with default parameters, you can simply launch: 

.. code-block:: R

    library(climb) ; library(Biobase)
    climb.res = climb(scRNA.es, bulk.es)
    celltype_proportions = climb.res$props

By default, CLIMB only predict cell-type proportions and do not predict cell-type expression. If you wish to include cell-type expression deconvolution, just add the parameter `predict_expression=TRUE` when you launch CLIMB : 

.. code-block:: R

    library(climb) ; library(Biobase)
    climb.res = climb(scRNA.es, bulk.es, predict_expression=TRUE)
    celltype_proportions = climb.res$props
    celltype_expression = climb.res$expr.pred



