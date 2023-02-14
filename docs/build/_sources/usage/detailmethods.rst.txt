.. _detailmethods:

CLIMB method
________________

**CLIMB deconvolution for cell-type abundance prediction**

CLIMB method proposes an innovative approach to solve a bulk deconvolution problem. Unlike conventional methods that fit coefficients for each cell-type based on a signature matrix (e.g., cell-type average expression), CLIMB fits coefficients, denoted as :math:`\alpha_{i}`, for each individual single-cell :math:`i`, with the goal of finding the optimal linear combination of single-cell expression vectors :math:`c_i` that can predict the bulk expression vector :math:`\mathbf{x}_n`. Thus, their relationship is defined by :math:`\mathbf{x}_n \sim \mathbf{C} \cdot \boldsymbol{\alpha}_n`, shown below in more details:

.. math:: \begin{bmatrix}
        x_{n1} \\
        \vdots \\
        x_{nG}
    \end{bmatrix}
    \propto
    \underbrace{ \begin{bmatrix}
        c_1^1 & ... & c_1^C \\
        \vdots & \ddots & \vdots \\
        c_G^1 & ... & c_G^C
    \end{bmatrix} }_{ \text{ $\mathbf{C}^{(G \times C)}$ } }
    \cdot
    \begin{bmatrix}
        \alpha_n^1 \\
        \vdots \\
        \alpha_n^C
    \end{bmatrix}

where the vector :math:`\mathbf{x}_n` represent the bulk expression matrix for a mixture :math:`n`, :math:`\mathbf{C}` is the raw counts matrix from scRNA-seq, and :math:`\boldsymbol{\alpha}_n` is the vector of coefficients fitted for each single-cell :math:`c` for a given mixture :math:`n`. The coefficients :math:`\alpha_i` are fitted with `glmnet` library with a non-negativity constraint on the coefficients, and no regularization. Thus, CLIMB loss function is defined by :

.. math:: L(\boldsymbol{\alpha}_n \vert \mathbf{C}, \mathbf{x}_n ) =  \min\limits_{\boldsymbol{\alpha}} \frac{1}{2G} \left\vert \mathbf{x}_n - \mathbf{C} \cdot \boldsymbol{\alpha}_n \right\vert^2 \;\;\text{, with}\;\; \alpha_{ni} \in [0,+\infty]

We then group and normalize :math:`\alpha_{ni}` coefficients at the cell-subtype level to obtain the cell-subtype proportions :math:`w_{nk}` for each cell-subtype :math:`k` and mixture :math:`n` :

.. math:: w_{nk} = \frac{\sum_{i=1}^{C} \tau_{ki} \cdot \alpha_{ni} / \theta_i }{\sum_{i=1}^{C} \alpha_{ni} / \theta_i }

where :math:`\boldsymbol{\tau}_k` is an indicator vector of length :math:`C` and equal to 1 whenever cell :math:`i` is of type :math:`k`, and :math:`\theta_i` denotes the total expression (e.g., the sum of RNA counts) in cell :math:`i`. 

**CLIMB deconvolution of bulk expression into cell-type expression** 

CLIMB allows the prediction of cancer cell expression in a two step manner. Let first define our target *high-resolution* expression :math:`\hat{\mathbf{S}}_{(ngk)}`, which is the cell-type-specific and bulk sample-specific expression for each cell-type :math:`k`. We assume that for a given group of samples, the *grouped* expression :math:`\hat{\mathbf{S}}_{(gk)}` relates to the *high-resolution* expression with the following formula:`

.. math:: s_{gk}^{(grouped)} = \frac{1}{N} \sum_{n=1}^N s_{ngk}^{(high-res)}

First, we take advantage of the bulk to single-cell mapping to directly derive the cell-type expression from the bulk-to-single-cell deconvolution or *mapping*. We can simply sum up the expression of cells mapped by CLIMB to get a *high-resolution* cell-type-specific and bulk-specific expression :math:`s_{ngk*^{(mapping)*`:

.. math:: s_{ngk}^{(mapping)} = \sum_{i=1}^C \alpha_{in} \cdot c_{ig} \cdot \tau_{ik}

with :math:`\alpha_{in}` being the coefficient fitted by GLMNET on each single-cell, linking each cell :math:`i` to each bulk :math:`n`. :math:`c_{ig}` denotes the single-cell expression of cell :math:`i` and gene :math:`g`, and  :math:`\tau_{ik}` is an indicator vector equal to 1 if cell :math:`i` is of type :math:`k`.

Although this is potentially sufficient in cases where scRNA-seq reference dataset maps all cell states present in bulk samples, it can be insufficient to detect genes deferentially-expressed (DE) at cell-type level. As DE genes will likely induce error on the cell-type expression deconvolution, we use the error itself to fit our model. We can also make the reasonable assumption that DE genes only originates in cancer cells. Thus, for a given gene :math:`g` and a given bulk :math:`n`, we have:

.. math:: y_{ng} \sim \hat{y}_{ng}^{(mapping)} + \epsilon_{ng} \;\;\; \rightarrow \;\;\; y_{ng} - \hat{y}_{ng}^{(mapping)} = \epsilon_{ng}

Thus, we use :math:`\epsilon_{ng} = y_{ng} - \hat{y}_{ng}^{(mapping)}` as input for a model to predict cancer cell-type expression. We decompose our error :math:`\epsilon_{ng}` between the error originated from DE genes :math:`\epsilon_{ng}^{(DE)}`, and an additional error term :math:`\epsilon_{ng}^*`. Similar as other methods [ref, ref], we use the information from all bulk samples :math:`N` to deconvolute the high-resolution cell-type specific expression :math:`\hat{s}_{ngk}^{(DE)}`. We consider that :math:`\epsilon_{ng}^{(DE)}` only originates from cancer cells :math:`C^{(C)}` of type :math:`K^{(C)}`, thus

.. math:: \epsilon_{ng} \sim \underbrace{\sum_{k=1}^{K^{(C)}} \left( w_{nk} \cdot \hat{s}_{ngk}^{(DE)} \right)}_{\epsilon_{ng}^{(DE)}} + \; \epsilon_{ng}^* = \sum_{i=1}^{C^{(C)}} \left( \alpha_{in} \cdot \hat{c}_{ig} \right) + \epsilon_{ng}^* = \mathbf{A} \cdot \hat{\mathbf{c}}_g + \epsilon_{ng}^*

Where :math:`\hat{\mathbf{c}}_g` is the vector of coefficients we wish to predict. Our final prediction of cell-type expresion for a cancer cell-type :math:`k`, a given bulk :math:`n`, and a given gene :math:`g` is:

.. math:: \hat{s}_{ngk}^{(high-res)} = \underbrace{\hat{s}_{ngk}^{(mapping)}}_{\text{from scRNA-seq}} + \underbrace{\hat{s}_{ngk}^{(DE)}}_{\text{learned}} = \sum_{i=1}^C \alpha_{in} \cdot ( c_{ig} + \hat{c}_{ig} ) \cdot \tau_{ik}

Following the assumptions made above, we set :math:`\hat{c}_{ig}` to zero for normal cells, which imply that :math:`\hat{s}_{ngk}^{(high-res)} = \hat{s}_{ngk}^{(mapping)}` for normal cell types.

**CLIMB-BA extension**

In the context of Acute Myeloid Leukemia, blast count is a commonly used metric in clinical practice and has been demonstrated to provide an accurate estimation of the proportion of cancer cells in a given mixture, as documented by Van Galen et al. In this work, we introduce a variable, denoted as :math:`b_n`, to quantify the blast count in sample :math:`n`. Additionally, we present an estimate of the cancer cell proportion through bulk deconvolution, defined as the sum of the proportions of each cancer cell subtype.

.. math:: \hat{b}_n = \sum_{k \in \kappa} w_{nk}

If we expect this proportion to be more precise than the estimation given by bulk deconvolution alone, cancer cell proportion can be normalized with the following approach :

.. math::
     w_{nk}^*=
        \begin{cases}
        w_{nk}  \; \frac{b_n + \epsilon}{\hat{b}_n + \epsilon},& \text{if } \;\; k \in \kappa\\
        w_{nk} \; \frac{1 - b_n + \epsilon}{1 - \hat{b}_n + \epsilon},  & \text{otherwise}
        \end{cases}

with :math:`b_n` being the blast count for a mixture :math:`n`, and :math:`\hat{b}_n` the cancer cell proportion estimated with CLIMB bulk deconvolution. :math:`\epsilon` is a constant factor that prevents division by zero. By modulating :math:`\epsilon`, the strength of the normalization can be controled, with the two following extreme situations :

.. math:: w^*_{nk} \to w_{nk} \frac{b_n}{\hat{b}_n} \;\;, \;\; \text{if } \;\; \epsilon \to 0

.. math:: w^*_{nk} \to w_{nk} \;\;, \;\; \text{if } \;\; \epsilon \to \infty

respectively :math:`w^*_{nk} \to w_{nk} (1 - b_n) / (1-\hat{b}_n)`, if :math:`\epsilon \to 0` and :math:`w^*_{nk} \to w_{nk}` if :math:`\epsilon \to \infty`, for normal cell subtypes. By increasing :math:`\epsilon` we are doing a smooth normalization that turned out to be beneficial when :math:`\epsilon=1`.

Note that we also tried to correct cancer cell-subtype proportions through a regularization penalty term but the latter was not able to reach better precision than the approach described here. As a simple normalization is easier to implement and faster to compute, we decided to choose this approach.


