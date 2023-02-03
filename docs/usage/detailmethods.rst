.. _detailmethods:

CLIMB method
________________

CLIMB method proposes an innovative approach to solve a bulk deconvolution problem. Instead of fitting coefficients for each cell-type using a signature matrix (e.g. cell-type average expression), CLIMB fits coefficients :math:`\alpha_{i}` for each single-cell :math:`i`, to find the best linear combination of single-cell expression vector :math:`c_i` to predict bulk expression vector :math:`\mathbf{x}_n`. Thus, we use the relation :math:`\mathbf{x}_n \sim \mathbf{C} \cdot \boldsymbol{\alpha}_n`, shown below in more details :

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

where the vector :math:`\mathbf{x}_n` represent the bulk expression matrix for a mixture :math:`n`, :math:`\mathbf{C}` is the raw counts matrix from scRNA-seq, and :math:`\boldsymbol{\alpha}_n` is the vector of coefficients fitted for each single-cell :math:`c` for a given mixture :math:`n` with GLMNET library. 

We then group and normalize these coefficients :math:`\alpha_{ni}` at the cell-subtype level to obtain the cell-subtype proportions :math:`w_{nk}` for each cell-subtype :math:`k` and mixture :math:`n` :

.. math:: w_{nk} = \frac{\sum_{i=1}^{C} \tau_{ki} \cdot \alpha_{ni} / \theta_i }{\sum_{i=1}^{C} \alpha_{ni} / \theta_i }

where :math:`\boldsymbol{\tau}_k` is an indicator vector of length :math:`C` and equal to 1 whenever cell :math:`i` is of type :math:`k`. 

Thus, CLIMB loss function is defined by :

.. math:: L(\boldsymbol{\alpha}_n \vert \mathbf{C}, \mathbf{x}_n ) =  \min\limits_{\boldsymbol{\alpha}} \frac{1}{2G} \left\vert \mathbf{x}_n - \mathbf{C} \cdot \boldsymbol{\alpha}_n \right\vert^2 \;\;\text{, with}\;\; \alpha_{ni} \in [0,+\infty]

Thus, we impose a non-negativity constraint on the coefficients, and no regularization.

**CLIMB-BA extension**

In the context of Acute Myeloid Leukemia, blast count is a routinely measured metric in the clinics. It provides an accurate estimation of cancer cell proportions in a given mixture, as shown by Van Galen et al. First, we introduce a variable for measure blast count in a sample :math:`n`, :math:`b_n`. Then, we introduce the estimated proportion of cancer cells through bulk deconvolution. The cancer cell proportion is the sum of the proportion of each cancer cell subtype :

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


