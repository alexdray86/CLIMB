Installation
____________

**Install CLIMB from R**

From R, you can use devtools to directly install **CLIMB** with the following command: 

.. code-block:: R

    devtools::install_github('alexdray86/CLIMB')

In R, you will need two libraries, GLMNET and Biobase. GLMNET is used to optimize our linear regression model, and Biobase is used for the `ExpressionSet` format, that is used as input for CLIMB (same as for MUSIC method). Both method should be installed automatically upon CLIMB installation. You will still need to load them in R before being able to use CLIMB. Thus, you will have to load the following libraries in R :

.. code-block:: R

    library(climb) ; library(glmnet) ; library(Biobase)

**Install CLIMB from GitHub**

Similarly, you can git clone our GitHub repository with the following command : 

.. code-block:: bash

    git clone https://github.com/alexdray86/CLIMB.git

For running CLIMB the first time, we advice to run the example notebook with mock data provided in the GitHub repository above. 

Enjoy !

