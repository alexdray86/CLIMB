# CLIMB
CLIMB performs bulk deconvolution with the full raw count matrix from scRNA-seq, fitting coefficients on each single-cells, aggregated at the cell-type level to produce cell-type proportions.

### Download CLIMB

CLIMB package can be installed directly with `devtools::install_github('alexdray86/CLIMB')` from R. It requies GLMNET and Biobase libraries to work. 

### Run the example notebook 

The notebook present in this folder can be used to run example data in the data/ folder.

### Easy usage : 

Go in your desired folder and git clone CLIMB repository:

`git clone https://github.com/alexdray86/CLIMB.git`

Launch the notebook with jupyter notebook:

`juypyter notebook climb_example.ipynb`

Follow the instructions to install CLIMB and load the data. 
