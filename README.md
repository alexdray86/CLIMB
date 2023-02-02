# CLIMB
CLIMB dissects the cellular composition of bulk samples by finding the best combination of single cells using a scRNA-seq dataset as a reference to reconstruct bulk expression from a bulk RNA sequencing target sample. The learned single-cell to bulk sample mapping can then be grouped to obtain cell-subtype proportion estimates, and cell-type gene expression prediction, for each target mixture.

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
