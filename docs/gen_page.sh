#!/bin/bash

# 1) Install sphinx - make sure to do it on a dedicated python environment
pip install -U Sphinx
pip install sphinx_rtd_theme # install good looking theme

# 2) Activate github pages on git, for the folder ./docs

# 3) rename all folder to remove trailing '_'. e.g. _build -> build 
