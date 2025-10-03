## Repository for "Metallicity fluctuation statistics II" ##

This repository contains the scripts required to carry out all the analysis and produce all the plots included in "Metallicity fluctuation statistics in the interstellar medium and young stars -- II. Elemental cross-correlations and the structure of chemical abundance space", by Mark Krumholz, Yuan-Sen Ting, Zefeng Li, Chuhan Zhang, Jennifer Mead, and Melissa Ness.

### Contents ###

The main script is the file `process_data.py`, in the main directory. This script calls all the remaining python code and processes the data; its output will be placed in a directory called `output`.

The remaining contents are:

* `asplund09_abd.py`: this file holds data on Solar abundances from [Asplund et al. (2009, ARA&A, 47, 481)](http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A)
* `seitenzahl13_yld.py`: this file hold data taken from the yield tables of [Seitenzahl et al. (2013, MNRAS, 429, 1156)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.429.1156S)
* `slug`: this directory contains output files produced by the [slug](https://bitbucket.org/krumholz/slug2) stellar population synthesis code; the `.param` files within this directory are the `slug` parameter files used to generate the data, and the `.fits` files are the outputs of `slug`.
* `stellar_data`: this directory contains tabulations of stellar element-element correlations extracted from [Ting & Weinberg (2022, ApJ, 927, 209)](https://ui.adsabs.harvard.edu/abs/2022ApJ...927..209T) and [Mead, De La Garza, & Ness (2025, arXiv:2504.18532)](https://ui.adsabs.harvard.edu/abs/2025arXiv250418532M). Descriptions of how these data are extracted are provided in the paper.

### Dependencies ###

Running the `process_data.py` script requires the following python packages:

* [numpy](https://numpy.org/)
* [scipy](https://scipy.org/)
* [matplotlib](https://matplotlib.org/)
* [astropy](https://www.astropy.org/) 
* [slugpy](https://bitbucket.org/krumholz/slug)

### License ###

This code is distributed under the terms of the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) version 3.0. The text of the license is included in the main directory of the repository as `LICENSE`.
