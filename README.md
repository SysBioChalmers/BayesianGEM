## Bayesian genome scale modeling reveals thermal determinants of yeast metabolism

#### Description
This repository contains code and data needed to replicate the analysis carried out in the manuscript Li G, et al. 2020.
* `code/` contains all scripts and detailed descrition can be found in `code/README.md`
* `data/` contains all input data needed, including experimental and estimated thermal parameters
* `models/` contains a list of yeast genome scale models with different settings used in this study.

All the results can be downloaded from *ZENODO*

#### Dependences
```
numpy                       1.15.0  
pandas                      0.23.4
scikit-learn                0.20.3
seaborn                     0.9.0
jupyter                     1.0.0
cobra                       0.15.3  
```
The repository was tested with Python 3.6.7.

#### Hardware
Since Bayesian approach is computational expsensive, all scripts except ones for visualizaiton have to be done with a computer cluster. Those scripts have been designed for parallel computation. 