## Using Bayesian statistical learning to integrate temperature dependence in enzyme-constrained GEMs
<p align="center">
  <img  src="figures/logo.png">
</p>

#### Description of folders
* `code/` contains all scripts and detailed descrition can be found in `code/README.md`.
* `data/` contains all input data needed, including experimental and estimated thermal parameters.
* `models/` contains a list of yeast genome scale models with different settings used in this study.
* `validate_smc_abc/` contains the scripts and simulation results for the validation of SMC-ABC apporach on toy models.

#### Dependences
```
numpy                       1.15.0  
pandas                      0.23.4
scikit-learn                0.20.3
seaborn                     0.9.0
jupyter                     1.0.0
cobra                       0.15.3  
Gurobi                      8.0.0
```
The repository was tested with Python 3.6.7.

#### Hardware
Since Bayesian approach is computational expsensive, all scripts except ones for visualizaiton have to be done with a computer cluster. Those scripts have been designed for parallel computation. The SMC-ABC approach takes around 3-5 days on a compute node with 32 cores (Intel Xeon Gold 6130 CPU). All visualization scripts take sveral seconds or minutes on a normal PC.

#### Reproduce the figures
(1) Clone this repository.  
(2) Install all required packages. This step takes at most several minutes.  
(2) Download the pre-computed results from Zenodo (https://zenodo.org/record/3996543#.X0J1BNP7S3I). Download the `results.tar.gz` file to the current directory and uncompress with 
```
tar -xzvf results.tar.gz
```
Then the figures in the manuscript can be reproduced by using jupyter notebooks
* `visualization.ipynb` : **Figs 2abcdefi, 3ab; 5abde; S5, S6, S7, S10, S11a**
* `simulate_with_start_params.ipynb` **Fig S2abc**
* `additional_plots.ipynb`: - additional plots (**Figs 2gh, S13**)
* `visualize_cv.ipynb`: **Fig S3, S4**
* `visualize_temperature_on_enzymes_posterior.ipynb`: **Fig 3cdef**
* `visualize_chemostat_metabolic_shift.ipynb`: **Fig 4abc, S9**
* `expdata.ipynb`: **Fig 5c, S11b**
* `single_enzyme.ipynb`: **Fig S8**
* `Case1.ipynb`~`Case5.ipynb`: **Fig S16-S20**

  
Results for experimental validation of ERG1 is available as csv files:
```
data/OD600_different_passages_40C.csv       
data/OD600_different_passages_42C.csv 
```

One can also recompute those results by following the introductions in `code/README.md` and the visualized by the above scripts.
