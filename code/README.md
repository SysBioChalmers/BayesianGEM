### Description
This folder contains the scripts that carry out all analysis reported in the paper.

#### Fundmental scripts
* `etcpy` - the script for incorporating temperature parameters into the enzyme constrained Yeast 7.6 model, the detail introduction can be found from `etcpy/README.md`.
* `GEMS.py` - the function for simulating the aerobic, anaerobic growth rate in batch cultivation and aerobic fluxes in chemostat cultivation, as well as a list of distance functions used for SMC-ABC approach.
* `abc_etc.py` - the function to perform SMC-ABC approach.

#### Simulations with initial parameters without considering uncertainties
* `simulate_with_start_params.ipynb` - script contains the function for visualization (**Fig S2**), this script can run run on PC or laptop.

#### SMC-ABC approach updates model parameters
Three experimental datasets under different temperatures are used in this section:
- `ExpGrowth.tsv` - the maximal specific growth rate in aerobic (Caspeta L., et al. Mbio, 2015) and anaerobic (Zakhartsev M., et al. J. Therm. Biol., 2015) batch cultivations

- `Chemostat_exp_data.txt` - fluxes of carbon dioxide (CO2), ethanol and glucose in chemostat cultivations (Postmus J., J. Biol. Chem., 2008)  

- `model_enzyme_params.csv` - 

Cross-validation scripts:
* `gem_smcabc_at_three_conditions_cv1.py`
* `gem_smcabc_at_three_conditions_cv2.py`
* `gem_smcabc_at_three_conditions_cv3.py`
The results can be visualized with `visualize_cv.ipynb`. (**Figs S3, S4**)  
* `gem_smcabc_at_three_conditions.py` - SMC-ABC approach update parameter space with all three observed datasets.

Above scripts need be run on high-performance cluster, and may take a few days.

* `visualization.ipynb` - the analysis of resulted Posterior models were analyzed with (**Figs 2abcdefi, 3abc; S5, S6, S7, S8**)

#### Machine learning applied to identify the most important thermal parameters
* `machine_learning_on_particles.py`


#### Study the effect of different processes individually on yeast cell growth rate
* `split_factors_population.py`


#### Simulate the metabolic shift in chemostat cultivation
* `simulate_chostat_metabolic_shift.py`
* `simulate_chemostat_double_protein_limt.ipynb`
* `simulate_chemostat_stable_mitochodria.ipynb`

#### Identify most rate-limiting enzymes at 42 °C based on flux control coefficients
* `fcc_population.py` Calculate the FCC for all enzymes at 42 °C
* `resucue_ERG1.py` Simulate the specific growth rate with Posterior models with/without a temperature-insensitive ERG1
* `simulate_down_regulation_of_ERG_genes.ipynb` Simulate the down-regulation of ERG pathway
* `remove_temperature_constraints_sequentially_based_on_ind_fcc_BS.py` At each step, calculate FCC for all enzymes and identify the one with the highest FCC, remove the temperature constraint for that enzyme. Then repeat until temperature contraint has been removed for all enzymes.


#### Visualization
* `visualization.ipynb` : **Figs 2abcdefi, 3ab; 5abde; S5, S6, S7, S10, S11a**
* `additional_plots.ipynb`: - additional plots (**Figs 2gh, S13**)
* `visualize_cv.ipynb`: **Fig S3, S4**
* `visualize_temperature_on_enzymes_posterior.ipynb`: **Fig 3cdef**
* `visualize_chemostat_metabolic_shift.ipynb`: **Fig 4abc, S9**
* `expdata.ipynb`: **Fig 5c, S11b**
* `single_enzyme.ipynb`: **Fig S8**