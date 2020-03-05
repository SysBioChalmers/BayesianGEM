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
The results can be visualized with `visualize_cv.ipynb`. (**Fig S3**)  
* `gem_smcabc_at_three_conditions.py` - SMC-ABC approach update parameter space with all three observed datasets.

Above scripts need be run on high-performance cluster, and may take a few days.

* `visualization.ipynb` - the analysis of resulted Posterior models were analyzed with  (**Fig 2bcdef, 3bc; Fig S4, S5, S6**)

#### Machine learning applied to identify the most important thermal parameters
* `machine_learning_on_particles.py`
The results can be visualized with `visualize_cv.ipynb`. (**Fig 2g**)

#### Study the effect of different processes individually on yeast cell growth rate
* `split_factors_population.py`
The results can be visualized with `visualize_cv.ipynb`. (**Fig 3a**)

#### Visualize the effect of temperature on individual enzymes in Posterior models.
* `single_enzyme.ipynb`
The script is for reproducing results in (**Fig S7**)

* `visualize_temperature_on_enzymes_posterior.ipynb`
The script is for reproducing results in (**Fig 3def**)

#### Simulate the metabolic shift in chemostat cultivation
* `simulate_chostat_metabolic_shift.py`
The results can be visualized with `visualize_chemostat_metabolic_shift.ipynb`. (**Fig 4**)

#### Identify most rate-limiting enzymes at 42 °C based on flux control coefficients
* `fcc_population.py` Calculate the FCC for all enzymes at 42 °C, the results are used for (**Fig 5a**)
* `resucue_ERG1.py` Simulate the specific growth rate with Posterior models with/without a temperature-insensitive ERG1 (**Fig 5b**)
* `remove_temperature_constraints_sequentially_based_on_ind_fcc_BS.py` At each step, calculate FCC for all enzymes and identify the one with the highest FCC, remove the temperature constraint for that enzyme. Then repeat until temperature contraint has been removed for all enzymes. (**Fig 5cd**)

Above results can be visualized with `visualize_cv.ipynb`
