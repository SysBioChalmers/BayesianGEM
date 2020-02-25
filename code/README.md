### Description
This folder contains all scripts to carry out all analysis reported in the paper. 

#### Fundmental scripts
* `etcpy` contains the script and introduction `etcpy/README.md` of how to include temperature as an addtional variable in the enzyme constrained yeast 7.6 model. 
* `GEMS.py` contains the functions for simulating the aerobic, anaerobic growth rate in batch cultivation and aerobic fluxes in chemostat clutivation. A list of distance functions used for SMC-ABC approach are also included. 
* `abc_etc.py` contains the function to perform SMC-ABC approach.

#### Simulations with start parameters (without considering undertainties).
* `simulate_with_start_params.ipynb`
The script also contains the function for visualization (**Fig S2**). The script can be done on PC.

#### SMC-ABC approach updates model parameters
Three experimental datasets will be used in this section: 
At different temperatures  
(1) the maximal specific growth rate in aerobic batch cultivations (Caspeta L., et al. Mbio, 2015)
(2) anaerobic  batch cultivations (Zakhartsev M., et al. J. Therm. Biol., 2015)
(3) fluxes of carbon dioxide (CO2), ethanol and glucose in chemostat cultivations (Postmus J., J. Biol. Chem., 2008)  
   
Cross-validation:
* `gem_smcabc_at_three_conditions_cv1.py`
* `gem_smcabc_at_three_conditions_cv2.py`
* `gem_smcabc_at_three_conditions_cv3.py`
The results can be visualized with `visualize_cv.ipynb`. (**Fig S3**)  
  
SMC-ABC approach update parameter space with all three observed datasets.
* `gem_smcabc_at_three_conditions.py`

Above scripts have to be performed with a computer cluster. It may take a few days.  
The analysis of resulted Posterior models were analyzed with `visualization.ipynb` (**Fig 2bcdef, 3bc; Fig S4, S5, S6**)

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

All above results can be visualized with `visualize_cv.ipynb`





