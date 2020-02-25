# etcpy
#### Important: this is now only for ecYeast7 only since the protein usage, NGAMT are different in other models. Will be adapted in the future if necessary.
This script provides several functions to map temperature effects on enzymes in
an enzyme constrained model. In the model, three temperature effects are considered:
* protein denaturation, which is described with a two-state denaturation model  
* enzyme kcat, which is described with a macromolecular theory macromolecular rate theory (MMRT): Hobbs J. et al., ACS Chem Biol. 2013  
* temperature dependent non-growth associated ATP maintenance (NGAM), which is from the experimentally measured NGAM at different temperatures under anaerobic conditions. The same dependence was assumed for aerobic growth  

## Installation
##### (1). Download tome package
##### (2). Open your terminal
##### (3). Change directory to the package directory
```linux
cd [directory to etcpy, where setup.py is]
```
##### (4). Run following command
```linux
pip install -e .
```

## Usage:  
The input files:
* an enzyme constrained model
* a dataframe containing the following parameters of all enzymes in the model
  * dHTH: the enthalpy change at convergence temperature TH (373.5 K), in J/mol
  * dSTS: the entropy change at convergence temperature TS (385 K), in J/mol/K
  * dCpu: the heat-capacity change of protein unfolding process, in J/mol/K
  * dCpt: the heat-capacity change between transition state and ground state in catalytic process, in J/mol/K
  * Topt: the optimal catalytic temperature of an enzymes, in K. At this temperature, the specific activity is maximal

```python
from etcpy import etc
```
1. temperature effect on protein denaturation  
```python
etc.map_fNT(model,T,df)
```
model is the cobra model object. T is temperature in K. df, the dataframe containing the following parameters of all enzymes in the model  

2. temperature effect on kcat  
```python
etc.map_kcatT(model,T,df)
```
model is the cobra model object. T is temperature in K. df, the dataframe containing the following parameters of all enzymes in the model  

3. temperature effect on NGAM
```python
etc.set_NGAMT(model,T)
```

4. set sigma
```python
etc.set_sigma(model,sigma)
```

5. simulate growth at different temperatures.
```python
etc.simulate_growth(model,Ts,sigma,df)
```
Ts is a list of tempertures in K.

6. calculate dHTH, dSTS and dCpu of the denaturation process. There are two functions provided for two difference scenarios.
```python 
dHTH,dSTS,dCpu = etc.get_dH_dS_dCpu_from_TmT90(Tm,T90)  #for a protein with experimental Tm and T90
dHTH,dSTS,dCpu = etc.get_dH_dS_dCpu_from_TmLength(Tm,proteinLength) # for protein with only Tm
```
Make sure that dCpu obtained from `etc.get_dH_dS_dCpu_from_TmT90(Tm,T90)` must be positive. If not, use `etc.get_dH_dS_dCpu_from_TmLength(Tm,proteinLength)`. 

7. sample the uncertainties in the thermal parameters Tm,Topt and dCpt. Given a dataframe containing following columns:Tm,Tm_std,T90,dCpt,dCpt_std,Topt,Topt_std,Length. Randomly generate a new value from a normal distribution N(Tm,Tm_std) taking Tm as an example. 
```python 
new_params = etc.sample_data_uncertainty(params) # One can also specify columns to be sampled. The default is to sample all columns: [Tm,dCpt,Topt]
thermalparams = etc.calculate_thermal_params(new_params)
```
Then ```thermalparams``` could be used to simulate growth rate ```etc.simulate_growth(model,Ts,sigma,thermalparams)```


8. simulate chemostat data.  
(1) fix growth rate, set objective function as minimizing glucose uptatke rate  
(2) map temperature parameters  
(3) get minimal glucose uptatke rate, then fix glucose uptatke rate to this minimal value (\*1.001 for simulation purpose)  
(4) minimize enzyme usage
(5) return solution = model.optimize()
```python 
%%time
params = pd.read_csv('./model_enzyme_params.csv',index_col=0) # contains at least Tm,T90,dCpt,Topt,Length
Ts = np.array([30,40,42,43,44,45,50])+273.15
growth_id = 'r_2111'
glc_up_id = 'r_1714_REV'
prot_pool_id = 'prot_pool_exchange'
solutions = etc.simulate_chomostat(model,0.1,params,Ts,0.5,growth_id,glc_up_id,prot_pool_id)
```
  
9. calculate flux control coefficients
```python
dfres = etc.do_fcc_at_Ts(model,Ts,params,sigma=0.5,delta=10)
```
Obtained `dfres` is a dataframe with enzymes as index and temperatures as columns. Values are flux control coefficients.


