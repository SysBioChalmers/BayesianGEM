#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import pickle
import os
from etcpy import etc
from multiprocessing import Process,cpu_count,Manager


# In[2]:


def format_input(thermalParams):
    # thermalParams: a dictionary with ids like uniprotid_Topt 
    params = pd.read_csv('../data/model_enzyme_params.csv',index_col=0)
    new_params = params.copy()
    for key,val in thermalParams.items():
        [ind,col] = key.split('_')
        new_params.loc[ind,col] = val
    
    # Update T90
    new_params['T90'] = params['T90']-params['Tm'] + new_params['Tm']

    df = etc.calculate_thermal_params(new_params)
    return df,new_params


# In[17]:


def simulate_growth(model,Ts,sigma,df,prot):
    '''
    # model, cobra model
    # Ts, a list of temperatures in K
    # sigma, enzyme saturation factor
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ensure that Topt is in K. Other parameters are in standard units.
    # prot: uniprotid, enzyme to be rescued
    '''
    rs0 = list()
    rs = list()
    
    met1 = model.metabolites.get_by_id('prot_{0}'.format(prot))
    old_kcat_coeffs = {rxn.id: rxn.get_coefficient(met1) for rxn in met1.reactions}
            
    met2 = model.metabolites.prot_pool
    rxn = model.reactions.get_by_id('draw_prot_{0}'.format(prot))
    old_fnt_coeff = rxn.get_coefficient(met2)
    
    for T in Ts:
        with model:
            
            # map temperature constraints
            etc.map_fNT(model,T,df)
            etc.map_kcatT(model,T,df)
            etc.set_NGAMT(model,T)
            etc.set_sigma(model,sigma)
            
            try: r = model.optimize().objective_value
            except:
                print('Failed to solve the problem')
                r = 0
            print(T-273.15,r)
            rs0.append(r)
            
            
            # rescue prot
            met1 = model.metabolites.get_by_id('prot_{0}'.format(prot))
            for rxn_id,old_coeff in old_kcat_coeffs.items():
                rxn = model.reactions.get_by_id(rxn_id)
                etc.change_rxn_coeff(rxn,met1,old_coeff)
            
            met2 = model.metabolites.prot_pool
            rxn = model.reactions.get_by_id('draw_prot_{0}'.format(prot))
            etc.change_rxn_coeff(rxn,met2,old_fnt_coeff)

            try: r = model.optimize().objective_value
            except:
                print('Failed to solve the problem')
                r = 0
            print(T-273.15,r)
            rs.append(r)
    return rs0,rs


# In[12]:


def worker(particle,index,Q,Ts):
    model = pickle.load(open('../models/aerobic.pkl','rb'))
    prot = 'P32476'
    df, _ = format_input(particle)
    rs = simulate_growth(model,Ts,0.5,df,prot)
    
    Q.put((index,rs))


# In[ ]:


outfile = '../results/rescue_ERG1_all_particles.pkl'


# In[ ]:


particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population
Ts = np.array([30, 35, 37, 38,39,40,41,42,45]) + 273.15
Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle,index,Q,Ts)) 
                               for index,particle in enumerate(particles)]

for p in jobs: p.start()
for p in jobs: p.join()


# In[ ]:


results_population = [None for _ in particles] 
for index,res in [Q.get(timeout=1) for p in jobs]: results_population[index] = res
pickle.dump(results_population,open(outfile,'wb'))


# In[ ]:




