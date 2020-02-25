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


outfile = '../results/fcc_population_all_particles.pkl'


# In[3]:


def do_fcc_at_T(T,params,sigma,delta):
    # Ts:     temperatures at which simulation will be performed
    # params: a dataframe with Topt, Tm, T90, dCpt of all enzymes
    # sigma:  enzyme saturation factor
    # delta:  the fold change to be made for kcat of each enzyme
    #
    # return a dataframe with temperature as column and enzyme id as index
    
    model = pickle.load(open('../models/aerobic.pkl','rb'))
    
    # 1. Calculate thermal parameters
    dfparam = etc.calculate_thermal_params(params)
    enzymes = list(set([met.id for met in model.metabolites 
                        if met.id.startswith('prot_') and met.id != 'prot_pool']))
    
    # 2. do fcc at each temperature
    data = dict()
    
    etc.map_fNT(model,T,dfparam)
    etc.map_kcatT(model,T,dfparam)
    etc.set_NGAMT(model,T)
    etc.set_sigma(model,sigma)

    # Get all enzymes
    try:    u0 = model.optimize().objective_value
    except: u0 = None
    print('Model Track u0:',T,u0)
    for enz_id in enzymes: 
        if u0 is None or u0 == 0: data[enz_id.split('_')[1]] = None
        else:
            with model as m:
                enz = m.metabolites.get_by_id(enz_id)
                # Perturbe kcat in all reactions of this enzyme 
                for rxn in enz.reactions:
                    if rxn.id.startswith('draw_'): continue
                    new_coeff = rxn.get_coefficient(enz)/(1+delta)
                    etc.change_rxn_coeff(rxn,enz,new_coeff)
                try:u = m.optimize().objective_value
                except: u = None
                
                if u is None: fc = None
                else: fc = (u-u0)/u0/delta

                data[enz.id.split('_')[1]] = fc

    # 3. Create a dataframe with temperature as column and enzyme ids as index
    lst, index = [], []
    for k,v in data.items():
        lst.append(v)
        index.append(k)
    dfres = pd.DataFrame(data={T:lst},index=index)
    return dfres        


# In[4]:


def do_fcc_at_Ts(Ts,params,sigma=0.5,delta=10):
    # Ts:     a list of temperatures at which simulation will be performed
    # params: a dataframe with Topt, Tm, T90, dCpt of all enzymes
    # sigma:  enzyme saturation factor
    # delta:  the fold change to be made for kcat of each enzyme
    #
    # return a dataframe with temperatures as columns and enzyme ids as index
    
    results = None
    for T in Ts:
        dfres = do_fcc_at_T(T,params,sigma,delta)
        if results is None: results = dfres
        else: results = pd.merge(results,dfres,left_index=True,right_index=True,how='inner')
                
    print('FCC at Ts:',results.shape)
    return results


# In[5]:


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


# In[6]:


def worker(particle,index,Q,Ts):
    df,new_params = format_input(particle)
    
    results = do_fcc_at_Ts(Ts,new_params,sigma=0.5,delta=10)
    
    Q.put((index,results))


# In[7]:


import time


# In[ ]:


particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population
Ts = np.array([5, 10, 15, 20, 25, 30, 35, 37, 40, 42, 45]) + 273.15
Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle,index,Q,Ts)) 
                               for index,particle in enumerate(particles)]

for p in jobs: p.start()
for p in jobs: p.join()


# In[ ]:


results_population = [None for _ in particles] 
for index,res in [Q.get(timeout=1) for p in jobs]: results_population[index] = res
pickle.dump(results_population,open(outfile,'wb'))

