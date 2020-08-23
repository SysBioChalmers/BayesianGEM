#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import pickle
import GEMS
import os
from multiprocessing import Process,cpu_count,Manager
from etcpy import etc


# In[3]:


def chemostat_double_protein_limt(thermalParams,Ts):
    df,new_params = GEMS.format_input(thermalParams)
    
    mae = pickle.load(open('../models/aerobic.pkl','rb'))
    mae.reactions.get_by_id('prot_pool_exchange').upper_bound *=2
    print(mae.reactions.get_by_id('prot_pool_exchange').upper_bound)
    growth_id = 'r_2111'
    glc_up_id = 'r_1714_REV'
    prot_pool_id = 'prot_pool_exchange'
    dilut = 0.1
    sigma = 1.0

    met_names = ['Glucose','CO2','Ethanol']
    
    try:
        solution = etc.simulate_chomostat(mae,dilut,new_params,Ts+273.15,
                                          sigma,growth_id,glc_up_id,prot_pool_id)
    
    except: solution = None

    return solution


# In[3]:


def worker(particle,index,Q,Ts):
    
    results = chemostat_double_protein_limt(particle,Ts)
    
    Q.put((index,results))


# In[4]:


outfile = '../results/chemostat_double_protein_limt.pkl'


# In[ ]:


particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population
Ts = np.arange(30,45,0.5)
Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle,index,Q,Ts)) 
                               for index,particle in enumerate(particles)]

for p in jobs: p.start()
for p in jobs: p.join()


# In[ ]:


results_population = [None for _ in particles] 
for index,res in [Q.get(timeout=1) for p in jobs]: results_population[index] = res
pickle.dump([Ts,results_population],open(outfile,'wb'))


# In[ ]:




