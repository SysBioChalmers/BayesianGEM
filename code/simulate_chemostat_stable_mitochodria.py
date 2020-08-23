#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import pickle
import GEMS
import os
from multiprocessing import Process,cpu_count,Manager
from etcpy import etc


# In[5]:


def load_gene_names():
    gene_names = dict()
    gene_ids = dict()
    for line in open('../data/enzyme_uniprot_gene_name.csv'):
        cont = line.strip().split(',')
        gene_names[cont[0]] = cont[1]
        gene_ids[cont[1]] = cont[0]
    return gene_names,gene_ids


# In[9]:


_,gene_ids = load_gene_names()


# In[2]:


def chemostat_double_protein_limt(thermalParams,Ts):
    mae = pickle.load(open('../models/aerobic.pkl','rb'))
    
    # stablize ATP1, HEM1 and PDB1, set Tm as 50 C
    for name in ['ATP1', 'HEM1', 'PDB1']:
        thermalParams['{0}_Tm'.format(gene_ids[name])] = 50 + 273.15
    
    df,new_params = GEMS.format_input(thermalParams)
    
    print( mae.reactions.get_by_id('prot_pool_exchange').upper_bound)
    growth_id = 'r_2111'
    glc_up_id = 'r_1714_REV'
    prot_pool_id = 'prot_pool_exchange'
    dilut = 0.1
    sigma = 0.5
    
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


outfile = '../results/chemostat_stable_mito.pkl'


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

