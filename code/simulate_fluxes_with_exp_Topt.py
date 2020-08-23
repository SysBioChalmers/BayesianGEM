#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import pandas as pd
import pickle
import GEMS
import os
from multiprocessing import Process,cpu_count,Manager
from etcpy import etc


# In[8]:


def get_target_rxns(model):
    topt_exp = {}
    for line in open('../data/1_temperature_optimum_data_filtered_assay_at.tsv'):
        if line.startswith('ec'): continue
        cont = line.split()
        topts = [float(item) for item in cont[2].split(',')]
        topt_exp[cont[-1]] = topt_exp.get(cont[-1],[]) + topts
        
    rxns = []
    for uni in topt_exp.keys(): 
        try:
            met_prot = model.metabolites.get_by_id('prot_{0}'.format(uni))
            lst = [rxn.id for rxn in met_prot.reactions if 'prot_pool' not in rxn.reaction]
            rxns.extend(lst)
        except: continue
    return rxns,topt_exp


# In[9]:


def set_temperature_constraints(model,particle,T):
    df,_ = GEMS.format_input(particle)
    
    etc.map_fNT(model,T,df)
    etc.map_kcatT(model,T,df)
    etc.set_NGAMT(model,T)
    etc.set_sigma(model,0.5)


# In[10]:


def get_flux_range_through_given_rxn(model,rxn_id,growth_id):

    with model:
        rmax = model.optimize().objective_value
        model.reactions.get_by_id(growth_id).upper_bound = rmax*0.999
        
        model.objective = rxn_id
        model.objective.direction = 'max'
        vmax = model.optimize().objective_value
        
        model.objective.direction = 'min'
        vmin = model.optimize().objective_value

    return vmin,vmax


# In[28]:


def simulate_one_particle(particle,T,target_rxns):
    growth_id = 'r_2111'
    print('Load model')
    mae = pickle.load(open('../models/aerobic.pkl','rb'))
    
    print('Set temperature constraits')
    set_temperature_constraints(mae,particle,T)
    
    print('Get vmax, vmin')
    # This is to avoid some unbounded results
    for rxn_id in target_rxns:
        mae.reactions.get_by_id(rxn_id).upper_bound = 1000
    
    fvas = {}
    for rxn_id in target_rxns:
        try:
            vmin,vmax = get_flux_range_through_given_rxn(mae,rxn_id,growth_id)
            fvas[rxn_id] = [vmin,vmax]
        except:
            fvas[rxn_id] = [None,None]
    return fvas


# In[23]:


def worker(particle,index,Q,target_rxns):
    results = {32:None,40:None}
    #try:
    results[32] = simulate_one_particle(particle,32+273.15,target_rxns)
    results[40] = simulate_one_particle(particle,40+273.15,target_rxns)
    #except: None
    Q.put((index,results))


# In[87]:


outfile = '../results/32_vs_40_topt_exp_rxns.pkl'
particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population

mae = pickle.load(open('../models/aerobic.pkl','rb'))
target_rxns,topt_exp = get_target_rxns(mae)

Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle,index,Q,target_rxns)) 
                               for index,particle in enumerate(particles)]

for p in jobs: p.start()
for p in jobs: p.join()


# In[ ]:


results_population = [None for _ in particles] 
for index,res in [Q.get(timeout=1) for p in jobs]: results_population[index] = res
pickle.dump([topt_exp,results_population],open(outfile,'wb'))

