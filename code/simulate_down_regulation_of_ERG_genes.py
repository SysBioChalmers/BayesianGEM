#!/usr/bin/env python
# coding: utf-8

# In[104]:


import numpy as np
import pandas as pd
import pickle
import GEMS
import os
from multiprocessing import Process,cpu_count,Manager
from etcpy import etc


# In[105]:


def load_gene_names():
    gene_names = dict()
    gene_ids = dict()
    for line in open('../data/enzyme_uniprot_gene_name.csv'):
        cont = line.strip().split(',')
        gene_names[cont[0]] = cont[1]
        gene_ids[cont[1]] = cont[0]
    return gene_names,gene_ids


# In[17]:





# In[4]:


gene_names,gene_ids = load_gene_names()


# In[19]:





# In[107]:


def get_erg_rxns(model):
    ERG_genes = ['ERG10','ERG13','HMG1','HMG2','ERG12','ERG8','MVD1','IDI1','ERG20',
             'ERG9','ERG7','ERG11','ERG24','ERG25','ERG26','ERG27','ERG6','ERG2',
             'ERG3','ERG5','ERG4'] # MVD1 is ERG19
    rxns = []
    for g in ERG_genes: 
        uni = gene_ids[g]
        met_prot = model.metabolites.get_by_id('prot_{0}'.format(uni))
        lst = [rxn.id for rxn in met_prot.reactions if 'prot_pool' not in rxn.reaction]
        rxns.extend(lst)
    return rxns


# In[108]:


def set_temperature_constraints(model,particle,T):
    df,_ = GEMS.format_input(particle)
    
    etc.map_fNT(model,T,df)
    etc.map_kcatT(model,T,df)
    etc.set_NGAMT(model,T)
    etc.set_sigma(model,0.5)


# In[109]:


def get_maximal_flux_through_given_rxn(model,rxn_id,growth_id):

    with model:
        rmax = model.optimize().objective_value
        model.reactions.get_by_id(growth_id).upper_bound = rmax*0.999
        
        model.objective = rxn_id
        model.objective.direction = 'max'
        vmax = model.optimize().objective_value

    return vmax


# In[110]:


def down_regulate_erg_genes(model,erg_vmax,p_down):
    # p_down, percentage of down regulation
    # erg_vmax, {rxn_id, vmax}
    for rxn_id, vmax in erg_vmax.items():
        model.reactions.get_by_id(rxn_id).upper_bound = vmax*p_down


# In[111]:


def simulate_one_particle(particle,T,ps):
    growth_id = 'r_2111'
    print('Load model')
    mae = pickle.load(open('../models/aerobic.pkl','rb'))
    
    # get ori coeffs for ERG1
    met1 = mae.metabolites.get_by_id('prot_{0}'.format('P32476'))
    old_kcat_coeffs = {rxn.id: rxn.get_coefficient(met1) for rxn in met1.reactions}
            
    met2 = mae.metabolites.prot_pool
    rxn = mae.reactions.get_by_id('draw_prot_{0}'.format('P32476'))
    old_fnt_coeff = rxn.get_coefficient(met2)
    
    print('Set temperature constraits')
    set_temperature_constraints(mae,particle,T)
    erg_rxns = get_erg_rxns(mae)
    
    print('Get vmax')
    # This is to avoid some unbounded results
    for rxn_id in erg_rxns:
        mae.reactions.get_by_id(rxn_id).upper_bound = 1000
    
    erg_vmax = {rxn_id:get_maximal_flux_through_given_rxn(mae,rxn_id,growth_id) for rxn_id in erg_rxns}
    
    rmax_vary_p = []
    rmax_vary_p_rescue = []
    
    print('Simulate down regulation')
    for p in ps:
        with mae:
            down_regulate_erg_genes(mae,erg_vmax,p)
            rmax_vary_p.append(mae.optimize().objective_value)
            
            # rescue prot
            for rxn_id,old_coeff in old_kcat_coeffs.items():
                rxn = mae.reactions.get_by_id(rxn_id)
                etc.change_rxn_coeff(rxn,met1,old_coeff)
            
            rxn = mae.reactions.get_by_id('draw_prot_{0}'.format('P32476'))
            etc.change_rxn_coeff(rxn,met2,old_fnt_coeff)
            rmax_vary_p_rescue.append(mae.optimize().objective_value)
    return rmax_vary_p, rmax_vary_p_rescue


# In[112]:


def worker(particle,index,Q,T,ps):
    
    try:results = simulate_one_particle(particle,T,ps)
    except: results = None
    Q.put((index,results))


# In[87]:


ps = np.arange(0.01,1.01,0.01)
T = 40 + 273.15
outfile = '../results/down_regulate_erg_pathway_{0}.pkl'.format(int(T-273.15))
particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population

Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle,index,Q,T,ps)) 
                               for index,particle in enumerate(particles)]

for p in jobs: p.start()
for p in jobs: p.join()


# In[ ]:


results_population = [None for _ in particles] 
for index,res in [Q.get(timeout=1) for p in jobs]: results_population[index] = res
pickle.dump([T,ps,results_population],open(outfile,'wb'))


# In[114]:


T = 40 + 273.15
Q = Manager().Queue()
particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population
worker(particles[0],0,Q,T,ps)


# In[115]:


s = Q.get()


# In[116]:




