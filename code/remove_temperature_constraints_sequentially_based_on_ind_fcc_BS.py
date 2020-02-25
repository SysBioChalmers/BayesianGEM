#!/usr/bin/env python
# coding: utf-8

# At 42 °C, 
# * Perform FCC at each step
# * Remove the temperature constraints one by one from the one with the biggest FCC score. 
# * Simulate the maximal growth rate.

# In[16]:


import numpy as np
import pandas as pd
import pickle
import os
from etcpy import etc
from multiprocessing import Process,cpu_count,Manager


# In[17]:


def sort_proteins_on_fcc(dfmean,T):
    proteins = dfmean.index[np.argsort(-dfmean[T])] # from biggest to smallest
    return proteins


# In[ ]:


def do_a_fcc(model,tested_enzymes,enzymes):
    #  tsted_enzymes: a list of enzyme already test. In the order of test iteration
    
    # 2. do fcc at each temperature
    fccs = dict()
    delta = 10
    # Get all enzymes
    try:    u0 = model.optimize().objective_value
    except: u0 = np.nan
    print('Model Track u0:',T,u0)
    for enz_id in enzymes: 
        if enz_id in tested_enzymes: continue
        if np.isnan(u0) or u0 == 0: fccs[enz_id] = np.nan
        else:
            with model as m:
                enz = m.metabolites.get_by_id(enz_id)
                # Perturbe kcat in all reactions of this enzyme 
                for rxn in enz.reactions:
                    if rxn.id.startswith('draw_'): continue
                    new_coeff = rxn.get_coefficient(enz)/(1+delta)
                    etc.change_rxn_coeff(rxn,enz,new_coeff)
                try:u = m.optimize().objective_value
                except: u = np.nan
                
                if u is None: fc = np.nan
                else: fc = (u-u0)/u0/delta

                fccs[enz.id] = fc
    return fccs


# In[1]:


def find_most_rate_limiting_enzyme(fccs):
    # return 
    lst = [(v,k) for k,v in fccs.items() if ~np.isnan(v)]
    lst.sort()
    try: return lst[-1][1]
    except: return None


# In[3]:


def remove_temperature_constraints(ecmodel,etcmodel,prot):
    # ecmodel: original model without temperature constraints
    # etcmodel: model which have temperature constraints for some/all enzymes
    # prot: uniprot id
    
    # original coeffes
    met1 = ecmodel.metabolites.get_by_id('prot_{0}'.format(prot))
    old_kcat_coeffs = {rxn.id: rxn.get_coefficient(met1) for rxn in met1.reactions}
            
    met2 = ecmodel.metabolites.prot_pool
    rxn = ecmodel.reactions.get_by_id('draw_prot_{0}'.format(prot))
    old_fnt_coeff = rxn.get_coefficient(met2)
    
    
    # rescue prot
    met1 = etcmodel.metabolites.get_by_id('prot_{0}'.format(prot))
    for rxn_id,old_coeff in old_kcat_coeffs.items():
        rxn = etcmodel.reactions.get_by_id(rxn_id)
        etc.change_rxn_coeff(rxn,met1,old_coeff)
            
    met2 = etcmodel.metabolites.prot_pool
    rxn = etcmodel.reactions.get_by_id('draw_prot_{0}'.format(prot))
    etc.change_rxn_coeff(rxn,met2,old_fnt_coeff)


# In[4]:


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


# In[5]:


def simulate_a_particle(particle,T):
    df,_ = format_input(particle)
    
    sigma = 0.5
    # map all temperature constrains on all enzymes
    ecmodel = pickle.load(open('../models/aerobic.pkl','rb'))
    etcmodel = pickle.load(open('../models/aerobic.pkl','rb'))
    
    # map temperature constraints
    etc.map_fNT(etcmodel,T,df)
    etc.map_kcatT(etcmodel,T,df)
    etc.set_NGAMT(etcmodel,T)
    etc.set_sigma(etcmodel,sigma)
    
    # get enzyme list
    enzymes = list(set([met.id for met in ecmodel.metabolites 
                        if met.id.startswith('prot_') and met.id != 'prot_pool']))
    
    tested_enzymes = []
    
    
    try: r = etcmodel.optimize().objective_value
    except: r = 0
    rs = [r]
    
    for i in range(len(enzymes)-1):
        fccs = do_a_fcc(etcmodel,tested_enzymes,enzymes)
        most_limit_enz_id = find_most_rate_limiting_enzyme(fccs)
        tested_enzymes.append(most_limit_enz_id)
        
        if most_limit_enz_id is not None: 
            remove_temperature_constraints(ecmodel,etcmodel,most_limit_enz_id.split('_')[1])
            try: r = etcmodel.optimize().objective_value
            except: r = np.nan
        else: 
            r = np.nan
            break
        print(most_limit_enz_id,r)
        rs.append(r)

    return rs, tested_enzymes


# In[6]:


def worker(particle,index,Q,T,dffcc):
    try:
        rs, tested_enzymes = simulate_a_particle(particle,T)
    except:
        rs = [None]+[None for _ in dffcc.index]
        tested_enzymes = [None for _ in dffcc.index]
    
    Q.put((index,rs,tested_enzymes))


# In[7]:


outfile = '../results/sequentially_remove_temperature_constraints_ind_BS_all_particles.pkl'


# In[ ]:


particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population
dffccs = pickle.load(open('../results/fcc_population_all_particles.pkl','rb')) # just used to get the number of enzymes
T = 42 + 273.15
Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle,index,Q,T,dffccs[index])) 
                               for index,particle in enumerate(particles)]

for p in jobs: p.start()
for p in jobs: p.join()


# In[ ]:


results_population = [None for _ in particles] 
tested_enzymes_population = [None for _ in particles] 
for index,res,tested_enzymes in [Q.get(timeout=1) for p in jobs]: 
    results_population[index] = res
    tested_enzymes_population[index] = tested_enzymes
    
pickle.dump([results_population,tested_enzymes_population],open(outfile,'wb'))

