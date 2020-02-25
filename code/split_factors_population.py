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


outfile = '../results/splited_factors_poppulation_all_particles.pkl'


# In[3]:


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


# In[4]:


def split_factor(particle, Ts, sigma):
    
    df,new_params = format_input(particle)
    dfparam = etc.calculate_thermal_params(new_params)
    
    results = dict() # results = {'etcYeast7':[]}
    
    # Case (1) No temperature constraints
    model = pickle.load(open('../models/aerobic.pkl','rb'))
    with model as m:
        etc.set_sigma(m,sigma)
        r = m.optimize().objective_value
        results['ori'] = [r for T in Ts]
        print('Case 1:',results['ori'])
    
    # Case (2) based on (1) Include the temperature dependent NGAM
    model = pickle.load(open('../models/aerobic.pkl','rb'))
    for T in Ts:
        with model as m:
            etc.set_NGAMT(m,T)
            etc.set_sigma(m,sigma)
            
            r = m.optimize().objective_value
            results['NGAM'] = results.get('NGAM', []) + [r]
    print('Case 2:',results['NGAM'])
    
    # Case (3) based on (1) Include the temperatuer dependent kcat
    model = pickle.load(open('../models/aerobic.pkl','rb'))
    for T in Ts:
        with model as m:
            etc.map_kcatT(m,T,dfparam)
            etc.set_sigma(m,sigma)
            
            r = m.optimize().objective_value
            results['kcat'] = results.get('kcat', []) + [r]
    print('Case 3:', results['kcat'])
    
    # Case (4) based on (1) Include the temperature dependent unfolding process
    model = pickle.load(open('../models/aerobic.pkl','rb'))
    for T in Ts:
        with model as m:
            etc.map_fNT(m,T,dfparam)
            etc.set_sigma(m,sigma)

            try: r = m.optimize().objective_value
            except: r = None
            results['unfolding'] = results.get('unfolding', []) + [r]
            
    print('Case 4:',results['unfolding'])
    
    # Case (5) with all constraints
    model = pickle.load(open('../models/aerobic.pkl','rb'))
    results['etc'] = etc.simulate_growth(model,Ts,sigma,dfparam)
    print('Case 5:',results['etc'])
    
    return results


# In[5]:


def worker(particle, Ts, sigma, index, Q):
    results = split_factor(particle, Ts, sigma)
    Q.put((index,results))


# In[6]:


particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population
Ts = np.arange(5,46,3) + 273.15
Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle, Ts, 0.5, index, Q)) 
                               for index,particle in enumerate(particles)]
        
for p in jobs: p.start()
for p in jobs: p.join()


# In[7]:


results_population = [None for _ in particles] 
for index,res in [Q.get(timeout=1) for p in jobs]: results_population[index] = res
pickle.dump(results_population,open(outfile,'wb'))


# In[ ]:




