#!/usr/bin/env python
# coding: utf-8

# In[1]:


# do chemostat simulation at dilution rate of 0.1, at various temperatures


# In[2]:


import numpy as np
import pandas as pd
import pickle
import os
from etcpy import etc
from multiprocessing import Process,cpu_count,Manager


# In[ ]:


path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)


# In[ ]:


def format_input(thermalParams):
    new_params = params.copy()
    for key,val in thermalParams.items():
        [ind,col] = key.split('_')
        new_params.loc[ind,col] = val
    
    # Update T90
    new_params['T90'] = params['T90']-params['Tm'] + new_params['Tm']

    df = etc.calculate_thermal_params(new_params)
    return df,new_params


# In[ ]:


def chemostat(thermalParams,Ts):
    df,new_params = format_input(thermalParams)
    mae = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
    
    
    growth_id = 'r_2111'
    glc_up_id = 'r_1714_REV'
    prot_pool_id = 'prot_pool_exchange'
    dilut = 0.1
    sigma = 0.5
    
    solutions = etc.simulate_chomostat(mae,dilut,new_params,Ts,
                                              sigma,growth_id,glc_up_id,prot_pool_id)

    return  solutions


# In[5]:


def worker(particle,index,Q,Ts):
    solutionTs = chemostat(particle,Ts)
    Q.put((index,solutionTs))


# In[ ]:


outfile = '../results/chemostat_solutions_all_particles.pkl'
particles = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb')).population
Ts = np.arange(30,40,0.2) + 273.15

Q = Manager().Queue()
jobs = [Process(target=worker,args=(particle,index,Q,Ts)) 
                               for index,particle in enumerate(particles)]


for p in jobs: p.start()
for p in jobs: p.join()
    
    
results_population = [None for _ in particles] 
for index,res in [Q.get(timeout=1) for p in jobs]: results_population[index] = res
pickle.dump(results_population,open(outfile,'wb'))


# In[ ]:




