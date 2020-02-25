#!/usr/bin/env python
# coding: utf-8

# In[1]:


import abc_etc as abc
import numpy as np
import GEMS
import os
import pandas as pd
import pickle


# In[2]:


Yobs_batch = GEMS.aerobic_exp_data()
Yobs_chemo = GEMS.chemostat_exp_data()
#Yobs_batch_an = GEMS.anaerobic_exp_data()
dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}


# In[5]:


Yobs = {'rae':Yobs_batch['data'],
        'chemostat':Yobs_chemo['data'],
        'ran':Yobs_batch_an['data']}


# In[ ]:


path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)


# #### Define priors

# In[ ]:


priors = dict()
for ind in params.index: 
    for col in ['Tm','Topt','dCpt']: 
        priors['{0}_{1}'.format(ind,col)] = abc.RV('normal',
                                                      loc=params.loc[ind,col],
                                                      scale=params.loc[ind,col+'_std'])


# #### Define model settings

# In[ ]:


min_epsilon = -1 # equivalent to r2 score of 1
population_size = 100
outfile = '../results/smcabc_gem_three_conditions_save_all_particles.pkl'


# In[ ]:


if not os.path.exists(outfile):
    print('Initialize model')
    model = abc.SMCABC(GEMS.simulate_at_three_conditions_2,
                        priors,
                        min_epsilon,
                        population_size,
                        GEMS.distance_2,
                        Yobs,
                        outfile,
                        generation_size=128)
else: model = pickle.load(open(outfile,'rb'))


# #### Run simulation

# In[ ]:


print('start simulations')
model.run_simulation()

