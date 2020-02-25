#!/usr/bin/env python
# coding: utf-8

# #### Objectives: Use a random forest regressor to regress R2 score on thermal parameters
# 1. figure out the performance of the model
# 2. use a three fold cross validation to calculate the feature importance

# In[6]:


import pickle
import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestRegressor
#from sklearn.model_selection import cross_val_score
#from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import r2_score
from scipy.stats import spearmanr,pearsonr
from scipy.stats import boxcox
#import matplotlib.pyplot as plt


# #### 1. prepare dataset
# * Column names are in the format of `P32476_Tm`, `P32476_Topt`,`P32476_dCpt`  
# * Indexes are numbers
# * The last column is the $R^2$ score of each particle  

# In[2]:


results = pickle.load(open('../results/smcabc_gem_three_conditions_save_all_particles.pkl','rb'))


# In[3]:


columns = list(results.all_particles[0].keys())
columns.sort()

data = list()
for p in results.all_particles:
    data.append([p[k] for k in columns])
df = pd.DataFrame(data=data,columns=columns)
df['r2'] = results.all_distances
print(df.shape)
print(df.columns)
print(df.index)
df.head(n=10)


# In[4]:


# Remove samples with a R2 score smaller than -3
df['r2'] = -df['r2']
sel_index = df.index[df['r2']>-3]
df = df.loc[sel_index,:]
print(df.shape)
print(df.columns)
print(df.index)
df.head(n=10)


# In[13]:


## apply boxcox transformation on r2
# df['r2'] = boxcox(df['r2']+4)[0]


# import matplotlib.pyplot as plt

# plt.hist(df['r2'],50)
# plt.show(9)

# #### 2. Split dataset into train-validation-test 80-10-10

# In[12]:


random.seed(1)
index = list(df.index)
random.shuffle(index)
df_sh = df.loc[index,:]
df_sh.head(n=10)


# X, y = df_sh.values[:,:-1], df_sh.values[:,-1].ravel()
# s = X.shape[0]
# X_train, y_train = X[:int(0.8*s),:], y[:int(0.8*s)]
# X_val,   y_val   = X[int(0.8*s):int(0.9*s),:], y[int(0.8*s):int(0.9*s)]
# X_test,  y_test  = X[int(0.9*s):,:], y[int(0.9*s):]
# print('train:',      X_train.shape, y_train.shape)
# print('validation:', X_val.shape, y_val.shape)
# print('test:',       X_test.shape, y_test.shape)

# In[ ]:


X, y = df_sh.values[:,:-1], df_sh.values[:,-1].ravel()
s = X.shape[0]
X_train, y_train = X[:int(0.8*s),:], y[:int(0.8*s)]
#X_val,   y_val   = X[int(0.8*s):int(0.9*s),:], y[int(0.8*s):int(0.9*s)]
X_test,  y_test  = X[int(0.8*s):,:], y[int(0.8*s):]
print('train:',      X_train.shape, y_train.shape)
#print('validation:', X_val.shape, y_val.shape)
print('test:',       X_test.shape, y_test.shape)


# #### 2. Optimize hyper-parameters 

# In[ ]:


report=open('../results/machine_learning_report.txt','a+',buffering=1)


# val_scores = list()
# mfs = np.arange(0.1,0.6,0.1)
# for mf in mfs:
#     rf = RandomForestRegressor(n_estimators=1000,n_jobs=-1,max_features=mf)
#     rf.fit(X_train,y_train)
#     train_score = rf.score(X_train,y_train)
#     val_score   = rf.score(X_val,y_val)
#     test_score  = rf.score(X_test,y_test)
#     print('max_features    :',mf)
#     print('Hyper-Opt, train:',train_score)
#     print('Hyper-Opt, val  :',val_score)
#     print('Hyper-Opt, test :',test_score)
#     print(' ')
#     
#     report.write('max_features    : ' + str(mf) + '\n')
#     report.write('Hyper-Opt, train: ' + str(train_score) + '\n')
#     report.write('Hyper-Opt, val  : ' + str(val_score) + '\n')
#     report.write('Hyper-Opt, test : ' + str(test_score) + '\n\n')
#     
#     val_scores.append(val_score)

# In[ ]:





# best_mf = mfs[np.argsort(val_scores)][-1]
# print('Best max_features',best_mf)
# report.write('Best max_features' + str(best_mf) + '\n\n')

# In[ ]:





# #### 3. Train final model, save feature importance

# In[ ]:





# In[ ]:


# create model with optimized params
#best_mf = 0.5
model = RandomForestRegressor(n_estimators=1000,n_jobs=-1)
model.fit(X_train,y_train)


# In[ ]:


features = df.columns[:-1]


# In[ ]:


# Model stats
p = model.predict(X)
rmse_cv = np.sqrt(MSE(p,y))
r2_cv = r2_score(y,p)
r_spearman = spearmanr(p,y)
r_pearson = pearsonr(p,y)
res = 'rmse:{:.4}\nr2:{:.4}\nspearmanr:{:.4}\np_value:{:.4}\npearonr:{:.4}\np_pearsonr:{:.4}'.format(rmse_cv,r2_cv,r_spearman[0],r_spearman[1],r_pearson[0],r_pearson[1])
print(res)
report.write(res+'\n')
report.write('test score:' + str(model.score(X_test,y_test)))


# In[ ]:


# feature importance
fhand = open('../results/machine_learning_feature_importance.csv','w')
fhand.write('feature,importance,std\n')

std = np.std([tree.feature_importances_ for tree in model.estimators_],axis=0)
for i in range(X.shape[1]):
    fhand.write('{0},{1},{2}\n'.format(features[i],model.feature_importances_[i],std[i]))
fhand.close()


# In[ ]:




