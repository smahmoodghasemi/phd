#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import time
from time import time
import pickle
import warnings
warnings.filterwarnings("ignore")
from matplotlib import pyplot as plt


# In[2]:


df=pd.read_csv(r'C:\Users\smahm\NFSIM\NF.csv')

df.columns=np.append(np.array([-1]),np.arange(0,165,15))
df=df.reset_index()
del df['index']


# In[3]:


#probability of Activation for a single allele ==> March Experiment (Non Flavopiridol)
M2=[0.13,0.22,0.32,0.41,0.46,0.51]
#probability of Activation for a single allele ==> September Experiment (Non Flavopiridol)
S2=[0.32,0.40,0.48,0.54,0.59,0.64]
#probability of Activation for a single allele ==> November Experiment (Non Flavopiridol)
N2=[0.26, 0.30, 0.36, 0.39, 0.39, 0.41]


# In[ ]:



bi=time()
df['MTD']=""
df['LIF']=""
df['A']=""
df['C']=""
df['B']=""
df['Mar15']=""
df['Sep15']=""
df['Nov15']=""

for i in range(df.shape[0]):

    df['LIF'][i]=float((np.array(df[-1])[i].split("LIF")[1]).split("A")[0])
    df['A'][i]=float((np.array(df[-1])[i].split("A")[1]).split("C")[0])
    df['C'][i]=float((np.array(df[-1])[i].split("C")[1]).split("B")[0])
    df['B'][i]=float((np.array(df[-1])[i].split("B")[1]).split("MTD")[0])
    df['MTD'][i]=float((np.array(df[-1])[i].split("B")[1]))
    df['Mar15'][i]=np.round(np.max([np.abs(M2[5]-df[90][i]),np.abs(M2[4]-df[75][i]),                                      np.abs(M2[3]-df[60][i]),np.abs(M2[2]-df[45][i]),                                   np.abs(M2[1]-df[30][i]),np.abs(M2[0]-df[15][i])]),2)
    
    df['Sep15'][i]=np.round(np.max([np.abs(S2[5]-df[90][i]),np.abs(S2[4]-df[75][i]),                                      np.abs(S2[3]-df[60][i]),np.abs(S2[2]-df[45][i]),                                   np.abs(S2[1]-df[30][i]),np.abs(S2[0]-df[15][i])]),2)
    
    df['Nov15'][i]=np.round(np.max([np.abs(N2[5]-df[90][i]),np.abs(N2[4]-df[75][i]),                                      np.abs(N2[3]-df[60][i]),np.abs(N2[2]-df[45][i]),                                   np.abs(N2[1]-df[30][i]),np.abs(N2[0]-df[15][i])]),2)
print(time()-bi)


# In[8]:


MAR=df.sort_values('Mar15')
MAR=MAR[MAR['Mar15']<0.05]
MAR= MAR.reset_index()
del MAR['index']

SEP=df.sort_values('Sep15')
SEP=SEP[SEP['Sep15']<0.05]
SEP= SEP.reset_index()
del SEP['index']

NOV=df.sort_values('Nov15')
NOV=NOV[NOV['Nov15']<0.05]
NOV= NOV.reset_index()
del NOV['index']

