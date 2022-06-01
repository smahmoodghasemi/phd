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


df=pd.read_csv(r'C:\Users\smahm\FVSIM\FVMTD'+str(i)+'.csv')

df.columns=np.append(np.array([-1]),np.arange(0,150,15))
df=df.reset_index()
del df['index']


# In[4]:


#probability of Activation for a single allele ==> December Experiment (Flavopiridol)
D2=[0.05,0.09,0.22,0.43,0.56,0.63] 
#probability of Activation for a single allele ==> January Experiment (Flavopiridol)
J2=[0.04,0.05,0.15,0.37,0.52,0.53]
#probability of Activation for a single allele ==> August Experiment (Flavopiridol)
A2=[0.03, 0.06, 0.14, 0.22, 0.30, 0.36]


# In[5]:



bi=time()
df['MTD']=""
df['A']=""
df['LIF']=""

df['Dec45']=""
df['Jan45']=""
df['Aug45']=""

for i in range(df.shape[0]):

    df['A'][i]=float((np.array(df[-1])[i].split("A")[1]).split("LIF")[0])
    df['LIF'][i]=float((np.array(df[-1])[i].split("LIF")[1]).split("MTD")[0])
    df['MTD'][i]=float((np.array(df[-1])[i].split("LIF")[1]))
    df['Dec45'][i]=np.round(np.max([np.abs(D2[5]-df[90][i]),np.abs(D2[4]-df[75][i]),                                      np.abs(D2[3]-df[60][i]),np.abs(D2[2]-df[45][i])]),2)
    
    df['Jan45'][i]=np.round(np.max([np.abs(J2[5]-df[90][i]),np.abs(J2[4]-df[75][i]),                                      np.abs(J2[3]-df[60][i]),np.abs(J2[2]-df[45][i])]),2)
    df['Aug45'][i]=np.round(np.max([np.abs(A2[5]-df[90][i]),np.abs(A2[4]-df[75][i]),                                      np.abs(A2[3]-df[60][i]),np.abs(A2[2]-df[45][i])]),2)
print(time()-bi)


# In[7]:


DEC=df.sort_values('Dec45')
DEC=DEC[DEC['Dec45']<0.04]
DEC= DEC.reset_index()
del DEC['index']

JAN=df.sort_values('Jan45')
JAN=JAN[JAN['Jan45']<0.04]
JAN= JAN.reset_index()
del JAN['index']

AUG=df.sort_values('Aug45')
AUG=AUG[AUG['Aug45']<0.04]
AUG= AUG.reset_index()
del AUG['index']


# In[ ]:




