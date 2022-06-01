#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import time
from time import time
import pickle


# In[10]:


def Allele(MTD=48,A=15,LIF=10,prop=0.75,T=np.array(range(0,195,15))): 
	'''
 	Input:
	MTD: Mean Transcription Duration (integer)
	A: Mean time between two successive transcriptions (float)
	LIF: Average lifetime for a nascent mRNA to be a mature mRNA (float)
	prop: proportion of the length of the gene to be transcribed for visibility (float)
	T: time points to record the number of mRNAs (list of integers)

	Output:
	incomp: list of incomplete transcription for each time point (dictionary)
	comp: list of complete transcription for each time point before degradation (dictionary)


	'''
    
    MTV=prop * MTD   
    TERM=np.max(T)+MTD+5
    size=int((TERM+MTD)/A)
    

    w=np.random.exponential(scale=A,size=size) ## Paragraph 6
    t=np.zeros(size)
    for i in range(size):
        t[i]=sum(w[0:i+1]) 

    STOP=np.count_nonzero(t <= TERM) 
    
    S=np.random.exponential(scale=LIF,size=STOP) ## Paragraph 11
    D=np.zeros(STOP)  
    EMERGE=np.zeros(STOP)             

    for j in range(STOP):
        EMERGE[j]=t[j]+MTV
        D[j]=t[j]+MTD+S[j] 
    
    incomp={}
    comp={}
    for k in T:
        incomp[str(k)]=[]
        comp[str(k)]=[]
        
        for j in range(STOP):
            if (t[j]<=k<EMERGE[j] ):
                incomp[str(k)].append(j) 
            elif (EMERGE[j]<=k<D[j]):
                comp[str(k)].append(j)             
    return incomp , comp  


# In[23]:


gridLIF=np.linspace(19, 22, 7,endpoint=True)  
gridA=np.linspace(5, 25, 41,endpoint=True)
gridMTD=np.arange(40,47,1)
VALUE=[(lif,a,mtd) for lif in gridLIF for a in gridA for mtd in gridMTD ] 
sim=np.arange(1,10001)
ST=np.arange(0,135,15)
MTD=np.arange(40,47,1)

viz=2


# In[24]:


ti=time()
df=pd.DataFrame()
for val in VALUE:
    name='A'+str(val[1])+'LIF'+str(val[0])+'MTD'+str(val[2])
    
    NAS1={}

    for i in sim:
        NAS1[i]={}
        _,NAS1[i]=Allele(MTD=val[2],A=val[1],LIF=val[0],prop=0.75,T=ST)    

    simulation=pd.DataFrame(0,index=[name],columns=np.arange(0,135,15))       
    for k in ST:
        for i in sim: 
            if (viz<= len(NAS1[i][str(k)])):
                simulation[k]=simulation[k]+1   
            elif (viz > len(NAS1[i][str(k)])):
                simulation[k]=simulation[k]+0
                    
    simulation=(simulation / (len(sim))).round(2)       
    df=df.append(simulation)
        
print(time()-ti)


# In[ ]:


df.to_csv("FV.csv")

