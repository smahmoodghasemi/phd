#!/usr/bin/env python
# coding: utf-8

# In[ ]:



# In[2]:


import numpy as np
import pandas as pd
import time
from time import time
import pickle


# In[14]:


def Allele(MTD=48,LIF=10,prop=0.75,C=25,B=120,A=15,T=np.array(range(0,195,15))): 
	'''
 	Input:
	MTD: Mean Transcription Duration (integer)
	A: Mean time between two successive transcriptions after treatment (float)
	LIF: Average lifetime for a nascent mRNA to be a mature mRNA (float)
	C: Mean time between two successive transcriptions after treatment (float)
	B: Time point that a transcription began before the begining of treatment (integer)
	prop: proportion of the length of the gene to be transcribed for visibility (float)
	T: time points to record the number of mRNAs (list of integers)

	Output:
	incomp: list of incomplete transcription for each time point (dictionary)
	comp: list of complete transcription for each time point before degradation (dictionary)


	'''
    
    MTV=prop * MTD   
    TERM=np.max(T)+MTD+5
    size=int((TERM+MTD)/A)
    
    w0=np.random.exponential(scale=C,size=size)
    t0=np.zeros(size)
    for i in range(size):
        t0[i]=sum(w0[0:i+1])
    t0=t0-B
    t0=t0[t0<1]
    w=np.random.exponential(scale=A,size=size) ## Paragraph 6
    t1=np.zeros(size)
    for i in range(size):
        t1[i]=sum(w[0:i+1]) ## Paragraph 7
        
    t=np.append(t0,t1)
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


# In[48]:


gridLIF=np.linspace(19, 22, 7,endpoint=True)  
gridA=np.linspace(5, 25, 41,endpoint=True)
gridC=np.arange(10,36,1)
gridB=np.arange(50,125,5)
gridMTD=np.arange(40,47,1)
VALUE=[(lif,a,c,b,mtd) for lif in gridLIF for a in gridA for c in gridC[gridC>a+5] for b in gridB for mtd in gridMTD] 
sim=np.arange(1,10001)
ST=np.arange(0,135,15)
viz=2


# In[24]:


ti=time()
df=pd.DataFrame()
for val in VALUE:
    name='LIF'+str(val[0])+'A'+str(val[1])+'C'+str(val[2])+'B'+str(val[3])+'MTD'+str(val[4])
    
    NAS1={}

    for i in sim:
        NAS1[i]={}
        _,NAS1[i]=Allele(MTD=val[4],LIF=val[0],A=val[1],C=val[2],B=val[3],prop=0.75,T=ST)    

    simulation=pd.DataFrame(0,index=[name],columns=np.arange(0,135,15))       
    for k in ST:
        for i in sim: 
            if (2<= len(NAS1[i][str(k)])):
                simulation[k]=simulation[k]+1   
            elif (2 > len(NAS1[i][str(k)])):
                simulation[k]=simulation[k]+0
                    
    simulation=(simulation / (len(sim))).round(2)       
    df=df.append(simulation)
        
print(time()-ti)


# In[ ]:


df.to_csv("NF.csv")

