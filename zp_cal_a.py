# -*- coding:utf-8 -*-
from bct import *
import pandas as pd
import numpy as np
import os
#os.chdir('I:/PAPER/干旱氮添加实验/R_data')



adja=pd.read_table("net.adja.csv",sep=',',index_col=0) 
adja_array = np.array(adja)
louvain=community_louvain(adja_array,gamma=0.9)
ci=np.array(louvain[0]) 
z=centrality.module_degree_zscore(adja_array, ci, flag=0)
p=participation_coef(adja_array, ci, degree='undirected')
sum(z>2.5)
sum(p>0.62)
np.max(ci)
#np.savetxt("zi.txt",z,fmt="%1.4f")
#np.savetxt("pi.txt",p,fmt="%1.4f")
#np.savetxt("module_ci.txt",ci,fmt="%1.4f") 

zipi=np.vstack((np.array(adja.index),z,p,ci))
np.savetxt("zipi.csv",np.transpose(zipi),fmt="%s,%1.4f,%1.4f,%d") 

## fin

