# -*- coding: utf-8 -*-
# Copyright (C) Dan He
#
# This file is part of the code accompanying:
# "Regional environmental heterogeneity under contrasting anthropogenic pressures
#  has differential effects on particle-attached than free-living bacteria
#  communities in coral reef waters."
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# For commercial or non-GPL licensing, please contact:
#   Dan He <hugh_20@163.com>


# -*- coding:utf-8 -*-
from bct import *
import pandas as pd
import numpy as np
import os
#os.chdir('I:/R_data')



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


