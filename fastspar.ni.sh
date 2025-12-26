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


#!/bin/bash  
Ni="FL PA"  
for item in $Ni; do  
  # 检查目录存在性  
  if [ -d "/home/ubt/fastspar_data/$item" ]; then  
    cd /home/ubt/fastspar_data/$item  
  else  
    echo "Directory /home/ubt/fastspar_data/$item does not exist!"  
    exit 1  
  fi  

  # 创建所需目录  
  mkdir -p bootstrap_counts  
  mkdir -p bootstrap_correlation  

  # 检查输入文件是否存在  
  if [ -f "/home/ubt/fastspar_data/${item}/${item}_otbtp2k.tsv" ]; then  
   fastspar --otu_table /home/ubt/fastspar_data/${item}/${item}_otbtp2k.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv --iterations 50 --exclude_iterations 20 --threshold 0.2 --threads 10
   fastspar_bootstrap --otu_table /home/ubt/fastspar_data/${item}/${item}_otbtp2k.tsv --number 1000 --prefix bootstrap_counts/${item}  
  else  
    echo "File /home/ubt/fastspar_data/${item}/${item}_otbtp2k.tsv does not exist!"  
    exit 1  
  fi  

  # 执行 parallel  
  parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 10 ::: bootstrap_counts/*  
  fastspar_pvalues --otu_table /home/ubt/fastspar_data/${item}/${item}_otbtp2k.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_${item}_ --permutations 1000 --outfile pvalues.tsv

  rm -fr /home/ubt/fastspar_data/${item}/bootstrap_c*
done
