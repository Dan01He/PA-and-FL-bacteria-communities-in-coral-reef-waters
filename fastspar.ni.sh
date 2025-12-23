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
