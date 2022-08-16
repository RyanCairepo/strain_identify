#!/bin/bash
chmod +x bwa-mem2-2.0pre2_x64-linux/bwa*

echo "start read extraction"

fL=$[$(awk '{if(NR>1)print}' ${1}|wc -m)-$(awk '{if(NR>1)print}' ${1}|wc -l)-100]
rL=$[$(awk '{if(NR==2)print$1}' ${2}|wc -m)-1]

#--- build index----
bwa-mem2-2.0pre1_x64-linux/bwa-mem2 index $1 

#--- map by BWA-MEM ---
bwa-mem2-2.0pre1_x64-linux/bwa-mem2 mem -T 90 $1 $2 > gene.sam;
awk 'BEGIN{tmp="'"${rL}M"'"}{if(NR>2 && $6==tmp) print $1" "$2 " "$4 " "$10}' gene.sam |sort -nk3 > extract.sam;

#--- run correction----
a="--gene_L=${fL}";
b="--read_L=${rL}";
echo "python3 0_error_cor.py $a $b";
python3 error_cor.py $a $b;