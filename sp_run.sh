#!/bin/bash
DIR="$(cd "$(dirname "$0")" && pwd)"
#chmod +x ${DIR}/bwa-mem2-2.0pre2_x64-linux/bwa*
#chmod +x ${DIR}/minimap2/minimap2*

echo "start read extraction"

#fL=$[$(awk '{if(NR>1)print}' ${1}|wc -m)-$(awk '{if(NR>1)print}' ${1}|wc -l)-100] #why fixed 100? the length of first line could vary
fL=$[$( awk 'BEGIN{RS="\r\n"}{if(NR>1) printf("%s", $NF);next}  ' $1  |wc -m)]
#rL=$[$(awk '{if(NR==2)print$1}' ${2}|wc -m)-1]

#$[$( awk '{if(NR%2==0)print length($1)}' ${2} | sort -n | tail  -1)]

#--- build index----
#${DIR}/bwa-mem2/bwa-mem2 index $1

#--- map by BWA-MEM ---
#${DIR}/bwa-mem2/bwa-mem2 mem -T 90 -t 10 $1 $2 > gene.sam;
if [ -z $3 ] || [ "$4" == "s" ] || [ -z $4 ]; then
    echo "align single read file"
    ${DIR}/minimap2/minimap2 -ax sr $1 $2 > gene.sam
elif [ "$4" == "c" ];then
    echo "align both read files together"
    ${DIR}/minimap2/minimap2 -ax sr $1 $2 $3 > gene.sam
fi

touch extract.sam

awk 'BEGIN{tmp="'"*"'"}{if(NR>2 && $6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' $f |sort -nk3 > extract_$count.sam  &

echo "first part extract.sam  with" $[$(wc -l extract.sam|cut -d' ' -f1)] "reads"
if [ ! -z $3 ] && [ $4 == "s" ] ;then

    ${DIR}/minimap2/minimap2 -ax sr $1 $3 > gene_2.sam


    touch extract2.sam
    cp /dev/null extract2.sam
 
    awk 'BEGIN{tmp="'"*"'"}{if(NR>2 && $6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' $f |sort -nk3 > extract2_$count.sam  &
 
    echo "second part extract2.sam with" $[$(wc -l extract2.sam|cut -d' ' -f1)] "reads"
    cat extract2.sam >> extract.sam
fi

if [ $[$(wc -l extract2.sam)] == 0 ] && [ $[$(wc -l extract.sam)] == 0 ];then
    echo "no matched reads"
    exit
fi



#awk 'BEGIN{tmp="'"*"'"}{if(NR>2 && $6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' gene.sam |sort -nk3 > extract.sam;
# awk 'match($12,/[1-9][0-9]*/) {print $0}' gene.sam quality score

rL=$(awk '{print length($4)}' extract.sam | sort -n | tail -1)
printf "\n\n produce extrat \n rl \"${rL}\" \n"


#--- run correction----
a="--gene_L=${fL}";
b="--read_L=${rL}";
#ref=$(<../${1})
#echo $ref
echo "python3 error_cor.py $a $b --ref=${1}";
python3 ${DIR}/error_cor.py $a $b --ref=${1};
