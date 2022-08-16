#!/bin/bash
DIR="$(cd "$(dirname "$0")" && pwd)"
#chmod +x ${DIR}/bwa-mem2-2.0pre2_x64-linux/bwa*
#chmod +x ${DIR}/minimap2/minimap2*

echo "start read extraction"

#fL=$[$( awk 'BEGIN{RS="\r\n"}{if(NR>1) printf("%s", $NF);next}  ' $1  |wc -m)]


#--- build index----
if [ -e rebuild_index=$6 ];then
	rebuild_index=$6
else
	rebuild_index=no
fi
#

fasta="f"
if [ "$#" -ge 5 ] ; then
    if [ "$5" == "bowtie2" ]; then
	if [ -d ${DIR}/s_bowtie2_index ]; then
        	if [ "$rebuild_index" == "yes" ];then
			echo "rebuild index on s_bowtie2_index, removing previous"
	        	#exit
			rm ${DIR}/s_bowtie2_index/*	
        		${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2-build --large-index $1 ${DIR}/s_bowtie2_index/reference_index	
		else
			echo " not rebuilding index, do nothing"
			#exit
		fi
	else
		echo "no previous index exist,exit"
        	#exit
		mkdir ${DIR}/s_bowtie2_index
        	#rm ${DIR}/s_bowtie2_index/*
        	${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2-build --large-index $1 ${DIR}/s_bowtie2_index/reference_index
	fi	
    fi
fi
core_count=$(($( grep -c ^processor /proc/cpuinfo)))

#--- map by BWA-MEM ---
#${DIR}/bwa-mem2/bwa-mem2 mem -T 90 -t 10 $1 $2 > gene.sam;
#if [ -z $3 ] || [ "$4" == "s" ] || [ -z $4 ]; then
if [ "$#" -eq 3 ] || [[ "$#" -eq 4  &&  "$4" == "s" ]]; then
    echo "align single read file"
    #${DIR}/minimap2/minimap2 -ax sr $1 $2 > gene.sam
#    if [ "$#" -eq 6 ] ; then
#        if [ "$6" == "bowtie2" ]; then
#            ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -x ${DIR}/bowtie2_index/reference_index -i $2 > gene.sam
 #           none
  #      fi
  #  fi


elif [ "$#" -ge 4 ] && [ "$4" == "c" ];then
    echo "align both read files together, $fasta"

    if [ "$#" -ge 5 ] ; then
        if [ "$5" == "bowtie2" ]; then
            if [ $fasta == "f" ]; then
           
               echo "s bowtie2 alignment local"
               ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -p ${core_count} -x ${DIR}/s_bowtie2_index/reference_index -1 $2 -2 $3 --local > gene.sam
            
            fi
        elif [ "$5" == "minimap2" ]; then
            ${DIR}/minimap2/minimap2 -t ${core_count} -ax sr $1 $2 $3 > gene.sam
        else
            echo "bwa-mem2 alignment"
            ${DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t ${core_count} $1 $2 $3 > gene.sam

        fi

    fi


fi
gene_size=$[$(stat --printf="%s" gene.sam)]
if [ "$gene_size" -eq 0 ]; then
    echo "error empty alignment result, exiting"
    exit
fi


#--- sample about 600 reads to estimate the lines of gene.sam useful for very large gene.sam ---
sample_size=300
head -$((sample_size+2)) gene.sam|tail -$sample_size > temp_size.sam && tail -$sample_size gene.sam >> temp_size.sam
#sample_size=$[$(wc -l temp_size.sam|cut -d' ' -f1 )]
avg=$[(($(ls -l temp_size.sam | awk '{print $5}')/(sample_size*2)))]
estimate=$[(($(ls -l gene.sam | awk '{print $5}')/avg))]
divide=$((estimate/core_count))

echo $divide $estimate $avg
if [ -z $divide ] || [ -z $estimate ] || [ -z $avg ]; then
    echo "error values, d $divide e $estimate a $avg c $core_count"
    exit
fi
time split -l $divide gene.sam split_gene
touch extract.sam
cp /dev/null extract.sam
cp /dev/null full_extract.sam
rm temp_size.sam
count=0
arr=()

for f in split_gene*
do
    touch extract_$count.sam
    cp /dev/null extract_$count.sam
    echo "split gene.sam into files $count"

    awk 'BEGIN{tmp="'"*"'"}{if(NF>10 && $6!=tmp) print $0}' $f  > human_extract_$count.sam  &
    awk 'BEGIN{tmp="'"*"'"}{if(NF>10 && $6==tmp) print $0}' $f  > full_extract_$count.sam  &

    arr+=($!)
    count=$((count+1))
done

for ((x=0;x<=$count;x++))
do

    wait ${arr[$x]}
done

arr=()
x=0
cp /dev/null full_extract.sam
for full1 in full_extract_*
do
    cat $full1 >> full_extract.sam
        awk '{print $1" "$2 " "$4 " "$10" "$6}' $full1  > extract_$x.sam  &
	arr+=($!)
	x=$((x+1))
done


for ((x1=0;x1<=$x;x1++))
do
    wait ${arr[$x1]}
done

for f1 in extract_*
do
    echo "$[$(grep -c '^' $f1 )] lines from $f1 into extract.sam"
    cat $f1 >>extract.sam
done

for f2 in split_gene*
do
	echo "remove $f2"    
	#rm $f2
done
echo "first part extract.sam  with" $[$(wc -l extract.sam|cut -d' ' -f1)] "reads"

if [   $[$(wc -l extract.sam|cut -d' ' -f1)] -eq 0 ];then
    echo "no matched reads"
    exit
fi

#cat full_extract.sam | sort -nk4 > s_full_extract.sam
#cp s_full_extract.sam full_extract.sam

#cat extract.sam | sort -nk3 > s_extract.sam
#cat s_extract.sam > extract.sam

#remove firstline of bowtie2 output that causes error
if [ "$#" -ge 5 ] ; then
    if [ "$5" == "bowtie2" ]; then
        head -2 extract.sam
	#sed -i 1d extract.sam
#        sed -i '$ d' extract.sam
    fi
fi

echo  $[$(wc -l extract.sam|cut -d' ' -f1)] " total reads"
echo "get non-human reads to non_human_extract.sam"

awk '{if($4 !~ /^[N]+$/) print $0}' extract.sam > non_human_extract.sam
#python ../get_ori.py non_human_extract.sam $2 $3
#wc -l non_human_extract.sam
exit

#awk 'BEGIN{tmp="'"*"'"}{if(NR>2 && $6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' gene.sam |sort -nk3 > extract.sam;
# awk 'match($12,/[1-9][0-9]*/) {print $0}' gene.sam quality score

rL=$(awk '{print length($4)}' extract.sam | sort -n | tail -1)
printf "\n\n produce extract \n rl \"${rL}\" \n"


#--- run correction----
a="--gene_L=${fL}";
b="--read_L=${rL}";
#ref=$(<../${1})
#echo $ref
echo "python3 error_cor.py $a  --ref=${1}--narrowing=True --match_l=1 --sam_file=extract.sam --r1_file=${2} --r2_file=${3}";
#python3 ${DIR}/strain_finder.py $a $b --ref=${1}  --narrowing=True --match_l=1 --sam_file=extract.sam --r1_file="$2" --r2_file="$3";
exitcode=$?
#${DIR}/megahit/build/megahit -1 half_real_R1.fastq -2 half_real_R2.fastq -o first_half_real_contig
