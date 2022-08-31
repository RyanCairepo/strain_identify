#!/usr/bin/bash
DIR="$(cd "$(dirname "$0")" && pwd)"
cp /dev/null combine_"$1"_extract.sam
for f in {56..64};do
    ${DIR}/find_sub.sh -r "$1" -1 ${DIR}/multi_support/${f}_filtered_R1.fastq -2 ${DIR}/multi_support/${f}_filtered_R2.fastq -m tog -a bowtie2 -c Y -d Y -o ${f}_out
    cat ${f}_out/narrowed_extract.sam >> combine_"$1"_extract.sam
    cat combine_"$1"_extract.sam > sorted_combine_"$1"_extract.sam
    cp sorted_combine_"$1"_extract.sam combine_"$1"_extract.sam
done
