#!/usr/bin/bash
DIR="$(cd "$(dirname "$0")" && pwd)"
cp /dev/null combined_extract.sam
while getopts ":s:p:" option; do
    case ${option} in
    p)
        echo "paired_end mode"
        for ((i = 3; i < (($#));i+=2)); do

            j=$((${i}+1))
            echo "$i" "${!i}" "$j" "${!j}"
            #echo "${@[$f]}" "${@[(($f+1))]}"

            ${DIR}/find_sub.sh -r "$2" -1 "${!i}" -2 "${!j}" -m tog -a bowtie2 -c Y -d Y -o ${i}_out
            cat ${i}_out/extract.sam >> combine_extract.sam

        done;;
    s)
        echo "single_end mode"
        for f in "${@:3}"; do
            echo $f
             ${DIR}/find_sub.sh -r "$2" -0 $f -c Y -d Y -o ${i}_out
             cat ${i}_out/extract.sam >> combine_extract.sam
        done;;
    \?)
        echo "-s for single end read files, -p for paired end read files"
        echo "when using paired end read files, place the paired files from the same sample together, eg. combine_align.sh reference.fasta SRR11092061_1.fastq SRR11092062_2.fastq SRR11092063_1.fastq SRR22092063_2.fastq"
        exit;;
    esac
done