#!/usr/bin/bash
#python $insec/countbase.py --get_gap_reads=yes files $1 $2 $3 $4 --strain_num $5 |tee  temp_log.txt ; rm -r batch_*

echo ${1}
strain=${1}
#rm -r batch_*
cp /dev/null strain_${strain}_combine_extract.sam 
echo "working on strain ${1}"

cp "subbed_reads_${1}.sam" "subbed_reads_${1}.sam_back"
cp "final_strain_${1}_reference.fa" "final_strain_${1}_reference.fa_back"

for f in {56..64};do
	#ls $insec/multi_support
	#ls $insec/multi_support/"${f}"_filtered_R*
	$insec/find_sub.sh -r final_strain_"${1}"_reference.fa -1 $insec/multi_support/"${f}"_filtered_R1.fastq -2 $insec/multi_support/"${f}"_filtered_R2.fastq -m tog -a bowtie2 -c Y -d Y -o "${f}"_out
	cat "${f}"_out/round_1_extract.sam >> strain_"${1}"_combine_extract.sam
done

python $insec/countbase.py --verify_misp=yes files $insec/402124/oneline_MN996528.1.fasta strain_"${1}"_combine_extract.sam $insec/spike_test/karect_n2_batch/subbed_reads_"${1}".sam $insec/spike_test/karect_n2_batch/final_strain_"${1}"_reference.fa true_single_valid_candidate.sam
