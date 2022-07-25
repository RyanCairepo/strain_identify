#!/bin/bash
#set -o errexit
DIR="$(cd "$(dirname "$0")" && pwd)"
#chmod +x ${DIR}/bwa-mem2-2.0pre2_x64-linux/bwa*
#chmod +x ${DIR}/minimap2/minimap2*
core_count=$(($( grep -c ^processor /proc/cpuinfo)-2))
echo "using $core_count threads"


#set -e
#set -o errexit

#---extract commandline arguments---
function  help {
    echo "Usage:"
    echo "Typical use for paired-end reads, align both read files (-1 read1 -2 read2) together (-m tog) with \$aligner (-a \$aligner):"
    echo "run.sh -r reference.fasta -1 first_read_file.fastq -2 second_read_file.fastq -m c -a \$aligner"
    echo "-a:     Aligner of choice: minimap2 bowtie2 bwa-mem2"
    echo "-m:     Alignment mode, tog means align both paired-end read files together, sep means align the
        two read files separately"
    echo "-f:    For fasta reads, add -f "
    echo "-0: -0 read_file.fastq for a single read file"
    exit 2
}
mode="tog"
aligner="bowtie2"
round=5
check_gap=False
out_dir="none"
#echo "args are:"
#echo "$@"
#echo "${#@}"
delete="n"
while getopts ":h:r:1:2:m:0:f:a:n:c:d:o:" option; do
    ((optnum++))
   case ${option} in
      h) # display Help

         help
         exit;;
     \?)
         echo "Invalid option, use run.sh -h for usage"
         exit;;
     r)
         ref=${OPTARG};;
     1)
         read1=${OPTARG};;
     2)
         read2=${OPTARG};;
     m)
         mode=${OPTARG};;
     0)
         read=${OPTARG};;
     f)
        fasta="t"
        echo "$fasta";;
    a)
        aligner=${OPTARG};;
    n)
        round=${OPTARG}
        echo "$round";;
    c)
        check_gap=${OPTARG}
        echo "$check_gap";;
    d)
        delete=${OPTARG};;
    o)
        out_dir=${OPTARG};;

   esac
done


if [ $OPTIND -eq 1 ]; then
    help
fi
shift $((OPTIND -1))






#--- build index----
function build_index {
if ! test -f $ref;then
    echo "reference file $ref not existed"
    exit
fi
echo "building index for $ref"


if [ "$aligner" == "bowtie2"  ]; then
    echo
    #mkdir ${DIR}/bowtie2_index
    #rm ${DIR}/bowtie2_index/*
    ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2-build --large-index $ref ${DIR}/bowtie2_index/reference_index
elif [ "$aligner" == "minimap2" ]; then

    echo "minimap2 no need index"
else
    ${DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index $ref

fi
}


#---alignment---

function  alignment {
        #${DIR}/bwa-mem2/bwa-mem2 mem -T 90 -t 10 $1 $2 > gene.sam;
    #if [ -z $3 ] || [ "$4" == "s" ] || [ -z $4 ]; then
    echo "start read extraction"
    if [ -e "$read"  ]; then
        echo "align single read file"
        if [ $aligner == "minimap2" ]; then
            ${DIR}/minimap2/minimap2 -ax sr "$ref" "$read" > gene.sam
        elif [ $aligner == "bowtie2" ]; then
            ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -p ${core_count} -x ${DIR}/bowtie2_index/reference_index -i "$read" > gene.sam
        else
            ${DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t $core_count "$ref" "$read"  > gene.sam
        fi

    elif [ -e "$read1" ] && [ -e "$read2" ];then
        if [ $mode != "sep" ]; then
            echo "align both read files together, $fasta"
            if [ "$aligner" == "bowtie2" ]; then
                if [ -e "$fasta" ] && [ $fasta == "t" ]; then
                   echo "bowtie2 fasta alignment"
                   ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -p ${core_count} -x ${DIR}/bowtie2_index/reference_index -f -1 "$read1" -2 "$read2" --local > gene.sam
                else
                    echo "bowtie2 alignment"
                    ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -p ${core_count} -x ${DIR}/bowtie2_index/reference_index -1 "$read1" -2 "$read2" --local > gene.sam
                fi
            elif [ "$aligner" == "minimap2" ]; then
                echo "minimap2 alignment"
                ${DIR}/minimap2/minimap2 -ax sr "$ref" "$read1" "$read2" > gene.sam
            else
                echo "bwa-mem2 alignment"
                ${DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t ${core_count}  "$ref" "$read1" "$read2" > gene.sam
            fi
        else
        echo "align read file $read1 and $read2 separately"
            if [ $aligner == "minimap2" ]; then
                ${DIR}/minimap2/minimap2 -t ${core_count} -ax sr "$ref" "$read1" > gene.sam
                ${DIR}/minimap2/minimap2 -t ${core_count} -ax sr "$ref" "$read2" > gene2.sam
            elif [ $aligner == "bowtie2" ]; then
                ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -p ${core_count} -x ${DIR}/bowtie2_index/reference_index -i "$read1" > gene.sam
                ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -p ${core_count} -x ${DIR}/bowtie2_index/reference_index -i "$read1" > gene2.sam
            else
                ${DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t $((core_count)) "$ref" "$read1"  > gene.sam
                ${DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t $((core_count)) "$ref" "$read1"  > gene2.sam
            fi
        fi
    else
        echo "Error getting read files, use -0 for single read file, use -1 read1.fq -2 read2.fq for paired-end reads"
        exit
    fi
}
function split_sam {
    gene_size=$[$(stat --printf="%s" $1)]
    if [ "$gene_size" -eq 0 ]; then
        echo "error empty alignment result, exiting"
        exit
    fi

    #--- sample about 600 reads to estimate the lines of gene.sam useful for very large gene.sam ---
    sample_size=300
    head -$((sample_size+2)) $1|tail -$sample_size > temp_size.sam && tail -$sample_size $1 >> temp_size.sam
    #sample_size=$[$(wc -l temp_size.sam|cut -d' ' -f1 )]
    avg=$[(($(ls -l temp_size.sam | awk '{print $5}')/(sample_size*2)))]
    estimate=$[(($(ls -l $1 | awk '{print $5}')/avg))]
    divide=$((estimate/core_count))

    echo $divide $estimate $avg
    if [ -z $divide ] || [ -z $estimate ] || [ -z $avg ]; then
        echo "error values, d $divide e $estimate a $avg c $core_count"
        exit
    fi
    time split -l $divide $1 split_gene

    rm temp_size.sam

}

function extract_sam {

    touch extract.sam
    cp /dev/null extract.sam
    cp /dev/null full_extract.sam
    count=0
    arr=()

    for f in split_gene*
    do
        touch extract_$count.sam
        cp /dev/null extract_$count.sam
        cp /dev/null full_extract_$count.sam
        echo "split gene.sam into files $count"

        awk 'BEGIN{tmp="'"*"'"}{if(NF>10 && $6!=tmp) print $0}' $f  > full_extract_$count.sam  &
        #awk 'BEGIN{tmp="'"*"'"}{if(NF>10 && $6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' $f  > extract_$count.sam  &

        arr+=($!)
        count=$((count+1))
    done

    for ((x=0;x<=$count;x++))
    do

        wait ${arr[$x]}
    done
    echo "full_extract_ ready"
    arr=()

    x=0
    for full1 in full_extract_*
    do
        cat $full1 >> full_extract.sam
        awk '{ print $1" "$2 " "$4 " "$10" "$6}' $full1  > extract_$x.sam  &
        arr+=($!)
        x=$((x+1))
    done

    for ((x=0;x<=$count;x++))
    do
        echo "wait, $(wc -l "full_extract_$x.sam") lines in full_extract_$x.sam "
        wait ${arr[$x]}
    done
    echo "extract_ ready"
    for f1 in extract_*
    do
        echo "$[$(grep -c '^' $f1 )] lines from $f1 into extract.sam"
        cat $f1 >>extract.sam
    done

    for f2 in split_gene*
    do
        echo "remove $f2"
        rm "$f2"
    done

    rm extract_*
    rm full_extract_*
    if [   $[$(wc -l extract.sam|cut -d' ' -f1)] -eq 0 ];then
        echo "no matched reads"
        exit
    fi
    cat full_extract.sam | sort -nk4 > s_full_extract.sam
    cp s_full_extract.sam full_extract.sam

    cat extract.sam | sort -nk3 > s_extract.sam
    cat s_extract.sam > extract.sam
    #remove firstline of bowtie2 output that causes error

    if [ "$aligner" != "minimap2" ]; then
        sed -i 1d extract.sam
#        sed -i '$ d' extract.sam
    fi

    echo  $[$(wc -l extract.sam|cut -d' ' -f1)] " total reads"
    #awk 'BEGIN{tmp="'"*"'"}{if(NR>2 && $6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' gene.sam |sort -nk3 > extract.sam;
    # awk 'match($12,/[1-9][0-9]*/) {print $0}' gene.sam quality score


    printf "\n\n produce extrat \n rl \"${rL}\" \n"
    echo " extract.sam  with" $[$(wc -l extract.sam|cut -d' ' -f1)] "reads"
}

cp /dev/null find_sub_log.txt
len_limit=$(bc<<<$( tr -d '\r' < "$ref" | awk 'BEGIN{RS="\n"}{if(substr($0,1,1) ~ /\>/ ) printf("%s", $NF);next}  '  | wc -m)\*0.7)
for i in {1..1};do
    echo "$i $ref $len_limit"
    #out_dir=round_"$i"_output
    #contig_dir=round_"$i"_contig
    if  [ "$out_dir" == "none" ]; then
        out_dir="${ref##*/}"_output
    fi

    echo output dir is  "$out_dir"
    #rm split_gene*

    fL=$( tr -d '\r' < "$ref" | awk 'BEGIN{RS="\n"}{if(substr($0,1,1) ~ /\>/ ) printf("%s", $NF);next}  '  | wc -m)
    if (( $(echo "$fL < $len_limit" | bc -l) )); then
        echo "gene length $fL too short"
        exit
    fi

    if [ -d "$out_dir" ]; then
        echo remove "$out_dir" ? Y/N

        if [ $delete == "Y" ]; then
            if [ -d "$out_dir" ]; then
                mv "$out_dir" "$out_dir"_del
                rm -r "$out_dir"
                mkdir "$out_dir"
            fi
#            if [ -d "$contig_dir" ]; then
#
#                rm -r "$contig_dir"
#            fi
        else
            exit
        fi
    else
        mkdir "$out_dir"
    fi

    echo "$ref" "$read1" "$read2" "$mode" "$check_gap" "$aligner" > "$out_dir"/args.txt


    build_index
    alignment

    if [ ! -e "$read" ] && [ "$mode" != "sep" ]; then
        echo "get matched reads"
        split_sam gene.sam
        extract_sam
    elif [ "$mode" == "sep" ] && [ -e "$read" ]; then
        extract_sam gene.sam
        cp extract.sam 1_extract.sam
        cp full_extract.sam 1_full_extract.sam
        split_sam gene2.sam
        extract_sam
        cat extract.sam >> 1_extract.sam
        cat full_extract.sam >> 1_full_extract.sam
        cp 1_extract.sam extract.sam
        cp 1_full_extract.sam full_extract.sam
    fi
    if [ ! -d intermit_out ]; then
        mkdir intermit_out
    fi
    cp extract.sam "$out_dir"/round_"$i"_extract.sam
    cp full_extract.sam "$out_dir"/round_"$i"_full_extract.sam
    #--- run correction----
    #fL=$( awk 'BEGIN{RS="\r\n"}{if(substr($0,1,1) ~ /\>/ ) printf("%s", $NF);next}  ' "$ref" | wc -m)
    rL=$(awk '{print length($4)}' extract.sam | sort -n | tail -1)
    a="--gene_L=${fL}";
    b="--read_L=${rL}";
    #ref=$(<../${1})
    #echo $ref
    match_limit=1 #$(echo "scale=2; ((100.0-$i+1)/100)" |bc -l)
    echo "python3 strain_finder.py $a  --ref=${ref} --narrowing=True --match_l=$match_limit --sam_file=extract.sam --r1_file=${read1} --r2_file=${read2} --round=$i --excluded_IDs=/dev/null --find_sub=True --bubble_mode=True --check_gap=${check_gap} --output_dir=${out_dir}";

    if [ $i -eq 1 ]; then
        cp /dev/null  "excluded_IDs.txt"
    fi

    if [ -e "$read1" ] && [ -e "$read2" ]; then
        python3 "${DIR}"/strain_finder.py $a $b --ref="${ref}"  --narrowing=True --match_l=${match_limit} --sam_file=extract.sam --r1_file="$read1" --r2_file="$read2" --round="$i" --excluded_IDs="excluded_IDs.txt" --find_sub=True --brute_force=True --check_gap="$check_gap" --output_dir="$out_dir";
    else

       python3 "${DIR}"/strain_finder.py $a $b --ref="${ref}"  --narrowing=True --match_l=${match_limit} --sam_file=extract.sam --r1_file="$read" --round="$i" --excluded_IDs="excluded_IDs.txt" --find_sub=True ;
    fi
    err=$?

    if [ $err -eq 5 ] ; then
        echo "unable to form contig"

    fi

    if [ $err -ne 0 ] && [ $err -ne 5 ] ; then
        echo " Error in python program $err"
        exit
    fi
    subamount=$(wc -l "$out_dir"/paired_real_narrowed_extract.sam)
    echo "$out_dir"/paired_real_narrowed_extract.sam $subamount >> find_sub_log.txt
    cp narrowed_cvg.txt "$out_dir"/narrowed_cvg.txt
    cp real_narrowed_cvg.txt "$out_dir"/
    cp nearly_real_narrowed_cvg.txt "$out_dir"/
    cp thin_* "$out_dir"/
    cp mutated_read_freq* "$out_dir"/


#    "${DIR}"/megahit/build/megahit -1 "$out_dir"/half_real_R1.fastq -2 "$out_dir"/half_real_R2.fastq -o "$contig_dir"
#
#    if [ "$(wc -l "$contig_dir"/final.contigs.fa | cut -d' ' -f1)" -gt 2 ]; then
#        echo "more than one contig formed"
#        longest=$(grep -oEi 'len=([0-9]+)' -A 1 "$contig_dir"/final.contigs.fa | cut -d '=' -f2 | sort -n | tail -1)
#        grep $longest -A 1 "$contig_dir"/final.contigs.fa > "$out_dir"/round_"$i"_contig.fa
#    else
#        cp "$contig_dir"/final.contigs.fa "$out_dir"/round_"$i"_contig.fa
#    fi
#    cL=$(tr -d '\r' < "$out_dir"/round_"$i"_contig.fa | awk 'BEGIN{RS="\n"}{if(substr($0,1,1) ~ /\>/ ) printf("%s", $NF);next}  '  | wc -m)
#    if [ $cL -lt 25000 ]; then
#        echo "contig formed too short, $cL"
#        exit
#    fi
#    ref="$out_dir"/round_"$i"_contig.fa


done
