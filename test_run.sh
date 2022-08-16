#!/bin/bash
#read_num=$[$(wc -l ./minimap2/gene.sam| cut -d' ' -f1)]
core_count=$[(($( grep -c ^processor /proc/cpuinfo)-2))]

sample_size=300
head -$((sample_size+2)) $1|tail -$sample_size > temp_size.sam && tail -$sample_size $1 >> temp_size.sam
avg=$[(($(ls -l temp_size.sam | awk '{print $5}')/(2*sample_size)))]
estimate=$[(($(ls -l $1 | awk '{print $5}')/$avg))]
divide=$((estimate/core_count))
echo $divide $estimate $avg
#convert fastq to fasta
divide=$((divide+(4-(divide%4))))
echo $divide

time split -l $divide $1 split_fastq
touch $1.fasta
cp /dev/null $1.fasta
rm temp_size.sam
count=0
arr=()

for f in split_fastq* 
do
    touch fasta_$count
    cp /dev/null  fasta_$count
    
    #if [ "$count" == 0 ];then
    #    echo "first file $f asxr_extract_$count"

        #awk 'BEGIN{tmp="'"*"'"}{if(NF>10 && $6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' $f  > multi_extract_$count.sam  &
        #awk 'BEGIN{tmp="'"*"'"}{if(NF>10 && $6!=tmp) print $0}' $f  > multi_extract_$count.sam  &
        #sed -n '1~4s/^@/>/p;2~4p' $f > fasta_$count &
        
        #convert fastq to fasta
        #awk '{if(NR%4==1) {sub("@",">"); printf("%s\n",$1);} else if(NR%4==2) print;}' $f > fasta_$count

        #fix ID
        awk '{if(NR%2==1) {print $1} else print $2}' $f > fasta_$count &
#    else
        echo "$f into fasta_$count" 
#        awk 'BEGIN{tmp="'"*"'"}{if($6!=tmp) print $1" "$2 " "$4 " "$10" "$6}' $f > axsr_extract_$count.sam  &
#    fi
    arr+=($!)
    count=$((count+1))
done

echo ${arr[*]}
echo ${#arr[@]}
running=1
prev=0

for ((x=0;x<$count;x++))
    do
        echo "wait for ${arr[$x]}"
        wait ${arr[$x]}
        echo $(date +"%T")
done
   
echo "all subproc finished"
: '
while [ "$running" -eq 1 ]
do
    for ((x=0;x<$count;x++))
    do
        
        #if [ $pid_run -ne 0 ]; then
        #if pid_run=$(ps -p ${arr[$x]})
     
        if [ -d "/proc/${arr[$x]}" ]
        then
            if [ $prev -ne ${arr[$x]} ]
            then
                prev=${arr[$x]}
                echo "wait for ${arr[$x]}"
            fi
            #wait ${arr[$x]}
            running=1
            break
        else
            #echo "${arr[$x]} has finished"
            #echo $[$(wc -l asxr_extract_"$x".sam| cut -d' ' -f1)]
            running=0
        fi
        

    done
    
done < /dev/null
'
for f1 in fasta_*
do
    
    echo "$[$(grep -c '^' $f1 )] lines from $f1 into $1.fasta"
    cat $f1 >> $1.fasta
done


#cp /dev/null multi_extract_s.sam
#cat multi_extract.sam | sort -nk4 > s_multi_extract.sam

#cp s_multi_extract.sam multi_extract.sam

#cat axsr_extract_* | sort -nk3 >> axsr_extract.sam

for f2 in split_fastq*
do 
    rm $f2
done



echo "$1.fasta with" $[$(grep -c '^' $1.fasta)] "reads"