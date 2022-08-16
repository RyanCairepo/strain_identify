#!/bin/sh

if [ ! -z $1 ]; then
	echo $1
fi
#if [ $# -ge 1 ]; then
#	echo $1
#fi

#lines=$(cat "top40.s3000-3315.txt")
#arr=()
#for line in $lines
#do
#	line1=$((line - 1))
#	tline=$((line+100))
#	echo $line1
	#echo $(awk '{if(NR==$line1){head -$line "naive_match.txt" | awk '{if(NR==$line1)print$3}'print$3};}' "naive_match.txt" )
#	echo $(head -$((line1)) naive_match.txt | tail -1 | awk '{print$3}')
#	arr+=$(head -$((line1)) naive_match.txt | tail -1 | awk '{print$3}')
#done

#echo $arr