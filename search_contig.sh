#!/bin/bash

while IFS= read -r line; do
	grep "len=$line" -A 1 $1 
done < "length_list.txt"