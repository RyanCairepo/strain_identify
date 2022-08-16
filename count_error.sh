#!/usr/bin/bash

for f in ./sra_downloads/*.fastq; do
	echo "$f" | grep -oh 'SRR11[^.]*' >> errors.txt
	./run.sh ../EPI_ISL_402124-ORF1ab.fasta "$f" | grep -oh '.*errors' >> errors.txt

	done
