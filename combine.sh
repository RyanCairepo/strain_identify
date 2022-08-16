#!/usr/bin/bash

for f in ./sra_downloads/*; do
	if echo "$f" | grep 'SRR11.*.fastq'; then
		cat "$f" >> combinesra.fastq
fi
done

