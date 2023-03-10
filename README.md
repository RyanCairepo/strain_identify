This program contains code for identifying within-host diversity.

The required python3 packages are:
collections, copy, pandas, numpy, scipy.sparse, typing, argparse, statistics

Steps to run the program:
1. obtain the read files, preferably in fastq format
2. obtain the reference sequence in fasta format
2.1. (optional) Correcting FASTQ read with small_RNA_propor, a compiled C++ program running on linux that does error correction for reads with command "./smallRNA_proper -f read_file". The output of correction is in "changed_list.txt".
3. run the strain identification process with /path/to/find_sub.sh -r reference.fa -1 read_1.fastq -2 read_2.fastq -m tog 
if the read files are in fasta format,
/path/to/find_sub.sh -r reference.fa -1 read_1.fasta -2 read_2.fasta -m tog -f
if the read file is single-ended:
find_sub.sh -r reference.fa -0 read_file.fastq

The output consists of the nucleotide sequences of detected strains, named in the format "final_strain_x_reference.fa", x is the numerical label of strains, and "subbed_read_x.fa", a set of reads that belong to this strains and are different from the reference sequence.

3.1. (optional) step of verification. Start by obtaining relevant samples. After that, run verify.py N original_reference.fa -p(paired end reads, for single end use -s) sample1_r1.fastq sample1_r2.fastq sample2_r1.fastq sample2_r2.fastq.
N is the numerical labelling of strain from step 3.

5. optional step for inferring synonymous state. It determines the changes in the nucleotide sequences are synonymous or non-synonymous. 
python synonymous_stat.py original_reference_sequence.fa subbed_read_N.sam translation_code.txt protein_pos.txt
transation_code.txt is the translation table from nucleotide bases to amino acid bases. protein_pos.txt is the position of proteins, in the form of protein_name:start..end




 protein position file, with protein start and end positions in the format protein:start1..end1, example in protein_pos.txt 
