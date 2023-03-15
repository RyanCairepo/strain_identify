
# Proportional correction workflow
- **step 0**(optional): compile smallRNA_proper.c with any C++ compiler if the smallRNAproper cannot run on your system, it should work on any linux system. 

- **Step 1**: Given the initial dataset $D_0$, use the command "./smallRNA_proper -f read_file" to perform a proportional correction and obtain the corrected dataset $D_1$. The list of changes made during the correction process will be stored in the "changed_list.txt" file. If you have multiple read files, make sure run smallRNAproper in DIFFERENT directories.

# Singletons correction workflow

- **Step 2**: We then split the corrected dataset $D_1$ into two sub-datasets: singleton records (i.e., records with a frequency of 1) which are stored in $D_2$, and non-singleton records which are stored in $D_3$.

- **Step 3**: Next, we correct errors in the singleton dataset $D_2$ using the program karect, but we keep only the corrected reads that differ from their original counterparts by one base, and save them as the corrected dataset $D_2^{\prime}$.

- **Step 4**: Finally, we merge the corrected dataset $D_2^{\prime}$ and the non-singleton dataset $D_3$ to obtain the final corrected dataset.

### Dependency

biopython (pip install biopython==1.79)

editdistance (pip install editdistance==0.6.0)

karect (https://github.com/aminallam/karect)

### run singleton correction
Install the dependency in your environment and set input files in the script of singleton_correction.py, then run 

`python singleton_correction.py -i <input.fastq> -o <output_dir> -k <./karect>`

or

`python singleton_correction.py --input <input.fastq> --output_dir <output_dir> --karect <./karect>`