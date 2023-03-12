# Dependency

biopython (pip install biopython==1.79)

editdistance (pip install editdistance==0.6.0)

networkx (pip install networkx==2.8.5)

mpire (pip install mpire==2.6.0)

tqdm (pip install tqdm==4.64.0)

karect (https://github.com/aminallam/karect)

# Isolates correction workflow

- **Step 1**: Begin with a smallRNA_propor corrected dataset $D$. To distinguish between high and low-frequency reads, we apply a predetermined frequency threshold, denoted as τ. Divide the reads into two groups accordingly: high-frequency and low-frequency.

- **Step 2**: Identify all 1-edit-distance relationships between the high-frequency and low-frequency reads. Construct a read graph with sequences as nodes and relationships as edges.

- **Step 3**: Based on the isolated sequences in the read graph, separate the original dataset $D$ into two sets: isolates contained dataset $D1$ and non-isolates contained dataset $D2$.

- **Step 4**: Correct dataset $D1$ using karect, resulting in the corrected dataset $D1'$. Merge the corrected dataset $D1'$ and $D2$ to obtain the final corrected dataset.

# run
Install the dependency in your environment and set parameters and input in the script of isolates_correction.py, then run 
`python isolates_correction.py`