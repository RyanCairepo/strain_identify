# update in Nov

import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_alg
import re
import argparse
import json

#set arguments
def get_arguments():
    parser = argparse.ArgumentParser(description='EC')
    parser.add_argument("--gene_L", type=int, default=4892)
    parser.add_argument("--read_L", type=int, default=100)
    parser.add_argument("--sam_file", type=str, default='extract.sam')
    parser.add_argument("--write_file", type=str, default='corrected.fa')
    parser.add_argument("--write_file_seq", type=str, default='contigs.fa')
    return parser.parse_args()

#convert to reverse_complement seq
def DNA_Re_complement(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    comp_s = sequence.upper()
    return comp_s[::-1]

#statistic proportion of bases in each column
def get_baseprop(tmp_count):
    base_porp = np.array([each/sum(tmp_count) for each in tmp_count])
    return base_porp

    
#%% main function is here-----------------------------------------------------

args = get_arguments()
print('Running parameters:\n')
print(json.dumps(vars(args), indent=4, separators=(',', ':')))


R_file = args.sam_file#chage next line 'cigar'
r = pd.read_csv(R_file, delimiter =' ', names = ['ID','strand','sta_p','sam_q', 'cigar'],encoding = 'unicode_escape')

dict_tonu = {'A':1,'C':2,'T':3,'G':4,'N':5}
dict_tole = dict(zip(dict_tonu.values(),dict_tonu.keys()))

#setting
gene_length = args.gene_L
read_length = args.read_L
read_number = r.shape[0]

row_l = []
col_l = []
val_l = []

#relocated seq to matrix A
for i in range(read_number):
    index = int(r['sta_p'].loc[i])
    sam_q = list(r['sam_q'].loc[i])
    sam_q_num = []
    for j in range(len(sam_q)):
        sam_q_num.append(dict_tonu[sam_q[j]])
    #belwo two lines added pad
    if len(sam_q_num)<read_length:
        sam_q_num.extend((read_length-len(sam_q_num)) * [0])

    val_l.extend(sam_q_num)
    row_tmp = [int(i) for n in range(read_length)]
    row_l.extend(row_tmp)
    col_tmp = [n for n in range(index,index+read_length)]
    col_l.extend(col_tmp) 

row = np.array(row_l)
col = np.array(col_l)
val = np.array(val_l)

csc = sp.coo_matrix((val, (row, col))).tocsc()
lil = sp.coo_matrix((val, (row, col))).tolil()

p_row_l = []
p_col_l = []
p_val_l = []

Seq_tru = []
Seq_let = []
Num_lowfre = 0
Num_col = 0
Num_base_porp = 0
num_col_list = []
#labeled reads with different bases with the dominant base,but unchanged
read_lab = np.full((read_number),0,dtype=int)
read_lab_num = 0

#-------------------------compute probabilities of each base
#count unknow bases in the updated sequence
non_update=0
#change 51 to 0
for i in range(0,gene_length+0):
    # go through each columntm
    tmp = np.squeeze(csc.getcol(i).toarray())
    tmp_count = np.bincount(tmp)[1:]
    #col_num store the value of total_coverage_col
    col_num = sum(tmp_count)
    #maximum letter's number is max_l + 1
    if(tmp_count.shape[0] == 0):
        #0 means this base is unknown
        non_update = non_update+1
        Seq_tru.append(5)
        Seq_let.append(dict_tole[5])
        continue
    else:      
        nonzero_num = np.nonzero(tmp_count)[0].shape[0]
        max_l = np.argmax(tmp_count)
        max_num = tmp_count[max_l]
        #update letter of sequence
        Seq_tru.append(max_l+1)
        Seq_let.append(dict_tole[max_l+1])
        
        #whether there are errors in the column
        if(nonzero_num == 1):
            #no error in this column, only one dominant base
             continue
        else:
            tmp_count[tmp_count == 0] = read_number
            if(col_num < 100):
                #according to base count to decide the dominant base
                if(np.min(tmp_count) < 3):
                    min_l = np.argmin(tmp_count)
                    letter = min_l + 1
                    for k in np.where(tmp == letter)[0]:
                        p_row_l.append(k)
                        p_col_l.append(i)
                        prob = tmp_count[min_l]/max_num * tmp_count[min_l]/col_num
                        p_val_l.append(prob)
                    tmp_count[min_l] = read_number
                #find whether there is other dominant base and label it
                min_l = np.argmin(tmp_count)
                while (tmp_count[min_l] < max_num):
                    letter = min_l + 1
                    for lab_index in np.where(tmp == letter)[0]:
                        read_lab[lab_index] = 1
                        read_lab_num +=1                        
                    tmp_count[min_l] = read_number
                    min_l = np.argmin(tmp_count)
            else:
                continue


#sparse matrix for storing probability matrix
p_row = np.array(p_row_l)
p_col = np.array(p_col_l)
p_val = np.array(p_val_l)
p_csr = sp.coo_matrix((p_val, (p_row, p_col))).tocsr()

#-----------------error correct---------------------------
Num_errors = 0
Num_oneerr = 0
minus = 0
num_perr = 0
sum_perr =0
aver_perr =0

for i in range(p_csr.shape[0]):
    tmp = p_csr.getrow(i).toarray()[0]
    #whether has the lowest frequency in a read
    #minus stores the number of -0.9
    minus = tmp[tmp<0].shape[0]
    #pretantial number of errors,stores in num_perr
    num_perr = np.nonzero(tmp)[0].shape[0]

    if(num_perr == 0):
        continue
    elif(num_perr == 1):
        #print(i,num_perr,sum_perr)
        #print(aver_perr)
        #correct the first error(which p is -0.9)
        col_index = np.where(tmp!=0)[0][0]
        seq_index = col_index - 50-1
        lil[i,col_index] = Seq_tru[seq_index]
        Num_oneerr += 1
        Num_errors += 1
        # more than 2 base is -0.45
    elif(num_perr == 2 and num_perr < 4):
        #print(i,num_perr)
        #correct error1
        col_index = np.where(tmp!=0)[0][0]
        seq_index = col_index - 50-1
        lil[i,col_index] = Seq_tru[seq_index]
        #correct error2
        col_index = np.where(tmp!=0)[0][1]
        seq_index = col_index - 50-1
        lil[i,col_index] = Seq_tru[seq_index]
        Num_errors += 2  
    sum_perr = 0
    aver_perr = 0
            
print("%d errors are corrected by InsEC" % Num_errors)
#corrected reads in csr
csr = lil.tocsr()

#store .fa
W_file = args.write_file
fa_file = open(W_file,'w')

#pdb.set_trace();
#store seq.txt
W_seq_file = args.write_file_seq
seq_file = open(W_seq_file,'w')
seq_file.write('>the updated nucleotide sequence\n') 
seq_file.write('%s\n' %("".join(Seq_let)))
seq_file.close()

#corrected_read & output

#i record the number of row
i=0
#the number of stored reads
read_sto=0 
for R in csr:
    if(read_lab[i] == 1):
        i += 1
        continue
    else:
        cor_seq = []
        read_ind = np.squeeze(R.toarray())
        read_letter = read_ind[np.nonzero(read_ind)]
        #make sure read_length is 100
        if(read_letter.shape[0] == args.read_L):
            #change number to letter
            for j in read_letter:
                cor_seq.append(dict_tole[j])
                seq_t = ''.join(cor_seq)      
            #transform reads map to opposite strand
            if(r['strand'].loc[i] == 16):
                seq = DNA_Re_complement(seq_t)
            else:
                seq = seq_t
                #write correct read into file
                
            fa_file.write('>%s\n%s\n' %(r['ID'].loc[i],seq))
        else:
            '''
            print(i,"read length are changed falsely")
            '''
        read_sto += 1
        i +=1
        
print("%d bases are updated in output sequence" % int(gene_length-non_update))
print("%d reads are stored in the corrected.fa" % read_sto)
print("error correction finished!")
fa_file.close()