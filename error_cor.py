# update in Nov
import os,random,time,sys

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
    parser.add_argument("--gene_L", type=int)
    parser.add_argument("--read_L", type=int, default=151)
    parser.add_argument("--sam_file", type=str, default='extract.sam')
    parser.add_argument("--write_file", type=str, default='corrected.fa')
    parser.add_argument("--write_file_seq", type=str, default='contigs.fa')
    parser.add_argument("--ref", type=str)
    parser.add_argument("--narrowing", type=str, default="False")
    parser.add_argument("--domin_l", type=float, default=0.5)
    parser.add_argument("--count_l", type=int, default=10)
    parser.add_argument("--match_l", type=float, default=0.975)
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

def correct_base(index, length):
    base_count = [A_count[index],T_count[index],C_count[index],G_count[index],N_count[index],de_count[index]]
    #[i[0] for i in sorted(enumerate(base_count), key=lambda x: x[1])]
    dominants = []
    letters = ["A", "T", "G", "C", "N", "-"]



    if length >= 100:
        for i in range(0,len(base_count)):
            if base_count[i]/sum(base_count) > threshold:
                dominants.append(i)
    else:
        for i  in range(0, len(base_count)):
            if base_count[i] > 3:
                dominants.append(i)
    domin_sum = 0
    for j in dominants:
        domin_sum += base_count[j]
    rand = random.randint(0, domin_sum)
    if len(dominants) == 1:
        base = letters[dominants[0]]
    elif len(dominants) == 2:
        if rand <= base_count[dominants[0]]:
            base = letters[dominants[0]]
        else:
            base = letters[dominants[1]]
    elif len(dominants) == 3:
        if rand <= base_count[dominants[0]]:
            base = letters[dominants[0]]
        elif base_count[dominants[0]] < rand <= base_count[dominants[1]] :
            base = letters[dominants[1]]
        else:
            base = letters[dominants[2]]
    else:
        base = Seq_let[index]
    return dict_tonu[base]

def get_dominant_type(index, length):

    base_count = [A_count[index], T_count[index], C_count[index], G_count[index], N_count[index], de_count[index]]
    # [i[0] for i in sorted(enumerate(base_count), key=lambda x: x[1])]
    dominants = []
    letters = ["A", "T", "G", "C", "N", "-"]

    if length >= 100:
        for i in range(0, len(base_count)):
            if base_count[i] / sum(base_count) > threshold:
                dominants.append(i)
    else:
        for i in range(0, len(base_count)):
            if base_count[i] > 3:
                dominants.append(i)
    return len(dominants)

if __name__=="__main__":
    #-------------------- main function is here-----------------------------------------------------
    start_time = time.time()

    args = get_arguments()
    print('Running parameters:\n')
    print(json.dumps(vars(args), indent=4, separators=(',', ':')))


    R_file = args.sam_file
    '''
    with open(R_file) as rf:
        for line in rf:
            print(line)
       '''

    r = pd.read_csv(R_file, delimiter =' ', names = ['ID','strand','sta_p','sam_q', 'cigar'],encoding = 'unicode_escape')


    dict_tonu = {'A':1,'C':2,'T':3,'G':4,'N':5, '-':6}
    dict_tole = dict(zip(dict_tonu.values(),dict_tonu.keys()))

    #setting
    gene_length = args.gene_L
    read_length = args.read_L
    read_number = r.shape[0]
    count_threshold = args.count_l
    start = 0
    domin_threshold = args.domin_l
    match_limit = args.match_l
    only_contigs = True
    threshold = 0.025


    row_l = []
    col_l = []
    val_l = []

    maxposition = gene_length #debug variable
    maxindex =0 #debug variable
    labindexes = {}#debug
    added_read = {} # read index with added N in the end and the original length
    # read reference genome
    ref = ""
    changed = 0#debug
    with open(args.ref,'r') as refg:
        for line in refg:
            if ">" not in line:
                ref += line.strip()

    #relocated seq to matrix A
    exclude_reads = set({})
    insertion_reads = {}
    real_len = [0] * read_number
    included_i = 0
    narrowed_read = []
    for i in range(read_number):
        exclude=False
        index = int(r['sta_p'].loc[i])-1
        sam_q = list(r['sam_q'].loc[i])
        cigar = str(r['cigar'].loc[i])
        cigar_str = re.findall(r"[0-9]+[MIDSH]",cigar)
        blk_pos = []
        blk_type = []
        ini = 0
        tmp_length = 0

        matched = 0
        for block in cigar_str:
            m = re.search(r'([0-9]+)([MIDSH])', block)
            bl = int(m.group(1)) + ini
            bt = str(m.group(2))
            #if bt == "S" or bt == "H":
            #    continue
            if bt=="M" :#or bt=="I" or bt=="D":
                matched += int(m.group(1))
            blk_type.append(bt)

            blk_pos.append(bl)
            ini = bl

            if bt!="D" and bt!="H":
                tmp_length += int(m.group(1)) # get read length without counting clipped bases

        real_len[i] = tmp_length

        if( matched/tmp_length < match_limit):
        #if(matched/read_length<match_limit):
            exclude_reads.add(i)
            #print(i,str(r.loc[i]), matched)
            exclude = True
            continue

        narrowed_read.append([str(r["ID"].loc[i]), int(r["strand"].loc[i]), index+1, str(r["sam_q"].loc[i]), cigar])
        #print("sam_q:", sam_q, "\n")
        sam_q_num = []

        #softclipping adjust index position
        s_start = 0
        s_end = 0


        c = 0
        inserts = []
        begin = 0
        reduce = 0 #deduct "non-existent" bases in total length, "D" and "H" not shown in reads
        for j in range(0,blk_pos[-1]): #change here to fill the blank with 0?


            if blk_type[c] == "M" :
                try:
                    sam_q_num.append(dict_tonu[sam_q[j-reduce]])
                except:
                    print(j-reduce, len(sam_q),i)
                    exit(1)


            elif blk_type[c] == "I" :
                inserts.append(dict_tonu[sam_q[j-reduce]])
            elif blk_type[c] == "S":
                sam_q_num.append(dict_tonu[sam_q[j - reduce]])
            elif  blk_type[c] == "H":
                sam_q_num.append(0)
            elif blk_type[c] == "D":
                sam_q_num.append(6)
            if blk_type[c] == "H" or blk_type[c] == "D":
                reduce += 1

            if j ==  blk_pos[c] -1 :   #update start and c, put inserts into hashtable

                if blk_type[c] == "I":
                    if i in insertion_reads.keys():

                        newinsert = insertion_reads.get(included_i).copy()
                        newinsert.append((index+begin, inserts))
                        insertion_reads.update({included_i:newinsert.copy()})
                    else:
                        insertion_reads.update({included_i:[(index+begin,inserts.copy())]})
                begin = blk_pos[c]
                inserts = []

                c += 1
                if c == len(blk_type):
                    break

            #else:
             #   sam_q_num.append(0)


        #index = index-s_start


        if blk_type[0]=="S":
            if index-blk_pos[0] < 0:
                start_pos = blk_pos[0] - index
                sam_q_num = sam_q_num[start_pos:]
            else:
                index = index - blk_pos[0]

        if len(sam_q_num) < read_length:
            added_read.update({i:len(sam_q_num)})
            sam_q_num += [0] * (read_length-len(sam_q_num))
            pad = 0
        else:
            pad =len(sam_q_num) -read_length
        #print(r["ID"].loc[i], cigar, index)
        #t = 0
        #for t0 in sam_q_num:
        #    print("("+str(t)+","+str(t0)+")",end=",")
         #   t+=1
        #print()

        val_l.extend(sam_q_num)
        row_tmp = [int(included_i) for n in range(read_length+pad)]
        row_l.extend(row_tmp)
        col_tmp = [n for n in range(index,index+read_length+pad)]
        col_l.extend(col_tmp)
        if len(sam_q_num) != len(row_tmp) or len(sam_q_num) != len(col_tmp):
            print(r["ID"].loc[i], i, index, len(row_tmp), len(col_tmp), len(sam_q_num))
            print(sam_q_num)
            exit(14)
        if (index + len(sam_q)) > maxposition:
            maxposition = index + read_length
            maxindex = index
        included_i += 1
    row = np.array(row_l)
    col = np.array(col_l)
    val = np.array(val_l)

    csc = sp.coo_matrix((val, (row, col))).tocsc() #matrix
    cor_record = {}
    #construct sparse matrix with insertion
    insertion_columns = set({})
    #add_matrix = sp.coo_matrix((read_number, maxposition+len(insertion_reads.keys())),dtype=np.int32).tolil()
    add_matrix = sp.coo_matrix((read_number-len(exclude_reads), maxposition+len(insertion_reads.keys())),dtype=np.int32).tocsc()

    print("csc shape", csc.shape, "add_matrix shape", add_matrix.shape)
    prev_time = time.time()
    print("stage 1 time ", time.time()-start_time)
    print(insertion_reads)

    if args.narrowing == "True":
        ID_count = {}
        print(len(narrowed_read), " 100% M reads")
        mated = []
        for rl in narrowed_read:
            if rl[0] in ID_count.keys():
                ID_count[rl[0]] += 1

            else:
                ID_count[rl[0]] = 1


        true_total_match = []
        for mate_rl in narrowed_read:
            index = mate_rl[2]-1
            if ref[index:index + len(mate_rl[3])] == mate_rl[3]:
                true_total_match.append(mate_rl)
        narrowed_read = true_total_match

        print(len(true_total_match)," truely matched reads")
        for mate_rl in narrowed_read:
            if ID_count[mate_rl[0]] >= 2:
                mated.append(mate_rl)
        narrowed_read = mated
        print(len(mated), "mated reads")
        with open("narrowed_extract.sam", "w+") as nf1:
            for line in narrowed_read:

                nf1.write(line[0]+" "+str(line[1])+" "+str(line[2])+" "+line[3]+" "+line[4]+"\n")
        with open("narrowed_read.fa", "w+") as nrf1:
            for line in narrowed_read:
                nrf1.write(">"+line[0]+"\n")
                nrf1.write(line[3]+"\n")
        #exit()

    print("r is",r, " readnumber ", read_number, "read length ", read_length)#debug
    print("last base at: " + str(maxposition) + "last read start: "+ str(maxindex) + " gene length is: " + str(gene_length))
    print(len(exclude_reads)," reads were excluded because matched bases less than limit", match_limit)



    for i in insertion_reads.keys():
        for j in insertion_reads[i]:
            index = 0
            index1 = 0
            while index < len(j[1]):
                add_matrix[i,j[0]+index] = j[1][index]
                insertion_columns.add(j[0]+index)
                cor_record.update({j[0]+index: "-" + str(index)+"*"+"|"+dict_tole[j[1][index]]})

                index += 1
    #print(sorted(np.unique(np.array(insertion_columns))))

    #------------------determine whether should be inserted----------

    remove_columns = []
    insertion_columns_list = list(insertion_columns)
    insertion_columns_list.sort()
    #print(insertion_columns_list)
    #print(sorted([x for x in cor_record.keys()]))

    for i in insertion_columns_list:
        tmp = np.squeeze(add_matrix.getcol(i).toarray())
        tmp_count = sum(np.bincount(tmp)[1:])

        if tmp_count < count_threshold:
            #print(i,"with", tmp_count,end=", ")

            remove_columns.append(i)


    #remove corresponding insertions in correct records
    #for i in remove_columns:
    #    cor_record.pop(i)

    remove_columns.sort()
    print("remove columns ",remove_columns)
    print(len(insertion_columns_list),len(remove_columns))
    print("phase 3 time ", time.time()-start_time)

    remove_copy = np.array(remove_columns)

    #remove columns with reads less than count_threshold in insertion_columns_list
    remove_set = set(remove_columns)
    #print(remove_copy)

    ir = 0
    n_icl = []
    icl_ori = []
    for ir in insertion_columns_list:
        if ir not in remove_set:


            move = np.where(remove_copy<ir)
            n_icl.append(ir-move[0].shape[0])
            icl_ori.append(ir)
            #print(ir, ir-move[0].shape[0], ir in remove_set, remove_copy[move])
        else:
            del cor_record[ir]
    #print(insertion_columns_list,n_icl)
    insertion_columns_list = n_icl
    #insertion_columns_list = n_icl
    #exit(2)
    '''
    ir = 0
    while(len(remove_columns)>0):
    
        #if ir >= len(remove_columns):
        #    break
    
            if remove_columns[ir] in insertion_columns_list:
                for ii, iv in enumerate(insertion_columns_list):
                    if iv >remove_columns[ir]:
                        insertion_columns_list[ii] = iv -1
                insertion_columns_list.remove(remove_columns[ir])
                for ir1 in range(0,len(remove_columns)):
                    if remove_columns[ir1] > remove_columns[ir]:
                        remove_columns[ir1] -= 1
                update_pairs = []
                #print("icl",insertion_columns_list)
                #print("before",cor_record)
                del cor_record[remove_columns[ir]]
                for ki in cor_record.keys():
                    #print(ki,end=",")
    
                    if ki > remove_columns[ir]:
                        update_pairs.append((ki-1,cor_record[ki]))
                        #kremove.append(ki)
                        #kv = cor_record[ki]
                        #print(ki in cor_record.keys(),cor_record[ki], end=",")
                #print(sorted([x[0] for x in update_pairs]))
                for kp in update_pairs:
    
                    del cor_record[kp[0]+1]
                for kp in update_pairs:
                    cor_record.update({kp[0]:kp[1]})
    
                    #try:
                        #del cor_record[ki]
                    #except:
                    #    print("error")
                    #    print(ir, remove_columns, insertion_columns_list)
                    #    exit(-5)
                        #cor_record.update({ki-1:kv})
    
                #print("ck",sorted([x for x in cor_record.keys()]))
                remove_columns.pop(ir)
    
            else:
                ir += 1
    '''

    #print(cor_record)
    #print(insertion_columns_list)
    print("\nset up insertion column list ",len(insertion_columns_list), insertion_columns_list)
    print("phase 4 time ", time.time()-start_time, "step time ", time.time()-prev_time)
    prev_time = time.time()

    #---------------move inserted columns in insertion matrix

    update_pairs=[]
    #print(remove_copy)
    del i
    for i2,i in enumerate(icl_ori):
        move = np.where(remove_copy<i)
        #print(move)
        move_len = move[0].shape[0]
        add_matrix[:,insertion_columns_list[i2]] = add_matrix[:,i]
        tmp = np.squeeze(add_matrix.getcol(i).toarray())
        tmp_1 = np.squeeze(add_matrix.getcol(insertion_columns_list[i2]).toarray())
        tmp_count = np.bincount(tmp)[1:]
        tmp_1_count = np.bincount(tmp_1)[1:]
        update_pairs.append((insertion_columns_list[i2],cor_record[i]))

        #print(i, move_len,insertion_columns_list[i2], np.array_equal(tmp_count,tmp_1_count),tmp_count, tmp_1_count)


    for up in update_pairs:
        cor_record.update({up[0]:up[1]})


    #combine insertions and other parts of csc
    #exit(-2)

    j0 = 0
    i = 0


    #---------------------determine all inserted columns-----------------
    '''
    con_list = []
    insert_col_length = {}
    while (i < len(insertion_columns_list)):
        prevtail = False
    
        j= i
        prevtail = False
        start = insertion_columns_list[i]
        count = 0
        while(insertion_columns_list[j+1]==insertion_columns_list[j]+1):
            if j+1==len(insertion_columns_list):
                break
            if count == 0:
                insert_col_length.update({start:2})
            else:
    
                insert_col_length.update({start:insert_col_length[start]+1})
            con_list.append(insertion_columns_list[j])
            count += 1
            prevtail = True
            j+=1
    
        if j==len(insertion_columns_list)-1 and insertion_columns_list[j] == insertion_columns_list[j-1]+1:
            con_list.append(insertion_columns_list[j])
            insert_col_length.update({start:count+1})
        elif j == i:
            insert_col_length.update({start: 1})
        i = j
    
      #'''
    if len(insertion_columns_list)>0:
        print("copyting csc to addmatrix if insertions are sure")
        del i
        for i in insertion_columns_list:
            print(i)
            test_matrix= sp.coo_matrix((val, (row, col))).tocsc()
            #print("before insert: ",np.array(test_matrix.getcol(i)))
            col[col >= i] += 1
            tmp = np.squeeze(add_matrix.getcol(i).toarray())
            tmp_base = tmp[np.nonzero(tmp)]

            row = np.append(row,np.array(np.nonzero(tmp)))
            col = np.append(col,np.array([i]*len(tmp_base)))

            val = np.append(val,np.array(tmp_base))


            #print(np.nonzero(tmp))
            #print("nprow",row.shape ,row[-len(tmp_base):])
            #print(np.array(tmp_base))
            #print("npval", val.shape, val[-len(tmp_base):])
            #print(np.array([i]*len(tmp_base)))
            #print("npcol", col.shape, col[-len(tmp_base):])
            #print( np.array(np.nonzero(tmp)).size,len(np.array([i] * len(tmp_base))),len(tmp_base))



            tmp_val=np.array(tmp_base)
            tmp_row = np.array(np.nonzero(tmp)[0])
            tmp_col = np.array([i]*len(tmp_base))
            #print(tmp_val.shape,tmp_row.shape,tmp_col.shape)
            test_matrix= sp.coo_matrix((val, (row, col))).tocsc()
            #print(np.array(test_matrix.getcol(i)))
            #print(np.array(test_matrix.getcol(i+1)))
            #exit(1)
        com_matrix = sp.coo_matrix((val, (row, col))).tocsc()
        test_icl = np.array(insertion_columns_list)
        #------------make sure copy process is correct---------

        for i in range(0,com_matrix.shape[1]):
            #print(i)
            if i not in insertion_columns_list:
                tmp = com_matrix.getcol(i).toarray()

                tmp1 = csc.getcol(i-len(test_icl[test_icl<=i])).toarray()
                if np.sum(tmp!=tmp1)!=0:
                    print(i, i-len(test_icl[test_icl<=i]))
                    print("com", tmp[np.nonzero(tmp)])
                    #print(np.nonzero(tmp))
                    print("csc",tmp1[np.nonzero(tmp1)])
                    exit(5)
                #print(np.nonzero(tmp1))

        add_matrix = com_matrix
        '''
        #---------------make sure forming matrix is equal to copy column---------
        for j1 in range(0,add_matrix.shape[1]):
    
            if j1 in insertion_columns_list:
                #tmp = np.squeeze(add_matrix.getcol(j1).toarray())
                #tmp_count = np.bincount(tmp)[1:]
                #print(j1, "with", sum(tmp_count), end=", ")
    
                continue
            if j0 >= csc.shape[1]:
                break
            #tmp = np.squeeze(csc.getcol(j0).toarray())
            #tmp_count = np.bincount(tmp)[1:]
            #if tmp_count.shape[0] == 0:
             #   j0 += 1
             #   continue
            #try:
            #if j1 % 100 == 0:
            #    print(time.time()-start_time)
            #print(j0,csc.getcol(j0).shape)
    
            add_matrix[:,j1] = csc[:,j0]
            #check if the copy process is correct
            add_arr = add_matrix.getcol(j1).toarray()
            csc_arr = csc.getcol(j0).toarray()
            com_arr = com_matrix.getcol(j1).toarray()
            j0 += 1
            if not np.array_equiv(add_arr,csc_arr) or not np.array_equiv(com_arr,csc_arr):
    
                print(j1,j0)
                print((add_arr ==csc_arr).all(),(com_arr==csc_arr).all())
                print("add",add_arr)
                print("csc",csc_arr)
                print("com",com_arr)
                exit(1)
    
    
    
        del i
        diff_col = []
        with open("copy_log.txt","w+") as cf:
            for i in range(0,add_matrix.shape[1]):
                #print(i)
                add_arr = add_matrix.getcol(i).toarray()
                tmpadd = add_arr[add_arr>0]
                if i >= com_matrix.shape[1] :
                    print("add still ", np.nonzero(add_arr))
                    if len(tmpadd) > 0:
                        exit(5)
                    else:
                        continue
    
                com_arr = com_matrix.getcol(i).toarray()
    
                if not np.array_equiv(add_arr, com_arr):
                    diff_col.append(i)
                    cf.write(str(i)+"\n")
                    cf.write("com"+str(com_arr[np.nonzero(com_arr)])+"\n")
                    cf.write("add"+str(add_arr[np.nonzero(add_arr)])+"\n")
        print(len(diff_col), diff_col)
        print(add_matrix.shape,com_matrix.shape)
                #exit(2)
        exit(1)
        '''


        print("copied csc to add_matrix", add_matrix.shape)
    else:
        print("no insertion")
        add_matrix = csc
    print("phase 5 time ",time.time()-start_time,"step ", time.time()-prev_time)
    prev_time = time.time()
    #find difference between contigs (without insertion) and reference
    newseq = []
    i = start
    for i in range(start, csc.shape[1]):

        tmp = np.squeeze(csc.getcol(i).toarray())#matrix

        tmp_count = np.bincount(tmp)[1:]

        #col_num store the value of total_coverage_col
        col_num = sum(tmp)
        #maximum letter's number is max_l + 1
        #print(i, tmp)


        if(sum(tmp_count) == 0):

            if i > len(ref)-1:
                continue
            newseq.append(ref[i-start])

            continue
        else:

            nonzero_num = np.nonzero(tmp_count)[0].shape[0]
            max_l = np.argmax(tmp_count)
            max_num = tmp_count[max_l]

            if i > len(ref)-1:
                newseq.append(dict_tole[max_l + 1])
                cor_record.update({i: "*|" + dict_tole[max_l + 1]})
                continue


            ref_le = ref[i-start]
            ref_pos = i -start
            if dict_tole[max_l+1] != ref[i-start]:
                #print(tmp_count)
                if sum(tmp_count) >= count_threshold:
                    newseq.append(dict_tole[max_l+1])
                    if max_l == 6:
                        cor_record.update({i: ref[i - start] + "|*"})
                    else:
                        cor_record.update({i: ref[i - start] + "|" + dict_tole[max_l + 1]})
                else:

                    newseq.append(ref[i-start])
            else:
                newseq.append(dict_tole[max_l + 1])

    print("cor_record set", cor_record)
    print("phase 6 time ",time.time()-start_time," step ", time.time()-prev_time)
    prev_time = time.time()
    #lil = sp.coo_matrix((val, (row, col))).tolil() #matrix
    lil = add_matrix.copy()
    #exit(2)
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


    A_count = [0] * add_matrix.shape[1]
    C_count = [0] * add_matrix.shape[1]
    T_count = [0] * add_matrix.shape[1]
    G_count = [0] * add_matrix.shape[1]
    N_count = [0] * add_matrix.shape[1]
    de_count = [0] * add_matrix.shape[1]



    #-------------------------compute probabilities of each base
    #count unknow bases in the updated sequence
    non_update=0
    insufficient = []
    emptycol = []
    sufficient = []
    domin_values  = {}
    for x1,x1v in enumerate(real_len):
        domin_values.update({x1:[1]*x1v})

    i=start
    i1 = start #i1 is cursor for ref
    ins_count = 0
    for i in range(start, add_matrix.shape[1]):
    #for i in range(1,gene_length+1):
        # go through each columntm

        try:
            tmp = np.squeeze(add_matrix.getcol(i).toarray())#matrix
        except:
            print(i)
            exit(2)
            continue
        tmp_count = np.bincount(tmp)[1:]

        #col_num store the value of total_coverage_col
        col_num = sum(tmp_count)
        #maximum letter's number is max_l + 1

        if(sum(tmp_count) == 0):
            #0 means this base is unknown
            # use original ref sequence if unknown
            emptycol.append(i)

            if i < len(ref) + len(insertion_columns_list) -2:

                if ins_count < len(insertion_columns_list):
                    if i <= insertion_columns_list[ins_count]:
                        ins_num = ins_count
                    else:
                        ins_count += 1
                        ins_num = ins_count
                else:
                    ins_num=len(insertion_columns_list)

                if(i in insertion_columns_list):
                    Seq_let.append(dict_tole[5])
                    Seq_tru.append(5)
                    print("shouldn't insert? ",i,insertion_columns_list)
                    print(i1)
                    print(tmp)
                    print(tmp_count)
                else:
                    non_update = non_update+1

                    Seq_tru.append(dict_tonu[ref[i1-start-ins_num]])
                    Seq_let.append(ref[i1-start-ins_num])
                    i1 += 1
            else:
                Seq_let.append(dict_tole[5])
                Seq_tru.append(5)
            continue
        else:

            if tmp_count.shape[0] >= 1:
                A_count[i-start] = tmp_count[0]
            if tmp_count.shape[0] >= 2:
                C_count[i-start] = tmp_count[1]
            if tmp_count.shape[0] >= 3:
                T_count[i-start] = tmp_count[2]
            if tmp_count.shape[0] >= 4:
                G_count[i-start] = tmp_count[3]
            if tmp_count.shape[0] >= 5:
                N_count[i-start] = tmp_count[4]
            if tmp_count.shape[0] >= 6:
                de_count[i-start] = tmp_count[5]

            nonzero_num = np.nonzero(tmp_count)[0].shape[0]
            max_l = np.argmax(tmp_count)
            max_num = tmp_count[max_l]
            #print(tmp)
            #print(tmp_count)
            #print(i, max_l, max_num, i in insertion_columns_list)

            #if(3175<=i<=3180):
            #    print(i, tmp_count, ref[i], dict_tole[max_l+1], max_l)
            #    print(tmp[tmp == max_l])

            if i >= len(ref)+len(insertion_columns_list)-2:
                Seq_tru.append(max_l + 1)
                Seq_let.append(dict_tole[max_l + 1])

            else:
                #update letter of sequence


                if i in insertion_columns_list:

                    Seq_tru.append(max_l + 1)
                    Seq_let.append(dict_tole[max_l + 1])
                else:

                    try:
                        ref_le = ref[i1-start]
                    except:
                        print("ref index error")
                        print(i1, len(insertion_columns_list))
                        print(i,len(ref))
                        exit(1)
                    ref_pos = i1 -start
                    if dict_tole[max_l+1] != ref[i1-start]:
                        #print(tmp_count, dict_tole[max_l+1], ref[i-start])
                        if sum(tmp_count) >= count_threshold:

                            Seq_tru.append(max_l+1)
                            Seq_let.append(dict_tole[max_l+1])
                        else:
                            non_update += 1
                            insufficient.append((i,sum(tmp_count)))
                            Seq_tru.append(dict_tonu[ref[i1-start]])
                            Seq_let.append(ref[i1-start])
                    else:
                        non_update += 1
                        sufficient.append(i)
                        Seq_tru.append(max_l + 1)
                        Seq_let.append(dict_tole[max_l + 1])

                    i1 += 1
            if only_contigs:
                continue
            #whether there are errors in the column
            if(nonzero_num == 1):
                #no error in this column, only one dominant base
                 continue
            else:


                tmp_count[tmp_count == 0] = read_number
                #if(col_num < 100):
                #according to base count to decide the dominant base
                #finding errors and mark with prob
                err_count = [0] * tmp.shape[0]


                if sum(tmp_count) < 100:
                    if(np.min(tmp_count) < 3):

                        min_l = np.argmin(tmp_count)
                        letter = min_l + 1
                        for k in np.where(tmp == letter)[0]:
                            p_row_l.append(k)
                            p_col_l.append(i)
                            #prob = tmp_count[min_l]/max_num * tmp_count[min_l]/col_num
                            prob = min_l/col_num + 0.1*get_dominant_type(i,real_len[i])
                            #domin_values.update({k:prob})

                            p_val_l.append(prob)
                            err_count[k] += 1
                        tmp_count[min_l] = read_number
                    else:
                        continue
                    #find whether there is other dominant base and label it
                    min_l = np.argmin(tmp_count)

                    while (tmp_count[min_l] < 3):
                        letter = min_l + 1

                        for lab_index in np.where(tmp == letter)[0]:
                            err_count[lab_index] += 1
                            p_row_l.append(lab_index)
                            p_col_l.append(i)
                            #prob = tmp_count[min_l]/max_num * tmp_count[min_l]/col_num
                            prob = min_l/col_num + 0.1*get_dominant_type(i,real_len[i])
                            p_val_l.append(prob)
                            #if err_count[lab_index] / len(str(r['sam_q'].loc[lab_index])) > threshold:
                            #    read_lab[lab_index] = 1
                            #    read_lab_num += 1

                                #labindexes.append(lab_index)

                        tmp_count[min_l] = read_number
                        min_l = np.argmin(tmp_count)
                else:
                    if (np.min(tmp_count)/col_num < threshold):
                        min_l = np.argmin(tmp_count)
                        letter = min_l + 1
                        for k in np.where(tmp == letter)[0]:
                            p_row_l.append(k)
                            p_col_l.append(i)
                            prob = min_l/col_num + 0.1*get_dominant_type(i,real_len[i])
                            p_val_l.append(prob)
                            err_count[k] += 1
                        tmp_count[min_l] = read_number
                    else:
                        continue
                        # find whether there is other dominant base and label it
                    min_l = np.argmin(tmp_count)

                    while (tmp_count[min_l]/col_num < threshold):
                        letter = min_l + 1

                        for lab_index in np.where(tmp == letter)[0]:
                            err_count[lab_index] += 1
                            p_row_l.append(lab_index)
                            p_col_l.append(i)
                            prob = min_l / col_num + 0.1 * get_dominant_type(i, real_len[i])
                            p_val_l.append(prob)
                            #if err_count[lab_index] / len(str(r['sam_q'].loc[lab_index])) > threshold:
                            #    read_lab[lab_index] = 1
                            #    read_lab_num += 1
                            # labindexes.append(lab_index)

                        tmp_count[min_l] = read_number
                        min_l = np.argmin(tmp_count)

    print("majortiy of contigs acquired ")
    print("time is ", time.time()-start_time,"step ",time.time()-prev_time)
    prev_time = time.time()
    #remove deleted sequence
    Seq_let = list(filter(lambda x: x != "-", Seq_let))
    seq_length = len(Seq_let)
    let = gene_length
    while let < add_matrix.shape[1] :
        if let >= seq_length:
            break
        tmp = np.squeeze(add_matrix.getcol(let).toarray())
        tmp_count = np.bincount(tmp)[1:]
        if let > len(ref)+len(insertion_columns_list) -2:

            Seq_let.pop(let)
            seq_length -= 1

        else:
            let += 1
    if add_matrix.shape[1] < gene_length or len(Seq_let) < gene_length:
        unknown_tail_len = gene_length-len(Seq_let)
        Seq_let.extend(ref[len(Seq_let)+len(insertion_columns_list):])

        print(unknown_tail_len,"bases at the end has no covering reads, use reference sequence's bases from",len(Seq_let)+len(insertion_columns_list) )

    print("got the contig , time ", time.time()-start_time, " step ", time.time()-prev_time, "ignore reads? ",only_contigs)
    print(non_update," bases remained the same as reference, ", len(insufficient)+len(emptycol)," positions not updated becasue of insufficient information")
    print(insufficient)
    print(len(sufficient))

    print("empty positions:",emptycol)
    print(len(emptycol))
    #corrected reads in csr
    csr = lil.tocsr()#matrix



    #pdb.set_trace();
    #store contigs.fa
    W_seq_file = args.write_file_seq
    seq_file = open(W_seq_file,'w')
    seq_file.write('>updated_sequence\n')
    seq_file.write('%s\n' %("".join(Seq_let)))
    seq_file.close()

    print("error correction finished!")
    print("correction records are in records.txt")
    #print(Seq_let)
    with open("records.txt",'w') as rf:
        c = 1
        sorted_record = sorted(cor_record.keys())
        for k in sorted_record:
            if k > len(ref)+len(insertion_columns_list) :
                continue
            eol = "\n" if c % 15 ==0 else "   "
            rf.write(str(k)+": "+cor_record.get(k)+eol)
            c+=1
    print("%d bases are updated in output sequence" % int(gene_length-non_update))
    #corrected_contig
    if only_contigs:
        exit(0)

    # -----------------read error correct---------------------------
    #sparse matrix for storing probability matrix
    p_row = np.array(p_row_l)
    p_col = np.array(p_col_l)
    p_val = np.array(p_val_l)

    p_csr = sp.coo_matrix((p_val, (p_row, p_col))).tocsr()#matrix
    #corrected.fa
    W_file = args.write_file
    fa_file = open(W_file,'w')
    Num_errors = 0
    Num_oneerr = 0
    minus = 0
    num_perr = 0
    sum_perr = 0
    aver_perr = 0

    print("labling time start")
    prev_time = time.time()
    i = 0
    for i in range(p_csr.shape[0]):  # matrix
        tmp = p_csr.getrow(i).toarray()[0]  # matrix
        # whether has the lowest frequency in a read
        # minus stores the number of -0.9
        minus = tmp[tmp < 0].shape[0]
        # pretantial number of errors,stores in num_perr
        num_perr = np.nonzero(tmp)[0].shape[0]

        if (num_perr == 0):
            continue
        else:
            count = 0

            np.ndarray.sort(tmp)
            if real_len[i] >= 100:
                err_est = int(real_len[i] * threshold)
            else:
                err_est = 3
            if num_perr >= err_est:
                nonzero_index = 0
                for t in tmp:
                    if t > 0:
                        nonzero_index += 1
                    if nonzero_index == err_est:
                        if t < domin_threshold:
                            read_lab[i] = 1
                            read_lab_num += 1
                        break
                    nonzero_index += 1
            if (read_lab[i] != 1):
                col_index = np.where(tmp != 0)[0][0]
                length = real_len[i]
                lil[i, col_index] = correct_base(col_index, length)
                Num_errors += num_perr

        sum_perr = 0
        aver_perr = 0

    print("finish labelling time", time.time()-start_time, "Step ", time.time()-prev_time)
    prev_time = time.time()
    #i record the number of row
    i=0
    #the number of stored reads
    read_sto=0
    for R in csr:
        if(read_lab[i] == 1 or i in exclude_reads):
            i += 1
            continue
        else:
            cor_seq = []
            read_ind = np.squeeze(R.toarray())
            read_letter = read_ind[np.nonzero(read_ind)]
            #make sure read_length is 100
            #if(read_letter.shape[0] == args.read_L):
                #change number to letter
            for j in read_letter:
                cor_seq.append(dict_tole[j])
                seq_t = ''.join(cor_seq)
            #transform reads map to opposite strand
            if(r['strand'].loc[i] & 16 != 0):
                seq = DNA_Re_complement(seq_t) #reverse complement of bwa
            else:
                seq = seq_t
                #write correct read into file
            #if i in added_read.keys():#delete added 0s at the end
            #    seq = seq[:added_read.get(i)]

            fa_file.write('>%s\n%s\n' %(r['ID'].loc[i],seq))
            #else:
            '''
            print(i,"read length are changed falsely")
            '''
            read_sto += 1
            i +=1

    fa_file.close()
    print("finish correcting reads", time.time()-start_time, "step ", time.time()-prev_time)

    print("%d reads are stored in the corrected.fa" % read_sto)
    print("%d errors are corrected by InsEC" % Num_errors)
    '''
    with open("difference.txt","w+") as dif:
        c = 0
        for i in range(0,len(ref)):
            if Seq_let[i] != ref[i]:
                eol = "\n" if c % 15 == 0 else "   "
                dif.write(str(i) + ": " + ref[i] + "|"+Seq_let[i] + eol)
                c += 1
    '''