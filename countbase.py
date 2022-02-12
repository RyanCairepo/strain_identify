import os
import random
import subprocess
import sys, re, statistics, itertools as its, multiprocessing as mp, numpy as np, scipy.sparse as sp, matplotlib.pyplot as pit
from typing import Counter
import  argparse,math
import editdistance
import matplotlib.pyplot as plt


def countbase():
    baselength = {}
    with open(sys.argv[1]) as f:
        for line in f:
            sl = line.split(" ")
            if len(sl) > 1:

                length = int(re.findall(r'[0-9]+', sl[2])[0])

                if(length in baselength.keys()):
                    baselength.update({length:baselength.get(length) + 1})
                else:
                    baselength.update({length:1})
    sortedlength = {k : v for k, v in sorted(baselength.items(), key=lambda item : item[1])}
    print(sortedlength)
def countmatch():
    count = 0
    lengths ={}
    with open("gene2.match", "w+") as fw:
        with open("SRR15096407.sam", "r+") as f:
            for line in f:
                sl = re.split(' |\t', line)
                #print(len(sl))
                if(len(sl)>3 and sl[5] != "*"):
                    fw.write(str(sl[5])+" "+ str(len(sl[9]))+ "\n")
                    if len(sl[9]) not in lengths.keys():
                        lengths.update({len(sl[9]):1})
                    else:
                        lengths.update({len(sl[9]):lengths.get(len(sl[9]))+1})
                    count+=1

    print(count)
    sortedlength = {k: v for k, v in sorted(lengths.items(), key=lambda item: item[1])}
    print(sortedlength)

def combine():
    os.chdir("./sra_downloads")
    with open("../combinedsra.fastq","a+") as wf:
        for file in os.listdir(os.getcwd()):
            with open(file,"r") as rf:
                #print(file)
                if re.findall(r'SRR11', file):

                    for line in rf.readlines():
                        print(line)
                        wf.write(line)

def countdiff():
    with open(sys.argv[1], "r") as f1:

        seq1 = f1.read()
        if seq1[0] == ">":
            f1.seek(0)
            lines = f1.readlines()[1:]
            seq1 = ""
            for l in lines:
                seq1 += l.replace("\n","").replace(" ","")
        else:
            seq1 = seq1.replace("\n","").replace(" ","")
    with open(sys.argv[2], "r") as f2:
        seq2 = f2.read()
        if seq2[0] == ">":
            f2.seek(0)
            lines = f2.readlines()[1:]
            seq2 = ""
            for l in lines:
                seq2 += l.replace("\n", "").replace(" ", "")
        else:
            seq2 = seq2.replace("\n", "").replace(" ", "")
    #print(len(seq1), len(seq2))
    diffnum = 0
    diff = {}
    if (len(seq1) < len(seq2)):
        for i in range(len(seq2)):
            if i < len(seq1):
                if seq1[i] != seq2[i]:
                    #print("#"+str(i)+": "+seq1[i]+"|"+seq2[i])
                    diffnum += 1
                    diff.update({i: seq1[i] + "|" + seq2[i]})
            else:
                diffnum+=1
                diff.update({i: "*|" + seq2[i]})
                #print("#"+str(i)+" seq2: "+seq2[i])

    else:
        for i in range(len(seq1)):
            if i < len(seq2):
                if seq1[i] != seq2[i]:
                    #print("#"+str(i)+": "+seq1[i]+"|"+seq2[i])
                    diffnum+=1
                    diff.update({i: seq1[i] + "|" + seq2[i]})
            else:
                diffnum+=0
                diff.update({i: seq1[i] + "|*" })
                #print("#"+str(i)+" seq1: "+seq1[i])
    #print(diffnum)
    if diffnum > 0:
        s_diff = sorted(diff.keys())
        c = 0
        for k in s_diff:
            eol = "\n" if c % 15 == 0 and c!= 0 else "   "
            #print(str(k) + ": " + diff.get(k), end=eol)
            c+=1
    print(editdistance.eval(seq1,seq2))

#def splitsequence():
def oneline_fasta():
    with open(sys.argv[1], "r") as f:
        seq = ""
        for line in f:
            if line[0] != ">":

                temp =   line.replace(" ","").replace("\n","").replace("\d","").upper()
                seq += re.sub('\d', "", temp)

        with open("oneline_"+sys.argv[1].replace("..","").replace("/",""), "w+") as fw:
            fw.write(seq)

def count_record():
    record1 = {}
    with open(sys.argv[1],"r") as f1:
        for line in f1:
            matches = re.findall('[0-9]+: [\-]*[A-Z*]+\|[A-Z*]+', line)
            for m in matches:
                parts = re.search('([0-9]+): ([\- ]*[A-Z*]+\|[A-Z*]+)', m)
                index = int(parts.group(1))
                cor = parts.group(2)
                record1.update({index:cor})
    record2 = {}
    with open(sys.argv[2],"r") as f2:
        for line in f2:
            matches = re.findall('[0-9]+: [\-]*[A-Z*]+\|[A-Z*]+', line)
            for m in matches:
                parts = re.search('([0-9]+): ([\- ]*[A-Z*]+\|[A-Z*]+)', m)
                index = int(parts.group(1))
                cor = parts.group(2)
                record2.update({index:cor})
    count = 0
    for k in record1.keys():
        if k in record2.keys():
            if record1[k] == record2[k]:

                count+=1
    print(count," corrections are the same")
def count_match(): #matched bases
    readlist = []
    matches = []
    tmp = 0
    narrowed =[]
    count = 0

    with open(sys.argv[1],"r") as f:

        for line in f:

            s =line.split(" ")
            readlist.append(s)
            m = 0

            if len(s) < 3:
                continue
            #print("her?")
            '''
            group = re.findall('([0-9]+)([MIDSH])',s[4])

            tmp = 0
            for i in group:
                #if i[1] != "D" or i[1]!="H":
                tmp += int(i[0])
                if i[1] == "M":
                    m += int(i[0])

            if m==0 or m/tmp < 1 :
                continue
            '''
            count += 1
            matches.append(m)

            narrowed.append(s)
            #print(line)


    ID_count = {}
    for read in narrowed:
        if read[0] not in ID_count.keys():
            ID_count.update({read[0]:1})
        else:
            ID_count[read[0]] += 1
    strand_count = {}
    single = []
    for id in ID_count.keys():
        if ID_count[id] < 2:
            single.append(id)
    print(len(readlist),len(narrowed), count, len(single))
    #print(single)
    for read in narrowed:

        if read[0] in single:
            continue
        #if not re.match('^[0-9]+[M]$', read[4]):
        if int(read[1]) not in strand_count.keys():
            strand_count.update({int(read[1]):1})
        else:
            strand_count[int(read[1])] += 1
    #print(strand_count)
    print({k: v for k, v in sorted(strand_count.items(), key=lambda item: item[0])})
    field_count = {}
    for strand in strand_count.keys():
        bin = "{0:b}".format(int(strand))

        i = 0
        for i in range(0,len(bin)):
            if bin[i] == "1":
                if i not in field_count:
                    field_count[i] = strand_count[strand]
                else:
                    field_count[i] += strand_count[strand]
    #print(field_count)
    print({k:v for k,v in sorted(field_count.items(),key=lambda item: item[0])})


    #with open("raw_high_match_extract.sam","w+") as f:
    #    for n in narrowed:
    #        f.write(line)


    print(len(matches),"reads met standard")
    print("average matches: ", statistics.mean(matches))
    print("median match: ",statistics.median(matches))
    matches.sort()
    print("30% ", matches[int(len(matches) * 0.3)], "10%", matches[int(len(matches) * 0.1)])


def hash_str():
    i = 0
    #hexp = hash(exp[0])
    #----------------dup_count is a dict, uses read as the key and a list
    # [frequency, starting position, cigar, match] as value----------------------
    with open(sys.argv[1],"r") as f:
        dup_count = {}
        for line in f:

            sl = line.split(" ")
            sl[4] = re.sub('\n','',sl[4])
            read_id = sl[0]

            read = sl[3]
            match_number = cigar_processing(sl[4])
            #if match_number < 145:# or match_number/len(sl[3]) < 0.97:
            if len(read) < 145:
                continue
            #hr = hash(read)
            #duplicate read found
            if(read in dup_count.keys()):
                #value_list = dup_count[read]
                #see if starting index is different
                #if value_list[1] != int(sl[2]):

                '''
                count = 1
                fre_count = 0
                index_exist = False
                while (read+str(count) in dup_count.keys()):
                    #fre_count += dup_count[read+str(count)][0]
                    #dup_count[read+str(count)][0] += 1
                    if int(sl[2]) == dup_count[read+str(count)][1]:
                        dup_count[read+str(count)][0] += 1
                        index_exist = True
                        break
                    else:
                        count+=1
                if not index_exist:
                    dup_count.update({read+str(count):[1,int(sl[2]), sl[4]]})
                '''
                #else:
                if dup_count[read][3] < match_number:
                    dup_count[read][1] = int(sl[2])
                    dup_count[read][2] = sl[4]
                    dup_count[read][3] = match_number

                dup_count[read][0] += 1

                    #dup_count.update({read:[dup_count[read][0]+1, dup_count[read][1],sl[4]]})
            else:
                #[frequency,index]
                dup_count.update({read:[1,int(sl[2]),sl[4],match_number]})
    '''
    # count neighbors
    for dr in dup_count.keys():
        count = 0
        for dr1 in dup_count.keys():
            if editdistance.eval(dr,dr1) == 1:
                count += dup_count[dr1][0]
        print(count,end=", ")
        dup_count[dr].append(count)
    '''
    #sort here by frequency
    sorted_distinct = {k:v for k,v in sorted(dup_count.items(),reverse=True,key=lambda x : x[1][0])}
    with open("distinct_report.txt", "w+") as dwf:
        for k,v in sorted_distinct.items():
            dwf.write(k+" "+str(v[0])+"\n")
            #print(k, "f:", v[0], "p:", v[1], "c:", v[2])#, "n", v[4])

    print("fasta count duplicate")
    f = open("extract.sam")
    lines = f.readlines()
    with open('extract.fasta', 'w+') as ff:
        for line in lines:
            name_seq = line.split(' ')
            ff.write(name_seq[0])
            ff.write('\n')
            ff.write(name_seq[3])
            ff.write('\n')
    f.close()
    seqs1 = []

    len_seqs1_dict = {}
    with open("extract.fasta") as f1:
        lines1 = f1.readlines()

        print('lines: ', len(lines1))
        for i in range(len(lines1)):
            if i % 2 != 0:
                seqs1.append(lines1[i])

                seq1_len = len(lines1[i])

                # if seq1_len in len_seqs1_dict.keys():
                len_seqs1_dict.setdefault(seq1_len, []).append(lines1[i])

            else:
                continue

    seq1_dict = Counter(seqs1)
    #print(seq1_dict)
    sorted_seq1_dict = {k:v for k,v in sorted(seq1_dict.items(), reverse=True, key=lambda x:x[1])}
    with open("pengyao_ferq.txt","w+") as pf, open("compare.txt","w+") as wf:
        if sorted_seq1_dict != sorted_distinct:
            for k in sorted_seq1_dict.keys():
                k1 = re.sub('\n','',k)
                pf.write(k+" "+str(sorted_seq1_dict[k])+"\n")
                if len(k1) >=145 and (k1 not in sorted_distinct or sorted_seq1_dict[k] != sorted_distinct[k1][0]):
                    if k1 in sorted_distinct.keys():
                        rv = sorted_distinct[k1]
                    else:
                        rv = "none"

                    wf.write("p: "+k1+" "+str(sorted_seq1_dict[k])+"\n"+"r:"+str(rv)+"\n")
            for k2 in sorted_distinct.keys():
                if k2 not in sorted_seq1_dict or (sorted_distinct[k2][0] -sorted_seq1_dict[k2+"\n"])!=0:
                    if k2+"\n" in sorted_seq1_dict.keys():
                        rv = str(sorted_seq1_dict[k2+"\n"])+"!="+str(sorted_distinct[k2][0])
                        print( (sorted_distinct[k2][0] -sorted_seq1_dict[k2+"\n"]),type(sorted_distinct[k2][0]), type(sorted_seq1_dict[k2+"\n"]))

                    else:
                        rv = "none"

                    wf.write("r: "+k2+" "+str(sorted_distinct[k2][0])+"\n"+"p:"+str(rv)+"\n")

    '''
    print_list = []
    for k in dup_count.keys():
        printed = set({})
        if dup_count[k][0] > 1 :
            if k not in printed:
                #print(k, dup_count[k])
                print_list.append(k+" "+str(dup_count[k]))
                printed.add(k)
        elif re.search('[0-9]',k):

            #print("different index")
            #print(k,dup_count[k])
            print_list.append("different index")
            print_list.append(k+" "+str(dup_count[k]))
            ori = re.sub('[0-9]', '', k)

            if dup_count[ori][0]  <2 :


                if dup_count[ori][1] != dup_count[k][1]:
                    #print(ori, dup_count[ori])

                    #printed.add(k)
                    if ori in printed:
                        print_list.remove(ori+" "+str(dup_count[ori]))
                        print_list.append(ori + " " + str(dup_count[ori]))
                    else:
                        printed.add(ori)
                        print_list.append(ori + " " + str(dup_count[ori]))
    for i in print_list:
        print(i)
    '''

    #print(dup_count)
            #print(read)

def cigar_processing(cigar):
    cigar_str = re.findall(r"[0-9]+[MIDSH]", cigar)
    blk_pos = []
    blk_type = []
    ini = 0
    tmp_length = 0

    matched = 0
    for block in cigar_str:
        m = re.search(r'([0-9]+)([MIDSH])', block)
        bl = int(m.group(1)) + ini
        bt = str(m.group(2))

        if bt == "M":  # or bt=="I" or bt=="D":
            matched += int(m.group(1))
        blk_type.append(bt)

        blk_pos.append(bl)
        ini = bl

        tmp_length += int(m.group(1))  # get read length without counting clipped bases

    return matched

def count_length():
    len_count = {}
    with open(sys.argv[1],"r") as f:
        for line in f:
            sl = line.split(" ")
            read_len = len(sl[3])
            if read_len in len_count.keys():
                len_count[read_len] += 1
            else:
                len_count[read_len] = 1
    s_count={k:v for k,v in sorted(len_count.items(), key=lambda i : i[1], reverse=True)}
    print(s_count)
    x = np.asarray([k for k in s_count.keys()])
    y = np.asarray([v for k, v in s_count.items()])
    px = 1 / plt.rcParams['figure.dpi']
    fig,ax = plt.subplots(figsize=(1900*px, 600*px))
    print(y)
    print(x)
    plt.bar(x,y, align='center')
    plt.xticks(x)
    plt.xlabel('read lengths')
    plt.ylabel('frequency')
    ax.set(xlim=[100,156])
    #plt.xticks([x for x in range(100,153)])
    #ax.xaxis.set_tick_params(width=5)
    #fig, axs = plt.subplots(1,2, sharey=True, tight_layout=True)


    print(fig)
    #axs[0].hist(x, bins=20)
    #axs[1].hist(y,bins=20)
    plt.show()

def getoverlap(f1,f2):
    read_set1 = {}
    overlapped = {}
    read_set2 = {}
    with open(f1, "r") as f1:
        for line in f1:
            sl = line.split(" ")
            if sl[0] not in read_set1.keys():
                read_set1.update({sl[0]: {sl[3]}})
            else:
                read_set1[sl[0]].add(sl[3])
    with open(f2, "r") as f2:
        for line in f2:
            sl = line.split(" ")
            if sl[0] not in read_set2.keys():
                read_set2.update({sl[0]: {sl[3]}})
            else:
                read_set2[sl[0]].add(sl[3])
    rs2_only = []
    for rid in read_set1.keys():

        if rid in read_set2.keys():
            if len(read_set1[rid]) >= len(read_set2[rid]):
                if len(read_set1[rid].intersection(read_set2[rid])) == len(read_set1[rid]) :
                    overlapped.update({rid:read_set1[rid].intersection(read_set2[rid])})
                else:
                    print(read_set1[rid],read_set2[rid])
            else:

                rs2_only.extend(list(read_set2[rid]))


    print(len(overlapped.keys()), "reads ID overlapped ", len(read_set1.keys()), " unique read ID in ", f1, " and ",len(read_set2.keys()),
          " in ", f2)
    overlap_ids=set(read_set1.keys()).intersection(set(read_set2.keys()))
    print(len(read_set1),len(read_set2))
    print(len(overlap_ids))
    print(len(rs2_only),"reads only in "+f2.name)
    exit()
    with open("overlap_ids.txt","w+")as wf:
        for ri in overlap_ids:
            wf.write(ri+",")
    #for k in overlapped.keys():
    #    print()
def get_sam_dup():
    read_set = set({})
    count = 0
    line_count = 0
    with open(sys.argv[1], "r") as f:
        for line in f:
            line_count += 1
            id = line.split(" ")[0]
            if id in read_set :
                count += 1
            else:
                read_set.add(id)
    print(len(read_set)," unique reads with ", count, "duplicates in total ", line_count )

def adjustindex(moveby):
    allreads = []
    with open(sys.argv[1], "r") as f:
        for line in f:
            fields = line.split(" ")
            allreads.append(fields.copy())
    with open("after_"+sys.argv[1],"w+") as f1:
        for i in allreads:
            c= 0
            for j in i:
                if c< 4 and c != 2:
                    f1.write(j+" ")
                elif c == 2:
                    f1.write(str(int(j)-moveby)+" ")
                elif c== 4:
                    f1.write(j)
                c += 1
            #f1.write("\n")

def sep_reads():
    with open(sys.argv[1], "r") as f1:
        for line in f1:
            fiels = line.split(" ")
            if "I" in fiels[4] or "D" in fiels[4]:
                print(line, end="")

def seq_depth(narrow=False):
    dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
    dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))
    unlikely=[]
    with open(sys.argv[1], "r") as f:
        read_length= 151
        row_l = []
        col_l= []
        val_l= []
        line_number = 0
        for line in f:
            sl = line.split()
            #exclude = False
            id = sl[0]
            index = int(sl[2]) - 1
            sam_q = sl[3]
            cigar = sl[4]
            cigar_str = re.findall(r"[0-9]+[MIDSH]", cigar)
            blk_pos = []
            blk_type = []
            ini = 0
            tmp_length = 0

            matched = 0
            for block in cigar_str:
                m = re.search(r'([0-9]+)([MIDSH])', block)
                bl = int(m.group(1)) + ini
                bt = str(m.group(2))

                if bt == "M":  # or bt=="I" or bt=="D":
                    matched += int(m.group(1))
                blk_type.append(bt)

                blk_pos.append(bl)
                ini = bl

                tmp_length += int(m.group(1))  # get read length without counting clipped bases

            if (matched / tmp_length < 0.6):
                unlikely.append([id,sam_q])
            #    exclude_reads.add(i)
                # print(i,str(r.loc[i]), matched)
           #     exclude = True
           #     continue
            if narrow:
                continue
            #narrowed_read.append([str(r["ID"].loc[i]), int(r["strand"].loc[i]), index + 1, str(r["sam_q"].loc[i]), cigar])
            # print("sam_q:", sam_q, "\n")
            sam_q_num = []

            # softclipping adjust index position
            s_start = 0
            s_end = 0

            c = 0
            inserts = []
            begin = 0
            reduce = 0  # deduct "non-existent" bases in total length, "D" and "H" not shown in reads
            for j in range(0, blk_pos[-1]):  # change here to fill the blank with 0?

                if blk_type[c] == "M":
                    try:
                        sam_q_num.append(dict_tonu[sam_q[j - reduce]])
                    except:
                        print(j - reduce, len(sam_q), j)
                        exit(1)


                elif blk_type[c] == "I":
                    inserts.append(dict_tonu[sam_q[j - reduce]])
                elif blk_type[c] == "S":
                    sam_q_num.append(dict_tonu[sam_q[j - reduce]])
                elif blk_type[c] == "H":
                    sam_q_num.append(0)
                elif blk_type[c] == "D":
                    sam_q_num.append(6)
                if blk_type[c] == "H" or blk_type[c] == "D":
                    reduce += 1

                if j == blk_pos[c] - 1:  # update start and c, put inserts into hashtable
                    '''
                    if blk_type[c] == "I":
                        if i in insertion_reads.keys():

                            newinsert = insertion_reads.get(included_i).copy()
                            newinsert.append((index + begin, inserts))
                            insertion_reads.update({included_i: newinsert.copy()})
                        else:
                            insertion_reads.update({included_i: [(index + begin, inserts.copy())]})
                    '''
                    begin = blk_pos[c]
                    inserts = []

                    c += 1
                    if c == len(blk_type):
                        break

            if blk_type[0] == "S":
                if index - blk_pos[0] < 0:
                    start_pos = blk_pos[0] - index
                    sam_q_num = sam_q_num[start_pos:]
                else:
                    index = index - blk_pos[0]

            if len(sam_q_num) < read_length:
                sam_q_num += [0] * (read_length - len(sam_q_num))
                pad = 0
            else:
                pad = len(sam_q_num) - read_length

            val_l.extend(sam_q_num)
            row_tmp = [line_number for n in range(read_length + pad)]
            row_l.extend(row_tmp)
            col_tmp = [n for n in range(index, index + read_length + pad)]
            col_l.extend(col_tmp)
            if len(sam_q_num) != len(row_tmp) or len(sam_q_num) != len(col_tmp):
                #print(r["ID"].loc[i], i, index, len(row_tmp), len(col_tmp), len(sam_q_num))
                print(sam_q_num)
                exit(14)

            line_number += 1
        if narrow:
            with open("unlikely.fa", "w+") as f:
                for r in unlikely:
                    f.write(">"+r[0]+"\n")
                    f.write(r[1]+"\n")
        else:
            row = np.array(row_l)
            col = np.array(col_l)
            val = np.array(val_l)

            csc = sp.coo_matrix((val, (row, col))).tocsc()
            depths = []
            for i in range(0,csc.shape[1]):
                tmp = np.squeeze(csc.getcol(i).toarray())
                tmp_count = sum(np.bincount(tmp)[1:])

                depths.append(tmp_count)
            print(min(depths), max(depths))
            print("average depth for ", sys.argv[1]," is ", statistics.mean(depths))
            print("median depth for ", sys.argv[1], "is ", statistics.median(depths))

def extract_to_read():
    with open(sys.argv[1],"r")as f:
        with open(sys.argv[2],"w+") as w:
            for line in f:
                sl = line.split(" ")
                w.write(">"+sl[0]+"\n")
                w.write(sl[3]+"\n")

def missing_narrow():
    match151 = {}
    with open("299_extract.sam","r") as ef:
        for line in ef:
            sl = line.split(" ")
            sl[4] = re.sub('\n','',sl[4])
            if sl[4] == "151M":
                if sl[0] not in match151.keys():
                    match151.update({sl[0]:1})
                else:
                    match151[sl[0]] += 1
    x = 0
    for k,v in match151.items():
        x += v

    print(x ," 151M read IDs")
    missed = set({})
    narrowed = {}
    with open("narrowed_extract.sam","r") as rf:
        for line in rf:
            sl = line.split(" ")
            sl[4] = re.sub('\n', '', sl[4])
            if sl[4] == "151M":
                if sl[0] not in narrowed.keys():
                    narrowed.update({sl[0]:1})
                else:
                    narrowed[sl[0]] += 1

    for k,v in match151.items():
        if k not in narrowed or v != narrowed[k]:
            missed.add(k)

    print(len(missed)," 151M reads missed in narrowing reads")
    print(missed)


def check_order():
    error_lines = {}
    count = 1
    with open(sys.argv[1],"r") as f1, open(sys.argv[2],"r") as f2:
        for line1, line2 in zip(f1,f2):
            if count %2 != 1:
                count += 1
                continue
            sl1 = line1.split(" ")[0]
            sl2 = line2.split(" ")[0]
            if sl1 != sl2:
                print("order error", sl1, sl2,count)
                error_lines.update({count-1:[sl1,sl2]})
            count += 1
    if len(error_lines) == 0:
        print("order is fine")
    else:
        with open("no match ID.txt", "w+") as errorf:
            for k,v in error_lines:
                errorf.write(str(k)+str(v))

def get_read_from_listf(id_file,readfile1,readfile2):
    if len(sys.argv) < 4:
        print("short of args")
        exit(-1)
    IDs = set({})
    with open(id_file,"r") as f:
        for line in f:
            #Eline = line.replace("[","")
            #line = line.replace("]","")
            #line = line.replace("'","")
            #line = line.replace(",","")
            IDs = set(line.split(" "))
    print("getting", len(IDs), " IDs of reads from " + readfile1 + " and " + readfile2)
    #exit(1)
    with mp.Pool(2) as pool:
        lparam = [(IDs,readfile1, "s_newest_reduced_R1.fastq"),(IDs, readfile2, "s_newest_reduced_R2.fastq")]
        pool.starmap(search_ID,lparam)
def get_ori_half(samfile,readfile1,readfile2):
    IDs = set({})
 #   if not re.match('.*\.sam',sys.argv[1]) or not re.match('.*\.fastq',sys.argv[2]) or not re.match('.*\.fastq',sys.argv[3]):
 #      print("accept sam as first arg")
 #      exit(1)
    with open(samfile,"r") as f:
        for line in f:
            #print(len(line.split(" ")),"IDs counted")
            IDs.add(line.split(" ")[0])
    #print(IDs)
    print("getting", len(IDs) ," IDs of reads from "+readfile1+" and "+readfile2)
    if len(sys.argv) < 4:
        print("get_ori_read() insufficient args")
        exit(2)
    with mp.Pool(2) as pool:
        lparam = [(IDs,readfile1, "half_real_R1.fastq"),(IDs, readfile2, "half_real_R2.fastq")]
        pool.starmap(search_ID,lparam)



def search_ID(id_set,readfile,outfile):
    ori_reads = []
    count = 0
    con_count = 0
    with open(readfile,"r") as r1:
        '''
        for block in iter(lambda: list(its.islice(r1,4)),[]):
            tmpid = re.sub('@','',block[0].split(" ")[0])
            #print(tmpid)
            if tmpid not in id_set:
                continue
            ori_reads.append(block)
            '''

        while True:
            block = list(its.islice(r1, 4))
            count += 1
            if not block:
                break
            tmpid = re.sub('@','',block[0].split(" ")[0])
            #print(tmpid)
            if tmpid not in id_set:
                con_count += 1
                continue
            ori_reads.append(block)

    print(count*4,"lines accessed",con_count*4,"reads skipped",len(ori_reads),"reads to write")

    with open(outfile, "w+") as wr1:
        for lines in ori_reads:
            for line in lines:
                wr1.write(line)


def remove_ID_fastq():
    ID_list = []
    with open(sys.argv[1], "r") as trf:
        for line in trf:
            ID_list.extend(line.split(" "))
    if len(ID_list) == 0:
        print("no reads to be removed")
        exit(2)
    id_set = set(ID_list)
    id_set.remove("")
    print(id_set)
    #exit(1)
    with mp.Pool(2) as pool:
        lparam = [(id_set,sys.argv[2], sys.argv[4]+"remained_reduced_R1.fastq"),(id_set, sys.argv[3],sys.argv[4]+ "remained_reduced_R2.fastq")]
        pool.starmap(remove_by_ID,lparam)

def remove_by_ID(id_set,readfile,outfile):
    remained_reads = []
    with open(readfile,"r") as r1f:
        while True:
            block = list(its.islice(r1f, 4))

            if not block:
                break
            tmpid = re.sub('@','',block[0].split(" ")[0])
            #print(tmpid)
            if tmpid in id_set:
                continue
            remained_reads.append(block)
    print(len(remained_reads),"reads remained")
    with open(outfile,"w+") as r1of:
        for lines in remained_reads:
            for line in lines:
                r1of.write(line)

def get_contigs():
    with open(sys.argv[1],"r") as f:
        for lines in iter(lambda: list(its.islice(f,2)),[]):
            #print(lines)
            filename = "contig_" + lines[0].split(" ")[3].strip()

            with open(filename,"w+") as wf:
                for line in lines:
                    wf.write(line)


def get_rev_comp(input_file,revcomp=False):
    print("getting reverse complement of ",input_file)
    with open(input_file, "r") as f1:

        seq1 = f1.read()
        if seq1[0] == ">":
            f1.seek(0)
            lines = f1.readlines()[1:]
            seq1 = ""
            for l in lines:
                seq1 += l.replace("\n","").replace(" ","")
        else:
            seq1 = seq1.replace("\n","").replace(" ","")
    '''
    seq_rev = seq1[::-1]
    #print(seq1)
    #print(seq_rev)
    with open("rev_"+sys.argv[1],"w+") as f2:
        f2.write(">reverse_of_"+sys.argv[1]+"\n")
        f2.write(seq_rev)
    '''

    seq_revc = seq1.upper()[::-1]
    seq_revc = seq_revc.replace('A', 't')
    seq_revc = seq_revc.replace('T', 'a')
    seq_revc = seq_revc.replace('C', 'g')
    seq_revc = seq_revc.replace('G', 'c')
    seq_rev_comp = seq_revc.upper()
    with open("rev_comp_"+os.path.basename(input_file),"w+") as f2:
        f2.write(">reverse_complement_of_"+input_file+"\n")
        f2.write(seq_rev_comp)
def rev_comp_read(seq):
    seq_revc = seq.upper()[::-1]
    seq_revc = seq_revc.replace('A', 't')
    seq_revc = seq_revc.replace('T', 'a')
    seq_revc = seq_revc.replace('C', 'g')
    seq_revc = seq_revc.replace('G', 'c')
    seq_rev_comp = seq_revc.upper()
    return seq_rev_comp
def fastq_rev_comp():
    count = 1
    with open(sys.argv[1], "r") as f1,open("rev_comp_"+os.path.basename(sys.argv[1]),"w+") as f2:
        for line in f1:
            if count % 2 == 0:
                seq = line.strip()

                f2.write(rev_comp_read(seq)+"\n")
            else:
                f2.write(line)
            count += 1
def dup_read_by_ID():
    ID_list = []
    with open(sys.argv[1], "r") as trf:
        for line in trf:
            ID_list.extend(line.split(" "))
    if len(ID_list) == 0:
        print("no reads to be removed")
        exit(2)
    id_set = set(ID_list)
    id_set.remove("")
    print(id_set)
    read_set = set({})
    rc_read_set = set({})
    dup_read_IDs = {}
    rc_dup_read_IDs = {}
    with open(sys.argv[2],"r") as rnf:
        for line in rnf:
            fields = line.split(" ")
            if fields[0] in id_set:
                flag = format(fields[1].strip(),'b')[::-1]
                if flag[4] == "1":
                    rc_read_set.add(rev_comp_read(fields[3].strip()))
                    rc_dup_read_IDs.update({rev_comp_read(fields[3].strip()): {line.split(" ")[0].strip()}})
                else:
                    read_set.add(fields[3].strip())
                    dup_read_IDs.update({fields[3].strip(): {line.split(" ")[0].strip()}})
    with mp.Pool(2) as pool:
        lparam = [(read_set,dup_read_IDs,sys.argv[2], "dup_reads_R1.txt"),(read_set,dup_read_IDs, sys.argv[3], "dup_reads_R2.txt")]
        pool.starmap(search_by_read,lparam)

def search_by_read(read_set,dup_read_IDs, readfile,outfile):

    with open(readfile,"r") as r1f:
        while True:
            block = list(its.islice(r1f, 4))

            if not block:
                break
            tmpid = re.sub('@','',block[0].split(" ")[0])
            tmpread = block[1].strip()
            #print(tmpid)
            if tmpread not in read_set:
                continue
            dup_read_IDs[tmpread].add(tmpid)
    #print(len(remained_reads),"reads remained")
    #return dup_read_IDs
        with open(outfile,"w+") as r1of:
            for lines in dup_read_IDs:
                for line in lines:
                    r1of.write(line)
def split_seq(input_file,length=9050,step=1000):
    seq = ""
    with open(input_file,"r") as f:
        for line in f:
            if line[0] != ">":
                seq += line.strip()
    for i in range(0,len(seq)-length,step):
        with open("slice_"+str(i)+"_"+str(i+length)+"_"+ os.path.basename(input_file),"w+") as wf:
            wf.write(">index"+str(i)+"_"+str(i+length)+"\n")
            wf.write(seq[i:i+length])

def count_pair_dist(input_file):
    dist = {}
    with open(input_file,"r") as f:
        for line in f:
            sl = line.split(" ")
            if sl[0] in dist.keys():
                dist[sl[0]].append(int(sl[2]))
            else:
                dist.update({sl[0]:[int(sl[2])]})
    dist_values = []
    for i in dist.keys():
        dist_values.append(abs(dist[i][0]-dist[i][1]))
    dist_val_count = Counter(dist_values)
    print(len(dist_values),len(dist),dist_val_count, max([x for x in dist_val_count.values()]))
    print(statistics.mean(dist_values),statistics.median(dist_values))

def test_shell(readfile1,readfile2):
    verify_sub_command = os.path.dirname(
        __file__) + "/find_sub.sh" + " " + "-r" + " "  + "oneline_NC045512.2.fasta_output/mt_read_subbed_0.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True'+ " -d Y"
    # print(verify_sub_command)
    # verify_sub_command= shlex.split(verify_sub_command)

    verify_proc = subprocess.Popen(verify_sub_command,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
    # verify_proc = subprocess.run(os.path.dirname(__file__)+" count_error.sh 1 2 3",shell=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    for line in verify_proc.stdout:
        print(line,end="")
    verify_proc.stdout.close()
    proc_code = verify_proc.wait()
    if proc_code:
        raise subprocess.CalledProcessError(proc_code,verify_sub_command)
    print(verify_proc.args)
    print(verify_proc.stdout)

def get_gap_reads(ref,strain,samfile ,readfile1,readfile2):

    #strain = 2
    if strain < 0:
        print("invalid strain number, error")
        exit(-2)
    prev_readset = set({})
    if strain > 0:
        for i in range(0,strain):
            with open("subbed_reads_"+str(0)+".sam","r") as subf:
                for line in subf:
                    sline = line.strip().split(" ")
                    prev_readset.add(sline[3])
    readlist=[]
    with open(samfile,"r") as samf:
        for line in samf:
            if line.split(" ")[3] not in prev_readset:
                readlist.append(line.strip().split(" "))
    ref_seq = ""
    with open(ref,"r") as rf:
        for line in rf:
            if line[0]!=">":
                ref_seq += line.strip()
    ref = ref_seq

    for read in readlist:
        tmp_misP = []
        read_index = int(read[2])-1
        for i, base in enumerate(read[3]):
            if ref[read_index + i] != base :
                tmp_misP.append(read_index + i)
        for mp in tmp_misP:
            read.append((int(mp),ref[mp],read[3][mp-read_index]))

    batch = 0


    while len(readlist) > 0:
        subbed_read = []
        misPs = []
        ref = ref_seq

        for batch,read in enumerate(readlist):
            overlap = False
            read_index = int(read[2])-1
            read_misPs = [int(x[0]) for x in read[5:]]
            if len(subbed_read) == 0:
                subbed_read.append(read)
                misPs.extend(read_misPs)
            else:
                if read_index > max(misPs) or read_index+len(read[3]) < min(misPs):
                    subbed_read.append(read)
                    misPs.extend(read_misPs)

                else:
                    for mp in misPs:
                        if read_index <= mp and read_index+len(read[3]) >= mp:
                            misp_found = False
                            misp_overlap = False
                            for r_misp in read[5:]:
                                if mp == r_misp[0]:
                                    for tmp_r in subbed_read:
                                        for subbed_r_misp in tmp_r[5:]:
                                            if r_misp[0] == subbed_r_misp[0]:
                                                misp_found = True
                                                #conflict by another misp
                                                if r_misp[2] != subbed_r_misp[2]:
                                                    misp_overlap = True
                                                    break
                                if misp_overlap:
                                    overlap = True
                                    break
                            # conflict by ref
                            if not misp_found:
                                overlap = True
                                break

                    if overlap:

                        batch += 1
                        continue
                    else:
                        subbed_read.append(read)
                        misPs.extend(read_misPs)

            ref = ref[:read_index] + read[3] + ref[read_index + len(read[3]):]
            with open("batch_" + str(batch) + "_reference.fa", "w+") as bf:
                bf.write(">batch_" + str(batch) + "\n")
                bf.write(ref)
            verify_sub_command = os.path.dirname(
                __file__) + "/find_sub.sh" + " " + "-r" + " "  +"batch_" + str(
                batch) + "_reference.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"

            verify_proc = subprocess.run(verify_sub_command,
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                         shell=True)
            print(verify_proc.stdout.split("\n")[-5:])
            gap_list = verify_proc.stdout.split("\n")[-3]
            if gap_list != "no gaps":

                print(gap_list)
                gap_list = re.sub('[\[\]]', '', gap_list)
                print(gap_list.split(",")[:-1])
                print(read)
                #gap_list = gap_list.split(",")[:-1]
                ref = ref[:read_index] + ref_seq[read_index:read_index+len(read[3])] + ref[read_index + len(read[3]):]
                subbed_read.remove(read)
            else:
                print(batch,"no gaps for ",subbed_read)

            batch += 1
        with open("final_strain_" + str(strain) + "_reference.fa", "w+") as bf:
                bf.write(">final_strain_" + str(strain) + "\n")
                bf.write(ref)
        with open("subbed_reads_"+str(strain)+".sam", "w+") as bf:
            for line in subbed_read:
                bf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
        #subprocess.run("rm -r batch_*",shell=True)
        print(subbed_read)
        for read in subbed_read:
            readlist.remove(read)

        strain += 1
        exit()

        '''
        for read in subbed_read:
            print(read)
            read_index = int(read[2]) -1
            ref = ref[:read_index] + read[3] + ref[read_index+len(read[3]):]
            readlist.remove(read)
        with open("batch_" +str(batch)+"_reference.fa","w+") as bf:
            bf.write(">batch_"+str(batch)+"\n")
            bf.write(ref)
        print(len(readlist),"reads left after batch",batch,"diff")
        verify_sub_command = os.path.dirname(
            __file__) + "/find_sub.sh" + " " + "-r" + " " + "batch_" +str(batch)+"_reference.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"
        #verify_sub_command = " python3 $insec/strain_finder.py --gene_L=29903  --ref=/mnt/c/Users/tan/Desktop/UTS/mafft-win/InsEC/spike_test/batch_"+str(batch)+"_reference.fa "+ "--narrowing=True --match_l=1 --sam_file=batch_"+str(batch)+"_reference.fa_output"  + "/round_1_extract.sam --r1_file=/mnt/c/Users/tan/Desktop/UTS/mafft-win/InsEC/corrected_shi_reduced_R1.fastq --r2_file=/mnt/c/Users/tan/Desktop/UTS/mafft-win/InsEC/corrected_shi_reduced_R2.fastq --round=1 --excluded_IDs=excluded_IDs.txt --bubble_mode=True --find_sub=True --check_gap=True"
        # print(verify_sub_command)
        # verify_sub_command= shlex.split(verify_sub_command)

        verify_proc = subprocess.run(verify_sub_command,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                       shell=True)
        # verify_proc = subprocess.run(os.path.dirname(__file__)+" count_error.sh 1 2 3",shell=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)

        #print(verify_proc.args)

        gap_list = verify_proc.stdout.split("\n")[-3]
        if gap_list != "no gaps":

            print(gap_list)
            gap_list = re.sub('[\[\]]','',gap_list)
            print(gap_list.split(",")[:-1])
            gap_list = gap_list.split(",")[:-1]
            gapped_reads = []
            for read in subbed_read:
                read_index = int(read[2]) - 1
                read_misPs = [int(x[0]) for x in read[5:]]
                for gap in gap_list:
                    gap = int(gap)
                    if gap > read_misPs[0]-150 and gap < read_misPs[-1]+150:
                        gapped_reads.append(read)
                        break

            if len(gapped_reads) > 0:
                print(gapped_reads)
                with open("gapped_reads"+str(batch)+".sam","w+") as grf:
                    for line in gapped_reads:
                        grf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
        
        else:
            batch += 1
            continue
        
        #exit()
        batch += 1
        '''
def get_output():
    out = ""
    with open("test_output.txt","r") as tf:
        for line in tf:
            out += line
    print(out.split("\n")[-2])
def count_fq_freq(fastqfile):
    readlist=[]
    with open(fastqfile,"r") as fq:
        count = 0
        while True:
            block = list(its.islice(fq, 4))

            if not block:
                break
            count += 1
            readlist.append(block[1].strip())
    print(len(readlist), count)
    read_freq = Counter(readlist)
    count = {}
    print(len(read_freq))

    for k,v in read_freq.items():
        if v not in count.keys():
            count.update({v:1})
        else:
            count[v] += 1

        #print(k,v)
    sum_freq = 0
    for k,v in {k:v for k,v in sorted(count.items(),key=lambda  x : x[1] )}.items():
        print(k,v)
        sum_freq += k*v
    print(sum_freq)

def remove_freq(samfile,freq):
    readfreq = {}
    readlist = []
    target_id = set({})
    with open(samfile,"r") as f:
        for line in f:
            fields = line.strip().split(" ")
            readlist.append(fields)
            #print(fields)
            if len(fields)>0 and fields[3] not in readfreq.keys():
                readfreq.update({fields[3]:1})
                #target_id.add(fields[0])
            else:
                readfreq[fields[3]] += 1
                #target_id.remove(fields[0])
    print(len(readlist),"reads in ",samfile, " uniq reads", len(readfreq))
    filtered_readlist = []
    for read in readlist:
        if readfreq[read[3]] <= freq:
            target_id.add(read[0])
        else:
            filtered_readlist.append(read)
    print(len(filtered_readlist),"reads with frequency > ",freq)
    lastpos = -1
    gapped_index = set({})
    for read in filtered_readlist:
        if lastpos > 0 and  lastpos < int(read[2]) :
            if int(read[2]) not in gapped_index:
                #gapped_index.add(int(read[2]))
                print("gap between ",prev, read)
                #print(read)
        if read[0] not in target_id:
            if lastpos == -1:
                print("first read", read)
            prev = read
            lastpos = int(read[2]) + len(read[3])
    #print(len(gapped_index), " reads have gaps between it and the read before ")
    #print(gapped_index)
    if lastpos < 29800:
        print("last read end at",lastpos)
    #print(count ,"reads have frequency 1")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="some tool functions",usage='%(prog)s [options]'+str(sys.argv))
    parser.add_argument('--rev_comp',type=str,default="no")
    parser.add_argument("--remove_by_ID",type=str,default="no")
    parser.add_argument("files",type=str,default="",nargs="+")
    parser.add_argument("--split_seq",type=str,default="no")
    parser.add_argument("--pair_dist", type=str,default="no")
    parser.add_argument("--overlap_id",type=str,default="no")
    parser.add_argument("--get_read_by_id",type=str,default="no")
    parser.add_argument("--test_shell",type=str,default="no")
    parser.add_argument("--get_gap_reads",type=str,default="no")
    parser.add_argument("--get_ori_half",type=str,default="no")
    parser.add_argument("--count_fq_freq",type=str,default="no")
    parser.add_argument("--strain_num",type=int,default=-1)
    parser.add_argument("--remove_freq", type=int,default=-1)
    args = parser.parse_args()

    input_files = args.files
    print(input_files[1:])
    # hash_str()
    #count_match()
    #countdiff()
    #count_record()
    #oneline_fasta()
    if args.overlap_id != "no":

        getoverlap(input_files[1],input_files[2])
    #get_sam_dup()

    #adjustindex(377)
    #sep_reads()

    #seq_depth()

    #extract_to_read()
    #seq_depth()
    #hash_str()
    #count_length()
    #missing_narrow()

    #check_order()
    #fix_ID()

    #get_contigs()
    if args.rev_comp != "no":
        get_rev_comp(input_files[1],True)
    if args.get_read_by_id != "no":
        get_read_from_listf(input_files[1],input_files[2],input_files[3])
    #rev_comp_read()
    if args.remove_by_ID != "no":
        remove_ID_fastq()
    if args.split_seq != "no":
        split_seq(input_files[1],length=2949,step=300)
    if args.pair_dist != "no":
        count_pair_dist(input_files[1])
    if args.test_shell != "no":
        test_shell(input_files[1],input_files[2])
    if args.get_gap_reads != "no":
        get_gap_reads(input_files[1],args.strain_num,input_files[2],input_files[3],input_files[4])
    if args.get_ori_half != "no":
        get_ori_half(input_files[1],input_files[2],input_files[3])
    if args.count_fq_freq != "no":
        count_fq_freq(input_files[1])
    if args.remove_freq > 0:
        remove_freq(input_files[1], args.remove_freq)
    corenum = mp.cpu_count()-2
    line_amount = 245216120

    manager = mp.Manager()
    lock = manager.Lock()
    dup_count = {}
    #print(dup_count)
    '''
    with mp.Pool(corenum) as pool:
        start = 0
        lparam = []
        while start < line_amount:
            lparam.append((start,start+int(line_amount/corenum),lock,dup_count))
            start += int(line_amount/corenum)+1

        results = pool.starmap(hash_str,lparam )
    '''
    #hash_str(0,line_amount,lock,dup_count)
    #print(dup_count)

    #with open("dup_reads_1","w+") as f:
    #    for i in dup_count.keys():
    #        #if(dup_count.get(i)[1]>1):
    #        f.write(i+str(dup_count.get(i))+" \n")
    #count_match()
