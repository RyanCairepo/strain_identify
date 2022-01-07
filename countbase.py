import os
import random
import sys, re, statistics, itertools as its, multiprocessing as mp, numpy as np, scipy.sparse as sp, matplotlib.pyplot as pit
from typing import Counter

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

def getoverlap():
    read_set1 = {}
    overlapped = {}
    read_set2 = {}
    with open(sys.argv[1], "r") as f1:
        for line in f1:
            sl = line.split(" ")
            if sl[0] not in read_set1.keys():
                read_set1.update({sl[0]: {sl[3]}})
            else:
                read_set1[sl[0]].add(sl[3])
    with open(sys.argv[2], "r") as f2:
        for line in f2:
            sl = line.split(" ")
            if sl[0] not in read_set2.keys():
                read_set2.update({sl[0]: {sl[3]}})
            else:
                read_set2[sl[0]].add(sl[3])
    for rid in read_set1.keys():

        if rid in read_set2.keys():
            if len(read_set1[rid].intersection(read_set2[rid]))>1:
                overlapped.update({rid:read_set1[rid].intersection(read_set2[rid])})

    print(len(overlapped.keys()), "reads ID overlapped ", len(read_set1.keys()), " unique read ID in ", sys.argv[1], " and ",len(read_set2.keys()),
          " in ", sys.argv[2])
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

def get_read_from_listf():
    if len(sys.argv) < 4:
        print("short of args")
        exit(-1)
    IDs = set({})
    with open(sys.argv[1],"r") as f:
        for line in f:
            #Eline = line.replace("[","")
            #line = line.replace("]","")
            #line = line.replace("'","")
            #line = line.replace(",","")
            IDs = set(line.split(" "))
    print("getting", len(IDs), " IDs of reads from " + sys.argv[2] + " and " + sys.argv[3])
    #exit(1)
    with mp.Pool(2) as pool:
        lparam = [(IDs,sys.argv[2], "shi_reduced_R1.fastq"),(IDs, sys.argv[3], "shi_reduced_R2.fastq")]
        pool.starmap(search_ID,lparam)
def get_ori_half():
    IDs = set({})
    if not re.match('.*\.sam',sys.argv[1]) or not re.match('.*\.fastq',sys.argv[2]) or not re.match('.*\.fastq',sys.argv[3]):
       print("accept sam as first arg")
       exit(1)
    with open(sys.argv[1],"r") as f:
        for line in f:
            IDs.add(line.split(" ")[0])
    #print(IDs)
    print("getting", len(IDs) ," IDs of reads from "+sys.argv[2]+" and "+sys.argv[3])
    if len(sys.argv) < 4:
        print("get_ori_read() insufficient args")
        exit(2)
    with mp.Pool(2) as pool:
        lparam = [(IDs,sys.argv[2], "half_real_R1.fastq"),(IDs, sys.argv[3], "half_real_R2.fastq")]
        pool.starmap(search_ID,lparam)



def search_ID(id_set,readfile,outfile):
    ori_reads = []
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

            if not block:
                break
            tmpid = re.sub('@','',block[0].split(" ")[0])
            #print(tmpid)
            if tmpid not in id_set:
                continue
            ori_reads.append(block)



    with open(outfile, "w+") as wr1:
        for lines in ori_reads:
            for line in lines:
                wr1.write(line)

def get_contigs():
    with open(sys.argv[1],"r") as f:
        for lines in iter(lambda: list(its.islice(f,2)),[]):
            #print(lines)
            filename = "contig_" + lines[0].split(" ")[3].strip()

            with open(filename,"w+") as wf:
                for line in lines:
                    wf.write(line)


def get_rev_comp(revcomp=False):
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
    with open("rev_comp_"+sys.argv[1],"w+") as f2:
        f2.write(">reverse_complement_of_"+sys.argv[1]+"\n")
        f2.write(seq_rev_comp)
def rev_comp_read():
    count = 1
    with open(sys.argv[1], "r") as f1,open("rev_comp_"+os.path.basename(sys.argv[1]),"w+") as f2:
        for line in f1:
            if count % 2 == 0:
                seq = line.strip()
                seq_revc = seq.upper()[::-1]
                seq_revc = seq_revc.replace('A', 't')
                seq_revc = seq_revc.replace('T', 'a')
                seq_revc = seq_revc.replace('C', 'g')
                seq_revc = seq_revc.replace('G', 'c')
                seq_rev_comp = seq_revc.upper()
                f2.write(seq_rev_comp+"\n")
            else:
                f2.write(line)
            count += 1
if __name__ == "__main__":
    # hash_str()
    #count_match()
    #countdiff()
    #count_record()
    #oneline_fasta()
    #getoverlap()
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
    #get_ori_half()
    #get_contigs()
    get_rev_comp(True)
    #get_read_from_listf()
    #rev_comp_read()


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
