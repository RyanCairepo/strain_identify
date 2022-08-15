import copy
import csv
import os
import random
import subprocess
import sys, re, statistics, itertools as its, multiprocessing as mp, numpy as np, scipy.sparse as sp, \
	matplotlib.pyplot as pit
from typing import Counter
import argparse, math
import editdistance
import matplotlib.pyplot as plt
import scipy.interpolate
import build_matrix as bm
import strain_finder as st_find


def countbase():
	baselength = {}
	with open(sys.argv[1]) as f:
		for line in f:
			sl = line.split(" ")
			if len(sl) > 1:

				length = int(re.findall(r'[0-9]+', sl[2])[0])

				if (length in baselength.keys()):
					baselength.update({length: baselength.get(length) + 1})
				else:
					baselength.update({length: 1})
	sortedlength = {k: v for k, v in sorted(baselength.items(), key=lambda item: item[1])}
	print(sortedlength)


def countmatch():
	count = 0
	lengths = {}
	with open("gene2.match", "w+") as fw:
		with open("SRR15096407.sam", "r+") as f:
			for line in f:
				sl = re.split(' |\t', line)
				# print(len(sl))
				if (len(sl) > 3 and sl[5] != "*"):
					fw.write(str(sl[5]) + " " + str(len(sl[9])) + "\n")
					if len(sl[9]) not in lengths.keys():
						lengths.update({len(sl[9]): 1})
					else:
						lengths.update({len(sl[9]): lengths.get(len(sl[9])) + 1})
					count += 1

	print(count)
	sortedlength = {k: v for k, v in sorted(lengths.items(), key=lambda item: item[1])}
	print(sortedlength)


def combine():
	os.chdir("./sra_downloads")
	with open("../combinedsra.fastq", "a+") as wf:
		for file in os.listdir(os.getcwd()):
			with open(file, "r") as rf:
				# print(file)
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
				seq1 += l.replace("\n", "").replace(" ", "")
		else:
			seq1 = seq1.replace("\n", "").replace(" ", "")
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
	# print(len(seq1), len(seq2))
	diffnum = 0
	diff = {}
	if (len(seq1) < len(seq2)):
		for i in range(len(seq2)):
			if i < len(seq1):
				if seq1[i] != seq2[i]:
					# print("#"+str(i)+": "+seq1[i]+"|"+seq2[i])
					diffnum += 1
					diff.update({i: seq1[i] + "|" + seq2[i]})
			else:
				diffnum += 1
				diff.update({i: "*|" + seq2[i]})
			# print("#"+str(i)+" seq2: "+seq2[i])

	else:
		for i in range(len(seq1)):
			if i < len(seq2):
				if seq1[i] != seq2[i]:
					# print("#"+str(i)+": "+seq1[i]+"|"+seq2[i])
					diffnum += 1
					diff.update({i: seq1[i] + "|" + seq2[i]})
			else:
				diffnum += 0
				diff.update({i: seq1[i] + "|*"})
			# print("#"+str(i)+" seq1: "+seq1[i])
	# print(diffnum)
	if diffnum > 0:
		s_diff = sorted(diff.keys())
		c = 0
		for k in s_diff:
			eol = "\n" if c % 15 == 0 and c != 0 else "   "
			# print(str(k) + ": " + diff.get(k), end=eol)
			c += 1
	print(editdistance.eval(seq1, seq2))


# def splitsequence():
def oneline_fasta(fasta_file):
	with open(fasta_file, "r") as f:
		seq = ""
		name = ""
		for line in f:
			if line[0] != ">":

				temp = line.replace(" ", "").replace("\n", "").replace("\d", "").upper()
				seq += re.sub('\d', "", temp)
			else:
				name = line
		with open("oneline_" + fasta_file.replace("..", "").replace("/", ""), "w+") as fw:
			fw.write(name)
			fw.write(seq)


def count_record():
	record1 = {}
	with open(sys.argv[1], "r") as f1:
		for line in f1:
			matches = re.findall('[0-9]+: [\-]*[A-Z*]+\|[A-Z*]+', line)
			for m in matches:
				parts = re.search('([0-9]+): ([\- ]*[A-Z*]+\|[A-Z*]+)', m)
				index = int(parts.group(1))
				cor = parts.group(2)
				record1.update({index: cor})
	record2 = {}
	with open(sys.argv[2], "r") as f2:
		for line in f2:
			matches = re.findall('[0-9]+: [\-]*[A-Z*]+\|[A-Z*]+', line)
			for m in matches:
				parts = re.search('([0-9]+): ([\- ]*[A-Z*]+\|[A-Z*]+)', m)
				index = int(parts.group(1))
				cor = parts.group(2)
				record2.update({index: cor})
	count = 0
	for k in record1.keys():
		if k in record2.keys():
			if record1[k] == record2[k]:
				count += 1
	print(count, " corrections are the same")


def count_match():  # matched bases
	readlist = []
	matches = []
	tmp = 0
	narrowed = []
	count = 0

	with open(sys.argv[1], "r") as f:

		for line in f:

			s = line.split(" ")
			readlist.append(s)
			m = 0

			if len(s) < 3:
				continue
			# print("her?")
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
		# print(line)

	ID_count = {}
	for read in narrowed:
		if read[0] not in ID_count.keys():
			ID_count.update({read[0]: 1})
		else:
			ID_count[read[0]] += 1
	strand_count = {}
	single = []
	for id in ID_count.keys():
		if ID_count[id] < 2:
			single.append(id)
	print(len(readlist), len(narrowed), count, len(single))
	# print(single)
	for read in narrowed:

		if read[0] in single:
			continue
		# if not re.match('^[0-9]+[M]$', read[4]):
		if int(read[1]) not in strand_count.keys():
			strand_count.update({int(read[1]): 1})
		else:
			strand_count[int(read[1])] += 1
	# print(strand_count)
	print({k: v for k, v in sorted(strand_count.items(), key=lambda item: item[0])})
	field_count = {}
	for strand in strand_count.keys():
		bin = "{0:b}".format(int(strand))

		i = 0
		for i in range(0, len(bin)):
			if bin[i] == "1":
				if i not in field_count:
					field_count[i] = strand_count[strand]
				else:
					field_count[i] += strand_count[strand]
	# print(field_count)
	print({k: v for k, v in sorted(field_count.items(), key=lambda item: item[0])})

	# with open("raw_high_match_extract.sam","w+") as f:
	#    for n in narrowed:
	#        f.write(line)

	print(len(matches), "reads met standard")
	print("average matches: ", statistics.mean(matches))
	print("median match: ", statistics.median(matches))
	matches.sort()
	print("30% ", matches[int(len(matches) * 0.3)], "10%", matches[int(len(matches) * 0.1)])


def hash_str():
	i = 0
	# hexp = hash(exp[0])
	# ----------------dup_count is a dict, uses read as the key and a list
	# [frequency, starting position, cigar, match] as value----------------------
	with open(sys.argv[1], "r") as f:
		dup_count = {}
		for line in f:

			sl = line.split(" ")
			sl[4] = re.sub('\n', '', sl[4])
			read_id = sl[0]

			read = sl[3]
			match_number = cigar_processing(sl[4])
			# if match_number < 145:# or match_number/len(sl[3]) < 0.97:
			if len(read) < 145:
				continue
			# hr = hash(read)
			# duplicate read found
			if (read in dup_count.keys()):
				# value_list = dup_count[read]
				# see if starting index is different
				# if value_list[1] != int(sl[2]):

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
				# else:
				if dup_count[read][3] < match_number:
					dup_count[read][1] = int(sl[2])
					dup_count[read][2] = sl[4]
					dup_count[read][3] = match_number

				dup_count[read][0] += 1

			# dup_count.update({read:[dup_count[read][0]+1, dup_count[read][1],sl[4]]})
			else:
				# [frequency,index]
				dup_count.update({read: [1, int(sl[2]), sl[4], match_number]})
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
	# sort here by frequency
	sorted_distinct = {k: v for k, v in sorted(dup_count.items(), reverse=True, key=lambda x: x[1][0])}
	with open("distinct_report.txt", "w+") as dwf:
		for k, v in sorted_distinct.items():
			dwf.write(k + " " + str(v[0]) + "\n")
		# print(k, "f:", v[0], "p:", v[1], "c:", v[2])#, "n", v[4])

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
	# print(seq1_dict)
	sorted_seq1_dict = {k: v for k, v in sorted(seq1_dict.items(), reverse=True, key=lambda x: x[1])}
	with open("pengyao_ferq.txt", "w+") as pf, open("compare.txt", "w+") as wf:
		if sorted_seq1_dict != sorted_distinct:
			for k in sorted_seq1_dict.keys():
				k1 = re.sub('\n', '', k)
				pf.write(k + " " + str(sorted_seq1_dict[k]) + "\n")
				if len(k1) >= 145 and (k1 not in sorted_distinct or sorted_seq1_dict[k] != sorted_distinct[k1][0]):
					if k1 in sorted_distinct.keys():
						rv = sorted_distinct[k1]
					else:
						rv = "none"

					wf.write("p: " + k1 + " " + str(sorted_seq1_dict[k]) + "\n" + "r:" + str(rv) + "\n")
			for k2 in sorted_distinct.keys():
				if k2 not in sorted_seq1_dict or (sorted_distinct[k2][0] - sorted_seq1_dict[k2 + "\n"]) != 0:
					if k2 + "\n" in sorted_seq1_dict.keys():
						rv = str(sorted_seq1_dict[k2 + "\n"]) + "!=" + str(sorted_distinct[k2][0])
						print((sorted_distinct[k2][0] - sorted_seq1_dict[k2 + "\n"]), type(sorted_distinct[k2][0]),
							  type(sorted_seq1_dict[k2 + "\n"]))

					else:
						rv = "none"

					wf.write("r: " + k2 + " " + str(sorted_distinct[k2][0]) + "\n" + "p:" + str(rv) + "\n")

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

# print(dup_count)
# print(read)


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
	with open(sys.argv[1], "r") as f:
		for line in f:
			sl = line.split(" ")
			read_len = len(sl[3])
			if read_len in len_count.keys():
				len_count[read_len] += 1
			else:
				len_count[read_len] = 1
	s_count = {k: v for k, v in sorted(len_count.items(), key=lambda i: i[1], reverse=True)}
	print(s_count)
	x = np.asarray([k for k in s_count.keys()])
	y = np.asarray([v for k, v in s_count.items()])
	px = 1 / plt.rcParams['figure.dpi']
	fig, ax = plt.subplots(figsize=(1900 * px, 600 * px))
	print(y)
	print(x)
	plt.bar(x, y, align='center')
	plt.xticks(x)
	plt.xlabel('read lengths')
	plt.ylabel('frequency')
	ax.set(xlim=[100, 156])
	# plt.xticks([x for x in range(100,153)])
	# ax.xaxis.set_tick_params(width=5)
	# fig, axs = plt.subplots(1,2, sharey=True, tight_layout=True)

	print(fig)
	# axs[0].hist(x, bins=20)
	# axs[1].hist(y,bins=20)
	plt.show()


def getoverlap(f1, f2):
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
				if len(read_set1[rid].intersection(read_set2[rid])) == len(read_set1[rid]):
					overlapped.update({rid: read_set1[rid].intersection(read_set2[rid])})
				else:
					print(read_set1[rid], read_set2[rid])
			else:

				rs2_only.extend(list(read_set2[rid]))

	print(len(overlapped.keys()), "reads ID overlapped ", len(read_set1.keys()), " unique read ID in ", f1, " and ",
		  len(read_set2.keys()),
		  " in ", f2)
	overlap_ids = set(read_set1.keys()).intersection(set(read_set2.keys()))
	print(len(read_set1), len(read_set2))
	print(len(overlap_ids))
	print(len(rs2_only), "reads only in " + f2.name)
	exit()
# with open("overlap_ids.txt","w+")as wf:
#    for ri in overlap_ids:
#        wf.write(ri+",")
# for k in overlapped.keys():
#    print()


def get_sam_dup():
	read_set = set({})
	count = 0
	line_count = 0
	with open(sys.argv[1], "r") as f:
		for line in f:
			line_count += 1
			id = line.split(" ")[0]
			if id in read_set:
				count += 1
			else:
				read_set.add(id)
	print(len(read_set), " unique reads with ", count, "duplicates in total ", line_count)


def adjustindex(moveby):
	"""needs work, plan use for read index adjust"""
	allreads = []
	with open(sys.argv[1], "r") as f:
		for line in f:
			fields = line.split(" ")
			allreads.append(fields.copy())
	with open("after_" + sys.argv[1], "w+") as f1:
		for i in allreads:
			c = 0
			for j in i:
				if c < 4 and c != 2:
					f1.write(j + " ")
				elif c == 2:
					f1.write(str(int(j) - moveby) + " ")
				elif c == 4:
					f1.write(j)
				c += 1
		# f1.write("\n")


def sep_reads():
	with open(sys.argv[1], "r") as f1:
		for line in f1:
			fiels = line.split(" ")
			if "I" in fiels[4] or "D" in fiels[4]:
				print(line, end="")


def seq_depth(narrow=False):
	dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
	dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))
	unlikely = []
	with open(sys.argv[1], "r") as f:
		read_length = 151
		row_l = []
		col_l = []
		val_l = []
		line_number = 0
		for line in f:
			sl = line.split()
			# exclude = False
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
				unlikely.append([id, sam_q])
			#    exclude_reads.add(i)
			# print(i,str(r.loc[i]), matched)
			#     exclude = True
			#     continue
			if narrow:
				continue
			# narrowed_read.append([str(r["ID"].loc[i]), int(r["strand"].loc[i]), index + 1, str(r["sam_q"].loc[i]), cigar])
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
				# print(r["ID"].loc[i], i, index, len(row_tmp), len(col_tmp), len(sam_q_num))
				print(sam_q_num)
				exit(14)

			line_number += 1
		if narrow:
			with open("unlikely.fa", "w+") as f:
				for r in unlikely:
					f.write(">" + r[0] + "\n")
					f.write(r[1] + "\n")
		else:
			row = np.array(row_l)
			col = np.array(col_l)
			val = np.array(val_l)

			csc = sp.coo_matrix((val, (row, col))).tocsc()
			depths = []
			for i in range(0, csc.shape[1]):
				tmp = np.squeeze(csc.getcol(i).toarray())
				tmp_count = sum(np.bincount(tmp)[1:])

				depths.append(tmp_count)
			print(min(depths), max(depths))
			print("average depth for ", sys.argv[1], " is ", statistics.mean(depths))
			print("median depth for ", sys.argv[1], "is ", statistics.median(depths))


def extract_to_read():
	with open(sys.argv[1], "r") as f:
		with open(sys.argv[2], "w+") as w:
			for line in f:
				sl = line.split(" ")
				w.write(">" + sl[0] + "\n")
				w.write(sl[3] + "\n")


def missing_narrow():
	match151 = {}
	with open("299_extract.sam", "r") as ef:
		for line in ef:
			sl = line.split(" ")
			sl[4] = re.sub('\n', '', sl[4])
			if sl[4] == "151M":
				if sl[0] not in match151.keys():
					match151.update({sl[0]: 1})
				else:
					match151[sl[0]] += 1
	x = 0
	for k, v in match151.items():
		x += v

	print(x, " 151M read IDs")
	missed = set({})
	narrowed = {}
	with open("narrowed_extract.sam", "r") as rf:
		for line in rf:
			sl = line.split(" ")
			sl[4] = re.sub('\n', '', sl[4])
			if sl[4] == "151M":
				if sl[0] not in narrowed.keys():
					narrowed.update({sl[0]: 1})
				else:
					narrowed[sl[0]] += 1

	for k, v in match151.items():
		if k not in narrowed or v != narrowed[k]:
			missed.add(k)

	print(len(missed), " 151M reads missed in narrowing reads")
	print(missed)


def check_order():
	error_lines = {}
	count = 1
	with open(sys.argv[1], "r") as f1, open(sys.argv[2], "r") as f2:
		for line1, line2 in zip(f1, f2):
			if count % 2 != 1:
				count += 1
				continue
			sl1 = line1.split(" ")[0]
			sl2 = line2.split(" ")[0]
			if sl1 != sl2:
				print("order error", sl1, sl2, count)
				error_lines.update({count - 1: [sl1, sl2]})
			count += 1
	if len(error_lines) == 0:
		print("order is fine")
	else:
		with open("no match ID.txt", "w+") as errorf:
			for k, v in error_lines:
				errorf.write(str(k) + str(v))


def get_read_from_listf(id_file, readfile1, readfile2):
	if len(sys.argv) < 4:
		print("short of args")
		exit(-1)
	IDs = set({})
	with open(id_file, "r") as f:
		for line in f:
			# Eline = line.replace("[","")
			# line = line.replace("]","")
			# line = line.replace("'","")
			# line = line.replace(",","")
			IDs = set(line.split(" "))
	print("getting", len(IDs), " IDs of reads from " + readfile1 + " and " + readfile2)
	# exit(1)
	with mp.Pool(2) as pool:
		lparam = [(IDs, readfile1, "s_newest_reduced_R1.fastq"), (IDs, readfile2, "s_newest_reduced_R2.fastq")]
		pool.starmap(search_ID, lparam)


def get_ori_half(samfile, readfile1, readfile2):
	IDs = set({})
	#   if not re.match('.*\.sam',sys.argv[1]) or not re.match('.*\.fastq',sys.argv[2]) or not re.match('.*\.fastq',sys.argv[3]):
	#      print("accept sam as first arg")
	#      exit(1)
	with open(samfile, "r") as f:
		for line in f:
			# print(len(line.split(" ")),"IDs counted")
			IDs.add(line.split(" ")[0])
	# print(IDs)
	print("getting", len(IDs), " IDs of reads from " + readfile1 + " and " + readfile2)
	if len(sys.argv) < 4:
		print("get_ori_read() insufficient args")
		exit(2)
	with mp.Pool(2) as pool:
		lparam = [(IDs, readfile1, "half_real_R1.fastq"), (IDs, readfile2, "half_real_R2.fastq")]
		pool.starmap(search_ID, lparam)


def search_ID(id_set, readfile, outfile):
	ori_reads = []
	count = 0
	con_count = 0
	with open(readfile, "r") as r1:
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
			tmpid = re.sub('@', '', block[0].split(" ")[0])
			# print(tmpid)
			if tmpid not in id_set:
				con_count += 1
				continue
			ori_reads.append(block)

	print(count * 4, "lines accessed", con_count * 4, "reads skipped", len(ori_reads), "reads to write")

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
	# exit(1)
	with mp.Pool(2) as pool:
		lparam = [(id_set, sys.argv[2], sys.argv[4] + "remained_reduced_R1.fastq"),
				  (id_set, sys.argv[3], sys.argv[4] + "remained_reduced_R2.fastq")]
		pool.starmap(remove_by_ID, lparam)


def remove_by_ID(id_set, readfile, outfile):
	remained_reads = []
	with open(readfile, "r") as r1f:
		while True:
			block = list(its.islice(r1f, 4))

			if not block:
				break
			tmpid = re.sub('@', '', block[0].split(" ")[0])
			# print(tmpid)
			if tmpid in id_set:
				continue
			remained_reads.append(block)
	print(len(remained_reads), "reads remained")
	with open(outfile, "w+") as r1of:
		for lines in remained_reads:
			for line in lines:
				r1of.write(line)


def get_contigs():
	with open(sys.argv[1], "r") as f:
		for lines in iter(lambda: list(its.islice(f, 2)), []):
			# print(lines)
			filename = "contig_" + lines[0].split(" ")[3].strip()

			with open(filename, "w+") as wf:
				for line in lines:
					wf.write(line)


def get_rev_comp(input_file, revcomp=False):
	print("getting reverse complement of ", input_file)
	with open(input_file, "r") as f1:

		seq1 = f1.read()
		if seq1[0] == ">":
			f1.seek(0)
			lines = f1.readlines()[1:]
			seq1 = ""
			for l in lines:
				seq1 += l.replace("\n", "").replace(" ", "")
		else:
			seq1 = seq1.replace("\n", "").replace(" ", "")
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
	with open("rev_comp_" + os.path.basename(input_file), "w+") as f2:
		f2.write(">reverse_complement_of_" + input_file + "\n")
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
	with open(sys.argv[1], "r") as f1, open("rev_comp_" + os.path.basename(sys.argv[1]), "w+") as f2:
		for line in f1:
			if count % 2 == 0:
				seq = line.strip()

				f2.write(rev_comp_read(seq) + "\n")
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
	with open(sys.argv[2], "r") as rnf:
		for line in rnf:
			fields = line.split(" ")
			if fields[0] in id_set:
				flag = format(fields[1].strip(), 'b')[::-1]
				if flag[4] == "1":
					rc_read_set.add(rev_comp_read(fields[3].strip()))
					rc_dup_read_IDs.update({rev_comp_read(fields[3].strip()): {line.split(" ")[0].strip()}})
				else:
					read_set.add(fields[3].strip())
					dup_read_IDs.update({fields[3].strip(): {line.split(" ")[0].strip()}})
	with mp.Pool(2) as pool:
		lparam = [(read_set, dup_read_IDs, sys.argv[2], "dup_reads_R1.txt"),
				  (read_set, dup_read_IDs, sys.argv[3], "dup_reads_R2.txt")]
		pool.starmap(search_by_read, lparam)


def search_by_read(read_set, dup_read_IDs, readfile, outfile):
	with open(readfile, "r") as r1f:
		while True:
			block = list(its.islice(r1f, 4))

			if not block:
				break
			tmpid = re.sub('@', '', block[0].split(" ")[0])
			tmpread = block[1].strip()
			# print(tmpid)
			if tmpread not in read_set:
				continue
			dup_read_IDs[tmpread].add(tmpid)
		# print(len(remained_reads),"reads remained")
		# return dup_read_IDs
		with open(outfile, "w+") as r1of:
			for lines in dup_read_IDs:
				for line in lines:
					r1of.write(line)


def split_seq(input_file, length=9050, step=1000):
	seq = ""
	with open(input_file, "r") as f:
		for line in f:
			if line[0] != ">":
				seq += line.strip()
	for i in range(0, len(seq) - length, step):
		with open("slice_" + str(i) + "_" + str(i + length) + "_" + os.path.basename(input_file), "w+") as wf:
			wf.write(">index" + str(i) + "_" + str(i + length) + "\n")
			wf.write(seq[i:i + length])


def count_pair_dist(input_file):
	dist = {}
	with open(input_file, "r") as f:
		for line in f:
			sl = line.split(" ")
			if sl[0] in dist.keys():
				dist[sl[0]].append(int(sl[2]))
			else:
				dist.update({sl[0]: [int(sl[2])]})
	dist_values = []
	for i in dist.keys():
		dist_values.append(abs(dist[i][0] - dist[i][1]))
	dist_val_count = Counter(dist_values)
	print(len(dist_values), len(dist), dist_val_count, max([x for x in dist_val_count.values()]))
	print(statistics.mean(dist_values), statistics.median(dist_values))


def test_shell(readfile1, readfile2):
	verify_sub_command = os.path.dirname(
		__file__) + "/find_sub.sh" + " " + "-r" + " " + "oneline_NC045512.2.fasta_output/mt_read_subbed_0.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"
	# print(verify_sub_command)
	# verify_sub_command= shlex.split(verify_sub_command)

	verify_proc = subprocess.Popen(verify_sub_command,
								   stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
	# verify_proc = subprocess.run(os.path.dirname(__file__)+" count_error.sh 1 2 3",shell=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
	for line in verify_proc.stdout:
		print(line, end="")
	verify_proc.stdout.close()
	proc_code = verify_proc.wait()
	if proc_code:
		raise subprocess.CalledProcessError(proc_code, verify_sub_command)
	print(verify_proc.args)
	print(verify_proc.stdout)


def get_gap_reads(ref_file, strain, samfile, readfile1, readfile2):
	"""produce and verify new strain
	@param ref_file fasta file of reference sequence
	@param strain integer number of strain
	"""

	# strain = 2
	if strain < 0:
		print("invalid strain number, error")
		exit(-2)
	prev_readset = set({})
	if strain > 0:
		for i in range(0, strain):
			print("getting reads in subbed_reads_" + str(i))
			with open("subbed_reads_" + str(i) + ".sam", "r") as subf:
				for line in subf:
					sline = line.strip().split(" ")
					prev_readset.add(sline[st_find.read_field])
			with open("rejected_subbed_reads_" + str(i) + ".sam", "r") as rejf:
				for line in rejf:
					fields = line.strip().split(" ")
					prev_readset.add(fields[st_find.read_field])
			print(len(prev_readset))
	readlist = []
	with open(samfile, "r") as samf:
		for line in samf:
			if line.split(" ")[st_find.read_field] not in prev_readset:
				fields = line.strip().split(" ")
				if "N" in fields[st_find.read_field]:
					continue
				readlist.append(fields)
	readlist = fix_s_pos(readlist)
	print(len(readlist), "candidate reads")
	ref_seq = ""
	with open(ref_file, "r") as rf:
		for line in rf:
			if line[0] != ">":
				ref_seq += line.strip()
	ref = ref_seq

	for ir, read in enumerate(readlist):
		tmp_misP = []
		read_index = int(read[st_find.index_field]) - 1
		for i, base in enumerate(read[st_find.read_field]):
			if ref[read_index + i] != base:
				tmp_misP.append(read_index + i)
		for mp in tmp_misP:
			read.append((int(mp), ref[mp], read[st_find.read_field][mp - read_index]))

	batch = 0

	# while len(readlist) > 0:
	subbed_read = []
	misPs = []
	ref = ref_seq
	covered_pos = {}
	for batch, read in enumerate(readlist):
		overlap = False
		read_index = int(read[2]) - 1
		read_misPs = [int(x[0]) for x in read[st_find.misp_field:]]
		'''   
        if len(subbed_read) == 0:
            subbed_read.append(read)
            if batch == 0:
                #misPs.extend(read_misPs)
                misPs.extend(read[5:])
            else:
                curr_index = 0
                for new_np_index, new_mp in enumerate(read_misPs):
                    if new_mp <= int(misPs[0][0]):
                        #misPs.insert(0, new_mp)
                        misPs.insert(0,read[5+new_np_index])
                    else:
                        while curr_index < len(misPs):
                            if new_mp <= curr_index:
                                #misPs.insert(curr_index, new_mp)
                                misPs.insert(curr_index,read[5+new_np_index])
                                curr_index += 1
                                break
                            curr_index += 1'''
		# else:
		print("curr read", batch, read, "bases covered", len(covered_pos))
		if len(subbed_read) > 0:
			# read clash detection, find all overlapping reads

			# mp-clash, compare all MPs, if the reads only itroduce new MP outise current mp ranges
			'''   
            if read_misPs[0] > misPs[-1][0] or read_misPs[-1] < misPs[0][0]:
                subbed_read.append(read)
                curr_index = 0
                for new_np_index, new_mp in enumerate(read_misPs):
                    if new_mp <= int(misPs[0][0]):
                        # misPs.insert(0, new_mp)
                        misPs.insert(0, read[5 + new_np_index])
                    else:
                        while curr_index < len(misPs):
                            if new_mp <= curr_index:
                                # misPs.insert(curr_index, new_mp)
                                misPs.insert(curr_index, read[5 + new_np_index])
                                curr_index += 1
                                break
                            curr_index += 1

                            '''

			clash = False
			for rei, base in enumerate(read[st_find.read_field]):
				r_start = int(read[st_find.index_field]) - 1
				if r_start + rei in covered_pos.keys():
					if read[st_find.read_field][rei] != covered_pos[r_start + rei]:
						clash = True
						print("clash at", r_start + rei, covered_pos[r_start + rei], read[st_find.read_field][rei])
						break
				# for sr in subbed_read:
				#    sr_start = int(sr[2])-1
				#    if sr_start < read_misp[0] < sr_start + len(sr[3]):
				#        if read_misp[2] != sr[3][int(read_misp[0])-sr_start]:
				#            clash = True
				#            break

			if clash:
				batch += 1
				continue

			'''   
                    if int(misPs[0][0]) <= int(read_misp[0]) <= int(misPs[-1][0]) :
                        #introducing new misP that clashes with subbed read, discard current read
                        #if read_misp[0] not in misPs:
                         #   overlap=True
                         #   break
                        #else:
                        index_overlap = False
                        misp_overlap = False
                        for subbed_misp in misPs:
                            if  read_misp[0] == subbed_misp[0]:
                                print("overlap misp",read_misp,subbed_misp)
                                if read_misp[2] != subbed_misp[2]:
                                    misp_overlap = True
                                    break
                                else:
                                    index_overlap = True
                        if misp_overlap or not index_overlap:
                            overlap = True
                            break

                if overlap:
                    batch += 1
                    continue'''

			# else:
			# subbed_read.append(read)
			# misPs.extend(read_misPs)
			''' 
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
                            break'''

		temp_ref = ref[:read_index] + read[3] + ref[read_index + len(read[3]):]
		# print(editdistance.eval(temp_ref,ref))
		with open("batch_" + str(batch) + "_reference.fa", "w+") as bf:
			bf.write(">batch_" + str(batch) + "\n")
			bf.write(temp_ref)
		verify_sub_command = os.path.dirname(
			__file__) + "/find_sub.sh" + " " + "-r" + " " + "batch_" + str(
			batch) + "_reference.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"

		verify_proc = subprocess.run(verify_sub_command,
									 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
									 shell=True)
		print(verify_proc.stdout.split("\n")[-5:])
		gap_list = verify_proc.stdout.split("\n")[-3]
		if gap_list != "no gaps":
			# print(read)
			# print(gap_list)
			gap_list = re.sub('[\[\]]', '', gap_list)
			print(gap_list.split(",")[:-1])

		# gap_list = gap_list.split(",")[:-1]
		# ref = ref[:read_index] + ref_seq[read_index:read_index+len(read[3])] + ref[read_index + len(read[3]):]
		# subbed_read.remove(read)
		else:

			subbed_read.append(read)
			# update covered positions with base in covered_pos
			read_start = int(read[2]) - 1
			for ind, base in enumerate(read[3]):

				if read_start + ind in covered_pos.keys():

					if covered_pos[read_start + ind] != base:
						print("clash happened for current read", ind, covered_pos[read_start + ind], base)
						exit()
				else:

					covered_pos.update({read_start + ind: base})
			if batch == 0:
				# misPs.extend(read_misPs)
				misPs.extend(read[st_find.misp_field:])
			else:
				misPs.extend(read[st_find.misp_field:])
				misPs = sorted(misPs, key=lambda x: int(x[0]))
				curr_index = 0
				for new_np_index, new_mp in enumerate(read_misPs):
					if new_mp <= int(misPs[0][0]):
						# misPs.insert(0, new_mp)
						misPs.insert(0, read[st_find.misp_field + new_np_index])
					elif new_mp >= int(misPs[-1][0]):
						misPs.insert(-1, read[st_find.misp_field + new_np_index])
					else:
						while curr_index < len(misPs):
							if new_mp <= curr_index:
								# misPs.insert(curr_index, new_mp)
								misPs.insert(curr_index, read[6 + new_np_index])
								curr_index += 1
								break
							curr_index += 1

			ref = temp_ref
			print(batch, "no gaps for curr read\n")
		# print(len(covered_pos),covered_pos.keys())
		# print(subbed_read,len(readlist))
		batch += 1

	with open("final_strain_" + str(strain) + "_reference.fa", "w+") as bf:
		bf.write(">final_strain_" + str(strain) + "\n")
		bf.write(ref)
	with open("subbed_reads_" + str(strain) + ".sam", "w+") as bf:
		for line in subbed_read:
			bf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
	# subprocess.run("rm -r batch_*",shell=True)
	print(subbed_read)
	for read in subbed_read:
		if read in readlist:
			readlist.remove(read)

	with open("strain_" + str(strain) + "_spike.fa", "w+") as sf:
		sf.write(">strain_" + str(strain) + "_spike" + "\n")
		sf.write(ref[21562:25383])
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
	with open("test_output.txt", "r") as tf:
		for line in tf:
			out += line
	print(out.split("\n")[-2])


def count_fq_freq(fastqfile):
	readlist = []
	with open(fastqfile, "r") as fq:
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

	for k, v in read_freq.items():
		if v not in count.keys():
			count.update({v: 1})
		else:
			count[v] += 1

	# print(k,v)
	sum_freq = 0
	for k, v in {k: v for k, v in sorted(count.items(), key=lambda x: x[1])}.items():
		print(k, v)
		sum_freq += k * v
	print(sum_freq)


def remove_freq(samfile, freq):
	readfreq = {}
	readlist = []
	target_id = set({})
	with open(samfile, "r") as f:
		for line in f:
			fields = line.strip().split(" ")
			readlist.append(fields)
			# print(fields)
			if len(fields) > 0 and fields[3] not in readfreq.keys():
				readfreq.update({fields[3]: 1})
			# target_id.add(fields[0])
			else:
				readfreq[fields[3]] += 1
			# target_id.remove(fields[0])
	print(len(readlist), "reads in ", samfile, " uniq reads", len(readfreq))
	filtered_readlist = []
	for read in readlist:
		if readfreq[read[3]] <= freq:
			target_id.add(read[0])
		else:
			filtered_readlist.append(read)
	print(len(filtered_readlist), "reads with frequency > ", freq)
	lastpos = -1
	gapped_index = set({})
	with open("freq" + str(freq) + "plus_" + os.path.basename(samfile), "w+") as wf:

		for read in filtered_readlist:
			wf.write(read[0] + " " + str(read[1]) + " " + str(read[2]) + " " + read[3] + " " + read[4] + "\n")
			if lastpos > 0 and lastpos < int(read[2]):
				if int(read[2]) not in gapped_index:
					# gapped_index.add(int(read[2]))
					print("gap between ", prev, read)
				# print(read)
			if read[0] not in target_id:
				if lastpos == -1:
					print("first read", read)
				prev = read
				lastpos = int(read[2]) + len(read[3])
	# print(len(gapped_index), " reads have gaps between it and the read before ")
	# print(gapped_index)
	if lastpos < 29800:
		print("last read end at", lastpos)
# print(count ,"reads have frequency 1")


def calc_coverage(ref_file, sam_file):
	'''
    verify_sub_command = os.path.dirname(
        __file__) + "/find_sub.sh" + " " + "-r" + " " + ref_file + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"

    verify_proc = subprocess.run(verify_sub_command,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                 shell=True)
    '''
	readlist = bm.read_sam(sam_file)
	matrix = bm.matrix_from_readlist(readlist, 0, set({}))

	ref = ""
	with open(ref_file, 'r') as refg:
		for line in refg:
			if ">" not in line:
				ref += line.strip()

	real_narrowed, paired_real_narrowed, nearly_real_narrowed, potential_mutated = bm.narrow_reads(ref,
																								   matrix.narrowed_read,
																								   "", True)

	rigor_matrix = bm.matrix_from_readlist(paired_real_narrowed, 0.975, set({}), False,
										   matrix, "real_narrowed")

	cvg_list = []
	rn_cvg_list = []

	for i in range(0, matrix.narrowed_matrix.shape[1]):
		tmp = np.squeeze(matrix.narrowed_matrix.getcol(i).toarray())
		tmp_count = np.bincount(tmp)[1:]
		cvg_list.append(sum(tmp_count))

	for i in range(0, rigor_matrix.real_narrowed_matrix.shape[1]):
		tmp = np.squeeze(rigor_matrix.real_narrowed_matrix.getcol(i).toarray())
		tmp_count = np.bincount(tmp)[1:]
		rn_cvg_list.append(sum(tmp_count))
	# with open("cvg.txt", "w+") as ncf:
	return cvg_list, rn_cvg_list


def get_coverage(ref_file, sam_file):
	cvg_list, rn_cvg_list = calc_coverage(ref_file, sam_file)
	with open("cvg.txt", "w+") as ncf, open("rn_cvg.txt", "w+") as rncf:
		for i, v in enumerate(cvg_list):
			ncf.write(str(i + 1) + ": " + str(v) + ", ")
		for i1, v1 in enumerate(rn_cvg_list):
			rncf.write(str(i1 + 1) + ": " + str(v1) + ", ")
# draw_coverage(["cvg.txt"])


def draw_coverage(files):
	# order = ["original SARS-CoV-2","discovered_strain_1","discovered_strain_2"]
	order = ["SARS-CoV-2/WIV04", "discovered strain 1", "discovered strain 2"]
	avg = 0
	for file, oname in zip(files, order):
		coverages = {}
		xaxis = []
		yaxis = []

		with open(file, "r") as cf:
			for line in cf:
				cov_list = line.split(", ")
				for cov in cov_list:
					nums = re.search('([0-9]+)\: ([0-9]+)', cov)
					if nums is None:
						break
					xaxis.append(int(nums.group(1)) + 1)
					yaxis.append(int(nums.group(2)))
					coverages.update({nums.group(1): nums.group(2)})
			avg = round(statistics.mean(yaxis), 2)
		x_y_spline = scipy.interpolate.make_interp_spline(xaxis, yaxis)
		x_ = np.linspace(min(xaxis), max(xaxis), 500)
		y_ = x_y_spline(x_)

		plt.plot(xaxis, yaxis, linewidth=2, label=oname + " Average:" + str(avg))
	legend = plt.legend(loc='upper center', shadow=True, fontsize='x-large')
	legend.get_frame()
	plt.xlabel("Position number ")
	plt.ylabel("Sequencing depth")
	plt.show()


def diff_coverage(cvg1, cvg2):
	cvg1_dict = {}
	cvg2_dict = {}
	with open(cvg1, "r") as f1:
		for line in f1:
			nums = line.split(", ")
			for pos in nums:
				if pos == '':
					continue
				kv = pos.split(":")
				cvg1_dict.update({int(kv[0]): int(kv[1])})

	with open(cvg2, "r") as f1:
		for line in f1:
			nums = line.split(", ")
			for pos in nums:
				if pos == '':
					continue
				kv = pos.split(":")
				cvg2_dict.update({int(kv[0]): int(kv[1])})
	print(len(cvg1_dict), len(cvg2_dict))
	diff_cvg = {}
	for k, v in cvg1_dict.items():
		if k in cvg2_dict.keys():
			diff_cvg.update({k: cvg1_dict[k] - cvg2_dict[k]})

	print(len(diff_cvg))

	with open("diff_cvg.txt", "w+") as wf:
		for k, v in diff_cvg.items():
			wf.write(str(k) + ":" + str(diff_cvg[k]) + ", ")


def count_subbed_freq(freq_file, sub_file):
	freq = {}
	with open(freq_file, "r") as ff:
		for line in ff:
			fields = line.strip().split(": ")
			freq.update({fields[0]: fields[1]})
	readlist = []
	with open(sub_file, "r") as sf:
		for line in sf:
			fields = line.strip().split(" ")
			fields.append(freq[fields[3]])
			readlist.append(fields)
	with open(sub_file + "freq.txt", "w+") as wf:
		for line in readlist:
			wf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + " " + str(
				line[5]) + "\n")


def filter_sam_freq(sam):
	freq = {}
	readlist = []
	with open(sam, "r") as f:
		for line in f:
			fields = line.strip().split(" ")
			readlist.append(fields)
			if fields[3] not in freq.keys():
				freq.update({fields[3]: 1})
			else:
				freq[fields[3]] += 1
	filtered_readlist = []
	for read in readlist:
		if freq[read[3]] > 1:
			filtered_readlist.append(read)

	with open("read_freq.txt", "w+") as ff:
		for k, v in sorted(freq.items(), key=lambda item: item[1]):
			ff.write(k + ": " + str(v) + "\n")

	with open("filtered_freq1_extract.sam", "w+") as wf:
		for line in filtered_readlist:
			wf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + " " + str(
				freq[line[3]]) + "\n")


def check_sam_gap(sam):
	readlist = []
	with open(sam, "r") as f:
		for line in f:
			fields = line.strip().split(" ")
			fields[2] = int(fields[2])
			readlist.append(fields)
	prev = int(readlist[0][2]) - 1
	prevlen = len(readlist[0][3])
	test_gaps = []
	for line in readlist:
		if int(line[2]) - prev > prevlen:
			# print("gaps at ", line, prev)
			test_gaps.extend([x for x in range(prev + prevlen - 1, int(line[2]) - 1)])
		prev = int(line[2])
		prevlen = len(line[3])
	return test_gaps
	matrix_info = bm.matrix_from_readlist(readlist, 1, set({}), True)
	gaps = []
	for i in range(20, matrix_info.narrowed_matrix.shape[1]):
		tmp = np.squeeze(matrix_info.narrowed_matrix.getcol(i).toarray())

		tmp_count = np.bincount(tmp)[1:]
		if sum(tmp_count) == 0:
			# print("gap col at ",i)

			# print("no gaps for ",sam)
			gaps.append(i)

	print(len(gaps), "gaps", "cols for checking gaps", 0, matrix_info.narrowed_matrix.shape[1])

	if len(gaps) > 0:

		with open("checked_gaps.txt", "w+") as gapf:
			for gap_col in gaps:
				gapf.write(str(gap_col) + ",")


def get_misp(ref_file, mreads_list, printing=True):
	"""

	:param ref_file: reference sequence file
	:param mreads_list: sam file list
	:param printing: True to print to stdout
	:return: misPs {position: ["ref_base|misp_base"]}, misP_source {position:[SRR_ID]}, misP_reads {position:[reads]}
	"""
	ref = ""
	misPs = {}
	misPs_count = {}
	misPs_source = {}
	misP_reads = {}
	pos = {}

	with open(ref_file, "r") as rf:
		for line in rf:
			if line[0] != ">":
				ref += line.strip()

	for mreads in mreads_list:
		with open(mreads, "r") as mf:

			for line in mf:
				fields = line.strip().split(" ")

				start = int(fields[2]) - 1

				for i, base in enumerate(fields[3]):
					if start + i >= len(ref):
						print(len(ref), start, i, fields)
						break
					if base != ref[start + i] and base != "N":

						if start + i not in misPs.keys():

							misPs.update({start + i: [ref[start + i] + "|" + base]})
							misPs_count.update({start + i: [1]})
							misPs_source.update({start + i: {base: [fields[0].split(".")[0]]}})
							misP_reads.update({start + i: {base: [fields]}})
						else:
							existed = False
							for i_mp, mp in enumerate(misPs[start + i]):
								if base == mp[-1]:
									existed = True
									misPs_count[start + i][i_mp] += 1
									misPs_source[start + i][base].append(fields[0].split(".")[0])
									misP_reads[start + i][base].append(fields)
									# print(start+i,mp,misP_reads[start+i][base])
									break

							if not existed:
								misPs[start + i].append(ref[start + i] + "|" + base)
								misPs_count[start + i].append(1)
								misPs_source[start + i].update({base: [fields[0].split(".")[0]]})

								misP_reads[start + i].update({base: [fields]})
						if start + i not in pos.keys():

							pos.update({start + i: [fields]})
						else:
							pos[start + i].append(fields)

	# with open("different_bases.txt","w+") as wf:

	for k, v in sorted(pos.items(), key=lambda x: x[0]):

		# if 21562<=k<=25384:

		for change in misPs[k]:
			support_list = [(k, v) for k, v in Counter(misPs_source[k][change[-1]]).items() if v > 1]
			multi = False
			for support_count in support_list:
				if support_count[1] > 1:
					multi = True
					break
			if len(support_list) > 1 or multi:
				if printing:
					print(k, ": ", end="")
					print(change, support_list, end=",    ")
					print()

			# wf.write(str(k)+":"+v+", ")

	for k, v in sorted(misPs.items(), key=lambda x: x[0]):
		print(k, v)
	return misPs, misPs_source, misP_reads


def overwrite_misp(ref, file62, mreads):
	misPs, misPs_source, misp_reads = get_misp(ref, mreads)
	nearly_narrowed_62 = []
	misp_count = 0
	with open(file62, "r") as f62:
		for line in f62:
			fields = line.strip().split(" ")
			fields.append([])
			start = int(fields[2]) - 1
			for i, base in enumerate(fields[3]):
				if start + i >= len(ref):
					print(start, i, fields)
				if base != ref[start + i] and base != "N":
					fields[-1].append((start + i, base))  # ref[start+i]+"|"+
					misp_count += 1
				nearly_narrowed_62.append(fields)
	# misPs_62,misPs_source_62 = get_misp(ref,file62)
	singleton = []
	multi = []
	print(len(nearly_narrowed_62), "reads and ", misp_count, "misPs in 62 file")
	for ri, read in enumerate(nearly_narrowed_62):
		found = False
		for tmp_misp in read[-1]:
			# overwritten
			if int(nearly_narrowed_62[ri + 1][2]) - 1 <= int(tmp_misp[0]) <= int(
					nearly_narrowed_62[ri + 1][2]) - 1 + len(nearly_narrowed_62[ri + 1][3]):
				print("overwrite", tmp_misp, read, nearly_narrowed_62[ri + 1])

			'''   
            if misPs_source[tmp_misp[0][tmp_misp[1]]] > 1:
                multi.append(tmp_misp)
            else:
                singleton.append(tmp_misp)
            '''


def fix_s_pos(subbed_reads):
	new_subbed_reads = []
	for ir, fields in enumerate(subbed_reads):
		index = int(fields[2]) - 1
		cigar = fields[4]
		cigar_str = re.findall(r"[0-9]+[MIDSH]", cigar)
		blk_pos = []
		blk_type = []
		ini = 0
		tmp_length = 0
		base_length = 0
		matched = 0
		for block in cigar_str:
			m = re.search(r'([0-9]+)([MIDSH])', block)
			bl = int(m.group(1)) + ini
			bt = str(m.group(2))
			# if bt == "S" or bt == "H":
			#    continue
			blk_type.append(bt)

			blk_pos.append(bl)
			ini = bl
		# if iread[4]=="140M2D10M":
		#	print(block, int(m.group(1)))

		c = 0
		inserts = []
		begin = 0
		reduce = 0  # deduct "non-existent" bases in total length, "D" and "H" not shown in reads
		curr_pos = 0

		if blk_type[0] == "S":
			if index - blk_pos[0] < 0:

				start_pos = blk_pos[0] - index
				tmp_length -= start_pos
				index = 0
			else:
				index = index - blk_pos[0]
			subbed_reads[ir][2] = index + 1
	return subbed_reads


dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))


def verify_misp(ref_file, samfile, mreads, new_strain_file, candidate_sam):
	ref = ""
	with open(new_strain_file, "r") as f:
		for line in f:
			if line[0] != ">":
				temp = line.replace(" ", "").replace("\n", "").replace("\d", "").upper()
				ref += re.sub('\d', "", temp)

	ori_ref = ""
	with open(ref_file, "r") as f:
		for line in f:
			if line[0] != ">":
				temp = line.replace(" ", "").replace("\n", "").replace("\d", "").upper()
				ori_ref += re.sub('\d', "", temp)

	strain = int(re.search('[0-9]+', os.path.basename(mreads))[0])

	test_misp = {}
	subbed_read = []
	with open(mreads, "r") as subf:
		for line in subf:
			fields = line.strip().split(" ")
			subbed_read.append(fields)
			start = int(fields[2]) - 1
			for ir, base in enumerate(fields[3]):
				if base != ori_ref[start + ir]:
					test_misp.update({start + ir: ori_ref[start + ir] + "|" + base})

	misPs, misPs_source, misP_reads = get_misp(ref_file, [mreads], False)
	if len(test_misp.keys() - misPs.keys()) > 0:
		print(test_misp)
		print(misPs)
		print(test_misp.keys() - misPs.keys())
		exit()
	# set up deletion handling

	total_misps = sum([len(x) for k, x in misPs.items()])
	print(total_misps, "misps")
	# print(sorted(misPs.items(),key=lambda x:x[0])[:10])

	read_freq = {}
	candidate_read = []
	with open(candidate_sam, "r") as canf:
		for line in canf:
			sline = line.strip().split(" ")
			candidate_read.append(sline)
			read_freq.update({sline[3]: int(sline[5])})
	initial_matrix_info = bm.matrix_from_readlist(bm.read_sam(samfile), 0, marked_id=set({}), target="raw")

	raw_matrix_info = bm.build_insertion(initial_matrix_info, 0)
	# not using insertion now,later change to raw_matrix_info.real_narrowed_matrix
	matrix = raw_matrix_info.narrowed_matrix

	found_segments = set({})
	rejected_read = []
	# accepted_misps = set({})
	accepted_misps = {}
	for pos, clist in misPs.items():
		reject = False
		if pos > matrix.shape[1] or pos in accepted_misps.keys():
			# print("max index exceeds at ", pos, clist)
			continue
		tmp = np.squeeze(matrix.getcol(pos).toarray())
		tmp_count = np.bincount(tmp)[1:]
		for change in clist:

			support_index = np.where(tmp == dict_tonu[change[-1]])
			support_count = support_index[0].shape[0]

			if support_count > 0:
				for sp in support_index[0]:
					temp_read = raw_matrix_info.narrowed_read[sp]
					# print(temp_read)

					# print(tmp_count)
					# print(ref[pos-5:pos],ref[pos],ref[pos:pos+5])
					# print(change,pos,sp,support_index,support_count)
					read_index = int(temp_read[2]) - 1
					window = 0
					tolerance = 1
					curr_mcount = 0
					curr = pos
					read_end = read_index + len(temp_read[3])
					reached = False
					l_reached = False
					l_window = 0
					l_temp_read = temp_read

					r_reached = False
					r_window = 0
					r_temp_read = temp_read
					encountered_misp = []

					while curr - 1 >= read_index:
						if ref[curr - 1] == temp_read[3][curr - 1 - read_index]:
							l_window += 1

						else:
							l_reached = True
							break
						curr -= 1
					curr_l_line = sp - 1
					if curr_l_line >= 0 and curr - 1 <= raw_matrix_info.narrowed_read[curr_l_line][2] - 1 + len(
							raw_matrix_info.narrowed_read[curr_l_line][3]):
						curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
					# print(curr_l_read)

					else:
						l_reached = True
					while not l_reached and curr_l_line >= 0:
						if curr - 1 < int(curr_l_read[2]) - 1:

							if curr_l_line - 1 < 0:
								l_reached = True
								break

							curr_l_line -= 1
							curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
							if curr - 1 > int(curr_l_read[2]) - 1 + len(curr_l_read[3]):
								l_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_l_read[2]) - 1 <= pos <= int(curr_l_read[2]) - 1 + len(curr_l_read[3]) - 1 and \
									curr_l_read[3][pos - (int(curr_l_read[2]) - 1)] != change[-1]:
								# print("read not supporting misp",pos,change)
								# print(curr_l_read)
								# l_reached= True
								curr_l_line -= 1

								continue
							# print(curr_l_read)
							continue

						# try:
						if ref[curr - 1] == curr_l_read[3][curr - 1 - (int(curr_l_read[2]) - 1)]:
							l_window += 1

						else:
							l_reached = True
							break

						# except:
						#    print(curr_l_read,curr-1-(int(curr_l_read[2]-1)))
						#    exit()
						curr -= 1
					l_seg = ref[curr:pos]
					curr = pos

					while curr + 1 <= read_index + len(temp_read[3]) - 1:

						if ref[curr + 1] == temp_read[3][curr + 1 - read_index]:
							r_window += 1
						else:
							r_reached = False
							break

						curr += 1
					curr_r_line = sp + 1
					if curr_r_line <= raw_matrix_info.narrowed_matrix.shape[0] and \
							raw_matrix_info.narrowed_read[curr_r_line][2] - 1 <= curr + 1:

						curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
					# print(curr_r_read)
					else:
						r_reached = True
					while not r_reached and curr_r_line <= raw_matrix_info.narrowed_matrix.shape[0]:
						if curr + 1 > int(curr_r_read[2]) - 1 + len(curr_r_read[3]) - 1:

							if curr_r_line + 1 >= raw_matrix_info.narrowed_matrix.shape[0]:
								r_reached = True
								break
							curr_r_line += 1
							curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
							# if current read does not cover the curr position (gaps), break
							if curr + 1 < int(curr_r_read[2]) - 1:
								r_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_r_read[2]) - 1 <= pos <= int(curr_r_read[2]) - 1 + len(curr_r_read[3]) and \
									curr_r_read[3][pos - (int(curr_r_read[2]) - 1)] != change[-1]:
								# print("read not supporting misp",pos,change, curr_r_read)
								# print(curr_r_read)
								curr_r_line += 1
								# r_reached= True
								continue

							# print(curr_r_read)
							continue
						if ref[curr + 1] == curr_r_read[3][curr + 1 - (int(curr_r_read[2]) - 1)]:
							r_window += 1
						else:
							r_reached = True
							break
						curr += 1
					seg = l_seg + ref[pos:curr + 1]
					window = l_window + r_window + 1
					if seg not in found_segments:
						found_segments.add(seg)
						print("target misp", pos, read_freq.get(misP_reads[pos][change[-1]][0][3]))

						print("seg_length", len(seg))

					# print("ref:",ori_ref[pos-l_window:pos+r_window+1])

					# print("seg:",seg)
					# reject this read candidate
					if len(seg) < 50 and len(misP_reads[pos][change[-1]]) < 2 and read_freq.get(
							misP_reads[pos][change[-1]][0][3]) == 1:
						rejected_read.extend(misP_reads[pos][change[-1]])
						# print("rejected ",misP_reads[pos][change[-1]])
						reject = True
						break
					# continue

					'''
                    DP find longest segment for a tolerance
                    left_mis = []
                    right_mis = []

                    while curr >= read_index + 1:

                        if ref[curr - 1] != temp_read[3][curr - 1 - read_index]:
                            left_mis.append(curr-1)

                        curr -= 1
                    if len(left_mis) == 0:
                        left_mis.append(-1)
                    curr = pos

                    while curr <= read_index + len(temp_read[3]) - 2:
                        if ref[curr + 1] != temp_read[3][curr + 1 - read_index]:
                            right_mis.append(curr+1)
                        curr += 1
                    if len(right_mis) == 0:
                        if left_mis[0] == -1:
                            window = len(temp_read[3])
                            reached = True
                        right_mis.append(-1)
                    print("misp position at ", pos,change)
                    #print(left_mis,right_mis)
                    if not reached:
                        #calculating longest window given a tolerance
                        length_matrix = [[-1]*len(right_mis) for i in range(len(left_mis))]
                        for i in range(0,len(left_mis)):
                            if left_mis[0] == -1:
                                length_matrix[i][0] = pos - read_index
                            else:
                                if curr_mcount <= tolerance:

                                    length_matrix[i][0] = pos-left_mis[i]
                                    curr_mcount += 1
                                else:

                                    length_matrix[i][0] = length_matrix[i-1][0]
                        i = 0
                        curr_mcount = 0
                        for i in range(0,len(right_mis)):
                            if right_mis[0] == -1:
                                length_matrix[0][i] = read_end-curr
                            else:
                                if curr_mcount <= tolerance:
                                    length_matrix[0][i] = right_mis[i] - pos
                                    curr_mcount += 1
                                else:

                                    length_matrix[0][i] = length_matrix[0][i-1]
                        curr_mcount = 0
                        window = 0
                        if len(left_mis)>1 or len(right_mis) > 1:
                            i = j =0
                            for i in range(1,len(left_mis)):
                                for j in range(1, len(right_mis)):
                                    if curr_mcount >= tolerance:

                                        length_matrix[i][j] = max(length_matrix[i-1][j],length_matrix[i][j-1])
                                    else:
                                        if curr_mcount + 1 <= tolerance:
                                            move_left = window + left_mis[i]-left_mis[i-1]
                                            move_right = window + right_mis[j]-right_mis[j-1]

                                            move = max(move_left,move_right)
                                            length_matrix[i][j] = move
                                            #if move == move_right or move == move_left:
                                            curr_mcount += 1
                        window = length_matrix[len(left_mis)-1][len(right_mis)-1]
                    '''

					print("\n segment size", window, "\n")
					if len(misP_reads[pos][change[-1]]) > 1:
						total_freq = 0
						for tr in misP_reads[pos][change[-1]]:
							total_freq += read_freq.get(tr[3])
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][change[-1]][0][3]), window))
						accepted_misps.update({pos: (change, total_freq, window)})
					else:
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][change[-1]][0][3]), window))
						accepted_misps.update({pos: (change, read_freq.get(misP_reads[pos][change[-1]][0][3]), window)})

			else:
				if len(misP_reads[pos][change[-1]]) < 2 and read_freq.get(
						misP_reads[pos][change[-1]][0][3]) == 1:
					rejected_read.extend(misP_reads[pos][change[-1]])
					# print("rejected ",misP_reads[pos][change[-1]])
					reject = True
					break
				else:
					if len(misP_reads[pos][change[-1]]) > 1:
						total_freq = 0
						for tr in misP_reads[pos][change[-1]]:
							total_freq += read_freq.get(tr[3])
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][change[-1]][0][3]), window))
						accepted_misps.update({pos: (change, total_freq, 0)})
					else:
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][change[-1]][0][3]), window))
						accepted_misps.update(
							{pos: (change, read_freq.get(misP_reads[pos][change[-1]][0][3]), 0)})

		if reject:
			continue

	if len(rejected_read) > 0:
		print(len(rejected_read), "reads rejected")
		write_sam(rejected_read, "rejected_" + os.path.basename(mreads), "a")
		rejected_set = set([x[3] for x in rejected_read])

		d_info = {}
		# for ir,read in enumerate(subbed_read):
		#    if "D" in read[4]:

		for ir, read in enumerate(subbed_read):
			if read[3] in rejected_set:
				print("revert", read)
				start = int(read[2]) - 1
				ref = ref[:start] + ori_ref[start:start + len(read[3])] + ref[start + len(read[3]):]
			# ref[start:start+len(read[3])] = ori_ref[start:start+len(read[3])]

			# print(pos,change,support_count,tmp_count)

		# remove rejected reads from  subbed_reads and file
		new_subbed_read = []
		for read in subbed_read:
			if read[3] not in rejected_set:
				new_subbed_read.append(read)
		write_sam(new_subbed_read, mreads)
		print(len(subbed_read), "reads originally in ", os.path.basename(mreads))
		subbed_read = new_subbed_read
		print(len(rejected_set), "reads reverted, removed from", os.path.basename(mreads), "now ", len(subbed_read))
		subbed_read_set = set([x[3] for x in subbed_read])
		reach_end = curr_gap_reads(ref, strain, subbed_read_set, candidate_read, ori_ref,
								   raw_matrix_info.narrowed_matrix, misP_reads, read_freq)
	# if reach_end:

	else:

		with open("misp_" + os.path.basename(mreads), "w+") as mf:
			for misp_pos in sorted(list(accepted_misps.keys())):
				val = accepted_misps[misp_pos]
				bases = val[0].split("|")
				print(str(misp_pos+1), "&", bases[0], "&", bases[1], "&", val[1], "&", val[2], "\\\\\n\hline")
				mf.write(str(misp_pos+1) + " & " + str(bases[0]) + " & " + bases[1] + "&" + str(val[1]) + " & " + str(
					val[2]) + "\\\\\n\hline\n")
	# return 0
	return rejected_read
	print()
	cvg_list = []
	for i in range(raw_matrix_info.real_narrowed_matrix.shape[1]):
		tmp = np.squeeze(raw_matrix_info.real_narrowed_matrix.getcol(i).toarray())
		tmp_count = np.bincount(tmp)[1:]
		cvg_list.append(sum(tmp_count))

	with open("cvg.txt", "w+") as near_rncf:
		for i, v in enumerate(cvg_list):
			near_rncf.write(str(i) + ": " + str(v) + ", ")


# facilitate verify_misp to stop at full verified reads
def fac_verify_misp(ref_file, samfile, mreads, new_strain_file, candidate_sam):
	strain = int(re.search('[0-9]+', os.path.basename(mreads))[0])
	with open("rejected_" + os.path.basename(mreads), "w+") as rejf:
		rejf.write("")
	rejected_reads = verify_misp(ref_file, samfile, mreads, new_strain_file, candidate_sam)
	tmp_rejected = []

	while (len(tmp_rejected) != len(rejected_reads)):
		tmp_rejected = copy.deepcopy(rejected_reads)

		rejected_reads = verify_misp(ref_file, samfile, mreads, new_strain_file, candidate_sam)

		verify_sub_command = "for f in {56..64}; do\n" + \
							 os.path.dirname(__file__) + "/find_sub.sh" + " " + "-r" + " " + "final_strain_" + str(
			strain) + "_reference.fa" + " " + "-1" + " " + "half_real_R1.fastq" + " " + "-2" + " " + "half_real_R2.fastq" + \
							 " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y;" + \
							 "cat " + os.path.dirname(
			__file__) + "/ mullti_support /${f}_out / round_1_extract.sam >> strain" + str(
			strain) + "_combine_extract.sam\ndone"

		verify_proc = subprocess.run(verify_sub_command,
									 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
									 shell=True)

#    print(rejected_reads)
# if rejected_reads is int:
#    break


'''
Fill more candidate reads after rejecting reads
'''


def curr_gap_reads(ref, strain, subbed_read_set, candidate_read, old_ref, matrix, misP_reads, read_freq):
	misp_bit = 6
	readfile1 = "half_real_R1.fastq"
	readfile2 = "half_real_R2.fastq"

	prev_readset = set({})
	covered_pos = {}
	for i in range(0, strain + 1):
		print("getting replacement reads in subbed_reads_" + str(i))
		with open("subbed_reads_" + str(i) + ".sam", "r") as subf:
			for line in subf:
				sline = line.strip().split(" ")
				prev_readset.add(sline[3])
		with open("rejected_subbed_reads_" + str(i) + ".sam", "r") as rejf:
			for line in rejf:
				fields = line.strip().split(" ")
				prev_readset.add(fields[3])
	with open("subbed_reads_" + str(strain) + ".sam", "r") as subf:
		for line in subf:
			fields = line.strip().split(" ")
			start = int(fields[2]) - 1
			for ib, base in enumerate(fields[3]):
				covered_pos.update({start + ib: base})

	print(len(prev_readset))

	readlist = candidate_read

	for ir, read in enumerate(readlist):
		tmp_misP = []
		read_index = int(read[2]) - 1
		for i, base in enumerate(read[3]):
			if old_ref[read_index + i] != base:
				tmp_misP.append(read_index + i)
		if len(tmp_misP) == 0:
			print(read)
			print(ref[read_index:read_index + len(read[3])])
			print(read[3])
			exit()
		for mp in tmp_misP:
			read.append((int(mp), ref[mp], read[3][mp - read_index]))

	subbed_read = []
	misPs = []
	for ir, read in enumerate(readlist):
		if read[3] in subbed_read_set:
			subbed_read.append(read)
			misPs.extend(read[misp_bit:])

	count = 0
	for batch, read in enumerate(readlist):
		reject = False
		count += 1
		if read[3] in prev_readset:
			# batch += 1
			continue
		read_index = int(read[2]) - 1
		read_misPs = [int(x[0]) for x in read[misp_bit:]]

		print("curr read", batch, read, "misp range", misPs[0], misPs[-1])
		clash = False
		for rei, base in enumerate(read[3]):
			r_start = int(read[2]) - 1
			if r_start + rei in covered_pos.keys():
				if read[3][rei] != covered_pos[r_start + rei]:
					clash = True
					print("clash at", r_start + rei, covered_pos[r_start + rei], read[3][rei])
					break

		if clash:
			batch += 1
			continue

		for rm in read[misp_bit:]:
			pos = rm[0]
			change = rm[2]
			tmp = np.squeeze(matrix.getcol(pos).toarray())
			tmp_count = np.bincount(tmp)[1:]

			support_index = np.where(tmp == dict_tonu[change[-1]])
			support_count = support_index[0].shape[0]
			if support_count == 0 and read_freq.get(read[3]) == 1 and pos not in misP_reads.keys():
				reject = True
				break
			# print("rejected ",misP_reads[pos][change[-1]])

		if reject:
			continue
		temp_ref = ref[:read_index] + read[3] + ref[read_index + len(read[3]):]
		# print(editdistance.eval(temp_ref,ref))
		with open("batch_" + str(batch) + "_reference.fa", "w+") as bf:
			bf.write(">batch_" + str(batch) + "\n")
			bf.write(temp_ref)
		verify_sub_command = os.path.dirname(
			__file__) + "/find_sub.sh" + " " + "-r" + " " + "batch_" + str(
			batch) + "_reference.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"

		verify_proc = subprocess.run(verify_sub_command,
									 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
									 shell=True)
		# print(verify_proc.stdout.split("\n")[-5:])
		gap_list = verify_proc.stdout.split("\n")[-3]
		if gap_list != "no gaps":

			print(gap_list)
			gap_list = re.sub('[\[\]]', '', gap_list)
			print(gap_list.split(",")[:-1])
			print(read)
		# gap_list = gap_list.split(",")[:-1]
		# ref = ref[:read_index] + ref_seq[read_index:read_index+len(read[3])] + ref[read_index + len(read[3]):]
		# subbed_read.remove(read)
		else:
			subbed_read.append(read)
			read_start = int(read[2]) - 1
			for ind, base in enumerate(read[3]):

				if read_start + ind in covered_pos.keys():

					if covered_pos[read_start + ind] != base:
						print("clash happened for current read", ind, covered_pos[read_start + ind], base)
						exit()
				else:

					covered_pos.update({read_start + ind: base})
			if batch == 0:
				# misPs.extend(read_misPs)
				misPs.extend(read[misp_bit:])
			else:
				curr_index = 0
				for new_np_index, new_mp in enumerate(read_misPs):
					if new_mp <= int(misPs[0][0]):
						# misPs.insert(0, new_mp)
						misPs.insert(0, read[misp_bit + new_np_index])
					else:
						while curr_index < len(misPs):
							if new_mp <= curr_index:
								# misPs.insert(curr_index, new_mp)
								misPs.insert(curr_index, read[misp_bit + new_np_index])
								curr_index += 1
								break
							curr_index += 1

			ref = temp_ref
			print(batch, "no gaps for ")
	# print(subbed_read,len(readlist))
	# batch += 1

	with open("final_strain_" + str(strain) + "_reference.fa", "w+") as bf:
		bf.write(">final_strain_" + str(strain) + "\n")
		bf.write(ref)
	with open("subbed_reads_" + str(strain) + ".sam", "w+") as bf:
		for line in subbed_read:
			bf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
	# subprocess.run("rm -r batch_*",shell=True)

	with open("fixed_strain_" + str(strain) + "_spike.fa", "w+") as sf:
		sf.write(">fixed_strain_" + str(strain) + "_spike" + "\n")
		sf.write(ref[21562:25383])
	strain += 1
	# remove temp files
	verify_sub_command = "rm" + " " + "-r" + " " + "./batch_*"

	verify_proc = subprocess.run(verify_sub_command,
								 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
								 shell=True)
	if count == len(readlist) - 1:
		return True
	else:
		return False


def sort_fastq(read1):
	count = 0
	ori_reads = []
	with open(read1, "r") as r1:
		while True:
			block = list(its.islice(r1, 4))
			count += 1
			if not block:
				break
			block[2] = block[2].strip() + re.sub('@', '', block[0])
			ori_reads.append(block)
	sorted_ori_reads = sorted(ori_reads, key=lambda x: x[0].split(" ")[0])
	with open("sorted_" + read1, "w+") as wf1:
		for block in sorted_ori_reads:
			for line in block:
				wf1.write(line)


def to_table(simfile):
	names = []
	contents = {}
	count = 0

	'''with open(simfile,"r") as sf:
        for line in sf:
            if count == 0:
                names = line.strip().split(" ")
                for n in names:

                    if len(n) > 0 and n[0] == ">":

                        contents.update({n[1:]:[]})

                names = list(contents.keys())
                print(names)
            if count >= 2:

                for n in line.strip().split(" ")[1:]:
                    if len(n)>0:
                        contents[names[count-2]].append(n)
            count += 1'''
	with open(simfile, "r") as simf:
		table = csv.reader(simf)
		count = 0
		for row in table:
			if count == 0:
				names = row[1:]
				for n in names:
					contents.update({n.replace("_", "\_"): []})

				names = list(contents.keys())
				print(names)
			else:
				for n in row[1:]:
					if len(n) > 0:
						contents[names[count - 1]].append(n)
			count += 1
	header = ""
	count = 0
	for n in names:
		if count == 0:
			header += "\\textbf{Virus name}     "
		else:
			if count >= len(names) - 3:
				header += "\\textbf{" + n + "}          "
		count += 1
	header += "\n           \\hline\n"
	print(header)
	count = 0
	for k in list(contents.keys()):
		print(k + " & ", end="")
		temp_c = 0
		length = len(contents[k])
		#print(contents[k])
		for n in contents[k]:
			if length-3 <= temp_c < length - 1:
				print(" " + n + " &", end="")
			elif temp_c == length -1:
				print(" " + n + " \\\\")
			temp_c += 1
		print("\n\\hline")


def search_seq(seq_file, ref1, ref2):
	seqs = []
	seq_names = []
	count = 0
	tmp_seq = ""
	ref1_seq = ""
	ref2_seq = ""
	with open(ref1, "r") as rf1:
		for line in rf1:
			if line[0] == ">":
				ref1_name = line.strip()
			else:
				ref1_seq += line.strip()
	with open(ref2, "r") as rf2:
		for line in rf2:
			if line[0] == ">":
				ref2_name = line.strip()
			else:
				ref2_seq += line.strip()
	with open(seq_file, "r") as f3:
		for line in f3:
			if ">" not in line:
				tmp_seq += line.strip()
			else:
				if count > 0:
					seqs.append(tmp_seq)
				tmp_seq = ""
				seq_names.append(line.strip())
				count += 1
		seqs.append(tmp_seq.strip())

	for i, target_seq in enumerate(seqs):
		dist1 = editdistance.eval(target_seq, ref1_seq)
		dist2 = editdistance.eval(target_seq, ref2_seq)
		if dist1 > dist2:
			print(seq_names[i], len(target_seq), dist1, dist2)


def get_changed_pos(ref1, ref2, sam_list):
	input1 = ""
	input2 = ""
	with open(ref1, "r") as rf1:
		for line in rf1:
			if line[0] != ">":
				input1 += line.strip()
			else:
				print("ref is ", line[1:].strip())
	with open(ref2, "r") as rf2:
		for line in rf2:
			if line[0] != ">":
				input2 += line.strip()
			else:
				print("change is ", line[1:].strip())
	pos = {}
	pos_match = {}
	for i in range(0, len(input1)):
		if input1[i] != input2[i]:
			if i >= 150 and i - 150 <= 3821:
				pos.update({i - 150: (input1[i], input2[i])})
				pos_match.update({i - 150: 0})

	for sam in sam_list:

		srr_num = re.search('[0-9]+', sam)
		# if int(srr_num.group(0)) == 59 or srr_num == 60:
		#    extended = True
		# else:
		#    extended = False
		extended = True
		print(srr_num.group(0))

		with open(sam, "r") as f:
			for line in f:
				fields = line.strip().split(" ")
				if extended:
					index = int(fields[2]) - 150 - 1
				else:
					index = int(fields[2]) - 1
				for k, v in pos.items():
					if index <= k <= index + len(fields[3]) - 1:
						read_base = fields[3][k - index]
						if read_base == v[1]:
							pos_match[k] += 1

		for k, v in pos_match.items():
			if v > 0:
				print(k, ":", v)

	with open("changed_nt_bases_spike.txt", "w+") as wf1:
		for k, v in pos.items():
			if pos_match[k] > 0:
				wf1.write(str(k + 151) + ":" + v[0] + "|" + v[1] + " " + str(pos_match[k]) + ", ")


def kmer_repeat(samfile, ref_file, readfile1, readfile2):
	ref = ""
	with open(ref_file, "r") as rf:
		ref = "".join([line for line in rf.readlines() if line[0] != ">"])
	# print(ref)

	all_read = bm.read_sam(samfile)
	read_freq = Counter([x[3] for x in all_read])
	freq2_read = [x for x in all_read if read_freq[x[3]] > 1]

	loop_read = []
	added_read_set = set({})
	for k in range(5, 149):
		curr = 0
		for read in freq2_read:
			if len(read[3]) < k:
				ktmp = len(read[3]) - 2
			else:
				ktmp = k
			if read[3][:ktmp] == read[3][-ktmp:]:
				if read[3] not in added_read_set:
					read.append(k)
					loop_read.append(read)
					added_read_set.add(read[3])

			# print("k is ",k,read)

	write_sam(loop_read, "")
	file_list = []
	for read in loop_read:
		new_ref = ref[:read[2] + len(read[3]) - 1] + read[3][read[-1]:] + ref[read[2] + len(read[3]):]
		# file_list.append(str(read[2])+"_loop_ref.fasta")
		# with open(str(read[2])+"_loop_ref.fasta","w+") as wf:
		#    wf.write(">"+str(read[2])+"_loop_read_reference\n")
		#    wf.write(new_ref)
		file_list.append(str(read[2]) + "_combined_str.fa")
		combined_str = read[3] + read[3][read[-1]:]
		with open(str(read[2]) + "_combined_str.fa", "w+") as wf:
			wf.write(">" + str(read[2]) + "_r_to_r[k:]" + "\n")
			wf.write(combined_str)

	for new_ref_file in file_list:
		verify_sub_command = os.path.dirname(
			__file__) + "/find_sub.sh" + " " + "-r" + " " + new_ref_file + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"

		verify_proc = subprocess.run(verify_sub_command,
									 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
									 shell=True)
		print(verify_proc.stdout.split("\n")[-5:])
		gap_list = verify_proc.stdout.split("\n")[-3]
		if gap_list != "no gaps":
			with open(new_ref_file + "_gaps", "w+") as gf:
				gf.write(gap_list)


def write_sam(readlist, filename, mode="w+"):
	if filename != '':
		with open(filename, mode) as wf:
			for line in readlist:
				s = ""
				for field in line:
					if field == line[-1]:
						if field == "False":
							s += "\n"
						else:
							s += str(field) + "\n"
					else:
						s += str(field) + " "
				wf.write(s)
			# wf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] +  "\n")
	else:
		for line in readlist:
			for field in line:
				if field == line[-1]:
					print(field)
				else:
					print(field, end=" ")

	# print(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")


def remove_cov(ref_file, samfile):
	readlist = bm.read_sam(samfile)
	matrix = bm.matrix_from_readlist(readlist, 0, set({}))

	ref = ""
	with open(ref_file, 'r') as refg:
		for line in refg:
			if ">" not in line:
				ref += line.strip()

	real_narrowed, paired_real_narrowed, nearly_real_narrowed, potential_mutated = bm.narrow_reads(ref,
																								   matrix.narrowed_read,
																								   "", True)

	rigor_matrix = bm.matrix_from_readlist(paired_real_narrowed, 0.975, set({}), False,
										   matrix, "real_narrowed")

	min_cvg = 1000000000
	min_index = -1
	all_cvg = {}
	for i in range(0, rigor_matrix.real_narrowed_matrix.shape[1]):
		tmp = np.squeeze(rigor_matrix.real_narrowed_matrix.getcol(i).toarray())
		tmp_count = np.bincount(tmp)[1:]
		all_cvg.update({i: sum(tmp_count)})
		if sum(tmp_count) < min_cvg and 1000 < i < 28000:
			min_cvg = sum(tmp_count)
			min_index = i
	all_cvg_sorted = [(k, v) for k, v in sorted(all_cvg.items(), key=lambda x: x[1]) if 40 < k < 29800]
	print("minumum cvg", min_cvg, "at", min_index, all_cvg_sorted[:10])

	tmp = np.squeeze(rigor_matrix.real_narrowed_matrix.getcol(min_index).toarray())
	read_nums = np.nonzero(tmp)
	print(read_nums)
	delete_list = [readlist[x] for x in read_nums[0]]
	delete_set = [readlist[x][3] for x in read_nums[0]]
	# print(delete_list)
	new_read_list = [x for x in readlist if x[3] not in delete_set]
	print(len(new_read_list), "reads left")
	new_matrix = bm.matrix_from_readlist(new_read_list, 1, set({}))
	print(new_matrix.narrowed_matrix.shape)
	gaps = []
	for i in range(new_matrix.narrowed_matrix.shape[1]):
		tmp = np.squeeze(new_matrix.narrowed_matrix.getcol(i).toarray())
		tmp_count = np.bincount(tmp)[1:]

		if i == min_index:
			print(tmp_count)
			print(np.nonzero(tmp))
		if sum(tmp_count) == 0:
			gaps.append(i)
	print("gaps at ", gaps)
	if len(gaps) > 0:
		with open("removed_gaps.txt", "w+") as wf:
			gaps_str = ",".join([str(i) for i in gaps])
			wf.write(gaps_str)


def get_batch_gap(read_file):
	readlist = []
	with open(read_file, "r") as sf:
		for line in sf:
			fields = line.strip().split(" ")
			if "N" in fields[3]:
				continue
			readlist.append(fields)
	print(len(readlist), "total reads")
	gapped_reads = []
	batches = []
	for ir, read in enumerate(readlist):

		gaps = check_sam_gap("batch_" + str(ir) + "_reference.fa_output/paired_real_narrowed_extract.sam")
		if len(gaps) > 0:
			batches.append(ir)
		else:
			gapped_reads.append(read)

	print(batches)
	print(len(gapped_reads), "single gap reads")
	write_sam(gapped_reads, "single_gap_read.sam", "w")


def single_gap_reads(ref, strain, samfile, readfile1, readfile2):
	misp_bit = 6
	# strain = 2
	if strain < 0:
		print("invalid strain number, error")
		exit(-2)
	prev_readset = set({})
	if strain > 0:
		for i in range(0, strain):
			print("getting reads in subbed_reads_" + str(i))
			with open("subbed_reads_" + str(i) + ".sam", "r") as subf:
				for line in subf:
					sline = line.strip().split(" ")
					prev_readset.add(sline[3])
			with open("rejected_subbed_reads_" + str(i) + ".sam", "r") as rejf:
				for line in rejf:
					fields = line.strip().split(" ")
					prev_readset.add(fields[3])
	print(len(prev_readset))
	readlist = []
	with open(samfile, "r") as samf:
		for line in samf:
			if line.split(" ")[3] not in prev_readset:
				fields = line.strip().split(" ")
				if "N" in fields[3]:
					continue
				readlist.append(fields)
	readlist = fix_s_pos(readlist)
	print(len(readlist), "candidate reads")
	ref_seq = ""
	with open(ref, "r") as rf:
		for line in rf:
			if line[0] != ">":
				ref_seq += line.strip()
	ref = ref_seq

	for ir, read in enumerate(readlist):
		tmp_misP = []
		read_index = int(read[2]) - 1
		for i, base in enumerate(read[3]):
			if ref[read_index + i] != base:
				tmp_misP.append(read_index + i)
		for mp in tmp_misP:
			read.append((int(mp), ref[mp], read[3][mp - read_index]))

	batch = 0

	# while len(readlist) > 0:
	subbed_read = []
	misPs = []
	ref = ref_seq
	covered_pos = {}
	for batch, read in enumerate(readlist):
		# if batch <= 1302:
		#   continue
		overlap = False
		read_index = int(read[2]) - 1
		read_misPs = [int(x[0]) for x in read[misp_bit:]]

		print("curr read", batch, read, "bases covered", len(covered_pos))
		if len(subbed_read) > 0:
			# read clash detection, find all overlapping reads

			# mp-clash, compare all MPs, if the reads only itroduce new MP outise current mp ranges

			clash = False
			for rei, base in enumerate(read[3]):
				r_start = int(read[2]) - 1
				if r_start + rei in covered_pos.keys():
					if read[3][rei] != covered_pos[r_start + rei]:
						clash = True
						print("clash at", r_start + rei, covered_pos[r_start + rei], read[3][rei])
						break

			if clash:
				batch += 1
				continue

		temp_ref = ref[:read_index] + read[3] + ref[read_index + len(read[3]):]
		# print(editdistance.eval(temp_ref,ref))
		with open("batch_" + str(batch) + "_reference.fa", "w+") as bf:
			bf.write(">batch_" + str(batch) + "\n")
			bf.write(temp_ref)
		verify_sub_command = os.path.dirname(
			__file__) + "/find_sub.sh" + " " + "-r" + " " + "batch_" + str(
			batch) + "_reference.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y"

		verify_proc = subprocess.run(verify_sub_command,
									 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
									 shell=True)
		print(verify_proc.stdout.split("\n")[-5:])
		gap_list = verify_proc.stdout.split("\n")[-3]
		if gap_list != "no gaps":
			# print(read)
			# print(gap_list)
			gap_list = re.sub('[\[\]]', '', gap_list)
			print(gap_list.split(",")[:-1])

		# gap_list = gap_list.split(",")[:-1]
		# ref = ref[:read_index] + ref_seq[read_index:read_index+len(read[3])] + ref[read_index + len(read[3]):]
		# subbed_read.remove(read)
		else:

			print(batch, "no gaps for curr read\n")
			subbed_read.append(read)
		# print(len(covered_pos),covered_pos.keys())
		# print(subbed_read,len(readlist))
		batch += 1

	with open("2nd_single_gap_read.sam", "w+") as bf:
		for line in subbed_read:
			bf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + " " + str(
				line[5]) + "\n")

	exit()


def freq_info(sub_sam, valid_sam):
	all_reads = []
	read_freq = {}
	with open(sub_sam, "r") as s1f:
		for line in s1f:
			fields = line.strip().split(" ")
			read_freq.update({fields[3]: int(fields[5])})
			all_reads.append(fields)
	val_reads = []
	with open(valid_sam, "r") as s2f:
		for line in s2f:
			fields = line.strip().split(" ")
			fields.append(read_freq[fields[3]])
			val_reads.append(fields)
	write_sam(val_reads, "freq_" + valid_sam, "w+")

def sequence_format(ref_file):
	ref = ""
	with open(ref_file,"r") as rf:
		for line in rf:
			if line[0] != ">":
				ref += line.strip()

	print(len(ref))
	i = 1
	space = 5
	#print("     ",1,end=" ")
	for b in ref:
		print(b,end="")
		if i % 20 == 0:
			if i % 100 == 0:
				space_str = (space- len(str(i)))*" "
				print("",i,"\\\\\n",space_str,i+1,end=" ")
			else:
				print(" ",end="")
		i += 1

def misp_syno(ref_file,sub_file,misp_file,code_file):


	endpoints = '''266..21555
	21563..25384
	25393..26220
	26245..26472
	26523..27191
	27202..27387
	27394..27759
	27756..27887
	27894..28259
	28274..29533
	29558..29674'''
	endpoints_list = endpoints.split("\n")
	protein_loc = []#(265,21555)
	for se in endpoints_list:
		start = se.split("..")[0]
		end = se.split("..")[1]
		protein_loc.append(range(int(start)-1,int(end)-1))
	translate = {}
	ref = ""
	spike_range = range(21562,25383)
	with open(code_file,"r") as cf:
		for line in cf:
			cols = line.strip().split(" ")
			pair = re.findall('([A-Z]{3}\ [A-Z*])',line)
			if len(pair) > 0:
				for p in pair:
					[nt, aa] = p.split(" ")
					translate.update({nt:aa})
	with open(ref_file,"r") as rf:
		for line in rf:
			if line[0] != ">":
				ref += line.strip()
	misP,misp_source,misp_reads =  get_misp(ref_file,[sub_file],False)
	syn_stat = {}
	move = []
	for pos,clist in misP.items():
		found = False
		for r in protein_loc:
			#move bases inside this protein region only

			if pos in r:
				found = True
				for change in clist:
					#location = int((pos-r.start)/3)
					index = (pos - r.start) % 3


					ref_codon = ref[pos-index:pos+(3-index)]
					if change[-1] != "-":
						mut_codon = ref_codon[0:index] + change[-1] + ref_codon[index+1:]
					else:
						move.append(pos)

						mut_codon = ref_codon[pos-index:pos] + ref[pos+1:pos+(3-index)+1]

					print(ref_codon,mut_codon,index, ref_codon[0:index],change,ref_codon[index+1:], translate[ref_codon],translate[mut_codon])
					if len(mut_codon) > 3:
						print(pos,clist)
						print(mut_codon)
						exit()
					del_move = len([x for x in move if x < pos])
					if translate[ref_codon] == translate[mut_codon]:
						print("synonymous ",pos,change)
						syn_stat.update({pos:"synonymous"})
					else:
						print("non-synonymous",pos,change)

						syn_stat.update({pos: "non-synonymous"})
		if not found:
			syn_stat.update({pos:"synonymous"})
	with open(misp_file,"r") as mif:
		for line in mif:
			m = re.search('^([0-9]+)',line)
			if m is None:
				continue
			index = int(m.group(1))
			new_line = line.replace('\\\\',' & '+syn_stat[index]+'\\\\')
			new_line = new_line.replace(str(index),str(index+1))
			print(new_line,end="")
			print("\hline")
def fix_position(sub_file,misp_file):
	insert = []
	delete = []
	read_list = []
	with open(sub_file,"r") as rf:
		for line in rf:
			fields = line.strip().split(" ")
			read_list.append(fields)
			index = int(fields[2])-1
			blocks = re.findall('[0-9]+[A-Z]',fields[4])
			start = index
			for blk in blocks:
				m = re.search('([0-9]+)([A-Z])',blk)
				num = int(m.group(1))
				type = m.group(2)
				if type == "I":
					insert.extend([ start + x for x in range(num)])
				elif type == "D":
					delete.extend([ start + x for x in range(num)])

				start += num

	print(insert,delete)
	new_misp_lines = []
	with open(misp_file,"r") as mf:
		for line in mf:

			m = re.search('^([0-9]+)',line)
			if m is None:
				continue
			index = m.group(1)

			del_pos =  len([x for x in delete if x < int(index)])
			ins_pos =  len([x for x in insert if x < int(index)])
			move_pos = ins_pos-del_pos
			new_l = re.sub(index,str(int(index)+move_pos+1),line)
			new_misp_lines.append(new_l)
	new_sub_lines = []
	with open(sub_file,"r") as mf:
		for line in mf:

			m = re.search('([0-9]+) [A-Z]',line)

			if m is None:
				continue
			index = m.group(1)

			del_pos =  len([x for x in delete if x < int(index)])
			ins_pos =  len([x for x in insert if x < int(index)])
			move_pos = ins_pos-del_pos

			new_l = re.sub(index,str(int(index)+move_pos),line)

			new_sub_lines.append(new_l)
			#print(new_l)
	with open("fixed_"+misp_file,"w+") as wf:
		for line in new_misp_lines:
			wf.write(line)
			wf.write("\hline\n")
	with open("fixed_"+sub_file,"w+") as w1f:
		for line in new_sub_lines:
			w1f.write(line)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="some tool functions", usage='%(prog)s [options]' + str(sys.argv))
	parser.add_argument('--rev_comp', type=str, default="no")
	parser.add_argument("--remove_by_ID", type=str, default="no")
	parser.add_argument("files", type=str, default="", nargs="+")
	parser.add_argument("--split_seq", type=str, default="no")
	parser.add_argument("--pair_dist", type=str, default="no")
	parser.add_argument("--overlap_id", type=str, default="no")
	parser.add_argument("--get_read_by_id", type=str, default="no")
	parser.add_argument("--test_shell", type=str, default="no")
	parser.add_argument("--get_gap_reads", type=str, default="no")
	parser.add_argument("--get_ori_half", type=str, default="no")
	parser.add_argument("--count_fq_freq", type=str, default="no")
	parser.add_argument("--strain_num", type=int, default=-1)
	parser.add_argument("--remove_freq", type=int, default=-1)
	parser.add_argument("--draw_coverage", type=str, default="no")
	parser.add_argument("--count_subbed_freq", type=str, default="no")
	parser.add_argument("--filter_sam_freq", type=str, default="no")
	parser.add_argument("--sam_gap", type=str, default="no")
	parser.add_argument("--get_misps", type=str, default="no")
	parser.add_argument("--sort_fastq", type=str, default="no")
	parser.add_argument("--oneline_fasta", type=str, default="no")
	parser.add_argument("--to_table", type=str, default="no")
	parser.add_argument("--search_seq", type=str, default="no")
	parser.add_argument("--get_changed_pos", type=str, default="no")
	parser.add_argument("--verify_misp", type=str, default="no")
	parser.add_argument("--get_coverage", type=str, default="no")
	parser.add_argument("--diff_coverage", type=str, default="no")
	parser.add_argument("--kmer_repeat", type=str, default="no")
	parser.add_argument("--remove_cov", type=str, default="no")
	parser.add_argument("--single_gap_reads", type=str, default="no")
	parser.add_argument("--get_batch_gap", type=str, default="no")
	parser.add_argument("--freq_info", type=str, default="no")
	parser.add_argument("--seq_format",type=str,default="no")
	parser.add_argument("--misp_syno",type=str, default="no")
	parser.add_argument("--fix_position", type=str, default="no")
	args = parser.parse_args()

	input_files = args.files
	print(input_files[1:])
	# hash_str()
	# count_match()
	# countdiff()
	# count_record()
	# oneline_fasta()
	if args.overlap_id != "no":
		getoverlap(input_files[1], input_files[2])
	# get_sam_dup()

	# adjustindex(377)
	# sep_reads()

	# seq_depth()

	# extract_to_read()
	# seq_depth()
	# hash_str()
	# count_length()
	# missing_narrow()

	# check_order()
	# fix_ID()

	# get_contigs()
	if args.rev_comp != "no":
		get_rev_comp(input_files[1], True)
	if args.get_read_by_id != "no":
		get_read_from_listf(input_files[1], input_files[2], input_files[3])
	# rev_comp_read()
	if args.remove_by_ID != "no":
		remove_ID_fastq()
	if args.split_seq != "no":
		split_seq(input_files[1], length=2949, step=300)
	if args.pair_dist != "no":
		count_pair_dist(input_files[1])
	if args.test_shell != "no":
		test_shell(input_files[1], input_files[2])
	if args.get_gap_reads != "no":
		get_gap_reads(input_files[1], args.strain_num, input_files[2], input_files[3], input_files[4])
	if args.get_ori_half != "no":
		get_ori_half(input_files[1], input_files[2], input_files[3])
	if args.count_fq_freq != "no":
		count_fq_freq(input_files[1])
	if args.get_coverage != "no":
		get_coverage(input_files[1], input_files[2])
	if args.draw_coverage != "no":
		draw_coverage(input_files[1:])
	if args.diff_coverage != "no":
		diff_coverage(input_files[1], input_files[2])
	if args.count_subbed_freq != "no":
		count_subbed_freq(input_files[1], input_files[2])
	if args.filter_sam_freq != "no":
		filter_sam_freq(input_files[1])
	if args.remove_freq > 0:
		remove_freq(input_files[1], args.remove_freq)
	if args.sam_gap != "no":
		check_sam_gap(input_files[1])
	if args.get_misps != "no":
		get_misp(input_files[1], input_files[2:])
	if args.sort_fastq != "no":
		sort_fastq(input_files[1])
	if args.oneline_fasta != "no":
		oneline_fasta(input_files[1])
	if args.to_table != "no":
		to_table(input_files[1])
	if args.search_seq != "no":
		search_seq(input_files[1], input_files[2], input_files[3])
	if args.get_changed_pos != "no":
		get_changed_pos(input_files[1], input_files[2], input_files[3:])
	if args.verify_misp != "no":
		fac_verify_misp(input_files[1], input_files[2], input_files[3], input_files[4], input_files[5])
	if args.kmer_repeat != "no":
		kmer_repeat(input_files[1], input_files[2], input_files[3], input_files[4])
	if args.remove_cov != "no":
		remove_cov(input_files[1], input_files[2])
	if args.single_gap_reads != "no":
		single_gap_reads(input_files[1], args.strain_num, input_files[2], input_files[3], input_files[4])
	if args.get_batch_gap != "no":
		get_batch_gap(input_files[1])
	if args.seq_format != "no":
		sequence_format(input_files[1])
	if args.freq_info != "no":
		freq_info(input_files[1], input_files[2])
	if args.misp_syno != "no":
		misp_syno(input_files[1],input_files[2],input_files[3], input_files[4])
	if args.fix_position != "no":
		fix_position(input_files[1], input_files[2])
	corenum = mp.cpu_count() - 2
	line_amount = 245216120

	manager = mp.Manager()
	lock = manager.Lock()
	dup_count = {}
	# print(dup_count)
	'''
    with mp.Pool(corenum) as pool:
        start = 0
        lparam = []
        while start < line_amount:
            lparam.append((start,start+int(line_amount/corenum),lock,dup_count))
            start += int(line_amount/corenum)+1

        results = pool.starmap(hash_str,lparam )
    '''
# hash_str(0,line_amount,lock,dup_count)
# print(dup_count)

# with open("dup_reads_1","w+") as f:
#    for i in dup_count.keys():
#        #if(dup_count.get(i)[1]>1):
#        f.write(i+str(dup_count.get(i))+" \n")
# count_match()
