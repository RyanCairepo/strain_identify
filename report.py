import copy
import os, random, time, sys
import multiprocessing as mp
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_alg
import re
import argparse
import json
from typing import Counter
import itertools as its
from dataclasses import dataclass
import tabulate
import csv
import editdistance


def compare(seq1, seq2):
	'''
	if len(first_seq) < len(second_seq):
		seq1 = first_seq
		seq2 = second_seq
	else:
		seq1 = second_seq
		seq2 = first_seq'''

	start = time.time()
	m = len(seq1)
	n = len(seq2)
	row = np.array([])
	col = np.array([])
	val = np.array([])

	# row = np.append(row,[range(m+1)])
	# col = np.append(col,[0]*(m+1))
	# val = np.append(val, row)

	# col = np.append(col,[range(n+1)])
	# row = np.append(row, [0]*(n+1))
	# val = np.append(val, [range(n+1)])

	row = list(range(m + 1))
	col = [0] * (m + 1)
	val = list(range(m + 1))

	col.extend(list(range(n + 1)))
	row.extend([0] * (n + 1))
	val.extend(list(range(n + 1)))
	print(n + 1, len(list(range(n + 1))), len([-100] * (n + 1)))
	print(len(row), len(col), len(val))
	# del i, j
	i = 0
	j = 0
	prev = {}
	for i in range(1, m + 1):
		# row = np.append(row, [i]*n)
		# col = np.append(col, range(1,n+1))
		# val = np.append(val, [-100] * n)
		row.extend([i] * n)
		col.extend(list(range(1, n + 1)))
		val.extend([-100] * n)
		prev[(i, 0)] = (i - 1, 0)

	del i
	for i in range(1, n + 1):
		prev[(0, i)] = (0, i - 1)
	del i
	n_row = np.array(row)
	n_col = np.array(col)
	n_val = np.array(val)
	# print(n_row.dtype,n_col.dtype,n_val.dtype)
	D = sp.coo_matrix((n_val, (n_row, n_col))).tocsr()
	print(D.shape)
	for i in range(1, D.shape[0]):
		for j in range(1, D.shape[1]):
			if seq1[i - 1] == seq2[j - 1]:
				sub_cost = 0
			else:
				sub_cost = 1
			# print(i,j)
			if int(D[i - 1, j]) < 0 or int(D[i - 1, j - 1]) < 0 or int(D[i, j - 1]) < 0:
				print("read on uninitialized cell", i, j, D[i - 1, j], D[i - 1, j - 1], D[i, j - 1])
				exit(-1)
			cost_cand = np.array([D[i - 1, j] + 1, D[i, j - 1] + 1, D[i - 1, j - 1] + sub_cost])
			min_cost = np.min(cost_cand)
			D[i, j] = min_cost
			operation = np.argmin(cost_cand)
			if operation == 0:
				prev_cell = (i - 1, j)
			elif operation == 1:
				prev_cell = (i, j - 1)
			else:
				prev_cell = (i - 1, j - 1)
			prev[(i, j)] = prev_cell

	print("editdistance between two sequences is ", D[m, n])
	comp = time.time()
	print(comp - start, "runnning for getting editdistance")

	edtdis = editdistance.eval(seq1, seq2)
	evalruntime = time.time()
	print(evalruntime - comp, "runnning for editdistance.eval")

	if D[m, n] != edtdis:
		print("error in counting distance,", D[m, n], edtdis)
		exit()
	i = m
	j = n
	rev_path = []
	while (i, j) in prev.keys():
		if prev[(i, j)] == (i - 1, j - 1):
			if seq1[i - 1] == seq2[j - 1]:
				rev_path.append("M")
			else:
				rev_path.append("m")
		elif prev[(i, j)] == (i - 1, j):
			rev_path.append("D")
		else:
			rev_path.append("I")
		(i, j) = prev[(i, j)]
	rev_path.reverse()
	print(len(rev_path), rev_path)
	count = 1
	type = rev_path[0]
	stat_list = ""
	for i in range(1, len(rev_path)):
		# print(i, stat_list)

		if type == rev_path[i] and i != len(rev_path) - 1:
			count += 1
		else:
			if i == len(rev_path) - 1:
				if type == rev_path[i]:
					count += 1
					if type == "I":
						type = "S"
					stat_list += str(count) + type
				else:
					if rev_path[i] == "I":
						rev_path[i] = "S"
					stat_list += str(count) + type + "1" + rev_path[i]

				break

			if (i - count == 0) and type == "I":
				stat_list += str(count) + "S"
			else:
				stat_list += str(count) + type
			count = 1
			type = rev_path[i]

	print("total run time", time.time() - start)
	# print(rev_path)
	print(stat_list)


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
			if count >= len(names) - 4:
				header += "\\textbf{" + n + "}          "
		count += 1
	header += "\n           \\hline\n"
	print(header)
	count = 0
	for k in list(contents.keys()):
		print(k + " & ", end="")
		temp_c = 0
		length = len(contents[k])
		for n in contents[k]:
			if temp_c < length - 1 and temp_c > length-4:
				print(" " + n + " &", end="")
			else:
				print(" " + n + " \\\\")
			temp_c += 1
		print("\n\\hline")


# input1 = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACT"
# input2 = "GGTATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTA"
# input1 = "CATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCG"
# input2 = "GTACATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCG"
input1 = ""
input2 = ""
changed = 0  # debug

input_seqs = []
seq_name = []
count = 0
with open(sys.argv[1], 'r') as refg:
	for line in refg:
		if ">" not in line:

			input1 += line.strip()
		else:
			if count > 0:
				input_seqs.append(input1)
			input1 = ""

			seq_name.append(line.strip().strip(">"))
			count += 1

	input_seqs.append(input1.strip())
print(len(input_seqs), len(seq_name))

with open(sys.argv[2], 'r') as refg2:
	for line in refg2:
		if ">" not in line:
			input2 += line.strip()
		else:
			input2_name = line.strip()
header = copy.deepcopy(seq_name)
header.insert(0, "Virus")
table = []

# seq_name.insert(0,input2_name)
# input_seqs.insert(0,input2)


for i in range(0, len(input_seqs)):
	tmp_row = []
	for j in range(0, len(input_seqs)):
		if j != i:
			dist = editdistance.eval(input_seqs[i], input_seqs[j])
			tmp_row.append(round((len(input_seqs[i]) - dist) / len(input_seqs[i]), 3))
		# tmp_row.append(0.1*j)
		# elif j < i:
		# dist = table[j][i]
		#	print(j+1,i)

		#	tmp_row.append(table[j+1][i])
		else:
			tmp_row.append(1.00)

	# print(len(tmp_row))
	table.append(tmp_row)
# print(len(table))
for i, t in enumerate(table):
	if i < len(seq_name):
		table[i].insert(0, seq_name[i])
with open("new_similarity.csv", "w+") as wf:
	writer = csv.writer(wf)
	writer.writerow(header)
	writer.writerows(table)

# tmp_row = [input2_name]
# for name,seq in zip(seq_name, input_seqs):
#	dist = editdistance.eval(seq,input2)
#	print(len(input1),len(input2),dist)
# print(name,dist,(len(seq)-dist)/len(input2))
# print(len(seq),seq[:10],seq[len(seq)-10:])
#	tmp_row.append(round((len(input2)-dist)/len(input2),3))

# table.append(tmp_row)

# break
# print(tabulate.tabulate(table,headers='firstrow'))
#exit()
print(len(input1), len(input2))
print(editdistance.eval(input1, input2))

# input1 = input1[:1000]
# input2 = input2[:1000]


print(dist, (len(input1) - dist) / len(input2))
# compare(input1,input2)
