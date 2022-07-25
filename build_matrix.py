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
from dataclasses import dataclass, field


@dataclass
class Matrix_info:
	max_shape: tuple = field()
	row: np.ndarray
	col: np.ndarray
	val: np.ndarray
	narrowed_read: list
	insertion_reads: dict
	marked_id: set
	# read_list : list
	nearly_real_narrowed_read : list = field(default="")
	nearly_real_narrowed_matrix : sp.coo_matrix = field(default=None)
	real_narrowed_read: list = field(default="")
	maxposition: int = field(default=-1)
	insertion_columns_list: list = field(default="")
	possible_insert: list = field(default="")
	narrowed_matrix: sp.coo_matrix = field(default="")
	real_narrowed_matrix: sp.coo_matrix = field(default="")


# narrowed_read = []
dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))
# maxposition = -1
# row = np.array([])
# col = np.array([])
# val = np.array([])
# cor_record = {}
# real_len = []
# gene_length =-1
# match_limit = -1
# read_list = []
markbit = 5


# marked_id = set({})
# marked_row_num = []
# all_marked_id = []
# half_real_reads = []
# half_real_ID = set({})

def read_sam(R_file):
	r = pd.read_csv(R_file, delimiter=' ', names=['ID', 'strand', 'sta_p', 'sam_q', 'cigar'], encoding='unicode_escape')
	read_number = r.shape[0]
	# real_len = [0] *read_number
	read_list = []
	for i in range(read_number):
		read_list.append([str(r["ID"].loc[i]), int(r["strand"].loc[i]), int(r["sta_p"].loc[i]), str(r["sam_q"].loc[i]),
						  str(r["cigar"].loc[i]), False])

	print("r is", r)

	return read_list


def matrix_from_readlist(all_read, match_limit, marked_id, initial=True, matrix_info=None, target="real_narrowed"):
	# global narrowed_read,maxposition, read_number,row,val,col

	read_number = 0
	insertion_reads = {}
	row_l = []
	col_l = []
	val_l = []
	narrowed = []
	maxposition = 0  # debug variable
	maxindex = 0  # debug variable
	labindexes = {}  # debug
	added_read = {}  # read index with added N in the end and the original length
	# read reference genome
	changed = 0  # debug
	included_i = 0
	if not initial:
		# building matrix for different narrowed reads
		if matrix_info is None or (target != "nearly_real_narrowed" and target != "real_narrowed" and target != "raw"):
			print("didn't provide matrix_info, or correct target exiting")
			exit(4)
		for i in range(len(all_read)):
			exclude = False
			iread = all_read[i]

			index = iread[2] - 1

			sam_q = iread[3]

			tmp_length = len(iread[3])


			# print("sam_q:", sam_q, "\n")

			# configuring sam_q_num for matrix
			sam_q_num = []

			# softclipping adjust index position
			s_start = 0
			s_end = 0

			c = 0
			inserts = []
			begin = 0
			reduce = 0  # deduct "non-existent" bases in total length, "D" and "H" not shown in reads
			for j in range(0,len(sam_q)):
				sam_q_num.append(dict_tonu[sam_q[j]])

			val_l.extend(sam_q_num)
			# row_tmp = [int(included_i) for n in range(tmp_length)]
			row_tmp = [int(included_i)] * tmp_length
			row_l.extend(row_tmp)
			col_tmp = [n for n in range(index, index + tmp_length)]
			col_l.extend(col_tmp)
			if len(sam_q_num) != len(row_tmp) or len(sam_q_num) != len(col_tmp):
				print(iread[0], i, index, len(row_tmp), len(col_tmp), len(sam_q_num), tmp_length, iread[4])
				print(sam_q_num)
				exit(14)
			if (index + len(sam_q)) > maxposition:
				maxposition = index + tmp_length
				maxindex = index
			included_i += 1

		# get possible max inserted columns
		# print(insertion_reads)
		target_matrix =	 sp.coo_matrix((val_l, (row_l,col_l))).tocsc()  # matrix
		if target == "real_narrowed":
			matrix_info.real_narrowed_matrix = target_matrix.copy()
		elif target == "nearly_real_narrowed":
			matrix_info.nearly_real_narrowed_matrix = target_matrix.copy()
		info_collection = matrix_info

	else:
		# initializing by selecting reads with match_limit
		reads = []
		for i0 in range(len(all_read)):
			reads.append(all_read[i0][3])
		freq_reads = Counter(reads)

		for i in range(len(all_read)):
			exclude = False
			iread = all_read[i]

			index = iread[2] - 1

			sam_q = iread[3]
			cigar = iread[4]
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
				if target == "raw" and bt == "M" or bt=="I" or bt=="D":
					matched += int(m.group(1))
				else:
					if bt == "M":  # or bt=="I" or bt=="D":
						matched += int(m.group(1))
				blk_type.append(bt)

				blk_pos.append(bl)
				ini = bl
				# if iread[4]=="140M2D10M":
				#	print(block, int(m.group(1)))
				base_length += int(m.group(1))
				if bt != "I" and bt != "H":
					tmp_length += int(m.group(1))  # get read length without counting clipped bases
			if target != "raw":
				#if (matched / base_length < match_limit) or iread[0] in marked_id:
				if (matched / base_length < match_limit) or iread[0] in marked_id:
					# if (len(cigar_str)>1 or not re.match('^[0-9]+[M]$',cigar_str[0])):
					# if(matched/read_length<match_limit):

					# print(i,str(r.loc[i]), matched)
					if  not (freq_reads[sam_q] > 1 and iread[0] not in marked_id):
						exclude = True
						continue

			narrowed.append(iread)
			# print("sam_q:", sam_q, "\n")

			# configuring sam_q_num for matrix
			sam_q_num = []

			# softclipping adjust index position
			s_start = 0
			s_end = 0

			c = 0
			inserts = []
			begin = 0
			reduce = 0  # deduct "non-existent" bases in total length, "D" and "H" not shown in reads
			curr_pos = 0
			for j in range(0, blk_pos[-1]):  # change here to fill the blank with 0?

				if blk_type[c] == "M":
					try:
						sam_q_num.append(dict_tonu[sam_q[j - reduce]])
					except:
						print(j - reduce, len(sam_q), i)
						exit(1)


				elif blk_type[c] == "I":
					inserts.append(dict_tonu[sam_q[j - reduce]])
				elif blk_type[c] == "S":

					sam_q_num.append(dict_tonu[sam_q[j - reduce]])
				elif blk_type[c] == "H":
					sam_q_num.append(0)
				elif blk_type[c] == "D":
					inserted_read =  iread[3][:curr_pos]+"-"+iread[3][curr_pos:]
					iread[3] = inserted_read
					sam_q_num.append(6)
				if blk_type[c] == "H" or blk_type[c] == "D":
					reduce += 1

				if j == blk_pos[c] - 1:  # update start and c, put inserts into hashtable

					if blk_type[c] == "I":

						if i in insertion_reads.keys():

							newinsert = copy.deepcopy(insertion_reads.get(included_i))
							newinsert.append((index + begin, inserts))
							insertion_reads.update({included_i: newinsert})
						else:
							insertion_reads.update({included_i: [(index + begin, copy.deepcopy(inserts))]})
					begin = blk_pos[c]
					inserts = []

					c += 1
					if c == len(blk_type):
						break
				curr_pos += 1
			if blk_type[0] == "S":
				if index - blk_pos[0] < 0:

					start_pos = blk_pos[0] - index
					tmp_length -= start_pos
					sam_q_num = sam_q_num[start_pos:]
					if i == 3:
						print("tmp_length", tmp_length, "base_length", base_length)
				else:
					index = index - blk_pos[0]
			'''
			if len(sam_q_num) < read_length:
				added_read.update({i: len(sam_q_num)})
				sam_q_num += [0] * (read_length - len(sam_q_num))
				pad = 0
			else:
				pad = len(sam_q_num) - read_length
			'''

			val_l.extend(sam_q_num)
			# row_tmp = [int(included_i) for n in range(tmp_length)]
			row_tmp = [int(included_i)] * tmp_length
			row_l.extend(row_tmp)
			col_tmp = [n for n in range(index, index + tmp_length)]
			col_l.extend(col_tmp)
			if len(sam_q_num) != len(row_tmp) or len(sam_q_num) != len(col_tmp):
				print(iread[0], i, index, len(row_tmp), len(col_tmp), len(sam_q_num), tmp_length, iread[4])
				print(sam_q_num)
				exit(14)
			if (index + len(sam_q)) > maxposition:
				maxposition = index + tmp_length
				maxindex = index
			included_i += 1

		# get possible max inserted columns
		# print(insertion_reads)
		insertion_lengths = {}
		for readnum, possible_insert in insertion_reads.items():
			for insertion_tuple in possible_insert:
				starting_index = insertion_tuple[0]
				if starting_index in insertion_lengths.keys():
					if len(insertion_tuple[1]) > insertion_lengths[starting_index]:
						insertion_lengths[starting_index] = len(insertion_tuple[1])
				else:
					insertion_lengths[starting_index] = len(insertion_tuple[1])
		extra_col_possible = sum(insertion_lengths.values())

		info_collection = Matrix_info(max_shape=(included_i, maxposition + extra_col_possible),
									  row=np.array(row_l), col=np.array(col_l), val=np.array(val_l),
									  narrowed_read=narrowed, insertion_reads=insertion_reads,
									  marked_id=marked_id)
		# print(np.bincount(info_collection.val))
		csc = sp.coo_matrix((info_collection.val, (info_collection.row, info_collection.col))).tocsc()  # matrix
		print("max position at", info_collection.maxposition, info_collection.col[-1], maxindex)
		print("insertion_reads", len(insertion_reads))
		print("csc shape", csc.shape)
		info_collection.narrowed_matrix = csc.copy()

	#	if initial:
	#		narrowed_read = narrowed.copy()
	return info_collection


def build_insertion(intermit_matrix_info, count_threshold):
	#print("currently not used")
	return intermit_matrix_info
	# global maxposition,exclude_reads,read_number,cor_record
	max_shape = intermit_matrix_info.max_shape
	insertion_reads = intermit_matrix_info.insertion_reads
	csc = intermit_matrix_info.narrowed_matrix
	insertion_columns = set({})
	print("add_matrix", max_shape)
	add_matrix = sp.coo_matrix(max_shape, dtype=np.int32).tocsc()

	for i in insertion_reads.keys():
		for j in insertion_reads[i]:
			index = 0
			index1 = 0
			while index < len(j[1]):
				add_matrix[i, j[0] + index] = j[1][index]
				insertion_columns.add(j[0] + index)

				index += 1

	remove_columns = []
	insertion_columns_list = list(insertion_columns)
	insertion_columns_list.sort()
	# print(insertion_columns_list)
	# print(sorted([x for x in cor_record.keys()]))

	for i in insertion_columns_list:
		tmp = np.squeeze(add_matrix.getcol(i).toarray())
		tmp_count = sum(np.bincount(tmp)[1:])

		if tmp_count < count_threshold:
			# print(i,"with", tmp_count,end=", ")

			remove_columns.append(i)

	# remove corresponding insertions in correct records
	# for i in remove_columns:
	#    cor_record.pop(i)

	remove_columns.sort()
	print("remove columns ", remove_columns)
	print(len(insertion_columns_list), len(remove_columns))

	remove_copy = np.array(remove_columns)

	# remove columns with reads less than count_threshold in insertion_columns_list
	remove_set = set(remove_columns)
	# print(remove_copy)

	ir = 0
	n_icl = []
	icl_ori_pos = []
	for ir in insertion_columns_list:
		if ir not in remove_set:

			move = np.where(remove_copy < ir)
			n_icl.append(ir - move[0].shape[0])
			icl_ori_pos.append(ir)
		# print(ir, ir-move[0].shape[0], ir in remove_set, remove_copy[move])
		else:
			continue
	# print(insertion_columns_list,n_icl)
	insertion_columns_list = n_icl
	# insertion_columns_list = n_icl
	# exit(2)

	# print(cor_record)
	# print(insertion_columns_list)
	print("\nset up insertion column list ", len(insertion_columns_list), insertion_columns_list)
	prev_time = time.time()

	# ---------------move inserted columns in insertion matrix

	update_pairs = []
	# print(remove_copy)
	# del i

	for i2, i in enumerate(icl_ori_pos):
		move = np.where(remove_copy < i)
		# print(move)
		add_matrix[:, insertion_columns_list[i2]] = add_matrix[:, i]
		# move_len = move[0].shape[0]
		# tmp = np.squeeze(add_matrix.getcol(i).toarray())
		# tmp_1 = np.squeeze(add_matrix.getcol(insertion_columns_list[i2]).toarray())
		# tmp_count = np.bincount(tmp)[1:]
		# tmp_1_count = np.bincount(tmp_1)[1:]
		# update_pairs.append((insertion_columns_list[i2], cor_record[i]))

	# print(i, move_len,insertion_columns_list[i2], np.array_equal(tmp_count,tmp_1_count),tmp_count, tmp_1_count)

	# for up in update_pairs:
	#	cor_record.update({up[0]: up[1]})

	# combine insertions and other parts of csc
	# exit(-2)

	j0 = 0
	i = 0

	# ---------------------determine all inserted columns-----------------

	if len(insertion_columns_list) > 0:
		# tmp_list = bm.get_rowcolval()
		row = intermit_matrix_info.row
		col = intermit_matrix_info.col
		val = intermit_matrix_info.val
		print("copyting csc to addmatrix if insertions are sure")
		del i
		for i in insertion_columns_list:
			# print(i)
			# test_matrix = sp.coo_matrix((val, (row, col))).tocsc()
			# print("before insert: ",np.array(test_matrix.getcol(i)))
			col[col >= i] += 1
			tmp = np.squeeze(add_matrix.getcol(i).toarray())
			tmp_base = tmp[np.nonzero(tmp)]

			row = np.append(row, np.array(np.nonzero(tmp)))
			col = np.append(col, np.array([i] * len(tmp_base)))

			val = np.append(val, np.array(tmp_base))

		# tmp_val = np.array(tmp_base)
		# tmp_row = np.array(np.nonzero(tmp)[0])
		# tmp_col = np.array([i] * len(tmp_base))
		# print(tmp_val.shape,tmp_row.shape,tmp_col.shape)
		# test_matrix = sp.coo_matrix((val, (row, col))).tocsc()
		# print(np.array(test_matrix.getcol(i)))
		# print(np.array(test_matrix.getcol(i+1)))
		# exit(1)
		com_matrix = sp.coo_matrix((val, (row, col))).tocsc()
		test_icl = np.array(insertion_columns_list)
		# ------------make sure copy process is correct---------

		for i in range(0, com_matrix.shape[1]):
			# print(i)
			if i not in insertion_columns_list:
				tmp = com_matrix.getcol(i).toarray()

				tmp1 = csc.getcol(i - len(test_icl[test_icl <= i])).toarray()
				if np.sum(tmp != tmp1) != 0:
					print(i, i - len(test_icl[test_icl <= i]))
					print("com", tmp[np.nonzero(tmp)])
					# print(np.nonzero(tmp))
					print("csc", tmp1[np.nonzero(tmp1)])
					exit(5)
		# print(np.nonzero(tmp1))

		add_matrix = com_matrix

		print("copied csc to add_matrix", add_matrix.shape)
	else:
		print("no insertion")
		add_matrix = csc
	intermit_matrix_info.real_narrowed_matrix = add_matrix
	intermit_matrix_info.insertion_columns_list = insertion_columns_list

	return intermit_matrix_info


def narrow_reads(ref, narrowed_read, out_dir, brute_force=True):
	# global narrowed_read,  half_real_reads, half_real_ID

	ID_count = {}
	'''
	for rl in narrowed_read:
		if rl[0] in ID_count.keys():
			ID_count[rl[0]] += 1

		else:
			ID_count[rl[0]] = 1
	errors = []

	for id in ID_count.keys():
		if ID_count[id] > 2:
			errors.append(id)
	if len(errors) > 0:
		print(str(len(errors))+" errors in "+str(errors))
		exit(-4)
	ID_count.clear()
	#'''
	loc_pair_narrowed = {}
	loc_pair_real_narrowed = {}
	loc_pair_paired_real_narrowed = {}
	print(len(narrowed_read), " 100% M reads in narrowed_extract.sam")
	count = 0
	half_real_ID = set({})
	# half_real_reads = []
	nearly_true_total = []
	true_total_match = []


	read_ferq = Counter([x[3] for x in narrowed_read])

	for ri,rl in enumerate(narrowed_read):
		# print(mate_rl)
		index = rl[2] - 1
		# print(len(ref[index:index + len(mate_rl[3])]),len(mate_rl[3]),end=" | ")

		if ref[index:index + len(rl[3])] == rl[3]:
			true_total_match.append(rl)
		# half_real_ID.add(rl[0])
		else:
			if re.match('^[0-9]+[M]$', rl[4]):

				nearly_true_total.append(copy.deepcopy(rl))
				narrowed_read[ri][4] = rl[4]+"*"
			else:
				if read_ferq[rl[3]] > 1:
					nearly_true_total.append(copy.deepcopy(rl))
					narrowed_read[ri][4] = rl[4] + "*"


	print(len(true_total_match), " truely matched reads in real_narrowed_extract.sam")
	print(len(nearly_true_total), "nearly real narrowed reads in nearly_narrowed_extract.sam")
	with open(out_dir + "real_narrowed_extract.sam", "w+") as nf1:
		for line in true_total_match:
			nf1.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
	with open(out_dir + "nearly_real_narrowed_extract.sam", "w+") as nrnf:
		for line in nearly_true_total:
			nrnf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
	with open(out_dir + "narrowed_extract.sam", "w+") as nf1:
		for line in narrowed_read:
			nf1.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")


	narrowed_read = true_total_match
	'''print(len(half_real_ID)*2," in paired_half_real_narrowed_extract.sam")
	with open(out_dir+"paired_half_real_narrowed_extract.sam", "w+") as hnf:
		for line in all_read:
			if line[0] in half_real_ID:
				half_real_reads.append(line)
				hnf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
	with open(out_dir+"paried_half_real_narrowed_read.fa", "w+") as nrf1:
		for line in half_real_reads:
			nrf1.write(">" + line[0] + "\n")
			nrf1.write(line[3] + "\n")'''

	mated = []

	for rl in narrowed_read:
		if rl[0] in ID_count.keys():
			ID_count[rl[0]] += 1

		else:
			ID_count[rl[0]] = 1

	for mate_rl in narrowed_read:
		if ID_count[mate_rl[0]] == 2:
			mated.append(mate_rl)

	narrowed_read = mated
	print(len(mated), len(narrowed_read), " reads in paired_real_narrowed_extract.sam", len(ID_count))
	with open(out_dir + "paired_real_narrowed_extract.sam", "w+") as nf1:
		for line in narrowed_read:
			nf1.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")

	if len(mated) == 0:
		print("no read satisfies condition, exiting")
		exit(-2)


	if brute_force:
		real_narrowed_ids = set([x[0] for x in true_total_match])
		count_nr_narrowed_ids = Counter([x[0] for x in nearly_true_total])
		potential_mutated_reads = []
		with open(out_dir + "potential_mutated_extract.sam", "w+") as pmf:
			for line in nearly_true_total:
				#if (line[0] in real_narrowed_ids or count_nr_narrowed_ids[line[0] == 2]) and 13 < line[2] < 29883 :
				if (line[0] in real_narrowed_ids or count_nr_narrowed_ids[line[0] == 2] or read_ferq[line[3]]>1) and 13 < line[2] < 29883:
					potential_mutated_reads.append(line)
					pmf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")




		return true_total_match,mated,nearly_true_total,potential_mutated_reads

	return (true_total_match, mated, nearly_true_total)


def marking(read_num, cvg, narrowed_read, col=0):
	# global marked_id,narrowed_read, marked_row_num

	marked_id_list = []
	if len(read_num) == 0:
		return marked_id_list, narrowed_read

	marked = 0
	for i in reversed(read_num):

		if narrowed_read[i][markbit]:
			marked += 1
		if marked >= cvg:
			return marked_id_list, narrowed_read
	if col != 0:
		print(marked, col, read_num)
	# marked_id = set({})
	for j in reversed(read_num):

		if not narrowed_read[j][markbit]:
			narrowed_read[j][markbit] = True

			# marked_id.add(narrowed_read[j][0])
			# marked_row_num.append(j)
			marked_id_list.append(narrowed_read[j][0])

			if marked == cvg:
				return marked_id_list, narrowed_read

			marked += 1
	return marked_id_list, narrowed_read


def marking_byid(read_num, cvg, narrowed_read, marked_id, col=0):
	# global marked_id,narrowed_read, marked_row_num

	if len(read_num) == 0:
		return marked_id

	marked = 0
	for i in reversed(read_num):

		if narrowed_read[i][0] in marked_id:
			marked += 1
		if marked >= cvg:
			return marked_id, narrowed_read

	# print(marked,col, read_num)
	# marked_id = set({})
	for j in reversed(read_num):

		if narrowed_read[j][0] not in marked_id:
			# narrowed_read[j][markbit] = True

			# marked_id.add(narrowed_read[j][0])
			# marked_row_num.append(j)
			marked_id.add(narrowed_read[j][0])

			if marked == cvg:
				return marked_id, narrowed_read

			marked += 1
	return marked_id, narrowed_read

def collecting_bubbles(read_num,read_list,brute_force=False):
	collected_reads = []
	if len(read_num) == 0:
		print("no reads coverd")
		return collected_reads, read_list

	for i in read_num:
		if brute_force:
			collected_reads.append(read_list[i])
		else:
			if not read_list[i][markbit]:
				collected_reads.append(read_list[i])
				read_list[i][markbit] = True

	return collected_reads,read_list

def get_bubble_reads(r1_file, r2_file, read_list, out_dir,rc_file):
	read_set = set({})
	rc_read_set = set({})
	'''	test_flag = format(read_list[0][1],'b')[::-1]
	if test_flag[4] == "1":
		if test_flag[6] == "1":
			rc_file = r1_file
		else:
			rc_file = r2_file
	else:
		if test_flag[6] == "1":
			rc_file = r2_file
		else:
			rc_file = r1_file'''
	reg_file = r1_file if rc_file == r2_file else r2_file
	for iread in read_list:
		flag = format(iread[1], 'b')[::-1]
		if flag[4] == "1":

			rc_read_set.add(rev_comp_read(iread[3]))
		else:
			read_set.add(iread[3])
	print("rev comp file",rc_file,"regular file",reg_file)
	#print(rc_read_set)
	print("forward", len(read_set),"reverse",len(rc_read_set))
	with mp.Pool(2) as pool:
		if rc_file == r1_file:
			lparam = [(rc_read_set, rc_file, out_dir + "side_bubble_reads_R1.fastq"), (read_set, reg_file, out_dir + "side_bubble_reads_R2.fastq")]
		else:
			lparam = [(read_set, reg_file, out_dir + "side_bubble_reads_R1.fastq"),
					  (rc_read_set, rc_file, out_dir + "side_bubble_reads_R2.fastq")]
		pool.starmap(extract_read_fastq, lparam)

def extract_read_fastq(read_set,read_file, outfile):
	ori_reads = []
	with open(read_file,"r") as f1:
		for block in iter(lambda: list(its.islice(f1, 4)), []):
			tmpread = block[1].strip()
			# print(tmpid)
			if tmpread not in read_set:
				continue
			ori_reads.append(block)
	with open(outfile,"w+") as wf:
		for read in ori_reads:
			for line in read:
				wf.write(line)

def write_new_extract(narrowed_read, marked_id, out_dir, round_num):
	linecount = 0
	with open(out_dir + "paired_reads_contig_round_" + str(round_num) + "_extract.sam", "w+") as prcf:
		for line in narrowed_read:
			if line[0] in marked_id:
				linecount += 1
				prcf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
	print(linecount, "reads in paired_reads_contig_round1_extract.sam")

	linecount = 0
	with open(out_dir + "next_round_extract.sam", "w+") as f:
		for line in narrowed_read:
			if line[0] not in marked_id:
				linecount += 1
				f.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
	print(linecount, "reads left in,", len(narrowed_read), "paired_real_narrowed_extract.sam", len(marked_id) * 2,
		  "reads marked in paired_real_narrowed_extract.sam")
	'''
	ID_count = Counter(marked_id_list)
	error_id = []
	for i in ID_count.keys():
		if ID_count[i] != 2:
			#print("error in number", i, ID_count[i])
			error_id.append(ID_count[i])
	'''

	marked = 0
	with open("excluded_IDs.txt", "w+") as exf:
		for mid in marked_id:
			exf.write(mid + ",")
	print(len(marked_id))


# narrowed_read[read_num][markbit] = True

def get_ori_half(r1_file, r2_file, marked_id, out_dir, read_list):
	id_set = set({})

	for iread in read_list:
		if iread[0] not in marked_id:
			id_set.add(iread[0])
	print("getting", len(id_set) * 2, "reads from " + r1_file + " and " + r2_file)
	with mp.Pool(2) as pool:
		lparam = [(id_set, r1_file, out_dir + "half_real_R1.fastq"), (id_set, r2_file, out_dir + "half_real_R2.fastq")]
		pool.starmap(extract_fastq_read, lparam)


def extract_fastq_read(id_set, readfile, outfile):
	ori_reads = []
	with open(readfile, "r") as r1:
		for block in iter(lambda: list(its.islice(r1, 4)), []):
			tmpid = re.sub('@', '', block[0].split(" ")[0])
			# print(tmpid)
			if tmpid not in id_set:
				continue
			ori_reads.append(block)
	with open(outfile, "w+") as wr1:
		for lines in ori_reads:
			for line in lines:
				wr1.write(line)

def rev_comp_read(seq):
	seq_revc = seq.upper()[::-1]
	seq_revc = seq_revc.replace('A', 't')
	seq_revc = seq_revc.replace('T', 'a')
	seq_revc = seq_revc.replace('C', 'g')
	seq_revc = seq_revc.replace('G', 'c')
	seq_rev_comp = seq_revc.upper()
	return seq_rev_comp



# def get_half_real_reads():
#	return half_real_reads
