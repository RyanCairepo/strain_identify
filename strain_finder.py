import os, random, time, sys
import statistics

import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_alg
import re
import argparse
import json
import build_matrix as bm
import copy
# import matplotlib; matplotlib.use("Qt5Agg")
# import matplotlib.pyplot as plt
from typing import Counter


# set arguments
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
	# parser.add_argument("--first_r", type=bool, default=True)
	parser.add_argument("--r1_file", type=str)
	parser.add_argument("--r2_file", type=str, default="")
	parser.add_argument("--round", type=int, default=1)
	parser.add_argument("--excluded_IDs", type=str, default="excluded_IDs.txt")
	parser.add_argument("--find_sub", type=str, default="no")
	parser.add_argument("--bubble_mode", type=str, default="False")
	return parser.parse_args()


# convert to reverse_complement seq
def DNA_Re_complement(sequence):
	sequence = sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	comp_s = sequence.upper()
	return comp_s[::-1]


dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))

if __name__ == "__main__":
	# -------------------- main function is here-----------------------------------------------------
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

	r = pd.read_csv(R_file, delimiter=' ', names=['ID', 'strand', 'sta_p', 'sam_q', 'cigar'], encoding='unicode_escape')

	dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
	dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))

	# setting
	round_num = args.round

	gene_length = args.gene_L
	read_length = args.read_L
	read_number = r.shape[0]
	count_threshold = args.count_l
	start = 0
	domin_threshold = args.domin_l
	match_limit = args.match_l
	threshold = 0.025
	# first_round = args.first_r
	readfile1 = args.r1_file
	readfile2 = args.r2_file
	find_sub = args.find_sub
	if find_sub != "no":
		out_dir = os.path.basename(args.ref) + "_output/"
	else:
		out_dir = "round_" + str(round_num) + "_output/"
	if round_num == 1:
		first_round = True
	else:
		first_round = False
	if args.bubble_mode == "True":
		bubble_mode = True
	else:
		bubble_mode = False

	row_l = []
	col_l = []
	val_l = []

	maxposition = 0  # debug variable
	maxindex = 0  # debug variable
	labindexes = {}  # debug
	added_read = {}  # read index with added N in the end and the original length
	# read reference genome
	ref = ""
	changed = 0  # debug
	with open(args.ref, 'r') as refg:
		for line in refg:
			if ">" not in line:
				ref += line.strip()

	# -------read marked IDs from last round
	marked_id = set({})
	if round_num > 1:
		with open(args.excluded_IDs, "r") as exr:
			for line in exr:
				marked_id.update(line.split(","))
	if round_num > 1 and len(marked_id) == 0:
		print("no marked read?")
	else:
		print(len(marked_id), "reads already marked")
	read_list = bm.read_sam(R_file)

	initial_matrix_info = bm.matrix_from_readlist(read_list, match_limit, marked_id, True)
	prev_time = time.time()
	real_narrowed, paired_real_narrowed, nearly_real_narrowed = bm.narrow_reads(ref, initial_matrix_info.narrowed_read, out_dir, read_list)
	# del initial_matrix_info
	intermit_matrix_info = bm.matrix_from_readlist(real_narrowed, match_limit, marked_id,False,initial_matrix_info,"real_narrowed")
	intermit_matrix_info = bm.matrix_from_readlist(nearly_real_narrowed,match_limit,marked_id,False,intermit_matrix_info,"nearly_real_narrowed")
	#intermit_matrix_info = bm.matrix_from_readlist(paired_real_narrowed, match_limit, marked_id)
	intermit_matrix_info.real_narrowed_read = real_narrowed
	intermit_matrix_info.nearly_real_narrowed_read = nearly_real_narrowed
	#intermit_matrix_info.narrowed_read = copy.deepcopy(initial_matrix_info.narrowed_read)
	#intermit_matrix_info.narrowed_matrix = initial_matrix_info.real_narrowed_matrix.copy()

	assert intermit_matrix_info.real_narrowed_matrix.shape[1] == intermit_matrix_info.narrowed_matrix.shape[1]
	print("narrowed_reads", len(intermit_matrix_info.narrowed_read), "real narrowed",
		  len(intermit_matrix_info.real_narrowed_read), "nearly real narrowed", len(intermit_matrix_info.nearly_real_narrowed_read))
	print("narrowed_shape", intermit_matrix_info.narrowed_matrix.shape, "real narrowed shape",
		  intermit_matrix_info.real_narrowed_matrix.shape, "nearly real narrowed shape", intermit_matrix_info.nearly_real_narrowed_matrix.shape)

	#del initial_matrix_info
	# print(intermit_matrix_info.narrowed_matrix.shape[1], intermit_matrix_info.real_narrowed_matrix.shape[1])

	maxposition = intermit_matrix_info.max_shape[1]
	intermit_matrix_info = bm.build_insertion(intermit_matrix_info, count_threshold)
	#print("real_narrowed matrix formed ", intermit_matrix_info.real_narrowed_matrix.shape, "time is ", time.time() - start_time)

	add_matrix = intermit_matrix_info.real_narrowed_matrix.copy()
	insertion_columns_list = intermit_matrix_info.insertion_columns_list
	print("insertion matrix set, time", time.time() - prev_time)
	# print(sorted(np.unique(np.array(insertion_columns))))
	print("add_matrix", add_matrix.shape)
	prev_time = time.time()

	# cor_record = bm.get_cor_record()

	# find difference between contigs (without insertion) and reference
	# skip for first round
	cor_record = {}
	if not first_round:
		csc = intermit_matrix_info.real_narrowed_matrix.copy()
		newseq = []
		i = start
		for i in range(start, csc.shape[1]):

			tmp = np.squeeze(csc.getcol(i).toarray())  # matrix

			tmp_count = np.bincount(tmp)[1:]

			# col_num store the value of total_coverage_col
			col_num = sum(tmp_count)

			# maximum letter's number is max_l + 1
			# print(i, tmp)

			if (sum(tmp_count) == 0):

				if i > len(ref) - 1:
					continue
				newseq.append(ref[i - start])

				continue
			else:

				nonzero_num = np.nonzero(tmp_count)[0].shape[0]
				max_l = np.argmax(tmp_count)
				max_num = tmp_count[max_l]

				if i > len(ref) - 1:
					newseq.append(dict_tole[max_l + 1])
					cor_record.update({i: "*|" + dict_tole[max_l + 1]})
					continue

				ref_le = ref[i - start]
				ref_pos = i - start
				if dict_tole[max_l + 1] != ref[i - start]:
					# print(tmp_count)
					if sum(tmp_count) >= count_threshold:
						newseq.append(dict_tole[max_l + 1])
						if max_l == 6:
							cor_record.update({i: ref[i - start] + "|*"})
						else:
							cor_record.update({i: ref[i - start] + "|" + dict_tole[max_l + 1]})
					else:

						newseq.append(ref[i - start])
				else:
					newseq.append(dict_tole[max_l + 1])

		print("cor_record set", cor_record)

	print("phase 3 time ", time.time() - start_time, " step ", time.time() - prev_time)
	prev_time = time.time()

	# exit(2)
	p_row_l = []
	p_col_l = []
	p_val_l = []

	Seq_tru = []
	Seq_let = []
	Num_lowfre = 0
	Num_col = 0
	Num_base_porp = 0
	num_col_list = []
	# labeled reads with different bases with the dominant base,but unchanged
	read_lab = np.full((read_number), 0, dtype=int)
	read_lab_num = 0

	# --------------------get cvg-------------------

	if bubble_mode:
		empty_count = 0
		i = 0
		narrowed_cvg_list = []
		nearly_narrowed_cvg_list = []
		real_narrowed_cvg_list = []
		bubbles = {}
		ratios = {}
		n_min = (0, 100000)
		rn_min = 100000
		rn_min_index = -1
		ratio_max = -1
		ratio_max_index = -1
		#for i in range(0, intermit_matrix_info.real_narrowed_matrix.shape[1]):
		for i in range(21562,25384):
			tmp = np.squeeze(intermit_matrix_info.real_narrowed_matrix.getcol(i).toarray())
			tmp_count = np.bincount(tmp)[1:]
			real_narrowed_cvg_list.append(sum(tmp_count))
			if rn_min > sum(tmp_count) > 0:
				rn_min = sum(tmp_count)
				rn_min_index = i

			tmp1 = np.squeeze(intermit_matrix_info.narrowed_matrix.getcol(i).toarray())
			tmp_count1 = np.bincount(tmp1)[1:]
			narrowed_cvg_list.append(sum(tmp_count1))
			if n_min[1] > sum(tmp_count1) > 0:
				n_min = (i, sum(tmp_count1))

			bub = narrowed_cvg_list[-1] - real_narrowed_cvg_list[-1]
			tmp2 = np.squeeze(intermit_matrix_info.nearly_real_narrowed_matrix.getcol(i).toarray())
			tmp_count2 = np.bincount(tmp1)[1:]
			nearly_narrowed_cvg_list.append(sum(tmp_count2))

			ratio = sum(tmp_count2)/sum(tmp_count)

			if ratio_max < ratio:
				ratio_max = ratio
				ratio_max_index = i
			bubbles.update({i: bub})
			ratios.update({i: ratio})
			if sum(tmp_count) ==0 or sum(tmp_count1) == 0:
				empty_count += 1

		#rn_min_index = 21542
		#rn_min = 486
		print(rn_min, rn_min_index, n_min)
		print(ratio_max,ratio_max_index,real_narrowed_cvg_list[ratio_max_index-21562], nearly_narrowed_cvg_list[ratio_max_index-21562])
		sg_range = list(range(ratio_max_index - 150, ratio_max_index + 150))
		rn_sg_extend = 75
		rn_covered_reads = []
		side_bubble_reads = []
		for rn_i in sg_range:
			rn_tmp = np.squeeze(intermit_matrix_info.real_narrowed_matrix.getcol(rn_i).toarray())
			nearly_rn_tmp = np.squeeze(intermit_matrix_info.nearly_real_narrowed_matrix.getcol(rn_i).toarray())
			# print("rn sg-i", np.bincount(rn_tmp)[1:])
			# print("n sg-i",np.bincount(n_tmp)[1:])
			rn_readnum = np.nonzero(rn_tmp)[0]

			# real_before = copy.deepcopy(intermit_matrix_info.real_narrowed_read)
			tmp_rn_read_list, intermit_matrix_info.real_narrowed_read = bm.collecting_bubbles(list(rn_readnum),
																				   intermit_matrix_info.real_narrowed_read)

			# test_thinread_set, test_narrowed_read = bm.marking_byid(list(thin_readnum),rn_min,real_before,set({}))

			nearly_rn_readnum = np.nonzero(nearly_rn_tmp)[0]

			tmp_nearly_rn_read_list, intermit_matrix_info.nearly_reaL_narrowed_read = bm.collecting_bubbles(list(nearly_rn_readnum),
																		 intermit_matrix_info.nearly_real_narrowed_read)

			rn_covered_reads.extend([x[3] for x in tmp_rn_read_list])
			rn_covered_reads.extend([x[3] for x in tmp_nearly_rn_read_list])
			#side_bubble_reads.extend(tmp_rn_read_list)
			side_bubble_reads.extend(tmp_nearly_rn_read_list)
			# print(len(whole_col), len(thin_readnum))
		print(len(rn_covered_reads)," reads from sg_range")
		del rn_i
		for rn_i in [x for x in range(sg_range[0]-rn_sg_extend, sg_range[-1]+rn_sg_extend) if x not in set(sg_range)]:
		#for rn_i in [x for x in range(21563,25384) if x not in set(sg_range)]:
			rn_tmp = np.squeeze(intermit_matrix_info.real_narrowed_matrix.getcol(rn_i).toarray())
			rn_readnum = np.nonzero(rn_tmp)[0]

			# real_before = copy.deepcopy(intermit_matrix_info.real_narrowed_read)
			tmp_rn_read_list, intermit_matrix_info.real_narrowed_read = bm.collecting_bubbles(list(rn_readnum),
												intermit_matrix_info.real_narrowed_read)
			rn_covered_reads.extend([x[3] for x in tmp_rn_read_list])

			side_bubble_reads.extend(tmp_rn_read_list)
		del rn_i
		print(len(rn_covered_reads)," reads from sg_range_extend")
		rn_read_set = set(rn_covered_reads)
		print( len(side_bubble_reads), "  reads for side bubbles ")
		test_flag = format(side_bubble_reads[0][1],'b')[::-1]
		print(test_flag)
		if test_flag[4] == "1":
			if test_flag[6] == "1":
				rc_file = readfile1
			else:
				rc_file = readfile2
		else:
			if test_flag[6] == "1":
				rc_file = readfile2
			else:
				rc_file = readfile1
		bm.get_bubble_reads(readfile1,readfile2,side_bubble_reads,out_dir, rc_file)

		'''
		with open("thin_segments_reads_IDs.txt", "w+") as tf:
			# for rn_num in thin_readnum:
			#	tri = intermit_matrix_info.real_narrowed_read[rn_num]
			#	tf.write(tri[0]+" "+str(tri[1])+" "+str(tri[2])+" "+tri[3]+"\n")
			for tr in thinread_set:
				tf.write(tr + " ")
		print(set(nearly_rn_covered_reads) - thinread_set)
		with open("thin_column_reads_ID.txt", "w+") as tcf:
			# for n_num in whole_col:
			#	tcri = intermit_matrix_info.narrowed_read[n_num]
			# if tcri[0] not in thinread_set:
			#	tcf.write(tcri[0]+" "+str(tcri[1])+" "+str(tcri[2])+" "+tcri[3]+"\n")
			for tcr in nearly_rn_covered_reads:
				if tcr not in thinread_set:
					tcf.write(tcr + " ")'''

		bubble_list = [v for v in bubbles.values()]
		print(len(bubble_list), len(narrowed_cvg_list), len(narrowed_cvg_list))
		print(statistics.median(bubble_list), statistics.median(bubble_list), max(bubble_list))
		with open("side_bubble_extract.sam","w+") as sbf:

			for line in side_bubble_reads:
				sbf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
		with open("bubbles.txt","w+") as bf:
			for k in bubbles.keys():
				bf.write(str(k)+": "+str(bubbles[k])+", ")
		exit()

		# ----------------------plotting bubbles----------
		'''x = np.asarray([k for k in bubbles.keys()])
		y = np.asarray([v for v in bubbles.values()])
		px = 1 / plt.rcParams['figure.dpi']
		#fig, ax = plt.subplots(figsize=(9000 * px, 900 * px))
		fig = plt.figure(figsize=(100,10))
		print(y)
		print(x)
		plt.bar(x, y, align='center')
		plt.xticks(x)
		plt.xlabel('position')
		plt.ylabel('bubble size')
		#ax.set(xlim=[3800, 20])
	
		plt.xlim(0, 3800)
		plt.ylim(0 ,80)
	
		# plt.xticks([x for x in range(100,153)])
		# ax.xaxis.set_tick_params(width=5)
		# fig, axs = plt.subplots(1,2, sharey=True, tight_layout=True)
		#plt.figure(figsize=(4000,200),dpi=20)
		print(fig)
		# axs[0].hist(x, bins=20)
		# axs[1].hist(y,bins=20)
		plt.show()'''
	else:
		cvg_list = []
		i= 0
		for i in range(0,intermit_matrix_info.real_narrowed_matrix.shape[1]):
			tmp = np.squeeze(intermit_matrix_info.real_narrowed_matrix.getcol(i).toarray())
			tmp_count = np.bincount(tmp)[1:]
			cvg_list.append(sum(tmp_count))

		# cvg = statistics.median(cvg_list)
		cvg_count = Counter(cvg_list)
		#s_cvg_list = sorted(cvg_list)
		#cvg = s_cvg_list[int(0.01 * len(cvg_list))]
		cvg = 26
		# cvg = s_cvg_list[0]
		if cvg == 0:
			cvg = 1
		# cvg = 1
		#print(cvg_count)
		print(len(cvg_list), "cvg is set to ", cvg, "medium coverage", statistics.median(cvg_list), "minimum", min(cvg_list),
			  "average ",
			  statistics.mean(cvg_list))

	# print({k: v for k, v in sorted(cvg_count.items(), key=lambda x: x[0])})

	# -------------------------compute probabilities of each base
	# count unknow bases in the updated sequence

	non_update = 0
	insufficient = []
	emptycol = []
	sufficient = []
	# domin_values = {}
	marked_id_list = []
	# real_len = bm.get_real_len()
	# for x1, x1v in enumerate(real_len):
	#	domin_values.update({x1: [1] * x1v})

	i = start
	i1 = start  # i1 is cursor for ref
	ins_count = 0
	for i in range(start, add_matrix.shape[1]):
		# for i in range(1,gene_length+1):
		# go through each columntm

		tmp = np.squeeze(add_matrix.getcol(i).toarray())  # matrix

		tmp_count = np.bincount(tmp)[1:]

		# col_num store the value of total_coverage_col
		col_num = sum(tmp_count)
		# maximum letter's number is max_l + 1

		if (sum(tmp_count) == 0):
			# 0 means this base is unknown
			# use original ref sequence if unknown
			emptycol.append(i)

			if i < len(ref) + len(insertion_columns_list) - 2:

				if ins_count < len(insertion_columns_list):
					if i <= insertion_columns_list[ins_count]:
						ins_num = ins_count
					else:
						ins_count += 1
						ins_num = ins_count
				else:
					ins_num = len(insertion_columns_list)

				if (i in insertion_columns_list):
					Seq_let.append(dict_tole[5])
					Seq_tru.append(5)
					print("shouldn't insert? ", i, insertion_columns_list)
					print(i1)
					print(tmp)
					print(tmp_count)
				else:
					non_update = non_update + 1

					Seq_tru.append(7)
					Seq_let.append(" ")
					i1 += 1
			else:
				Seq_let.append(dict_tole[5])
				Seq_tru.append(5)
			continue
		else:
			covered_reads = list(np.nonzero(tmp)[0])
			if i % 1000 == 0:
				tmp_marked, intermit_matrix_info.real_narrowed_read = bm.marking(covered_reads, cvg,
																				 intermit_matrix_info.real_narrowed_read, col=i)
			else:
				tmp_marked, intermit_matrix_info.real_narrowed_read = bm.marking(covered_reads, cvg,
																			 intermit_matrix_info.real_narrowed_read)
			marked_id_list.extend(tmp_marked)

			nonzero_num = np.nonzero(tmp_count)[0].shape[0]
			max_l = np.argmax(tmp_count)
			max_num = tmp_count[max_l]
			# print(tmp)
			# print(tmp_count)
			# print(i, max_l, max_num, i in insertion_columns_list)

			# if(3175<=i<=3180):
			#    print(i, tmp_count, ref[i], dict_tole[max_l+1], max_l)
			#    print(tmp[tmp == max_l])

			if i >= len(ref) + len(insertion_columns_list) - 2:
				Seq_tru.append(max_l + 1)
				Seq_let.append(dict_tole[max_l + 1])

			else:
				# update letter of sequence

				if i in insertion_columns_list:

					Seq_tru.append(max_l + 1)
					Seq_let.append(dict_tole[max_l + 1])
					cor_record.update({i: "*|" + dict_tole[max_l + 1]})
				else:

					try:
						ref_le = ref[i1 - start]
					except:
						print("ref index error")
						print(i1, len(insertion_columns_list))
						print(i, len(ref))
						exit(1)
					ref_pos = i1 - start
					if dict_tole[max_l + 1] != ref[i1 - start]:
						# print(tmp_count, dict_tole[max_l+1], ref[i-start])
						if sum(tmp_count) >= count_threshold:

							Seq_tru.append(max_l + 1)
							Seq_let.append(dict_tole[max_l + 1])
						else:
							non_update += 1
							insufficient.append((i, sum(tmp_count)))
							Seq_tru.append(dict_tonu[ref[i1 - start]])
							Seq_let.append(ref[i1 - start])
					else:
						non_update += 1
						sufficient.append(i)
						Seq_tru.append(max_l + 1)
						Seq_let.append(dict_tole[max_l + 1])

					i1 += 1

	print("majortiy of contigs acquired ")
	print("time is ", time.time() - start_time, "step ", time.time() - prev_time)
	prev_time = time.time()
	# remove deleted sequence

	seq_length = len(Seq_let)
	let = gene_length
	while let < add_matrix.shape[1]:
		if let >= seq_length:
			break
		tmp = np.squeeze(add_matrix.getcol(let).toarray())
		tmp_count = np.bincount(tmp)[1:]
		if let > len(ref) + len(insertion_columns_list) - 2:

			Seq_let.pop(let)
			seq_length -= 1

		else:
			let += 1
	# if add_matrix.shape[1] < gene_length or len(Seq_let) < gene_length:
	#	unknown_tail_len = gene_length - len(Seq_let)
	#	Seq_let.extend(ref[maxposition + 1:])

	#	print(len(emptycol), "columns before the end of coverage have gaps")
	#	emptycol.extend([c for c in range(add_matrix.shape[1], len(ref))])
	#	print(unknown_tail_len, "bases at the end has no covering reads, use reference sequence's bases from",
	#		  maxposition + 1)

	Seq_let = list(filter(lambda x: x != "-" or x != " ", Seq_let))
	Seq_let_str = "".join(Seq_let)
	Seq_let_str = Seq_let_str.replace("-", "")
	Seq_let_str = Seq_let_str.replace(" ", "")
	# run megahit here

	print("got the contig , time ", time.time() - start_time, " step ", time.time() - prev_time)
	print(non_update, " bases remained the same as reference, ", len(insufficient) + len(emptycol),
		  " positions not updated becasue of insufficient information")
	# print(insufficient, "columns with insufficient depths (column number, base amount)")
	# print(len(sufficient))
	print(len(emptycol), "cols are empty")
	with open(out_dir + "gaps.txt", "w+") as ef:
		ef.write(str(emptycol))
	# corrected reads in csr

	# pdb.set_trace();
	# store contigs.fa
	W_seq_file = args.write_file_seq
	seq_file = open(out_dir + W_seq_file, 'w')
	seq_file.write('>updated_sequence\n')
	# seq_file.write('%s\n' % ("".join(Seq_let)))
	seq_file.write(Seq_let_str + "\n")
	seq_file.close()
	if len(Seq_let_str) < 25000:
		print("too much gaps to form a sequence", len(emptycol), len(Seq_let), gene_length)
		if len(emptycol) < 3000:
			print(emptycol)
		bm.get_ori_half(readfile1, readfile2, marked_id, out_dir, real_narrowed)

		exit(5)
	print("error correction finished!")
	print("correction records are in records.txt")
	# print(Seq_let)
	if not first_round:
		with open(out_dir + "records.txt", 'w') as rf:
			c = 1
			sorted_record = sorted(cor_record.keys())
			for k in sorted_record:
				if k > len(ref) + len(insertion_columns_list):
					continue
				eol = "\n" if c % 15 == 0 else "   "
				rf.write(str(k) + ": " + cor_record.get(k) + eol)
				c += 1
		print("%d bases are updated in output sequence" % int(gene_length - non_update))
	if find_sub != "no":
		exit()
	print("prepare for next round, save all marked ID in excluded_IDs.txt")
	marked_id.update(marked_id_list)
	bm.write_new_extract(intermit_matrix_info.real_narrowed_read, marked_id, out_dir, round_num)
	bm.get_ori_half(readfile1, readfile2, marked_id, out_dir, intermit_matrix_info.narrowed_read)
# corrected_contig
