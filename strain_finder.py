import os, random, time, sys
import shlex
import statistics
import subprocess

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

import identify_verify

id_field = 0
flag_field = 1
index_field = 2
read_field = 3
cigar_field = 4
freq_field = 5
misp_field = 6
gene_length = 29891
'''
positions for each protein region
'''
endpoints = '''21563..25384
25393..26220
26245..26472
26523..27191
27202..27387
27394..27759
27756..27887
27894..28259
28274..29533
29558..29674'''
dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))
protein_loc = []  # (265,21555)

def global_init():
	"""
	initializing global variables for handling sam file fields
	"""
	global id_field, flag_field, index_field, read_field, cigar_field, freq_field, \
		misp_field, protein_loc, dict_tonu, dict_tole

	endpoints_list = endpoints.split("\n")
	for se in endpoints_list:
		start = se.split("..")[0]
		end = se.split("..")[1]
		protein_loc.append((int(start) - 1, int(end) - 1))

	#print(protein_loc)



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
	parser.add_argument("--match_l", type=float, default=0.7)
	# parser.add_argument("--first_r", type=bool, default=True)
	parser.add_argument("--r1_file", type=str)
	parser.add_argument("--r2_file", type=str, default="")
	parser.add_argument("--round", type=int, default=1)
	parser.add_argument("--excluded_IDs", type=str, default="excluded_IDs.txt")
	parser.add_argument("--find_sub", type=str, default="no")
	parser.add_argument("--bubble_mode", type=str, default="False")
	parser.add_argument("--brute_force", type=str,default="True")
	parser.add_argument("--misp_mode", type=str, default="False")
	parser.add_argument("--check_gap", type=str,default="False")
	parser.add_argument("--output_dir", type=str,default="")
	parser.add_argument("--region_break",type=str,default="")
	parser.add_argument("--gap_threshold",type=int,default=1)
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


if __name__ == "__main__":
	# -------------------- main function is here-----------------------------------------------------
	start_time = time.time()

	global_init()
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



	# setting
	round_num = args.round

	gene_length = args.gene_L
	read_length = args.read_L
	read_number = r.shape[0]
	count_threshold = args.count_l
	start = 0
	domin_threshold = args.domin_l
	match_limit = float(args.match_l)
	threshold = 0.025
	# first_round = args.first_r
	readfile1 = args.r1_file
	readfile2 = args.r2_file
	find_sub = args.find_sub


	if find_sub != "no":
		if args.output_dir == "":
			out_dir = os.path.basename(args.ref) + "_output/"
		else:
			out_dir = args.output_dir +"/"
	else:
		out_dir = "round_" + str(round_num) + "_output/"
	print(out_dir)

	if round_num == 1:
		first_round = True
	else:
		first_round = False

	if args.misp_mode == "True":
		misp_mode = True
	else:
		misp_mode = False

	row_l = []
	col_l = []
	val_l = []

	maxposition = 0  # debug variable
	maxindex = 0  # debug variable
	labindexes = {}  # debug
	# read reference genome
	ref = ""
	changed = 0  # debug
	with open(args.ref, 'r') as refg:
		for line in refg:
			if ">" not in line:
				ref += line.strip().strip(" ")

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

	start_time = time.time()

	if args.check_gap == "False":
		gaps = []
		cols = []
		read_list,read_freq,real_narrowed = bm.dup_read_sam(R_file,ref)
		initial_matrix_info = bm.matrix_from_readlist(real_narrowed, match_limit, marked_id, True,target="raw")
		for column in range(400, initial_matrix_info.narrowed_matrix.shape[1]):
			tmp = np.squeeze(initial_matrix_info.narrowed_matrix.getcol(column).toarray())
			tmp_count = np.bincount(tmp)[1:]
			cols.append(column)
			if sum(tmp_count) < args.gap_threshold:
				gaps.append(column)
		print("cols for checking gaps", min(cols), max(cols))
		if len(gaps) > 0:
			if args.check_gap != "false":
				print("only checking gaps")
				with open("gaps.txt", "w+") as gapf:
					for gap_col in gaps:
						gapf.write(str(gap_col) + ",")

			for gap in gaps:
				print(gap, end=",")
			print()
			# exit()
			exit(-4)
		else:
			print("no gaps")


	else:

		read_list = bm.read_sam(R_file)
		#read_list,freq = bm.dup_read_sam(R_file)
		print("read_sam takes", time.time()-start_time)
		initial_matrix_info = bm.matrix_from_readlist(read_list, match_limit, marked_id)

	print("matrix_initialization", time.time()-start_time)



	#if args.brute_force == "True":

	real_narrowed, paired_real_narrowed, nearly_real_narrowed, potential_mutated = bm.narrow_reads(ref,
		initial_matrix_info.narrowed_read,out_dir, True)

	#real_narrowed, paired_real_narrowed, nearly_real_narrowed = bm.narrow_reads(ref, initial_matrix_info.narrowed_read, out_dir, False)
	print("narrowed_read",time.time()-start_time)

	intermit_matrix_info = bm.matrix_from_readlist(real_narrowed, match_limit, marked_id,False,initial_matrix_info,"real_narrowed")
	print("real_narrowed_matrix",time.time()-start_time)
	#intermit_matrix_info = bm.matrix_from_readlist(paired_real_narrowed, match_limit, marked_id,False,initial_matrix_info,"real_narrowed")
	intermit_matrix_info = bm.matrix_from_readlist(nearly_real_narrowed,match_limit,marked_id,False,intermit_matrix_info,"nearly_real_narrowed")
	print("nearly_real_narrowed_matrix",time.time()-start_time)
	#intermit_matrix_info = bm.matrix_from_readlist(paired_real_narrowed, match_limit, marked_id)
	intermit_matrix_info.real_narrowed_read = real_narrowed
	#intermit_matrix_info.real_narrowed_read = paired_real_narrowed
	intermit_matrix_info.nearly_real_narrowed_read = nearly_real_narrowed
	#intermit_matrix_info.narrowed_read = copy.deepcopy(initial_matrix_info.narrowed_read)
	#intermit_matrix_info.narrowed_matrix = initial_matrix_info.real_narrowed_matrix.copy()


	#------------------------------------------brute force mode----------

	del initial_matrix_info


	if args.brute_force == "True":
		print("brute_force mode")
		gap_range = range(100, 28000)
		print(len(potential_mutated),"reads not fully matched")


		mutated_read_freq = Counter([ x[read_field] for x in potential_mutated])

		sorted_mutated_read = sorted(potential_mutated,key=lambda x : (mutated_read_freq.get(x[read_field],int(x[index_field]))),reverse=True)


		with open("mutated_read_freq.txt","w+") as mrfreq:
			for k,v in sorted(mutated_read_freq.items(),key= lambda x:x[1]):
				if v > 0:
					mrfreq.write(k+": "+str(v)+"\n")
		with open(out_dir+"mutated_read_freq.txt","w+") as mrfreq:
			for k,v in sorted(mutated_read_freq.items(),key= lambda x:x[1]):
				if v > 0:
					mrfreq.write(k+": "+str(v)+"\n")
		for tr in sorted_mutated_read:
			tr.append(mutated_read_freq[tr[read_field]])


		#print({k:v for k,v in mutated_read_freq.items() if v > 1 })
		#print(sorted_mutated_read)
		smr_index = 0
		misP =[]
		added_read = set({})
		reduced_sorted_mutated_read = []
		if args.region_break != "":
			misPs, misP_source, misP_reads, sorted_mutated_read = identify_verify.get_misp(ref, sorted_mutated_read)
			endpoints_list = []
			with open(args.region_break,"r") as ppf:
				for line in ppf:
					fields = line.split(":")
					endpoints_list.append(fields[1].strip())

			#endpoints_list = endpoints.split("\n")
			protein_loc = []  # (265,21555)
			for se in endpoints_list:
				start = se.split("..")[0]
				end = se.split("..")[1]
				protein_loc.append(range(int(start) - 1, int(end) - 1))

			for smr_index in range(len(sorted_mutated_read)):
				if sorted_mutated_read[smr_index][read_field] in added_read:

					continue
				else:
					added_read.add(sorted_mutated_read[smr_index][read_field])
				#print(smr_index,sorted_mutated_read[smr_index])
				tmp_misP = []
				smr_read_index = sorted_mutated_read[smr_index][index_field]-1
				insertion = False
				if "I" not in sorted_mutated_read[smr_index][cigar_field]:
					#get misp list for this read
					tmp_misP = [list(x.items())[0][0] for x in sorted_mutated_read[smr_index][misp_field]]

				else:
					#currently not handle insertion reads
					insertion = True

					#continue
				if not insertion:
					tmp_misP_set = set(tmp_misP)
					stop_con = False
					for mp in tmp_misP:
						found = False

						#for protein_pos in protein_loc:
						for r in protein_loc:
							# move bases inside this protein region only

							if mp in r and mp < r.stop-2:
								found = True
								# location = int((pos-r.start)/3)
								index = (mp - r.start) % 3

								ref_codon = ref[mp - index:mp + (3 - index)]
								if sorted_mutated_read[smr_index][read_field][mp-smr_read_index] != "-":
									#need to consider multiple misp together
									mut_codon = ""

									# mut_codon = ref_codon[0:index] + sorted_mutated_read[smr_index][read_field][mp-smr_read_index] + ref_codon[index + 1:]
									for i in range(mp-index,mp+(3-index)):
										# if the adjacent spots are also misPs, will use the read's adjecent spot instead of ref

										if i < len(ref):
											if i in tmp_misP_set and 0 <= i-smr_read_index < len(sorted_mutated_read[smr_index][read_field]):
												mut_codon += sorted_mutated_read[smr_index][read_field][i-smr_read_index]
											else:

												mut_codon += ref[i]


									if mut_codon == "TAA" or mut_codon == "TGA" or mut_codon == "TAG":
										print("stop codon produced by",sorted_mutated_read[smr_index], tmp_misP)
										stop_con = True
										break
								#else:

									#mut_codon = ref_codon[mp - index:mp] + ref[mp + 1:pos + (3 - index) + 1]

								# avoid stop condon produced by substitution
							if stop_con or found:
								break
					if stop_con:
						continue
					reduced_sorted_mutated_read.append(sorted_mutated_read[smr_index])
				else:
					#insertion read included

					reduced_sorted_mutated_read.append(sorted_mutated_read[smr_index])


			del smr_index
			del smr_read_index
		else:
			for read in sorted_mutated_read:
				if read[read_field] not in added_read:
					reduced_sorted_mutated_read.append(read)
					added_read.add(read[read_field])



		#assemble

		pos_sorted_mutated_read =copy.deepcopy(reduced_sorted_mutated_read)
		pos_sorted_mutated_read = identify_verify.fix_s_pos(pos_sorted_mutated_read,restore=True)
		#exclude reads with N
		with open(out_dir+"sub_read_candidate.sam","w+") as candf:
			for line in pos_sorted_mutated_read:
				#print(line)
				if "N" not in line[read_field]:
					candf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + " " + str(mutated_read_freq[line[3]]) + "\n")

		rn_cvg_list= []
		nearly_rn_cvg_list = []
		print(intermit_matrix_info.real_narrowed_matrix.shape,intermit_matrix_info.nearly_real_narrowed_matrix.shape)
		for i in gap_range:
			if i >= intermit_matrix_info.real_narrowed_matrix.shape[1]:
				break
			tmp = np.squeeze(intermit_matrix_info.real_narrowed_matrix.getcol(i).toarray())
			tmp_count = np.bincount(tmp)[1:]
			rn_cvg_list.append(sum(tmp_count))
			if i < intermit_matrix_info.nearly_real_narrowed_matrix.shape[1]:
				tmp2 = np.squeeze(intermit_matrix_info.nearly_real_narrowed_matrix.getcol(i).toarray())
				tmp_count2 = np.bincount(tmp2)[1:]
				nearly_rn_cvg_list.append(sum(tmp_count2))

		with open("nearly_real_narrowed_cvg.txt","w+") as near_rncf, open("real_narrowed_cvg.txt", "w+") as rncf:
			for i,v in enumerate(nearly_rn_cvg_list):
				near_rncf.write(str(i)+": "+str(v)+", ")
				rncf.write(str(i)+": "+str(rn_cvg_list[i])+", ")


		exit()



	else:
		del potential_mutated

	#assert intermit_matrix_info.real_narrowed_matrix.shape[1] == intermit_matrix_info.narrowed_matrix.shape[1]
	print("narrowed_reads", len(intermit_matrix_info.narrowed_read), "real narrowed",
		  len(intermit_matrix_info.real_narrowed_read), "nearly real narrowed", len(intermit_matrix_info.nearly_real_narrowed_read))
	print("narrowed_shape", intermit_matrix_info.narrowed_matrix.shape, "real narrowed shape",
		  intermit_matrix_info.real_narrowed_matrix.shape, "nearly real narrowed shape", intermit_matrix_info.nearly_real_narrowed_matrix.shape)

	#del initial_matrix_info
	# print(intermit_matrix_info.narrowed_matrix.shape[1], intermit_matrix_info.real_narrowed_matrix.shape[1])

	maxposition = intermit_matrix_info.max_shape[1]
	intermit_matrix_info = bm.build_insertion(intermit_matrix_info, count_threshold)
	#print("real_narrowed matrix formed ", intermit_matrix_info.real_narrowed_matrix.shape, "time is ", time.time() - start_time)

	add_matrix = intermit_matrix_info.narrowed_matrix.copy()
	insertion_columns_list = intermit_matrix_info.insertion_columns_list
	print("insertion matrix set, time", time.time() - start_time)
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
