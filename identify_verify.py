import collections
import copy
import csv
import os
import random
import subprocess
import sys, re, statistics, itertools as its, multiprocessing as mp, numpy as np, scipy.sparse as sp, \
	matplotlib.pyplot as pit
import argparse, math
import editdistance
import matplotlib.pyplot as plt
import scipy.interpolate
import build_matrix as bm
import strain_finder as st_find

# the misp base in change_str starts from misp_str[2]
change_base_bit = 2
def fix_s_pos(subbed_reads,restore=False):
	new_subbed_reads = []
	for ir, fields in enumerate(subbed_reads):
		tmp_read = copy.deepcopy(fields)

		index = int(fields[st_find.index_field]) - 1
		cigar = fields[st_find.cigar_field]
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
			if not restore:
				if index - blk_pos[0] < 0:

					start_pos = blk_pos[0] - index
					tmp_length -= start_pos
					index = 0
				else:
					index = index - blk_pos[0]

				tmp_read[st_find.index_field] = index + 1
			else:
				if index > 0:

					tmp_read[st_find.index_field] += blk_pos[0]
			#subbed_reads[ir][st_find.index_field] = index + 1
		new_subbed_reads.append(tmp_read)


	return new_subbed_reads

def if_clash(read,covered_pos):
	"""
	determine if the read has conflict with existing substitution
	:param read:
	:param covered_pos: dict for bases at pos
	:return: bool
	"""
	for rei, base in enumerate(read[st_find.read_field]):
		r_start = int(read[st_find.index_field]) - 1
		if r_start + rei in covered_pos.keys():
			if read[st_find.read_field][rei] != covered_pos[r_start + rei]:

				print("clash at", r_start + rei, covered_pos[r_start + rei], read[st_find.read_field][rei])
				return True
	return False

def if_reject(read,matrix,read_freq,misP_reads):
	"""
	find if a read should be rejected for lacking misp support
	:param misp_list: misps of a given read {pos:["ref|change"]
	:param matrix: csc_matrix for fining support bases
	:param read_freq: frequency of reads
	:param misP_reads: reads that support a misp at a position
	:return: bool
	"""
	misp_list = read[st_find.misp_field]
	for rm in misp_list:
		rm_tuple = list(rm.items())[0]
		pos = rm_tuple[0]
		change = rm_tuple[1][0].split("|")[1]
		tmp = np.squeeze(matrix.getcol(pos).toarray())
		tmp_count = np.bincount(tmp)[1:]

		support_index = np.where(tmp == st_find.dict_tonu[change[0]])
		support_count = support_index[0].shape[0]
		if support_count == 0 and read_freq.get(read[st_find.read_field]) == 1 and pos not in misP_reads.keys():
			return True
	return False

def verify_seq_support(temp_ref,batch,readfile1, readfile2):
	"""
	Find if the reads support a sequence with no gaps by running alignment
	:param temp_ref: sequence to be tested
	:param batch: batch number
	:param readfile1:
	:param readfile2:
	:return: gap list
	"""
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

		return gap_list
	# gap_list = gap_list.split(",")[:-1]
	# ref = ref[:read_index] + ref_seq[read_index:read_index+len(read[3])] + ref[read_index + len(read[3]):]
	# subbed_read.remove(read)
	else:
		return []

def identify_strain(ref_file, strain, candidate_sam, readfile1, readfile2):
	"""produce and verify new strain, also replace rejected reads
	@param ref_file fasta file of reference sequence
	@param strain integer number of strain
	This function identifies new strains and verify , it iteratively goes through the list
	of candidate reads, check if substitution will result in sequence not supported by sequencing reads
	or cause stop codon that separates a protein from the middle
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
	with open(candidate_sam, "r") as samf:
		for line in samf:
			if line.split(" ")[st_find.read_field] not in prev_readset:
				fields = line.strip().split(" ")
				if "N" in fields[st_find.read_field]:
					continue
				readlist.append(fields)
	# readlist = fix_s_pos(readlist)
	print(len(readlist), "candidate reads")
	ref_seq = ""
	with open(ref_file, "r") as rf:
		for line in rf:
			if line[0] != ">":
				ref_seq += line.strip()
		ref = ref_seq

	misPs,misP_source,misP_reads,new_readlist = get_misp(ref, readlist)
	original_readlist = copy.deepcopy(readlist)
	readlist = new_readlist
	batch = 0

	# while len(readlist) > 0:
	subbed_read = []
	misPs = []
	ref = ref_seq
	covered_pos = {}
	# remove reads that cause early stop condon


	for batch, read in enumerate(readlist):
		overlap = False
		read_index = int(read[st_find.index_field]) - 1
		read_misPs = [list(x.items())[0][0] for x in read[st_find.misp_field]]

		print("curr read", batch, read, "bases covered", len(covered_pos))
		if len(subbed_read) > 0:
			# read clash detection, find all overlapping reads
			# mp-clash, compare all MPs, if the reads only itroduce new MP outise current mp ranges
			clash = if_clash(read,covered_pos)
			if clash:
				batch += 1
				continue

		temp_ref = ref[:read_index] + read[st_find.read_field] + ref[read_index + len(read[st_find.read_field]):]

		if verify_seq_support(temp_ref,batch,readfile1,readfile2):

			subbed_read.append(original_readlist[batch])
			# update covered positions with base in covered_pos
			read_start = int(read[st_find.index_field]) - 1
			for ind, base in enumerate(read[st_find.read_field]):

				if read_start + ind in covered_pos.keys():

					if covered_pos[read_start + ind] != base:
						print("clash happened for current read", ind, covered_pos[read_start + ind], base)
						exit()
				else:

					covered_pos.update({read_start + ind: base})

			ref = temp_ref
			print(batch, "no gaps for curr read\n")

		batch += 1

	with open("final_strain_" + str(strain) + "_reference.fa", "w+") as bf:
		bf.write(">final_strain_" + str(strain) + "\n")
		bf.write(ref)
	with open("subbed_reads_" + str(strain) + ".sam", "w+") as bf:
		for line in subbed_read:
			bf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + " "+str(line[5])+"\n")
	# subprocess.run("rm -r batch_*",shell=True)
	print(subbed_read)
	for read in subbed_read:
		if read in readlist:
			readlist.remove(read)

	with open("strain_" + str(strain) + "_spike.fa", "w+") as sf:
		sf.write(">strain_" + str(strain) + "_spike" + "\n")
		sf.write(ref[21562:25383])
	strain += 1



def base_compare(base1,base2):
	"""
	compare two bases
	:param base1: base from ref
	:param base2: base from read
	:return: misp format if not equal, else empty
	"""
	if base1 != base2:
		return [base1+"|"+base2]
	else:
		return []

def misp_equal(misp1, misp2):
	"""

	:param misp1: (position(int), [change_str1,change_str2...]) existed misp at position misp[0]
	:param misp2: new misp, only one change_str (position ,[change_str])
	:return: index in of change_str list in {position:[]}
	"""
	#print(misp1)
	for i_mispstr,existed_mispstr in enumerate(misp1[1]):
		# for insertion
		if existed_mispstr[0] == "-" == misp2[1][0]:
			#inserting bases are the same
			if existed_mispstr == misp2[0]:
				return i_mispstr
		# for non-insertion
		else:

			if existed_mispstr == misp2[1][0]:
				return i_mispstr
	return None


def get_misp(ref, sub_read, printing=True, fix_pos=True):
	"""
	handle getting misp in this func
	:param ref: reference sequence file
	:param sub_read: list of reads
	:param printing: True to print to stdout
	:return: misPs {position: ["ref_base|misp_base"]}, misP_source {position:{change_str:[SRR_ID]}}, misP_reads {position:{change_str:[reads]}}
			read end with [{position: ["ref_base|misp_base"]}]
	"""
	if fix_pos:
		sub_read = fix_s_pos(sub_read)
	read_with_misp = []
	misPs = {}
	misPs_count = {}
	misPs_source = {}
	misP_reads = {}
	pos = {}

	read_freq = {}
	with open("mutated_read_freq.txt", "r") as mf:
		for line in mf:
			freq_fields = line.strip().split(": ")
			read_freq.update({freq_fields[0]: int(freq_fields[1])})

	insertion_reads = {}
	included_read_num = 0
	for read_num,read in enumerate(sub_read):
		temp_misp = []
		read_index = int(read[st_find.index_field]) - 1
		read_str = read[st_find.read_field]
		cigar = read[st_find.cigar_field]
		cigar_str = re.findall(r"[0-9]+[MIDSH]", cigar)
		blk_pos = []
		blk_type = []
		blk_length = []
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
			blk_length.append(int(m.group(1)))
			blk_pos.append(bl)
			ini = bl
			# if iread[4]=="140M2D10M":
			#	print(block, int(m.group(1)))
			base_length += int(m.group(1))
			if bt != "I" and bt != "H":
				tmp_length += int(m.group(1))  # get read length without counting clipped bases

		#correct read_index for soft-clipping at beginning, if index already handled in build_matrix, not doing ,
		#if blk_type[0] == "S":
		#	if read_index - blk_pos[0] < 0:
		#		start_pos = blk_pos[0] - read_index
		#		tmp_length -= start_pos
		#	else:
		#		read_index = read_index - blk_pos[0]

		c = 0
		inserts = []
		begin = 0 # insertion begin

		for j in range(0, blk_pos[-1]):  # change here to fill the blank with 0?
			# compare base directly, soft clipping at beginning has been accounted
			if blk_type[c] == "M" or blk_type[c] == "S":

				curr_misp= base_compare(ref[read_index + j],read_str[j])

				if len(curr_misp)>0:
					temp_misp.append({read_index + j:curr_misp})

			# handle insertion
			elif blk_type[c] == "I":
				inserts.append(st_find.dict_tonu[read_str[j]])

			# should have been inserted "-" in missing position
			elif blk_type[c] == "D" or blk_type[c] == "H":
				temp_misp.append({read_index + j: [ref[read_index  + j] + "|" + read_str[ j]]})
			if j == blk_pos[c] - 1:  # update start and c, put inserts into hashtable

				if blk_type[c] == "I":
					insert_str = ""
					for ins_base in inserts:
						insert_str += ins_base

					temp_misp.append({read_index + j:["-|"+insert_str]})
					if read_num in insertion_reads.keys():

						newinsert = copy.deepcopy(insertion_reads.get(included_read_num))
						newinsert.append((read_index+ begin, inserts))
						insertion_reads.update({included_read_num: newinsert})
					else:
						insertion_reads.update({included_read_num: [(read_index + begin, copy.deepcopy(inserts))]})

				begin = blk_pos[c]
				inserts = []

				c += 1
				if c == len(blk_type):
					break

		new_read = copy.deepcopy(read)
		if len(new_read) > st_find.freq_field+1 and type(new_read[st_find.freq_field]) is not int:
			#handle "False" appended to fields
			del new_read[st_find.freq_field]
		elif len(new_read) < st_find.freq_field+1:

			new_read.append(read_freq[new_read[st_find.read_field]])

		new_read.append([])

		# record misp in this blk
		for misp_kv in temp_misp:
			misp = [(k,v) for k,v in misp_kv.items()][0]
			#print(new_read)
			new_read[st_find.misp_field].append(misp_kv)
			if misp[0] not in misPs.keys():
				misPs.update(misp_kv)
				misPs_count.update({misp[0]: [1]})
				misPs_source.update({misp[0]: {misp[1][0].split("|")[1]: [read[st_find.id_field].split(".")[0]]}})
				misP_reads.update({misp[0]: {misp[1][0].split("|")[1]: [read]}})
			else:
				existed_change_index = misp_equal((misp[0],misPs[misp[0]]), misp)
				# found same misp
				if existed_change_index is not None:
					misPs_count[misp[0]][existed_change_index] += 1
					misPs_source[misp[0]][misp[1][0].split("|")[1]].append(read[st_find.id_field].split(".")[0])
					misP_reads[misp[0]][misp[1][0].split("|")[1]].append(read)
				#no same misp
				else:
					misPs[misp[0]].extend(misp[1])
					misPs_count[misp[0]].append(1)
					misPs_source[misp[0]].update({misp[1][0].split("|")[1]: [read[st_find.id_field].split(".")[0]]})
					misP_reads[misp[0]].update( {misp[1][0].split("|")[1]: [read]})


			if misp[0] not in pos.keys():

				pos.update({misp[0]: [read]})
			else:
				pos[misp[0]].append(read)

		read_with_misp.append(new_read)

	# with open("different_bases.txt","w+") as wf:

	for k, v in sorted(pos.items(), key=lambda x: x[0]):

		# if 21562<=k<=25384:
		if k not in misPs.keys():
			print(k,misPs.keys())
		for change in misPs[k]:
			changed_base = change.split("|")[1]
			support_list = [(k, v) for k, v in collections.Counter(misPs_source[k][changed_base]).items() if v > 1]
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

	return misPs, misPs_source, misP_reads,read_with_misp


def verify_misp(ref_file, support_matrix_info, sub_read_file, new_strain_file, candidate_sam):
	"""

	:param ref_file: original ref fa
	:param support_matrix_info: other SRR sam
	:param sub_read_file: subbed_read
	:param new_strain_file: generated strain fa
	:param candidate_sam: sub_read_candidates
	:return: rejected reads, updated misp_conflict, no_candidate to fill
	"""
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

	strain = int(re.search('[0-9]+', os.path.basename(sub_read_file))[0])

	test_misp = {}
	subbed_read = []
	with open(sub_read_file, "r") as subf:
		for line in subf:
			fields = line.strip().split(" ")
			subbed_read.append(fields)

	misPs, misPs_source, misP_reads,subbed_read = get_misp(ori_ref, subbed_read)
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
			read_freq.update({sline[st_find.read_field]: int(sline[st_find.freq_field])})
	candidate_read = fix_s_pos(candidate_read)
	#initial_matrix_info = bm.matrix_from_readlist(bm.read_sam(support_matrix), 0, marked_id=set({}), target="raw")

	#raw_matrix_info = bm.build_insertion(initial_matrix_info, 0)
	# not using insertion now,later change to raw_matrix_info.real_narrowed_matrix
	raw_matrix_info = support_matrix_info
	matrix = support_matrix_info.narrowed_matrix
	#raw_matrix_info.narrowed_read = fix_s_pos(raw_matrix_info.narrowed_read)
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
			changed_base = change.split("|")[1]
			support_index = np.where(tmp == st_find.dict_tonu[change.split("|")[1]])
			support_count = support_index[0].shape[0]
			visited_sindex = set({})
			if support_count > 0:
				for sp in support_index[0]:
					temp_read = raw_matrix_info.narrowed_read[sp]
					# print(temp_read)

					# print(tmp_count)
					# print(ref[pos-5:pos],ref[pos],ref[pos:pos+5])
					# print(change,pos,sp,support_index,support_count)
					read_index = int(temp_read[st_find.index_field]) - 1
					window = 0
					tolerance = 1
					curr_mcount = 0
					curr = pos
					read_end = read_index + len(temp_read[st_find.read_field])
					reached = False
					l_reached = False
					l_window = 0
					l_temp_read = temp_read

					r_reached = False
					r_window = 0
					r_temp_read = temp_read
					encountered_misp = []

					while curr - 1 >= read_index:
						try:
							if ref[curr - 1] == temp_read[st_find.read_field][curr - 1 - read_index]:
								l_window += 1

							else:
								l_reached = True
								break
						except:
							print(curr,temp_read,len(ref))
							raise "index error"
						curr -= 1
					curr_l_line = sp - 1
					if curr_l_line >= 0 and curr - 1 <= raw_matrix_info.narrowed_read[curr_l_line][st_find.index_field] - 1 + len(
							raw_matrix_info.narrowed_read[curr_l_line][st_find.read_field]):
						curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
					# print(curr_l_read)

					else:
						l_reached = True
					#while not l_reached and curr_l_line >= 0:
					for curr_l_line in range(sp-1,0,-1):
						if l_reached:
							break
						if curr_l_line in visited_sindex:
							continue
						else:
							visited_sindex.add(curr_l_line)
						curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
						if curr - 1 < int(curr_l_read[st_find.index_field]) - 1 or curr-1 >=int(curr_l_read[st_find.index_field]) - 1 +len(curr_l_read[st_find.read_field]) :

							if curr_l_line - 1 < 0:
								l_reached = True
								break

							#curr_l_line -= 1
							#curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
							if curr - 1 > int(curr_l_read[st_find.index_field]) - 1 + len(curr_l_read[st_find.read_field]):
								l_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_l_read[st_find.index_field]) - 1 <= pos <= int(curr_l_read[st_find.index_field]) - 1 + len(curr_l_read[st_find.read_field]) - 1 and \
									curr_l_read[st_find.read_field][pos - (int(curr_l_read[st_find.index_field]) - 1)] != change.split("|")[1]:
								# print("read not supporting misp",pos,change)
								# print(curr_l_read)
								# l_reached= True
								#curr_l_line -= 1
								#if curr_l_line < 0:
								#	l_reached = True
								#	break
								#else:
								#	curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
								continue

							# print(curr_l_read)
							continue

						try:
							if ref[curr - 1] == curr_l_read[st_find.read_field][curr - 1 - (int(curr_l_read[st_find.index_field]) - 1)]:
								l_window += 1

							else:
								l_reached = True
								break
						except:
							print(curr-1,curr - 1 - (int(curr_l_read[st_find.index_field]) - 1),curr_l_line,curr_l_read,raw_matrix_info.narrowed_read[curr_l_line])
							raise "index errror"
						# except:
						#    print(curr_l_read,curr-1-(int(curr_l_read[2]-1)))
						#    exit()
						curr -= 1
					l_seg = ref[curr:pos]
					curr = pos

					while curr + 1 <= read_index + len(temp_read[st_find.read_field]) - 1:

						if ref[curr + 1] == temp_read[st_find.read_field][curr + 1 - read_index]:
							r_window += 1
						else:
							r_reached = False
							break

						curr += 1
					curr_r_line = sp + 1
					if curr_r_line <= raw_matrix_info.narrowed_matrix.shape[0] and \
							raw_matrix_info.narrowed_read[curr_r_line][st_find.index_field] - 1 <= curr + 1:

						curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
					# print(curr_r_read)
					else:
						r_reached = True
					#while not r_reached and curr_r_line <= raw_matrix_info.narrowed_matrix.shape[0]:
					for curr_r_line in range(sp+1,raw_matrix_info.narrowed_matrix.shape[0]):
						if r_reached:
							break
						if curr_r_line in visited_sindex:
							continue
						else:
							visited_sindex.add(curr_r_line)
						curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
						if curr + 1 > int(curr_r_read[st_find.index_field]) - 1 + len(curr_r_read[st_find.read_field]) - 1 or\
								curr + 1 < int(curr_r_read[st_find.index_field])-1:

							if curr_r_line + 1 >= raw_matrix_info.narrowed_matrix.shape[0]:
								r_reached = True
								break
							#curr_r_line += 1
							#curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
							# if current read does not cover the curr position (gaps), break
							if curr + 1 < int(curr_r_read[st_find.index_field]) - 1:
								r_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_r_read[st_find.index_field]) - 1 <= pos <= int(curr_r_read[st_find.index_field]) - 1 + len(curr_r_read[st_find.read_field])-1 and \
									curr_r_read[st_find.read_field][pos - (int(curr_r_read[st_find.index_field]) - 1)] != changed_base:
								# print("read not supporting misp",pos,change, curr_r_read)
								# print(curr_r_read)
								#curr_r_line += 1
								#if curr_r_line + 1 >= raw_matrix_info.narrowed_matrix.shape[0]:
								#	r_reached = True
								#	break
								#else:
								#	curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
								# r_reached= True
								continue

							# print(curr_r_read)
							continue
						if ref[curr + 1] == curr_r_read[st_find.read_field][curr + 1 - (int(curr_r_read[st_find.index_field]) - 1)]:
							r_window += 1
						else:
							r_reached = True
							break
						curr += 1
					seg = l_seg + ref[pos:curr + 1]
					window = l_window + r_window + 1
					if seg not in found_segments:
						found_segments.add(seg)
						print("target misp", pos, read_freq.get(misP_reads[pos][changed_base][0][st_find.read_field]))

						print("seg_length", len(seg))

					# print("ref:",ori_ref[pos-l_window:pos+r_window+1])

					# print("seg:",seg)
					# reject this read candidate
					if len(seg) < 50 and len(misP_reads[pos][changed_base]) < 2 and read_freq.get(
							misP_reads[pos][changed_base][0][st_find.read_field]) == 1:
						rejected_read.extend(misP_reads[pos][changed_base])
						# print("rejected ",misP_reads[pos][changed_base])
						reject = True
						break
					# continue

					#print("\n segment size", window, "\n")
					if len(misP_reads[pos][changed_base]) > 1:
						total_freq = 0
						for tr in misP_reads[pos][changed_base]:
							total_freq += read_freq.get(tr[st_find.read_field])
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][changed_base][0][3]), window))
						accepted_misps.update({pos: (change, total_freq, window)})
					else:
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][changed_base][0][3]), window))
						accepted_misps.update({pos: (change, read_freq.get(misP_reads[pos][changed_base][0][st_find.read_field]), window)})

			else:
				if len(misP_reads[pos][changed_base]) < 2 and read_freq.get(
						misP_reads[pos][changed_base][0][st_find.read_field]) == 1:
					rejected_read.extend(misP_reads[pos][changed_base])
					# print("rejected ",misP_reads[pos][changed_base])
					reject = True
					break
				else:
					if len(misP_reads[pos][changed_base]) > 1:
						total_freq = 0
						for tr in misP_reads[pos][changed_base]:
							total_freq += read_freq.get(tr[3])
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][changed_base][0][3]), window))
						accepted_misps.update({pos: (change, total_freq, 0)})
					else:
						# accepted_misps.add((pos, change, read_freq.get(misP_reads[pos][changed_base][0][3]), window))
						accepted_misps.update(
							{pos: (change, read_freq.get(misP_reads[pos][changed_base][0][st_find.read_field]), 0)})

		if reject:
			continue

	if len(rejected_read) > 0:
		print(len(rejected_read), "reads rejected")
		#restore index pos for writing files
		write_sam(rejected_read, "rejected_" + os.path.basename(sub_read_file), "a",True)
		rejected_set = set([x[st_find.read_field] for x in rejected_read])

		d_info = {}
		# for ir,read in enumerate(subbed_read):
		#    if "D" in read[4]:

		for ir, read in enumerate(subbed_read):
			if read[st_find.read_field] in rejected_set:
				print("revert", read)
				start = int(read[st_find.index_field]) - 1
				ref = ref[:start] + ori_ref[start:start + len(read[st_find.read_field])] + ref[start + len(read[st_find.read_field]):]
			# ref[start:start+len(read[3])] = ori_ref[start:start+len(read[3])]

			# print(pos,change,support_count,tmp_count)

		# remove rejected reads from  subbed_reads and file
		new_subbed_read = []
		for read in subbed_read:
			if read[st_find.read_field] not in rejected_set:
				new_subbed_read.append(read)
		# restore index pos for writing files
		write_sam(new_subbed_read, sub_read_file,restore_pos=True)
		print(len(subbed_read), "reads originally in ", os.path.basename(sub_read_file))
		subbed_read = new_subbed_read
		print(len(rejected_set), "reads reverted, removed from", os.path.basename(sub_read_file), "now ", len(subbed_read))
		subbed_read_set = set([x[st_find.read_field] for x in subbed_read])
		subprocess.run("rm -r batch_*",stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

		#new_misp_conflict,no_fill = curr_gap_reads(ref, strain, subbed_read_set, candidate_read,
		#						   raw_matrix_info.narrowed_matrix, read_freq,misp_conflict)
		#misp_conflict = new_misp_conflict
	# if reach_end:

	'''else:

		with open("misp_" + os.path.basename(sub_read_file), "w+") as mf:
			for misp_pos in sorted(list(accepted_misps.keys())):
				val = accepted_misps[misp_pos]
				bases = val[0].split("|")
				print(str(misp_pos+1), "&", bases[0], "&", bases[1], "&", val[1], "&", val[2], "\\\\\n\hline")
				mf.write(str(misp_pos+1) + " & " + str(bases[0]) + " & " + bases[1] + "&" + str(val[1]) + " & " + str(
					val[2]) + "\\\\\n\hline\n")'''

	# return 0
	return subbed_read,accepted_misps
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
def fac_verify_misp(ref_file, support_samfile, sub_read_file, new_strain_file, candidate_sam):
	strain = int(re.search('[0-9]+', os.path.basename(sub_read_file))[0])
	with open("rejected_" + os.path.basename(sub_read_file), "w+") as rejf:
		rejf.write("")
	initial_matrix_info = bm.matrix_from_readlist(bm.read_sam(support_samfile), 0, marked_id=set({}), target="raw")
	matrix_info = initial_matrix_info
	verified_read, accpeted_misps = verify_misp(ref_file, matrix_info, sub_read_file, new_strain_file, candidate_sam)

	tmp_rejected = []

	tmp_verified_read, tmp_accepted_misps = verify_misp(ref_file, matrix_info, sub_read_file, new_strain_file, candidate_sam)
	misp_conflict = {}
	verified_misp,verified_misp_source,verified_misp_read, tmp_reads = get_misp(get_ref_seq(ref_file), verified_read,
																				fix_pos=False)
	for k in verified_misp.keys():
		misp_conflict.update({k:[]})
	strain = int(re.search('[0-9]+', os.path.basename(sub_read_file))[0])
	read_freq = {}
	candidate_read = []
	with open(candidate_sam, "r") as canf:
		for line in canf:
			sline = line.strip().split(" ")
			candidate_read.append(sline)
			read_freq.update({sline[st_find.read_field]: int(sline[st_find.freq_field])})
	candidate_read = fix_s_pos(candidate_read)


	while len(verified_read) != len(tmp_verified_read):
		verified_read = tmp_verified_read
		accpeted_misps = tmp_accepted_misps

		tmp_subbed_read,misp_conflict = curr_gap_reads(get_ref_seq(new_strain_file),strain,verified_read,candidate_read,matrix,read_freq,misp_conflict)
		tmp_verified_read, tmp_accepted_misps = verify_misp(ref_file, matrix_info, sub_read_file, new_strain_file,
															candidate_sam)


	#write misp
	accepted_misps = tmp_accepted_misps
	with open("misp_" + os.path.basename(sub_read_file), "w+") as mf:
		for misp_pos in sorted(list(accepted_misps.keys())):
			val = accepted_misps[misp_pos]
			bases = val[0].split("|")
			print(str(misp_pos + 1), "&", bases[0], "&", bases[1], "&", val[1], "&", val[2], "\\\\\n\hline")
			mf.write(str(misp_pos + 1) + " & " + str(bases[0]) + " & " + bases[1] + "&" + str(val[1]) + " & " + str(
				val[2]) + "\\\\\n\hline\n")

	'''while (len(tmp_rejected) != len(rejected_reads)):
	
		tmp_rejected = copy.deepcopy(rejected_reads)

		rejected_reads,misp_conflict, no_candidate = verify_misp(ref_file, support_samfile, sub_read_file, new_strain_file, candidate_sam, misp_conflict)

		verify_sub_command = "for f in {56..64}; do\n" + \
							 os.path.dirname(__file__) + "/find_sub.sh" + " " + "-r" + " " + "final_strain_" + str(
			strain) + "_reference.fa" + " " + "-1" + " " + "half_real_R1.fastq" + " " + "-2" + " " + "half_real_R2.fastq" + \
							 " " + "-m" + " " + "tog" + " " + "-a" + " " + "bowtie2" + " " + "-c" + " " + 'True' + " -d Y;" + \
							 "cat " + os.path.dirname(
			__file__) + "/ mullti_support /${f}_out / round_1_extract.sam >> strain" + str(
			strain) + "_combine_extract.sam\ndone"

		verify_proc = subprocess.run(verify_sub_command,
									 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
									 shell=True)'''

#    print(rejected_reads)
# if rejected_reads is int:
#    break
def get_ref_seq(ref_file):
	ref  =""
	with open(ref_file,"r") as rf:
		for line in rf:
			if line[0] != ">":
				ref += line.strip()
	return ref


def curr_gap_reads(ref, strain, subbed_read_set, candidate_read,  matrix, read_freq,misp_conflict):
	"""
	Fill more candidate reads after rejecting reads
	:param ref: reference sequence
	:param strain: strain number
	:param subbed_read_set: set of remained subbed reads
	:param candidate_read:
	:param matrix: matrix of other SRR samples
	:param read_freq: frequency of each read
	:param misp_conflict: misp that causes gaps when placed together, format {pos:{pos_misp:[changed_base]}}
	:return:
	"""
	readfile1 = "half_real_R1.fastq"
	readfile2 = "half_real_R2.fastq"

	prev_readset = set({})
	covered_pos = {}
	for i in range(0, strain + 1):
		print("getting replacement reads in subbed_reads_" + str(i))
		with open("subbed_reads_" + str(i) + ".sam", "r") as subf:
			for line in subf:
				sline = line.strip().split(" ")
				prev_readset.add(sline[st_find.read_field])
		with open("rejected_subbed_reads_" + str(i) + ".sam", "r") as rejf:
			for line in rejf:
				fields = line.strip().split(" ")
				prev_readset.add(fields[st_find.read_field])
	with open("subbed_reads_" + str(strain) + ".sam", "r") as subf:
		for line in subf:
			fields = line.strip().split(" ")
			start = int(fields[st_find.index_field]) - 1
			for ib, base in enumerate(fields[st_find.read_field]):
				covered_pos.update({start + ib: base})

	print(len(prev_readset))

	readlist = candidate_read

	subbed_read = []
	#misPs = []
	for ir, read in enumerate(readlist):
		if read[st_find.read_field] in subbed_read_set:
			subbed_read.append(read)
			#misPs.extend(read[st_find.misp_field])
	original_readlist = copy.deepcopy(readlist)
	misPs, misP_source, misP_reads, new_readlist = get_misp(ref, readlist, False)
	readlist = new_readlist

	count = 0
	for batch, read in enumerate(readlist):

		count += 1
		if read[st_find.read_field] in prev_readset:
			# batch += 1
			continue
		read_index = int(read[st_find.index_field]) - 1
		#read_misPs = [list(x.items())[0][0] for x in read[st_find.misp_field]]

		print("curr read", batch, read)

		if if_clash(read, covered_pos) or if_reject(read,matrix,read_freq,misP_reads):
			continue

		temp_ref = ref[:read_index] + read[st_find.read_field] + ref[read_index + len(read[st_find.read_field]):]
		# print(editdistance.eval(temp_ref,ref))
		gap_list = verify_seq_support(temp_ref,batch,readfile1,readfile2)
		if len(gap_list) == 0:
			subbed_read.append(read)
			read_start = int(read[2]) - 1
			for ind, base in enumerate(read[st_find.read_field]):

				if read_start + ind in covered_pos.keys():

					if covered_pos[read_start + ind] != base:
						print("clash happened for current read", ind, covered_pos[read_start + ind], base)
						exit()
				else:
					covered_pos.update({read_start + ind: base})


			ref = temp_ref
			print(batch, "no gaps for ")
		else:
			#construct relational gaps
			print("gaps at ",gap_list)
			with open("gaps.txt", "a+") as gf:
				gf.write("gaps for read," + str(read) + " is " + str(gap_list) + "\n")
				for tmp_read in subbed_read:
					tmp_read_index= int(tmp_read[st_find.index_field]) -1
					if abs(tmp_read_index-read_index) <= 300:
						gf.write(str(tmp_read))
				gf.write("\n")


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
	return misp_conflict



def write_sam(readlist, filename, mode="w+",restore_pos=False):
	if restore_pos:
		readlist = fix_s_pos(readlist,restore=True)
	if filename != '':
		with open(filename, mode) as wf:

			for line in readlist:
				s = ""
				field_count = 0
				for field in line:
					if field_count == st_find.misp_field-1:
						s += str(field) + "\n"
					elif field_count>=st_find.misp_field:
						s += ""
					else:
						s += str(field) + " "
					field_count += 1
				wf.write(s)
			# wf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] +  "\n")
	else:
		for line in readlist:
			for field in line:
				if field == line[-1]:
					print(field)
				else:
					print(field, end=" ")

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
	subbed_read = bm.read_sam(sub_file,True)
	misP,misp_source,misp_reads,subbed_read = get_misp(ref, subbed_read, False)
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
					changed_base = change.split("|")[1]

					ref_codon = ref[pos-index:pos+(3-index)]
					if changed_base != "-":
						# need to consider multiple misp together
						mut_codon = ""
						tmp_read = []
						for read in subbed_read:
							if misp_reads[pos][changed_base][0][st_find.read_field] == read[st_find.read_field]:
								tmp_read = read
						assert len(tmp_read) != 0
						#print([list(x.items()) for x in tmp_read[st_find.misp_field]])
						tmp_misP_set = set([list(x.items())[0][0] for x in tmp_read[st_find.misp_field]])
						read_index = int(tmp_read[2])-1
						# mut_codon = ref_codon[0:index] + sorted_mutated_read[smr_index][read_field][mp-smr_read_index] + ref_codon[index + 1:]
						for read_misp in tmp_read[st_find.misp_field]:
							mp = list(read_misp.items())[0][0]

							if mp != pos:
								#print(mp,pos)
								continue
							for i in range(mp - index, mp + (3 - index)):
								#print(tmp_read[st_find.read_field][i - read_index])
								# if the adjacent spots are also misPs, will use the read's adjecent spot instead of ref
								if i in tmp_misP_set and 0 <= i - read_index < len(tmp_read[st_find.read_field]):
									mut_codon += tmp_read[st_find.read_field][i - read_index]
								else:
									mut_codon += ref[i]
					assert len(mut_codon) ==3
					#print(ref_codon,mut_codon,index, ref[pos-index],change,ref[pos+(3-index)-1], translate[ref_codon],translate[mut_codon])

					del_move = len([x for x in move if x < pos])

					if translate[ref_codon] == translate[mut_codon]:
						print("synonymous ",pos+1,change)
						syn_stat.update({pos:"synonymous"})
					else:
						print("non-synonymous",pos+1,change)

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
if __name__ == "__main__":
	"""
	fix numerical order of arguments
	fix verify_misp combining process
	"""

	#ref_file, candidate_sam, readfile1, readfile2
	#single_gap_reads(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
	#verify_sub_command = "rm -r batch_*"

	#verify_proc = subprocess.run(verify_sub_command,
	#							 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
	#							 shell=True)
	# ref_file, strain, samfile, readfile1, readfile2


	strain_max = int(sys.argv[5])
	for strain in range(0,strain_max+1):
		#identify_strain(sys.argv[1], strain, sys.argv[2], sys.argv[3], sys.argv[4])
		#subprocess.run("rm -r batch_*",stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,shell=True)
		#other_SRR_command = os.path.dirname(__file__) + "/get_combined_extract.sh final_strain_"+str(strain)+"_reference.fa"
		#other_SRR_proc = subprocess.run(other_SRR_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)


		backup_command = "cp subbed_reads_"+str(strain)+".sam subbed_reads_"+str(strain)+".sam.original; cp final_strain_"+str(strain) \
			+ "_reference.fa final_strain_"+str(strain)+"_reference.fa.original"
		#backup_proc = subprocess.run(backup_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,shell=True)

		fac_verify_misp(sys.argv[1], "combine_final_strain_" + str(strain) + "_reference.fa_extract.sam",
						"subbed_reads_" + str(strain) + ".sam", "final_strain_" + str(strain) + "_reference.fa",
						sys.argv[2])
		subprocess.run("rm -r batch_* *_out", stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,shell=True)