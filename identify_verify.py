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
def fix_s_pos(subbed_reads):
	new_subbed_reads = []
	for ir, fields in enumerate(subbed_reads):
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
			if index - blk_pos[0] < 0:

				start_pos = blk_pos[0] - index
				tmp_length -= start_pos
				index = 0
			else:
				index = index - blk_pos[0]
			subbed_reads[ir][st_find.index_field] = index + 1
	return subbed_reads


def invalid_stop(read_list):
	"""
	use misp to derive stop codon, not handle indel
	:param read_list: the list of candidate reads
	:return: list of reads after removing problematic ones
	"""
	valid_reaads = []

	return valid_reaads


def identify_strain(ref_file, strain, samfile, readfile1, readfile2):
	"""produce and verify new strain
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

	#getting misp
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
	# remove reads that cause early stop condon


	for batch, read in enumerate(readlist):
		overlap = False
		read_index = int(read[st_find.index_field]) - 1
		read_misPs = [int(x[0]) for x in read[st_find.misp_field:]]

		print("curr read", batch, read, "bases covered", len(covered_pos))
		if len(subbed_read) > 0:
			# read clash detection, find all overlapping reads

			# mp-clash, compare all MPs, if the reads only itroduce new MP outise current mp ranges
			clash = False
			for rei, base in enumerate(read[st_find.read_field]):
				r_start = int(read[st_find.index_field]) - 1
				if r_start + rei in covered_pos.keys():
					if read[st_find.read_field][rei] != covered_pos[r_start + rei]:
						clash = True
						print("clash at", r_start + rei, covered_pos[r_start + rei], read[st_find.read_field][rei])
						break

			if clash:
				batch += 1
				continue

		temp_ref = ref[:read_index] + read[st_find.read_field] + ref[read_index + len(read[st_find.read_field]):]

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

			gap_list = re.sub('[\[\]]', '', gap_list)
			print(gap_list.split(",")[:-1])

		else:

			subbed_read.append(read)
			# update covered positions with base in covered_pos
			read_start = int(read[st_find.index_field]) - 1
			for ind, base in enumerate(read[st_find.read_field]):

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
	print(misp1)
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


def get_misp(ref_file, sub_read_file_list, printing=True):
	"""

	:param ref_file: reference sequence file
	:param sub_read_file_list: sam file list
	:param printing: True to print to stdout
	:return: misPs {position: ["ref_base|misp_base"]}, misP_source {position:{change_str:[SRR_ID]}}, misP_reads {position:{change_str:[reads]}}
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

	for sub_read_file in sub_read_file_list:
		with open(sub_read_file, "r") as mf:
			insertion_reads = {}
			included_read_num = 0
			for read_num,read_str in enumerate(mf):
				read = read_str.strip().split(" ")
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

				#correct read_index for soft-clipping at beginning
				if blk_type[0] == "S":
					if read_index - blk_pos[0] < 0:
						start_pos = blk_pos[0] - read_index
						tmp_length -= start_pos
					else:
						read_index = read_index - blk_pos[0]

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


				# record misp in this blk
				for misp_kv in temp_misp:
					misp = [(k,v) for k,v in misp_kv.items()][0]
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
							misP_reads.update({misp[0]: {misp[1][0].split("|")[1]: [read]}})

					if misp[0] not in pos.keys():

						pos.update({misp[0]: [read]})
					else:
						pos[misp[0]].append(read)



	# with open("different_bases.txt","w+") as wf:

	for k, v in sorted(pos.items(), key=lambda x: x[0]):

		# if 21562<=k<=25384:
		if k not in misPs.keys():
			print(k,misPs.keys())
		for change in misPs[k]:
			support_list = [(k, v) for k, v in collections.Counter(misPs_source[k][change[-1]]).items() if v > 1]
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


def verify_misp(ref_file, samfile, sub_read_file_list, new_strain_file, candidate_sam):
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

	strain = int(re.search('[0-9]+', os.path.basename(sub_read_file_list))[0])

	test_misp = {}
	subbed_read = []
	with open(sub_read_file_list, "r") as subf:
		for line in subf:
			fields = line.strip().split(" ")
			subbed_read.append(fields)

	misPs, misPs_source, misP_reads = get_misp(ref_file, [sub_read_file_list], False)
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

			support_index = np.where(tmp == st_find.dict_tonu[change.split("|")[1]])
			support_count = support_index[0].shape[0]

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
						if ref[curr - 1] == temp_read[st_find.read_field][curr - 1 - read_index]:
							l_window += 1

						else:
							l_reached = True
							break
						curr -= 1
					curr_l_line = sp - 1
					if curr_l_line >= 0 and curr - 1 <= raw_matrix_info.narrowed_read[curr_l_line][st_find.index_field] - 1 + len(
							raw_matrix_info.narrowed_read[curr_l_line][st_find.read_field]):
						curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
					# print(curr_l_read)

					else:
						l_reached = True
					while not l_reached and curr_l_line >= 0:
						if curr - 1 < int(curr_l_read[st_find.index_field]) - 1:

							if curr_l_line - 1 < 0:
								l_reached = True
								break

							curr_l_line -= 1
							curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
							if curr - 1 > int(curr_l_read[st_find.index_field]) - 1 + len(curr_l_read[st_find.read_field]):
								l_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_l_read[st_find.index_field]) - 1 <= pos <= int(curr_l_read[st_find.index_field]) - 1 + len(curr_l_read[st_find.read_field]) - 1 and \
									curr_l_read[st_find.read_field][pos - (int(curr_l_read[st_find.index_field]) - 1)] != change.split("|")[1]:
								# print("read not supporting misp",pos,change)
								# print(curr_l_read)
								# l_reached= True
								curr_l_line -= 1

								continue
							# print(curr_l_read)
							continue

						# try:
						if ref[curr - 1] == curr_l_read[st_find.read_field][curr - 1 - (int(curr_l_read[st_find.index_field]) - 1)]:
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
					while not r_reached and curr_r_line <= raw_matrix_info.narrowed_matrix.shape[0]:
						if curr + 1 > int(curr_r_read[st_find.index_field]) - 1 + len(curr_r_read[st_find.read_field]) - 1:

							if curr_r_line + 1 >= raw_matrix_info.narrowed_matrix.shape[0]:
								r_reached = True
								break
							curr_r_line += 1
							curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
							# if current read does not cover the curr position (gaps), break
							if curr + 1 < int(curr_r_read[st_find.index_field]) - 1:
								r_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_r_read[st_find.index_field]) - 1 <= pos <= int(curr_r_read[st_find.index_field]) - 1 + len(curr_r_read[st_find.read_field]) and \
									curr_r_read[st_find.read_field][pos - (int(curr_r_read[st_find.index_field]) - 1)] != change[-1]:
								# print("read not supporting misp",pos,change, curr_r_read)
								# print(curr_r_read)
								curr_r_line += 1
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
						print("target misp", pos, read_freq.get(misP_reads[pos][change[-1]][0][st_find.read_field]))

						print("seg_length", len(seg))

					# print("ref:",ori_ref[pos-l_window:pos+r_window+1])

					# print("seg:",seg)
					# reject this read candidate
					if len(seg) < 50 and len(misP_reads[pos][change[-1]]) < 2 and read_freq.get(
							misP_reads[pos][change[-1]][0][st_find.read_field]) == 1:
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
						accepted_misps.update({pos: (change, read_freq.get(misP_reads[pos][change[-1]][0][st_find.read_field]), window)})

			else:
				if len(misP_reads[pos][change[-1]]) < 2 and read_freq.get(
						misP_reads[pos][change[-1]][0][st_find.read_field]) == 1:
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
							{pos: (change, read_freq.get(misP_reads[pos][change[-1]][0][st_find.read_field]), 0)})

		if reject:
			continue

	if len(rejected_read) > 0:
		print(len(rejected_read), "reads rejected")
		write_sam(rejected_read, "rejected_" + os.path.basename(sub_read_file_list), "a")
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
		write_sam(new_subbed_read, sub_read_file_list)
		print(len(subbed_read), "reads originally in ", os.path.basename(sub_read_file_list))
		subbed_read = new_subbed_read
		print(len(rejected_set), "reads reverted, removed from", os.path.basename(sub_read_file_list), "now ", len(subbed_read))
		subbed_read_set = set([x[st_find.read_field] for x in subbed_read])
		reach_end = curr_gap_reads(ref, strain, subbed_read_set, candidate_read, ori_ref,
								   raw_matrix_info.narrowed_matrix, misP_reads, read_freq)
	# if reach_end:

	else:

		with open("misp_" + os.path.basename(sub_read_file_list), "w+") as mf:
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

	for ir, read in enumerate(readlist):
		tmp_misP = []
		read_index = int(read[st_find.index_field]) - 1
		for i, base in enumerate(read[st_find.read_field]):
			if old_ref[read_index + i] != base:
				tmp_misP.append(read_index + i)
		if len(tmp_misP) == 0:
			print(read)
			print(ref[read_index:read_index + len(read[st_find.read_field])])
			print(read[st_find.read_field])
			exit()
		for mp in tmp_misP:
			read.append((int(mp), ref[mp], read[st_find.read_field][mp - read_index]))

	subbed_read = []
	misPs = []
	for ir, read in enumerate(readlist):
		if read[st_find.read_field] in subbed_read_set:
			subbed_read.append(read)
			misPs.extend(read[misp_bit:])

	count = 0
	for batch, read in enumerate(readlist):
		reject = False
		count += 1
		if read[st_find.read_field] in prev_readset:
			# batch += 1
			continue
		read_index = int(read[st_find.index_field]) - 1
		read_misPs = [int(x[0]) for x in read[misp_bit:]]

		print("curr read", batch, read, "misp range", misPs[0], misPs[-1])
		clash = False
		for rei, base in enumerate(read[st_find.read_field]):
			r_start = int(read[st_find.index_field]) - 1
			if r_start + rei in covered_pos.keys():
				if read[st_find.read_field][rei] != covered_pos[r_start + rei]:
					clash = True
					print("clash at", r_start + rei, covered_pos[r_start + rei], read[st_find.read_field][rei])
					break

		if clash:
			batch += 1
			continue

		for rm in read[misp_bit:]:
			pos = rm[0]
			change = rm[2]
			tmp = np.squeeze(matrix.getcol(pos).toarray())
			tmp_count = np.bincount(tmp)[1:]

			support_index = np.where(tmp == st_find.dict_tonu[change[-1]])
			support_count = support_index[0].shape[0]
			if support_count == 0 and read_freq.get(read[st_find.read_field]) == 1 and pos not in misP_reads.keys():
				reject = True
				break
			# print("rejected ",misP_reads[pos][change[-1]])

		if reject:
			continue
		temp_ref = ref[:read_index] + read[st_find.read_field] + ref[read_index + len(read[st_find.read_field]):]
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
			for ind, base in enumerate(read[st_find.read_field]):

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

def single_gap_reads(ref_file, candidate_sam, readfile1, readfile2):
	"""
	Find each read that substituted independently with no support
	:param ref:
	:param strain:
	:param samfile:
	:param readfile1:
	:param readfile2:
	:return:
	"""
	misp_bit = 6
	# strain = 2

	readlist = []
	with open(candidate_sam, "r") as samf:
		for line in samf:

			fields = line.strip().split(" ")
			if "N" in fields[3]:
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

	with open("single_valid_read.sam", "w+") as bf:
		for line in subbed_read:
			bf.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + " " + str(
				line[5]) + "\n")


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


if __name__ == "__main__":
	"""
	fix numerical order of arguments
	fix verify_misp combining process
	"""

	#ref_file, candidate_sam, readfile1, readfile2
	single_gap_reads(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
	verify_sub_command = "rm -r batch_*"

	verify_proc = subprocess.run(verify_sub_command,
								 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
								 shell=True)
	# ref_file, strain, samfile, readfile1, readfile2
	identify_strain(sys.argv[1],sys.argv[5],"single_valid_read.sam",sys.argv[3],sys.argv[4])
	#ref_file, samfile, mreads, new_strain_file, candidate_sam
