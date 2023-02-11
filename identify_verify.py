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
import verify

# the misp base in change_str starts from misp_str[2]
change_base_bit = 2
def fix_s_pos(subbed_reads,restore=False):
	new_subbed_reads = []
	for ir, fields in enumerate(subbed_reads):
		tmp_read = copy.deepcopy(fields)
		tmp_read[st_find.index_field] = int(tmp_read[st_find.index_field])


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
def restore_insertion_pos(subbed_read,inserted_pos):
	new_subbed_read = []
	sorted_inserted_poskey = sorted(inserted_pos.keys())
	index_offsets = [0] * len(subbed_read)
	curr_offset = 0
	#adjust iPos after each insertion restored
	for iPos in sorted_inserted_poskey:

		for ir,read in enumerate(subbed_read):
			read_index = int(read[st_find.index_field]) - 1

			if read_index >= iPos-curr_offset:

				index_offsets[ir] += inserted_pos[iPos]

				#assert read_index != iPos -curr_offset + inserted_pos[iPos]-1 and read_index != iPos-curr_offset
		curr_offset += inserted_pos[iPos]
	for ir,read in enumerate(subbed_read):
		new_read = copy.deepcopy(read)
		new_read_index = int(read[st_find.index_field]) - index_offsets[ir]
		new_read[st_find.index_field] = str(new_read_index)
		new_subbed_read.append(new_read)

	return new_subbed_read

def fix_insertion_pos(curr_read,subbed_read,candidate_read,covered_pos, inserted_pos):
	"""
	After an insertion read is replaced, fix subbed read index and candidate read index, also fix covered_pos position
	inserted_pos is also updated. Duplicate insertion will be detected with inserted_pos and no change will happen
	:param curr_read: read with insertion
	:param subbed_read: list of subbed read
	:param candidate_read: list of candidate read
	:param covered_pos: covered bases
	:param inserted_pos: posiiton of insertion and length dict
	:return: new_subbed_read, new_candidate_read,new_covered_pos,new_inserted_pos
	"""

	r_start = int(curr_read[st_find.index_field])-1
	cigar = curr_read[st_find.cigar_field]
	cigar_str = re.findall(r"[0-9]+[MIDSH]", cigar)
	blk_pos = []
	blk_type = []
	blk_length = []
	ini = 0
	fixed_ins_ini  = 0
	tmp_length = 0
	base_length = 0
	matched = 0
	currInsLength = {}


	for block in cigar_str:
		m = re.search(r'([0-9]+)([MIDSH])', block)
		bl = int(m.group(1)) + ini
		bLength = int(m.group(1))
		bt = str(m.group(2))

		blk_type.append(bt)

		blk_pos.append(bl)
		blk_length.append(int(m.group(1)))
		if bt == "I":
			fixed_ins_pos = r_start + ini - fixed_ins_ini
			# deduct insertion not recorded yet in inserted_pos and covered_pos
			if fixed_ins_pos not in inserted_pos.keys():

				currInsLength.update({r_start + ini:int(m.group(1))})
				fixed_ins_ini += bLength
		ini = bl
		# if iread[4]=="140M2D10M":
		#	print(block, int(m.group(1)))
		base_length += int(m.group(1))
		if bt != "I" and bt != "H":
			tmp_length += int(m.group(1))  # get read length without counting clipped bases

	new_inserted_pos = {}
	new_subbed_read = []
	for ir,read in enumerate(subbed_read):
		new_read = copy.deepcopy(read)
		for iPos in currInsLength.keys():
			if int(read[st_find.index_field])-1 >= iPos:

				new_read[st_find.index_field] = int(read[st_find.index_field])+ currInsLength[iPos]
		new_subbed_read.append(new_read)
	new_candidate_read = []
	for ir,read in enumerate(candidate_read):
		new_read = copy.deepcopy(read)
		for iPos in currInsLength.keys():
			if int(read[st_find.index_field])-1 >= iPos:
				#print(read, "add ",currInsLength[iPos],"at",iPos)
				new_read[st_find.index_field] = int(read[st_find.index_field])+ currInsLength[iPos]

		new_candidate_read.append(new_read)
	new_covered_pos = {}


	for pos in covered_pos.keys():
		posUpdLength = 0
		ins_length = 0
		for iPos in currInsLength.keys():
			if pos >= iPos - ins_length:
				posUpdLength += currInsLength[iPos]
			ins_length += currInsLength[iPos]
		new_covered_pos.update({pos + posUpdLength:covered_pos[pos]})

	for insPos in inserted_pos:
		posUpdLength = 0
		ins_length = 0
		for iPos in currInsLength.keys():
		# if insPos with fixed_ins_offset is larger than iPos then adjust pos
		# is fixed_ins_offset the insPOslength?
			if insPos > iPos - ins_length:
				posUpdLength += currInsLength[iPos]
			if insPos == iPos -ins_length:
				raise "duplicate insertion included? " +str(insPos) + str(iPos) + str(ins_length) + "\n"+str(inserted_pos)
			ins_length += currInsLength[iPos]
		new_inserted_pos.update({insPos+posUpdLength:inserted_pos[insPos]})
	new_inserted_pos.update(currInsLength)
	return new_subbed_read, new_candidate_read,new_covered_pos,new_inserted_pos

def if_clash(read,covered_pos,subbed_read, inserted_pos=None):
	"""
	determine if the read has conflict with existing substitution
	:param read:
	:param covered_pos: dict for bases at pos
	:param inserted_pos: positions for insertion, used to determine insertion reads, default None
	:param subbed_read:
	:return: bool
	"""
	r_start = int(read[st_find.index_field]) - 1
	if inserted_pos is None:
		for rei, base in enumerate(read[st_find.read_field]):

			if r_start + rei in covered_pos.keys():
				if read[st_find.read_field][rei] != covered_pos[r_start + rei]:

					print("clash at", r_start + rei, covered_pos[r_start + rei], read[st_find.read_field][rei])
					return True
		return False

	else:
		#insertion clash can happen either an insertion have different length/bases at the position of existed insertions
		# or the insertion is at a position in the middle of subbed reads with no insertion at such position
		# an important issue is the potential change of positions by the insertion of current read when doing comparison
		# if an insertion is in inserted_pos, include it in index computation, otherwise not

		cigar = read[st_find.cigar_field]
		cigar_insert_str = re.findall('[0-9]+[MIDSH]', cigar)

		clash = False
		ini = 0
		fixed_ins_ini = 0
		blk_type = []
		blk_pos = []
		currInsPosLength = {}
		ins_occur = 0
		for blk in cigar_insert_str:
			included_ins = False
			m = re.search(r'([0-9]+)([MIDSH])', blk)
			bl = int(m.group(1)) + ini
			bLength = int(m.group(1))
			bt = m.group(2)
			blk_type.append(bt)
			blk_pos.append(bl)
			if bt == "I":
				insStartPos = r_start + ini - fixed_ins_ini
				if insStartPos in inserted_pos.keys():
					included_ins = True
					if bLength != inserted_pos[insStartPos]:
						clash = True
						print("different length insertion clash at",insStartPos,r_start+ini, inserted_pos[r_start+ini],bLength)
						break
					else:
						for j in range(ini - fixed_ins_ini, bl-fixed_ins_ini):
							if r_start + j not in covered_pos:
								raise "Insertion not in covered_pos"
							if read[st_find.read_field][j+fixed_ins_ini] != covered_pos[r_start+j]:
								clash = True
								print("insertion base clash", r_start+j+fixed_ins_ini,fixed_ins_ini,read[st_find.read_field][j], covered_pos[r_start+j])
								break
						if clash:
							break

				else:
					#insertion not in inserted_pos, check whether other reads cover it

					if insStartPos in covered_pos:
						for test_read in subbed_read:
							if int(test_read[st_find.index_field])-1 < insStartPos < int(test_read[st_find.index_field])-1 + len(test_read[st_find.read_field])-1:
								print("insertion not existed in subbedread",test_read ,insStartPos,fixed_ins_ini)
								clash = True
								break
						if clash:
							break
						else:
							currInsPosLength.update({r_start + ini: bLength})
					else:
						currInsPosLength.update({r_start + ini : bLength})
					fixed_ins_ini += bLength

			ini = bl
		c = 0
		if not clash:
			insertion_offset = 0
			sortedInsPos = sorted(currInsPosLength.keys())
			print(currInsPosLength)
			index_sortedInsPos = 0

			for j in range(0, len(read[st_find.read_field])):  # change here to fill the blank with 0?

				if len(currInsPosLength) > 0:
					if r_start + j == sortedInsPos[index_sortedInsPos]:
						insertion_offset += currInsPosLength[sortedInsPos[index_sortedInsPos]]
						if index_sortedInsPos < len(sortedInsPos) - 1:

							index_sortedInsPos += 1


				if blk_type[c] == "M" or blk_type[c] == "S" or blk_type[c] == "D":
					if j + r_start - insertion_offset in covered_pos:
						if read[st_find.read_field][j  ] != covered_pos[j+r_start-insertion_offset]:
							clash = True
							print("base clash at read ",j + insertion_offset,read[st_find.read_field][j - insertion_offset],insertion_offset,covered_pos[j+r_start] )
							break
		
		return clash

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
		batch) + "_reference.fa" + " " + "-1" + " " + readfile1 + " " + "-2" + " " + readfile2 + " " + "-m" + " " + "tog" + " " + "-c" + " " + "True" + " -d Y"

	verify_proc = subprocess.run(verify_sub_command,
								 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
								 shell=True)
	# print(verify_proc.stdout.split("\n")[-5:])

	gap_list = verify_proc.stdout.split("\n")[-3]
	read_count_output = []
	for line in verify_proc.stdout.split("\n"):
		m = re.match(".*reads in.*extract.sam",line)
		if m is not None:
			read_count_output.append(m.group(0))
	if gap_list != "no gaps":


		gap_list = re.sub('[\[\]]', '', gap_list)
		gap_list= gap_list.split(",")[:-1]
		print(gap_list)
		for line in read_count_output:
			print(line)

		return gap_list
	# gap_list = gap_list.split(",")[:-1]
	# ref = ref[:read_index] + ref_seq[read_index:read_index+len(read[3])] + ref[read_index + len(read[3]):]
	# subbed_read.remove(read)
	else:
		return []

def write_result(strain, ref_seq, subbed_read, stage=0, restore_pos=False, curr=False):
	write_seq(ref_seq,"final_strain_" + str(strain) + "_reference.fa")
	#with open("final_strain_" + str(strain) + "_reference.fa", "w+") as bf:
	#	bf.write(">final_strain_" + str(strain) + "\n")
	#	bf.write(ref_seq)
	#with open("strain_" + str(strain) + "_spike.fa", "w+") as sf:
	#	sf.write(">strain_" + str(strain) + "_spike" + "\n")
	#	sf.write(ref_seq[21562:25384])
	if stage > 0:
		if curr:
			write_sam(subbed_read, "curr_stage_"+str(stage)+"_subbed_reads_" + str(strain) + ".sam",restore_pos=restore_pos)
			write_seq(ref_seq,"curr_stage_"+str(stage)+"_final_strain_" + str(strain) + "_reference.fa")

	write_sam(subbed_read,"subbed_reads_" + str(strain) + ".sam",restore_pos=restore_pos)
	subprocess.run("rm -r batch_*",shell=True)
	print(subbed_read)



def identify_strain(ref_file, strain, candidate_sam, readfile1, readfile2):
	"""produce and verify new strain, also replace rejected reads
	@param ref_file fasta file of reference sequence
	@param strain integer number of strain
	This function identifies new strains and verify , it iteratively goes through the list
	of candidate reads, check if substitution will result in sequence not supported by sequencing reads
	or cause stop codon that separates a protein from the middle
	"""

	prev_readset = set({})
	if strain > 0:
		for i in range(0, strain):
			print("getting reads in subbed_reads_" + str(i))
			with open("subbed_reads_" + str(i) + ".sam", "r") as subf:
				for line in subf:
					sline = line.strip().split(" ")
					prev_readset.add(sline[st_find.read_field])
			if os.path.exists("rejected_subbed_reads_" + str(i) + ".sam"):
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
	insertedPos = {}
	# remove reads that cause early stop condon


	i = 0
	while i < len(readlist):
		read = readlist[i]
		i += 1
		overlap = False
		read_index = int(read[st_find.index_field]) - 1
		read_misPs = [list(x.items())[0][0] for x in read[st_find.misp_field]]
		cigar = read[st_find.cigar_field]
		print("curr read", batch, read, "bases covered", len(covered_pos))
		if len(subbed_read) > 0:
			# read clash detection, find all overlapping reads
			# mp-clash, compare all MPs, if the reads only itroduce new MP outise current mp ranges
			if "I" in cigar:
				clash = if_clash(read, covered_pos, subbed_read, insertedPos)
			else:
				clash = if_clash(read, covered_pos, subbed_read)

			if clash:
				batch += 1
				continue
		if "I" in read[st_find.cigar_field]:
			insertLength = 0

			cigar_insert_str = re.findall('[0-9]+[MIDSH]',cigar)
			ini = 0
			fixed_ins_ini = 0
			for blk in cigar_insert_str:
				m = re.search(r'([0-9]+)([MIDSH])', blk)
				bLength = int(m.group(1))
				bl = ini + bLength
				bt = m.group(2)
				if bt == "I":
					if read_index + ini - fixed_ins_ini not in insertedPos:
						insertLength += int(m.group(1))
						#currinsPosLength.update({read_index + ini: int(m.group(1))})

						fixed_ins_ini += bLength

				ini = bl
			temp_ref = ref[:read_index] + read[st_find.read_field] + ref[read_index + len(read[st_find.read_field])-insertLength:]
		else:
			temp_ref = ref[:read_index] + read[st_find.read_field] + ref[read_index + len(read[st_find.read_field]):]
		gap_list = verify_seq_support(temp_ref,batch,readfile1,readfile2)
		if len(gap_list) == 0:
			# update covered positions with base in covered_pos
			if "I" in cigar:
				changed_subbed_read, readlist, covered_pos,insertedPos = fix_insertion_pos(read,subbed_read,readlist,covered_pos,insertedPos)


			subbed_read.append(copy.deepcopy(original_readlist[batch]))
			read_start = int(read[st_find.index_field]) - 1
			for ind, base in enumerate(read[st_find.read_field]):

				if read_start + ind in covered_pos.keys():

					if covered_pos[read_start + ind] != base:
						print("clash happened for current read", ind, covered_pos[read_start + ind], base)
						raise "if_clash() not working properly"
				else:

					covered_pos.update({read_start + ind: base})

			ref = temp_ref
			print(batch, "no gaps for curr read\n")
			#temp_subbed_read = restore_insertion_pos(subbed_read,insertedPos)
			write_sam(subbed_read,"batch_"+str(batch)+"_subbed_read.sam","w+",True)

		else:

			for gap_pos in gap_list:

				if abs(read_index-int(gap_pos)) > 300:
					print("gap 300 bases away")
					print(read, gap_pos)
					exit()

		batch += 1

	#restored_subbed_read = restore_insertion_pos(subbed_read,insertedPos)
	write_result(strain, ref,subbed_read)
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
		insertion_len = [0]
		insertion_total_len = 0

		for block in cigar_str:
			m = re.search(r'([0-9]+)([MIDSH])', block)

			bl = int(m.group(1)) + ini
			bLength  = int(m.group(1))
			bt = str(m.group(2))
			# if bt == "S" or bt == "H":
			#    continue

			blk_type.append(bt)
			blk_length.append(int(m.group(1)))
			blk_pos.append(bl)

			ini = bl
			if bt == "I":


				insertion_len.append(insertion_len[-1]+bLength)
				insertion_total_len += bLength
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

		insertion_count = 0

		for j in range(0, blk_pos[-1]):  # change here to fill the blank with 0?
			# compare base directly, soft clipping at beginning has been accounted

			if blk_type[c] == "M" or blk_type[c] == "S":
				#if len(insertion_len) > 0:

				if ref[read_index + j-insertion_len[insertion_count]]!=read_str[j]:
					temp_misp.append({read_index + j-insertion_len[insertion_count]:[ref[read_index + j-insertion_len[insertion_count]]+"|"+read_str[j]]})
				#else:
					#if ref[read_index + j]!=read_str[j]:
						#temp_misp.append({read_index + j:[ref[read_index + j]+"|"+read_str[j]]})
			# handle insertion
			elif blk_type[c] == "I":
				inserts.append(st_find.dict_tonu[read_str[j]])

			# should have been inserted "-" in missing position
			elif blk_type[c] == "D" or blk_type[c] == "H":
				#if len(insertion_len) > 0:
				temp_misp.append({read_index + j-insertion_len[insertion_count]: [ref[read_index  + j-insertion_len[insertion_count]] + "|" + read_str[ j]]})
				#else:
				#	temp_misp.append({read_index + j: [ref[read_index + j] + "|" + read_str[j]]})
			if j == blk_pos[c] - 1:  # update start and c, put inserts into hashtable
				#print(read[0])

				if blk_type[c] == "I":
					insert_str = ""
					for ins_base in inserts:
						insert_str += st_find.dict_tole[ins_base]

					temp_misp.append({read_index + j-insertion_len[insertion_count]:["-|"+insert_str]})
					if read_num in insertion_reads.keys():

						newinsert = copy.deepcopy(insertion_reads.get(included_read_num))
						newinsert.append((read_index+ begin, inserts))
						insertion_reads.update({included_read_num: newinsert})
					else:
						insertion_reads.update({included_read_num: [(read_index + begin, copy.deepcopy(inserts))]})
					if insertion_count < len(insertion_len) - 1:
						insertion_count += 1

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


def readlist_comapre(readlist1,readlist2):
	"""

	:param readlist1:
	:param readlist2:
	:return: count of total unique reads
	"""

	print(len(readlist1), len(readlist2))
	r1_set = set([x[3] for x in readlist1])
	#r1_set = set([x[3].replace("-","") for x in readlist1])
	r2_set = set([x[3] for x in readlist2])
	#r2_set = set([x[3].replace("-","") for x in readlist2])
	print("readlist1 unique")
	count = 0
	for read in readlist1:
		if read[3] not in r2_set:
			print(read)
			count +=1
	print("readlist2 unique")
	for read in readlist2:
		if read[3] not in r1_set:
			print(read)
			count += 1
	return count

def get_ref_seq(ref_file):
	ref  =""
	with open(ref_file,"r") as rf:
		for line in rf:
			if line[0] != ">":
				ref += line.strip()
	return ref



def write_misp(accepted_misps,filename):
	with open(filename, "w+") as mf:
		for misp_pos in sorted(list(accepted_misps.keys())):
			val = accepted_misps[misp_pos]
			bases = val[0].split("|")
			print(str(misp_pos + 1), "&", bases[0], "&", bases[1], "&", val[1], "&", val[2], "\\\\\n\hline")
			mf.write(str(misp_pos + 1) + " & " + str(bases[0]) + " & " + bases[1] + "&" + str(val[1]) + " & " + str(
				val[2]) + "\\\\\n\hline\n")
def write_seq(ref,filename):
	match = re.match('strain_([0-9]+)',filename)
	strain = ""
	if  match is not None:
		strain = match.group(1)
	with open(filename,"w+") as wf:
		wf.write(">strain_"+strain+"_reference\n")
		wf.write(ref)
	with open("spike_"+filename,"w+") as swf:
		swf.write(">strain_"+strain+"_spike\n")
		swf.write(ref[21562:25384])
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
	if len(sys.argv) > 6:
		strain_min = int(sys.argv[6])
	else:
		strain_min = 0

	for strain in range(strain_min,strain_max+1):
		identify_strain(sys.argv[1], strain, sys.argv[2], sys.argv[3], sys.argv[4])

		subprocess.run("mv batch_* strain_"+str(strain)+"intermit_out",stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,shell=True)

		#verify.verify_strain(strain,sys.argv[1],"subbed_reads_"+str(strain)+".sam")