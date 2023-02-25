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
import identify_verify
import strain_finder as st_find
import identify_verify as idv


def new_verify_misp(ref_file, support_matrix_info, sub_read_file, new_strain_file, candidate_sam,stage=0):
	"""

	:param ref_file: original ref fa
	:param support_matrix_info: other SRR sam
	:param sub_read_file: subbed_read
	:param new_strain_file: generated strain fa
	:param candidate_sam: sub_read_candidates
	:return: rejected reads, updated misp_conflict, no_candidate to fill
	if multiple reads contain same misp,the misp will be accpeted.
	 if one read contains rejected misp, the read will be rejected, removed from misp_reads
	 misp with no read support will also be removed from misPs, misp_reads and misp_source

	"""
	ref = ""
	with open(new_strain_file, "r") as f:
		for line in f:
			if line[0] != ">":
				temp = line.replace(" ", "").replace("\n", "").replace("\d", "").upper()
				ref += re.sub('\d', "", temp)
	unvarified_ref = ref
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
	original_subbed_read = copy.deepcopy(subbed_read)
	misPs, misPs_source, misP_reads,new_subbed_read = idv.get_misp(ori_ref, subbed_read)

	subbed_read = new_subbed_read


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
	candidate_read = idv.fix_s_pos(candidate_read)

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
	rejected_misps = {}

def verify_misp(ref_file, support_matrix_info, sub_read_file, new_strain_file, candidate_sam,stage=0):
	"""

	:param ref_file: original ref fa
	:param support_matrix_info: other SRR sam
	:param sub_read_file: subbed_read
	:param new_strain_file: generated strain fa
	:param candidate_sam: sub_read_candidates
	:return: rejected reads, updated misp_conflict, no_candidate to fill
	if multiple reads contain same misp,the misp will be accpeted.
	 if one read contains rejected misp, the read will be rejected, removed from misp_reads
	 misp with no read support will also be removed from misPs, misp_reads and misp_source

	"""
	ref = ""
	with open(new_strain_file, "r") as f:
		for line in f:
			if line[0] != ">":
				temp = line.replace(" ", "").replace("\n", "").replace("\d", "").upper()
				ref += re.sub('\d', "", temp)
	unvarified_ref = ref
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
	original_subbed_read = copy.deepcopy(subbed_read)
	misPs, misPs_source, misP_reads,new_subbed_read = idv.get_misp(ori_ref, subbed_read)


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
	candidate_read = idv.fix_s_pos(candidate_read)

	inserted_pos = {}
	i = 0
	while (i < len(new_subbed_read)):
		temp_read = new_subbed_read[i]
		new_subbed_read,candidate_read,covered_pos,inserted_pos = identify_verify.fix_insertion_pos(new_subbed_read[i],new_subbed_read,candidate_read,{},inserted_pos)
		i += 1
	misp_fixed_subbed_read = []
	fix_to_ori_mPos = {}
	i = 0
	while (i < len(new_subbed_read)):
		temp_read = new_subbed_read[i]
		misp_fixed_read = copy.deepcopy(temp_read[:st_find.misp_field])
		misp_fixed_read.append([])
		for misp in temp_read[st_find.misp_field]:
			pos = list(misp.keys())[0]
			ori_pos = pos
			change = misp[pos]
			for insPos in inserted_pos.keys():
				if pos >= insPos:
					pos += inserted_pos[insPos]
			misp_fixed_read[st_find.misp_field].append({pos:change})
			fix_to_ori_mPos.update({pos:ori_pos})
		print(misp_fixed_read)
		i += 1
		misp_fixed_subbed_read.append(misp_fixed_read)
	subbed_read = misp_fixed_subbed_read
	# fix insertion pos of subbed read


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
	rejected_misps = {}
	for ir,read in enumerate(subbed_read):
		#if read_freq[read[st_find.read_field]] > 1:
		#	rejected = False
		#else:
		#	rejected = True
		rejected = False
		tmp_accepted_misps = {}
		for misp in read[st_find.misp_field]:
			pos = list(misp.keys())[0]
			change = misp[pos][0]
			if pos in accepted_misps.keys() :
				accepted_misps.update({pos:(accepted_misps[pos][0],accepted_misps[pos][1]+read_freq[read[st_find.read_field]],accepted_misps[pos][2])})
				print(pos,accepted_misps[pos],"already accpted, add frequency, continue")
				continue
			if pos in rejected_misps.keys():
				if read_freq[read[st_find.read_field]] == 1 and change in set(rejected_misps[pos]):
					rejected_read.append(read)
					print(read,pos,change,"already rejected, continue",rejected_misps[pos])
					continue
			if pos > matrix.shape[1]:
				print( misp,"out of range")
				print(matrix.shape)
				exit()
			tmp = np.squeeze(matrix.getcol(pos).toarray())
			tmp_count = np.bincount(tmp)[1:]

			changed_base = change.split("|")[1]
			support_index = np.where(tmp == st_find.dict_tonu[changed_base])
			support_count = support_index[0].shape[0]
			visited_sindex = set({})
			if support_count > 0:
				for sp in support_index[0]:
					temp_read = raw_matrix_info.narrowed_read[sp]
					read_index = int(temp_read[st_find.index_field]) - 1
					curr = pos
					read_end = read_index + len(temp_read[st_find.read_field])
					reached = False
					l_reached = False
					l_window = 0


					r_reached = False
					r_window = 0


					while curr - 1 >= read_index:
						try:
							if ref[curr - 1] == temp_read[st_find.read_field][curr - 1 - read_index]:
								l_window += 1

							else:
								l_reached = True
								break
						except:
							print(curr, temp_read, len(ref))
							raise "index error"
						curr -= 1
					curr_l_line = sp - 1
					if curr_l_line >= 0 and curr - 1 <= raw_matrix_info.narrowed_read[curr_l_line][
						st_find.index_field] - 1 + len(
							raw_matrix_info.narrowed_read[curr_l_line][st_find.read_field]):
						curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
					# print(curr_l_read)

					else:
						l_reached = True
					# while not l_reached and curr_l_line >= 0:
					for curr_l_line in range(sp - 1, 0, -1):
						if l_reached:
							break
						if curr_l_line in visited_sindex:
							continue
						else:
							visited_sindex.add(curr_l_line)
						curr_l_read = raw_matrix_info.narrowed_read[curr_l_line]
						if curr - 1 < int(curr_l_read[st_find.index_field]) - 1 or curr - 1 >= int(
								curr_l_read[st_find.index_field]) - 1 + len(curr_l_read[st_find.read_field]):

							if curr_l_line - 1 < 0:
								l_reached = True
								break

							if curr - 1 > int(curr_l_read[st_find.index_field]) - 1 + len(
									curr_l_read[st_find.read_field]):
								l_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_l_read[st_find.index_field]) - 1 <= pos <= int(
									curr_l_read[st_find.index_field]) - 1 + len(curr_l_read[st_find.read_field]) - 1 and \
									curr_l_read[st_find.read_field][
										pos - (int(curr_l_read[st_find.index_field]) - 1)] != change.split("|")[1]:

								continue

							continue

						try:
							if ref[curr - 1] == curr_l_read[st_find.read_field][
								curr - 1 - (int(curr_l_read[st_find.index_field]) - 1)]:
								l_window += 1

							else:
								l_reached = True
								break
						except:
							print(curr - 1, curr - 1 - (int(curr_l_read[st_find.index_field]) - 1), curr_l_line,
								  curr_l_read, raw_matrix_info.narrowed_read[curr_l_line])
							raise "index errror"
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
					if curr_r_line < raw_matrix_info.narrowed_matrix.shape[0] and \
							raw_matrix_info.narrowed_read[curr_r_line][st_find.index_field] - 1 <= curr + 1:

						curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
					else:
						r_reached = True
					# while not r_reached and curr_r_line <= raw_matrix_info.narrowed_matrix.shape[0]:
					for curr_r_line in range(sp + 1, raw_matrix_info.narrowed_matrix.shape[0]):
						if r_reached:
							break
						if curr_r_line in visited_sindex:
							continue
						else:
							visited_sindex.add(curr_r_line)
						curr_r_read = raw_matrix_info.narrowed_read[curr_r_line]
						if curr + 1 > int(curr_r_read[st_find.index_field]) - 1 + len(
								curr_r_read[st_find.read_field]) - 1 or \
								curr + 1 < int(curr_r_read[st_find.index_field]) - 1:

							if curr_r_line + 1 >= raw_matrix_info.narrowed_matrix.shape[0]:
								r_reached = True
								break

							if curr + 1 < int(curr_r_read[st_find.index_field]) - 1:
								r_reached = True
								break
							# if current read does not support the current misp, break
							if int(curr_r_read[st_find.index_field]) - 1 <= pos <= int(
									curr_r_read[st_find.index_field]) - 1 + len(curr_r_read[st_find.read_field]) - 1 and \
									curr_r_read[st_find.read_field][
										pos - (int(curr_r_read[st_find.index_field]) - 1)] != changed_base:

								continue

							# print(curr_r_read)
							continue
						if ref[curr + 1] == curr_r_read[st_find.read_field][
							curr + 1 - (int(curr_r_read[st_find.index_field]) - 1)]:
							r_window += 1
						else:
							r_reached = True
							break
						curr += 1
					seg = l_seg + ref[pos:curr + 1]
					window = l_window + r_window + 1
					if seg not in found_segments:
						found_segments.add(seg)
					print("target misp", pos, read_freq.get(read[st_find.read_field]))

					print("seg_length", len(seg))

					'''if len(seg)<50 and read_freq.get(read[st_find.read_field]) == 1 :#and len(misP_reads[pos][changed_base])==1:
						rejected_read.append(read)
						#if misP_reads[pos][changed_base][0][st_find.read_field] == read[st_find.read_field]:
						#	misP_reads[pos].pop(changed_base)
						#	misPs[pos].remove(change)
						rejected = True
						if pos not in rejected_misps.keys():
							rejected_misps.update({pos:[change]})
						else:
							rejected_misps.get(pos).append(change)
						break'''
					if window < 50 and read_freq[read[st_find.read_field]]==1:
						rejected = True
						break
				#if not rejected:
				if rejected:  # and len(misP_reads[pos][changed_base])==1:

					rejected_read.append(read)
					# if misP_reads[pos][changed_base][0][st_find.read_field] == read[st_find.read_field]:
					#	misP_reads[pos].pop(changed_base)
					#	misPs[pos].remove(change)
					rejected = True
					if pos not in rejected_misps.keys():
						rejected_misps.update({pos: [change]})
					else:
						rejected_misps.get(pos).append(change)
					break
				else:
					tmp_accepted_misps.update({pos: (change, read_freq[read[st_find.read_field]], window)})

			else:
				if read_freq.get(read[st_find.read_field]) == 1:# and len(misP_reads[pos][changed_base]) < 2:
					#rejected_read.append(read)
					rejected_read.append(original_subbed_read[ir])
					rejected = True
					#if misP_reads[pos][changed_base][0][st_find.read_field] == read[st_find.read_field]:
					#	misP_reads[pos].pop(changed_base)
					#	misPs[pos].remove(change)
					if pos not in rejected_misps.keys():
						rejected_misps.update({pos: [change]})
					else:
						rejected_misps.get(pos).append(change)
					break
				else:
					tmp_accepted_misps.update({pos: (change, read_freq[read[st_find.read_field]], 0)})
		print("accept?",read,tmp_accepted_misps,len(tmp_accepted_misps),len(read[st_find.misp_field]))
		if len(tmp_accepted_misps) == len(read[st_find.misp_field]):

			for kpos in tmp_accepted_misps.keys():
				if kpos not in accepted_misps:
					accepted_misps.update({kpos:tmp_accepted_misps[kpos]})

	restored_accepted_misp = {}
	for aMisp in accepted_misps:
		assert fix_to_ori_mPos[aMisp] in misPs.keys()
		restored_accepted_misp.update({fix_to_ori_mPos[aMisp]:accepted_misps[aMisp]})
	fixed_accepted_misp = accepted_misps
	accepted_misps = restored_accepted_misp
	if len(rejected_read) > 0:
		print(len(rejected_read), "reads rejected")

		rejected_set = set([x[st_find.read_field] for x in rejected_read])

		d_info = {}

		'''if rejected read has more than 1 misp, check if all misp should be reversed'''
		for ir, read in enumerate(new_subbed_read):
			start = int(read[st_find.index_field]) - 1
			cigar = read[st_find.cigar_field]
			#if re.search('^[0-9]+S', cigar) is not None:
			#	assert start != int(original_subbed_read[ir][st_find.index_field]) - 1
			if read[st_find.read_field] in rejected_set:
				#print("revert", read)

				for tmp_misp in read[st_find.misp_field]:
					tmp_misp_pos = list(tmp_misp.keys())[0]
					if tmp_misp_pos not in accepted_misps:
						#assert ref[tmp_misp_pos] != ori_ref[tmp_misp_pos]
						ref = ref[:tmp_misp_pos] + ori_ref[tmp_misp_pos] + ref[tmp_misp_pos+1:]

		idv.write_seq(ref,"verified_"+new_strain_file)


				#ref = ref[:start] + ori_ref[start:start + len(read[st_find.read_field])] + ref[start + len(read[st_find.read_field]):]
			# ref[start:start+len(read[3])] = ori_ref[start:start+len(read[3])]

			# print(pos,change,support_count,tmp_count)

		# remove rejected reads from  subbed_reads and file
		verified_subbed_read = []
		for read in original_subbed_read:
			if read[st_find.read_field] not in rejected_set:
				verified_subbed_read.append(read)

		print(len(subbed_read), "reads originally in ", os.path.basename(sub_read_file))
		subbed_read = verified_subbed_read
		print(len(rejected_set), "reads reverted, removed from", os.path.basename(sub_read_file), "now ", len(subbed_read))
		subbed_read_set = set([x[st_find.read_field] for x in subbed_read])
		subprocess.run("rm -r batch_*",stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

		#new_misp_conflict,no_fill = curr_gap_reads(ref, strain, subbed_read_set, candidate_read,
		#						   raw_matrix_info.narrowed_matrix, read_freq,misp_conflict)
		#misp_conflict = new_misp_conflict
		'''after removing rejected reads, go through accepted_misp again
			writing restored index to rejected and subbed_read file'''

		if stage > 0:
			idv.write_sam(rejected_read, "stage_"+str(stage)+"_rejected_" + os.path.basename(sub_read_file), "w+", False)
			idv.write_sam(verified_subbed_read, "stage_" + str(stage) + "_verified_"+sub_read_file, restore_pos=False)
			idv.write_seq(ref,"stage_"+str(stage)+"_verified_"+new_strain_file)
		idv.write_seq(ref,new_strain_file)
		idv.write_sam(rejected_read, "rejected_" + os.path.basename(sub_read_file), "a+",True)
		idv.write_sam(subbed_read, sub_read_file, restore_pos=True)
	# if reach_end:

	# update misp total freq
	else:
		print("all reads are already verified",len(rejected_read),"reads rejected")
		if ref != unvarified_ref:
			idv.write_seq(ref,"error_ref_"+new_strain_file)
		'''
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
	#verified_read, accpeted_misps = verify_misp(ref_file, matrix_info, sub_read_file, new_strain_file, candidate_sam)
	initial_sub_read = bm.read_sam(sub_read_file)
	tmp_rejected = []
	stage = 1
	tmp_verified_read, tmp_accepted_misps = verify_misp(ref_file, matrix_info, sub_read_file, new_strain_file, candidate_sam,stage)
	misp_conflict = {}
	verified_misp,verified_misp_source,verified_misp_read, tmp_reads = idv.get_misp(idv.get_ref_seq(ref_file), tmp_verified_read,
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
	fixed_candidate_read = idv.fix_s_pos(candidate_read)
	for ir,read in enumerate(candidate_read):
		if re.search('^[0-9]+S',read[st_find.cigar_field]) is not None:
			if int(read[st_find.index_field]) == int(fixed_candidate_read[ir][st_find.index_field]):
				print(read,candidate_read[ir])
				raise "didin't fix s begin"
	candidate_read = fixed_candidate_read
	one_iter = True
	if not one_iter and len(tmp_verified_read) < len(initial_sub_read):
		print("run curr_gap_read")
		exit()
		verified_read = []
		tmp_subbed_read, misp_conflict = curr_gap_reads(idv.get_ref_seq(new_strain_file), strain, tmp_verified_read,
														candidate_read, matrix_info.narrowed_matrix, read_freq,
														misp_conflict, stage)
		if idv.readlist_comapre(tmp_subbed_read,tmp_verified_read) == 0:
			print("no new reads found, terminate")
			verified_read = tmp_subbed_read

	else:
		verified_read = tmp_verified_read

	#keep a global visited read set?
	while not one_iter and len(verified_read) != len(tmp_verified_read):
		print("iterateive curr gap read")
		break
		verified_read = tmp_verified_read
		accepted_misps = tmp_accepted_misps

		stage += 1

		tmp_verified_read, tmp_accepted_misps = verify_misp(ref_file, matrix_info, sub_read_file, new_strain_file,
															candidate_sam,stage)

		test_tmp_verified_read = bm.read_sam(sub_read_file)
		readlist_comapre(tmp_verified_read,verified_read)
		if  len(tmp_verified_read) != len(test_tmp_verified_read):

			readlist_comapre(tmp_verified_read,test_tmp_verified_read)
			raise "verify_misp introduce new reads or remove reads or write reads error"
		curr_verified_read = copy.deepcopy(tmp_verified_read)
		tmp_subbed_read,misp_conflict = curr_gap_reads(get_ref_seq(new_strain_file), strain, curr_verified_read,
													   candidate_read, matrix_info.narrowed_matrix, read_freq,
													   misp_conflict, stage)

		if readlist_comapre(tmp_subbed_read,tmp_verified_read) == 0:
			print("no new found reads, finish verification")
			break
		if len(tmp_subbed_read) < len(tmp_verified_read):
			raise "curr_gap_reads removes verirified reads"
		other_SRR_command = os.path.dirname(__file__) + "/get_combined_extract.sh final_strain_"+str(strain)+"_reference.fa"
		other_SRR_proc = subprocess.run(other_SRR_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
		subprocess.run("cp combine_final_strain_"+str(strain)+"_reference.fa_extract.sam stage_"+str(stage)+"_combine_final_strain_"+str(strain)+"_reference.fa_extract.sam",stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)



	#write misp
	verified_read = tmp_verified_read
	accepted_misps = tmp_accepted_misps
	with open("misp_" + os.path.basename(sub_read_file), "w+") as mf:
		for misp_pos in sorted(list(accepted_misps.keys())):
			val = accepted_misps[misp_pos]
			bases = val[0].split("|")
			print(str(misp_pos + 1), "&", bases[0], "&", bases[1], "&", val[1], "&", val[2], "\\\\\n\hline")
			mf.write(str(misp_pos + 1) + " & " + str(bases[0]) + " & " + bases[1] + "&" + str(val[1]) + " & " + str(
				val[2]) + "\\\\\n\hline\n")


def curr_gap_reads(ref, strain, subbed_read, candidate_read, matrix, read_freq, misp_conflict, stage):
	"""
	Fill more candidate reads after rejecting reads
	:param stage:
	:param ref: reference sequence
	:param strain: strain number
	:param subbed_read: set of remained subbed reads
	:param candidate_read: fixed position
	:param matrix: matrix of other SRR samples
	:param read_freq: frequency of each read
	:param misp_conflict: misp that causes gaps when placed together, format {pos:{pos_misp:[changed_base]}}
	:return:
	"""


	if os.path.exists("half_real_r1.fastq") and os.path.exists("half_real_r2.fastq"):
		readfile1 = "half_real_R1.fastq"
		readfile2 = "half_real_R2.fastq"
	else:
		if os.path.exists("half_real.fastq"):
			readfile1 = "half_real.fastq"
			readfile2 = "none"
		else:
			raise "Error: reduced version reads error"
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

	print(len(prev_readset))

	readlist = candidate_read

	file_subbed_read = bm.read_sam("subbed_reads_" + str(strain) + ".sam")
	fixed_subbed_read = idv.fix_s_pos(file_subbed_read)
	file_subbed_read = fixed_subbed_read
	assert idv.readlist_comapre(file_subbed_read,subbed_read)== 0
	for isr,sread in enumerate(file_subbed_read):
		if int(sread[st_find.index_field]) != int(subbed_read[isr][st_find.index_field]):
			print(sread, subbed_read[isr])
			raise "index error in file_subbed_read"
		start  = int(sread[st_find.index_field])-1
		for ib, base in enumerate(sread[st_find.read_field]):
			covered_pos.update({start + ib: base})
	#misPs = []
	subbed_read_set = set([x[st_find.read_field] for x in subbed_read])
	new_subbed_read = copy.deepcopy(subbed_read)
	subbed_read = new_subbed_read
	#for ir, read in enumerate(readlist):
		#if read[st_find.read_field] in subbed_read_set:
			#subbed_read.append(read)
			#misPs.extend(read[st_find.misp_field])
	original_readlist = copy.deepcopy(readlist)
	misPs, misP_source, misP_reads, new_readlist = idv.get_misp(ref, readlist, fix_pos=False)
	for tmp_read in new_readlist:
		if len(tmp_read[st_find.misp_field]) > 70 :
			print(tmp_read)
			raise "position error"

	readlist = new_readlist

	count = 0
	for batch, read in enumerate(readlist):

		count += 1
		if read[st_find.read_field] in prev_readset:
			# batch += 1
			continue
		if read[st_find.read_field] in subbed_read_set:
			continue
		read_index = int(read[st_find.index_field]) - 1
		#read_misPs = [list(x.items())[0][0] for x in read[st_find.misp_field]]

		print("curr read", batch, read)

		if idv.if_clash(read, covered_pos) or idv.if_reject(read,matrix,read_freq,misP_reads):
			continue

		temp_ref = ref[:read_index] + read[st_find.read_field] + ref[read_index + len(read[st_find.read_field]):]
		# print(editdistance.eval(temp_ref,ref))

		gap_list = identify_verify.verify_seq_support(temp_ref,batch,readfile1,readfile2)
		if len(gap_list) == 0:
			subbed_read.append(candidate_read[batch])
			subbed_read_set.add(read[st_find.read_field])
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
	if len(subbed_read) > len(file_subbed_read):
		idv.\
			write_result(strain, ref, subbed_read, stage, restore_pos=True,curr=True)
	else:
		print("no more reads can fit")
	strain += 1
	# remove temp files
	verify_sub_command = "rm" + " " + "-r" + " " + "./batch_*"

	verify_proc = subprocess.run(verify_sub_command,
								 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
								 shell=True)
	return subbed_read,misp_conflict

def verify_strain(strain,ref_file,candidate_sam):

	backup_command = "cp subbed_reads_" + str(strain) + ".sam subbed_reads_" + str(
		strain) + ".sam.original; cp final_strain_" + str(strain) \
					 + "_reference.fa final_strain_" + str(strain) + "_reference.fa.original"
	backup_proc = subprocess.run(backup_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
								 universal_newlines=True, shell=True)

	other_SRR_command = os.path.dirname(__file__) + "/get_combined_extract.sh final_strain_" + str(
		strain) + "_reference.fa"
	other_SRR_proc = subprocess.run(other_SRR_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

	fac_verify_misp(ref_file, "combine_final_strain_" + str(strain) + "_reference.fa_extract.sam",
					"subbed_reads_" + str(strain) + ".sam", "final_strain_" + str(strain) + "_reference.fa",
					candidate_sam)
	subprocess.run("rm -r batch_* *_out", stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
				   shell=True)

if __name__ == "__main__":
	verify_strain(int(sys.argv[1]),sys.argv[2],sys.argv[3])