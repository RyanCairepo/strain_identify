import os.path
import subprocess
import sys
import numpy as np
import identify_strain,re
import build_matrix as bm
import strain_finder as st_find

def misp_syno(ref_file,sub_file,code_file,pos_file, read1,read2,offset=0):
	endpoints_list = []
	region_list = []
	with open(pos_file, "r") as ppf:
		for line in ppf:
			fields = line.split(":")
			endpoints_list.append(fields[1].strip())
			region_list.append(fields[0])
	protein_loc = []#(265,21555)
	for se in endpoints_list:
		start = se.split("..")[0]
		end = se.split("..")[1]
		protein_loc.append(range(int(start)-1-offset,int(end)-1-offset))

	assert len(protein_loc) == len(region_list),str(len(protein_loc))+" "+str(protein_loc)
	translate = {}
	ref = ""
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

	misP,misp_source,misp_reads,subbed_read = identify_strain.get_misp(ref, subbed_read, False)
	print(len(misP))

	cov_matrix = bm.matrix_from_readlist(subbed_read, 0.9, set(), target="real_narrowed",fix_s_pos=False)
	syn_stat = {}
	move = []
	syn_count = 0
	for pos,clist in sorted([ (k,v) for k,v in misP.items()],key=lambda x:x[0]):
		found = False
		for pi,r in enumerate(protein_loc):
			#move bases inside this protein region only

			if pos in r:
				print("region:",region_list[pi])
				found = True
				for change in clist:
					index = (pos - r.start) % 3
					changed_base = change.split("|")[1]
					original_base = change.split("|")[0]

					ref_codon = ref[pos-index:pos+(3-index)]
					if changed_base != "-" and original_base != "-":
						# need to consider multiple misp together
						mut_codon = ""
						tmp_read = []
						for read in subbed_read:
							if misp_reads[pos][changed_base][0][st_find.read_field] == read[st_find.read_field]:
								tmp_read = read
						assert len(tmp_read) != 0
						#print(tmp_read)
						#print([list(x.items()) for x in tmp_read[st_find.misp_field]])
						tmp_misP_set = set([list(x.items())[0][0] for x in tmp_read[st_find.misp_field]])
						read_index = int(tmp_read[2])-1
						# mut_codon = ref_codon[0:index] + sorted_mutated_read[smr_index][read_field][mp-smr_read_index] + ref_codon[index + 1:]
						test_ref_codon = ""
						for read_misp in tmp_read[st_find.misp_field]:
							mp = list(read_misp.items())[0][0]

							if mp != pos:
								#print(mp,pos)
								continue

							for i in range(mp - index, mp + (3 - index)):
								#print(i)
								#print(tmp_read[st_find.read_field][i - read_index])
								# if the adjacent spots are also misPs, will use the read's adjecent spot instead of ref
								if i < len(ref):
									if i in tmp_misP_set and 0 <= i - read_index < len(tmp_read[st_find.read_field]):
										mut_codon += tmp_read[st_find.read_field][i - read_index]
										test_ref_codon+=ref[i]
									else:
										mut_codon += ref[i]
										test_ref_codon += ref[i]
								else:
									mut_codon = "XXX"
						assert len(mut_codon) ==3
						assert ref_codon == test_ref_codon,ref_codon+" "+test_ref_codon



						#print(ref_codon,mut_codon,index, ref[pos-index],change,ref[pos+(3-index)-1], translate[ref_codon],translate[mut_codon])

						del_move = len([x for x in move if x < pos])

						if mut_codon not in translate.keys() or translate[ref_codon] != translate[mut_codon]:
							print("non-synonymous", pos + 1 + offset, change)
							syn_stat.update({pos: "non-synonymous"})
						else:
							print("synonymous ", pos + 1 + offset, change)
							syn_stat.update({pos: "synonymous"})
							syn_count += 1
					else:
						print("non-synonymous",pos+1 + offset,change)
						syn_stat.update({pos :"non-synonymous"})
		if not found:
			syn_stat.update({pos:"unavailable"})
	print("non-synonymous/synonymous ratio", syn_count,len(syn_stat),(len(syn_stat)-syn_count)/syn_count)

	for pos in misP:
		tmp = np.squeeze(cov_matrix.narrowed_matrix.getcol(pos).toarray())
		tmp_count = np.bincount(tmp)[1:]
		change = misP[pos]
		changed_base = misP[pos][0].split("|")[1]
		base_num = st_find.dict_tonu[changed_base]
		read_num = np.where(tmp==base_num)

		for matrix_val in tmp[tmp!=0]:
			assert matrix_val == base_num,str(matrix_val)+" "+str(base_num)+" "+str(misp_reads[pos])
		total_freq = sum([subbed_read[x][st_find.freq_field] for x in read_num[0]])
		print(change,misp_source[pos],total_freq)
		assert total_freq > 1, str(total_freq)+" "+str(pos)+" "+str(misp_reads[pos])
		misP[pos][0] += "|" + str(total_freq)

	identify_strain.write_misp(misP,"misp_"+os.path.basename(sub_file))

	with open("misp_"+os.path.basename(sub_file),"r") as mif:
		for line in mif:
			m = re.search('^([0-9]+)',line)
			if m is None:
				continue
			index = int(m.group(1)) -1
			new_line = line.replace('\\\\',' & '+syn_stat[index]+'\\\\')
			new_line = new_line.replace(str(index),str(index+1))
			print(new_line,end="")
			print("\hline")

if __name__ == "__main__":
	if len(sys.argv) > 7:
		misp_syno(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],int(sys.argv[7]))
	else:
		misp_syno(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],sys.argv[6])