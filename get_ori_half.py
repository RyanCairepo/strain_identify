import itertools as its
import multiprocessing as mp
import re
import sys


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


get_ori_half(sys.argv[1],sys.argv[2],sys.argv[3])