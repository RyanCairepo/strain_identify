import itertools
import multiprocessing
import os
import sys,re,editdistance,time, multiprocessing as mp, itertools as its

class Param(object):




	def __init__(self,index,splice,start_line,size, step):
		self.index,  self.splice, self.start_line, self.size, self.step = index,splice,start_line,size, step


#def inspect(index,  count,splice, start_line, size, w_reads, w_window, w_lines, lock ):
def inspect(param):
	#index, count, splice, start_line, size, w_reads, w_window, w_lines, lock = \
	#	param.index, param.count, param.splice, param.start_line,param.size, param.w_reads, param.w_window, param.w_lines, param.lock
	#index, splice,start_line,size = param[0].index, param[0].splice, param[0].start_line, param[0].size
	index, splice, start_line, size = param.index, param.splice, param.start_line, param.size
	#file = param[1]
	local_count = 0
	#count.value+=size
	lines = []
	window = []
	local_reads = []
	dis = []
	step = param.step


	while( local_count % step !=1):
		local_count += 1

	while(start_line % step != 1):
		start_line += 1

	with open(sys.argv[1],"r") as file:
		#print(os.getpid()," reading start ",start_line)
		for read in its.islice(file, start_line, start_line + size,step):
			if local_count % step == 1:
				#rread = read.split(' ')[3]
				#print(read, local_count)
				tmp =  editdistance.eval(splice, read)
				if tmp <= distance:

					#lock.acquire()
					#try:
					#	print(local_count+start_line)
					#	print(splice)
					#	print(read)
					#	print(local_count)
					#finally:
					#	lock.release()
					lines.append(local_count+start_line)
					window.append(index)
					local_reads.append(read)
					dis.append(tmp)
				local_count+=step
	'''
	lock.acquire()
	try:
		w_reads.extend(local_reads)
		w_window.extend(window)
		w_lines.extend(lines)
	finally:
		lock.release()
	'''
	return local_reads,window,lines,dis


if __name__ == "__main__":
	splice = ""
	whole_gene=""
	length=300
	step = 4
	distance=130
	#count = mp.Value(int, 0)
	with open("spike.fasta","r") as f:
		for line in f:
			if re.search('^[ATGCUatcgu]+', line):
				whole_gene = line
				break
	with open("naive_match_615.txt","w+") as f:
		f.write("")
	if len(whole_gene) < length:
		print(len(whole_gene))
		print(whole_gene)
		exit(1)
	l = mp.Lock()
	line_count = mp.Value('I',0)
	core_num = mp.cpu_count()-2
	line_amount = 245216120
	#line_amount = 10099368

	index = 0
	matched_reads = []
	matched_window = []
	matched_lines = []
	match_count = 0

	n_lines = []
	n_window = []
	n_reads = []
	n_dis = []
	for index in range(1, 315):
		line_count.value = 0
		splice=whole_gene[index:index+length]
		#print("start searching for matches for ",splice," at location ",str(index))
		start = 0
		plist=[]
		#inspect(Param(index, splice,start,int(line_amount)))


		with mp.Pool(core_num) as pool:
			lparam = []
			#with open(sys.argv[1], "r") as f:


			while start < line_amount:
				#iterfile = itertools.islice(f, start, start + int(line_amount / core_num) + 1, 4)
				lparam.append((Param(index, splice, start, int(line_amount/core_num),step)))
				start += int(line_amount/core_num)+1

			results=pool.map(inspect,lparam)
		
		#print(results)
		#print(results[1])
		#print(results[2])

		for i in results:
			#print(i[0])
			#print(i[1])
			#print(i[2])

			n_reads.extend(i[0])
			n_window.extend(i[1])
			n_lines.extend(i[2])
			n_dis.extend(i[3])
			match_count += len(n_lines)


		print("found ", len(n_lines), " reads matches")
		if(len(n_lines)==0) :
			continue
		n_dis, n_lines, n_reads, n_window = zip(*sorted(zip(n_dis, n_lines, n_reads, n_window)))
		with open("naive_match_615.txt", "a") as log:
			for i in range(len(n_lines)):
				log.write("match at " + str(n_lines[i]) + " index at " + str(n_window[i]) + " with distance " + str(
					n_dis[i]) + "\nreads:" + n_reads[i] + " current found matches: " + str(match_count)+"\n")
		print(len(n_lines)," matches found at location ",index)
		n_reads = []
		n_window = []
		n_lines = []
		n_dis =[]
		#print(len(n_reads))
		#print(len(n_window))
		#exit(1)

	#n_dis, n_lines, n_reads, n_window = zip(*sorted(zip(n_dis, n_lines, n_reads, n_window)))
	print("found ",len(n_lines)," reads matches, recorded in naive_match_615.txt")
	#with open("naive_match_615.txt", "w+") as log:
	#	for i in range(len(n_lines)):
	#		log.write("match at " + str(n_lines[i])+" index at "+str(n_window[i])+" with distance "+ str(n_dis[i])+"\nreads:"+n_reads[i]+"\n")



	#print(len(lines)," of matched reads found, distance <= ",distance)
	#print(lines)