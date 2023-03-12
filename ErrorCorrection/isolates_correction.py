# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2022-09-08 22:37:36
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-12 15:33:58

from Bio import SeqIO
import os
import operator
from collections import Counter
import editdistance
import networkx as nx
from tqdm import tqdm
from mpire import WorkerPool
import sys
import gzip
import copy
import time

def sub_base(base):
    if base == 'A':
        return 'TCG'
    elif base == 'T':
        return 'ACG'
    elif base == 'C':
        return 'ATG'
    elif base == 'G':
        return 'ACT'
    elif base == 'N':
        return 'N'

def replace_char(seq, char, index):
    seq[index] = char
    return ''.join(seq)

def seq2substitution(read):
    """
    enumerate all the substitutions for one read 
    Args:
        read (str): a sequence
    Returns:
        set: a set contains reads
    """
    editdis1_list = []
    # substitution
    raw_seq = list(read)
    n = len(raw_seq)
    for i in range(n):
        seq = copy.deepcopy(raw_seq)
        temp_base = seq[i]
        #print(temp_base)
        temp_sub_base = list(sub_base(temp_base))
        #print(temp_sub_base)
        for b in temp_sub_base:
            if b != 'N':
                sub_seq = replace_char(seq, b, i)
                editdis1_list.append(sub_seq)
    return set(editdis1_list)

def seq2deletion(read):
    """
    enumerate all the deletions for one read 
    Args:
        read (str): a sequence
    Returns:
        set: a set contains reads
    """
    editdis1_list = []
    seq = list(read)
    n = len(seq)
    # deletion
    for i in range(n):
        del_seq = read[:i] + read[(i+1):]
        editdis1_list.append(del_seq)
    return set(editdis1_list)

def seq2insertion(read):
    """
    enumerate all the insertions for one read 
    Args:
        read (str): a sequence
    Returns:
        set: a set contains reads
    """
    editdis1_list = []
    seq = list(read)
    n = len(seq)
    # insertion
    bases = ['A', 'G', 'C', 'T']
    for i in range(n+1):
        for b in bases:
            raw_seq = copy.deepcopy(seq)
            raw_seq.insert(i, b)
            editdis1_list.append(''.join(raw_seq))
    return set(editdis1_list)

def enumerate_ed1_seqs(read):
    possible_ed1 = []
    possible_ed1.extend(seq2deletion(read))
    possible_ed1.extend(seq2substitution(read))
    possible_ed1.extend(seq2insertion(read))
    return set(possible_ed1)

def parse_file_type(data_set):
    items = data_set.split(".")
    ext = items[-1]
    if ext == 'fa' or ext == 'fasta':
        f_type = 'fasta'
    elif ext == 'fq' or ext == 'fastq':
        f_type = 'fastq'
    elif ext == 'gz':
        f_type = items[-2] + '.' + ext
    return f_type

def parse_data(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        ff_type = file_type.split('.')[0]
        handle = gzip.open(data_set, 'rt')
        record_iterator = SeqIO.parse(handle, ff_type)
        return record_iterator, ff_type
    else:
        record_iterator = SeqIO.parse(data_set, file_type) 
        return record_iterator, file_type

def real_ed1_seqs(total_seqs, read):
    """
    given a read, generate all its 1nt-edit-distance read counterparts existing in the dataset to form the edges  

    Args:
        total_seqs (list): The list consisting of all reads in the sequencing dataset.
        read (str): A DNA/RNA sequence.

    Returns:
        list: list of tuples of read pairs with only one base different
    """
    possible_ed1 = enumerate_ed1_seqs(read)
    real_seqs =  total_seqs.intersection(possible_ed1)
    read_lst = [read] * len(real_seqs)
    edges = list(zip(read_lst, real_seqs))
    return edges
    
def extract_isolates(input_file, file_type, graph, unique_seqs, seqs2id_dict, result_dir):
    """
    split the source file into two files of isolates and non-isolates based on constructed graph

    Args:
        graph (object): A graph constructed using NetworkX.
        unique_seqs (list): All the unique reads in the dataset
        seqs2id_dict (dict): key: read, id: sequencing id

    Returns:
        files: two files of isolates and non-isolates based on constructed graph
    """
    if not nx.is_connected(graph):
        print("G is a connected graph: {}".format(nx.is_connected(graph)))
        isolates = set(list(nx.isolates(graph)))
    else:
        print("G is a connected graph: {}".format(nx.is_connected(graph)) )

    print("Isolated Nodes Number: {}".format(len(isolates)))

    name_lst = []
    non_name_lst = []
    nonisolated_seqs = unique_seqs - isolates 
    # save isolated nodes as negative samples
    for k in isolates:
        name_lst.extend(seqs2id_dict[k])                  
            
    for s in nonisolated_seqs:
        non_name_lst.extend(seqs2id_dict[s])

    print("isolates name list: {}".format(len(name_lst)))
    print("non-isolates name list: {}".format(len(non_name_lst)))

    bases = input_file.split('/')[-1]
    base = bases.split('.')

    if file_type == 'fastq' or file_type == 'fq' or file_type == 'fastq.gz' or file_type == 'fq.gz':
        isolates_file = result_dir + base[0] + '_isolates.fastq'
        non_isolates_file = result_dir + base[0] + '_non_isolates.fastq'
    elif file_type == 'fasta' or file_type == 'fa' or file_type == 'fasta.gz' or file_type == 'fas.gz':
        isolates_file = result_dir + base[0] + '_isolates.fasta'
        non_isolates_file = result_dir + base[0] + '_non_isolates.fasta' 

    extract_records2(result_dir, name_lst, input_file, isolates_file)
    extract_records2(result_dir, non_name_lst, input_file, non_isolates_file)

    print("Isolated nodes extraction completed.")
    return isolates_file, non_isolates_file

def generate_graph(data_set, high_freq_thre, num_workers):
    """
    construct 1nt- or 2nt-edit-distance-based read graph

    Args:
        data_set (str): The filename including path to be corrected.

    Returns:
        MultiVariables: Multi Variables after constructing edit-distance-based read graph
    """
    record_iterator, file_type = parse_data(data_set)
    seqs2id_dict = {}
    total_seqs = []
    # seq_lens_set = set()
    for item in record_iterator:
        seq = str(item.seq)
        # ll = len(seq)
        # seq_lens_set.add(ll)
        total_seqs.append(seq)
        seqs2id_dict.setdefault(seq, []).append(str(item.id))
    unique_seqs = set(total_seqs)

    graph = nx.Graph()
    read_count = Counter(total_seqs)
    high_freq = []
    low_freq = []
    for read, frequency in tqdm(read_count.items(), miniters=int(len(read_count)/1000)):
        if not graph.has_node(read):
            graph.add_node(read, count = frequency, flag=False)  
        if frequency >= high_freq_thre:
            high_freq.append(read)
        else:
            low_freq.append(read)
    if len(high_freq) == 0:
        print("Error Correction Failed as no high-frequency reads detected.")
        sys.exit(1)
    print(len(read_count))
    ######################################################
    edges_lst = []
    shared_unique_seqs = unique_seqs
    with WorkerPool(num_workers, shared_objects=shared_unique_seqs, start_method='fork') as pool:
        for edge_lst in pool.imap(real_ed1_seqs, high_freq, progress_bar=True):
            edges_lst.extend(edge_lst)

    if len(edges_lst) > 0:
        print(len(edges_lst))
        print(edges_lst[0])
        graph.add_edges_from(edges_lst)

    # return graph, seq_lens_set, seqs2id_dict, unique_seqs, file_type
    return graph, seqs2id_dict, unique_seqs, file_type

def extract_records2(working_dir, name_lst, data_set, sub_dataset):
    f_name_lst = os.path.join(working_dir, 'tmp_name.lst')
    with open(f_name_lst, "w") as outfile:
        outfile.write("\n".join(name_lst))
    os.system("seqtk subseq %s %s > %s" % (data_set, f_name_lst, sub_dataset))
    os.system("rm %s" % f_name_lst)
    return

class IsolatesErrorCorrection():
    def __init__(self, 
            karect, 
            threads, 
            # matchtype, 
            # celltype, 
            minoverlap, 
            kmer, 
            kmererrors, 
            isolates, 
            file_type, 
            output_dir):
            
        self.karect = karect
        self.threads = threads
        # self.matchtype = matchtype
        # self.celltype = celltype
        self.minoverlap = minoverlap
        self.kmer = kmer
        self.kmererrors = kmererrors

        self.isolates = isolates
        self.file_type = file_type
        self.output_dir = output_dir
        bases = self.isolates.split('/')[-1]
        self.base = bases.split('.')

    def select_correction(self, f_original, f_karect, prefix, frequency_file):
        ori_seqs_lst = []
        cor_seqs_lst = []
        ori_seq2name_dict = {}
        karect_seq2name_dict = {}
        for ori_rec, cor_rec in zip(SeqIO.parse(f_original, self.file_type), SeqIO.parse(f_karect, self.file_type)):
            ori_seq = str(ori_rec.seq)
            cor_seq = str(cor_rec.seq)
            ori_name = ori_rec.id
            cor_name = cor_rec.id
            ori_seqs_lst.append(ori_seq)
            cor_seqs_lst.append(cor_seq)
            
            ori_seq2name_dict.setdefault(ori_seq,[]).append(ori_name)
            karect_seq2name_dict.setdefault(cor_seq,[]).append(cor_name)

        ori_seqs_set = set(ori_seqs_lst)
        cor_seqs_set = set(cor_seqs_lst)

        print("No. of sequences(isolates) before karect correction: {}".format(len(ori_seqs_set)))
        if len(ori_seqs_set) == len(ori_seqs_lst):
            print("Each the isolates' frequency equals to 1 before karect correction!")
        else:
            print("Before karect correction, not all the isolates' frequency equals to 1!")
        
        print("No. of unique sequences after karect correction: {}".format(len(cor_seqs_set)))

        intersect = ori_seqs_set & cor_seqs_set
        print("No. of the sequences' intersection before and after karect correction: {}".format(len(intersect)))

        union = cor_seqs_set | ori_seqs_set
        print("No. of the sequences' union before and after karect correction: {}".format(len(union)))

        new_sequences = cor_seqs_set - ori_seqs_set
        print("No. of newly generated sequences before and after karect correction: {}".format(len(new_sequences)))

        cor_seqs_lst_dict = Counter(cor_seqs_lst)
        
        # new_sequences_freq_lst = []
        # intersect_sequences_freq_lst = []
        # new_sequences_lst = []
        # intersect_sequences_lst = []
        # for key, v in cor_seqs_lst_dict.items():
        #     if key in intersect and v != 1:
        #         # intersect_sequences_freq_lst.append(v)
        #         intersect_sequences_lst.append(key)
        #     if key in new_sequences and v != 1:
        #         # new_sequences_freq_lst.append(v)
        #         new_sequences_lst.append(key)
        
        # new_sequences_lst = []
        # intersect_sequences_lst = []
        seqs_not_1 = []
        for key, v in cor_seqs_lst_dict.items():
            if v != 1:
                seqs_not_1.append(key)
        seqs_not_1_set = set(seqs_not_1)
        intersect_sequences_set = seqs_not_1_set & intersect
        new_sequences_set = seqs_not_1_set & new_sequences

        del cor_seqs_lst_dict
        del union
        del new_sequences
        del intersect
        del ori_seqs_lst
        del ori_seqs_set
        del cor_seqs_lst
        del cor_seqs_set
        ############################################################################
        karect_records = SeqIO.index(f_karect, self.file_type)
        original_records = SeqIO.index(f_original, self.file_type)

        keep_inter_karect_name_lst = []
        keep_inter_original_name_lst = []

        for seq in intersect_sequences_set:
            karect_names = karect_seq2name_dict[seq]
            temp_karect_name_lst = []
            temp_original_name_lst = []
            for name in karect_names:
                ori_seq = original_records[name]
                if editdistance.eval(seq, ori_seq) == 1:
                    temp_karect_name_lst.append(name)
                else:
                    temp_original_name_lst.append(name)                    
                    
            if len(temp_karect_name_lst) == 1:
                keep_inter_original_name_lst.extend(temp_karect_name_lst)
            if len(temp_karect_name_lst) >= 2:
                keep_inter_karect_name_lst.extend(temp_karect_name_lst)
            if len(temp_original_name_lst) >=1:
                keep_inter_original_name_lst.extend(temp_original_name_lst)
            del temp_karect_name_lst
            del temp_original_name_lst
            del karect_names
        ######################################################################
        keep_new_seq_karect_name_lst = []
        keep_new_seq_original_name_lst = []

        for new_seq in new_sequences_set:
            karect_names = karect_seq2name_dict[seq]
            temp_karect_name_lst = []
            temp_original_name_lst = []
            for name in karect_names:
                ori_seq = original_records[name].seq
                if editdistance.eval(new_seq, ori_seq) <= 1:
                    temp_karect_name_lst.append(name)
                else:
                    temp_original_name_lst.append(name)
                    
            if len(temp_karect_name_lst) == 1:
                keep_new_seq_original_name_lst.extend(temp_karect_name_lst)
            if len(temp_karect_name_lst) >= 2:
                keep_new_seq_karect_name_lst.extend(temp_karect_name_lst)
            if len(temp_original_name_lst) >=1:
                keep_new_seq_original_name_lst.extend(temp_original_name_lst)
            del karect_names
            del temp_karect_name_lst
            del temp_original_name_lst
        #################################################################################
        # save these frquency changing reads before and after correction
        original_karect_inter_only_fastq_file = prefix + '_original_karect_inter_only.' + self.file_type            
        self.extract_records(keep_inter_karect_name_lst, f_original, original_karect_inter_only_fastq_file)
        
        karect_inter_only_fastq_file = prefix + '_karect_inter_only.' + self.file_type            
        self.extract_records(keep_inter_karect_name_lst, f_karect, karect_inter_only_fastq_file)
        
        with open(frequency_file, 'w') as f:
            f.write('===========================================')
            f.write('\n')
            f.write('Intersection sequences: ')
            f.write('\n')
            f.write('===========================================')
        self.freqency_seq_extraction(self.output_dir + original_karect_inter_only_fastq_file, self.output_dir + karect_inter_only_fastq_file, frequency_file)
        # print("333333333333")
        #############################################################################
        # save these frquency changing reads before and after correction
        original_karect_new_seq_only_fastq_file = prefix + '_original_karect_new_seq_only.' + self.file_type            
        self.extract_records(keep_new_seq_karect_name_lst, f_original, original_karect_new_seq_only_fastq_file)
        karect_new_seq_only_fastq_file = prefix + '_karect_new_seq_only.' + self.file_type            
        self.extract_records(keep_new_seq_karect_name_lst, f_karect, karect_new_seq_only_fastq_file)

        with open(frequency_file, 'a') as f:
            f.write('===========================================')
            f.write('\n')
            f.write('Newly generated sequences: ')
            f.write('\n')
            f.write('===========================================')
        self.freqency_seq_extraction(self.output_dir + original_karect_new_seq_only_fastq_file, self.output_dir + karect_new_seq_only_fastq_file, frequency_file)
        ##################################################################################################################
        keep_original_name_lst = list(set(keep_new_seq_original_name_lst).union(set(keep_inter_original_name_lst)))
        keep_original_f = prefix + '_keep_original_records.' + self.file_type   
        keep_original_fastq_file = self.extract_records(keep_original_name_lst, f_original, keep_original_f)

        keep_karect_name_lst = list(set(keep_new_seq_karect_name_lst).union(set(keep_inter_karect_name_lst)))
        keep_correct_f = prefix + '_keep_corrected_records.' + self.file_type            
        keep_correct_fastq_file = self.extract_records(keep_karect_name_lst, f_karect, keep_correct_f)
        #################################################################################
        f_output_name = 'no_change_isolates.' + self.file_type
        # extract sequences without change before and after correction 
        total_name_lst = list(karect_records)
        rest_name_lst = list(set(total_name_lst) - set(keep_original_name_lst).union(set(keep_karect_name_lst)))
        no_change_fastq_file = self.extract_records(rest_name_lst, f_original, f_output_name)
        print("rest_name_lst:{}".format(len(rest_name_lst)))

        del rest_name_lst
        del keep_karect_name_lst
        del keep_original_name_lst
        del keep_new_seq_original_name_lst
        del keep_inter_original_name_lst
        del keep_new_seq_karect_name_lst
        del keep_inter_karect_name_lst

        karect_corrected = self.output_dir + 'karect.' + self.file_type
        os.system("cat %s %s > %s" % (keep_original_fastq_file, keep_correct_fastq_file, karect_corrected))

        corrected_isolates = self.output_dir + 'corrected_isolates.' + self.file_type
        os.system("cat %s %s > %s" % (karect_corrected, no_change_fastq_file, corrected_isolates))
        return corrected_isolates

    def karect_correct_isolates(self):
        # # karect correction
        t1 = time.time()
        # print(type(self.matchtype), type(self.minoverlap))
        # print(self.matchtype, self.celltype, self.minoverlap, self.kmer, self.kmererrors)
        # os.system("%s -correct -threads=%s -matchtype=%s -celltype=%s  -minoverlap=%s -kmer=%s -kmererrors=%s -inputfile=%s -resultdir=%s" % (self.karect, self.threads, self.matchtype, self.celltype, self.minoverlap, self.kmer, self.kmererrors, self.isolates, self.output_dir))
        try:
            os.system("%s -correct -threads=%s -matchtype=edit -celltype=haploid  -minoverlap=%s -kmer=%s -kmererrors=%s -inputfile=%s -resultdir=%s" % (self.karect, self.threads, self.minoverlap, self.kmer, self.kmererrors, self.isolates, self.output_dir))
        except:
            print("Please check wheter Karect exists and make sure its parameters setting correct!")
        # os.system("rm ./res_graph* ./temp_res*")
        # shutil.move('./res_graph_a.txt', self.output_dir)
        # shutil.move('./res_graph_b.txt', self.output_dir)
        # shutil.move('./input_file.txt', self.output_dir)
        # shutil.move('./*.fastq', self.output_dir)

        # os.remove('./res_graph_a.txt')
        # os.remove('./res_graph_b.txt')
        # os.remove('./input_file.txt')
        # os.remove('./*.fastq')        
        t2 = time.time()
        print("Karect time: {}".format(t2 - t1))
        karect_isolates = self.output_dir + 'karect_' + self.base[0] + '.' + self.file_type
        frequency_file = self.output_dir + self.base[0] + '_frequency.txt'
        corrected_isolates = self.select_correction(self.isolates, karect_isolates, self.base[0], frequency_file)
        t4 = time.time()
        print("Select Correction Time: {}".format(t4-t2))
        print("isolates Correction finished!")
        return corrected_isolates

    def freqency_seq_extraction(self, f_original, f_correct, frequency_file):
        correct_seq2name_dict = {}
        original_name2seq_dict = {}
        for ori_rec, cor_rec in zip(SeqIO.parse(f_original, self.file_type), SeqIO.parse(f_correct, self.file_type)):
            ori_seq = str(ori_rec.seq)
            cor_seq = str(cor_rec.seq)

            ori_id = ori_rec.id
            cor_id = cor_rec.id

            original_name2seq_dict.setdefault(ori_id, []).append(ori_seq)
            correct_seq2name_dict.setdefault(cor_seq, []).append(cor_id)

        seq_freq = {}
        for seq, name_lst in correct_seq2name_dict.items():
            seq_freq[seq] = len(name_lst)

        seq_freq_order = dict(sorted(seq_freq.items(), key=operator.itemgetter(1), reverse=True))

        with open(frequency_file, 'a') as f:        
            f.write('Top frequency: ')
            f.write('\n')
            f.write(str(list(seq_freq_order.values())))
            f.write('\n')
            for seq, freq in seq_freq_order.items():
                f.write('##########################################')
                f.write('\n')
                f.write(str(freq))
                f.write('\n')
                f.write(seq)
                f.write('\n')
                f.write(str(correct_seq2name_dict[seq]))
                f.write('\n')
                for name in correct_seq2name_dict[seq]:
                    # print(name)
                    f.write(name)
                    f.write('\n')
                    f.write(str(original_name2seq_dict[name]))
                    f.write('\n')
                f.write('##########################################')
        return

    def extract_records(self, name_lst, f_input, f_output_name):
        name_f_out = self.output_dir + 'temp_name.lst'
        if os.path.exists(name_f_out):
            os.system("rm %s" % name_f_out)
        with open(name_f_out, 'w') as f1:
            for item in name_lst:
                f1.write(str(item))
                f1.write('\n')

        new_fastq = self.output_dir + f_output_name
        if os.path.exists(new_fastq):
            os.system("rm %s" % new_fastq)
        # os.system("wc %s" % name_f_out)
        os.system("seqtk subseq %s %s > %s" % (f_input, name_f_out, new_fastq))
        return new_fastq    

def ErrorCorrection(input_file, high_freq_thre, num_workers, result_dir, karect, threads, minoverlap, kmer, kmererrors):
    graph, seqs2id_dict, unique_seqs, file_type = generate_graph(input_file, high_freq_thre, num_workers)
    isolates_file, non_isolates_file = extract_isolates(input_file, file_type, graph, unique_seqs, seqs2id_dict, result_dir)
    IEC = IsolatesErrorCorrection(karect, threads, minoverlap, kmer, kmererrors, isolates_file, file_type, result_dir)

    corrected_isolates = IEC.karect_correct_isolates() 
    bases = input_file.split('/')[-1]
    base = bases.split('.')
    corrected_file = result_dir + base[0] + '_corrected.' + file_type
    os.system("cat %s %s > %s" % (corrected_isolates, non_isolates_file, corrected_file))
    return corrected_file

if __name__ == '__main__':
    input_file1 = "/share/rycai/original_reduced_r1.fastq"
    high_freq_thre = 5
    num_workers = 60
    result_dir = "./result/"
    karect = "./karect"
    threads = 12
    minoverlap = 120
    kmer = 14
    kmererrors = 2
    if os.path.exists(result_dir):
        print("Warning: Directory '% s' already exists, we will use it." % result_dir)
    else:
        os.makedirs(result_dir)
        print("Directory '% s' created" % result_dir)
    ErrorCorrection(input_file1, high_freq_thre, num_workers, result_dir, karect, threads, minoverlap, kmer, kmererrors)
    input_file2 = "/share/rycai/original_reduced_r2.fastq"
    ErrorCorrection(input_file2, high_freq_thre, num_workers, result_dir, karect, threads, minoverlap, kmer, kmererrors)

    