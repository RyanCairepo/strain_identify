# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-03-14 17:41:08
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-15 12:51:18

from typing import Counter
from Bio import SeqIO
import os
import editdistance
import operator
import gzip
import argparse

class SingletonErrorCorrection():
    def __init__(self, karect, singleton, file_type, output_dir):
        self.karect = karect
        self.singleton = singleton
        self.file_type = file_type
        self.output_dir = output_dir
        bases = self.singleton.split('/')[-1]
        self.base = bases.split('.')

    def frequency(self, f_input1, f_input2):
        record_iterator1 = SeqIO.parse(f_input1, self.file_type)
        record_iterator2 = SeqIO.parse(f_input2, self.file_type)

        seqs_lst1 = []
        for item1 in record_iterator1:
            seqs_lst1.append(str(item1.seq))
        
        seqs_lst2 = []
        for item2 in record_iterator2:
            seqs_lst2.append(str(item2.seq))
        seqs_set2 = set(seqs_lst2)
        seqs_set1 = set(seqs_lst1)

        if len(seqs_lst1) == len(seqs_set1):
            print("No. of sequences(singletons) before karect correction: {}".format(len(seqs_lst1)))
        
        print("No. of sequences after karect correction: {}".format(len(seqs_set2)))

        intersect = seqs_set1 & seqs_set2
        print("No. of the sequences' intersection before and after karect correction: {}".format(len(intersect)))

        union = seqs_set1 | seqs_set2
        print("No. of the sequences' union before and after karect correction: {}".format(len(union)))

        new_sequences = seqs_set2 - seqs_set1
        print("No. of newly generated sequences before and after karect correction: {}".format(len(new_sequences)))

        seqs_lst2_dict = Counter(seqs_lst2)

        # new_sequences_freq_lst = []
        # intersect_sequences_freq_lst = []
        new_sequences_lst = []
        intersect_sequences_lst = []
        for key, v in seqs_lst2_dict.items():
            if key in intersect and v != 1:
                # intersect_sequences_freq_lst.append(v)
                intersect_sequences_lst.append(key)
            if key in new_sequences and v != 1:
                # new_sequences_freq_lst.append(v)
                new_sequences_lst.append(key)
        
        del seqs_lst2_dict
        del union
        del new_sequences
        del intersect
        del seqs_set1
        del seqs_set2
        del seqs_lst1
        del seqs_lst2
        return intersect_sequences_lst, new_sequences_lst

    def second_correction(self, f_input1, f_input2, intersect_sequences_lst, new_sequences_lst, prefix, frequency_file):
        ############################################################################
        records2 = SeqIO.index(f_input2, self.file_type)

        inter_name_lst = []
        new_name_lst = []
        rest_name_lst = []
        for record in list(records2):
            # print(record)
            # name = records2[record]
            seq = records2[record].seq
            if seq in intersect_sequences_lst:
                inter_name_lst.append(record)
            elif seq in new_sequences_lst:
                new_name_lst.append(record)
            else:
                rest_name_lst.append(record)
        f_output_name = self.output_dir + 'no_change_singletons.' + self.file_type
        no_change_fastq_file = self.extract_records(rest_name_lst, f_input1, f_output_name)
        print("inter_name_lst:{}".format(len(inter_name_lst)))
        print("new_name_lst:{}".format(len(new_name_lst)))
        print("rest_name_lst:{}".format(len(rest_name_lst)))
        del rest_name_lst

        records1 = SeqIO.index(f_input1, self.file_type)

        seq2name_inter_dict = {}
        for name in inter_name_lst:
            seq = records2[name].seq
            seq2name_inter_dict.setdefault(seq,[]).append(name)

        del inter_name_lst
        
        keep_inter_karect_name_lst = []
        keep_inter_original_name_lst = []

        for seq, names in seq2name_inter_dict.items():
            i = 0
            temp_karect_name_lst = []
            temp_original_name_lst = []
            for name in names:
                seq1 = records1[name].seq
                if editdistance.eval(seq, seq1) <= 1:
                    i = i + 1
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
        del seq2name_inter_dict
        # save these frquency changing reads before and after correction
        original_karect_inter_only_fastq_file = self.output_dir + prefix + '_original_karect_inter_only.' + self.file_type            
        self.extract_records(keep_inter_karect_name_lst, f_input1, original_karect_inter_only_fastq_file)
        
        karect_inter_only_fastq_file = self.output_dir + prefix + '_karect_inter_only.' + self.file_type            
        self.extract_records(keep_inter_karect_name_lst, f_input2, karect_inter_only_fastq_file)
        
        with open(frequency_file, 'w') as f:
            f.write('===========================================')
            f.write('\n')
            f.write('Intersection sequences: ')
            f.write('\n')
            f.write('===========================================')
        self.freqency_seq_extraction(original_karect_inter_only_fastq_file, karect_inter_only_fastq_file, frequency_file)

    #############################################################################
        seq2name_new_seq_dict = {}
        for name in new_name_lst:
            seq = records2[name].seq
            seq2name_new_seq_dict.setdefault(seq,[]).append(name)
        del new_name_lst
        keep_new_seq_karect_name_lst = []
        keep_new_seq_original_name_lst = []

        for seq, names in seq2name_new_seq_dict.items():
            i = 0
            temp_karect_name_lst = []
            temp_original_name_lst = []
            for name in names:
                seq1 = records1[name].seq
                if editdistance.eval(seq, seq1) <= 1:
                    i = i + 1
                    temp_karect_name_lst.append(name)
                else:
                    temp_original_name_lst.append(name)
            if len(temp_karect_name_lst) == 1:
                keep_new_seq_original_name_lst.extend(temp_karect_name_lst)
            if len(temp_karect_name_lst) >= 2:
                keep_new_seq_karect_name_lst.extend(temp_karect_name_lst)
            if len(temp_original_name_lst) >=1:
                keep_new_seq_original_name_lst.extend(temp_original_name_lst)
            
            del temp_karect_name_lst
            del temp_original_name_lst
        del seq2name_new_seq_dict
        # save these frquency changing reads before and after correction
        original_karect_new_seq_only_fastq_file = self.output_dir + prefix + '_original_karect_new_seq_only.' + self.file_type            
        self.extract_records(keep_new_seq_karect_name_lst, f_input1, original_karect_new_seq_only_fastq_file)
        karect_new_seq_only_fastq_file = self.output_dir + prefix + '_karect_new_seq_only.' + self.file_type            
        self.extract_records(keep_new_seq_karect_name_lst, f_input2, karect_new_seq_only_fastq_file)
        with open(frequency_file, 'a') as f:
            f.write('===========================================')
            f.write('\n')
            f.write('Newly generated sequences: ')
            f.write('\n')
            f.write('===========================================')
        self.freqency_seq_extraction(original_karect_new_seq_only_fastq_file, karect_new_seq_only_fastq_file, frequency_file)

        # del records1
        # del records2

    ##################################################################################################################
        keep_original_name_lst = list(set(keep_new_seq_original_name_lst).union(set(keep_inter_original_name_lst)))
        keep_original_f = self.output_dir + prefix + '_keep_original_records.' + self.file_type   
        keep_original_fastq_file = self.extract_records(keep_original_name_lst, f_input1, keep_original_f)

        keep_karect_name_lst = list(set(keep_new_seq_karect_name_lst).union(set(keep_inter_karect_name_lst)))
        keep_correct_f = self.output_dir + prefix + '_keep_corrected_records.' + self.file_type            
        keep_correct_fastq_file = self.extract_records(keep_karect_name_lst, f_input2, keep_correct_f)

        del keep_karect_name_lst
        del keep_original_name_lst
        del keep_new_seq_original_name_lst
        del keep_inter_original_name_lst
        del keep_new_seq_karect_name_lst
        del keep_inter_karect_name_lst

        karect_corrected = self.output_dir + 'karect.' + self.file_type
        os.system("cat %s %s > %s" % (keep_original_fastq_file, keep_correct_fastq_file, karect_corrected))

        corrected_singleton = self.output_dir + 'corrected_singleton.' + self.file_type
        os.system("cat %s %s > %s" % (karect_corrected, no_change_fastq_file, corrected_singleton))

        # os.system("rm %s" % (f_output_name))
        os.system("rm %s" % (original_karect_inter_only_fastq_file)) 
        os.system("rm %s" % (karect_inter_only_fastq_file))
        os.system("rm %s" % (original_karect_new_seq_only_fastq_file)) 
        # os.system("rm %s" % (keep_original_f))
        os.system("rm %s" % (keep_original_fastq_file))
        os.system("rm %s" % (keep_correct_fastq_file))
        # os.system("rm %s" % (keep_correct_f))
        os.system("rm %s" % (karect_new_seq_only_fastq_file))
        os.system("rm %s" % (karect_corrected))
        os.system("rm %s" % (no_change_fastq_file))

        return corrected_singleton

    def extract_records(self, name_lst, f_input, f_output_name):
        name_f_out = self.output_dir + 'temp_name.lst'
        if os.path.exists(name_f_out):
            os.system("rm %s" % name_f_out)
        with open(name_f_out, 'w') as f1:
            for item in name_lst:
                f1.write(str(item))
                f1.write('\n')

        if os.path.exists(f_output_name):
            os.system("rm %s" % f_output_name)
        # os.system("wc %s" % name_f_out)
        os.system("seqtk subseq %s %s > %s" % (f_input, name_f_out, f_output_name))
        os.system("rm %s" % (name_f_out))
        return f_output_name

    def karect_correct_singleton(self):
        # # karect correction
        os.system("%s -correct -threads=16 -matchtype=edit -celltype=haploid  -minoverlap=120 -kmer=14 -kmererrors=2 -inputfile=%s -resultdir=%s" % (self.karect, self.singleton, self.output_dir))
        # # os.system("rm ./res_graph* ./temp_res*")
        karect_singleton = self.output_dir + 'karect_' + self.base[0] + '.' + self.file_type
        frequency_file = self.output_dir + self.base[0] + '_frequency.txt'
        intersect_sequences_lst, new_sequences_lst = self.frequency(self.singleton, karect_singleton)

        corrected_singleton = self.second_correction(self.singleton, karect_singleton, intersect_sequences_lst, new_sequences_lst, self.base[0], frequency_file)
        print("Singleton Correction finished!")
        return corrected_singleton

    def seq2name(self, record_iterator):
        seq2name_dict = {}
        for item in record_iterator:
            seq = str(item.seq)
            name = item.id
            seq2name_dict.setdefault(seq, []).append(name)
        return seq2name_dict

    def name2seq(self, record_iterator):
        name2seq_dict = {}
        for item in record_iterator:
            seq = str(item.seq)
            name = item.id
            name2seq_dict.setdefault(name, []).append(seq)
        return name2seq_dict

    def freqency_seq_extraction(self, f_original, f_corrected, frequency_file):
        record_iterator1 = SeqIO.parse(f_original, self.file_type)
        record_iterator2 = SeqIO.parse(f_corrected, self.file_type)
        # seq2name_dict1 = seq2name(record_iterator1)
        name2seq_dict1 = self.name2seq(record_iterator1)
        seq2name_dict2 = self.seq2name(record_iterator2)
        # name2seq_dict2 = name2seq(record_iterator2)

        seq_freq = {}
        for seq, name_lst in seq2name_dict2.items():
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
                f.write(str(seq2name_dict2[seq]))
                f.write('\n')
                for name in seq2name_dict2[seq]:
                    # print(name)
                    f.write(name)
                    f.write('\n')
                    f.write(str(name2seq_dict1[name]))
                    f.write('\n')
                f.write('##########################################')
        return

    def count(self, f_input, file_type):
        record_iterator = SeqIO.parse(f_input, self.file_type)
        seqs_lst = []
        for item in record_iterator:
            seqs_lst.append(str(item.seq))
        print("No. of sequences: {}".format(len(seqs_lst)))
        seqs_lst_dict = Counter(seqs_lst)
        i = 0
        for key, v in seqs_lst_dict.items():
            if v == 1:
                i = i + 1

        print("No. of singleton: {}".format(i))
        return

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

def splitIsolates(fin, result_dir):
    record_iterator, file_type = parse_data(fin)
    seqs2id_dict = {}
    total_seqs = []

    for item in record_iterator:
        seq = str(item.seq)
        total_seqs.append(seq)
        seqs2id_dict.setdefault(seq, []).append(str(item.id))
    # unique_seqs = set(total_seqs)

    isolates_id_lst = []
    non_isolates_id_lst = []
    for seq, id_lst in seqs2id_dict.items():
        if len(id_lst) == 1:
            isolates_id_lst.append(id_lst[0])
        else:
            non_isolates_id_lst.extend(id_lst)

    bases = fin.split('/')[-1]
    base = bases.split('.')
    
    if file_type == 'fastq' or file_type == 'fq' or file_type == 'fastq.gz' or file_type == 'fq.gz':
        isolates_file = result_dir + base[0] + '_isolates.fastq'
        non_isolates_file = result_dir + base[0] + '_non_isolates.fastq'
    elif file_type == 'fasta' or file_type == 'fa' or file_type == 'fasta.gz' or file_type == 'fas.gz':
        isolates_file = result_dir + base[0] + '_isolates.fasta'
        non_isolates_file = result_dir + base[0] + '_non_isolates.fasta' 

    extract_records2(result_dir, isolates_id_lst, fin, isolates_file)
    extract_records2(result_dir, non_isolates_id_lst, fin, non_isolates_file)

    print("Singletons extraction completed.")
    return isolates_file, non_isolates_file, file_type

def extract_records2(working_dir, name_lst, data_set, sub_dataset):
    f_name_lst = os.path.join(working_dir, 'tmp_name.lst')
    with open(f_name_lst, "w") as outfile:
        outfile.write("\n".join(name_lst))
    os.system("seqtk subseq %s %s > %s" % (data_set, f_name_lst, sub_dataset))
    os.system("rm %s" % f_name_lst)
    return

def errorCorrection(karect, f_in, result_dir):
    if os.path.exists(result_dir):
        print("Warning: Directory '% s' already exists, we will use it." % result_dir)
    else:
        os.makedirs(result_dir)
        print("Directory '% s' created" % result_dir)
    isolates_file, non_isolates_file, file_type = splitIsolates(f_in, result_dir)
    SEC = SingletonErrorCorrection(karect, isolates_file, file_type, result_dir)
    corrected_isolates = SEC.karect_correct_singleton()
    bases = f_in.split('/')[-1]
    base = bases.split('.')
    corrected_dataset = result_dir + base[0] + '.corrected.' + file_type
    os.system("cat %s %s > %s" % (corrected_isolates, non_isolates_file, corrected_dataset))
    os.system("rm %s" % (corrected_isolates))
    os.system("rm %s" % (non_isolates_file))
    os.system("rm %s" % (isolates_file))
    os.system("rm *.txt")
    os.system("rm *.fastq")
    print("All done! Have a good day.")
    return corrected_dataset

if __name__ == '__main__':

    # input_file1 = "../../data/propor_corrected_reduced_r1.fastq"
    # result_dir = "./new_result/"
    # karect = "./karect"
    
    parser = argparse.ArgumentParser(description='singletions correction')

    parser.add_argument('-i', '--input', type=str, help='Input file name')
    parser.add_argument('-o', '--output_dir', type=str, default="./result/", help='Output directory')
    parser.add_argument('-k', '--karect', type=str, default="./karect", help='karect full path')


    args = parser.parse_args()

    fin = args.input
    karect = args.karect
    result_dir = args.output_dir

    errorCorrection(karect, fin, result_dir)
    # input_file2 = "../../data/propor_corrected_reduced_r2.fastq"
