import argparse
import re
import sys
import unittest

import build_matrix
import identify_verify
import strain_finder as st
import build_matrix as bm
import countbase
from unittest.mock import patch
import pytest
filelist = []
with open("test_param.txt","r") as tf:
    for line in tf:
        files = line.strip().split(" ")
        filelist.extend(files)


class TestSubmit(unittest.TestCase):
    skipgetmisp = re.search('sam$', filelist[0]) is not None
    skipmoveback = re.search('subbed', filelist[0]) is None

    '''
    def __init__(self,filelist):
        super(TestSubmit, self).__init__()
        self.filelist = filelist
        if re.search('sam$', self.filelist[0]) is not None:
            self.skipgetmisp = True
        else:
            self.skipgetmisp = False
        self.skipmoveback = re.search('subbed', self.filelist[0]) is None
    
    def setUp(self):
       
        super(TestSubmit, self).setUp()

    '''
    #@unittest.skipIf(skipgetmisp, "skiped as target is " + sys.argv[1])
    def test_GetMisp_1(self):
        if self.skipgetmisp:
            return
        ref_file = filelist[0]
        sub_read_file_list = filelist[1:]
        sub_read = []
        with open(sub_read_file_list[0],"r") as sf:
            for line in sf:
                fields = line.strip().split(" ")
                fields.append(0)
                sub_read.append(fields)
        ref = ""
        with open(ref_file,"r") as rf:
            for line in rf:
                if line[0] != ">":
                    ref += line.strip()

        fix_sub_read = identify_verify.fix_s_pos(sub_read)
        for i,read in enumerate(fix_sub_read):
            begin_s = re.match('^([0-9]+)S' ,read[st.cigar_field])
            if begin_s is not None:
                #print(begin_s.group(1),read,sub_read[i])
                move = int(begin_s.group(1))
                assert int(sub_read[i][st.index_field])- int(read[st.index_field]) == move,str(move)+" "+str(read)+" "+str(sub_read[i])

        cmisPs, cmisPs_source, cmisP_reads = countbase.get_misp(ref_file, sub_read_file_list, False)
        imisPs, imisPs_source, imisp_reads, new_sub_read = identify_verify.get_misp(ref, sub_read)
        assert cmisPs.keys()==imisPs.keys(),str(sorted(list(set.difference(set(imisPs.keys()),set(cmisPs.keys())))))
        for read in new_sub_read:
            if len(new_sub_read[st.misp_field]) > 3:
                print(read)
        for pos in list(imisPs.keys()):
            ichanges = set(imisPs[pos])
            cchanges = set(cmisPs[pos])
            assert len(set.difference(ichanges,cchanges)) == 0
        for i,read in enumerate(new_sub_read):
            begin_s = re.match('^([0-9]+)S' ,read[st.cigar_field])
            if begin_s is not None:
                #print(read,sub_read[i])
                move = int(begin_s.group(1))
                assert int(sub_read[i][st.index_field])- int(read[st.index_field]) == move

    #@unittest.skipIf(skipmoveback,"skipped")
    def test_moveBack(self):

        if self.skipmoveback:
            return
        fixed_list = countbase.fix_sub_pos(filelist[0])
        candidate_list = bm.read_sam(filelist[0])
        candidate_index = {}
        for read in candidate_list:
            candidate_index.update({read[st.read_field]:read[st.index_field]})
        for ir,read in enumerate(fixed_list):
            assert read[st.index_field] == candidate_index[read[st.read_field]]
            if re.search('^[0-9]+S',read[st.cigar_field]) is not None:
                assert int(read[st.index_field]) != int(candidate_list[ir][st.index_field]), str(read[st.index_field])+" "+str(candidate_list[ir][st.index_field])

    def test_Narrowed(self):
        match_limit = 0.95
        extract_reads = build_matrix.read_sam(filelist[0])
        matrix_info = build_matrix.matrix_from_readlist(extract_reads,match_limit,set({}),target="raw")
        print(len(extract_reads),len(matrix_info.narrowed_read))
        for read in matrix_info.narrowed_read:
            matched = re.findall('[0-9]+M',read[st.cigar_field])
            assert matched is not None
            total = 0
            for blk in matched:
                m = re.match('([0-9]+)M',blk)
                if m is not None:
                    total += int(m.group(1))
            #self.assertEqual( total/len(read[st.read_field]) , 0.7)
            if total/len(read[st.read_field]) < match_limit:
                print(total,read)
                raise "narrow error"

    def test_Getmisp(self):
        readlist = []
        ref = ""
        with open(filelist[0], "r") as rf:
            for line in rf:
                if line[0] != ">":
                    ref += line.strip()
        readlist = bm.read_sam(filelist[1],True)

        old_misPs = {}
        old_readlist = []
        with open(filelist[1],"r") as rf:
            for line in rf:
                fields  = line.strip().split(" ")
                fields.append([])
                old_readlist.append(fields)
        old_readlist = identify_verify.fix_s_pos(old_readlist)
        misPs, misP_source, misP_reads, new_readlist = identify_verify.get_misp(ref, old_readlist, False, False)
        for fields in old_readlist:
            read_index = int(fields[st.index_field])-1
            for ir,base in enumerate(fields[st.read_field]):
                if base != ref[read_index+ir]:
                    old_misPs.update({read_index+ir:[ref[read_index+ir]+"|"+base]})
                    fields[st.misp_field].append({read_index+ir:[ref[read_index+ir]+"|"+base]})
        #assert len(misPs) == len(old_misPs), str(len(misPs))+"\n\n"+str(len(old_misPs))
        for ir,read in enumerate(new_readlist):
            assert len(read[st.misp_field]) == len(old_readlist[ir][st.misp_field]),str(read)+"\n"+str(old_readlist[ir])
            assert len(read[st.misp_field]) > 0
            for im,mp in enumerate(read[st.misp_field]):
                #print(mp,list(mp.items()))

                assert int(list(mp.keys())[0]) == int(list(old_readlist[ir][st.misp_field][im].keys())[0]),str(list(mp.keys())[0])+"\n"+str(list(old_readlist[ir][st.misp_field][im].keys())[0])
        #for key in misPs:
            #assert key in old_misPs.keys()
            #assert misPs[key][0] == old_misPs[key][0]
    def test_misp_read(self):
        ref = identify_verify.get_ref_seq(filelist[0])

        subbed_read = bm.read_sam(filelist[1],True)

        misPs,misP_source,misP_reads,new_subbed = identify_verify.get_misp(ref,subbed_read)
        new_read_dict = {}
        for read in new_subbed:
            assert len(read[st.misp_field])>0
            new_read_dict.update({read[st.read_field]:read})

        for pos in misP_reads:
            for change in misP_reads[pos].keys():
                if len(misP_reads[pos][change]) > 1:
                    print("========")
                    print(pos,change)
                    for read in misP_reads[pos][change]:
                        print(new_read_dict[read[st.read_field]])


    def test_seqMisp(self):
        ref1 = identify_verify.get_ref_seq(filelist[0])
        ref2 = identify_verify.get_ref_seq(filelist[2])
        readlist = bm.read_sam(filelist[1],True)
        misPs, misP_source, misP_reads, new_readlist = identify_verify.get_misp(ref1,readlist,False,True)
        for i,base in enumerate(ref1):
            if base != ref2[i]:
                assert i in misPs.keys() and ref2[i] == misPs[i][0].split("|")[1]

    def test_consecMisp(self):
        readlist = bm.read_sam(filelist[1])
        ref = identify_verify.get_ref_seq(filelist[0])
        misP,misP_source,misP_reads, new_readlist = identify_verify.get_misp(ref,readlist)
        read_dict = {}
        for read in new_readlist:
            read_dict.update({read[3]:read})
        consec_misp_reads = []
        last_pos = 0
        curr_count = 0
        sorted_pos = sorted(list(misP_reads.keys()))
        added_read = set({})
        for ip,pos in enumerate(sorted_pos):
            if ip>0:
                if pos - last_pos <= 3:
                    curr_count += 1
                else:
                    if curr_count > 10:
                        #print(curr_count,ip,[sorted_pos[x] for x in range(ip-curr_count,ip)],misP_reads[sorted_pos[ip-curr_count+1]])
                        if list(misP_reads[sorted_pos[ip-1]].values())[0][0][3] not in added_read:

                            added_read.add(list(misP_reads[sorted_pos[ip-1]].values())[0][0][3])
                            consec_misp_reads.extend([list(misP_reads[sorted_pos[x]].values())[0][0] for x in range(ip-curr_count,ip) ])
                    curr_count = 0

            last_pos = pos
        new_consec_reads = []

        #for i in consec_misp_reads:
            #print(read_dict[i[3]])
    def test_misp_syno(self):
        identify_verify.misp_syno(filelist[0],filelist[1],filelist[2],filelist[3])
    def test_verify(self):
        ref = identify_verify.get_ref_seq(filelist[0])
        sub_read = bm.read_sam(filelist[1])
        support_read = bm.read_sam(filelist[2])
        support_matrix_info = bm.matrix_from_readlist(support_read,0.7,set({}))

        misP, misP_source, misP_reads, new_readlist = identify_verify.get_misp(ref, sub_read)
        ver_read, ver_misp = identify_verify.verify_misp(filelist[0],support_matrix_info,filelist[1],filelist[3],filelist[4],stage=2)

    def test_misp_seq(self):
        ref = identify_verify.get_ref_seq(filelist[0])
        sub_read = bm.read_sam(filelist[1])
        new_ref = identify_verify.get_ref_seq(filelist[2])
        assert len(ref) == len(new_ref),str(len(ref))+" "+str(len(new_ref))
        #sub_read = identify_verify.fix_s_pos(sub_read)
        misPs, misP_source,misP_read,new_sub_read = identify_verify.get_misp(ref,sub_read,fix_pos=True)
        error = False
        for i,base in enumerate(ref):
            if base != new_ref[i]:
                if i not in misPs.keys() or misPs[i][0].split("|")[1] != new_ref[i]:
                    error = True
                    print(i,ref[i],new_ref[i])
                    if i in misPs.keys():
                        print(misPs[i])
                    else:
                        print(i,"not in misP")
        if error:
            raise "missing misp"
    def test_read_compare(self):
        readlist1 = bm.read_sam(filelist[0])
        readlist2 = bm.read_sam(filelist[1])
        assert identify_verify.readlist_comapre(readlist1,readlist2) == 0
    def test_curr_gap_read(self):
        curr_verified_read = bm.read_sam(filelist[1])
        candidate_read = bm.read_sam(filelist[2])
        support_read = bm.read_sam(filelist[3])
        matrix_info = bm.matrix_from_readlist(support_read,0.7,set({}))
        read_freq = {}
        with open("mutated_read_freq.txt","r") as mf:
            for line in mf:
                fields = line.strip().split(": ")
                read_freq.update({fields[0]:int(fields[1])})
        misp_conflict = {}
        stage = 0
        strain = 0
        curr_verified_read = identify_verify.fix_s_pos(curr_verified_read)
        print("originally ",len(curr_verified_read),"verified reads")
        candidate_read = identify_verify.fix_s_pos(candidate_read)
        tmp_subbed_read, misp_conflict = identify_verify.curr_gap_reads(identify_verify.get_ref_seq(filelist[0]), strain, curr_verified_read,
                                                        candidate_read, matrix_info.narrowed_matrix, read_freq,
                                                        misp_conflict, stage)
        identify_verify.write_sam(tmp_subbed_read,"curr_subbed_read_0.sam",restore_pos=True)
        assert identify_verify.readlist_comapre(curr_verified_read,tmp_subbed_read)== 0
if __name__ == '__main__':
    #parser = argparse.ArgumentParser()

    #parser.add_argument("files",type=str,default="",nargs="+")
    #parser.add_argument("--test_args",nargs="*")
    #args = parser.parse_args()
    #filelist = args.files

    #sys.argv[1:] = args.test_args


    #suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestSubmit)
    #suite.addTest(TestSubmit("testNarrowed"))
    #result = unittest.TextTestRunner().run(suite)

    unittest.main()