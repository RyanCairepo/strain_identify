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
    def test_MoveBack(self):

        if self.skipmoveback:
            return
        fixed_list = countbase.fix_sub_pos(filelist[0])
        candidate_list = bm.read_sam(filelist[1])
        candidate_index = {}
        for read in candidate_list:
            candidate_index.update({read[st.read_field]:read[st.index_field]})
        for read in fixed_list:
            assert read[st.index_field] == candidate_index[read[st.read_field]]

    def test_Narrowed(self):
        extract_reads = build_matrix.read_sam(filelist[0])
        matrix_info = build_matrix.matrix_from_readlist(extract_reads,0.7,set({}),target="read_narrowed")
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
            assert total/len(read[st.read_field]) >= 0.7
    def test_Getmisp(self):
        readlist = []
        ref = ""
        with open(filelist[0], "r") as rf:
            for line in rf:
                if line[0] != ">":
                    ref += line.strip()
        readlist = bm.read_sam(filelist[1],True)
        misPs, misP_source, misP_reads, new_readlist = identify_verify.get_misp(ref, readlist, False, True)
        old_misPs = {}
        old_readlist = []
        with open(filelist[1],"r") as rf:
            for line in rf:
                fields  = line.strip().split(" ")
                old_readlist.append(fields)
        old_readlist = identify_verify.fix_s_pos(old_readlist)
        for fields in old_readlist:
            read_index = int(fields[st.index_field])-1
            for ir,base in enumerate(fields[st.read_field]):
                if base != ref[read_index+ir]:
                    old_misPs.update({read_index+ir:[ref[read_index+ir]+"|"+base]})
        assert len(misPs) == len(old_misPs), str(misPs)+"\n"+str(old_misPs)
        for read in new_readlist:
            for mp in read[st.misp_field]:
                print(mp,list(mp.items()))
                assert len(list(mp.items())[0][1]) == 1
        for key in misPs:
            assert key in old_misPs.keys()
            assert misPs[key][0] == old_misPs[key][0]
    def test_consecMisp(self):
        readlist = bm.read_sam(filelist[1])
        ref = identify_verify.get_ref_seq(filelist[0])
        misP,misP_source,misP_reads, new_readlist = identify_verify.get_misp(ref,readlist)
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
        for i in consec_misp_reads:
            print(i)
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