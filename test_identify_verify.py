import sys
import unittest
import identify_verify
import countbase
import pytest

class TestGetMisp(unittest.TestCase):
	ref_file = sys.argv[1]
	sub_read_file_list = sys.argv[2:]
	cmisPs, cmisPs_source, cmisP_reads = countbase.get_misp(ref_file, sub_read_file_list)
	imisPs, imisPs_source, imisp_reads = identify_verify.get_misp(ref_file, sub_read_file_list)
	assert len(set.difference(set(cmisPs.keys()),set(imisPs.keys()))) == 0,str(set.difference(set(cmisPs.keys()),set(imisPs.keys())))
	for pos in list(imisPs.keys()):
		ichanges = set(imisPs[pos])
		cchanges = set(cmisPs[pos])
		assert len(set.difference(ichanges,cchanges)) == 0


class Test


