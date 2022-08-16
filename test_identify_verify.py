import sys
import unittest
import identify_verify
import strain_finder as st
import build_matrix as bm
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


class TestPotentialMutated(unittest.TestCase):
	readlist = bm.read_sam(sys.argv[1])
	initial_matrix_info = bm.matrix_from_readlist(readlist, 1, set({}), True)

	# if args.brute_force == "True":
	real_narrowed, paired_real_narrowed, nearly_real_narrowed, potential_mutated = bm.narrow_reads(sys.argv[2],
																								   initial_matrix_info.narrowed_read,
																								   "./", True)

	intermit_matrix_info = bm.matrix_from_readlist(paired_real_narrowed, 1, set({}),False,initial_matrix_info,"real_narrowed")
	intermit_matrix_info = bm.matrix_from_readlist(nearly_real_narrowed,1,set({}),False,intermit_matrix_info,"nearly_real_narrowed")



