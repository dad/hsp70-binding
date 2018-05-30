import math, string
import unittest
import scipy as sp
import motif

########################
# TEST CASES
########################
	
class test_scoreRegion(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		s = motif._scoreRegion(seq, motif.rudiger_hsp70_weight_matrix)
		self.assertTrue(abs(s+11.84) < 0.001)

class test_scoreRegion_randomAA(unittest.TestCase):
	def test_run(self):
		seq = 'DALLAN@ANDREI'
		# Don't throw an error if there's a character you don't know, just
		# score it as zero.
		s = motif.score(seq, motif.rudiger_hsp70_weight_matrix)

class test_scoreWindows(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		matrix = motif.hsp70_weight_matrix
		window_size = len(matrix['A'])

		s = motif._scoreWindows(seq, matrix)
		self.assertTrue(s.window_size == window_size)
		self.assertTrue(len(s) == len(seq)+window_size-1)
		self.assertTrue(abs(s[window_size-1]-11.84) < 0.001)

class test_scoreWindows(unittest.TestCase):
	def test_run(self):
		seq = 'MNNAANTGTTNESNV'
		matrix = motif.rudiger_hsp70_weight_matrix
		window_size = len(matrix['A'])

		score_res = motif.score(seq, matrix, return_windows=True)
		self.assertTrue(score_res._windows[0]=='-'*(window_size-1)+seq[0])
		self.assertTrue(score_res._windows[-1]==seq[-1]+'-'*(window_size-1))
		#for (s,r) in zip(scores,regions):
		#	print s,r

class test_scoreResidues(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		matrix = motif.rudiger_hsp70_weight_matrix
		window_size = len(matrix['A'])
		res = motif.score(seq, matrix)
		rs = [s for s in res.residue_scores]
		self.assertTrue(len(rs)==len(seq))

class test_scoreResidues_hsf(unittest.TestCase):
	def test_run(self):
		seq = 'MNNAANTGTTNESNV'
		matrix = motif.rudiger_hsp70_weight_matrix
		window_size = len(matrix['A'])
		res = motif.score(seq, matrix)
		rs = [s for s in res.residue_scores]
		self.assertTrue(len(rs)==len(seq))

class test_scoreWindows_individual(unittest.TestCase):
	def test_run(self):
		matrix = motif.rudiger_hsp70_weight_matrix
		window_size = len(matrix['A'])
		for aa in 'ACDEFGHIKLMNPQRSTVWY':
			res = motif.score(aa, matrix)
			self.assertTrue(res._scores[0] == matrix[aa][-1]) # score in first window should always be the score for the first residue in the last position.
			#print scores

if __name__=="__main__":
	# Run tests
	unittest.main(verbosity=3)