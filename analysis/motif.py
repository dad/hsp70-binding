import math, string, random
import unittest
import scipy as sp


class ScoreResult:
	def __init__(self):
		self.scores = []
		self.regions = None
		self.window_size = None
	
	def len(self):
		return len(self.scores)
	
	def __len__(self):
		return len(self.scores)

class ScoreSummary:
	def __init__(self):
		self.run_frequency = None
		self.min_score = None
		
def summarizeScores(residue_scores, threshold, num_freqs):
	score_sum = ScoreSummary()
	# Retrieve histogram of lengths
	score_string = ''
	good_char = 'Y'
	for s in residue_scores:
		if s <= threshold:
			score_string += good_char
		else:
			score_string += ' '
	runs = score_string.split()
	run_frequency = [0]*num_freqs
	for i in range(num_freqs):
		ny = good_char*(i+1)
		run_frequency[i] = runs.count(ny)
	score_sum.run_frequency = run_frequency
	# Check to see if we've covered all sites
	if sum(run_frequency) < len(runs):
		score_sum.num_longer_runs = len(runs) - sum(run_frequency)
	else:
		assert sum(run_frequency) == len(runs)
		score_sum.num_longer_runs = 0
	
	score_sum.min_score = min(residue_scores)
	score_sum.num_motifs = len(runs)
	score_sum.num_sites = score_string.count(good_char)
	return score_sum

hsp70_weight_matrix = {
	'A':[-0.02,-0.05,-0.07,-0.11,0.79,0.79,0.79,0.79,0.79,0.69,0.46,0.3,0.15],
	'C':[1.61,3.21,4.87,7.3,6.35,6.35,6.35,6.35,6.35,0.37,0.25,0.16,0.08],
	'D':[0.14,0.29,0.44,0.65,4.91,4.91,4.91,4.91,4.91,0.53,0.35,0.23,0.12],
	'E':[0.49,0.98,1.48,2.22,5.14,5.14,5.14,5.14,5.14,2.47,1.65,1.09,0.54],
	'F':[0.05,0.09,0.14,0.21,-1.17,-1.17,-1.17,-1.17,-1.17,0.79,0.53,0.35,0.17],
	'G':[-0.11,-0.22,-0.33,-0.5,1.95,1.95,1.95,1.95,1.95,0.05,0.03,0.02,0.01],
	'H':[-0.08,-0.16,-0.24,-0.37,1.74,1.74,1.74,1.74,1.74,0.13,0.09,0.06,0.03],
	'I':[0.34,0.69,1.04,1.56,-2.05,-2.05,-2.05,-2.05,-2.05,0.17,0.11,0.08,0.04],
	'K':[-0.28,-0.56,-0.85,-1.28,0.4,0.4,0.4,0.4,0.4,-1.62,-1.08,-0.71,-0.36],
	'L':[0.56,1.12,1.7,2.54,-3.62,-3.62,-3.62,-3.62,-3.62,-0.03,-0.02,-0.01,-0.01],
	'M':[0.03,0.05,0.08,0.12,1.1,1.1,1.1,1.1,1.1,0.26,0.17,0.11,0.06],
	'N':[0.25,0.49,0.74,1.12,2.36,2.36,2.36,2.36,2.36,-0.44,-0.29,-0.19,-0.1],
	'P':[0.05,0.1,0.15,0.22,1.63,1.63,1.63,1.63,1.63,-0.42,-0.28,-0.18,-0.09],
	'Q':[-0.37,-0.74,-1.13,-1.69,1.6,1.6,1.6,1.6,1.6,-0.22,-0.15,-0.1,-0.05],
	'R':[-0.39,-0.78,-1.19,-1.78,-0.79,-0.79,-0.79,-0.79,-0.79,-2.58,-1.72,-1.14,-0.57],
	'S':[-0.04,-0.09,-0.13,-0.2,1.27,1.27,1.27,1.27,1.27,-0.34,-0.23,-0.15,-0.07],
	'T':[-0.3,-0.6,-0.91,-1.36,0.27,0.27,0.27,0.27,0.27,-0.73,-0.48,-0.32,-0.16],
	'V':[-0.08,-0.17,-0.26,-0.39,-1.75,-1.75,-1.75,-1.75,-1.75,1.05,0.7,0.46,0.23],
	'W':[-0.14,-0.29,-0.43,-0.65,3.49,3.49,3.49,3.49,3.49,0.17,0.12,0.08,0.04],
	'X':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'-':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'Y':[0.06,0.13,0.19,0.29,-1.88,-1.88,-1.88,-1.88,-1.88,1.73,1.15,0.76,0.38]
}

# For testing
flat_weight_matrix = {
	'A':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'C':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'D':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'E':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'F':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'G':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'H':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'I':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'K':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'L':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'M':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'N':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'P':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'Q':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'R':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'S':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'T':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'V':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'W':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'X':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'-':[1,1,1,1,1,1,1,1,1,1,1,1,1],
	'Y':[1,1,1,1,1,1,1,1,1,1,1,1,1]
}

def scoreRegion(seq, matrix): #, exclude_from_left=0, exclude_from_right=0):
	"""Score a region of the protein. This sequence assumes that the input sequence is no longer than the window-size
	specified by the scoring matrix.
	
		@return a number, the sum of the scores.
	"""
	n = len(seq)
	assert(n <= len(matrix.values()[0]))
	try:
		scores = [matrix[x][i] for (x,i) in zip([a for a in seq], range(n))]
	except KeyError, ke:
		# Go through and be more careful if necessary
		scores = [matrix[x][i] for (x,i) in zip([a for a in seq], range(n)) if matrix.has_key(x)]
	return sum(scores)
	
def scoreWindows(seq, matrix, return_regions=False):
	# Scan sequence
	window_size = len(matrix.values()[0])
	# Pad the sequence so that sliding window does not require any special-casing
	pad = '-'
	assert(matrix[pad][0] == 0) # Check to make sure padding does not alter the score.
	padded_seq = pad*window_size + seq + pad*window_size
	n = len(seq)
	res = ScoreResult()
	res.window_size = window_size
	#scores = []
	if return_regions:
		res.regions = []
	for i in range(len(seq)+window_size-1):
		region = padded_seq[(i+1):(i+window_size+1)]
		assert len(region) == window_size
		region_score = scoreRegion(region, matrix)
		res.scores.append(region_score)
		if return_regions:
			res.regions.append(region)
	return res

def scoreResidues(seq, scoreresult, window_size): #, summaryFunction=sp.mean):
	# Assumption is that scores is a ScoreResult generated by scoreProtein
	# score for a residue is defined as the score for the window centered on this residue
	assert(len(scoreresult) == len(seq)+window_size-1)
	mid_window = int(math.floor(window_size/2.0))
	residue_scores = [0]*len(seq)
	for aai in range(len(seq)):
		pos = aai + mid_window
		residue_scores[aai] = scoreresult.scores[pos]
	return residue_scores

def getResidueScores(seq, matrix):
	score_result = scoreWindows(seq, matrix, return_regions=False)
	residue_scores = scoreResidues(seq, score_result, window_size=len(matrix.values()[0]))
	return residue_scores
	

########################
# TEST CASES
########################
	
class test_scoreRegion(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		s = scoreRegion(seq, hsp70_weight_matrix)
		self.assertTrue(abs(s-11.84) < 0.001)

class test_scoreRegion_randomAA(unittest.TestCase):
	def test_run(self):
		seq = 'DALLAN@ANDREI'
		# Don't throw an error if there's a character you don't know, just
		# score it as zero.
		s = scoreRegion(seq, hsp70_weight_matrix)

class test_scoreWindows(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		matrix = hsp70_weight_matrix
		window_size = len(matrix.values()[0])

		s = scoreWindows(seq, matrix)
		self.assertTrue(len(s) == len(seq)+window_size-1)
		self.assertTrue(abs(s[window_size-1]-11.84) < 0.001)

class test_scoreWindows(unittest.TestCase):
	def test_run(self):
		seq = 'MNNAANTGTTNESNV'
		matrix = hsp70_weight_matrix
		window_size = len(matrix.values()[0])

		scores = scoreWindows(seq, matrix, return_regions=True)
		self.assertTrue(scores.regions[0]=='-'*(window_size-1)+seq[0])
		self.assertTrue(scores.regions[-1]==seq[-1]+'-'*(window_size-1))
		#for (s,r) in zip(scores,regions):
		#	print s,r

class test_scoreResidues(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		matrix = hsp70_weight_matrix
		window_size = len(matrix.values()[0])
		scores = scoreWindows(seq, matrix)
		resscores = scoreResidues(seq, scores, window_size)
		self.assertTrue(len(resscores)==len(seq))

class test_scoreResidues_hsf(unittest.TestCase):
	def test_run(self):
		seq = 'MNNAANTGTTNESNV'
		matrix = hsp70_weight_matrix
		window_size = len(matrix.values()[0])
		scores = scoreWindows(seq, matrix)
		resscores = scoreResidues(seq, scores, window_size)
		print scores
		print resscores
		self.assertTrue(len(resscores)==len(seq))

class test_scoreWindows_individual(unittest.TestCase):
	def test_run(self):
		matrix = hsp70_weight_matrix
		window_size = len(matrix.values()[0])
		for aa in 'ACDEFGHIKLMNPQRSTVWY':
			scores = scoreWindows(aa, matrix)
			self.assertTrue(scores.scores[0] == matrix[aa][-1]) # score in first window should always be the score for the first residue in the last position.
			#print scores

if __name__=="__main__":
	# Run tests
	unittest.main(verbosity=3)
	
	