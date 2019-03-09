import math, string
import unittest
import scipy as sp


class ScoreEntry:
	def __init__(self, pos, score, residue, window):
		self.score = score
		self.pos = pos # 1-based index!
		self.window = window
		self.residue = residue

class ScoreSummary:
	def __init__(self):
		self.run_frequency = None
		self.min_score = None

class ScoreResult:
	def __init__(self, seq, matrix):
		self._scores = []
		self._sequence = seq
		self._matrix = matrix
		self._windows = []
		self._window_size = len(list(matrix.values())[0])
		self._mid_window = int(math.floor(self._window_size/2.0))
	
	def len(self):
		return len(self._scores)
	
	def positionScores(self):
		for s in self._scores:
			yield s
	
	def __len__(self):
		return len(self._scores)
	
	def addScore(self, score):
		self._scores.append(score)
	
	def addWindow(self, window):
		self._windows.append(window)

	def summary(self, threshold, num_freqs):
		score_sum = ScoreSummary()
		# Retrieve histogram of lengths
		score_string = ''
		good_char = 'Y'
		for s in self.results:
			if s.pos>0 & s.pos<=len(self._sequence):
				if s.score >= threshold:
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

		score_sum.min_score = min(self._scores)
		score_sum.max_score = max(self._scores)
		score_sum.num_motifs = len(runs)
		score_sum.num_sites = score_string.count(good_char)
		return score_sum
	
	def maskSequence(self, threshold, mask_char='x', cmp=lambda x,t: x<t):
		assert len(mask_char)==1
		mask = [cmp(s,threshold) for s in self.residue_scores]
		masked_seq = ''
		seq = self._sequence
		for i in range(len(seq)):
			if mask[i]:
				masked_seq += mask_char
			else:
				masked_seq += seq[i]
		return masked_seq
	
	def __getitem__(self, id):
		return self._scores[id]
	
	@property
	def scores(self):
		for s in self._scores:
			yield s

	@property
	def residue_scores(self):
		seqstart = self._mid_window
		seqend = self._mid_window + len(self._sequence)
		for s in self._scores[seqstart:seqend]:
			yield s

	@property
	def window_size(self):
		return self._window_size

	@property
	def mid_window(self):
		return self._mid_window

	@property
	def windows(self):
		for s in self._windows:
			yield s
	
	@property
	def results(self):
		for pos in range(len(self._scores)):
			seqpos = pos-self._mid_window
			window = None
			if len(self._windows) > 0:
				window = self._windows[pos]
			residue = '-'
			if seqpos >= 0 and seqpos < len(self._sequence):
				residue = self._sequence[seqpos]
			entry = ScoreEntry(seqpos, self._scores[pos], residue, window)
			yield entry
	
	@property
	def residue_results(self):
		for sentry in self.results:
			seqpos = sentry.pos-self._mid_window
			if seqpos >= 0 and seqpos < len(self._sequence):
				yield sentry



# From Rudiger S, Germeroth L, Schneider-Mergener J, Bukau B (1997) Substrate specificity of the DnaK chaperone determined by screening cellulose-bound peptide libraries. EMBO J 16: 1501-1507. 
rudiger_hsp70_weight_matrix = {
	'A':[0.02,0.05,0.07,0.11,-0.79,-0.79,-0.79,-0.79,-0.79,-0.69,-0.46,-0.3,-0.15],
	'C':[-1.61,-3.21,-4.87,-7.3,-6.35,-6.35,-6.35,-6.35,-6.35,-0.37,-0.25,-0.16,-0.08],
	'D':[-0.14,-0.29,-0.44,-0.65,-4.91,-4.91,-4.91,-4.91,-4.91,-0.53,-0.35,-0.23,-0.12],
	'E':[-0.49,-0.98,-1.48,-2.22,-5.14,-5.14,-5.14,-5.14,-5.14,-2.47,-1.65,-1.09,-0.54],
	'F':[-0.05,-0.09,-0.14,-0.21,1.17,1.17,1.17,1.17,1.17,-0.79,-0.53,-0.35,-0.17],
	'G':[0.11,0.22,0.33,0.5,-1.95,-1.95,-1.95,-1.95,-1.95,-0.05,-0.03,-0.02,-0.01],
	'H':[0.08,0.16,0.24,0.37,-1.74,-1.74,-1.74,-1.74,-1.74,-0.13,-0.09,-0.06,-0.03],
	'I':[-0.34,-0.69,-1.04,-1.56,2.05,2.05,2.05,2.05,2.05,-0.17,-0.11,-0.08,-0.04],
	'K':[0.28,0.56,0.85,1.28,-0.4,-0.4,-0.4,-0.4,-0.4,1.62,1.08,0.71,0.36],
	'L':[-0.56,-1.12,-1.7,-2.54,3.62,3.62,3.62,3.62,3.62,0.03,0.02,0.01,0.01],
	'M':[-0.03,-0.05,-0.08,-0.12,-1.1,-1.1,-1.1,-1.1,-1.1,-0.26,-0.17,-0.11,-0.06],
	'N':[-0.25,-0.49,-0.74,-1.12,-2.36,-2.36,-2.36,-2.36,-2.36,0.44,0.29,0.19,0.1],
	'P':[-0.05,-0.1,-0.15,-0.22,-1.63,-1.63,-1.63,-1.63,-1.63,0.42,0.28,0.18,0.09],
	'Q':[0.37,0.74,1.13,1.69,-1.6,-1.6,-1.6,-1.6,-1.6,0.22,0.15,0.1,0.05],
	'R':[0.39,0.78,1.19,1.78,0.79,0.79,0.79,0.79,0.79,2.58,1.72,1.14,0.57],
	'S':[0.04,0.09,0.13,0.2,-1.27,-1.27,-1.27,-1.27,-1.27,0.34,0.23,0.15,0.07],
	'T':[0.3,0.6,0.91,1.36,-0.27,-0.27,-0.27,-0.27,-0.27,0.73,0.48,0.32,0.16],
	'V':[0.08,0.17,0.26,0.39,1.75,1.75,1.75,1.75,1.75,-1.05,-0.7,-0.46,-0.23],
	'W':[0.14,0.29,0.43,0.65,-3.49,-3.49,-3.49,-3.49,-3.49,-0.17,-0.12,-0.08,-0.04],
	'X':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'-':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'Y':[-0.06,-0.13,-0.19,-0.29,1.88,1.88,1.88,1.88,1.88,-1.73,-1.15,-0.76,-0.38]
}

# From http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000475
# van Durme et al. PLoS Computational Biology 2009
van_durme_hsp70_weight_matrix = {
	'A':[-1.60,-8.43,-0.67,-8.43,-0.37,0.11,-0.14],
	'C':[-6.13,-6.23,-5.81,-5.30,-4.75,-6.35,-6.26],
	'D':[0.05,-1.23,-8.71,-9.61,-1.24,-8.60,-1.30],
	'E':[1.15,-8.12,0.35,-9.61,-0.37,-2.24,-0.54],
	'F':[0.04,2.87,4.03,-0.11,2.78,3.36,2.03],
	'G':[-1.25,-2.06,-8.85,-9.15,-0.10,-8.87,0.43],
	'H':[0.22,0.08,0.51,-6.31,-5.78,0.19,0.72],
	'I':[0.54,2.80,0.47,4.60,0.63,1.59,-0.70],
	'K':[1.54,2.00,-6.92,0.16,0.48,1.25,0.25],
	'L':[0.37,3.00,4.51,5.67,2.37,2.12,-0.36],
	'M':[-0.17,-5.63,-4.55,5.22,2.68,1.41,0.44],
	'N':[0.58,0.30,-7.35,-7.15,-0.02,0.39,1.22],
	'P':[-7.76,-7.51,-6.70,-8.39,1.83,1.43,1.14],
	'Q':[1.47,-7.02,0.79,-5.60,0.50,1.02,-0.32],
	'R':[0.15,1.83,1.66,0.87,1.99,-0.08,2.34],
	'S':[0.34,-8.36,-8.61,-8.65,-0.86,-1.54,-0.33],
	'T':[0.47,-7.91,-0.49,1.46,0.22,0.88,0.72],
	'V':[-0.44,0.86,1.68,1.79,1.45,-0.45,-1.87],
	'W':[10.07,-2.05,2.89,-10.15,-2.19,10.45,-4.31],
	'Y':[1.46,3.39,5.08,-12.11,3.42,2.63,1.29],
	'-':[0,0,0,0,0,0,0],
	'X':[0,0,0,0,0,0,0]
}

# For testing
flat_weight_matrix = {
	'A':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'C':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'D':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'E':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'F':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'G':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'H':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'I':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'K':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'L':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'M':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'N':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'P':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'Q':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'R':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'S':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'T':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'V':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'W':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'X':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'-':[0,0,0,0,0,0,0,0,0,0,0,0,0],
	'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0]
}

weight_matrices = {'rudiger':rudiger_hsp70_weight_matrix, 'vandurme':van_durme_hsp70_weight_matrix, 'flat':flat_weight_matrix}


def _scoreRegion(seq, matrix): #, exclude_from_left=0, exclude_from_right=0):
	"""Score a region of the protein. This sequence assumes that the input sequence is no longer than the window-size
	specified by the scoring matrix.
	
		@return a number, the sum of the scores.
	"""
	n = len(seq)
	assert(n == len(list(matrix.values())[0]))
	scores = [matrix[x][i] for (x,i) in zip([a for a in seq], range(n)) if x in matrix]
	return sum(scores)
	
def _scoreWindows(seq, matrix, return_windows=False):
	# Scan sequence
	window_size = len(list(matrix.values())[0])
	# Pad the sequence so that sliding window does not require any special-casing
	pad = '-'
	assert(matrix[pad][0] == 0) # Check to make sure padding does not alter the score.
	padded_seq = pad*window_size + seq + pad*window_size
	n = len(seq)
	res = ScoreResult(seq, matrix)
	window_size = res.window_size
	for i in range(len(seq)+window_size-1):
		window = padded_seq[(i+1):(i+window_size+1)]
		assert len(window) == window_size
		window_score = _scoreRegion(window, matrix)
		res.addScore(window_score)
		if return_windows:
			res.addWindow(window)
	return res

def score(seq, matrix, return_windows=False):
	return _scoreWindows(seq, matrix, return_windows)

########################
# TEST CASES
########################
	
class test_scoreRegion(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		s = _scoreRegion(seq, rudiger_hsp70_weight_matrix)
		self.assertTrue(abs(s+11.84) < 0.001)

class test_scoreRegion_randomAA(unittest.TestCase):
	def test_run(self):
		seq = 'DALLAN@ANDREI'
		# Don't throw an error if there's a character you don't know, just
		# score it as zero.
		s = score(seq, rudiger_hsp70_weight_matrix)

class test_scoreWindows(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		matrix = hsp70_weight_matrix
		window_size1 = len(list(matrix.values())[0])

		s = _scoreWindows(seq, matrix)
		self.assertTrue(s.window_size == window_size1)
		self.assertTrue(len(s) == len(seq)+window_size1-1)
		self.assertTrue(abs(s[window_size-1]-11.84) < 0.001)

class test_scoreWindows(unittest.TestCase):
	def test_run(self):
		seq = 'MNNAANTGTTNESNV'
		matrix = rudiger_hsp70_weight_matrix
		window_size = len(list(matrix.values())[0])

		score_res = score(seq, matrix, return_windows=True)
		self.assertTrue(score_res._windows[0]=='-'*(window_size-1)+seq[0])
		self.assertTrue(score_res._windows[-1]==seq[-1]+'-'*(window_size-1))
		#for (s,r) in zip(scores,regions):
		#	print s,r

class test_scoreResidues(unittest.TestCase):
	def test_run(self):
		seq = 'DALLANSANDREI'
		matrix = rudiger_hsp70_weight_matrix
		window_size = len(list(matrix.values())[0])
		res = score(seq, matrix)
		rs = [s for s in res.residue_scores]
		self.assertTrue(len(rs)==len(seq))

class test_scoreResidues_hsf(unittest.TestCase):
	def test_run(self):
		seq = 'MNNAANTGTTNESNV'
		matrix = rudiger_hsp70_weight_matrix
		window_size = len(list(matrix.values())[0])
		res = score(seq, matrix)
		rs = [s for s in res.residue_scores]
		self.assertTrue(len(rs)==len(seq))

class test_scoreWindows_individual(unittest.TestCase):
	def test_run(self):
		matrix = rudiger_hsp70_weight_matrix
		window_size = len(list(matrix.values())[0])
		for aa in 'ACDEFGHIKLMNPQRSTVWY':
			res = score(aa, matrix)
			self.assertTrue(res._scores[0] == matrix[aa][-1]) # score in first window should always be the score for the first residue in the last position.
			#print scores

if __name__=="__main__":
	# Run tests
	unittest.main(verbosity=3)
	
	