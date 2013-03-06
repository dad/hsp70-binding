import sys, os, math, string, argparse, re
import stats, util, biofile, na, translate
import motif
import scipy as sp

# Goal:
#	Read in FASTA database and energy threshold
# 	Scan protein sequence and 
#	Write out, for each protein:
# 1) For each length n = 1..., number of regions in protein with n contiguous sites below threshold

def maskSequence(seq, mask, mask_char):
	masked_seq = ''
	for i in range(len(seq)):
		if mask[i]:
			masked_seq += mask_char
		else:
			masked_seq += seq[i]
	return masked_seq

def realignSequence(seq, aligned_seq, gap='-'):
	alseq_pos = 0
	seq_pos = 0
	realigned_seq = ''
	while alseq_pos < len(aligned_seq):
		aa = aligned_seq[alseq_pos]
		if aa == gap:
			realigned_seq += gap
		else:
			realigned_seq += seq[seq_pos]
			seq_pos += 1
		alseq_pos += 1
	assert len(realigned_seq) == len(aligned_seq)
	return realigned_seq

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Extraction of evidence from MaxQuant evidence files")
	parser.add_argument("--d", "--database", dest="fasta_fname", default=None, help="filename of FASTA database containing proteins to score")
	parser.add_argument("--id", dest="protein_id", default=None, help="particular protein identifier to score")
	parser.add_argument("--s", "--sequence", dest="sequence", default=None, help="particular protein sequence to score")
	parser.add_argument("--translate", dest="translate",action="store_true", default=False, help="translate incoming sequences?")
	parser.add_argument("-t", "--threshold", type=float, dest="score_threshold", default=-5.0, help="score threshold defining a putative binding site")
	parser.add_argument("--max-frequency-bin", dest="maximum_frequency_bin", default=10, help="maximum number of sequential binding sites to count")
	parser.add_argument("-r", "--report", dest="write_report", action="store_true", default=False, help="write out specific report for each protein?")
	parser.add_argument("-m", "--mask", dest="mask_sequences", action="store_true", default=False, help="mask input sequences?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output (summary) filename")
	options = parser.parse_args()
	
	# Set up some output
	info_outs = util.OutStreams(sys.stdout)
	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)

	orf_dict = None
	if not options.fasta_fname is None:
		fname = os.path.expanduser(options.fasta_fname)
		(headers, sequences) = biofile.readFASTA(fname)
		orf_dict = zip([biofile.firstField(h) for h in headers], sequences)
	
	# Set the weight matrix
	matrix = motif.hsp70_weight_matrix
	window_size = len(matrix.values()[0])
	# for associating windows with residues, center them
	mid_window = int(math.floor(window_size/2.0))
	
	# Write out
	for (k,v) in vars(options).items():
		outs.write("# {} = {}\n".format(k,v))
	
	seq = None
	if not options.protein_id is None and not orf_dict is None: # analyze specific protein pulled from DB
		header = 'pos\taa\taa.score\taa.below.threshold\tregion\twindow.score\twindow.below.threshold\n'
		outs.write(header)
		try:
			seq = orf_dict[options.protein_id]
		except KeyError, ke:
			raise KeyError, "# No protein found for ID {}".format(options.protein_id)
	
	if not options.sequence is None:
		seq = options.sequence
	
	if options.write_report and not seq is None:
		score_res = motif.scoreWindows(seq, matrix, return_regions=True)
		residue_scores = motif.scoreResidues(seq, score_res, window_size) #, min)
		for pos in range(len(seq)+window_size-1):
			aa = '-'
			resscore_out = 'NA'
			resthresh_out = ' '
			aai = pos-mid_window
			if aai >= 0 and aai < len(seq):
				aa = seq[aai]
				resscore_out = residue_scores[aai]
				if resscore_out <= options.score_threshold:
					resthresh_out = "*"
			# Window results
			winthresh_out = ' '
			window = score_res.regions[pos]
			winscore_out = score_res.scores[pos]
			if winscore_out <= options.score_threshold:
				winthresh_out = "*"
			line = "{pos}\t{aa}\t{resscore}\t{resthresh}\t{win}\t{winscore}\t{winthresh}\n".format(
				#os=pos+1, 
				pos=aai+1,
				aa=aa, resscore=na.formatNA(resscore_out,"{:1.2f}"), resthresh=resthresh_out,
				win=window, winscore=na.formatNA(winscore_out,"{:1.2f}"), winthresh=winthresh_out)
			outs.write(line)
		sys.exit()
		
	if options.write_report and not orf_dict is None:
		header = "orf\tnum.sites\tnum.motifs\tprop.sites\tmin.score\t" + \
			'\t'.join(['num.motifs.{:d}'.format(i+1) for i in range(options.maximum_frequency_bin)]) + '\tnum.motifs.longer\n'
		outs.write(header)
		for (hdr,rawseq) in zip(headers,sequences):
			orf = biofile.firstField(hdr)
			if options.translate:
				seq = translate.translate(rawseq)
				if seq is None:
					outs.write("# Skipping {} -- bad translation\n".format(orf))
					continue
			else:
				seq = rawseq
			if seq[-1] == '*':
				seq = seq[0:-1]
			residue_scores = motif.getResidueScores(seq, matrix)
			score_summary = motif.summarizeScores(residue_scores, options.score_threshold, options.maximum_frequency_bin)
			line = "{orf}\t{nsites}\t{nmotifs}\t{propsite}\t{ms}\t{freq}\t{nlonger}\n".format(
				orf = orf, ms = score_summary.min_score,
				nsites = score_summary.num_sites, nmotifs = score_summary.num_motifs, propsite=score_summary.num_sites/float(len(seq)),
				freq = '\t'.join(["{:d}".format(score_summary.run_frequency[i]) for i in range(options.maximum_frequency_bin)]),
				nlonger = score_summary.num_longer_runs
				)
			outs.write(line)
			
	if options.mask_sequences:
		# To accept gapped alignments;
		# Make shadow alignment
		# Detect sites
		# Restore alignment, masking all but binding sites
		for (hdr, rawseq) in zip(headers,sequences):
			seq = rawseq.replace('-','')
			if options.translate:
				seq = translate.translate(seq)
				if seq is None:
					outs.write("# Skipping due to bad translation: {}\n".format(hdr))
					continue
			residue_scores = motif.getResidueScores(seq, matrix)
			masked_seq = maskSequence(seq, [x>options.score_threshold for x in residue_scores], mask_char='x')
			realigned_seq = realignSequence(masked_seq, rawseq)
			line = ">{}\n{}\n".format(hdr, realigned_seq)
			outs.write(line)
			
	

