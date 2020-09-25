import sys, os, math, argparse, random
import stats, util, biofile, na, translate
import motif
import scipy as sp

# Goal:
#	Read in FASTA database and energy threshold
# 	Scan protein sequence and 
#	Write out, for each protein:
# 1) For each length n = 1..., number of regions in protein with n contiguous sites below threshold


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

# DAD: need more comments here!

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Detection of motifs using position-specific scoring matrices")
	parser.add_argument("-d", "--database", dest="fasta_fname", default=None, help="filename of FASTA database containing proteins to score")
	parser.add_argument("--id", dest="protein_id", default=None, help="particular protein identifier to score")
	parser.add_argument("-s", "--sequence", dest="sequence", default=None, help="particular protein sequence to score")
	parser.add_argument("--pssm", dest="pssm_name", default='rudiger', help="position-specific scoring matrix name [rudiger, vandurme, flat]")
	parser.add_argument("--translate", dest="translate",action="store_true", default=False, help="translate incoming sequences?")
	parser.add_argument("-t", "--threshold", type=float, dest="score_threshold", default=5.0, help="score threshold defining a putative binding site")
	parser.add_argument("--max-frequency-bin", dest="maximum_frequency_bin", default=10, help="maximum number of sequential binding sites to count")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove gaps from input FASTA file?")
	parser.add_argument("--report", dest="write_report", action="store_true", default=False, help="write out specific report for each protein?")
	parser.add_argument("--distribution", dest="write_distribution", action="store_true", default=False, help="write out complete distribution of scores for all windows?")
	parser.add_argument("--randomize", dest="randomize", action="store_true", default=False, help="shuffle sequences before calculation?")
	parser.add_argument("-m", "--mask", dest="mask_sequences", action="store_true", default=False, help="mask input sequences?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output (summary) filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = open(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in list(optdict.items()):
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	orf_dict = None
	gene_orf_map = None
	if not options.fasta_fname is None:
		fname = os.path.expanduser(options.fasta_fname)
		(headers, sequences) = biofile.readFASTA(fname)
		if options.randomize:
			shuffled_seqs = []
			for seq in sequences:
				shuf = [x for x in seq]
				random.shuffle(shuf)
				shuffled_seqs.append(''.join(shuf))
			sequences = shuffled_seqs
		orf_dict = dict(zip([biofile.firstField(h) for h in headers], sequences))
		gene_orf_map = dict([(biofile.secondField(h), biofile.firstField(h)) for h in headers])

	# Set the weight matrix
	try:
		matrix = motif.weight_matrices[options.pssm_name]
	except KeyError as ke:
		outs.write("# Unable to find weight matrix {}; try one of {}\n".format(options.pssm_name, ','.join(motif.weight_matrices.keys())))
	window_size = len(matrix['A']) #len(matrix.values()[0])
	# for associating windows with residues, center them
	mid_window = int(math.floor(window_size/2.0))
	
	seq = None
	if not orf_dict is None and not options.protein_id is None:
		try:
			seq = orf_dict[options.protein_id]
			data_outs.write("# Found ORF {}\n".format(options.protein_id))
		except KeyError as ke:
			# Maybe the user passed in a gene name.
			try:
				orf = gene_orf_map[options.protein_id]
				data_outs.write("# Found gene {} = {}\n".format(options.protein_id, orf))
				seq = orf_dict[orf]
			except KeyError:
				raise KeyError("# No protein found for ID {}".format(options.protein_id))
	
	if not options.sequence is None:
		seq = options.sequence
	
	if options.degap and not seq is None:
		seq = seq.replace('-','')

	if options.write_report and not seq is None:
		score_res = motif.score(seq, matrix, return_windows=True)
		header = 'pos\tresidue\tscore\tabove.threshold\twindow\n'
		data_outs.write(header)
		for sentry in score_res.results:
			resthresh_out = ' '
			if sentry.score >= options.score_threshold and sentry.residue != '-':
				resthresh_out = "*"
			line = "{pos}\t{aa}\t{resscore}\t{resthresh}\t{win}\n".format(
				pos=sentry.pos+1, aa=sentry.residue, resscore=na.formatNA(sentry.score,"{:1.2f}"), resthresh=resthresh_out, win=sentry.window)
			data_outs.write(line)
		sys.exit()

	# Write output
	dout_dist = util.DelimitedOutput()
	dout_dist.addHeader('orf','systematic ORF identifier','s')
	dout_dist.addHeader('pos','position (1-based) of window','d')
	dout_dist.addHeader('score','score of windowed sequence','f')

	n_written = 0
	if options.write_distribution:
		# Write the header descriptions
		dout_dist.describeHeader(data_outs)
		# Write the header fields
		dout_dist.writeHeader(data_outs)

		# Write entire distribution of scores
		for (hdr,rawseq) in zip(headers,sequences):
			orf = biofile.firstField(hdr)
			if options.translate:
				seq = translate.translate(rawseq)
				if seq is None:
					outs.write("# Skipping {} -- bad translation\n".format(orf))
					continue
			else:
				seq = rawseq
			# Remove trailing '*' (stop codon) if present
			if seq[-1] == '*':
				seq = seq[0:-1]
			score_res = motif.score(seq, matrix)
			L = len(seq)

			result = dout_dist.createResult(default=None)
			for (aai, score) in enumerate(score_res.scores):
				if aai>=window_size and aai<L-window_size: # Select for complete windows only.
					result['pos'] = aai+1 # 1-based indexing
					result['score'] = score
					result['orf'] = orf
					line = dout_dist.formatLine(result)
					data_outs.write(line)
					n_written += 1
		
	if options.write_report and not orf_dict is None:
		header = "orf\tnum.sites\tnum.motifs\tprop.sites\tmax.score\t" + \
			'\t'.join(['num.motifs.{:d}'.format(i+1) for i in range(options.maximum_frequency_bin)]) + '\tnum.motifs.longer\n'
		data_outs.write(header)
		for (hdr,rawseq) in zip(headers,sequences):
			orf = biofile.firstField(hdr)
			if options.translate:
				seq = translate.translate(rawseq)
				if seq is None:
					data_outs.write("# Skipping {} -- bad translation\n".format(orf))
					continue
			else:
				seq = rawseq
			# Remove trailing '*' (stop codon) if present
			if seq[-1] == '*':
				seq = seq[0:-1]
			score_res = motif.score(seq, matrix)
			score_summary = score_res.summary(options.score_threshold, options.maximum_frequency_bin)
			line = "{orf}\t{nsites}\t{nmotifs}\t{propsite:1.3f}\t{ms:1.3f}\t{freq}\t{nlonger}\n".format(
				orf = orf, ms = score_summary.max_score,
				nsites = score_summary.num_sites, nmotifs = score_summary.num_motifs, propsite=score_summary.num_sites/float(len(seq)),
				freq = '\t'.join(["{:d}".format(score_summary.run_frequency[i]) for i in range(options.maximum_frequency_bin)]),
				nlonger = score_summary.num_longer_runs
				)
			data_outs.write(line)
			
	if options.mask_sequences:
		# To accept gapped alignments;
		# Make shadow alignment
		# Detect sites
		# Restore alignment, masking all but binding sites
		for (hdr, rawseq) in zip(headers,sequences):
			if options.degap:
				seq = rawseq.replace('-','')
			if options.translate:
				seq = translate.translate(seq)
				if seq is None:
					data_outs.write("# Skipping due to bad translation: {}\n".format(hdr))
					continue
			score_res = motif.score(seq, matrix)
			masked_seq = score_res.maskSequence(options.score_threshold, mask_char='x')
			realigned_seq = realignSequence(masked_seq, rawseq)
			line = ">{}\n{}\n".format(hdr, realigned_seq)
			data_outs.write(line)
			
	

