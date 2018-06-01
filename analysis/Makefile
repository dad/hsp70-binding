PY = python3

DROPBOX = "/Users/dad/Dropbox (Drummond Lab)"
DATA = ../data
GENOMES = $(DROPBOX)/bio/genomes
CUPR = ../../cupr/data

all: motif

motif:
	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --pssm rudiger --threshold 10.0 --id YGL073W --report --out $(DATA)/hsf1-motifs-r10.txt
	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --pssm rudiger --threshold 10.0 --id YER165W --report --out $(DATA)/pab1-motifs-r10.txt
#	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --pssm vandurme --threshold 10.0 --id YGL073W --report --out $(DATA)/hsf1-motifs-vd5.txt
	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --pssm rudiger --threshold 10.0 --id YAL005C --report --out $(DATA)/ssa1-motifs-r10.txt
	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --pssm rudiger --threshold 10.0 --id YER006W --report --out $(DATA)/nug1-motifs-r10.txt

sum:
#	$(PY) motif-finder.py --database test.fa --threshold 0.0 --pssm vandurme --report --out test.out
#	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_coding.fasta --translate --report --threshold 0.0 --out $(DATA)/scer-motifs-summary-0.txt
#	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --report --threshold 3.0 --out $(DATA)/scer-motifs-summary-3.txt
	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --pssm rudiger --report --threshold 12.15 --out $(DATA)/scer-motifs-summary-r12_15.txt
#	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --pssm vandurme --report --threshold 0.0 --out $(DATA)/scer-motifs-summary-vd0.txt
#	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --pssm vandurme --report --threshold 5.0 --out $(DATA)/scer-motifs-summary-vd5.txt
#	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --report --threshold 7.0 --out $(DATA)/scer-motifs-summary-7.txt
#	$(PY) motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --report --threshold 10.0 --out $(DATA)/scer-motifs-summary-10.txt
	
mask:
	$(PY) motif-finder.py --database $(CUPR)/hsf-alignment/hsf-all-manual.fa --pssm rudiger --mask --threshold 5.0 --out $(DATA)/hsf-hsp70-binding-sites-r5.fa
#	$(PY) motif-finder.py --database $(CUPR)/hsf-alignment/hsf-all-manual.fa --pssm vandurme --mask --threshold 0.0 --out $(DATA)/hsf-hsp70-binding-sites-vd0.fa
#	$(PY) motif-finder.py --database $(CUPR)/hsf-alignment/hsf-all-manual.fa --pssm vandurme --mask --threshold 2.0 --out $(DATA)/hsf-hsp70-binding-sites-vd2.fa
#	$(PY) motif-finder.py --database $(CUPR)/hsp70-alignment.fa --pssm rudiger --mask --threshold 5.0 --out $(DATA)/hsp70-hsp70-binding-sites.fa	
#	$(PY) motif-finder.py --database $(CUPR)/hsf-alignment/sigma32-alignment.fa --pssm rudiger --mask --threshold 5.0 --out $(DATA)/sigma32-hsp70-binding-sites.fa	

seq:
	$(PY) motif-finder.py --sequence RYLLKNRANSSTSS --pssm vandurme --threshold 0.0 --report
	$(PY) motif-finder.py --sequence RYLLKNRANSSTSS --pssm rudiger --threshold 5.0 --report
#	$(PY) motif-finder.py --sequence RYLLKNRANDSTSS --report

test:
	nosetests --verbosity=3 motif_test.py