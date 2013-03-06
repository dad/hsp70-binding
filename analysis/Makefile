PY = $$PYTHONEXE

DATA = ../data
GENOMES = ../../../bio/genomes
CUPR = ../../cupr/data/hsf-alignment

all: motif

motif:
	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --id YGL073W --report --out hsf1.out

sum:
#	$(PY) hsp70-motif-finder.py --database test.fa --report --out test.out
#	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --report --threshold -5.0 --out scer-hsp70-motifs-summary-5.txt
#	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --report --threshold -10.0 --out scer-hsp70-motifs-summary-10.txt
	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_trans.fasta --report --threshold -3.0 --out scer-hsp70-motifs-summary-3.txt
#	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_coding.fasta --translate --report --threshold 0.0 --out scer-hsp70-motifs-summary-0.txt
	
hsf:
	$(PY) hsp70-motif-finder.py --database $(CUPR)/hsf-all-manual.fa --mask --threshold -5.0 --out $(DATA)/hsf-hsp70-binding-sites.fa

seq:
	$(PY) hsp70-motif-finder.py --sequence RYLLKNRANSSTSS --threshold -10.0 --report
#	$(PY) hsp70-motif-finder.py --sequence RYLLKNRANDSTSS --report

test:
	nosetests --verbosity=3 motif.py