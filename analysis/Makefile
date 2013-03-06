PY = $$PYTHONEXE

GENOMES = ../../../bio/genomes

all: motif

motif:
	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --id YGL073W --report --out hsf1.out

sum:
#	$(PY) hsp70-motif-finder.py --database test.fa --report --out test.out
#	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --report --threshold -5.0 --out scer-hsp70-motifs-summary-5.txt
	$(PY) hsp70-motif-finder.py --database $(GENOMES)/scer_orf_trans_all.fasta --report --threshold -10.0 --out scer-hsp70-motifs-summary-10.txt
	

seq:
	$(PY) hsp70-motif-finder.py --sequence RYLLKNRANSSTSS --threshold -5.0 --report
#	$(PY) hsp70-motif-finder.py --sequence RYLLKNRANDSTSS --report

test:
	nosetests --verbosity=3 motif.py