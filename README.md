hsp70-binding
=============

To get a full profile for a particular gene:
python hsp70-motif-finder.py --database scer_orf_trans_all.fasta --pssm rudiger --threshold 5.0 --id YAL005C --report --out ssa1-motifs.txt

(You can use --pssm vandurme to use the van Durme 2009 paper's scoring matrix.)

To get short profiles for all genes in a FASTA file:
python hsp70-motif-finder.py --database scer_orf_trans_all.fasta --pssm rudiger --report --threshold 5.0 --out scer-hsp70-motifs-summary.txt

To profile a sequence you just made up:
python hsp70-motif-finder.py --sequence RYLLKNRANSSTSS --pssm rudiger --threshold 5.0 --report

To mask an alignment you pass in, covering up all sites that are not scored as the center of an Hsp70 binding motif:
python hsp70-motif-finder.py --database hsf-all-manual.fa --pssm rudiger --mask --threshold 5.0 --out hsf-masked-hsp70-binding-sites.fa

To get some information about these options:
python hsp70-motif-finder.py --help
