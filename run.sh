# flye/tests/data/ecoli_500kb_reads.fastq.gz
reads_file=flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# reads_file=flye/tests/data/ecoli_500kb.fasta
out_dir=o
bin/flye --pacbio-corr $reads_file -g 500k -o $out_dir -t 8 -m 1000
