# flye/tests/data/ecoli_500kb_reads.fastq.gz
# reads_file=flye/tests/data/ecoli_500kb.fasta
reads_file=flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
out_dir=o
rm -rf ${out_dir}
bin/flye \
	--pacbio-corr $reads_file \
	--debug -g 500k -o $out_dir -t 8 -m 1000
exit

flye_root=/home/goshng/all/polap/Flye
out_dir=${flye_root}/o
rm -rf ${out_dir}
mkdir -p \
	${out_dir}/00-assembly \
	${out_dir}/10-consensus \
	${out_dir}/20-repeat \
	${out_dir}/30-contigger \
	${out_dir}/40-polishing

# Stage: assembly
# input1: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# output: ${out_dir}/00-assembly/draft_assembly.fasta
#
bin/flye-modules assemble \
	--reads ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--out-asm ${out_dir}/00-assembly/draft_assembly.fasta \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--genome-size 500000 \
	--min-ovlp 1000

exit

flye_root=/home/goshng/all/polap/Flye
out_dir=${flye_root}/l
mkdir -p \
	${out_dir}/00-assembly \
	${out_dir}/10-consensus \
	${out_dir}/20-repeat \
	${out_dir}/30-contigger \
	${out_dir}/40-polishing

bin/flye-modules assemble \
	--reads ${flye_root}/l.fq \
	--out-asm ${out_dir}/00-assembly/draft_assembly.fasta \
	--config ${flye_root}/flye/config/bin_cfg/asm_raw_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--genome-size 500000 \
	--min-ovlp 1000

exit

rm -rf ${out_dir}
bin/flye \
	--pacbio-corr $reads_file \
	--debug -g 500k -o $out_dir -t 8 -m 1000
exit

flye_root=/home/goshng/all/polap/Flye
out_dir=${flye_root}/o
mkdir -p \
	${out_dir}/00-assembly \
	${out_dir}/10-consensus \
	${out_dir}/20-repeat \
	${out_dir}/30-contigger \
	${out_dir}/40-polishing

# Stage: assembly
# input1: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# output: ${out_dir}/00-assembly/draft_assembly.fasta
#
bin/flye-modules assemble \
	--reads ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--out-asm ${out_dir}/00-assembly/draft_assembly.fasta \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--genome-size 500000 \
	--min-ovlp 1000

# Stage: consensus 1
# input: ${out_dir}/00-assembly/draft_assembly.fasta
# output: ${out_dir}/10-consensus/minimap.bam
#
flye-minimap2 \
	${out_dir}/00-assembly/draft_assembly.fasta \
	${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz -x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	flye-samtools view -T ${out_dir}/00-assembly/draft_assembly.fasta -u - |
	samtools view -h - |                                # dflye:
	awk '{if ($1 ~ /^@/ || ($2 & 16) == 0) print $0}' | # dflye: Filter for forward-forward reads
	samtools view -bS - |                               # dflye: Convert filtered SAM back to BAM
	flye-samtools sort -T ${out_dir}/10-consensus/sort_241207_112403 -O bam -@ 4 -l 1 -m 1G -o \
		${out_dir}/10-consensus/minimap.bam

# Stage: consensus 2
# Compute consensus -> hard
# minimap.bam -> consensus.fasta
# input1: ${out_dir}/10-consensus/minimap.bam
# input1: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# output: ${out_dir}/10-consensus/consensus.fasta

# Stage: repeat
# input1: ${out_dir}/10-consensus/consensus.fasta
# output: ${out_dir}/20-repeat
#
flye-modules repeat \
	--disjointigs ${out_dir}/10-consensus/consensus.fasta \
	--reads ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--out-dir ${out_dir}/20-repeat \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--min-ovlp 1000

# Stage: contigger
# input1: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# input2: ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg
# input3: ${out_dir}/20-repeat/repeat_graph_edges.fasta
# input4: ${out_dir}/20-repeat/repeat_graph_dump
# input5: ${out_dir}/20-repeat/read_alignment_dump
# output: ${out_dir}/30-contigger
#
flye-modules contigger \
	--graph-edges ${out_dir}/20-repeat/repeat_graph_edges.fasta \
	--reads ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--out-dir ${out_dir}/30-contigger \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--repeat-graph ${out_dir}/20-repeat/repeat_graph_dump \
	--graph-aln ${out_dir}/20-repeat/read_alignment_dump \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--min-ovlp 1000

# Stage: polishing 1
# input1: ${out_dir}/30-contigger/contigs.fasta
# input2: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# output: ${out_dir}/40-polishing/minimap_1.bam
#
flye-minimap2 \
	${out_dir}/30-contigger/contigs.fasta \
	${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz -x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	flye-samtools view -T ${out_dir}/30-contigger/contigs.fasta -u - |
	flye-samtools sort -T ${out_dir}/40-polishing/sort_241207_112408 -O bam -@ 4 -l 1 -m 1G \
		-o ${out_dir}/40-polishing/minimap_1.bam

# Stage: polishing 2
# input1: ${out_dir}/40-polishing/filtered_contigs.fasta
# input1: ${out_dir}/30-contigger/graph_final.fasta
# output: ${out_dir}/40-polishing/edges_aln.bam
#
flye-minimap2 \
	${out_dir}/40-polishing/filtered_contigs.fasta \
	${out_dir}/30-contigger/graph_final.fasta -x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	flye-samtools view -T ${out_dir}/40-polishing/filtered_contigs.fasta -u - |
	flye-samtools sort -T ${out_dir}/40-polishing/sort_241207_112414 -O bam -@ 4 -l 1 -m 1G \
		-o ${out_dir}/40-polishing/edges_aln.bam

# Stage: finalize
#
