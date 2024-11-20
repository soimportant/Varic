#! /bin/bash

source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio

ROOT_DIR="/mnt/ec/ness/yolkee/thesis"
DATA_DIR="$ROOT_DIR/data"
REF_DIR="$DATA_DIR/ref"

ref="$REF_DIR/human/GRCh38/chr21.fasta"
# ref="$REF_DIR/Ecoli/K12/GCF_000005845.2/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
# ref="$REF_DIR/yeast/S288C/GCF_000146045.2/ncbi_dataset/data/GCF_000146045.2/whole_genome.fasta"
# ref="$REF_DIR/Drosophila/ncbi_dataset/data/GCF_000001215.4/whole_genome.fasta"


haplotype=2
ratio=0.002
thread=40
min_read_length=0
method="QS"
# method="err"
output_dir="$DATA_DIR/human/GRCh38_$haplotype"_v2
# output_dir="$DATA_DIR/Ecoli/K12_"$haplotype"_v2"
# output_dir="$DATA_DIR/yeast/S288C_"$haplotype"_v3"
# output_dir="$DATA_DIR/"


python3 /mnt/ec/ness/yolkee/thesis/scripts/simulators/simulate_haplotype_and_read \
    -r $ref \
    -c $haplotype \
    --ratio $ratio \
    -t $thread \
    -m $method \
    -o $output_dir \
    -l $min_read_length
