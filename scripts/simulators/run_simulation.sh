#! /bin/bash

source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio

# ref="/mnt/ec/ness/yolkee/thesis/data/ref/human/GRCh38/chr1.fasta"
ref="/mnt/ec/ness/yolkee/thesis/data/ref/Ecoli//K12/GCF_000005845.2/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
haplotype=2
sz=4641652
snp_count=$(( $sz/500 ))
indel_count=$(( $sz/5000 ))
thread=6
method=qs
output_dir="/mnt/ec/ness/yolkee/thesis/data/human/GRCh38/"

python3 /mnt/ec/ness/yolkee/thesis/scripts/simulators/simulate_haplotype_and_read \
    -r $ref \
    -c $haplotype \
    --snp_count $snp_count \
    --indel_count $indel_count \
    -t $thread \
    -m $method \
    -o $output_dir
