#! /bin/bash

set -e

# Ecoli K12 -> 4641652 bp
# variant count = 4000 -> 0.08%

variant_count=4000
snp_rate=0.8
indel_rate=0.2

snp_count= `expr $variant_count \* $snp_rate`
indel_count= `expr $variant_count \* $indel_rate`

output_prefix="/mnt/ec/ness/yolkee/thesis/data/snp_ref/Ecoli/K12/"
ref="/mnt/ec/ness/yolkee/thesis/data/ref/Ecoli/K12/GCF_000005845.2/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"

simuG="/mnt/ec/ness/yolkee/thesis/tools/simuG/simuG.pl"

if [ ! -f $ref ]; then
    echo "Reference file not found: $ref"
    exit 1
fi

perl $simuG -r $ref \
            -snp_count $snp_count \
            -indel_count $indel_count \
            -prefix $output_prefix
