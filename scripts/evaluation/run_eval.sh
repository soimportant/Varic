#! /bin/bash

set -e
set -x

root="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12"

exe="/mnt/ec/ness/yolkee/thesis/build/eval"
origin_read="$root/reads/ONT/merged_reads_D20.fastq"
correct_read="/mnt/ec/ness/yolkee/thesis/result/vechat/Ecoli/a.fasta"
plodiy=2

$exe -o           $origin_read \
     -n           $plodiy \
     -c           $correct_read \
     -r           "$root/snp_ref/h1.simseq.genome.fa" \
                  "$root/snp_ref/h2.simseq.genome.fa" \
     --maf        "$root/reads/ONT/h1_D20_0001.maf" \
                  "$root/reads/ONT/h2_D20_0001.maf" \
     --snp_vcf    "$root/snp_ref/h1.refseq2simseq.SNP.vcf" \
                  "$root/snp_ref/h2.refseq2simseq.SNP.vcf" \
     --indel_vcf  "$root/snp_ref/h1.refseq2simseq.INDEL.vcf" \
                  "$root/snp_ref/h2.refseq2simseq.INDEL.vcf"
