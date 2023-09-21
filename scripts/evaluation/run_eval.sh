#! /bin/bash

set -e

source "/mnt/ec/ness/yolkee/miniconda3/bin/activate" bio

exe="/mnt/ec/ness/yolkee/thesis/scripts/evaluation/run_eval.py"
 # ./run_eval.py -d ~/thesis/data/Ecoli/K12 -c ~/thesis/result/Ecoli/K12/ONT/D50/racon.fasta -p ONT --depth 50 -t 8
platform="ONT"
depth=20
tool="vechat"
threads=40

data_root="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12"
result_root="/mnt/ec/ness/yolkee/thesis/result/Ecoli/K12/$platform/D$depth/$tool"  
correct_read="$result_root/corrected.fasta"

python3 $exe -d           $data_root \
             -c           $correct_read \
             -p           $platform \
             --depth      $depth \
             -o           $result_root \
             -t           $threads


# root="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K13"
# exe="/mnt/ec/ness/yolkee/thesis/fbuild/eval"
# origin_read="$root/reads/ONT/merged_reads_D20.fastq"
# correct_read="/mnt/ec/ness/yolkee/thesis/result/Ecoli/K13/canu.fasta"
# plodiy=2

# $exe -o           $origin_read \
#      -n           $plodiy \
#      -c           $correct_read \
#      -r           "$root/snp_ref/h1.simseq.genome.fa" \
#                   "$root/snp_ref/h2.simseq.genome.fa" \
#      --maf        "$root/reads/ONT/h1_D20_0001.maf" \
#                   "$root/reads/ONT/h2_D20_0001.maf" \
#      --snp_vcf    "$root/snp_ref/h1.refseq2simseq.SNP.vcf" \
#                   "$root/snp_ref/h2.refseq2simseq.SNP.vcf" \
#      --indel_vcf  "$root/snp_ref/h1.refseq2simseq.INDEL.vcf" \
#                   "$root/snp_ref/h2.refseq2simseq.INDEL.vcf"
