#! /bin/bash

source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio

exe="/mnt/ec/ness/yolkee/thesis/scripts/correction_tools/run_correction_tools.py"
read="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/reads/ONT/D10/merged_reads.fastq"
ref="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/snp_ref/h1.simseq.genome.fa"
platform="ONT"
out="/mnt/ec/ness/yolkee/thesis/result/Ecoli/K12/ONT/D10"
thread=80

$exe -r $read \
     --ref $ref \
     -p $platform \
     -o $out \
     -t $thread
