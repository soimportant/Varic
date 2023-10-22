#! /bin/bash

source "/mnt/ec/ness/yolkee/miniconda3/bin/activate" bio

exe="/mnt/ec/ness/yolkee/thesis/scripts/correction_tools/run_correction_tools.py"
platform="ONT"
depth=20
threads=80
read="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/reads/$platform/D$depth/merged_reads.fastq"
ref="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/snp_ref/h1.simseq.genome.fa"
out="/mnt/ec/ness/yolkee/thesis/results/Ecoli/K12/$platform/D$depth"

$exe -r $read \
     --ref $ref \
     -p $platform \
     -o $out \
     -t $threads
