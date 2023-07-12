#! /bin/bash

source /mnt/ec/ness/yolkee/miniconda3/bin/activate base

seq_file="/mnt/ec/mammoth/yolkee/thesis/data/Ecoli/K12/reads/ONT/merged_reads_D20.fastq"
thread=40
output_file="/mnt/ec/mammoth/yolkee/thesis/result/vechat/Ecoli/a.fastq"

/mnt/ec/ness/yolkee/thesis/tools/vechat/scripts/vechat \
  $seq_file \
  -t $thread \
  -o $output_file
