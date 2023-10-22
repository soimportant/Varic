#! /usr/bin/bash

set -e

source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio

exe="/mnt/ec/ness/yolkee/thesis/build/main"
platform="ONT"
depth=20
raw_read="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/reads/$platform/D$depth/merged_reads.fastq"
overlap_paf="/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/reads/$platform/D$depth/merged_reads.overlap.paf"
output="/tmp/a.fasta"
thread=64

massif_output="/mnt/ec/ness/yolkee/thesis/scripts/massif.txt"

$exe -r $raw_read -c $overlap_paf -o $output -t $thread -p $platform --debug
