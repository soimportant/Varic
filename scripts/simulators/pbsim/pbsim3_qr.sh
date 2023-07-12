#! /bin/bash
set -e

q:pbsim_dir=/mnt/ec/ness/yolkee/thesis/tools/pbsim3

if [ $# < 5 ]; then
  echo "Usage: $0 <ref_genome> <qs_model> <depth> <output>"
  exit 1
fi

ref_genome=$1
qs_model=$2
depth=$3
output=$4

exe=$pbsim_dir/src/pbsim

$exe  --strategy wgs \
      --method qshmm \
      --qshmm $err_model \
      --depth $depth \
      --genome $ref_genome \
      --prefix $output