#! /usr/bin/bash

set -e
source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio
ulimit -s unlimited


# argparse platform input output and thread

if [ $# -ne 4 ]; then
  echo "Usage: cal_overlap.sh <platform> <input> <output> <thread>"
  exit 1
fi

platform=$1
input=$2
output=$3
thread=$4

if [ "$platform" == "PacBio" ]; then
  minimap2 -x ava-pb --dual=yes $input $input -t $thread -o $output
elif [ "$platform" == "ONT" ]; then
  minimap2 -x ava-ont --dual=yes $input $input -t $thread -o $output
else
  echo "Invalid platform: $platform"
  exit 1
fi