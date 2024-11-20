#! /bin/bash

set -e
source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio
ulimit -s unlimited

ROOT_DIR="/mnt/ec/ness/yolkee/thesis"

# read args and set default values
while getopts ":t:d" opt; do
  case $opt in
    t)
      thread=$OPTARG
      ;;
    d)
      debug=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

if [ -z "$thread" ]; then
  echo "Thread is not set"
  exit 1
fi

if [ -z "$debug" ]; then
  debug=0
fi

echo "Debug mode: $debug"

# TODO: change this in run_correction_tools.py
species="Ecoli"
strain="K12_4_v2"

DATA_ROOT="/mnt/ec/mammoth/yolkee/thesis/data"
RESULT_ROOT="/mnt/ec/mammoth/yolkee/thesis/results"

# iterate seq_platform("PacBio", "ONT") and seq_depth(20, 50, 100)
for seq_platform in "ONT"; do
  for seq_depth in 10; do 
    
    DATA_DIR=$DATA_ROOT/$species/$strain
    RESULT_DIR=$RESULT_ROOT/$species/$strain

    EXE="$ROOT_DIR/scripts/correction_tools/run_correction_tools.py"
    ARGS="-r $DATA_DIR -p $seq_platform -d $seq_depth -t $thread -o $RESULT_DIR"
    echo "Running $seq_platform $seq_depth"
    python3 $EXE $ARGS
  done
done
