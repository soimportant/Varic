#! /usr/bin/bash

# set -x
set -e
source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio
ulimit -s unlimited

ROOT_DIR="/mnt/ec/ness/yolkee/thesis"

# read and parse args, if not set, use default values
while getopts ":a:t:p:d" opt; do
  case $opt in
    a)
      depth=$OPTARG
      ;;
    t)
      thread=$OPTARG
      ;;
    p)
      pruned=$OPTARG
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

if [ -z "$depth" ]; then
  depth=-1
fi

if [ -z "$pruned" ]; then
  pruned=0.95
fi

echo "Debug mode: $debug, max_depth = $depth, pruned = $pruned"

# massif_output="/mnt/ec/ness/yolkee/thesis/scripts/massif.txt"
#heap_profile="/mnt/ec/ness/yolkee/thesis/tmp/profile/heap"
#export HEAPPROFILE=$heap_profile

species="Ecoli"
strain="K12_v2"

DATA_DIR="$ROOT_DIR/data/$species/$strain"
RESULT_DIR="$ROOT_DIR/results/$species/$strain"

# iterate seq_platform("PacBio", "ONT") and seq_depth(10, 20, 30, 50)

# for seq_platform in "ONT" "ONT_HQ" "PacBio-RSII" "PacBio-SEQUEL"; do

for seq_depth in 10; do
  for seq_platform in "ONT"; do

    raw_read="$DATA_DIR/reads/$seq_platform/D$seq_depth/merged_reads.fastq"
    overlap_paf="$DATA_DIR/reads/$seq_platform/D$seq_depth/merged_reads.overlap.paf"
    OUTPUT_DIR="$RESULT_DIR/$seq_platform/D$seq_depth/me_$pruned"
    if [ $depth -ne -1 ]; then
      OUTPUT_DIR="$OUTPUT_DIR"_d"$depth"_early
    fi

    corrected_read="$OUTPUT_DIR/corrected.fasta"

    # create directory
    mkdir -p $OUTPUT_DIR
    echo "Running $seq_platform $seq_depth"

    # Run correction
    echo "Running correction"
    CORRECTOR_EXE="$ROOT_DIR/build/varic"
    # CORRECTOR_ARGS="-r $raw_read -c $overlap_paf -o $corrected_read -t $thread -p $seq_platform --match $match --mismatch $mismatch --gap $gap --extend $gap"
    CORRECTOR_ARGS="-r $raw_read -c $overlap_paf -o $corrected_read -t $thread -p $seq_platform -d $depth --prune $pruned"
    if [ $debug -eq 1 ]; then
      CORRECTOR_ARGS="$CORRECTOR_ARGS --debug"
    fi
    time $CORRECTOR_EXE $CORRECTOR_ARGS
    
    
    # Run evaluation
    echo "Running evaluation"
    EVAL_EXE="$ROOT_DIR/eval/scripts/run_eval.py"
    EVAL_ARGS="-d $DATA_DIR -c $corrected_read -p $seq_platform --depth $seq_depth -o $OUTPUT_DIR -t $thread"
    if [ $debug -eq 1 ]; then
      EVAL_ARGS="$EVAL_ARGS --debug"
    fi
    time $EVAL_EXE $EVAL_ARGS
    
  done
done



