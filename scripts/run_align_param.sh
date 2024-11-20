#! /usr/bin/bash

set -e
source /mnt/ec/ness/yolkee/miniconda3/bin/activate bio
ulimit -s unlimited

ROOT_DIR="/mnt/ec/ness/yolkee/thesis"

# read and parse args, if not set, use default values
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

species="Ecoli"
strain="K12"
seq_platform="ONT"
seq_depth=20

DATA_DIR="$ROOT_DIR/data/$species/$strain"
RESULT_DIR="$ROOT_DIR/results/align_params_try/$species/$strain"

# Create result directory
mkdir -p $RESULT_DIR/$seq_platform/D$seq_depth/me

raw_read="$DATA_DIR/reads/$seq_platform/D$seq_depth/merged_reads.fastq"
overlap_paf="$DATA_DIR/reads/$seq_platform/D$seq_depth/merged_reads.overlap.paf"


# massif_output="/mnt/ec/ness/yolkee/thesis/scripts/massif.txt"
#heap_profile="/mnt/ec/ness/yolkee/thesis/tmp/profile/heap"
#export HEAPPROFILE=$heap_profile
# Run correction

# match from 3 to 10
# mismatch from 1 to 5
# gap from 1 to 5

for match in {3..9}; do
  for mismatch in {2..6}; do
    for gap in {2..6}; do
      
      output="$RESULT_DIR/$seq_platform/D$seq_depth/me/corrected.fasta"
      echo "Running correction on match = $match, mismatch = $mismatch, gap = $gap"
      echo "Output: $output"

      CORRECTOR_EXE="$ROOT_DIR/build/main2"
      CORRECTOR_ARGS="-r $raw_read -c $overlap_paf -o $output -t $thread -p $seq_platform --match $match --mismatch -$mismatch --gap -$gap --extend -$gap"
      if [ $debug -eq 1 ]; then
        CORRECTOR_ARGS="$CORRECTOR_ARGS --debug"
      fi
      $CORRECTOR_EXE $CORRECTOR_ARGS

      # Run evaluation
      echo "Running evaluation"
      EVAL_EXE="$ROOT_DIR/eval/scripts/run_eval.py"
      EVAL_ARGS="-d $DATA_DIR -c $output -p $seq_platform --depth $seq_depth -o $RESULT_DIR/$seq_platform/D$seq_depth/me -t $thread"
      if [ $debug -eq 1 ]; then
        EVAL_ARGS="$EVAL_ARGS --debug"
      fi

      $EVAL_EXE $EVAL_ARGS

      # Draw evaluation plots"
      echo "Drawing plots"
      PLOT_EXE="$ROOT_DIR/eval/scripts/draw/draw_eval.py"
      PLOT_ARGS="-r $RESULT_DIR -p $seq_platform -d $seq_depth"
      $PLOT_EXE $PLOT_ARGS

      
      GET_EXE="$ROOT_DIR/eval/scripts/get_quality.py"
      GET_ARGS="-f $RESULT_DIR/$seq_platform/D$seq_depth/me/result.csv --match $match --mismatch $mismatch --gap $gap"
      python3 $GET_EXE $GET_ARGS
    done
  done
done
