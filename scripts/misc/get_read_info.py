#%%
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys
import argparse
import itertools
import edlib

DATA_ROOT="/mnt/ec/ness/yolkee/thesis/data"
RESULT_ROOT="/mnt/ec/ness/yolkee/thesis/results"

def read_fasta(fasta_file):
  if not os.path.isfile(fasta_file):
    print(f"File {fasta_file} not found")
    sys.exit(1)
  return [record for record in SeqIO.parse(fasta_file, "fasta")]


def read_fastq(fastq_file):
  if not os.path.isfile(fastq_file):
    print(f"File {fastq_file} not found")
    sys.exit(1)
  return [record for record in SeqIO.parse(fastq_file, "fastq")]


def check_read_exists(reads, id):
  for read in reads:
    if read.id == id:
      return True
  return False


def main():
  parser = argparse.ArgumentParser(description="Compare the corrected reads with the perfect reads")
  parser.add_argument("-a", "--species", help="The species name")
  parser.add_argument("-b", "--strain", help="The strain name")
  parser.add_argument("-p", "--platform", help="The sequencing platform")
  parser.add_argument("-d", "--depth", help="The sequencing depth")
  parser.add_argument("-t", "--tool", help="The correction tool")
  parser.add_argument("-r", "--read", help="The read id")
  args = parser.parse_args()

  species = args.species
  strain = args.strain
  platform = args.platform
  depth = args.depth
  tool = args.tool
  read_id = args.read

  data_dir = f"{DATA_ROOT}/{species}/{strain}/reads/{platform}/D{depth}"
  result_dir = f"{RESULT_ROOT}/{species}/{strain}/{platform}/D{depth}/{tool}"

  raw_reads = read_fastq(f"{data_dir}/merged_reads.fastq")
  perfect_reads = read_fasta(f"{data_dir}/merged_perfect_reads.fasta")
  corrected_reads = read_fasta(f"{result_dir}/corrected.fasta")

  if check_read_exists(raw_reads, read_id) == False or check_read_exists(perfect_reads, read_id) == False or check_read_exists(corrected_reads, read_id) == False:
    print("Read id not found")
    sys.exit(1)

  tmp_dir = f"/tmp/{read_id}"
  os.makedirs(tmp_dir, exist_ok=True)

  raw_read = [read for read in raw_reads if read.id == read_id][0]
  perfect_read = [read for read in perfect_reads if read.id == read_id][0]
  corrected_read = [read for read in corrected_reads if read.id == read_id][0]

  with open(f"{tmp_dir}/raw_read.fasta", "w") as f:
    SeqIO.write(raw_read, f, "fasta")
  with open(f"{tmp_dir}/perfect_read.fasta", "w") as f:
    SeqIO.write(perfect_read, f, "fasta")
  with open(f"{tmp_dir}/corrected_read.fasta", "w") as f:
    SeqIO.write(corrected_read, f, "fasta")

  seqs = {
    "raw": raw_read.seq,
    "perfect": perfect_read.seq,
    "corrected": corrected_read.seq
  }
  # take permutation of seqs

  os.makedirs(f"{tmp_dir}/align", exist_ok=True)
  for perm in itertools.permutations(seqs.keys(), 2):
    seq1 = seqs[perm[0]]
    seq2 = seqs[perm[1]]
    result = edlib.align(seq1, seq2, mode="NW", task="path")
    nice = edlib.getNiceAlignment(result, seq1, seq2, show_end=True)
    file = f"{tmp_dir}/align/{perm[0]}_{perm[1]}.txt"
    with open(file, "w") as f:
      f.write("\n".join(nice.values()))


if __name__ == "__main__":
  main()
