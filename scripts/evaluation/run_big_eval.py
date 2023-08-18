#! /usr/bin/python3

# run big eval on all platform and all depth

from subprocess import *
import os
import sys
import argparse
import logging

def check_dir(path):
  if not os.path.exists(path) or not os.path.isdir(path):
    print(f"{path} not exists or not a directory.")
    exit(-1)

def run_cmd(cmd : list):
  logging.info(f"Run command: {' '.join(cmd)}")
  p = run(cmd, stdout=PIPE)
  if p.returncode != 0:
    logging.error(f"cmd = {' '.join(cmd)} failed!")
    exit(p.returncode)

def main():
  logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s %(levelname)s] %(message)s")
  parser = argparse.ArgumentParser(prog="eval")
  parser.add_argument("-d",
                      dest="data_dir", type=str,
                      required=True,
                      help="directory containing reads and haplotype")
  parser.add_argument("-r", "--read",
                      dest="res_dir", type=str,
                      required=True,
                      help="directory containing corrected reads")
  args, unknown = parser.parse_known_args(sys.argv[1:])
  
  data_dir = args.data_dir
  res_dir = args.res_dir
  snp_ref_dir = data_dir + "/snp_ref"
  reads_dir = data_dir + "/reads"

  check_dir(data_dir)
  check_dir(res_dir)
  check_dir(snp_ref_dir)
  check_dir(reads_dir)

  plodiy = len([s for s in os.listdir(snp_ref_dir) if s.endswith(".fa")])
  if plodiy == 0:
    print(f"No haplotype found in {snp_ref_dir}")
    exit(-1)
  
  eval_exe = "/mnt/ec/ness/yolkee/thesis/build/eval"
  eval2_exe = "/mnt/ec/ness/yolkee/thesis/build/eval2"
  eval2_args = [eval2_exe]
  
  raw_reads = f"{reads_dir}/merged_reads.fastq"

  for platform in ["ONT", "PacBio"]:
    for depth in [10, 20, 30, 50, 100]:
      eval_args = [eval_exe, "-n", f"{plodiy}"]
      raw_reads = f"{reads_dir}/{platform}/D{depth}/merged_reads.fastq"
      # raw reads
      eval_args += ["-o", raw_reads]
      res = f"{res_dir}/{platform}/D{depth}"

      # corrected reads
      eval_args += ["-c"]
      for corrected_reads in os.listdir(res):
        if corrected_reads.endswith(".fasta") or corrected_reads.endswith(".fa"):
          eval_args += [f"{res}/{corrected_reads}"]
      
      # haplotype (reference genome)
      eval_args += ["-r"]
      for p in range(1, plodiy + 1):        
        eval_args += [f"{snp_ref_dir}/h{p}.simseq.genome.fa"]

      # vcf contains snp information
      eval_args += ["--snp_vcf"]
      for p in range(1, plodiy + 1):
        eval_args += [f"{snp_ref_dir}/h{p}.refseq2simseq.SNP.vcf"]

      # vcf contains indel information
      eval_args += ["--indel_vcf"]
      for p in range(1, plodiy + 1):
        eval_args += [f"{snp_ref_dir}/h{p}.refseq2simseq.INDEL.vcf"]

      # maf contains the mapping information of raw reads to haplotype 
      eval_args += ["--maf"]
      for p in range(1, plodiy + 1):        
        eval_args += [f"{reads_dir}/h{p}.simseq.genome.fa"]
        
      print(' '.join(eval_args))  

if __name__ == "__main__":
  main()
