#! /usr/bin/python3

# TODO: run small eval on single directory, not all platform and all depth

from subprocess import *
import os
import sys
import argparse
import logging

def check_dir(path):
  if not os.path.exists(path) or not os.path.isdir(path):
    print(f"{path} not exists or not a directory.")
    exit(-1)


def check_argument(args):
  data_dir = args.data_dir
  corrected_read = args.corrected_read
  snp_ref_dir = data_dir + "/snp_ref"
  reads_dir = f"{data_dir}/reads/{args.platform}/D{args.depth}"

  check_dir(data_dir)
  check_dir(snp_ref_dir)
  check_dir(reads_dir)
  for cread in corrected_read:
    if not os.path.exists(cread):
      print(f"Corrected reads {cread} not exists.")
      exit(-1)

  plodiy = len([s for s in os.listdir(snp_ref_dir) if s.endswith(".fa")])
  if plodiy == 0:
    print(f"No haplotype found in {snp_ref_dir}")
    exit(-1)
  
  if args.threads < 1:
    print(f"Invalid thread number {args.thread}")
    exit(-1)


def run_cmd(cmd : list):
  logging.info(f"Run command: {' '.join(cmd)}")
  p = run(cmd)
  if p.returncode != 0:
    logging.error(f"cmd = {' '.join(cmd)} failed!")
    exit(p.returncode)


def main():
  logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s %(levelname)s] %(message)s")
  parser = argparse.ArgumentParser(prog="eval")

  # should be a data dir
  parser.add_argument("-d", "--data_dir", 
                      dest="data_dir",
                      type=str, required=True,
                      help="data directory contains snp_ref and reads")
  # could be multiple corrected reads
  parser.add_argument("-c", "--corrected_read", 
                      dest="corrected_read",
                      type=str, required=True, nargs='+',
                      help="path of corrected reads")
  parser.add_argument("-p", "--platform",
                      dest="platform", type=str,
                      required=True,
                      choices=["PacBio", "ONT"],
                      help="Sequencing platform of raw reads")
  parser.add_argument("--depth",
                      dest="depth", type=int,
                      required=True,
                      help="Sequencing depth of raw reads")
  parser.add_argument("-t", "--threads",
                      dest="threads", type=int,
                      default=1,
                      help="Number of threads to use")
  args, unknown = parser.parse_known_args(sys.argv[1:])
  check_argument(args)

  data_dir = args.data_dir
  corrected_read = args.corrected_read
  snp_ref_dir = data_dir + "/snp_ref"
  reads_dir = f"{data_dir}/reads/{args.platform}/D{args.depth}"
  plodiy = len([s for s in os.listdir(snp_ref_dir) if s.endswith(".fa")])

  # run eval on single haplotype
  eval_exe = "/mnt/ec/ness/yolkee/thesis/build/eval2"
  
  # raw reads
  raw_reads = f"{reads_dir}/merged_reads.fastq"
  eval_args = [eval_exe, "-n", f"{plodiy}"]
  eval_args += ["-o", raw_reads]

  # corrected reads
  eval_args += ["-c"] + corrected_read

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
    eval_args += [f"{reads_dir}/h{p}_0001.maf"]

  # threads
  eval_args += ["-t", f"{args.threads}"]

  # platform
  eval_args += ["-p", args.platform]

  print(' '.join(eval_args))
  run_cmd(eval_args)

if __name__ == "__main__":
  main()
