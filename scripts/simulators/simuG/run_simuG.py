#! /usr/bin/python3

from subprocess import *
import os
import sys
import argparse
import logging

def run_cmd(cmd : list):
  logging.info(f"Run command: {' '.join(cmd)}")
  p = run(cmd, stdout=PIPE)
  if p.returncode != 0:
    logging.error(f"cmd = {' '.join(cmd)} failed!")
    exit(-1)


def main():
  logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s %(levelname)s] %(message)s")
  parser = argparse.ArgumentParser(prog="run_simuG.py",
                                   description="Simulate SNP genome by simuG and pbsim")
  parser.add_argument("-r", "--ref",
                      dest="ref", type=str,
                      required=True,
                      help="Path of reference genome that used to introduce SNPs")
  parser.add_argument("-c", "--haplotype",
                      dest="haplotype", type=int,
                      default=1,
                      help="Number of haplotypes, indicate how many runs for simuG")
  parser.add_argument("--snp_count",
                      dest="snp_count", type=int,
                      default=0,
                      help="Number of SNPs to introduce for simuG")
  parser.add_argument("--indel_count",
                      dest="indel_count", type=int,
                      default=0,
                      help="Number of indels to introduce for simuG")
  parser.add_argument("-o",
                      dest="output_dir", type=str,
                      required=True,
                      help="directory of output file")

  args = parser.parse_args(sys.argv[1:])


  # simuG paramter
  ref = args.ref
  haplotype = args.haplotype
  snp_count = args.snp_count
  indel_count = args.indel_count
  output_dir = args.output_dir
  
  simuG = f"/mnt/ec/ness/yolkee/thesis/tools/simuG/simuG.pl"

  if not os.path.exists(ref) or not os.path.isfile(ref):
    logging.error(f"Reference genome {ref} not exists or not a file.")
    exit(-1)
  if haplotype < 1:
    logging.error(f"Haplotype count must greater than 1.")
    exit(-1)
  if snp_count < 0:
    logging.error(f"SNP count must greater than 0.")
    exit(-1)
  if indel_count < 0:
    logging.error(f"Indel count must greater than 0.")
    exit(-1)

  # run simuG for each haplotype
  snp_ref_dir = f"{output_dir}/snp_ref"
  os.makedirs(snp_ref_dir, exist_ok=True)

  for h in range(haplotype):
    pref = f"{snp_ref_dir}/h{h + 1}"
    simuG_cmd = [
      "perl",
      simuG,
      "-r", f"{ref}",
      "-snp_count", f"{snp_count}",
      "-indel_count", f"{indel_count}",
      "-prefix", f"{pref}"
    ]
    run_cmd(simuG_cmd)
  
if __name__ == "__main__":
  main()
