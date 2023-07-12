#! /usr/bin/python3

from subprocess import *
import os
import sys
import argparse
import logging

"""
For haplotype 1 corrected reads, use following steps for evaluating miss
corrected variants between haplotype 1 and haplotype 2.

1. Align corrected reads to haplotype 1
  - bwa mem -t 8 -x ont2d -p haplotype_1.fa corrected_reads.fastq > corrected_reads_haplotype_1.sam
2. Convert sam coordinate to haplotype 2
  - python3 src/convert_sam_coordinate.py corrected_reads_haplotype_1.sam haplotype_1.fa haplotype_2.fa corrected_reads_haplotype_2.sam
3. For each variant on haplotype 2, check if it exists on alignment result
  - 
"""

def run_cmd(cmd : list):
  logging.info(f"Run command: {' '.join(cmd)}")
  p = run(cmd, stdout=PIPE)
  if p.returncode != 0:
    logging.error(f"cmd = {' '.join(cmd)} failed!")
    exit(p.returncode)

def main():
  logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s %(levelname)s] %(message)s")

  parser = argparse.ArgumentParser(prog="eval_between_haplotype",
                                    description="Evaluate wrongly corrected variants between haplotype 1 and haplotype 2")
  parser.add_argument("--h1", dest="h1", type=str, required=True, help="Path of haplotype 1")
  parser.add_argument("--h2", dest="h2", type=str, required=True, help="Path of haplotype 2")
  parser.add_argument("--reads", dest="reads", type=str, required=True, help="Path of corrected reads of haplotype 1")
  parser.add_argument("--vcf", dest="vcf", type=str, required=True, help="Path of vcf file of haplotype 2")
  parser.add_argument("-p", "--platform", dest="platform", type=str, required=True, choices=["PacBio", "ONT"], help="Platform of reads")

  args = parser.parse_args(sys.argv[1:])
  h1 = args.h1
  h2 = args.h2
  reads = args.read
  vcf = args.vcf

  # 1. Align corrected reads to haplotype 1
  # minimap2 command
  
  minimap2_cmd = [
    "minimap2",
    "-x", 

  ]


  simuG = f"/mnt/ec/ness/yolkee/thesis/tools/simuG/simuG.pl"
  pbsim_dir = "/mnt/ec/ness/yolkee/thesis/tools/pbsim3"
  pbsim = f"{pbsim_dir}/src/pbsim"

  # simuG paramter
  ref = args.ref
  haplotype = args.haplotype
  snp_count = args.snp_count
  indel_count = args.indel_count

  # pbsim parameter
  depth = args.depth
  output_dir = args.output_dir
  model_dir = f"{pbsim_dir}/data"
  models = {
    "qs": {
      "PacBio": "QSHMM-RSII.model",
      "ONT": "QSHMM-ONT.model"
    },
    "err": {
      "PacBio": "ERRHMM-RSII.model",
      "ONT": "ERRHMM-ONT.model"
    }
  }

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
  if depth < 1:
    logging.error(f"Depth must greater than 1.")
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
  
  # run pbsim for each haplotype
  reads_dir = f"{output_dir}/reads"
  os.makedirs(reads_dir, exist_ok=True)
  for platform, model in models[args.method].items():
    method = f"{args.method}hmm"
    pref = f"{reads_dir}/{platform}"
    os.makedirs(pref, exist_ok=True)
    for h in range(haplotype):
      genome = f"{snp_ref_dir}/h{h + 1}.simseq.genome.fa"
      pbsim_cmd = [
        pbsim,
        "--strategy", "wgs",
        "--method", f"{method}",
        f"--{method}", f"{model_dir}/{model}",
        "--depth", f"{depth}",
        "--genome", f"{genome}",
        "--id-prefix", f"{h}",
        "--prefix", f"{pref}/h{h + 1}_D{depth}"
      ]
      run_cmd(pbsim_cmd)

  # use seqkit to merge reads
  for platform in models[args.method].keys():
    pref = f"{reads_dir}/{platform}"
    merged_reads = f"{pref}/merged_reads_D{depth}.fastq"
    concat = [
      "seqkit",
      "concat",
      "-f",
      "-o", merged_reads,
    ]
    concat += [f"{pref}/h{h + 1}_D{depth}_0001.fastq" for h in range(haplotype)]
    run_cmd(concat)
  
  # use minimap2 for calculating overlap on merged read
  for platform in models[args.method].keys():
    pref = f"{reads_dir}/{platform}"
    merged_reads = f"{pref}/merged_reads_D{depth}.fastq"
    overlap = f"{pref}/merged_reads_D{depth}.overlap.paf"
    if platform == "PacBio":
      platform_str = "ava-pb"
    elif platform == "ONT":
      platform_str = "ava-ont"
    minimap2_cmd = [
      "minimap2",
      "-x", f"{platform_str}",
      "-X",
      f"{merged_reads}",
      f"{merged_reads}",
      "-o", f"{overlap}"
    ]
    run_cmd(minimap2_cmd)
  
  
if __name__ == "__main__":
  main()

