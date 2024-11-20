#! python3

# TODO: run small eval on single directory, not all platform and all depth

from subprocess import *
import os
import sys
import argparse
import logging

# library can read fasta
from Bio import SeqIO


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


def run_cmd(cmd : list, exe_dir = None):
  cwd = os.getcwd()
  if exe_dir != None:
    if not os.path.exists(exe_dir):
      os.mkdir(exe_dir)
    os.chdir(exe_dir)
  logging.info(f"Run command: {' '.join(cmd)}")
  p = run(cmd)
  if p.returncode != 0:
    logging.error(f"cmd = {' '.join(cmd)} failed!")
    exit(p.returncode)
  os.chdir(cwd)  


def run_variant_eval(data_dir, corrected_read, threads, depth, platform, output):

  snp_ref_dir = data_dir + "/snp_ref"
  reads_dir = f"{data_dir}/reads/{platform}/D{depth}"
  plodiy = len([s for s in os.listdir(snp_ref_dir) if s.endswith(".fa")])

  # run eval on single haplotype
  eval_exe = "/mnt/ec/ness/yolkee/thesis/build/eval"
  
  raw_reads = f"{reads_dir}/merged_reads.fastq"
  eval_args = [eval_exe, "-n", f"{plodiy}"]

  # raw reads
  eval_args += ["-a", raw_reads]

  # corrected reads -> list
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
  eval_args += ["-t", f"{threads}"]

  # platform
  eval_args += ["-p", platform]

  # output_path
  eval_args += ["-o", f"{output}"]
  run_cmd(eval_args)


# run quast evaluation(reference-needed)
# since we only have haplotype
def run_quast_eval(data_dir, corrected_read, thread, depth, platform):

  snp_ref_dir = data_dir + "/snp_ref"
  plodiy = len([s for s in os.listdir(snp_ref_dir) if s.endswith(".fa")])

  for cread in corrected_read:
    read_dir = os.path.dirname(cread)
    for p in range(0, plodiy):
      reads = [r for r in SeqIO.parse(cread, "fasta") if r.id.startswith(f"{p}")]
      SeqIO.write(reads, f"{read_dir}/h{p + 1}.cread.fasta", "fasta")

    # TODO: split cread into to different haplotype
    for p in range(1, plodiy + 1):
      quast_cmd = ["quast.py"]
      quast_cmd += ["-r", f"{snp_ref_dir}/h{p}.simseq.genome.fa"]
      quast_cmd += ["-o", f"{os.path.dirname(cread)}/quast_h{p}"]
      quast_cmd += ["--threads", f"{thread}"]
      quast_cmd += ["--min-contig", "500"]
      quast_cmd += ["--ambiguity-usage", "one"]
      quast_cmd += ["--fast", f"{read_dir}/h{p}.cread.fasta"]
      run_cmd(quast_cmd)
    
    # ? can we combine two result from quast?
    
# run Merqury evaluation(reference-free)
def run_assemble_eval(corrected_read, threads, platform):
  for cread in corrected_read:
    logging.info(f"Run Merqury assemble evaluation on {cread}")
    cwd = os.path.dirname(cread)
    flye_cmd = ["flye", "--threads", f"{threads}", "--out-dir", cwd + "/flye"]
    if platform == "PacBio":
      flye_cmd += ["--pacbio-corr"]
    elif platform == "ONT":
      flye_cmd += ["--nano-corr"]
    flye_cmd += [cread]
    run_cmd(flye_cmd)
    meryl_cmd = f"meryl count k=16 {cread} output {cwd}/reads.meryl threads={threads}"
    run_cmd(meryl_cmd.split())
    merqury_cmd = f"merqury.sh {cwd}/reads.meryl {cwd}/flye/assembly.fasta merqury"
    run_cmd(merqury_cmd.split(), cwd + "/merqury")


# TODO: parse result from Quast and Merqury and our evaluation
def parse_result():
  pass


def draw_fig():
  pass



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
  parser.add_argument("-D", "--depth",
                      dest="depth", type=int,
                      required=True,
                      help="Sequencing depth of raw reads")
  parser.add_argument("-t", "--threads",
                      dest="threads", type=int,
                      default=1,
                      help="Number of threads to use")
  parser.add_argument("-o", "--output",
                      dest="output", type=str,
                      default="result",
                      help="Output directory to use")
  parser.add_argument("--tmp",
                      dest="tmp", type=str,
                      default="/tmp",
                      help="Temporary directory to use")
  args, unknown = parser.parse_known_args(sys.argv[1:])
  check_argument(args)

  data_dir = args.data_dir
  corrected_read = args.corrected_read
  threads = args.threads
  depth = args.depth
  platform = args.platform
  output = args.output

  snp_ref_dir = data_dir + "/snp_ref"
  reads_dir = f"{data_dir}/reads/{platform}/D{depth}"
  run_variant_eval(data_dir, corrected_read, threads, depth, platform, output)
  # run_assemble_eval(corrected_read, threads, platform)
  # run_quast_eval(data_dir, corrected_read, threads, depth, platform)

if __name__ == "__main__":
  main()
