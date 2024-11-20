#! python3

from subprocess import *
from resource import *
import sys
import os
import argparse
import shutil

# empty for running all
wanted_tools = ['canu', 'racon', 'vechat']

def run_cmd(cmd):
  proc = Popen(cmd)
  proc.wait()
  if proc.returncode != 0:
    print(f"Error: {cmd} failed with exit code {proc.returncode}")
    raise Exception(f"Error: {cmd} failed with exit code {proc.returncode}")


def get_genome_size(ref):
  with open(ref, "r") as f:
    for line in f:
      if line.startswith(">"):
        continue
      else:
        return len(line.strip())


ROOT_DIR="/mnt/ec/ness/yolkee/thesis"
EVAL_ROOT_DIR="/mnt/ec/ness/yolkee/thesis/eval"

def main():
  parser = argparse.ArgumentParser(prog="run_correction_tool")
  parser.add_argument("-r", dest="data_dir", type=str, required=True, help="Root directory of the data")
  parser.add_argument("-p", dest="platform", type=str, required=True, help="Platform")
  parser.add_argument("-d", dest="depth", type=int, required=True, help="Depth")
  parser.add_argument("-t", dest="threads", type=int, required=True, help="Number of threads")
  parser.add_argument("-o", dest="result_dir", type=str, required=True, help="Output directory")
  args, unknown = parser.parse_known_args(sys.argv[1:])

  threads = args.threads  
  data_dir = args.data_dir
  platform = args.platform
  depth = args.depth
  result_dir = f"{args.result_dir}/{platform}/D{depth}"
  
  raw_reads = f"{data_dir}/reads/{platform}/D{depth}/merged_reads.fastq"
  overlap_paf = f"{data_dir}/reads/{platform}/D{depth}/merged_reads.overlap.paf"
  ref = f"{data_dir}/snp_ref/h1.simseq.genome.fa"

  if not os.path.exists(raw_reads) or not os.path.isfile(raw_reads):
    print(f"Raw reads {raw_reads} not exists or not a file.")
    exit(-1)
  if not os.path.exists(overlap_paf) or not os.path.isfile(overlap_paf):
    print(f"Overlap {overlap_paf} not exists or not a file.")
    exit(-1)
  if not os.path.exists(ref) or not os.path.isfile(ref):
    print(f"Reference genome {ref} not exists or not a file.")
    exit(-1)
  if not os.path.exists(result_dir):
    tool = os.path.dirname(result_dir)
    os.makedirs(tool, exist_ok=True)

  print(f"Read Path = {raw_reads}")
  print(f"Reference Path = {ref}")
  print(f"Platform = {platform}")
  print(f"Output Path = {result_dir}")
  print(f"Threads = {threads}")
  
  genome_size = get_genome_size(ref)
  cwd = os.path.dirname(__file__)

  for tool in os.listdir(cwd):
    if not os.path.isdir(cwd + "/" + tool):
      continue

    if len(wanted_tools) > 0 and tool not in wanted_tools:
      continue

    tool_tmp_dir = f"{cwd}/{tool}/tmp"
    shutil.rmtree(tool_tmp_dir, ignore_errors=True)
    if not os.path.exists(tool_tmp_dir):
      os.mkdir(tool_tmp_dir)
    os.chdir(tool_tmp_dir)

    corrected_reads = f"{result_dir}/{tool}/corrected.fasta"

    tool_cmd = [
      "python3",
      f"{cwd}/{tool}/run_{tool}.py",
      "-r", raw_reads,
      "--ref", ref,
      "-p", platform,
      "-o", corrected_reads,
      "-t", str(threads),
      "--size", str(genome_size),
    ]
    eval_cmd = [
      "python3",
      f"{EVAL_ROOT_DIR}/scripts/run_eval.py",
      "-d", data_dir,
      "-c", corrected_reads,
      "-p", platform,
      "--depth", str(args.depth),
      "-o", f"{result_dir}/{tool}",
      "-t", str(threads),
    ]
    plot_cmd = [
      "python3",
      f"{EVAL_ROOT_DIR}/scripts/draw/draw_eval.py",
      "-r", f"{result_dir}/{tool}",
      "-p", platform,
      "-d", str(args.depth),
    ]
    
    print("Start running", tool)
    try:
      # run_cmd(tool_cmd)
      run_cmd(eval_cmd)
    #  run_cmd(plot_cmd)
    except Exception as e:
      print(f"Error: {tool} failed")
      print(e)
    
    print(f"{tool} done", end="\n===========================================\n")

if __name__ == "__main__":
  main()
