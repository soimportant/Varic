#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse
import shutil

def get_genome_size(ref):
  with open(ref, "r") as f:
    for line in f:
      if line.startswith(">"):
        continue
      else:
        return len(line.strip())

def main():
  parser = argparse.ArgumentParser(prog="run_correction_tool")
  parser.add_argument("-r", dest="read", type=str, required=True, help="Path of raw reads")
  parser.add_argument("--ref", dest="ref", type=str, required=True, help="Path of reference genome of raw reads(Used by Canu)")
  parser.add_argument("-p", dest="platform", type=str, required=True, choices=["PacBio", "ONT", "HiFi"], help="Platform of raw reads")
  parser.add_argument("-o", dest="output", type=str, required=True, help="Path of output directory")
  parser.add_argument("-t", dest="threads", type=int, required=True, help="Number of threads")
  args, unknown = parser.parse_known_args(sys.argv[1:])

  read = os.path.abspath(args.read)
  ref = os.path.abspath(args.ref)
  platform = args.platform
  output = os.path.abspath(args.output)

  # check argument
  if not os.path.exists(read) or not os.path.isfile(read):
    print(f"Raw reads {read} not exists or not a file.")
    exit(-1)
  if not os.path.exists(ref) or not os.path.isfile(ref):
    print(f"Reference genome {ref} not exists or not a file.")
    exit(-1)
  if not os.path.exists(output):
    tool = os.path.dirname(output)
    os.makedirs(tool, exist_ok=True)

  print(f"Read Path = {read}")
  print(f"Reference Path = {ref}")
  print(f"Platform = {platform}")
  print(f"Output Path = {output}")
  print(f"Threads = {args.threads}")
  
  size = get_genome_size(ref)
  cwd = os.path.dirname(__file__)
  
  tmp = cwd + "/tmp"
  if not os.path.exists(tmp):
    os.mkdir(tmp)
  os.chdir(tmp)

  # TODO: add time and memory usage at output section
  for tool in os.listdir(cwd):
    if not os.path.isdir(cwd + "/" + tool) or tool == "tmp":
      continue
    cmd = [
      f"{cwd}/{tool}/run_{tool}.py",
      "-r", read,
      "--ref", ref,
      "-p", platform,
      "-o", f"{output}/{tool}/corrected.fasta",
      "-t", str(args.threads),
      "--size", str(size),
    ]
    print(f"Start run {tool}")
    print(' '.join(cmd))

    try:
      proc = Popen(cmd)
      proc.wait()
      if proc.returncode != 0:
        print(f"Error: {tool} failed with exit code {proc.returncode}")
    except:
      print(f"Error: {tool} failed")
    
    print(f"{tool} done", end="\n===========================================\n")
  
  # delete all file in current folder exclude this file and directory
  # shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
  main()
