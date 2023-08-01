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
  parser.add_argument("-d",
                      dest="data_dir", type=str,
                      required=True,
                      help="directory containing reads and haplotype")
  parser.add_argument("-r", "--read",
                      dest="res_dir", type=str,
                      required=True,
                      help="directory containing corrected reads")  

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

  for tool in os.listdir(cwd):
    if not os.path.isdir(cwd + "/" + tool) or tool == "tmp":
      continue
    if tool != "lorma":
      continue
    cmd = [
      f"{cwd}/{tool}/run_{tool}.py",
      "-r", read,
      "--ref", ref,
      "-p", platform,
      "-o", f"{output}/{tool}.fasta",
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
  shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
  main()