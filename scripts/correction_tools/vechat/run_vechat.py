#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse
import shutil
from Bio import SeqIO

TOOL="vechat"

def main():
  parser = argparse.ArgumentParser(prog="run_vechat")
  parser.add_argument("-r", dest="read", type=str, required=True, help="Path of raw reads")
  parser.add_argument("-p", dest="platform", type=str, required=True, choices=["PacBio", "ONT", "HiFi"], help="Platform of raw reads")
  parser.add_argument("-o", dest="output", type=str, required=True, help="Path of output file")
  parser.add_argument("-t", dest="threads", type=int, required=False, help="Number of threads")
  args, unknown = parser.parse_known_args(sys.argv[1:])

  read = os.path.abspath(args.read)
  platform = args.platform
  output = os.path.abspath(args.output)
  threads = args.threads

  # check argument
  if not os.path.exists(read) or not os.path.isfile(read):
    print(f"Raw reads {read} not exists or not a file.")
    exit(-1)
  if not os.path.exists(output):
    dir = os.path.dirname(output)
    os.makedirs(dir, exist_ok=True)

  if platform == "PacBio":
    platform = "pb"
  elif platform == "ONT":
    platform = "ont"
  elif platform == "HiFi":
    platform = "pb"


  exe="/mnt/ec/ness/yolkee/thesis/tools/correction_tools/vechat/scripts/vechat"
  opts=[
    "-t", str(threads),
    "--platform", platform,
    "-o", output,
    read
  ]
  cmd = [exe] + opts

  print(f"{TOOL} command =", ' '.join(cmd))
  proc = Popen(cmd)
  proc.wait()
  if proc.returncode != 0:
    print(f"Error: {TOOL} failed with exit code {proc.returncode}")
    exit(-1)

  # read result and remove suffix "rr" in read name 
  records = [r for r in SeqIO.parse(output, "fasta")]
  for r in records:
    r.id = r.id[:-2]
  SeqIO.write(records, output, "fasta")

if __name__ == "__main__":
  main()
