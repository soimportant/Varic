#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse
import shutil

TOOL="Canu"

def main():
  parser = argparse.ArgumentParser(prog="run_canu")
  parser.add_argument("-r", dest="read", type=str, required=True, help="Path of raw reads")
  parser.add_argument("-p", dest="platform", type=str, required=True, choices=["PacBio", "ONT", "HiFi"], help="Platform of raw reads")
  parser.add_argument("-o", dest="output", type=str, required=True, help="Path of output file")
  parser.add_argument("-t", dest="threads", type=int, required=False, default=2, help="Number of threads")
  parser.add_argument("--size", dest="size", type=str, required=True, default="1g", help="Genome size")
  args, unknown = parser.parse_known_args(sys.argv[1:])

  read = os.path.abspath(args.read)
  platform = args.platform
  output = os.path.abspath(args.output)
  size = args.size
  size = 4641851

  # check argument
  if not os.path.exists(read) or not os.path.isfile(read):
    print(f"Raw reads {read} not exists or not a file.")
    exit(-1)
  if not os.path.exists(output):
    dir = os.path.dirname(output)
    os.makedirs(dir, exist_ok=True)
  if platform == "PacBio":
    platform = "-pacbio"
  elif platform == "ONT":
    platform = "-nanopore"
  elif platform == "HiFi":
    platform = "-pacbio-hifi"

  exe="/mnt/ec/ness/yolkee/thesis/tools/canu-2.2/bin/canu"
  opts=[
    "-correct",
    f"genomeSize={size}",
    "-p", "tmp"
  ]
  cmd = [exe] + opts + [platform, read]

  print(f"{TOOL} command =", ' '.join(cmd))
  proc = Popen(cmd)
  proc.wait()
  if proc.returncode != 0:
    print(f"Error: {TOOL} failed with exit code {proc.returncode}")
    exit(-1)

  # copy result to correct directory
  shutil.copy("./tmp.correctedReads.fasta.gz", output)

if __name__ == "__main__":
  main()