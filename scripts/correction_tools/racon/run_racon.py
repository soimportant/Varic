#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse
import shutil
import tempfile

TOOL="racon"

def main():
  parser = argparse.ArgumentParser(prog="run_racon")
  parser.add_argument("-r", dest="read", type=str, required=True, help="Path of raw reads")
  parser.add_argument("-p", dest="platform", type=str, required=True, choices=["PacBio", "ONT", "HiFi"], help="Platform of raw reads")
  parser.add_argument("-o", dest="output", type=str, required=True, help="Path of output file")
  parser.add_argument("-t", dest="threads", type=int, required=False, default=2, help="Number of threads")
  args, unknown = parser.parse_known_args(sys.argv[1:])

  read = os.path.abspath(args.read)
  platform = args.platform
  output = os.path.abspath(args.output)

  # check argument
  if not os.path.exists(read) or not os.path.isfile(read):
    print(f"Raw reads {read} not exists or not a file.")
    exit(-1)
  if not os.path.exists(output):
    dir = os.path.dirname(output)
    os.makedirs(dir, exist_ok=True)

  if platform == "PacBio" or platform == "HiFi":
    platform = "ava-pb"
  elif platform == "ONT":
    platform = "ava-ont"

  # use minimap generate temporary overlap file
  # minimap2 -x ava-pb -X -t $threads $target_reads $quary_reads

  tmp = tempfile.mktemp() + ".paf"
  minimap2_cmd = [
    "minimap2",
    "-x", platform,
    "-X",
    "-t", str(args.threads),
    "-o", tmp,
    read,
    read
  ]
  minimap2_proc = Popen(minimap2_cmd)
  minimap2_proc.wait()

  if minimap2_proc.returncode != 0:
    print(f"Error: minimap2 failed with exit code {minimap2_proc.returncode}")
    exit(-1)

  exe="/mnt/ec/ness/yolkee/thesis/tools/racon/build/bin/racon"
  opts=[
    "-t", str(args.threads),
    read,
    tmp,
    read
  ]
  cmd = [exe] + opts

  print(f"{TOOL} command =", ' '.join(cmd))
  with open(output, "w") as f:
    proc = Popen(cmd, stdout=f)
    proc.wait()
  
  # clean temporary file or do something like that
  os.remove(tmp)

  if proc.returncode != 0:
    print(f"Error: {TOOL} failed with exit code {proc.returncode}")
    exit(-1)

if __name__ == "__main__":
  main()