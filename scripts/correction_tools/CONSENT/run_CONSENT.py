#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse

TOOL="CONSENT"

def main():
  parser = argparse.ArgumentParser(prog="run_consent")
  parser.add_argument("-r", dest="read", type=str, required=True, help="Path of raw reads")
  parser.add_argument("-p", dest="platform", type=str, required=True, choices=["PacBio", "ONT", "HiFi"], help="Platform of raw reads")
  parser.add_argument("-t", dest="threads", type=int, required=True, help="Number of threads")
  parser.add_argument("-o", dest="output", type=str, required=True, help="Path of output file")
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
    platform = "PB"
  
  exe="/mnt/ec/ness/yolkee/thesis/tools/correction_tools/CONSENT/CONSENT-correct"
  opts=[
    "--in", read,
    "--out", output,
    "--type", platform,
    "-j", str(args.threads)
  ]
  
  # cmd = ["CONSENT-correction", "--in", read, "--out", output, "--type", platform] + opts
  cmd = [exe] + opts
  
  print(f"{TOOL} cmd =", ' '.join(cmd))
  proc = Popen(cmd)
  proc.wait()
  
  # remove paf file
  
  if proc.returncode != 0:
    print(f"Error: CONSENT failed with exit code {proc.returncode}")
    exit(-1)

if __name__ == "__main__":
  main()