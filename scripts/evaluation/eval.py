#! /usr/bin/python3

from subprocess import *
import os
import sys
import argparse


def main():
  parser = argparse.ArgumentParser(prog="eval")
  parser.add_argument("-d",
                      dest="dir", type=str,
                      required=True,
                      help="directory containing reads and snp genome")
  parser.add_argument("-n",
                      dest="plodiy", type=int,
                      required=True,
                      help="ploidy of the genome")
  parser.add_argument("-p", "--platform",
                      dest="platform", type=str,
                      required=True,
                      choices=["PacBio", "ONT"])
  parser.add_argument("-r", "--read",
                      dest="read", type=str,
                      required=True,
                      help="corrected read file")

                    
  args = parser.parse_args(sys.argv[1:])
  dir = args.dir
  snp_ref_dir = dir + "/snp_ref"
  reads_dir = dir + "/reads"

  if not os.path.exists(snp_ref_dir) or not os.path.isdir(snp_ref_dir):
    print(f"{snp_ref_dir} not exists or not a directory.")
    exit(-1)
  if not os.path.exists(reads_dir) or not os.path.isdir(reads_dir):
    print(f"{reads_dir} not exists or not a directory.")
    exit(-1)

  plodiy = args.plodiy
  eval_exe = "/mnt/ec/ness/yolkee/thesis/build/eval"
  eval_args = []
  for i in range(1, args.plodiy + 1):



if __name__ == "__main__":
  main()
