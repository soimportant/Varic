#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse
import shutil
from resource import *
import time
from Bio import SeqIO

TOOL="vechat"


def time_wrapper(t):
  unit = {
    1000: "s",
    1: "ms"
  }
#for v, s in unit.items():
#   if t >= v:
#     return f"{t / v:.4f} {s}"
  return f"{t:.4f} s"


def memory_wrapper(m):
  unit = {
    1024*1024: "GB",
    1024: "MB",
    1: "KB",
  }
  for v, s in unit.items():
    if m >= v:
      return f"{m / v:.4f} {s}"
  return f"{m:.4f} KB"


def run_cmd(cmd, stdout = sys.stdout, stderr = sys.stderr):
  if type(cmd) == str:
    cmd = cmd.split()
  p = Popen(cmd, stdout=stdout, stderr=stderr)
  try:
    p.wait()
  except KeyboardInterrupt:
    p.kill()

  if p.returncode != 0:
    print(f"\nYour job is failed with returncode = {abs(p.returncode)}")
  else:
    print("\nYour job is successfully finished\n")

  src = getrusage(RUSAGE_CHILDREN)
  print("==========================")
  print("Resource usage of your job")
  print("==========================")

  m = [
    ("User time", src.ru_utime, time_wrapper),
    ("System time", src.ru_stime, time_wrapper),
    ("Peak Memory usage", src.ru_maxrss, memory_wrapper),
    ("Shared memory size", src.ru_ixrss, memory_wrapper),
    ("Unshared memory size", src.ru_idrss, memory_wrapper),
    ("Unshared stack size", src.ru_isrss, memory_wrapper),
    ("Page faults not requiring I/O", src.ru_minflt, None),
    ("Page faults requiring I/O", src.ru_majflt, None),
    ("Number of swap out", src.ru_nswap, None),
    ("Block input operations", src.ru_inblock, None),
    ("Block output operations", src.ru_oublock, None),
    ("Messages sent", src.ru_msgsnd, None),
    ("Messages received", src.ru_msgrcv, None),
    ("Signal received", src.ru_nsignals, None),
    ("Voluntary context switches", src.ru_nvcsw, None),
    ("Involunatary context switches", src.ru_nivcsw, None)
  ]
  for a, b, f in m:
    if f != None:
      print(f"{a:{33}}: {f(b)}")
    else:
      print(f"{a:{33}}: {b}")

  return p.returncode
  


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

  # print(f"{TOOL} command =", ' '.join(cmd))
  returncode = run_cmd(cmd)
  if returncode != 0:
    print(f"Error: {TOOL} failed with exit code {returncode}")
    exit(-1)

  # read result and remove suffix "rr" in read name 
  records = [r for r in SeqIO.parse(output, "fasta")]
  for r in records:
    r.id = r.id[:-2]
  SeqIO.write(records, output, "fasta")

if __name__ == "__main__":
  main()
