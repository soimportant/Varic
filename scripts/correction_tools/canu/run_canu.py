#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse
import shutil
import gzip
from resource import *
import time

TOOL="Canu"

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


def run_cmd(cmd, stdout=sys.stdout, stderr=sys.stderr):
  if type(cmd) == str:
    cmd = cmd.split()
    
  start = time.time()
  p = Popen(cmd)
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
  print(f"Time elapsed: {time.time() - start:.4f} s")

  return p.returncode

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
  threads = args.threads
  size = args.size

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

  exe="/mnt/ec/ness/yolkee/thesis/tools/correction_tools/canu-2.2/bin/canu"
  opts=[
    "-correct",
    f"genomeSize={size}",
    f"maxThreads={threads}",
    "-p", "tmp"
  ]
  cmd = [exe] + opts + [platform, read]

  running_log="running.log"
  with open(running_log, "w") as f:
    print(f"{TOOL} command =", ' '.join(cmd))
    returncode = run_cmd(cmd, stdout=f, stderr=f)
    if returncode != 0:
      print(f"Error: {TOOL} failed with exit code {returncode}")
      exit(-1)

  # copy result to correct directory
  canu_output = "./tmp.correctedReads.fasta.gz"
  with gzip.open(canu_output, "rb") as fin:
    with open(output, "wb") as fout:
      shutil.copyfileobj(fin, fout)

if __name__ == "__main__":
  main()