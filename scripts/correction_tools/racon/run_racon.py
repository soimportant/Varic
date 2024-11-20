#! /usr/bin/python3

from subprocess import *
import sys
import os
import argparse
import shutil
import tempfile
from resource import *
import time

from Bio import SeqIO

TOOL="racon"


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
  print(f"Elapsed time: {time.time() - start:.4f} s")

  return p.returncode
  

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
    "--dual=yes",
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

  exe="/mnt/ec/ness/yolkee/thesis/tools/correction_tools/racon/build/bin/racon"
  opts=[
    "-t", str(args.threads),
    "-f",
    read,
    tmp,
    read
  ]
  cmd = [exe] + opts

  print(f"{TOOL} command =", ' '.join(cmd))
  with open(output, "w") as f:
    returncode = run_cmd(cmd, f)
  
  if returncode != 0:
    print(f"Error: {TOOL} failed with exit code {returncode}")
    exit(-1)

  # racon output will be this, trim it
  # >01_2157r LN:i:13058 RC:i:86 XC:f:1.000000
  records = [r for r in SeqIO.parse(output, "fasta")]
  for r in records:
    r.id = r.id[0 : r.id.find('r')]
  SeqIO.write(records, output, "fasta")

if __name__ == "__main__":
  main()
