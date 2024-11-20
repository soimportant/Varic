from Bio import SeqIO
import argparse
from collections import Counter, defaultdict
import re
import readline

def read_fragment(frag_file):
  fragments = []
  with open(frag_file) as f:
    for line in f:
      fragments.append(line.strip().split())
  return fragments

def main():
  # take perfect_read as position paramter
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--perfect_read", help="perfect read file")
  args = parser.parse_args()

  # read perfect read
  reads = {r.id : r for r in SeqIO.parse(args.perfect_read, "fasta")}

  while True:
    read_id = input("Enter read id: ")
    if read_id not in reads:
      print("Read id not found")
      continue
    read = reads[read_id]
    print(f"Read len = {len(read.seq)}")
    fragments = read_fragment(f"/mnt/ec/ness/yolkee/thesis/tests/fragments/{read_id}.txt")
    len_sum = sum([len(f[2]) for f in fragments])
    print(f"Average coverage = {len_sum / len(read.seq)}")

    while True:
      kmer = input("Enter kmer: ").strip()
      if kmer == "exit" or kmer == "q":
        break
      seq = ''.join([s for s in read.seq])
      pos = [(m.start(0), m.end(0)) for m in re.finditer(kmer, seq)]
      if len(pos) == 0:
        print("Origin kmer not found")
      else:
        prev = Counter({seq[p[0] - 1] for p in pos if p[0] > 0})
        nxt = Counter({seq[p[1]] for p in pos if p[1] < len(seq)})
        if len(prev) == 1:
          prev = prev.most_common(1)[0][0]
        if len(nxt) == 1:
          nxt = nxt.most_common(1)[0][0]
        print(f"Found {kmer} at {pos}, Answer: prev, next = ({prev}, {nxt})")

      prev = Counter()
      nxt = Counter()
      max_pos = 0
      min_pos = 1e9
      for frag in fragments:
        left, right, seq = frag
        sz = len(kmer)
        for m in re.finditer(kmer, seq):
          start = m.start(0)
          end = m.end(0)
          max_pos = max(max_pos, start + int(left))
          min_pos = min(min_pos, start + int(left))
          if start != 0:
            prev.update([seq[start - 1]])
          if end < len(seq):
            nxt.update([seq[end]])
      print("({}, {}) Prev: {}, Next: {}".format(min_pos, max_pos, prev, nxt))

if __name__ == '__main__':
  main()

