import math
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

from Bio import SeqIO

def draw_kmer_profile(kmer_profile : dict, kmer_size, output_file):
  # sort the kmer profile by kmer size


  # plot the kmer profile
  kmer_count = [0] * (max(kmer_profile.values()) + 1)
  for v in kmer_profile.values():
    kmer_count[v] += 1

  # sns.set()  
  # sns.histplot(kmer_count, log_scale=(True, False))
  # plt.title("kmer profile (k = {})".format(kmer_size))
  # plt.xlabel('freq')
  # plt.ylabel('count')
  # plt.savefig(output_file)

  plt.plot(kmer_count)
  plt.title("kmer profile (k = {})".format(kmer_size))
  plt.yscale('log')
  plt.xlabel('freq')
  plt.ylabel('count')
  plt.savefig(output_file)


def main():
  # read the sequence file and kmer size for making kmer profile
  print(sys.argv)

  perfect_read_file = sys.argv[1]
  fragments_dir = sys.argv[2]
  output_dir = sys.argv[3]

  # read the raw read file fastq
  perfect_reads = list(SeqIO.parse(perfect_read_file, 'fasta'))

  for read in perfect_reads:
    read_name = read.name
    read_len = len(read.seq)
    kmer_size = math.ceil(math.log2(read_len))
    
    # read the fragments
    path = fragments_dir + '/' + read_name + '.txt'
    if os.path.exists(path) == False:
      continue
    fragments = []
    with open(fragments_dir + '/' + read_name + '.txt', 'r') as f:
      for line in f.readlines():
        s = line.split()
        s[0] = int(s[0])
        s[1] = int(s[1])
        fragments.append((s[0], s[1], s[2]))
    
    # make the kmer profile
    kmer_profile = dict()
    head_kmer_profile = dict()
    tail_kmer_profile = dict()
    boundary = 100

    for left, right, fragment in fragments:
      for i in range(len(fragment) - kmer_size + 1):
        kmer = fragment[i:i+kmer_size]
        # normal
        if kmer not in kmer_profile:
          kmer_profile[kmer] = 0
        kmer_profile[kmer] += 1

        # head
        if i + left + kmer_size < boundary:
          if kmer not in head_kmer_profile:
            head_kmer_profile[kmer] = 0
          head_kmer_profile[kmer] += 1
          
        # tail
        if i + left > len(fragment) - boundary:
          if kmer not in tail_kmer_profile:
            tail_kmer_profile[kmer] = 0
          tail_kmer_profile[kmer] += 1
    
    # check the profile of kmer in perfect read 
    perfect_kmer_profile = list()
    for i in range(len(read.seq) - kmer_size + 1):
      kmer = read.seq[i:i+kmer_size]
      if kmer not in kmer_profile:
        perfect_kmer_profile.append(0)
      else:
        perfect_kmer_profile.append(kmer_profile[kmer])
    
    # draw the kmer profile
    plt.title("Read {} (len = {})".format(read_name, read_len))
    plt.plot(perfect_kmer_profile)
    plt.savefig(output_dir + '/' + read_name + '.png')
    plt.clf()

    # draw the head profile
    head_perfect_kmer_profile = list()
    for i in range(boundary - kmer_size + 1):
      kmer = read.seq[i:i+kmer_size]
      if kmer not in head_kmer_profile:
        head_perfect_kmer_profile.append(0)
      else:
        head_perfect_kmer_profile.append(head_kmer_profile[kmer])
    plt.title("Read {} head (len = {})".format(read_name, read_len))
    plt.plot(head_perfect_kmer_profile)
    plt.savefig(output_dir + '/' + read_name + '_head.png')
    plt.clf()

    # draw the tail profile
    tail_perfect_kmer_profile = list()
    for i in range(len(read.seq) - boundary, len(read.seq) - kmer_size + 1):
      kmer = read.seq[i:i+kmer_size]
      if kmer not in tail_kmer_profile:
        tail_perfect_kmer_profile.append(0)
      else:
        tail_perfect_kmer_profile.append(tail_kmer_profile[kmer])
    plt.title("Read {} tail (len = {})".format(read_name, read_len))
    plt.plot(tail_perfect_kmer_profile)
    plt.savefig(output_dir + '/' + read_name + '_tail.png')
    plt.clf()

  exit(0)

  sequence_file = sys.argv[1]
  kmer_size = int(sys.argv[2])
  output_file = sys.argv[3]
  with open(sequence_file, 'r') as f:
    sequence = [s.split()[2] for s in f.readlines() if len(s) > 0]
  
  # make the kmer profile
  kmer_profile = dict()
  for read in sequence:
    for i in range(len(read) - kmer_size + 1):
      kmer = read[i:i+kmer_size]
      if kmer not in kmer_profile:
        kmer_profile[kmer] = 0
      kmer_profile[kmer] += 1

  # draw the kmer profile
  draw_kmer_profile(kmer_profile, kmer_size, output_file)

if __name__ == '__main__':
  main()