import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

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
  sequence_file = sys.argv[1]
  kmer_size = int(sys.argv[2])
  output_file = sys.argv[3]
  with open(sequence_file, 'r') as f:
    sequence = [s for s in f.readlines() if len(s) > 0]
  
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