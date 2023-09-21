# %%

import os

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

platform="ONT"
depth=20
tool="vechat"
f1 = f"/mnt/ec/ness/yolkee/thesis/result/Ecoli/K12/{platform}/D{depth}/{tool}/h0_result.csv"
f2 = f"/mnt/ec/ness/yolkee/thesis/result/Ecoli/K12/{platform}/D{depth}/{tool}/h1_result.csv"
df1 = pd.read_csv(f1)
df2 = pd.read_csv(f2)
df = pd.concat([df1, df2], ignore_index=True)

# %%

def draw_variant_acc(mismatch, out_of_range, total):
  acc = 1 - (mismatch + out_of_range) / total
  plot = acc.plot.hist(bins=20)
  plot.set_title("variant accuracy")
  plot.set_xlabel(f"mean: {acc.mean():.3f}")
  plot.set_ylabel("count")

# print histogram of variant corrected accuracy
variation_read = df[df["covered_variants"] > 0]
draw_variant_acc(variation_read["v_mismatch"], variation_read["v_out_of_range"], variation_read["covered_variants"])

# %%

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11, 11))
def plot_end_acc(match, mismatch, insertion, deletion, front=True):
  end_len = match + mismatch + insertion + deletion
  plot = end_len.clip(0, 400).plot.hist(bins=20, ax=axes[0, 0 if front else 1])
  plot.set_title("front end length" if front else "rear end length")
  plot.set_xlabel(f"mean: {end_len.mean():.3f}")
  plot.set_ylabel("count")

  # plot deletion distribution
  deletion_portion = deletion / end_len
  plot = deletion_portion.plot.hist(bins=100, ax=axes[1, 0 if front else 1])
  plot.set_title("front deletion ratio" if front else "rear deletion ratio")
  plot.set_xlabel(f"mean length: {deletion.mean():.3f}")

plot_end_acc(df["f_match"], df["f_mismatch"], df["f_insertion"], df["f_deletion"], True)
plot_end_acc(df["b_match"], df["b_mismatch"], df["b_insertion"], df["b_deletion"], False)

# %%

# print histogram of error profile in middle part
# def plot_error_profile(readlen, mismatch, insertion, deletion):
#   mismatch_ratio = mismatch / readlen
#   insertion_ratio = insertion / readlen
#   deletion_ratio = deletion / readlen

#   # create bar plot for mismatch, insertion, deletion and stack them
#   plot = pd.DataFrame({
#     "insertion": insertion_ratio,
#     "deletion": deletion_ratio
#   }).plot.hist(bins=50)
#   plot.set_yscale("log")
#   plot.set_title("error profile in middle part")
#   plot.set_xlabel(f"mean length: {readlen.mean():.3f}")
#   plot.set_ylabel("count")

# plot_error_profile(df["read_len"], df["e_mismatch"], df["e_insertion"], df["e_deletion"])

# %%

def plot_avg_match_len(avg_len):
  plot = avg_len.plot.hist(bins=50)
  plot.set_title("average match or mismatch length when meet a deletion")
  plot.set_xlabel(f"mean length: {avg_len.mean():.3f}")
  plot.set_ylabel("count")

plot_avg_match_len(df[df["avg_del_len"] > 0]["avg_del_len"])

# %%
