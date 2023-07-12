#! /bin/bash

set -x

h1=`realpath ~/thesis/data/Ecoli/K12/snp_ref/h1.simseq.genome.fa`
h2=`realpath ~/thesis/data/Ecoli/K12/snp_ref/h2.simseq.genome.fa`
chain_file="/tmp/h1.h2.chain"
reads=`realpath ~/thesis/result/vechat/Ecoli/polished/iter2/diploid_correct.fasta`

# check h1 and h2 and read exists
if [ ! -f $h1 ]; then
    echo "h1 not exists"
    exit 1
fi
if [ ! -f $h2 ]; then
    echo "h2 not exists"
    exit 1
fi
if [ ! -f $reads ]; then
    echo "reads not exists"
    exit 1
fi

# align read from h1_cread to h1
minimap2 -ax map-ont $h1 $reads -o h1.sam

# create chain file
flo $h1 $h2 $chain_file

# samtools view -bo h1.bam h1.sam

# liftover coordinate from h1 coordinate to h2
CrossMap.py bam $chain_file h1.sam h2

# # convert bam to sam
# samtools view -h -o h2.sam h2.bam
