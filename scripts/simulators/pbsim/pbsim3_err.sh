#! /bin/bash
set -e

pbsim_dir=/mnt/ec/ness/yolkee/thesis/tools/pbsim3
data_dir="/mnt/ec/ness/yolkee/thesis/data"
ref_dir=$data_dir/ref
reads_dir=$data_dir/reads
species="Ecoli"
# "Ecoli" "human" "yeast"
strains="Sakai"

ref_genome="/mnt/ec/ness/yolkee/thesis/data/ref/Ecoli/Sakai/GCF_000008865.2/ncbi_dataset/data/GCF_000008865.2/GCF_000008865.2_ASM886v2_genomic.fna"
# ref_genome="/mnt/ec/ness/yolkee/thesis/data/ref/Ecoli/K12/GCF_000005845.2/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
err_model="/mnt/ec/ness/yolkee/thesis/tools/pbsim3/data/ERRHMM-ONT.model"
depth=20

# output="$reads_dir/$species/$strains/PacBio/D$depth/PacBio_D$depth"
output="$reads_dir/$species/$strains/ONT/D$depth/ONT_D$depth"

echo $ref_genome
echo $err_model
echo $depth
echo $output

exe=$pbsim_dir/src/pbsim

$exe  --strategy wgs \
      --method errhmm \
      --errhmm $err_model \
      --depth $depth \
      --genome $ref_genome \
      --prefix $output