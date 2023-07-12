#! /bin/bash

set -e 
set -x 

canu="/mnt/ec/ness/yolkee/thesis/tools/canu-2.2/bin/canu"
genomeSize="4641605"
# platform="pacbio"
platform="nanopore"
# platform="pacbio-hifi"


$canu genomeSize=$genomeSize \
      -$platform \
      -p ecoli_K12
