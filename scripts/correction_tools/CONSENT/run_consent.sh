
set -x

exe="/mnt/ec/ness/yolkee/thesis/tools/CONSENT/CONSENT-correct"
reads="~/thesis/data/Ecoli/K13/reads/ONT/merged_reads_D20.fastq"
output="~/thesis/result/CONSENT/Ecoli/K13/ONT/result.fasta"

readsTechnology="ONT"

opt=""

$exe $opt --in $reads --out $output --type $readsTechnology