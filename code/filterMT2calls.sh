module load gatk
input=$1
output=$2

gatk FilterMutectCalls \
   -V $input \
   -O $output.filtered.vcf.gz
