input=$1
output=$2

module load gatk

gatk VariantsToTable \
     -V $1 \
     -F CHROM -F POS -F TYPE -F REF -F ALT -GF AD -GF AF \
     -O $2
