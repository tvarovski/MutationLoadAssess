#$ -q UI,TELOMER
#$ -pe smp 16
#$ -cwd
#$ -r y
# -j y
# -o /dev/null
#$ -l datacenter=LC


REFERENCE=hg19/hg19.fa

SAMPLE='SRR4047707.bam'
OUTPUT=SRR4047707_haplotypeSP.vcf.gz

module load gatk


gatk HaplotypeCallerSpark  \
   -R $REFERENCE \
   -I $SAMPLE \
   -O $OUTPUT
