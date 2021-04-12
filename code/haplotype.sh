#$ -q UI,TELOMER
#$ -pe smp 16
#$ -cwd
#$ -r y
# -j y
# -o /dev/null
#$ -l datacenter=LC


REFERENCE=hg19/hg19.fa

SAMPLE='SRR4047715.bam'
OUTPUT=SRR4047715_haplotype.g.vcf.gz

module load gatk

gatk --java-options "-Xmx4g" HaplotypeCaller --native-pair-hmm-threads 16 \
   -R $REFERENCE \
   -I $SAMPLE \
   -O $OUTPUT
