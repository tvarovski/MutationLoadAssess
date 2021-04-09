#$ -q UI,TELOMER
#$ -pe smp 16
#$ -cwd
#$ -r y
# -j y
# -o /dev/null
#$ -l datacenter=LC

REFERENCE=hg19/hg19.fa

BLOOD='SRR4047715.bam'
BLOOD_name='h4.lib1'

FIBROBLAST='SRR4047707.bam'
FIBROBLAST_name='DAG_H26'

OUTPUT='file_mutect.vcf.gz'
module load gatk

gatk Mutect2 -R $REFERENCE -I $FIBROBLAST -I $BLOOD -normal $BLOOD_name -tumor $FIBROBLAST_name -O $OUTPUT
