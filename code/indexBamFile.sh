#$ -q TELOMER,UI
#$ -pe smp 56
#$ -cwd
#$ -r y
# -j y
# -o /dev/null
#$ -l datacenter=LC

module load stack/2020.1
module load samtools/1.10_intel-19.0.5.281

samtools index -@ 55 SRR4047722.bam
