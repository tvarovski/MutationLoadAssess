## Figure to Reproduce:
![Figure3 B, ref1](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/references/figureToReproduce.PNG)

>The number of somatic mutations detected in each clone and the rate of accumulation of mutations per year are provided. (B) The spectra of base changes in the clones. For each base change the reverse complements are also included.

---

## Materials And Methods For Figure Reproduction:

1. Skin fibroblast biopsy taken from a healthy individual.
2. Isolation of a single fibroblast cell.
3. Cell culture growth of single-cell derived lineages.
4. DNA extraction from the clonal cells and blood (control), and sequencing.
5. Sequencing read quality assessment and filtering.
6. Read alignment to the reference [human genome (GRCh37)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) using [BWA-MEM-0.7.10](https://sourceforge.net/projects/bio-bwa/files/).
7. Deduplication of reads by using [Picard Tools](https://broadinstitute.github.io/picard/) MarkDuplicates.
8. BAM files processed according to the [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us) [best practices pipeline](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).
9. SNVs were called by using three independent tools (only variants common between all three were considered): haplotype caller from [GATK](https://gatk.broadinstitute.org/hc/en-us), [VarScan2](https://github.com/dkoboldt/varscan) and [MuTect](https://github.com/broadinstitute/mutect). Variant calling was limited to genomic regions with 10X+ coverage with 3+ reads supporting the call. 
10. Only the SNVs that were present in fibroblasts but not blood (within individual) were called as somatic.
11. Known dbSNPs ([version 138](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi?view+summary=view+summary&build_id=138)) and known simple repeat regions (UCSC Genome Browser, [hg19 build](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)) were removed from the variant list.
12. Variant calls with allelic frequencies different than 45-55% (heterozygous) and 90%+ (homozygous) were discarted.
13. Resulting calls were analysed and classified according to the spectra of base changes within the clones resulting in the figure 3-B.

