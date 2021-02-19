# Mutation Load Assessment in Human Whole Genome Sequencing Data

## Introduction

Cancer and aging are believed to arise through mutations within the genomes. These genetic changes constantly arise within the somatic cells during the lifespan and can be caused by both endogenous and exogenous mutagenic factors. The ability to distinguish such newly arising mutations from lineage-inherited variants is a critical for linking potential mutagenic factors with a change in somatic mutation rates.

In this project I want to focus on developing a way to capture *de novo* mutations within individual genomes from clonal populations of cells. The research articles linked and described below examine one possible strategy for finding new mutations arising within single genomes (cells) that I eventually want to use for analysis of WGS data from human cell cultures that were subjected to different mutagenic treatments (e.g. hydroxyurea, mitomycin-C ) to directly study their effects in various mutant backgrounds.

---
### Project Reference Articles
This project is based on techniques and methods described in the two seminal research projects conducted under the leadership of Natalie Saini linked and briefly described below.

[**Research Article Reference #1**]()
>Natalie Saini, Steven A. Roberts,Leszek J. Klimczak, Kin Chan, Sara A. Grimm, Shuangshuang Dai, David C. Fargo, Jayne C. Boyer, William K. Kaufmann, Jack A. Taylor, Eunjung Lee,Isidro Cortes-Ciriano, Peter J. Park, Shepherd H. Schurman, Ewa P. Malc, Piotr A. Mieczkowski, Dmitry A. Gordenin, "[The Impact of Environmental and Endogenous Damage on Somatic Mutation Load in Human Skin Fibroblasts](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006385)", *PLOS Genetics* October 27, 2016

The authors of this study proposed a way to examine the link between UV exposure and mutation within genomes of single cells (human skin fibroblasts) of health individuals. This analysis was enabled by the specific mutatation signuture of UV-induced DNA damage which results in C→T changes and CpC→TpT dinucleotide changes. The authors of the study observed higher rates of these signatures within individuals in fibroblasts taken from forearm as compared to the hip.

[**Research Article Reference #2**]()
>Natalie Saini, Camille K. Giacobone, Leszek J. Klimczak, Brian N. Papas, Adam B. Burkholder, Jian-Liang Li, David C. Fargo, Re Bai, Kevin Gerrish, Cynthia L. Innes, Shepherd H. Schurman, Dmitry A. Gordenin, "[UV-exposure, endogenous DNA damage, and DNA replication errors shape the spectra of genome changes in human skin](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009302)", *PLOS Genetics* January 14, 2021

The second, more recent paper expands on the paper above and examines more individuals from more diverse backgrounds, among 21 healthy volunteers, ranging in ages from 25 to 79 years. The authors didn't find a connection between the age and sex of the donor, however, skin cells from darker-skined individuals had a lower median mutation load by ~2.5x compared to the skin cells from lighter-skinned individuals. This difference was attributed to the difference in UV-induced mutation signatures which suggests that melanin is protective against UV DNA damage.

---
### Figure to Reproduce
The following *Figure 3* taken from the ```Research Article #1``` shows in *A* a mutation load in clones derived from skin fibroblasts biopsied from two individuals (D1 & D2), from either the left (L) or right (R) side of the body, from hip (H) or forearm (F). Each sample shows how many new mutations were aquired throughout the the individual's life in that particular clone, and estimated yearly mutation rated based on the individuals's age. *B*, breaks down the clone-specific mutations by substitution class and shows their relative proportion for each clone. Investigating the changes in the proportion of base changes spectra might give insight into the mechanism by which these mutations were aquired. In this project, I want to reproduce the outcome shown in *Figure 3-B*.

![Figure3 B, ref1](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/references/figureToReproduce.PNG)
*Figure Legend:*
>The number of somatic mutations detected in each clone and the rate of accumulation of mutations per year are provided. (B) The spectra of base changes in the clones. For each base change the reverse complements are also included.

---
## Materials And Methods For Figure Reproduction

The following is the ordered list of steps and materials required to conduct the type of analysis summarized in *Figure 3-B* from the ```Research Article #1```:

1. Skin fibroblast biopsy taken from a healthy individual.
2. Isolation of a single fibroblast cell.
3. Cell culture growth of single-cell derived lineages.
4. DNA extraction from the clonal cells and blood (control), and sequencing.
5. Sequencing read quality assessment and filtering.
6. Read alignment to the reference [human genome (GRCh37)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) using [BWA-MEM-0.7.10](https://sourceforge.net/projects/bio-bwa/files/).
7. Deduplication of reads by using [Picard Tools](https://broadinstitute.github.io/picard/) MarkDuplicates.
8. Processing BAM files according to the [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us) [best practices pipeline](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).
9. Calling the SNVs by using three independent tools: haplotype caller from [GATK](https://gatk.broadinstitute.org/hc/en-us), [VarScan2](https://github.com/dkoboldt/varscan) and [MuTect](https://github.com/broadinstitute/mutect). Variant calling was limited to genomic regions with 10X+ coverage with 3+ reads supporting the call. Only variants common between outputs of all three tools were used in the analysis.
10. Selecting for the somatic variants that are present only in the individuals' fibroblasts, but not blood. This step removes all common mutations between the clonal fibroblast sample and blood.
11. Removal of the known Human variants - dbSNPs ([version 138](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi?view+summary=view+summary&build_id=138)), and known Human simple repeat regions (UCSC Genome Browser, [hg19 build](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)) from the variant list.
12. Removal of variant calls with allelic frequencies different than 45-55% (heterozygous) and 90%+ (homozygous).
13. Analysis and classification of the resulting variant calls according to the spectra of base changes within the clones.
14. Drawing the final analysis based on the base change spectra results in the reproduced figure.


---
## Results

To be updated when the project results are available

---
## Discussion

To be updated when the project results are available

---
## Conclusions

To be updated when the project results are available
