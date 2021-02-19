# Introduction to Scientific Computing Project Goals

In this project I want to focus on developing a way to capture *de novo* mutations within individual genomes in populations of cells. The research articles linked and described below examine one possible strategy for finding new mutations arising within single genomes (cells) that I eventually want to use for analysis of WGS data from human cell cultures that were subjected to different mutagenic treatments (e.g. hydroxyurea, mitomycin C) to directly study their effects in various mutant backgrounds.

---

## Figure to Reproduce:
![Figure3 B, ref1](https://storage.googleapis.com/plos-corpus-prod/10.1371/journal.pgen.1006385/1/pgen.1006385.g003.PNG_L?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=wombat-sa%40plos-prod.iam.gserviceaccount.com%2F20210212%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20210212T043652Z&X-Goog-Expires=3600&X-Goog-SignedHeaders=host&X-Goog-Signature=b650ea279fdf3c0cd70b6ea27703e2816fe86a4f48c919b024b21a4755046bdfef1dd7e29e0b33df4130184c93e792c19b7da6a614c769f20e0a2b8244e9af9daffba1eecd38b2f3b79d2c48cc5a3be2b934c2f6dd96c85fac286edd0e9c8b5b4c547c4bebdf527bf2de99b8ffdae64565c2dcf5d2be10dd02455be795d7824ceb2b29e3f0464a50af5258e52aaecb2015804d89288171705141f9bcbaca8dd3954b358f5a4b76a994f85272d0c6b3db31e656dd5223e5e1dc2732ab7561c0f464d4b0cebce637a546ed6016771db49278018f462334f864e40e09cec52619e3fbb01cb337d99493c02ba34fe3974d13b38049b340583aa252a311b582283652)

*Figure Legend:*

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
13. Resulting calls were analysed and classified according to the spectra of base changes within the clones resulting in the figure 3-B below.


#### Research Article Reference #1

>Natalie Saini, Steven A. Roberts,Leszek J. Klimczak, Kin Chan, Sara A. Grimm, Shuangshuang Dai, David C. Fargo, Jayne C. Boyer, William K. Kaufmann, Jack A. Taylor, Eunjung Lee,Isidro Cortes-Ciriano, Peter J. Park, Shepherd H. Schurman, Ewa P. Malc, Piotr A. Mieczkowski, Dmitry A. Gordenin, "[The Impact of Environmental and Endogenous Damage on Somatic Mutation Load in Human Skin Fibroblasts](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006385)", *PLOS Genetics* October 27, 2016

Cancer and aging are believed to arise through mutations within the genomes. These genetic changes constantly arise within the somatic cells during the individuals lifespan and can be caused by both endogenous and exogenous factors. The authors of this study proposed a way to examine the link between UV exposure and mutation within genomes of single cells (human skin fibroblasts) of health individuals. This analysis was enabled by the specific mutatation signuture of UV-induced DNA damage which results in C→T changes and CpC→TpT dinucleotide changes. The authors of the study observed higher rates of these signatures within individuals in fibroblasts taken from forearm as compared to the hip.

#### Research Article Reference #2

>Natalie Saini, Camille K. Giacobone, Leszek J. Klimczak, Brian N. Papas, Adam B. Burkholder, Jian-Liang Li, David C. Fargo, Re Bai, Kevin Gerrish, Cynthia L. Innes, Shepherd H. Schurman, Dmitry A. Gordenin, "[UV-exposure, endogenous DNA damage, and DNA replication errors shape the spectra of genome changes in human skin](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009302)", *PLOS Genetics* January 14, 2021

The second, more recent paper exands on the paper above and examines more individuals from more diverse backgrounds among 21 healthy volunteers, ranging in ages from 25 to 79 years. The authors didn't find a connection between the age and sex of the donor, however, the skin cells from Black individuals had a lower median mutation load by ~2.5x compared to the skin cells from White individuals. This difference was attributed to the difference in UV-induced mutation signatures which suggests that melanin is protective against UV DNA damage.

---

## Results

To be updated when the project results are available

---

## Discussion

To be updated when the project results are available

---

## Conclusions

To be updated when the project results are available
