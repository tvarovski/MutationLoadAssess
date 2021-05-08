# Mutation Load Assessment in Human Whole Genome Sequencing Data

## Introduction

*In vivo* DNA is under the constant stress of endogenous and exogenous damaging factors that can lead to changes in the genetic sequence and genome instability. In somatic cells, nucleotide substitutions associated with DNA damage are implicated in cancer and aging. The mechanistic role the DNA damaging factors play in somatic mutagenesis is not well understood and the ability to accurately find and characterize nucleotide substitution events is a critical step for linking the potential effects of the mutagenic factors with genetic sequence changes within somatic cells. 

The authors of "The Impact of Environmental and Endogenous Damage on Somatic Mutation Load in Human Skin Fibroblasts", the study that I am replicating in this document, proposed a way to examine the link between UV exposure, a known DNA damaging agent, and mutation within genomes of single cells (human skin fibroblasts) of healthy individuals. Since UV-induced DNA damage results in C→T changes and CpC→TpT dinucleotide changes, enrichment in C→T substitutions allowed the authors of the study to establish a causal relationship between UV exposure and the skin fibroblast mutations, and distinguish those from other mutation sources. The authors of the study observed higher rates of C→T signatures in fibroblasts biopsied from the forearm as compared to the hip of individual donors. The difference could be explained by the observation that skin around the hip is generally more protected from the UV by clothing, as compared to the skin on the forearm which is more exposed. 

To assess the extent to which UV-induced DNA damage impacts somatic cells, Natalie Saini and her collaborators proposed a strategy for finding and characterizing *de novo* mutations arising within single genomes (cells), by creating fibroblast-derived clonal cell lineages and comparing their WGS mutational signature with the donor's blood sample. Since skin and blood have the same embryonic origin, comparing genomes of these two tissues and removing variants common between them allows for capturing the genetic changes that have been acquired after the developmental divergence of the two tissues, changes that have been accumulating throughout the life of the individual and have occurred by the means of UV damage.

A similar method can be employed for the study of the impacts of other mutagenic treatments (e.g. hydroxyurea, mitomycin-C, ionizing radiation) to directly study their effects in various mutant backgrounds of subclonal human cell cultures by comparison of WGS data from treated and control populations of cells. Since single-cell WGS is expensive, single-cell-derived clonal populations of treated and control cells can be used. I am particularly interested in studying the mutagenic effects of such treatments in knockout mutants of genes involved in the DNA repair mechanisms which would be my future direction if I can replicate the outcomes of the original studies.

---
### Project Reference Articles
This project is based on techniques and methods described in the two seminal research projects conducted under the leadership of Natalie Saini.

Similar to the first paper that I briefly described in the introduction and by using the same methods, the second article examines more individuals from more diverse backgrounds, among 21 healthy volunteers, ranging in ages from 25 to 79 years. The authors didn't find a connection between the age and sex of the donor, however, skin cells from darker-skinned individuals had a lower median mutation load by ~2.5x compared to the skin cells from lighter-skinned individuals. This difference was attributed to the difference in UV-induced mutation signatures which suggests that melanin is protective against UV DNA damage.

[**Research Article Reference #1**](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/references/Research%20Article%201.pdf)
>Natalie Saini, Steven A. Roberts,Leszek J. Klimczak, Kin Chan, Sara A. Grimm, Shuangshuang Dai, David C. Fargo, Jayne C. Boyer, William K. Kaufmann, Jack A. Taylor, Eunjung Lee,Isidro Cortes-Ciriano, Peter J. Park, Shepherd H. Schurman, Ewa P. Malc, Piotr A. Mieczkowski, Dmitry A. Gordenin, "[The Impact of Environmental and Endogenous Damage on Somatic Mutation Load in Human Skin Fibroblasts](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006385)", *PLOS Genetics* October 27, 2016

[**Research Article Reference #2**](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/references/Research%20Article%202.pdf)
>Natalie Saini, Camille K. Giacobone, Leszek J. Klimczak, Brian N. Papas, Adam B. Burkholder, Jian-Liang Li, David C. Fargo, Re Bai, Kevin Gerrish, Cynthia L. Innes, Shepherd H. Schurman, Dmitry A. Gordenin, "[UV-exposure, endogenous DNA damage, and DNA replication errors shape the spectra of genome changes in human skin](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009302)", *PLOS Genetics* January 14, 2021

---
### Figure to Reproduce
The following *Figure 3* taken from the `Research Article #1` shows in *A* a mutation load in clones derived from skin fibroblasts biopsied from two individuals (D1 & D2), from either the left (L) or right (R) side of the body, from the hip (H) or forearm (F). Each sample shows how many new mutations were acquired throughout the individual's life in that particular clone, and estimated yearly mutation rate based on the individuals' age. *B* breaks down the clone-specific mutations by substitution class and shows their relative proportion for each clone. Investigating the changes in the proportion of base changes spectra might give insight into the mechanism by which these mutations were acquired. In this project, I want to reproduce the outcome shown in *Figure 3-B*.

![Figure3 B, ref1](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/references/figureToReproduce.PNG)
*Figure Legend:*
>The number of somatic mutations detected in each clone and the rate of accumulation of mutations per year are provided. (B) The spectra of base changes in the clones. For each base change the reverse complements are also included.

---

## Materials And Methods For Figure Reproduction

The following is the ordered list of steps and materials required to conduct the type of analysis summarized in *Figure 3-B* from the `Research Article #1`. Steps 1-8 are presented here to show the entire process of the sample preparation and analysis. Steps 1-8 were conducted by the authors of the study and step 9 is a starting point for my efforts to reproduce the original outcomes:

1. Skin fibroblast biopsy taken from a healthy individual.
2. Isolation of a single fibroblast cell.
3. Cell culture growth of single-cell derived lineages.
4. DNA extraction from the clonal cells and blood (control), and sequencing.
5. Sequencing read the quality assessment and filtering.
6. Read alignment to the reference [human genome (GRCh37)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) using [BWA-MEM-0.7.10](https://sourceforge.net/projects/bio-bwa/files/).
7. Deduplication of reads by using [Picard Tools](https://broadinstitute.github.io/picard/) MarkDuplicates.
8. Processing BAM files according to the [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us) [best practices pipeline](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

>The preprocessed (steps 1-8) WGS data used in this study was uploaded into dbGAP under accession number ```phs001182.v1.p1``` in a BAM format.

---
9. Calling the SNVs by using three independent tools: haplotype caller from [GATK](https://gatk.broadinstitute.org/hc/en-us), [VarScan2](https://github.com/dkoboldt/varscan) and [MuTect](https://github.com/broadinstitute/mutect). Variant calling was limited to genomic regions with 10X+ coverage with 3+ reads supporting the call. Only variants common between outputs of all three tools were used in the analysis.
10. Selecting for the somatic variants that are present only in the individuals' fibroblasts, but not blood. This step removes all common mutations between the clonal fibroblast sample and blood.
11. Removal of the known Human variants - dbSNPs ([version 138](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi?view+summary=view+summary&build_id=138)), and known Human simple repeat regions (UCSC Genome Browser, [hg19 build](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)) from the variant list.
12. Removal of variant calls with allelic frequencies different than 45-55% (heterozygous) and 90%+ (homozygous).
13. Analysis and classification of the resulting variant calls according to the spectra of base changes within the clones.
14. Drawing the final analysis based on the base change spectra results in the reproduced figure.
### Data Aquisition and Pre-Processing

The data for this project was acquired by downloading it from the [dbGaP](https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?page=list_wishlists) repository. I've included a summary of this multi-step process [here](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/input/README.md#input-data). 

### Analysis Workflow

The workflow below is a high-level depiction of the steps performed to reproduce the final figure.

![Flowchart](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/figures/flowchart.png)

## Variant Detection
For finding somatic variants, GATK v4.1.8.1 and Picard 2.23.0 were used. I used two variant callers to detect mutations: Mutect2 and HaplotypeCaller. The instructions on how I used Mutect2 can be found [here](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/tree/main/code#mutect2), and instructions for HaplotypeCaller variant calling and subsequent filtering can be found [here](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/tree/main/code#haplotype-caller). After the varaints were called, all VCF outputs were exported into a TSV tables as detailed in instructions above. Due unexpected problems, my analysis does not include Varscan2, additional program that was used in the original analysis ([problem details](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/tree/main/code#varscan2)).

## Intersecting Caller Common Variants And Removing Known SNPs
Only common variants between all callers were taken into an account and filtered based on the various quality metrics. To do that I decided to write my program as I didn't find anything else that would satisfy the needs of this project. Next, I took the resulting tables and removed variants that matched positions of known SNPs in dbSNPs ([version 138](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi?view+summary=view+summary&build_id=138)). The resulting dataset was used for further analysis.

The first step is to standardize all of the outputs from the callers into a simple table. To make all of the filterings I am using a [PySpark](https://spark.apache.org/docs/latest/api/python/index.html) library for python. The code is available in the repository-attached [Jupyter Notebook file](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/data_parser.ipynb).

### Dealing With [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format) (VCF) Files
The outputs of the above variant calling programs are in a VCF format, a format that is quite difficult to work with as the columns are not always of standard input. Therefore, I converted the VCF files into a CSV file for it to be much easier to use during the filtering. To do that I've found a relevant python library with bioinformatics tools called `scikit-allel` that seems to be able to convert VCF files to CSV format. To acquire the library on a Linux machine alongside all of the required dependencies for full functionality follow my instructions [here](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/tree/main/code#dealing-with-variant-call-format-vcf-files)

## Plotting The Data

The resulting datasets were then combined into one table with one additional column containing the sample information. This table was used for making the final figure showing stacked 100% percent bar chart in MS Excel with the use of Power Pivot. The tables and figures can be viewed in the repository-attached file [here](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/figures/SpectraSummarySpreadsheet.xlsx).

---
## Results

Results from my analysis follow the mutation spectra expectations with all proportions of the events conserved and with C->T mutations constituting the majority of the detected variants. Interestingly, in concordance with the results described by research article #2 there seems to be no correlation between the location of the biopsied sample and the observed signature which suggests that even short exposure intervals can have a measurable outcome.

Mutation Spectra of My Results
![Results](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/figures/FinalSpectra.PNG)
>The spectra of SNP base changes in the analyzed clones as identified by intersecting results Mutect2 and HaplotypeCaller variant callers. Known dbSNPs ([version 138](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi?view+summary=view+summary&build_id=138)) variants were removed from the results. Only variants with allelic frequencies between 0.45-0.55 or > 0.9 were taken into account. For each base change the reverse complements are also included. Legend with spectra color codes is to the right of the figure.

![ResultsTable](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/figures/FinalSpectraTable.PNG)
>The number of somatic mutations (SNPs) detected in each analysed clone with breakdown by base change spectra. The C→T mutations constituting the majority of substitutions was highlighted in red.

Some discrepancies between the original and my reproduced figure can be observed. For one sample I have ~100 fewer calls already than the final figure in the paper and this will likely decrease more when intersected with calls made by other programs that I'm still working on. This may be due to the original summary being conducted for all mutations and not just SNPs as I have done.

## Discussion

Irrespective of minor discrepancies, my reproduced pipeline faithfully reproduces the conclusions from the original paper with trends and proportions of the mutation spectra conserved. As observed in the original study, the endogenous and exogenous factors have a comparable impact on the mutagenesis. Additionally, skin areas that are more sun-exposed acquired more C->T mutations supporting the claim that it is UV radiation that is responsible for these nucleotide changes. With the kind of analysis undertaken in this project, it is important to remember that it is likely to have a high right of false positives in the variant discovery and filtering stage due to the very stringent requirements that were established to ensure the mutation spectra were identified correctly. In the future studies it might be usefull to revisit some of those strict parameters. Finally, this pipeline can be used now to study the effects of other somatic mutagenic factors in mammalian tissue culture systems.

## Unexpected Contingencies

Some challenges I've encountered were mostly centered around how to acquire data in the first place. For one, it is not publicly available so I had to create an account for the database it was held in, I had to be added as a downloader for that study through my PI, I had to download proprietary software for securely downloading the data from the servers (which hasn't been without problems and errors and troubleshooting to get it to work) and later processing. Even before downloading, I wasn't prepared for acquiring 400GB of data, I had to make space for that on my Argon account where I will be continuing the project. After downloading it turned out that the data is in an unfamiliar file format (SRA) which has to be converted to BAM for it to work in bioinformatic pipelines, and it also took me a while to figure out how to do it... Also, working with large datasets takes a proportionately more amount of time for downloading and processing which only adds to the frustration when something doesn't work.

Trying to understand the different file formats and their structure (VCF was particularly frustrating) was time intensive and required deep engagement to be able to use these files in my pipeline. Additionally, figuring out what given options for each tool and program from GATK mean, when they are used and what for was also tedious. While GATK's tools are very well documented, information on how these tools should be used, and in which contexts is spares and confusing. Unfortunately, for most of the settings I had to go by feeling or just sticking with default settings and "hoping for the best" since none of that is described in the methods sections of my research articles.

Overall, I was very confused but in the end it worked out so I consider it time well spent.
