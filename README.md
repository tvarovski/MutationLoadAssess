# Mutation Load Assessment in Human Whole Genome Sequencing Data

## Introduction

Genomic mutations are believed to be one of the main contributing factors implicated in cancer and aging. These genetic changes constantly arise within the somatic cells during the lifespan and can be caused by both endogenous and exogenous mutagenic factors. The ability to find such newly arising mutations from lineage-inherited variants is a critical step for linking the potential mutagenic factors with the change in somatic mutation rates. In this project, I want to focus on replicating a way of finding these *de novo* mutations that can be captured within individual genomes by sequencing DNA from single-cell-derived, clonal populations of cells.

The research articles linked and described below examine one possible strategy for finding new mutations arising within single genomes (cells), by creating fibroblast-derived clonal cell lineages and comparing their WGS mutational signature with the donor's blood sample. Since skin and blood have the same embryonic origin, comparing genomes of these two tissues to one another allows for capturing the genetic changes that have been acquired throughout the life of the individual, that occured after the developmental divergence of the two tissues (for example UV damage).

A similar method can be employed for the study of the impact of different mutagenic treatments (e.g. hydroxyurea, mitomycin-C, ionizing radiation) to directly study their effects in various mutant backgrounds by the analysis of WGS data from treated, subclonal human cell cultures. I am particularly interested in studying the mutagenic effects of such treatments in knockout mutants of genes involved in the DNA repair mechanisms which would be my future direction if I can replicate the outcomes of the original studies.

---
### Project Reference Articles
This project is based on techniques and methods described in the two seminal research projects conducted under the leadership of Natalie Saini linked and briefly described below.

The authors of the first study proposed a way to examine the link between UV exposure and mutation within genomes of single cells (human skin fibroblasts) of healthy individuals. This analysis was enabled by the specific mutatation signuture of UV-induced DNA damage which results in C→T changes and CpC→TpT dinucleotide changes. The authors of the study observed higher rates of these signatures within individuals in fibroblasts taken from forearm as compared to the hip.

The second, more recent paper expands on the paper above by using the same methods and examines more individuals from more diverse backgrounds, among 21 healthy volunteers, ranging in ages from 25 to 79 years. The authors didn't find a connection between the age and sex of the donor, however, skin cells from darker-skined individuals had a lower median mutation load by ~2.5x compared to the skin cells from lighter-skinned individuals. This difference was attributed to the difference in UV-induced mutation signatures which suggests that melanin is protective against UV DNA damage.

[**Research Article Reference #1**](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/references/Research%20Article%201.pdf)
>Natalie Saini, Steven A. Roberts,Leszek J. Klimczak, Kin Chan, Sara A. Grimm, Shuangshuang Dai, David C. Fargo, Jayne C. Boyer, William K. Kaufmann, Jack A. Taylor, Eunjung Lee,Isidro Cortes-Ciriano, Peter J. Park, Shepherd H. Schurman, Ewa P. Malc, Piotr A. Mieczkowski, Dmitry A. Gordenin, "[The Impact of Environmental and Endogenous Damage on Somatic Mutation Load in Human Skin Fibroblasts](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006385)", *PLOS Genetics* October 27, 2016

[**Research Article Reference #2**](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/references/Research%20Article%202.pdf)
>Natalie Saini, Camille K. Giacobone, Leszek J. Klimczak, Brian N. Papas, Adam B. Burkholder, Jian-Liang Li, David C. Fargo, Re Bai, Kevin Gerrish, Cynthia L. Innes, Shepherd H. Schurman, Dmitry A. Gordenin, "[UV-exposure, endogenous DNA damage, and DNA replication errors shape the spectra of genome changes in human skin](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009302)", *PLOS Genetics* January 14, 2021

---
### Figure to Reproduce
The following *Figure 3* taken from the `Research Article #1` shows in *A* a mutation load in clones derived from skin fibroblasts biopsied from two individuals (D1 & D2), from either the left (L) or right (R) side of the body, from hip (H) or forearm (F). Each sample shows how many new mutations were aquired throughout the the individual's life in that particular clone, and estimated yearly mutation rate based on the individuals's age. *B*, breaks down the clone-specific mutations by substitution class and shows their relative proportion for each clone. Investigating the changes in the proportion of base changes spectra might give insight into the mechanism by which these mutations were aquired. In this project, I want to reproduce the outcome shown in *Figure 3-B*.

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
5. Sequencing read quality assessment and filtering.
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

The data for this project was aquired by downloading it from the [dbGaP](https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?page=list_wishlists) repository. 

Unfortunately, the dataset used is 1) not publicly available since the sequencing data are from human subjects which qualifies it is confidential medical data, and 2) it is vastly too large to be uploaded to the GitHub repository, which prevents me to be able to share this data publicly.

In my efforts, I needed to make sure there is enough storage space to hold the data (~400GB) before downloading. I am using University of Iowa's Argon cluster resources for holding and processing the data. To download the data, one must obtain a permission form the authors of the study (I got mine through my PI) and create a unique accession key (.ngc file) by following the steps described [here](https://www.ncbi.nlm.nih.gov/sra/docs/dbgap-cloud-access/).

Since the repository holds all different kinds of data, it is important to download only the data that one needs. The particular datasets (details in table below) can be selected for downloading through the dbGaP web interface which will organize your specific download request into a downloadable cart file (.krt) that will be used later for downloading all datasets at once.

Before downloading, one needs to install software for the linux environment that allows for securly downloading the datasets. After the installation, download the data by using your key and SRA or cart file for prefeching and decryption. Detailed procedure is available [here](https://www.ncbi.nlm.nih.gov/books/NBK36439/#Download.Aspera_Connect). 

Due to the space and time constraints, only the following datasets, collected from one individual, were selected for the further analysis:

| Accession # | Sample Name | Type |
| --- | --- | --- |
| SRR4047707 | D1-L-H | Skin, left hip |
| SRR4047715 | D1-blood | Blood |
| SRR4047717 | D1-L-F1 | Skin, left forearm |
| SRR4047722 | D1-R-F | Skin, right forearm |
| SRR4047723 | D1-R-H1 | Skin, right hip |

---
>Detailed procedure for the aquisition of the above datasets is shown below.

Download the latest version of the NCBI SRA Toolkit. Untar or unzip downloaded toolkit file.

```bash
$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-ubuntu64.tar.gz
```

Check the checksum. The output should be ```e21a5ba21196328e7bd1a417055f5b32```
```bash
$ md5sum -b sratoolkit.2.10.9-ubuntu64.tar.gz
```

Before running the download commands below, make sure the dbGaP repository key (.ngc) and the cart files are ready.

Download a fresh dbGaP repository key (.ngc) file and re-config the toolkit with the command below.

```bash
$ /path-to-your-sratoolkit-installation-dir/bin/vdb-config -i
```
From the sratoolkit GUI interface, import the repository key
Download dbGaP data files
Run the command below to download the files specified in the cart file.
```bash
$ /path-to-your-sratoolkit-installation-dir/bin/prefetch --ngc /path-to-ngc-file-dir/xxxxx.ngc /path-to-your-cart-file/xxxxx.krt
```
Please make sure the sratoolkit, ngc, and cart files are on the same disk drive.

The downloaded dbGaP non-SRA files need to be decrypted before use/ Run the command below to decrypt the files.
```bash
$ /path-to-your-sratoolkit-installation-dir/bin/vdb-decrypt --ngc /path-to-ngc-file-dir/xxxxx.ngc /path-to-top-level-download-dir/
```

After downloading and decrypting, the datasets are in the .sra format. To change the format to the bam file run the following command:
```bash
$ sam-dump SRRnumber | samtools view -bS - > SRRnumber.bam
```

Downloading the reference genome
```bash
$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
$ gunzip hg19.fa.gz
```
Creating an index file
```bash
$ samtools indx hg19.fa
```

Creating a dictionary file
```bash
$ java -jar picard.jar CreateSequenceDictionary \ 
      R=reference.fasta \ 
      O=reference.dict
```



---
## Results

To be updated when the project results are available

---
## Discussion

To be updated when the project results are available

---
## Unexpected Contingencies

Some challenges I've encountered were mostly centered around how to acquire data in the first place. For one, it is not publicly available so I had to create an account for the database it was held in, I had to be added as a downloader for that study through my PI, I had to download proprietary software for securely downloading the data from the servers (which hasn't been without problems and errors and troubleshooting to get it to work) and later processing. Even before downloading, I wasn't prepared for acquiring 400GB of data, I had to make space for that on my Argon account where I will be continuing the project. After downloading it turned out that the data is in an unfamiliar file format (SRA) which has to be converted to BAM for it to work in bioinformatic pipelines, and it also took me a while to figure out how to do it... Also, working with latge datasets takes proportionately more amount of time for downloading and processingwhich only adds to the frustration when something doesn't work.
