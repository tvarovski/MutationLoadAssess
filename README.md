# Mutation Load Assessment in Human Whole Genome Sequencing Data

## Introduction

*In vivo* DNA is under a constant stress of endogenous and exogenous damaging factors that can lead to changes in the genetic sequence and genome instability, that in somatic cells, are implicated in cancer and aging. However, the factors and the mechanistic role they play in causing mutations is not well understood and the ability to find and characterisze such events is a critical step for linking the potential effects of the mutagenic factors with changes in somatic cells on a genetic sequence level. The authors of the study that I want to replicate proposed a way to examine the link between UV exposure, known DNA damaging agent, and mutation within genomes of single cells (human skin fibroblasts) of healthy individuals. UV-induced DNA damage results in C→T changes and CpC→TpT dinucleotide changes, which can be used as a marker for the UV exposure. The authors of the study observed higher rates of such signatures in fibroblasts biopsied from forearm as compared to the hip of individual donors. The difference can be explained by the observation that skin around the hip is generally more protected from the UV by clothing, as compared to the skin on the forearm which is more exposed. 

Natalie Saini and her collaborators propose a strategy for finding and characterising *de novo* mutations arising within single genomes (cells), by creating fibroblast-derived clonal cell lineages and comparing their WGS mutational signature with the donor's blood sample. Since skin and blood have the same embryonic origin, comparing genomes of these two tissues and removing variants common between them allows for capturing the genetic changes that have been acquired after the developmental divergence of the two tissues, changes that have been accumulating throughout the life of the individual and have occured by the means of UV damage.

A similar method can be employed for the study of the impacts of other mutagenic treatments (e.g. hydroxyurea, mitomycin-C, ionizing radiation) to directly study their effects in various mutant backgrounds of subclonal human cell cultures by compairson of WGS data from treated and control populations of cells. Since single cell WGS is expensive, single-cell-derived clonal populations of treated and control cells can be used. I am particularly interested in studying the mutagenic effects of such treatments in knockout mutants of genes involved in the DNA repair mechanisms which would be my future direction if I can replicate the outcomes of the original studies.

---
### Project Reference Articles
This project is based on techniques and methods described in the two seminal research projects conducted under the leadership of Natalie Saini.

Similarily to the first paper that I briefly described in the introduction and by using the same methods, the second article examines more individuals from more diverse backgrounds, among 21 healthy volunteers, ranging in ages from 25 to 79 years. The authors didn't find a connection between the age and sex of the donor, however, skin cells from darker-skined individuals had a lower median mutation load by ~2.5x compared to the skin cells from lighter-skinned individuals. This difference was attributed to the difference in UV-induced mutation signatures which suggests that melanin is protective against UV DNA damage.

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

#### Reference Genome
The reference genome can be aquired by downloading it using the command below
```bash
$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
$ gunzip hg19.fa.gz
```
Next, the reference genome needs to be indexed using [samtools](http://www.htslib.org/) and dictionary needs to be created by using [PicardTools](https://broadinstitute.github.io/picard/)
```bash
$ samtools faidx hg19.fa
```
```bash
$ java -jar picard.jar CreateSequenceDictionary R=hg19.fa O=hg19.dict
```
#### Indexing Samples
Before the samples can be used by GATK they must be indexed as well. To do that you can use samtools as well.
```bash
$ samtools index SampleName.bam
```
To save time you can also use multi-threaded processing.
```bash
$ samtools index -@ <number_of_threads> SampleName.bam
```
## Variant Detection
For finding somatic variants, GATK v4.1.8.1 and Picard 2.23.0 was used.

### Mutect2
To run the Mutect2 variant detection program, first one needs to find the `<matched_blood_sample_name>` and `<fibroblast_sample_name>`. These can be extracted from the BAM read group headers by using the following command:
```bash
$ gatk GetSampleName -I <input.bam> -O <output.txt>
$ cat output.txt
```
```bash
$ gatk Mutect2 -R <reference_genome> \
    -I <fibroblast_sample> \
    -I <matched_blood_sample> \
    -normal <matched_blood_sample_name> \
    -tumor <fibroblast_sample_name> \
    -O <outputname_name.vfc.gz>
```

After the variant detection has finished (this takes a long time) it is needed to filter the calls based on probability. I used default settings. To do that the following command was used:

```bash
$ gatk FilterMutectCalls \
   -V <input.vcf.gz> \
   -O <output.filtered.vcf.gz>
```
Next, I exported the table into a TSV file. This was done by using another GATK program by taking the official GATK advice: "No, really, do not write your own parser if you can avoid it. This is not a comment on how smart or how competent we think you are -- it is a comment on how annoyingly obtuse and convoluted the VCF format is.", and "Why are we sticking with it [VCF format] anyway? Because, as Winston Churchill famously put it, VCF is the worst variant call representation, except for all the others."

Needless to say, I had my owns struggles to understand how this file format follows logic.

Options are available to select particular fields for analysis. I stuck with `CHROM`, `TYPE`, `REF`, `ALT`, `AD`, `AF`. To convert to a table:

```bash
$ gunzip <input.vcf.gz>

$ gatk VariantsToTable \
     -V <input.vcf> \
     -F CHROM -F POS -F TYPE -F REF -F ALT -GF AD -GF AF \
     -O <oOutput.tsv>
```

Then I downloaded the resulting table files and used them for further filtering and analysis with my custom python script / Jupyter Notebook.

### Haplotype Caller
First, I decided to run a Beta Spark version of Haplotype Caller for distributed computation([HaplotypeCallerSpark](https://gatk.broadinstitute.org/hc/en-us/articles/360037433931-HaplotypeCallerSpark-BETA-)) to save on processing time since the production version of the [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360036452392-HaplotypeCaller) doesn't have such functionality. The results should be nontheless comparable. To run, one can use the command below:

```bash
$ gatk HaplotypeCallerSpark  \
   -R <reference_genome> \
   -I <sample_name> \
   -O <output_name.vcf.gz>
```

Unfortunately, this method didn't work for me on the argon cluster. Instead I used the production version of HaplotypeCaller:

```bash
$ gatk --java-options "-Xmx4g" HaplotypeCaller --native-pair-hmm-threads 16 \
   -R $REFERENCE \
   -I $SAMPLE \
   -O $OUTPUT
```

HaplotypeCaller doesn't have a functionality of filtering variants from a matched normal like Mutect2 does therefore I will be writing custom python code to resolve this. Additionally, HaplotypeCaller's output (GVCF) is different from Mutect2 (VCF), so I need to find out how to use the GVCF format and how to call/filter variants based on this output file.

### Varscan2
This program is somewhat problematic. It requires a use of `samtools mpileup` to create a mpileup file. This step takes a really long time and creates enormous in size files... Next these file need to be piped into varscan's `mpileup2snp` for variant calling. I have not been able to perform this step yet.


## Removing Common Variants
Common variants between all three callers need to be removed and then filtered based on the various quality metrics. To do that I decided to write my own program as I didn't see anything that would satisfy the needs of this project. The first step is to standardize all of the outputs from the callers into a simple table.

### Dealing With [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format) (VCF) Files
The outputs of the above variatnt calling programs are in a VCF format, format that is quite difficult to work with as the columns are not always of a standard input. Therefore, I want to convert the VCF files into a CSV file that will be much easier to use during the filtering. To do that I've found a relevant python library with bioinformatics tools called `scikit-allel` that seems to be able to convert VCF files to CSV format. To aquire the library on a linux machine alongside with all of the required dependencies for full functionality type:
```bash
$ pip install scikit-allel[full]
```
Now, to convert a file from VCF to CSV a following code snippet should suffice:
```python
import allel
allel.vcf_to_csv('example.vcf', 'example.csv', fields=['CHROM', 'POS', 'DP', 'REF', 'ALT', 'QUAL'])
```
Where `fields` are names of the relevant positional and quality metrics for the variant calls as outlined by the [VCF file encoding standards](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

### Plotting The Data

The resulting datasets are then combined into one table with additional column containing the sample information. Such table can be used for making the final figure (stacked 100% percent bar chart) in MS Excel.

### Putting Everything Together.
To make make all of the filtering I am using a [PySpark](https://spark.apache.org/docs/latest/api/python/index.html) library for python. The code is available in the repository-attached [Jupyter Notebook file](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/data_parser.ipynb).

---
## Results

Preliminary Results from Mutect2 analysis follow the expectations with C->T mutations constituting the majority of the new variants. Interestingly, in corcordance with the results described by the research article #2 there seems to be no corellation between the location of the biopsied sample and observed signature which suggests that the even a short exposure intervals can have a measurable outcome. 

Preliminary summary of the Mutect Output

![Preliminary_Results](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/figures/Preliminary_Result_Mutect.png?raw=true)
>The spectra of SNP base changes in the clones as preliminairly identified by Mutect2 variant caller. For each base change the reverse complements are also included.


![Preliminary_Results_Table](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/figures/Preliminary_Result_Mutect_table.png?raw=true)
>The number of somatic mutations (SNPs) detected in each analysed clone with breakdown by base change spectra.

Already seen, some discrepancies are visible. For one sample I have ~100 less calls already than the final figure in the paper and this will likely decrease more when intersected with calls made by other programs that I'm still working on. It is possible that this is due the original summary being conducted for all mutations and not just SNPs like I have done.

---
## Discussion

To Be Continued

---
## Unexpected Contingencies

Some challenges I've encountered were mostly centered around how to acquire data in the first place. For one, it is not publicly available so I had to create an account for the database it was held in, I had to be added as a downloader for that study through my PI, I had to download proprietary software for securely downloading the data from the servers (which hasn't been without problems and errors and troubleshooting to get it to work) and later processing. Even before downloading, I wasn't prepared for acquiring 400GB of data, I had to make space for that on my Argon account where I will be continuing the project. After downloading it turned out that the data is in an unfamiliar file format (SRA) which has to be converted to BAM for it to work in bioinformatic pipelines, and it also took me a while to figure out how to do it... Also, working with latge datasets takes proportionately more amount of time for downloading and processingwhich only adds to the frustration when something doesn't work.

More challanges with using the right tool for the job, understanding the different file formats (VCF in particular), figuring out what given options for each program mean and when they are used etc. Unfortunately, none of that is described in the methods sections of my research articles so I have to go by feeling or just sticking with default settings and hoping for the best.
