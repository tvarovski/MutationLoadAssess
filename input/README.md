# Data

## Figure to Reproduce:
![Figure3 B, ref1](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/references/figureToReproduce.PNG)

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

To be updated when the project results are available

---

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
