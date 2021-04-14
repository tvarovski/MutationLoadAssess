# Code

---

## Index of Bash Code Snippets, ARGON Cluster Submission Scripts, And Jupyter Notebooks:

| File Name | Short Description |
| --- | --- |
| [Mutect2GATK.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/Mutect2GATK.sh) | Takes input bam files of a treated/tumor (fibroblast) sample and a control/normal sample (Blood) and runs Mutect2 to call variants |
| [haplotype.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/haplotype.sh) | Takes a bam file and a reference genome and runs production version of the Haplotype Caller. Returns a .vcf.gz output file |
| [haplotypeSpark.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/haplotypeSpark.sh) | Takes a bam file and a reference genome and runs beta version of the Spark Haplotype Caller. Doesn't function in my environment. Returns a .vcf.gz output file |
| [indexBamFile.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/indexBamFile.sh) | Indexes a bamfile using multithreading. User can Specify number of threads to be used.
| [filterMT2calls.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/filterMT2calls.sh) | Used for filtering the output of the Mutect2 caller |
| [vcf_toTable.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/vcf_toTable.sh) | Used for converting the standard vcf output of the filterMT2calls.sh script to a table that can be parsed by my Jupyter Notebook program |
| [data_parser.ipynb](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/data_parser.ipynb) | This is a Jupyter Notebook that includes all pyton code used for analysis of data from different variant callers.|

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
