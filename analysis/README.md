# Data Analysis

## Figure to Reproduce:
![Figure3 B, ref1](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/references/figureToReproduce.PNG)

*Figure Legend:*

>The number of somatic mutations detected in each clone and the rate of accumulation of mutations per year are provided. (B) The spectra of base changes in the clones. For each base change the reverse complements are also included.

---

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
gatk --java-options "-Xmx4g" HaplotypeCaller --native-pair-hmm-threads 16 \
   -R $REFERENCE \
   -I $SAMPLE \
   -O $OUTPUT
```

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

The resulting datasets will be combined into one table with additional column containing the sample information. Such table can be used for making the final figure.

I imagine the final table to have the following structure:

| Sample | chromosome | position | reference_allele | sample_allele |
| --- | --- | --- | --- | --- |
| SRR4047707 | chr1 | 10133 | A | T |

### Putting Everything Together.
To make make all of the filtering I am using a [PySpark](https://spark.apache.org/docs/latest/api/python/index.html) library for python. The code is available in the repository attached Jupyter Notebook file.
