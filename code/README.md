# Code

---

## Figure to Reproduce:
![Figure3 B, ref1](https://raw.githubusercontent.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/main/references/figureToReproduce.PNG)

*Figure Legend:*

>The number of somatic mutations detected in each clone and the rate of accumulation of mutations per year are provided. (B) The spectra of base changes in the clones. For each base change the reverse complements are also included.

---

## Bash Submission Scripts for Argon Cluster Processing:

| File Name | Short Description |
| --- | --- |
| [Mutect2GATK.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/Mutect2GATK.sh) | Takes input bam files of a treated/tumor (fibroblast) sample and a control/normal sample (Blood) and runs Mutect2 to call variants |
| [haplotype.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/haplotype.sh) | Takes a bam file and a reference genome and runs production version of the Haplotype Caller. Returns a .vcf.gz output file |
| [haplotypeSpark.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/haplotypeSpark.sh) | Takes a bam file and a reference genome and runs beta version of the Spark Haplotype Caller. Doesn't function in my environment. Returns a .vcf.gz output file |
| [indexBamFile.sh](https://github.com/Intro-Sci-Comp-UIowa/biol-4386-course-project-tvarovski/blob/main/code/indexBamFile.sh) | Indexes a bamfile using multithreading. User can Specify number of threads to be used.


