#  Introduction to DNA-Seq analysis 

 
by Mathieu Bourgey, _Ph.D_



## Introduction

This workshop will show you how to launch the various steps of a DNA-Seq pipeline in order to call quality variants.

We will be working on a 1000 genome sample, NA12878. You can find the whole raw data on the 1000 genome website:
<http://www.1000genomes.org/data>

NA12878 is the child of the trio while NA12891 and NA12892 are her parents.

| Mother | Father | Child |
| --- | --- | --- |
| NA12892 | NA12891 | NA12878 |


If you finish early, feel free to perform the same steps on the other two individuals: NA12891 & NA12892. 


For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.
We're going to focus on the reads extracted from a 300 kbp stretch of chromosome 1

| Chromosome | Start | End |
| --- | --- | --- |
| chr1 | 17704860 | 18004860 |


This analysis for exome sequencing or Whole genome sequencing. In both the cases the analysis is very similar and the main differences are the use of restricted area of search during the different steps. 



## Original Setup

### Acces to the server

Read these directions for information on how to log into the server.

#### To connect to the server

##### Port connection
depending on the network you are connected you may use one of this to port to connect to the workshop server:

    -- TCP connections to the port 22 from eduroam
    -- TCP connections to the port 20022 from networks other than eduroam

##### linux / Mac user
 
    * Open a terminal
    * ssh -p <PORT> <USERNAME>@gw.genome.med.kyoto-u.ac.jp

##### Windows

    * Open your terminal software (TeraTerm or PuTTy)
    * create a new connection using the following information
        .  Server : gw.genome.med.kyoto-u.ac.jp
        .  Port :   22 (from eduroam)  20022 (from the network other than eduroam)
        .  USERNAME, passwd :  provided on the first day of the course




### Software requirements

These are all already installed, but here are the original links.

  * [BVATools](http://bitbucket.org/mugqic/bvatools/downloads/)
  * [SAMTools](http://sourceforge.net/projects/samtools/)
  * [BWA](http://bio-bwa.sourceforge.net/)
  * [Genome Analysis Toolkit](http://www.broadinstitute.org/gatk/)
  * [Picard](http://broadinstitute.github.io/picard/)
  * [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  * [SnpEff](http://snpeff.sourceforge.net/)



### Environment setup

```
#set up
export SOFT_DIR=/usr/local/bin
export WORK_DIR=$HOME/workshop/VariantCalling

export TRIMMOMATIC_JAR=$SOFT_DIR/trimmomatic.jar
export PICARD_JAR=$SOFT_DIR/picard.jar
export GATK_JAR=$SOFT_DIR/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$SOFT_DIR/bvatools-full.jar
export SNPEFF_HOME=/usr/local/src/snpEff/
export SNPEFF_JAR=$SNPEFF_HOME/snpEff.jar
export REF=$WORK_DIR/reference/


rm -rf $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR
ln -s /home/mathieu/cleanCopy/* .
```

### Data files

The initial structure of your folders should look like this:   


```
ROOT  
    |-- raw_reads/               # fastqs from the center (down sampled)  
        `-- NA12878/             # Child sample directory  
        `-- NA12891/             # Father sample directory  
        `-- NA12892/             # Mother sample directory  
    |-- scripts/                 # command lines scripts         
    |-- saved_results/           # precomputed final files  
```

### Cheat sheets

* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)
* [commands file of this module](https://github.com/mbourgey/Kyoto_DNAseq_workshop/blob/master/scripts/commands.sh)


## First data glance

So you've just received an email saying that your data is ready for download from the sequencing center of your choice.

**What should you do?**  

____________________
____________________
____________________

### Fastq files

Let's first explore the fastq file.

Try these commands

```
zless -S raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz

```

These are fastq file. 

**Could you descride the fastq format?**   

____________________
____________________
____________________


```
zcat raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz | head -n4
zcat raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz | head -n4
```

**What was special about the output and why was it like that?**   

____________________
____________________
____________________


### Quality

We can't look at all the reads. Especially when working with whole genome 30x data. You could easilly have Billions of reads.

Tools like FastQC and BVATools readsqc can be used to plot many metrics from these data sets.

Let's look at the data:

```
mkdir -p originalQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  --threads 1 --regionName ACTL8 --output originalQC/
```

Copy the images from the originalQC folder to your desktop and open the images following these directions

#### To copy data from the server

##### Port connection
depending on the network you are connected you may use one of this to port to connect to the workshop server:

    -- TCP connections to the port 22 from eduroam
    -- TCP connections to the port 20022 from networks other than eduroam

##### linux / Mac user
 
    * Open a terminal
    * scp -r -P <PORT> <USERNAME>@gw.genome.med.kyoto-u.ac.jp:<DISTANT_PATH> <LOCAL_PATH>

##### Windows

    * Open your file  transfert software (WinSCP or FilleZilla)
    * create a new connection using the following information
        .  Server : gw.genome.med.kyoto-u.ac.jp
        .  Port :   22 (from eduroam)  20022 (from the network other than eduroam)
        .  USERNAME, passwd :  provided on the first day of the course
        .  Navigate to the <DISTANT_PATH>
        .  Copy the data into your <LOCAL_PATH>


**What stands out in the graphs?**  

____________________
____________________
____________________


All the generated graphics have their uses. This being said 2 of them are particularly useful to get an overal picture of how good or bad a run went. 


These are the Quality box plots 

![Quality box plots](QualityBoxPlot.png)


and the nucleotide content graphs.

![Nucleotide content](nucleotide_content.png)


The Box plot shows the quality distribution of your data. The Graph goes > 100 because both ends are appended one after the other.


The quality of a base is computated using the Phread quality score.

![Phred quality score formula](phred_formula.png)

The formula outputs an integer that is encoded using an [ASCII](http://en.wikipedia.org/wiki/ASCII) table. 

![ASCII table](ascii_table.png)


The way the lookup is done is by taking the the phred score adding 33 and using this number as a lookup in the table. The Wikipedia entry for the [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format) has a summary of the varying values.

Older illumina runs were using phred+64 instead of phred+33 to encode their fastq files.



### Trimming

After this careful analysis of the raw data we see that

- Some reads have bad 3' ends.
- no read has adapter sequences in it.

Although nowadays this doesn't happen often, it does still happen. In some cases, miRNA, it is expected to have adapters. Since they are not part of the genome of interest they should be removed if enough reads have them.


To be able to remove adapters and low qualtity beses we will use Trimmomatic. 

The adapter file is already in your reference  folder.

We can look at the adapters

```
cat $REF/adapters.fa
```

**Why are there 2 different ones ?**   

____________________
____________________
____________________


Let's try removing them and see what happens.

```
mkdir -p reads/NA12878/

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:$REF/adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12878/NA12878.trim.out

cat reads/NA12878/NA12878.trim.out
```

**What does Trimmomatic says it did ?**   

____________________
____________________
____________________

Let's look at the graphs now

```
mkdir -p postTrimQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 1 --regionName ACTL8 --output postTrimQC/
```

**How does it look now?**   

____________________
____________________
____________________



**Could we have done a better job?**   

____________________
____________________
____________________



## Alignment

The raw reads are now cleaned up of artefacts we can align the read to the reference.

In case you have multiple readsets or library you should align them separatly !

**Why should this be done separatly?**   

____________________
____________________
____________________


```
mkdir -p alignment/NA12878/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPU:runNA12878_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/hg19.fa\
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx8G -jar ${PICARD_JAR} SortSam \
  INPUT=/dev/stdin \
  OUTPUT=alignment/NA12878/NA12878.sorted.bam \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000
```
\pagebreak

**Why is it important to set Read Group information?**   

____________________
____________________
____________________



The details of the fields can be found in the SAM/BAM specifications [Here](http://samtools.sourceforge.net/SAM1.pdf)
For most cases, only the sample name, platform unit and library one are important. 


### Lane merging (optional)

In case we generate multiple lane of sequencing or mutliple library. It is not practical to keep the data splited and all the reads should be merge into one massive file. 

Since we identified the reads in the BAM with read groups, even after the merging, we can still identify the origin of each read.


### SAM/BAM

Let's spend some time to explore bam files.

try

```
samtools view alignment/NA12878/NA12878.sorted.bam | head -n4
```

Here you have examples of alignment results.
A full description of the flags can be found in the [SAM specification](http://samtools.sourceforge.net/SAM1.pdf)

Try using [picards explain flag site](https://broadinstitute.github.io/picard/explain-flags.html) to understand what is going on with your reads

The flag is the 2nd column.

**What do the flags of the first 4 reads mean?**   

____________________
____________________
____________________



You can use samtools to filter reads as well.

```
# Say you want to count the *un-aligned* reads, you can use
samtools view -c -f4 alignment/NA12878/NA12878.sorted.bam

# Or you want to count the *aligned* reads you, can use
samtools view -c -F4 alignment/NA12878/NA12878.sorted.bam

```

**How many reads mapped and unmapped were there?**   

____________________
____________________
____________________



Another useful bit of information in the SAM is the CIGAR string.
It's the 6th column in the file. This column explains how the alignment was achieved.

 * M == base aligns *but doesn't have to be a match*. A SNP will have an M even if it disagrees with the reference.
 * I == Insertion
 * D == Deletion
 * S == soft-clips. These are handy to find un removed adapters, viral insertions, etc.

An in depth explanation of the CIGAR can be found [here](http://genome.sph.umich.edu/wiki/SAM)
The exact details of the cigar string can be found in the SAM spec as well.
Another good site


## Cleaning up alignments

We started by cleaning up the raw reads. Now we need to fix and clean some alignments.

### Indel realignment (optional if you use GATK HaplotypeCaller variant caller)

The first step for this is to realign around indels and snp dense regions.
The Genome Analysis toolkit has a tool for this called IndelRealigner.

It basically runs in 2 steps

1- Find the targets
2- Realign them.

```
java -Xmx8G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/hg19.fa\
  -o alignment/NA12878/realign.intervals \
  -I alignment/NA12878/NA12878.sorted.bam \
  -L chr1

java -Xmx8G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/hg19.fa\
  -targetIntervals alignment/NA12878/realign.intervals \
  -o alignment/NA12878/NA12878.realigned.sorted.bam \
  -I alignment/NA12878/NA12878.sorted.bam

```

**How could we make this go faster?**   

____________________
____________________
____________________



**How many regions did it think needed cleaning?**   

____________________
____________________
____________________



### Mark duplicates

As the step says, this is to mark duplicate reads.

**What are duplicate reads ?**   

____________________
____________________
____________________

\pagebreak

**What are they caused by ?**   

____________________
____________________
____________________



Here we will use picards for duplicates marking:

```
java -Xmx8G -jar ${PICARD_JAR} MarkDuplicates \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/NA12878/NA12878.realigned.sorted.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.bam \
  METRICS_FILE=alignment/NA12878/NA12878.sorted.dup.metrics
```

We can look in the metrics output to see what happened.

```
less alignment/NA12878/NA12878.sorted.dup.metrics
```
**How many duplicates were there ?**   

____________________
____________________
____________________


This is very low, we expect in general <2%.


### Recalibration

This is the last BAM cleaning-up step.

The goal for this step is to try to recalibrate base quality scores. The vendors tend to inflate the values of the bases in the reads.
Also, this step tries to lower the scores of some biased motifs for some technologies.

It runs in 2 steps, 
1- Build covariates based on context and known snp sites
2- Correct the reads based on these metrics

```
java -Xmx8G -jar ${GATK_JAR} \
  -T BaseRecalibrator \
  -R ${REF}/hg19.fa\
  -knownSites ${REF}/dbSNP_135_chr1.vcf.gz \
  -L chr1:17700000-18100000 \
  -o alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -I alignment/NA12878/NA12878.sorted.dup.bam

java -Xmx8G -jar ${GATK_JAR} \
  -T PrintReads \
  -R ${REF}/hg19.fa\
  -BQSR alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -o alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -I alignment/NA12878/NA12878.sorted.dup.bam
```


## Extract Metrics (optional in the workshop)

Once your whole bam is generated, it's always a good thing to check the data again to see if everything makes sens.


### Compute coverage (optional in the workshop)

If you have data from a capture kit, you should see how well your targets worked

Both GATK and BVATools have depth of coverage tools. We wrote our own in BVAtools because
- GATK was deprecating theirs, but they changed their mind
- GATK's is very slow
- We were missing some output that we wanted from the GATK's one (GC per interval, valid pairs, etc)

Here we'll use the GATK one

```
java  -Xmx8G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/hg19.fa\
  -o alignment/NA12878/NA12878.sorted.dup.recal.coverage \
  -I alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -L chr1:17700000-18100000

#### Look at the coverage
less -S alignment/NA12878/NA12878.sorted.dup.recal.coverage.sample_interval_summary
```

Coverage is the expected ~30x. 

summaryCoverageThreshold is a usefull function to see if your coverage is uniform.
Another way is to compare the mean to the median. If both are almost equal, your coverage is pretty flat. If both are quite different, that means something is wrong in your coverage. A mix of WGS and WES would show very different mean and median values.

### Insert Size (optional in the workshop)

```
java -Xmx8G -jar ${PICARD_JAR} CollectInsertSizeMetrics \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/hg19.fa\
  INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv \
  HISTOGRAM_FILE=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.histo.pdf \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#look at the output
less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv
```

**What is the insert size and the corresponding standard deviation ?**   

____________________
____________________
____________________



**Is the insert-size important ?**   

____________________
____________________
____________________



### Alignment metrics (optional in the workshop)

For the alignment metrics, we used to use ```samtools flagstat``` but with bwa mem since some reads get broken into pieces, the numbers are a bit confusing.
You can try it if you want.

We prefer the Picard way of computing metrics

```
java -Xmx8G -jar ${PICARD_JAR} CollectAlignmentSummaryMetrics  \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/hg19.fa\
  INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#### explore the results

less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv

```

**What is the percent of aligned reads ?**   

____________________
____________________
____________________




## Calling variants with GATK
<a name="variants"></a>

We mapped the reads to hg19 and then we removed duplicate reads and realigned the reads around the indels. Let's now call SNPs in NA12878 :

```
mkdir -p variants/
java -Xmx8g -jar $GATK_JAR -T HaplotypeCaller -l INFO -R ${REF}/hg19.fa\
-I alignment/NA12878/NA12878.sorted.dup.recal.bam  --variant_index_type LINEAR --variant_index_parameter 128000 -dt none \
-o variants/NA12878.hc.vcf -L chr1:17704860-18004860

```

`-Xmx8g` instructs java to allow up 2 GB of RAM to be used for GATK.
 
`-l` INFO specifies the minimum level of logging. 
 
`-R` specifies which reference sequence to use. 
 
`-I` specifies the input BAM files. 

`--variant_index_type` LINEAR  specifies the indexing strategy to use for VCFs. LINEAR creates a Linear index with bins of equal width, specified by the Bin Width parameter. 

`--variant_index_parameter` 128000 specifies the bin width.

`-dt` NONE specifies to do not downsample the data.

`-L` indicates the reference region where SNP calling should take place  



HERE WE ARE !! variant have been called

\pagebreak

**Do you think this is enough ?**   

____________________
____________________
____________________



## Filter the variants
<a name="filter"></a>

Typically variant callers will only perform a minimal amount of filtering when presenting variant calls. 

To perform more rigorous filtering, another program must be used. In our case, we will use the *VariantFiltration* tool in GATK.

**NOTE:** The best practice when using GATK is to use the *VariantRecalibrator*. In our data set, we had too few variants to accurately use the variant recalibrator and therefore we used the *VariantFiltration* tool instead. 

```
java -Xmx8g -jar $GATK_JAR -T VariantFiltration \
-R ${REF}/hg19.fa--variant variants/NA12878.hc.vcf -o variants/NA12878.hc.filter.vcf --filterExpression "QD < 2.0" \
--filterExpression "FS > 200.0" \
--filterExpression "MQ < 40.0" \
--filterName QDFilter \
--filterName FSFilter \
--filterName MQFilter
```


`-Xmx8g` instructs java to allow up 8 GB of RAM to be used for GATK. 

`-R` specifies which reference sequence to use. 

`--variant` specifies the input vcf file. 

`-o` specifies the output vcf file. 

`--filterExpression` defines an expression using the vcf INFO and genotype variables. 

`--filterName` defines what the filter field should display if that filter is true. 

**What is QD, FS, and MQ?**    

____________________
____________________
____________________




## Adding functional consequence
<a name="function"></a>

The next step in trying to make sense of the variant calls is to assign functional consequence to each variant.

At the most basic level, this involves using gene annotations to determine if variants are sense, missense, or nonsense. 

We typically use SnpEff but many use Annovar and VEP as well.

Let's run snpEff
```
java -Xmx8G -jar $SNPEFF_JAR eff  -v -no-intergenic \
-i vcf -o vcf hg19 variants/NA12878.hc.filter.vcf >  variants/NA12878.hc.filter.snpeff.vcf
```


`-Xmx8g` instructs java to allow up 4 GB of RAM to be used for snpEff. 

`-c` specifies the path to the snpEff configuration file 

`-v` specifies verbose output. 

`-no-intergenic` specifies that we want to skip functional consequence testing in intergenic regions. 

`-i` and `-o` specify the input and output file format respectively. In this case, we specify vcf for both. 

`hg19` specifies that we want to use the hg19 annotation database. 

`variants/NA12878.hc.filter.vcf` specifies our input vcf filename 

`variants/NA12878..hc.filter.snpeff.vcf` specifies our output vcf filename 


## Investigating the functional consequence of variants
<a name="consequence"></a>

You can learn more about the meaning of snpEff annotations [here](http://snpeff.sourceforge.net/SnpEff_manual.html#output).

Use less to look at the new vcf file: 

```
less -S variants/NA12878.hc.filter.snpeff.vcf
```

We can see in the vcf that snpEff added few sections. These are hard to decipher directly from the VCF other tools or scripts, need to be used to make sens of this.


The annotation is presented in the INFO field using the new ANN format. For more information on this field see [here](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf). Typically, we have: 


`ANN=Allele|Annotation|Putative impact|Gene name|Gene ID|Feature type|Feature ID|Transcript biotype|Rank Total|HGVS.c|...`

Here's an example of a typical annotation: 

`ANN=C|intron_variant|MODIFIER|PADI6|PADI6|transcript|NM_207421.4|Coding|5/16|c.553+80T>C||||||`

**What does the example annotation actually mean?**   

____________________
____________________
____________________



additionally, you could download and look at the report generated by snpEff `snpEff_summary.html`.




**How many variants had a high impact?**   

____________________
____________________
____________________




**How many variants had a moderate impact?**   

____________________
____________________
____________________


## Adding dbSNP annotations
<a name="dbSNP"></a>

Go back to look at your last vcf file:

```
less -S variants/NA12878.hc.filter.snpeff.vcf
```

\pagebreak

**What do you see in the third column?**   

____________________
____________________
____________________



The third column in the vcf file is reserved for identifiers. Perhaps the most common identifier is the dbSNP rsID.

Use the following command to generate dbSNP rsIDs for our vcf file: 

```
java -Xmx8g -jar $GATK_JAR -T VariantAnnotator -R ${REF}/hg19.fa \
--dbsnp $REF/dbSNP_135_chr1.vcf.gz --variant variants/NA12878.hc.filter.snpeff.vcf \
-o variants/NA12878.hc.filter.snpeff.dbsnp.vcf -L chr1:17704860-18004860

```


`-Xmx2g` instructs java to allow up 2 GB of RAM to be used for GATK. 

`-R` specifies which reference sequence to use. 

`--dbsnp` specifies the input dbSNP vcf file. This is used as the source for the annotations. 

`--variant` specifies the input vcf file. 

`-o` specifies the output vcf file. 

`-L` defines which regions we should annotate. In this case, I chose the chromosomes that contain the regions we are investigating. 



**Can you find a variant that passed and wasn't in dbSNP?**   

____________________
____________________
____________________


## (Optional 1) Use IGV to investigate the SNPs

The best way to see and understand the differences between the two vcf files will be to look at them in IGV.

If you need, the IGV color codes can be found here: [IGV color code by insert size](http://software.broadinstitute.org/software/igv/interpreting_insert_size) and [IGV color code by pair orientation](http://software.broadinstitute.org/software/igv/interpreting_pair_orientations).


You can download all the NA12878.* files in the current directory to your local computer:


To do this you can use the procedure that was described previously.

After that you need to follow the steps as in Option 1 except that you need to load the files in IGV using (File->Load from File...).


Finally, go to a region on chromsome 1 with reads (chr1:17704860-18004860) and spend some time SNP gazing...

**Do the SNPs look believable?**   

____________________
____________________
____________________



**Are there any positions that you think should have been called as a SNP, but weren't?**   

____________________
____________________
____________________



## (Optional 2) Investigating the trio
<a name="trio"></a>

At this point we have aligned and called variants in one individual. However, we actually have FASTQ and BAM files for three family members!

As additional practice, perform the same steps for the other two individuals (her parents): NA12891 and NA12892. Here are some additional things that you might want to look at:

\pagebreak

 1. **If you load up all three realigned BAM files and all three final vcf files into IGV, do the variants look plausible?** Use a [Punnett square](https://en.wikipedia.org/wiki/Punnett_square) to help evaluate this. i.e. if both parents have a homozygous reference call and the child has a homozygous variant call at that locus, this might indicate a trio conflict.   

____________________
____________________
____________________



 2. **Do you find any additional high or moderate impact variants in either of the parents?**   

____________________
____________________
____________________



 3. **Do all three family members have the same genotype for Rs7538876 and Rs2254135?**   

____________________
____________________
____________________



 4. GATK produces even better variant calling results if all three BAM files are specified at the same time (i.e. specifying multiple `-I filename` options). Try this and then perform the rest of module 5 on the trio vcf file. **Does this seem to improve your variant calling results? Does it seem to reduce the trio conflict rate?**   

____________________
____________________
____________________



## Summary

In this lab, we aligned reads from the sample NA12878 to the reference genome `hg19`:

- We became familiar with FASTQ and SAM/BAM formats. 
- We checked read QC with BVAtools.
- We trimmed unreliable bases from the read ends using Trimmomatic. 
- We aligned the reads to the reference using BWA. 
- We sorted the alignments by chromosome position using PICARD. 
- We realigned short indels using GATK. 
- We recalibrate the Base Quality using GATK.
- We call and filter variants using GATK
- We annotate variant with SnpEff


## Acknowledgments


I would like to thank and acknowledge Louis Letourneau, Michael Stromberg and Guillaume Bourque for this help and for sharing his material. The format of the tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI.


-----------------------

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

-----------------------

