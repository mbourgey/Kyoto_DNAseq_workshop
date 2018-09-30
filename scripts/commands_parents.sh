#set up
export SOFT_DIR=/usr/local/
export WORK_DIR=$HOME/workshop/VariantCalling

export TRIMMOMATIC_JAR=$SOFT_DIR/Trimmomatic-0.38/trimmomatic-0.38.jar
export PICARD_JAR=$SOFT_DIR/picard/picard.jar
#export PICARD_JAR=/home/mathieu/workshop/VariantCalling/picard.jar
export GATK_JAR=$SOFT_DIR/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$SOFT_DIR/bvatools-1.6/bvatools-1.6-full.jar
export SNPEFF_HOME=$SOFT_DIR/snpEff/
export SNPEFF_JAR=$SNPEFF_HOME/snpEff.jar
export REF=$WORK_DIR/reference/



#########
##Father's command
#########

#zless -S raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz

zcat raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz | head -n4
zcat raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz | head -n4

mkdir -p originalQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz \
  --threads 1 --regionName ACTL8 --output originalQC/
  
#scp -r <USER>@www.genome.med.kyoto-u.ac.jp:$WORK_DIR/originalQC/ ./

cat $REF/adapters.fa
mkdir -p reads/NA12891/

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:$REF/adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12891/NA12891.trim.out

cat reads/NA12891/NA12891.trim.out
mkdir -p postTrimQC/

java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 1 --regionName ACTL8 --output postTrimQC/

mkdir -p alignment/NA12891/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12891\tSM:NA12891\tLB:NA12891\tPU:runNA12891_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/hg19.fa \
  reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx8G -jar ${PICARD_JAR} SortSam \
  INPUT=/dev/stdin \
  OUTPUT=alignment/NA12891/NA12891.sorted.bam \
  CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000

samtools view alignment/NA12891/NA12891.sorted.bam | head -n4
samtools view alignment/NA12891/NA12891.sorted.bam | grep "1313:19317:61840"

# Say you want to count the *un-aligned* reads, you can use
samtools view -c -f4 alignment/NA12891/NA12891.sorted.bam

# Or you want to count the *aligned* reads you, can use
samtools view -c -F4 alignment/NA12891/NA12891.sorted.bam

java -Xmx8G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/hg19.fa \
  -o alignment/NA12891/realign.intervals \
  -I alignment/NA12891/NA12891.sorted.bam \
  -L chr1

java -Xmx8G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/hg19.fa \
  -targetIntervals alignment/NA12891/realign.intervals \
  -o alignment/NA12891/NA12891.realigned.sorted.bam \
  -I alignment/NA12891/NA12891.sorted.bam

java -Xmx8G -jar ${PICARD_JAR} MarkDuplicates \
  REMOVE_DUPLICATES=false CREATE_INDEX=true \
  INPUT=alignment/NA12891/NA12891.realigned.sorted.bam \
  OUTPUT=alignment/NA12891/NA12891.sorted.dup.bam \
  METRICS_FILE=alignment/NA12891/NA12891.sorted.dup.metrics

#less alignment/NA12891/NA12891.sorted.dup.metrics
java -Xmx8G -jar ${GATK_JAR} \
  -T BaseRecalibrator \
  -R ${REF}/hg19.fa \
  -knownSites ${REF}/dbSNP_135_chr1.vcf.gz \
  -L chr1:17700000-18100000 \
  -o alignment/NA12891/NA12891.sorted.dup.recalibration_report.grp \
  -I alignment/NA12891/NA12891.sorted.dup.bam

java -Xmx8G -jar ${GATK_JAR} \
  -T PrintReads \
  -R ${REF}/hg19.fa \
  -BQSR alignment/NA12891/NA12891.sorted.dup.recalibration_report.grp \
  -o alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -I alignment/NA12891/NA12891.sorted.dup.bam

java  -Xmx8G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/hg19.fa \
  -o alignment/NA12891/NA12891.sorted.dup.recal.coverage \
  -I alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -L chr1:17700000-18100000

#### Look at the coverage
#less -S alignment/NA12891/NA12891.sorted.dup.recal.coverage.sample_interval_summary
java -Xmx8G -jar ${PICARD_JAR} CollectInsertSizeMetrics \
  REFERENCE_SEQUENCE=${REF}/hg19.fa \
  INPUT=alignment/NA12891/NA12891.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.tsv \
  HISTOGRAM_FILE=alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.histo.pdf \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#look at the output
#less -S alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.tsv
java -Xmx8G -jar ${PICARD_JAR} CollectAlignmentSummaryMetrics \
  REFERENCE_SEQUENCE=${REF}/hg19.fa \
  INPUT=alignment/NA12891/NA12891.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12891/NA12891.sorted.dup.recal.metric.alignment.tsv \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#### explore the results

#less -S alignment/NA12891/NA12891.sorted.dup.recal.metric.alignment.tsv

mkdir -p variants/
java -Xmx8g -jar $GATK_JAR -T HaplotypeCaller -l INFO -R ${REF}/hg19.fa \
-I alignment/NA12891/NA12891.sorted.dup.recal.bam  --variant_index_type LINEAR --variant_index_parameter 128000 -dt none \
-o variants/NA12891.hc.vcf -L chr1:17704860-18004860

java -Xmx8g -jar $GATK_JAR -T VariantFiltration \
-R ${REF}/hg19.fa --variant variants/NA12891.hc.vcf -o variants/NA12891.hc.filter.vcf --filterExpression "QD < 2.0" \
--filterExpression "FS > 200.0" \
--filterExpression "MQ < 40.0" \
--filterName QDFilter \
--filterName FSFilter \
--filterName MQFilter

java -Xmx8G -jar $SNPEFF_JAR eff -c ${REF}/snpEff.config -v -no-intergenic \
-i vcf -o vcf hg19 variants/NA12891.hc.filter.vcf >  variants/NA12891.hc.filter.snpeff.vcf
#less -S variants/NA12891.hc.filter.snpeff.vcf
#less -S variants/NA12891.hc.filter.snpeff.vcf

java -Xmx8g -jar $GATK_JAR -T VariantAnnotator -R ${REF}/hg19.fa \
--dbsnp $REF/dbSNP_135_chr1.vcf.gz --variant variants/NA12891.hc.filter.snpeff.vcf \
-o variants/NA12891.hc.filter.snpeff.dbsnp.vcf -L chr1:17704860-18004860





#########
##Mother's command
#########


#zless -S raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz

zcat raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz | head -n4
zcat raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz | head -n4

mkdir -p originalQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz \
  --threads 1 --regionName ACTL8 --output originalQC/
  
#scp -r <USER>@www.genome.med.kyoto-u.ac.jp:$WORK_DIR/originalQC/ ./

cat $REF/adapters.fa
mkdir -p reads/NA12892/

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:$REF/adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12892/NA12892.trim.out

cat reads/NA12892/NA12892.trim.out
mkdir -p postTrimQC/

java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 1 --regionName ACTL8 --output postTrimQC/

mkdir -p alignment/NA12892/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12892\tSM:NA12892\tLB:NA12892\tPU:runNA12892_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/hg19.fa \
  reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx8G -jar ${PICARD_JAR} SortSam \
  INPUT=/dev/stdin \
  OUTPUT=alignment/NA12892/NA12892.sorted.bam \
  CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000

samtools view alignment/NA12892/NA12892.sorted.bam | head -n4
samtools view alignment/NA12892/NA12892.sorted.bam | grep "1313:19317:61840"

# Say you want to count the *un-aligned* reads, you can use
samtools view -c -f4 alignment/NA12892/NA12892.sorted.bam

# Or you want to count the *aligned* reads you, can use
samtools view -c -F4 alignment/NA12892/NA12892.sorted.bam

java -Xmx8G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/hg19.fa \
  -o alignment/NA12892/realign.intervals \
  -I alignment/NA12892/NA12892.sorted.bam \
  -L chr1

java -Xmx8G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/hg19.fa \
  -targetIntervals alignment/NA12892/realign.intervals \
  -o alignment/NA12892/NA12892.realigned.sorted.bam \
  -I alignment/NA12892/NA12892.sorted.bam

java -Xmx8G -jar ${PICARD_JAR} MarkDuplicates \
  REMOVE_DUPLICATES=false CREATE_INDEX=true \
  INPUT=alignment/NA12892/NA12892.realigned.sorted.bam \
  OUTPUT=alignment/NA12892/NA12892.sorted.dup.bam \
  METRICS_FILE=alignment/NA12892/NA12892.sorted.dup.metrics

#less alignment/NA12892/NA12892.sorted.dup.metrics
java -Xmx8G -jar ${GATK_JAR} \
  -T BaseRecalibrator \
  -R ${REF}/hg19.fa \
  -knownSites ${REF}/dbSNP_135_chr1.vcf.gz \
  -L chr1:17700000-18100000 \
  -o alignment/NA12892/NA12892.sorted.dup.recalibration_report.grp \
  -I alignment/NA12892/NA12892.sorted.dup.bam

java -Xmx8G -jar ${GATK_JAR} \
  -T PrintReads \
  -R ${REF}/hg19.fa \
  -BQSR alignment/NA12892/NA12892.sorted.dup.recalibration_report.grp \
  -o alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -I alignment/NA12892/NA12892.sorted.dup.bam

java  -Xmx8G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/hg19.fa \
  -o alignment/NA12892/NA12892.sorted.dup.recal.coverage \
  -I alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -L chr1:17700000-18100000

#### Look at the coverage
#less -S alignment/NA12892/NA12892.sorted.dup.recal.coverage.sample_interval_summary
java -Xmx8G -jar ${PICARD_JAR} CollectInsertSizeMetrics \
  REFERENCE_SEQUENCE=${REF}/hg19.fa \
  INPUT=alignment/NA12892/NA12892.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.tsv \
  HISTOGRAM_FILE=alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.histo.pdf \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#look at the output
#less -S alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.tsv
java -Xmx8G -jar ${PICARD_JAR} CollectAlignmentSummaryMetrics \
  REFERENCE_SEQUENCE=${REF}/hg19.fa \
  INPUT=alignment/NA12892/NA12892.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12892/NA12892.sorted.dup.recal.metric.alignment.tsv \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#### explore the results

#less -S alignment/NA12892/NA12892.sorted.dup.recal.metric.alignment.tsv

mkdir -p variants/
java -Xmx8g -jar $GATK_JAR -T HaplotypeCaller -l INFO -R ${REF}/hg19.fa \
-I alignment/NA12892/NA12892.sorted.dup.recal.bam  --variant_index_type LINEAR --variant_index_parameter 128000 -dt none \
-o variants/NA12892.hc.vcf -L chr1:17704860-18004860

java -Xmx8g -jar $GATK_JAR -T VariantFiltration \
-R ${REF}/hg19.fa --variant variants/NA12892.hc.vcf -o variants/NA12892.hc.filter.vcf --filterExpression "QD < 2.0" \
--filterExpression "FS > 200.0" \
--filterExpression "MQ < 40.0" \
--filterName QDFilter \
--filterName FSFilter \
--filterName MQFilter

java -Xmx8G -jar $SNPEFF_JAR eff -c ${REF}/snpEff.config -v -no-intergenic \
-i vcf -o vcf hg19 variants/NA12892.hc.filter.vcf >  variants/NA12892.hc.filter.snpeff.vcf
#less -S variants/NA12892.hc.filter.snpeff.vcf
#less -S variants/NA12892.hc.filter.snpeff.vcf

java -Xmx8g -jar $GATK_JAR -T VariantAnnotator -R ${REF}/hg19.fa \
--dbsnp $REF/dbSNP_135_chr1.vcf.gz --variant variants/NA12892.hc.filter.snpeff.vcf \
-o variants/NA12892.hc.filter.snpeff.dbsnp.vcf -L chr1:17704860-18004860


