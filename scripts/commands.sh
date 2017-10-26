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
ln -s /home/mBourgey/cleanCopy/* . ##TO UPDATE

#zless -S raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz

zcat raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz | head -n4
zcat raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz | head -n4

mkdir -p originalQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  --threads 1 --regionName ACTL8 --output originalQC/
  
#scp -r <USER>@www.genome.med.kyoto-u.ac.jp:$WORK_DIR/originalQC/ ./

cat $REF/adapters.fa
mkdir -p reads/NA12878/

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12878/NA12878.trim.out

cat reads/NA12878/NA12878.trim.out
mkdir -p postTrimQC/

java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 1 --regionName ACTL8 --output postTrimQC/

mkdir -p alignment/NA12878/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPU:runNA12878_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/genome/bwa_index/Homo_sapiens.GRCh37.fa \
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx8G -jar ${PICARD_HOME}/SortSam.jar \
  INPUT=/dev/stdin \
  OUTPUT=alignment/NA12878/NA12878.sorted.bam \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000

samtools view alignment/NA12878/NA12878.sorted.bam | head -n4
samtools view alignment/NA12878/NA12878.sorted.bam | grep "1313:19317:61840"

# Say you want to count the *un-aligned* reads, you can use
samtools view -c -f4 alignment/NA12878/NA12878.sorted.bam

# Or you want to count the *aligned* reads you, can use
samtools view -c -F4 alignment/NA12878/NA12878.sorted.bam

java -Xmx8G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12878/realign.intervals \
  -I alignment/NA12878/NA12878.sorted.bam \
  -L 1

java -Xmx8G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -targetIntervals alignment/NA12878/realign.intervals \
  -o alignment/NA12878/NA12878.realigned.sorted.bam \
  -I alignment/NA12878/NA12878.sorted.bam

java -Xmx8G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/NA12878/NA12878.realigned.sorted.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.bam \
  METRICS_FILE=alignment/NA12878/NA12878.sorted.dup.metrics

#less alignment/NA12878/NA12878.sorted.dup.metrics
java -Xmx8G -jar ${GATK_JAR} \
  -T BaseRecalibrator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -knownSites ${REF}/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
  -L 1:17700000-18100000 \
  -o alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -I alignment/NA12878/NA12878.sorted.dup.bam

java -Xmx8G -jar ${GATK_JAR} \
  -T PrintReads \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -BQSR alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -o alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -I alignment/NA12878/NA12878.sorted.dup.bam

java  -Xmx8G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12878/NA12878.sorted.dup.recal.coverage \
  -I alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -L  1:17700000-18100000

#### Look at the coverage
#less -S alignment/NA12878/NA12878.sorted.dup.recal.coverage.sample_interval_summary
java -Xmx8G -jar ${PICARD_HOME}/CollectInsertSizeMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv \
  HISTOGRAM_FILE=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.histo.pdf \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#look at the output
#less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv
java -Xmx8G -jar ${PICARD_HOME}/CollectAlignmentSummaryMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#### explore the results

#less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv

mkdir -p variants/
java -Xmx8g -jar $GATK_JAR -T HaplotypeCaller -l INFO -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
-I alignment/NA12878/NA12878.sorted.dup.recal.bam  --variant_index_type LINEAR --variant_index_parameter 128000 -dt none \
-o variants/NA12878.hc.vcf -L 1:17704860-18004860

java -Xmx8g -jar $GATK_JAR -T VariantFiltration \
-R ${REF}/genome/Homo_sapiens.GRCh37.fa --variant variants/NA12878.hc.vcf -o variants/NA12878.hc.filter.vcf --filterExpression "QD < 2.0" \
--filterExpression "FS > 200.0" \
--filterExpression "MQ < 40.0" \
--filterName QDFilter \
--filterName FSFilter \
--filterName MQFilter

java -Xmx8G -jar $SNPEFF_JAR eff  -v -no-intergenic \
-i vcf -o vcf GRCh37.75 variants/NA12878.hc.filter.vcf >  variants/NA12878.hc.filter.snpeff.vcf
#less -S variants/NA12878.rmdup.realign.hc.filter.snpeff.vcf
#less -S variants/NA12878.rmdup.realign.hc.filter.snpeff.vcf

java -Xmx8g -jar $GATK_JAR -T VariantAnnotator -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
--dbsnp $REF/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz --variant variants/NA12878.rmdup.realign.hc.filter.snpeff.vcf \
-o variants/NA12878.rmdup.realign.hc.filter.snpeff.dbsnp.vcf -L 1:17704860-18004860

