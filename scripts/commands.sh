cd /lb/project/mugqic/projects/kyotoWorkshop_2015

m mugqic/tabix/0.2.6 mugqic/igvtools/2.3.14 mugqic/snpEff/4.0 mugqic/GenomeAnalysisTK/3.2-2 mugqic/bvatools/1.4 mugqic/trimmomatic/0.32 mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.10 mugqic/picard/1.123

export REF=$(pwd)/references/

zless -S raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz

zcat raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz | head -n4
zcat raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz | head -n4

zgrep -c "^@SRR" raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz

zgrep -c "^@" raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz

mkdir originalQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz   --read2 raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz   --threads 2 --regionName SRR --output originalQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair1.fastq.gz   --read2 raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair2.fastq.gz   --threads 2 --regionName ERR --output originalQC/

cat adapters.fa

mkdir -p reads/NA12878/runSRR_1/

mkdir -p reads/NA12878/runERR_1/

java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33   raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair1.fastq.gz   raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair2.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.single1.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.single2.fastq.gz   ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32   2> reads/NA12878/runERR_1/NA12878.ERR.trim.out
java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33   raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz   raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.single1.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.single2.fastq.gz   ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32   2> reads/NA12878/runSRR_1/NA12878.SRR.trim.out

cat reads/NA12878/runERR_1/NA12878.ERR.trim.out reads/NA12878/runSRR_1/NA12878.SRR.trim.out

mkdir postTrimQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz   --read2 reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz   --threads 2 --regionName ERR --output postTrimQC/

java -Xmx1G -jar ${BVATOOLS_JAR} readsqc   --read1 reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz   --read2 reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz   --threads 2 --regionName SRR --output postTrimQC/

mkdir -p alignment/NA12878/runERR_1
mkdir -p alignment/NA12878/runSRR_1

bwa mem -M -t 2   -R '@RG\tID:ERR_ERR_1\tSM:NA12878\tLB:ERR\tPU:runERR_1\tCN:Broad Institute\tPL:ILLUMINA'   ${REF}/b37.fasta   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz   reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz   | java -Xmx2G -jar ${PICARD_HOME}/SortSam.jar   INPUT=/dev/stdin   OUTPUT=alignment/NA12878/runERR_1/NA12878.ERR.sorted.bam   CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000

bwa mem -M -t 2   -R '@RG\tID:SRR_SRR_1\tSM:NA12878\tLB:SRR\tPU:runSRR_1\tCN:Broad Institute\tPL:ILLUMINA'   ${REF}/b37.fasta   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz   reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz   | java -Xmx2G -jar ${PICARD_HOME}/SortSam.jar   INPUT=/dev/stdin   OUTPUT=alignment/NA12878/runSRR_1/NA12878.SRR.sorted.bam   CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000






