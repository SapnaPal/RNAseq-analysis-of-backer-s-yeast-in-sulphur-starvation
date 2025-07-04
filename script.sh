module load FastQC/0.12.1
fastqc -t 48 SRR19383407_1.fastq.gz
fastqc -t 48 SRR19383407_2.fastq.gz
module load fastp/0.23.4
fastp -i SRR19383407_1.fastq.gz -I SRR19383407_2.fastq.gz -o SRR19383407_1.fq.gz -O SRR19383407_2.fq.gz -w 48 -h SRR19383407.html
module load FastQC/0.12.1
fastqc -t 48 SRR19383408_1.fastq.gz
fastqc -t 48 SRR19383408_2.fastq.gz
module load fastp/0.23.4
fastp -i SRR19383408_1.fastq.gz -I SRR19383408_2.fastq.gz -o SRR19383408_1.fq.gz -O SRR19383408_2.fq.gz -w 48 -h SRR19383408.html
module load FastQC/0.12.1
fastqc -t 48 SRR19383409_1.fastq.gz
fastqc -t 48 SRR19383409_2.fastq.gz
module load fastp/0.23.4
fastp -i SRR19383409_1.fastq.gz -I SRR19383409_2.fastq.gz -o SRR19383409_1.fq.gz -O SRR19383409_2.fq.gz -w 48 -h SRR19383409.html
module load FastQC/0.12.1
fastqc -t 48 SRR19383416_1.fastq.gz
fastqc -t 48 SRR19383416_2.fastq.gz
module load fastp/0.23.4
fastp -i SRR19383416_1.fastq.gz -I SRR19383416_2.fastq.gz -o SRR19383416_1.fq.gz -O SRR19383416_2.fq.gz -w 48 -h SRR19383416.html
module load FastQC/0.12.1
fastqc -t 48 SRR19383417_1.fastq.gz
fastqc -t 48 SRR19383417_2.fastq.gz
module load fastp/0.23.4
fastp -i SRR19383417_1.fastq.gz -I SRR19383417_2.fastq.gz -o SRR19383417_1.fq.gz -O SRR19383417_2.fq.gz -w 48 -h SRR19383417.html
module load FastQC/0.12.1
fastqc -t 48 SRR19383418_1.fastq.gz
fastqc -t 48 SRR19383418_2.fastq.gz
module load fastp/0.23.4
fastp -i SRR19383418_1.fastq.gz -I SRR19383418_2.fastq.gz -o SRR19383418_1.fq.gz -O SRR19383418_2.fq.gz -w 48 -h SRR19383418.html
module load  hisat2/2.2.1
#hisat2_extract_splice_sites.py genome.gtf > genome.ss
#hisat2_extract_exons.py genome.gtf > genome.exon
#hisat2-build -p 48 --ss genome.ss --exon genome.exon genome.fa genome
hisat2 -p 48 -x ref/Saccharomyces_cerevisiae -1 SRR19383407_1.fq.gz -2 SRR19383407_2.fq.gz -S SRR19383407.sam
hisat2 -p 48 -x ref/Saccharomyces_cerevisiae -1 SRR19383408_1.fq.gz -2 SRR19383408_2.fq.gz -S SRR19383408.sam
hisat2 -p 48 -x ref/Saccharomyces_cerevisiae -1 SRR19383409_1.fq.gz -2 SRR19383409_2.fq.gz -S SRR19383409.sam
hisat2 -p 48 -x ref/Saccharomyces_cerevisiae -1 SRR19383416_1.fq.gz -2 SRR19383416_2.fq.gz -S SRR19383416.sam
hisat2 -p 48 -x ref/Saccharomyces_cerevisiae -1 SRR19383417_1.fq.gz -2 SRR19383417_2.fq.gz -S SRR19383417.sam
hisat2 -p 48 -x ref/Saccharomyces_cerevisiae -1 SRR19383418_1.fq.gz -2 SRR19383418_2.fq.gz -S SRR19383418.sam
module load samtools/1.19.2
samtools sort -@ 48 -o SRR19383407.bam SRR19383407.sam
samtools index SRR19383407.bam
samtools sort -@ 48 -o SRR19383408.bam SRR19383408.sam
samtools index SRR19383408.bam
samtools sort -@ 48 -o SRR19383409.bam SRR19383409.sam
samtools index SRR19383409.bam
samtools sort -@ 48 -o SRR19383416.bam SRR19383416.sam
samtools index SRR19383416.bam
samtools sort -@ 48 -o SRR19383417.bam SRR19383417.sam
samtools index SRR19383417.bam
samtools sort -@ 48 -o SRR19383418.bam SRR19383418.sam
samtools index SRR19383418.bam
module load qualimap/v2.3
qualimap bamqc -bam SRR19383407.bam -outfile SRR19383407_qualimap.html
qualimap bamqc -bam SRR19383408.bam -outfile SRR19383408_qualimap.html
qualimap bamqc -bam SRR19383409.bam -outfile SRR19383409_qualimap.html
qualimap bamqc -bam SRR19383416.bam -outfile SRR19383416_qualimap.html
qualimap bamqc -bam SRR19383417.bam -outfile SRR19383417_qualimap.html
qualimap bamqc -bam SRR19383418.bam -outfile SRR19383418_qualimap.html
