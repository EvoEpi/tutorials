# Pre-processing short reads

[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) performs a variety of useful trimming tasks for Illumina paired-end (PE) and single ended (SE) data. It is always a good idea to clean-up gDNA-seq and RNA-seq prior to mapping.

PE:

```bash
NP="" #number of processors
FASTQ="" #basename of fastq file; i.e., string before [_1|_2].fastq|[_1|_2].fq
PATH="" #path to adapter and other Illumina-specific sequences
FASTA="" #fasta file name of adapter and PCR primer sequences to trim
#other parameters are kept as default

java -jar trimmomatic-0.36.jar PE \
-phred33 \
-threads ${NP} \
-trimlog ${FASTQ}_trimmomatic.log \
${FASTQ}_1.fastq \
${FASTQ}_2.fastq \
${FASTQ}_1_P.fastq \
${FASTQ}_1_U.fastq \
${FASTQ}_2_P.fastq \
${FASTQ}_2_U.fastq \
ILLUMINACLIP:${PATH}/${FASTA}:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36
```

SE:

```bash
NP="" #number of processors
FASTQ="" #basename of fastq file; i.e., string before .fastq|.fq
PATH="" #path to adapter and other Illumina-specific sequences
FASTA="" #fasta file name of adapter and PCR primer sequences to trim
#other parameters are kept as default

java -jar trimmomatic-0.36.jar SE \
-phred33 \
-threads ${NP} \
-trimlog ${FASTQ}_trimmomatic.log \
${FASTQ}.fastq \
${FASTQ}_T.fastq \
ILLUMINACLIP:${PATH}/${FASTA}:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36
```
