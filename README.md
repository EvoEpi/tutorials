# tutorials
Example workflows of various programs ranging from processing ChIP-seq, DNA-seq, RNA-seq (mRNA and sRNA), and WGBS data, and estimating gene and species phylogenies, testing for selection, and more.

## DNA-seq
Next-generation sequencing (NGS) methods provide cheap and reliable large-scale DNA sequencing. They are used extensively for _de novo_ sequencing, for disease mapping, and in population genetic studies.

__Chromosome conformation capture__. Chromosome conformation capture techniques are a set of molecular biology [methods](https://en.wikipedia.org/wiki/Chromosome_conformation_capture) used to analyze the spatial organization of chromatin in a cell. HiC, is an extension of 3C (quantifies interactions between a single pair of genomic loci) that is capable of identifying long range interactions in an unbiased, genome-wide fashion.

[juicer](https://github.com/aidenlab/juicer) is a platform for analyzing kilobase resolution HiC data. `juicer` is kind of a pain at times because it requires specific locations of the fastq and genome fasta files (i.e., in fastq/ and genome/ directories, respectively). But, hey, it works and is able to get you a chromosome-level assembly down-the-line.

Step 1. Trim reads using [HOMER](http://homer.ucsd.edu/homer/ngs/homerTools.html).

__-3__ trim adapter sequence  
__-matchStart__ don't start searching for adapter until this position, default: 0  
__-min__ remove sequence that are shorter than this after trimming  

```bash
RES="" #restriction enzyme sequence
R1="" #read 1 fastq
R2="" #read 2 fastq

homerTools trim -3 ${RES} -matchStart 20 -min 20 ${R1}.fastq > ${R1}_T.fastq
homerTools trim -3 ${RES} -matchStart 20 -min 20 ${R2}.fastq > ${R2}_T.fastq

#nuance of juicer
mkdir fastq
mv ${R1}_T.fastq ${R2}_T.fastq fastq/
```

Step 2. Prepare genome for mapping with `juicer`. This involves generating a file containing a huge list of RE cut site positions, a chromosome sizes file, and indexing the genome using `bwa`.

```bash
RE="" #restriction enzyme name; e.g., MboI
GENOME="" #genome identifier; e.g., myGenome_v1.0
REF="" #reference fasta file; e.g., myGenome_v1.0.fasta
#pay attention to the string for GENOME and the REF fasta file name

#nuance of juicer
mkdir genome
cp ${REF} genome/
cd genome/

python2 generate_site_positions.py \
${RE} \
${GENOME} \
${REF}

awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${GENOME}_${RE}.txt > ${GENOME}.chrom.sizes

bwa index ${REF}
```

Step 4. Map HiC reads using `juicer`.

__-d__ path to the top level directory  
__-s__ restriction enzyme  
__-z__ path for reference sequence file, bwa index files must be in same directory  
__-y__ path for restriction site file  
__-t__ threads  

```bash
GENOME="" #genome identifier; e.g., myGenome_v1.0
TOP="" #path of top directory
RE="" #restriction enzyme name; e.g., MboI
DIR="" #path to directory one above genome/
NP="" #number of processors

juicer.sh \
-g ${GENOME} \
-d ${TOP} \
-s ${RE} \
-z ${DIR}/${REF} \
-p ${DIR}/${GENOME}.chrom.sizes \
-y ${DIR}/${GENOME}_${RE}.txt \
-t ${NP}
```

The `.hic` file in directory `aligned/` produced by `juicer` can be loaded into [Juicebox](https://github.com/aidenlab/Juicebox/wiki) for visualization.

The `merged_nodups.txt` in directory `aligned/` along with the reference fasta file are used by [3D-DNA](https://github.com/theaidenlab/3d-dna) for genome assembling. `3D-DNA` pipeline consists of one bash wrapper script `run-asm-pipeline.sh` that calls individual modules to assemble a genome. Running with default parameters, and HiC sequencing depth of ~20X, I went from a genome assembly of 5745 scaffolds with N50=4.9Mb to 5505 scaffolds with N50=707Mb.

```bash
REF="" #reference fasta file HiC reads are mapped to

run-asm-pipeline.sh genome/${REF} aligned/merged_nodups.txt
```

__Single Nucleotide Polymorphisms (SNP) calling__. Disease mapping and population genetic studies typically involve aligning short sequencing reads from one or more individuals to a reference genome to identify variable sites or Single Nucleotide Polymorphisms (SNPs). The frequency of variants within a population can provide information regarding the genetic basis of diseases, loci under natural selection, and loci segregating in a population. In this example I go over mapping short-read re-sequencing data to a reference genome assembly followed by SNP calling.

Step 1. Index genome using [bwa](http://bio-bwa.sourceforge.net/).

```bash
REF="" #reference fasta file

bwa index ${REF}
```

Step 2. Clean-up reads using `Trimmomatic`. See `Trimmomatic` example below.

Step 3. Map reads using `bwa`. In this example we will be merging multiple bam files downstream, so we beed to add a read group header to distinguish alignments in the merged bam. We will pipe the output [SAMtools](http://www.htslib.org/) to create a sorted bam file. In this example we are mapping paired-end data.

```bash
NP="" #number of processors/threads
SAMPLE="" #unique identifier for sample
REF="" #reference fasta file
R1="" #read 1 fastq file
R2="" #read 2 fastq file

bwa mem \
-t ${NP} \
-R \"@RG\tID:${SAMPLE}\tSM:${SAMPLE}\" \
${REF} \
${R1} \
${R2} | \
samtools sort \
-o ${SAMPLE}_aln-pe_sorted.bam \
-@ ${NP} \
-
```

Step 4. Merge bam files using [bamtools](https://github.com/pezmaster31/bamtools).

```bash
BAM="" #space-separated list of bam files
OUT="" #name of merged bam file ending in .bam

bamtools merge \
-in ${BAM} \
-out ${OUT}
```

Step 5. We will be using [freebayes](https://github.com/ekg/freebayes/blob/master/README.md) to call Single Nucleotide Polymorphisms (SNPs). Specifically, `freebayes-parallel` to speed up SNP calling by breaking our reference genome assembly into regions. The default parameters of `freebayes` are very liberal, so I recommend increasing their stringency. Below are parameters I typically change (->) from their default value, but you should choose values that best suits your data:

__-C__ Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 2->5  
__-3__ Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 0->300  
__-p__ Sets the default ploidy for the analysis to N. default: 2->2  
__-m__ Exclude alignments from analysis if they have a mapping quality less than Q. default: 1->30  
__-q__ Exclude alleles from analysis if their supporting base quality is less than Q. default: 0->30  
__-n__ Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. (Set to 0 to use all; default: all)->4   

```bash
REF="" #reference fasta file
REGIONS="" #name of regions outfile
NP="" #number of processors
BAM="" #merged bam file
OUT="" #name of vcf outfile

samtools faidx \
${REF}

fasta_generate_regions.py \
${REF}.fai \
100000 > \
${REGIONS}

freebayes-parallel \
${REGIONS} \
${NP} \
-f ${REF} \
-C 5 \
-3 300 \
-p 2 \
-m 30 \
-q 30 \
-n 4 \
${BAM} \
> ${OUT}
```

## Pre-processing short reads
[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) performs a variety of useful trimming tasks for Illumina paired-end (PE) and single ended (SE) data. It is always a good idea to clean-up gDNA-seq and RNA-seq prior to mapping.

PE:

```bash
NP="" #number of processors
FASTQ="" #basename of fastq file; i.e., string before [_1|_2].fastq|[_1|_2].fq
PATH="" #path to adapter and other Illumina-specific sequences
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
ILLUMINACLIP:${PATH}/TruSeq2-PE.fa:2:30:10 \
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
#other parameters are kept as default

java -jar trimmomatic-0.36.jar SE \
-phred33 \
-threads ${NP} \
-trimlog ${FASTQ}_trimmomatic.log \
${FASTQ}.fastq \
${FASTQ}_T.fastq \
ILLUMINACLIP:${PATH}/TruSeq2-SE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36
```

## RNA-seq (mRNA)
A popular toolset used for analyzing RNA-seq data is the tuxedo suite, which consists of [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) and [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/). The suite provided a start to finish pipeline that allowed users to map reads, assemble transcripts, and perform differential expression analyses. A newer "tuxedo suite" has been developed and is made up of three tools: `HISAT`, `StringTie`, and `Ballgown`. [Pertea et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27560171) provides a summary of the new suite as well as a tutorial.

Step 1. Mapping is performed using `HISAT2` and usually the first step, prior to mapping, is to create an index of the reference genome. Caution: This step takes a lot of memory. _If you use --snp, --ss, and/or --exon, [hisat2-build](https://ccb.jhu.edu/software/hisat2/manual.shtml) will need about **200GB RAM** for the **human genome** size as index building involves a graph construction. Otherwise, you will be able to build an index on your desktop with 8GB RAM_.

```bash
SAMPLE="" #identifier associated with fastq
SPECIES="" #species
GFF="" #gff3 name
REF="" #genome reference fasta file
INDEX="" #basename of the index files to write
NP="" #number of processors
FASTQ="" #name of fastq file

gffread ${GFF} -T -o ${SPECIES}.gtf
python extract_splice_sites.py ${SPECIES}.gtf > ${SPECIES}.ss
python extract_exons.py ${SPECIES}.gtf > ${SPECIES}.exon
hisat2-build \
-p ${NP} \
--ss ${SPECIES}.ss \
--exon ${SPECIES}.exon \
${REF} \
${INDEX}

mkdir ${SAMPLE}_map

hisat2 -p ${NP} \
--dta \
-x ${INDEX} \
-U ${FASTQ} \
-S ${SAMPLE}_map/${SAMPLE}.sam

samtools sort \
-@ ${NP} \
-o ${SAMPLE}_map/${SAMPLE}.bam \
${SAMPLE}_map/${SAMPLE}.sam

rm ${SAMPLE}_map/*.sam
```

Step 2. Now we need to assemble the mapped reads into transcripts. `StringTie` can assemble transcripts with or without annotation. With annotation:

```bash
SAMPLE="" #identifier associated with fastq
SPECIES="" #species
NP="" #number of processors

mkdir ${SAMPLE}_assembly

stringtie ${SAMPLE}_map/${SAMPLE}.bam \
-l ${SAMPLE} \
-p ${NP} \
-G ${SAMPLE}.gtf \
-o ${SAMPLE}_assembly/${SPECIES}.gtf \
-e \
-A ${SAMPLE}_abundance.out
```

## WGBS
WGBS allows the interrogation of the methylation status at a single cytosine ([Cokus et al. 2008](https://www.ncbi.nlm.nih.gov/pubmed/18278030); [Lister et al. 2008](https://www.ncbi.nlm.nih.gov/pubmed/18423832)). This process uses sodium bisulfite to convert unmethylated cytosine to uracil and ultimately thymine via PCR ([Clark et al. 1994](https://www.ncbi.nlm.nih.gov/pubmed/8065911)). These can then be detected by sequencing the converted product and mapping the data to a reference genome. Reads that contain a thymine where the reference genome contains a cytosine indicate that the reference cytosine is unmethylated, whereas reads that still retain a cytosine indicate that the reference cytosine is methylated.

While there are many programs available to map WGBS reads ([Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/), [BSseeker2](https://github.com/BSSeeker/BSseeker2), [Methylkit](https://bioconductor.org/packages/release/bioc/html/methylKit.html), etc.), my preferred program is [methylpy](https://github.com/yupenghe/methylpy). `methylpy` is versatile and offers a number of statistical and analytical advantages over other programs. Specifically, discretely binning sites as methylated or unmethylated using binomial statistics and a Differentially Methylated Region (DMR) finder.

Step 1. Index your reference genome assembly.

```bash
REF="" #genome reference fasta file
INDEX="" #basename of the index files to write

methylpy build-reference \
--input-files ${REF} \
--output-prefix ${INDEX} \
--bowtie2 True
```

Step 2. Map reads and generate an allc file with output from a binomial test.

```bash
FASTQ="" #WGBS/methylC-seq fastq file
SAMPLE="" #identifier associated with fastq
FOR="" #basename of the index reference followed by _f
REV="" #basename of the index reference followed by _r
REF="" #genome reference fasta file
OUT="" #path to output directory
NP="" #number of processors
MEM="" #gb for sorting; e.g., '24G'
PICARD="" #path to picard
CONTROL="" #seqeunce in fasta file for sodium bisulfite non-conversion rate estimation

methylpy single-end-pipeline \
--read-files ${FASTQ} \
--sample ${SAMPLE} \
--forward-ref ${FOR} \
--reverse-ref ${REV} \
--ref-fasta ${REF} \
--path-to-output ${OUT} \
--num-procs ${NP} \
--sort-mem "${MEM}" \
--generate-allc-file False \
--trim-reads True \
--remove-clonal True \
--path-to-picard ${PICARD} \
--min-qual-score 10 \
--min-read-len 30

methylpy call-methylation-state \
--input-file ${SAMPLE}_processed_reads_no_clonal.bam \
--sample ${SAMPLE} \
--ref-fasta ${REF} \
--paired-end False \
--unmethylated-control ${CONTROL} \
--binom-test True \
--sig-cutoff 0.01 \
--min-cov 3
```

Optional. Run DMRfind. DMRfind is ideal for finding differences between treatments and control, tissues, time points, etc., and can handle ≥2 allc files. Caution: ensure allc file names match sample identifiers and chromosomes in allc file are strictly numerical starting at 1.

```bash
ALLC_CON="" #experimental control allc file; e.g., allc_control.tsv
ALLC_TRE="" #experimental treatment allc file; e.g., allc_treatment.tsv
CON="" #identifier for control that matches allc file name; e.g., control
TRE="" #identifier for treatment that matches allc file name; e.g., treatment
MC="" #sequence context; e.g., CGN, CHG, or CHH
CHR="" #space separated list of chromosomes; names must match allc file; e.g., 1 2 3 4 5 
NP="" #number of processors
OUT="" #output prefix

methylpy DMRfind \
--allc-files ${ALLC_CON} ${ALLC_TRE} \
--samples ${CON} ${TRE} \
--mc-type "${MC}" \
--chroms ${CHR} \
--num-procs ${NP} \
--output-prefix ${OUT}
```

## Phylogenetics

[RAxML (Randomized Axelerated Maximum Likelihood)](https://cme.h-its.org/exelixis/web/software/raxml/index.html) is a program for sequential and parallel Maximum Likelihood (ML) based inference of large phylogenetic trees. It can also be used for postanalyses of sets of phylogenetic trees, analyses of alignments and, evolutionary placement of short reads.

Step 1. Generate ${NR} ML trees on distinct starting trees and also print the tree with the best likelihood to a file called RAxML_bestTree.${ID}.

```bash
NP="" #number of processors
ALIGN="" #fasta or phylip aligned sequences
ID="" #identifier
MODEL="" #model of substitution
NR="" #number of replicates

raxmlHPC-PTHREADS-AVX \
-T ${NP} \
-s ${ALIGN} \
-n ${ID}-ML \
-m ${MODEL} \
-p 1234 \
-#${NR} \
> ${ID}-ML.log
```

Step 2. Get bootstrap support values for the best ML tree.

```bash
NP="" #number of processors
ALIGN="" #fasta or phylip aligned sequences
ID="" #identifier
MODEL="" #model of substitution
NR="" #number of replicates

raxmlHPC-PTHREADS-AVX \
-T ${NP} \
-s ${ALIGN} \
-n ${ID}-BS \
-m ${MODEL} \
-p 1234 \
-b 1234 \
-#${NR} \
> ${ID}-BS.log
```

Step 3. Having computed the bootstrap replicate trees that will be printed to a file called RAxML_bootstrap.${ID} we can now use them to draw bipartitions on the best ML tree as follows.

```bash
NP="" #number of processors
ALIGN="" #fasta or phylip aligned sequences
ID="" #identifier
MODEL="" #model of substitution
NR="" #number of replicates

raxmlHPC-PTHREADS-AVX \
-T ${NP} \
-m ${MODEL} \
-p 1234 \
-f b \
-t RAxML_bestTree.${ID}-ML \
-z RAxML_bootstrap.${ID}-BS \
-n ${ID}-TR
```
