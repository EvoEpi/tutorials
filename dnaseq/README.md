# DNA-seq

Next-generation sequencing (NGS) methods provide cheap and reliable large-scale DNA sequencing. They are used extensively for _de novo_ sequencing, for disease mapping, and in population genetic studies.

__Chromosome conformation capture__. Chromosome conformation capture techniques are a set of molecular biology [methods](https://en.wikipedia.org/wiki/Chromosome_conformation_capture) used to analyze the spatial organization of chromatin in a cell. HiC, is an extension of 3C (quantifies interactions between a single pair of genomic loci) that is capable of identifying long range interactions in an unbiased, genome-wide fashion.

[juicer](https://github.com/aidenlab/juicer) is a platform for analyzing kilobase resolution HiC data. `juicer` is kind of a pain at times because it requires specific locations of the fastq and genome fasta files (i.e., in fastq/ and genome/ directories, respectively). But, hey, it works and is able to get you a chromosome-level assembly down-the-line.

__Step 1__. Trim reads using [HOMER](http://homer.ucsd.edu/homer/ngs/homerTools.html).

__-3__ Trim adapter sequence  
__-matchStart__ Don't start searching for adapter until this position, default: 0  
__-min__ Remove sequence that are shorter than this after trimming  

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

__Step 2__. Prepare genome for mapping with `juicer`. This involves generating a file containing a huge list of RE cut site positions, a chromosome sizes file, and indexing the genome using `bwa`.

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

__Step 4__. Map HiC reads using `juicer`.

__-d__ Path to the top level directory  
__-s__ Restriction enzyme  
__-z__ Path for reference sequence file, bwa index files must be in same directory  
__-y__ Path for restriction site file  
__-t__ Threads  

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

__Step 1__. Index genome using [bwa](http://bio-bwa.sourceforge.net/).

```bash
REF="" #reference fasta file

bwa index ${REF}
```

__Step 2__. Clean-up reads using `Trimmomatic`. See `Trimmomatic` example in "Pre-processing reads" directory.

Step 3. Map reads using `bwa`. In this example we will be merging multiple bam files downstream, so we beed to add a read group header to distinguish alignments in the merged bam. We will pipe the output [SAMtools](http://www.htslib.org/) to create a sorted bam file. In this example we are mapping paired-end data.

__-t__ Number of threads  
__-R__ Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’  

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

__Step 4__. Merge bam files using [bamtools](https://github.com/pezmaster31/bamtools).

```bash
BAM="" #space-separated list of bam files
OUT="" #name of merged bam file ending in .bam

bamtools merge \
-in ${BAM} \
-out ${OUT}
```

__Step 5__. We will be using [freebayes](https://github.com/ekg/freebayes/blob/master/README.md) to call Single Nucleotide Polymorphisms (SNPs). Specifically, `freebayes-parallel` to speed up SNP calling by breaking our reference genome assembly into regions. The default parameters of `freebayes` are very liberal, so I recommend increasing their stringency. Below are parameters I typically change (->) from their default value, but you should choose values that best suits your data:

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
