# RNA-seq (mRNA)

A popular toolset used for analyzing RNA-seq data is the tuxedo suite, which consists of [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) and [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/). The suite provided a start to finish pipeline that allowed users to map reads, assemble transcripts, and perform differential expression analyses. A newer "tuxedo suite" has been developed and is made up of three tools: `HISAT`, `StringTie`, and `Ballgown`. [Pertea et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27560171) provides a summary of the new suite as well as a tutorial.

__Step 1__. Mapping is performed using `HISAT2` and usually the first step, prior to mapping, is to create an index of the reference genome. Caution: This step takes a lot of memory. _If you use --snp, --ss, and/or --exon, [hisat2-build](https://ccb.jhu.edu/software/hisat2/manual.shtml) will need about **200GB RAM** for the **human genome** size as index building involves a graph construction. Otherwise, you will be able to build an index on your desktop with 8GB RAM_.

`hisat2-build`:  
__-p__ Number of threads  
__-ss__ Note this option should be used with the following --exon option. Provide a list of splice sites  
__-exon__ Note this option should be used with the above --ss option. Provide a list of exons  

`hisat2`:  
__-p__ Number of threads  
__--dta__ Report alignments tailored for transcript assemblers including `StringTie`  
__-x__ The basename of the index for the reference genome  
__-U__ Comma-separated list of files containing unpaired reads to be aligned. Alternatively, use `-1` (comma-separated list of files containing mate 1s) and `-2` (comma-separated list of files containing mate 2s) for paired-end reads.  
__-S__ File to write SAM alignments to  

`samtools`:  
__-@__ Number of threads  
__-o__ Output to FILE [stdout]  

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

hisat2 \
-p ${NP} \
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

__Step 2__. Now we need to assemble the mapped reads into transcripts. `StringTie` can assemble transcripts with or without annotation. With annotation:

__-l__ Sets <label> as the prefix for the name of the output transcripts  
__-p__ Number of threads
__-G__ Use the reference annotation file (in GTF or GFF3 format) to guide the assembly process  
__-o__ Sets the name of the output GTF file where StringTie will write the assembled transcripts  
__-e__ Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option (requires -G, recommended for -B/-b)  
__-A__ Gene abundances will be reported (tab delimited format) in the output file with the given name  

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
