# tutorials
Example workflows of various programs ranging from processing ChIP-seq, DNA-seq, RNA-seq (mRNA and sRNA), and WGBS data, and estimating gene and species phylogenies, testing for selection, and more.

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
