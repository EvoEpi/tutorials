# ChIP-seq

ChIP-sequencing, also known as ChIP-seq, is a method used to analyze protein interactions with DNA. ChIP-seq combines chromatin immunoprecipitation (ChIP) with DNA sequencing to identify the binding sites of DNA-associated proteins. An important part of designing ChIP-Seq experiments is determining what controls to use for the experiment. Chromatin inputs serve as a good control for bias in chromatin fragmentation and variations in sequencing efficiency; additionally, they provide greater and more evenly distributed coverage of the genome ([Kidder et al. 2011](https://www.ncbi.nlm.nih.gov/pubmed/21934668)). Hence, the processing of ChIP-seq involves cleaning and mapping reads from a DNA interaction of interest to a control, and identifying binding sites (peaks) relative to a control.

__Step 1__. Clean-up reads using `Trimmomatic`. See `Trimmomatic` example below.

__Step 2__. Map reads using [bowtie](http://bowtie-bio.sourceforge.net/index.shtml).

__-S__ Print alignments in SAM format  
__-t__ Print the amount of wall-clock time taken by each phase  
__-p__ Launch \<int\> parallel search threads (default: 1)  
__-v__ Print verbose output (for debugging)  
__--best__ Make Bowtie guarantee that reported singleton alignments are "best" in terms of stratum and in terms of the quality values at the mismatched position(s)  
__--strata__ If many valid alignments exist and are reportable and they fall into more than one alignment "stratum", report only those alignments that fall into the best stratum  
__-m__ Suppress all alignments for a particular read or pair if more than \<int\> reportable alignments exist for it  

```bash
REF="" #genome reference fasta file
INDEX="" #basename of the index files to write
MOD_FQ="" #fastq filename for modification
MOD="" #modification identifier
INPUT_FQ="" #fastq filename for input control
INPUT="" #input identifier
NP="" #number of processers

bowtie-build \
${REF} \
${INDEX}

bowtie \
${INDEX} \
${MOD_FQ} \
-S ${MOD}.sam \
-t \
-p ${NP} \
-v 2 \
--best \
--strata \
-m 1

bowtie \
${INDEX} \
${INPUT_FQ} \
-S ${INPUT}.sam \
-t \
-p ${NP} \
-v 2 \
--best \
--strata \
-m 1
```

__Step 3__. Sort bam file using [SAMtools](http://www.htslib.org/).

__-O__ Format  
__-o__ Output to FILE \[stdout\]  
__-T__ Tmpprefix  
__-@__ Threads

```bash
samtools \
sort \
-O 'bam' \
-o ${MOD}_sorted.bam \
-T tmp \
-@ ${NP} \
${MOD}.sam

samtools \
sort \
-O 'bam' \
-o ${INPUT}_sorted.bam \
-T tmp \
-@ ${NP} \
${INPUT}.sam
```

__Step 4__. Remove/mark PCR duplicates using [Picard](https://broadinstitute.github.io/picard/).

```bash
java -jar picard.jar \
MarkDuplicates \
INPUT=${MOD}_sorted.bam \
OUTPUT=${MOD}_clean.bam \
METRICS_FILE=${MOD}_METRICS_FILE.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT

java -jar picard.jar \
MarkDuplicates \
INPUT=${INPUT}_sorted.bam \
OUTPUT=${INPUT}_clean.bam \
METRICS_FILE=${INPUT}_METRICS_FILE.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT
```

__Step 5__. Convert bam to bed using [bedtools](https://bedtools.readthedocs.io/en/latest/).

```bash
bedtools bamtobed \
-i ${MOD}_clean.bam \
> ${MOD}.bed

bedtools bamtobed \
-i ${INPUT}_clean.bam \
> ${INPUT}.bed
```

__Step 6__. Call peaks using [macs2](https://github.com/taoliu/MACS).

__-t__ Treatment FILENAME  
__-c__ Control FILENAME  
__-g__ It's the mappable genome size or effective genome size which is defined as the genome size which can be sequenced  
__--keep-dup__ It controls the MACS behavior towards duplicate tags at the exact same location -- the same coordination and the same strand  
__-n__ The name string of the experiment  
__--broad__ When this flag is on, MACS will try to composite broad regions in BED12 ( a gene-model-like format ) by putting nearby highly enriched regions into a broad region with loose cutoff  

```bash
GSIZE=""

macs2 callpeak \
-t ${MOD}.bed \
-c ${INPUT}.bed \
-g ${GSIZE} \
--keep-dup all \
-n ${MOD} \
--broad
```
