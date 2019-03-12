# Step 1. Generate [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) _de novo_ and _genome-guided_ assemblies

__Step 1.1__. `Trinity` _de novo_ assembly.

__--seqType__ Type of reads (fa or fq)  
__--max_memory__ Suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc) provided in Gb of RAM  
__--CPU__ Number of CPUs to use  
__--no_version_check__ Don't run a network check to determine if software updates are available  
__--full_cleanup__ Only retain the Trinity fasta file  
__--normalize_reads__ Run in silico normalization of reads  
__--left__ Left reads, one or more file names (separated by commas, not spaces)  
__--right__ Right reads, one or more file names (separated by commas, not spaces)

```bash
SEQTYPE="" #fa or fq
MEM="" #suggested max memory
CPU="" #number of CPUs
LEFT="" #comma separated list of left (1) reads
RIGHT="" #comma separated list of right (2) reads

Trinity \
--seqType ${SEQTYPE} \
--max_memory ${MEM} \
--CPU ${CPU} \
--no_version_check \
--full_cleanup \
--normalize_reads \
--left ${LEFT} \
--right ${RIGHT}
```

__Step 1.2__. `Trinity` _genome-guided_ assembly.

`hisat2`:  
__-p__ Number of threads  
__--dta__ Report alignments tailored for transcript assemblers including `StringTie`  
__-x__ The basename of the index for the reference genome  
__-1__ Comma-separated list of files containing mate 1s (Alternatively, -U, comma-separated list of files containing unpaired reads)  
__-2__ Comma-separated list of files containing mate 2s  
__-S__ Out SAM file

`samtools`:  
__-@__ Number of threads  
__-o__ Output to FILE [stdout]

`Trinity`:  
__genome_guided_bam__ Provide path to coordinate-sorted bam file  
__genome_guided_max_intron__ Maximum allowed intron length (also maximum fragment span on genome)  
__--max_memory__ Suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc) provided in Gb of RAM  
__CPU__ Number of CPUs to use

```bash
CPU="" #number of processors
INDEX="" #basename of the index files to write
LEFT="" #comma separated list of left (1) reads
RIGHT="" #comma separated list of right (2) reads
SAM="" #output sam file
BAM="" #output bam file
BP="" #max allowed intron length
MEM="" #max memory

hisat2 \
-p ${CPU} \
--dta \
-x ${INDEX} \
-1 ${LEFT} \
-2 ${RIGHT} \
-S ${SAM}

samtools \
sort \
-@ ${CPU} \
-o ${BAM} \
${SAM}

Trinity \
--genome_guided_bam ${BAM} \
--genome_guided_max_intron ${BP} \
--max_memory ${MEM} \
--CPU ${CPU}
```
