# Step 3. Identifying repetitive content using [RepeatMasker](http://www.repeatmasker.org/)

__Step 3.1__. _de novo_ repeat identification using `RepeatModeler` from `RepeatMasker`.

```bash
SAMPLE="" #sample name
GENOME="" #genome assembly fasta file
CPU="" #number of CPUs
LOG="" #log filename

BuildDatabase \
-name ${SAMPLE} \
-engine ncbi \
${GENOME}

RepeatModeler \
-pa ${CPU} \
-engine ncbi \
-database ${SAMPLE} \
2>&1 | \
tee ${LOG}
```

__Step 3.2__. Full repeat annotation.

Identify repeats using `RepeatMasker` and the appropriate library from `Repbase`.

```bash
REPBASE="" #name of repbase library
CPU="" #number of CPUs
LIB="" #repbase fasta file
ABBR="" #species abbreviation, typically first three letters of genus and species
GENOME="" #genome assembly fasta file

mkdir ${REPBASE}_mask
RepeatMasker \
-pa ${CPU} \
-e ncbi \
-lib ${LIB} \
-dir ${ABBR}_mask \
${GENOME}
```

Multiple repeat annotations can be performed using different libraries. It is a good idea to rename the outputs after each round of repeat annotation so they are more representative of what they contain.

Results from each round of annotation must be analyzed together to produce the final/full repeat annotation

```bash
mkdir full_mask
gunzip ${REPBASE}_mask/*.cat.gz rep_mask1/*.cat.gz rep_mask2/*.cat.gz
cat ${REPBASE}_mask/*.cat rep_mask1/*.cat rep_mask2/*.cat > full_mask/${GENOME}.full_mask.cat
cd full_mask
ProcessRepeats \
-species ${REPBASE} \
full_mask/${GENOME}.full_mask.cat
```

In order to feed these repeats into `MAKER` properly, we must separate out the complex repeats.

```bash
rmOutToGFF3.pl \
full_mask/${GENOME}.full_mask.out \
> full_mask/${GENOME}.full_mask.out.gff3

#isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" ${GENOME}.full_mask.gff3 \
> ${GENOME}.full_mask.complex.gff3

#reformat to work with MAKER
cat ${GENOME}.full_mask.complex.gff3 | \
perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
> ${GENOME}.full_mask.complex.reformat.gff3
```
