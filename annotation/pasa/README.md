# Step 2. Use [PASA](https://github.com/PASApipeline/PASApipeline) to build a comprehensive transcriptome database using _de novo_ and _genome-guided_ RNA-seq assemblies

__Step 2.1__. Concatenate the _de novo_ and _genome-guided_ `Trinity` assemblies into a single 'transcripts.fasta' file.

```bash
cat <de novo assembly> <genome-guided assembly > <transcripts.fasta>
```

__Step 2.2__. Create a file containing the list of transcript accessions that correspond to the `Trinity` _de novo_ assembly.

```bash
.../pasa/.../misc_utilities/accession_extractor.pl < <de novo assembly> > <tdn.accs>
```

__Step 2.3__. Run `PASA` using RNA-seq related options with '--TDN tdn.accs'.

```bash
CPU="" #number of CPUs
GENOME="" #genome assembly fasta file
TDN="" #Trinity de novo assembly
ACC="" #list of transcript accessions

Launch_PASA_pipeline.pl \
-c pasa.alignAssembly.config \
-C \
-r \
--PASACONF pasa.config \
--CPU ${CPU} \
-R \
-g ${GENOME} \
--ALIGNERS blat \
-t ${TDN} \
--TDN ${ACC} \
--transcribed_is_aligned_orient
```

__Step 2.4__. After completing the `PASA` alignment assembly, generate the comprehensive transcriptome database.

```bash
TSCRPT="" #concatenated assembly fasta
PID="" #min percent identity
PAL="" #min per aligned

.../pasa/.../scripts/build_comprehensive_transcriptome.dbi \
-c pasa.alignAssembly.config \
-t ${TSCRPT} \
--min_per_ID 95 \
--min_per_aligned 30
```
