# Data repositories and downloading data

## Genome repositories

__[CoGe](https://genomevolution.org/coge/).__ "_CoGe is a platform for performing Comparative Genomics research. It provides an open-ended network of interconnected tools to manage, analyze, and visualize next-gen data."_

__[ENSEMBL](https://useast.ensembl.org/index.html).__ "_Ensembl is a genome browser for vertebrate genomes that supports research in comparative genomics, evolution, sequence variation and transcriptional regulation. Ensembl annotate genes, computes multiple alignments, predicts regulatory function and collects disease data. Ensembl tools include BLAST, BLAT, BioMart and the Variant Effect Predictor (VEP) for all supported species._"

__[FlyBase](https://flybase.org/).__ A database for drosophila genetics and molecular biology.

__[i5k Workspace](https://i5k.nal.usda.gov/content/data-downloads).__ "_If you use genomics to study arthropods, you are an i5k member!_" Goals of the i5k are to: (i) "_Organize the sequencing and analysis of the genomes of 5,000 arthropod species_", (ii) "_Provide guidelines and best practices for arthropod genome projects and their data management_", (iii) "_Help exisiting and new arthropod genome projects to find the most appropriate repository for their needs_", and (iv) "_Grow a community around arthropod genomes that works towards improved sequencing, assembly, annotation, and data management standards._"

__[NCBI - Genome](https://www.ncbi.nlm.nih.gov/genome/).__ "_This resource organizes information on genomes including sequences, maps, chromosomes, assemblies, and annotations._"

__[Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html).__ "_Phytozome, the Plant Comparative Genomics portal of the Department of Energy's Joint Genome Institute, provides JGI users and the broader plant science community a hub for accessing, visualizing and analyzing JGI-sequenced plant genomes, as well as selected genomes and datasets that have been sequenced elsewhere. As of release v12.1.6, Phytozome hosts 93 assembled and annotated genomes, from 82 Viridiplantae species._"

__[Mycocosm](https://genome.jgi.doe.gov/mycocosm/home).__ Similar to Phytozome, but for fungi. Check if genome is publicly available ([green highlighted](https://genome.jgi.doe.gov/fungi/fungi.info.html))

__[Vectorbase](https://www.vectorbase.org/).__ "_VectorBase is a National Institute of Allergy and Infectious Diseases (NIAID) Bioinformatics Resource Center (BRC) providing genomic, phenotypic and population-centric data to the scientific community for invertebrate vectors of human pathogens._"

## Sequence respositories

__[Gene Expression Ombibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/).__ GEO is a public functional genomics data repository. It is very similar to the SRA (below) and often links between the two repositories are available.

All GEO records and raw data files are freely available for bulk download from the [FTP site](ftp://ftp.ncbi.nlm.nih.gov/geo/). The root directory (FTP site link) ftp://ftp.ncbi.nlm.nih.gov/geo/ has 4 main subdirectories corresponding to the types of GEO records:

datasets/
platforms/
samples/
series/

Each of those has range subdirectories to avoid browser timeouts. Range subdirectory name is created by replacing the three last digits of the accession with letters "nnn". For example:

ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/ contains record-specific subdirectories for GDS1001, GDS1002, ..., and GDS1995.

Individual record directory in turn has subdirectories named according to the type of data they contain. Detailed description for each type of record follows. Note that not all records have the data of all described types.

Please see the [README](ftp://ftp.ncbi.nlm.nih.gov/geo/README.txt) for details on directory structure and file formats.

__[Short Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra).__ The SRA stores raw sequencing data and alignment information from high-throughput sequencing platforms.

Downloading data from the SRA requires some _a priori_ knowledge on the accession and location. Beginning with a list of desired SRA data sets (e.g., a list of SRA Run accessions, “SRRs”), the exact download location for that data file can be determined as follows:

wget/FTP root: ftp://ftp-trace.ncbi.nih.gov

Remainder of path:

/sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra

Where

{SRR|ERR|DRR} should be either ‘SRR’, ‘ERR’, or ‘DRR’ and should match the prefix of the target .sra file

Examples:

Downloading SRR304976 by wget or FTP:

```bash
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR304/SRR304976/SRR304976.sra
```

[Phil Ewels](http://phil.ewels.co.uk/) has written a lovely program ([SRA explorer](https://ewels.github.io/sra-explorer/)) to collect SRA datasets and get a quick bash download script for either SRA files or FastQ files.

Extracting fastq file(s) from .sra file(s) requires using [SRA-toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software):

__I__ Append read id after spot id as 'accession.spot.readid' on defline

```bash
SRA="" sra file

fastq-dump -I ${SRA}
```
