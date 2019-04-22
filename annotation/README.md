# Annotation

## Genome

Genome annotation is a multi-data, -program, and -step process. Other than a __genome assembly__ to annotate, __repeat annotations__ (easy enough to perform), __RNA-seq data__ (preferrably from a diversity of tissues and experimental conditions), and __protein fasta files__ from several closely related species (these are easier to come by, e.g., [ENSEMBL](https://useast.ensembl.org/index.html), [NCBI](https://www.ncbi.nlm.nih.gov/genome/), and [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). Four of the major players (in order of operation) in genome annotation are (1) `Trinity`, (2) `PASA`, (3) `RepeatMasker`, and (4) `MAKER`.

## Gene

Annotating protein coding genes typically invovled classifying them into families and predicting domains and important sites. A commonly used, and my preferred, program is `InterProScan`.
