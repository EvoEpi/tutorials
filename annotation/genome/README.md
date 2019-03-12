# Annotation

## Genome

Genome annotation is a multi-data, -program, and -step process. Other than a __genome assembly__ to annotate, __repeat annotations__ (easy enough to perform), __RNA-seq data__ (preferrably from a diversity of tissues and experimental conditions), and __protein fasta files__ from several closely related species (these are easier to come by, e.g., [ENSEMBL](https://useast.ensembl.org/index.html), [NCBI](https://www.ncbi.nlm.nih.gov/genome/), and [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). Four of the major players (in order of operation) in genome annotation are (1) `Trinity`, (2) `PASA`, (3) `RepeatMasker`, and (4) `MAKER`.

...under construction.

## Gene

__InterProScan.__ "_InterProScan is a tool that combines different protein signature recognition methods into one resource. The number of signature databases and their associated scanning tools, as well as the further refinement procedures, increases the complexity of the problem._" More details are at [EMBL-EBI](http://www.ebi.ac.uk/interpro/).

__-appl__ Optional, comma separated list of analyses. If this option is not set, ALL analyses will be run.  
__-cpu__ Optional, number of cores for inteproscan.  
__-f__ Optional, case-insensitive, comma separated list of output formats. Supported formats are TSV, XML, JSON, GFF3, HTML and SVG. Default for protein sequences are TSV, XML and GFF3, or for nucleotide
__-goterms__ Optional, switch on lookup of corresponding Gene Ontology annotation (IMPLIES -iprlookup option)  
__-i__ Optional, path to fasta file that should be loaded on Master startup. Alternatively, in CONVERT mode, the InterProScan 5 XML file to convert.  
__-iprlookup__ Also include lookup of corresponding InterPro annotation in the TSV and GFF3 output formats.  
__-o__ Optional explicit output file name (relative or absolute path)

```bash
APPL="" #comma separated list of analyses, e.g., TIGRFAM,CDD,PANTHER,Pfam
CPU="" #number of cpus
F="" #comma separated list of output formats
IN="" #in aa fasta file
OUT="" #filename of outfile

sh interproscan.sh \
-appl ${APPL} \
-cpu ${CPU} \
-f ${F} \
-goterms \
-i ${IN} \
-iprlookup \
-o ${OUT}
```
