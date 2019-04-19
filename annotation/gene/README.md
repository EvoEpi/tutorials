# InterProScan

"_InterProScan is a tool that combines different protein signature recognition methods into one resource. The number of signature databases and their associated scanning tools, as well as the further refinement procedures, increases the complexity of the problem._" More details are at [EMBL-EBI](http://www.ebi.ac.uk/interpro/).

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
