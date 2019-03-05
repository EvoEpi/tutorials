# Multiple Sequence Alignment (MSA)

## MUSCLE (MUltiple Sequence Comparison by Log- Expectation)

[MUSCLE](https://www.ncbi.nlm.nih.gov/pubmed/15034147) is computer software for multiple sequence alignment of protein and nucleotide sequences. `MUSCLE` is faster and more accurate than T-Coffee, MAFFT and CLUSTALW.

```bash
IN="" fasta file of sequences to align
OUT="" fasta file of aligned sequences

muscle \
-in ${IN} \
-out ${OUT}
```

## PASTA (Practical Alignment using SATe and Transitivity)

[PASTA](https://www.ncbi.nlm.nih.gov/pubmed/25549288) is ideal for large-scale multiple sequence alignment estimation.

```bash
FASTA="" #fasta file of sequences to align
MEM="" #memory in mb
CPU="" #number of threads
LOG="" #log file

run_pasta.py \
-i ${FASTA} \
--max-mem-mb=${MEM} \
--num-cpus=${CPU} \
> ${LOG}
```
## PRANK

_[PRANK](http://wasabiapp.org/software/prank/) is a probabilistic multiple alignment program for DNA, codon and amino-acid sequences. It’s based on a novel algorithm that treats insertions correctly and avoids over-estimation of the number of deletion events. In addition, PRANK borrows ideas from maximum likelihood methods used in phylogenetics and correctly takes into account the evolutionary distances between sequences. Lastly, PRANK allows for defining a potential structure for sequences to be aligned and then, simultaneously with the alignment, predicts the locations of structural units in the sequences._

_`PRANK` can do translated alignments of protein-coding DNA sequences or align them using the codon model. Translation is selected with the options -translate (standard code) or -mttranslate (mitochondrial code), and the codon alignment with the option -codon. Using the example data input_dna.fas, the following commands make a translated alignment:_

```
prank -d=drm2.fna -o=drm2_translated -translate -F
```

_...and a codon alignment:_

```bash
prank -d=drm2.fna -o=drm2_codon -codon -F
```
