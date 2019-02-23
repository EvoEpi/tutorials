# MUSCLE (MUltiple Sequence Comparison by Log- Expectation)

[MUSCLE](https://www.ncbi.nlm.nih.gov/pubmed/15034147) is computer software for multiple sequence alignment of protein and nucleotide sequences. `MUSCLE` is faster and more accurate than T-Coffee, MAFFT and CLUSTALW.

```bash
IN="" fasta file of sequences to align
OUT="" fasta file of aligned sequences

muscle \
-in ${IN} \
-out ${OUT}
```

# PASTA (Practical Alignment using SATe and Transitivity)

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
