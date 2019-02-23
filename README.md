# tutorials
Example workflows of various programs ranging from processing ChIP-seq, DNA-seq, RNA-seq (mRNA and sRNA), and WGBS data, and estimating gene and species phylogenies, testing for selection, and more.

## Phylogenetics

[RAxML (Randomized Axelerated Maximum Likelihood)](https://cme.h-its.org/exelixis/web/software/raxml/index.html) is a program for sequential and parallel Maximum Likelihood (ML) based inference of large phylogenetic trees. It can also be used for postanalyses of sets of phylogenetic trees, analyses of alignments and, evolutionary placement of short reads.

Step 1. Generate ${NR} ML trees on distinct starting trees and also print the tree with the best likelihood to a file called RAxML_bestTree.${ID}.

```bash
NP="" #number of processors
ALIGN="" #fasta or phylip aligned sequences
ID="" #identifier
MODEL="" #model of substitution
NR="" #number of replicates

raxmlHPC-PTHREADS-AVX \
-T ${NP} \
-s ${ALIGN} \
-n ${ID}-ML \
-m ${MODEL} \
-p 1234 \
-#${NR} \
> ${ID}-ML.log
```

Step 2. Get bootstrap support values for the best ML tree.

```bash
NP="" #number of processors
ALIGN="" #fasta or phylip aligned sequences
ID="" #identifier
MODEL="" #model of substitution
NR="" #number of replicates

raxmlHPC-PTHREADS-AVX \
-T ${NP} \
-s ${ALIGN} \
-n ${ID}-BS \
-m ${MODEL} \
-p 1234 \
-b 1234 \
-#${NR} \
> ${ID}-BS.log
```

Step 3. Having computed the bootstrap replicate trees that will be printed to a file called RAxML_bootstrap.${ID} we can now use them to draw bipartitions on the best ML tree as follows.

```bash
NP="" #number of processors
ALIGN="" #fasta or phylip aligned sequences
ID="" #identifier
MODEL="" #model of substitution
NR="" #number of replicates

raxmlHPC-PTHREADS-AVX \
-T ${NP} \
-m ${MODEL} \
-p 1234 \
-f b \
-t RAxML_bestTree.${ID}-ML \
-z RAxML_bootstrap.${ID}-BS \
-n ${ID}-TR
```
