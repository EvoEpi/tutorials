# Species tree estimation using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) and [ASTRAL](https://github.com/smirarab/ASTRAL)

Many thanks to [Karolina Heyduk](https://karohey.wixsite.com/home) and [Shawn Thomas](https://jlmlab.wixsite.com/jlmlab/shawn-thomas) for sharing portions of their pipeline.

__Step 1__. Estimate gene trees using `RAxML`.

Conserved, single-copy genes tend to work well for ultimately estimating a species tree. At least in flowering plants, there has been concerted effort to develop a universal probe set for several hundred conserved single-copy nuclear genes ([Johnson et al. 2018](https://www.ncbi.nlm.nih.gov/pubmed/30535394)).

__-f__ Select algorithm (a: rapid Bootstrap analysis and search for bestscoring ML tree in one program run)  
__-x__ Specify an integer number (random seed) and turn on rapid bootstrapping  
__-p__ Specify a random number seed for the parsimony inferences  
__-N__ Specify the number of alternative runs on distinct starting trees. In combination with the "b" option, this will invoke a multiple bootstrap analysis  
__-m__ Model of Binary (Morphological), Nucleotide, MultiState, or Amino Acid Substitution  
__-s__ Specify the name of the alignment data file in PHYLIP or FASTA format  
__-n__ Specifies the name of the output file  
__-w__ FULL (!) path to the directory into which RAxML shall write its output files

```bash
x="" #random seed for rapid bootstrapping (1237)
p="" #random seed for parsimony inferences (12345)
N="" #number of alternative runs (500)
m="" #model (GTRGAMMA)
s="" #sequence file
n="" #output file
w="" #full out directory path

mpirun raxmlHPC-MPI-AVX \
-f a \
-x ${x} \
-p ${p} \
-N ${N} \
-m ${m} \
-s ${s} \
-n ${n} \
-w ${w}
```

__Step 2__. Separately concatenate `bipartitions` files and `bootstrap` files.

```bash
cat RAxML_bipartitions.* > biparts.tre
realpath RAxML_bootstrap.* > boots.txt 
```

__Step 3__. Estimate a species tree using `ASTRAL`.

__-i__ Concatenated bipartitions files  
__-b__ Concatenated boostrap files  
__-r__ Number of bootstrap replicates to consider (needs to be less than number of bootstraps generated i.e. if 500 bootstraps generated in RAxML use 400 in ASTRAL  
__-o__ Filename of output tree

```bash
i="" #concatenated bipartitions files
b="" #concatenated boostrap files
r="" #number of bootstrap replicates
o="" #filename of output tree

astral.5.6.1.jar -i ${i} -b ${b} -r ${r} -o ${o}
```
