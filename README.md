## Introduction

Phycfg is a **proof-of-concept** tool that demonstrates a new way to model
**rooted** phylogenetic trees in the maximum likelihood (ML) framework. It
applies [Stochastic Context-Free Grammar][scfg-wiki] (SCFG) to trees and
enables flexible parameter estimation with the EM algorithm. I derived the SCFG
formulation in 2006 but did not publish it. To the best of my knowledge, this
work is still novel.

You can find the theory of SCFG in the [tex](tex) directory. Older notes in
[tex-old](tex-old) provide more background.

## Functionality

In its current implementation, phycfg considers branches are independent of
each other in that there are no shared parameters between branches. It
internally infers the posterior joint count matrix at each branch and uses this
matrix for EM-based parameter estimation and local topology search.

Although the implementation of SCFG is alphabet agnostic, **only nucleotide
sequences are fully supported** for now. Phycfg implements three [nucleotide
substitution models][sub-model]:
* `FULL`: time irreversible model without restrictions (12 free parameters per branch)
* `GTR`: [General Time Reversible][GTR] model
* `TN93`: by [Tamura and Nei (1993)][TN93]
Adding the support of binary alphabet should be trivial but I do not have test
data to evaluate. Phycfg would not work well with amino acid sequences due to
the large number of free parameters.

### Re-estimate branch lengths with TN93

Different from conventional algorithms, phycfg assumes each branch is
associated with a different transition/transversion ratio and different
nucleotide frequencies. It uses a moment estimator to infer branch lengths. For
example, with
```sh
phycfg scfg -e100 test/SMC1.nhx.gz test/SMC1.mfa.gz
```
phycfg will apply 100 EM iterations, print log likelihood to stderr and output
a tree with re-estimated branch lengths under the TN93 model. You can see the
likelihood increasing with each iteration. Note that the moment estimator may
fail under long evolutionary distance.

### Nearest neighbor interchange (NNI)

NNI is a foundamental operation in ML tree building. Phycfg does not provide
the full procedure to find the globally optimal tree toplogy, which is highly
challenging, but it can perform NNIs to search optimal trees locally around
input trees:
```sh
phycfg scfg -n100 test/SMC1.nhx.gz test/SMC1.mfa.gz
```
In each round, phycfg computes the approximate likelihood of all topologies
within one NNI, selects the topology giving the largest likelihood increase and
applies multiple rounds of EM to update other branches. This process stops if
there is no locally better topology.

### Model comparison at each branch

Phycfg can compare two models **at each branch**:
```sh
phycfg scfg -t TN93 -m GTR test/CCNE.nhx.gz test/CCNE.mfa.gz
```
On `CD`-lines, it reports the log likelihood ratio, P-value of the likelihood
ratio test and BIC difference at the last three columns.

## History

[scfg-wiki]: https://en.wikipedia.org/wiki/Probabilistic_context-free_grammar
[sub-model]: https://iqtree.github.io/doc/Substitution-Models
[TN93]: https://en.wikipedia.org/wiki/Models_of_DNA_evolution#TN93_model_(Tamura_and_Nei_1993)
[GTR]: https://en.wikipedia.org/wiki/Models_of_DNA_evolution#GTR_model_(Tavar%C3%A9_1986)
