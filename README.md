## Introduction

Phycfg is a ***proof-of-concept*** tool that demonstrates a novel approach to
modeling **rooted** phylogenetic trees within the maximum likelihood (ML)
framework. By applying [Stochastic Context-Free Grammar][scfg-wiki] (SCFG) to
tree structures, it enables flexible parameter estimation via the EM
algorithm. I derived the SCFG formulation in 2006 but did not publish it. To
the best of my knowledge, this work is still novel.

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

I derived the SCFG formation for phylogenetic trees when [TreeFam][tf] was
still my primary project at the Sanger Institute. The [date I put][old-date] in
the initial notes on SCFG was 2006-11-26. A couple of months later, I started
working on next-generation sequencing data analysis. Although I [gave a
talk][youtube] at a [workshop][workshop] in September 2007, I still did not
write any code by then.

I first [implemented][phycfg-old] the theory in 2018 to model cell lineage
trees using single-cell SNVs produced for the [META-CS paper][meta-cs]. The
code worked but it was sloppy and the result is non-interesting, so we did not
include it in the paper. Half a month after my last code commit, I started a
faculty job and left the SCFG formulation behind again.

I came back to SCFG in 2025 when [Kevin Hu][kevin], a brilliant sophomore in
Harvard at the time, worked with me to provide a more careful implementation.
His code demonstrated the convergence of EM and could reconstruct ancestral
sequences that are near identical to the [IQ-TREE implementation][iqtree-asr],
confirming the theoretical correctness of SCFG. Encouraged by the result, I
implemented the theory again. This is phycfg.

Richard Durbin told me that he maintains a list of unpublished work that is not
abandoned yet. The SCFG formulation is the oldest in this list. I hope we can
publish it some day.

[scfg-wiki]: https://en.wikipedia.org/wiki/Probabilistic_context-free_grammar
[sub-model]: https://iqtree.github.io/doc/Substitution-Models
[TN93]: https://en.wikipedia.org/wiki/Models_of_DNA_evolution#TN93_model_(Tamura_and_Nei_1993)
[GTR]: https://en.wikipedia.org/wiki/Models_of_DNA_evolution#GTR_model_(Tavar%C3%A9_1986)
[tf]: https://www.treefam.org/
[old-date]: https://github.com/lh3/phycfg/blob/master/tex-old/scfg.tex#L28-L30
[youtube]: https://www.youtube.com/watch?v=Vm0EwVghNlQ
[workshop]: https://www.newton.ac.uk/seminar/10904/
[meta-cs]: https://www.pnas.org/doi/abs/10.1073/pnas.2013106118
[phycfg-old]: https://github.com/lh3/phycfg-old
[kevin]: https://www.linkedin.com/in/kevin-hu-a35249234/
[iqtree-asr]: https://iqtree.github.io/doc/Command-Reference#ancestral-sequence-reconstruction
