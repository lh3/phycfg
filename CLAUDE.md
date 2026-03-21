# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git Rules

- **Never stage or commit untracked files** unless explicitly instructed by the user.

## Build Commands

```sh
make          # build the phycfg binary
make clean    # remove build artifacts
make asan=1   # build with AddressSanitizer enabled
```

The binary is built as `./phycfg`. There are no automated tests at this stage; test data lives in `test/`.

Input files are gzip-compressed — phycfg reads gzip'd formats directly (do not decompress test files).

```sh
./phycfg view test/CCNE.nhx.gz   # example usage
```

## Architecture

**phycfg** is a C99 command-line tool for phylogenetic tree analysis (version 0.0, early development).

### Command dispatch

`main.c` is the entry point and contains `main_view()`, `main_msaflt()`, `main_reroot()`, and `main_scfg()`. It dispatches subcommands:
- `view` → `main_view()` at the bottom of `main.c`; accepts `-l STR` (comma/space-separated leaf names or `@file`) to extract and print the minimal induced subtree over those leaves
- `msaflt` → `main_msaflt()`; reads a gzip'd FASTA MSA, infers residue type, encodes, filters columns, and writes decoded FASTA to stdout; accepts `-m INT` (min non-gap/non-ambiguous residues per column, default 1) and `-1`/`-2`/`-3` to select codon positions
- `reroot` → `main_reroot()`; reroots a tree and writes Newick to stdout; by default uses global midpoint rooting; with `-l STR` roots at the midpoint of the branch leading to the LCA of the listed leaves
- `scfg` → `main_scfg()` in `scfg.c`; reads tree and MSA, encodes, matches sequences to leaves; with `-n INT` runs EM then up to INT rounds of NNI topology search, prints final Newick to stdout; `-m INT` sets EM iterations per round (default 100), `-b INT` sets branch EM iterations (default 50)
- `version` → prints `PC_VERSION` from `phycfg.h`

When `kom_verbose >= 3` and the command succeeds, timing/resource info is printed to stderr.

### Library (`libphycfg.a`)

Object files archived: `kommon.o knhx.o tree.o io.o msa.o model.o sfunc.o scfg.o`. Linked with `main.o`.

- **`kommon.c`/`kommon.h`** — general-purpose utilities:
  - Memory macros: `kom_malloc`, `kom_calloc`, `kom_realloc`, `kom_grow` (dynamic array growth with 1.5× expansion)
  - String: `kstring_t` (length-tracked), `kom_strdup`, `kom_strndup`, `kom_sprintf_lite` (limited format: `%d`, `%ld`, `%u`, `%s`, `%c`)
  - Numerics: `kom_parse_num` (parses K/M/G suffixes), `kom_splitmix64` (fast RNG), `kom_u64todbl`
  - Timing/RSS: `kom_realtime()`, `kom_cputime()`, `kom_peakrss()`, `kom_percent_cpu()` (cross-platform: POSIX + Win32)
  - Residues: `kom_nt4_table` (A/C/G/T→0–3, else 4), `kom_comp_table`, `kom_revcomp()`, `kom_aa20_table` (20 AAs + X + \*→0–21, else 22), `kom_aa_i2c` (index-to-char string `"ARNDCQEGHILKMFPSTWYV*X"`)
  - Globals: `kom_verbose` (default 3), `kom_dbg`, `kom_dbg_flag`

- **`knhx.c`/`knhx.h`** — Newick/NHX parser (parsing only; no formatting or tree operations):
  - `kn_parse(nhx, &n, &max, &error, &end)` — parses a Newick/NHX string into a flat `calloc`'d array of `knhx1_t` nodes; error bits defined as `KNERR_*` constants
  - `kn_destroy(n, a)` — frees node array (name strings, child index arrays, the array itself)
  - `knhx1_t` fields: `parent`, `n` (child count), `child[]`, `name`, `d` (branch length, -1.0 if absent), `aux`
  - **Gotcha**: `kn_parse` does not initialise `child` for leaf nodes. Always guard with `if (a[i].n > 0)` before freeing or dereferencing `child`; `kn_destroy` handles this correctly.
  - Used only by `pc_tree_parse` in `tree.c`; not part of the public API

- **`tree.c`** — primary tree data structures and operations (declared in `phycfg.h`):
  - `pc_node_t` fields: `n_child`, `ftime` (post-order index), `seq_id` (index into `pc_msa_t::name`, -1 if unmatched), `d` (branch length, -1.0 if absent), `tmp` (int32_t temporary field), `name`, `parent`, `q` (pointer to `pc_scfg_data_t`, NULL until `pc_scfg_alloc`), `child[]` (flexible array)
  - `pc_tree_t` fields: `n_node`, `m` (alphabet size, set by `pc_tree_match_msa`), `root`, `node` (pointer array in post-order; `node[i]->ftime == i`)
  - `pc_tree_parse(str, &end)` — parse Newick/NHX string via `kn_parse`, convert to `pc_tree_t`
  - `pc_tree_sync(t)` — rebuild `t->node[]` and `ftime` fields from scratch (always frees and reallocates)
  - `pc_tree_mark_leaf(t, n, leaf, mark)` — set `mark[i]=1` for leaves whose names appear in `leaf[0..n-1]`; caller allocates `mark`
  - `pc_tree_reduce(t, mark)` — return a new `pc_tree_t` spanning only leaves with `mark != 0`; suppresses unary internals, accumulates branch lengths; caller must `pc_tree_destroy` result
  - `pc_tree_lca(t, mark)` — return ftime of the LCA of all marked nodes via post-order count propagation; returns -1 if no nodes marked
  - `pc_tree_mid_longest(t, &dist_to_mid)` — find node `p` whose incoming branch contains the diameter midpoint (O(n) post-order); writes distance from `p` to midpoint into `*dist_to_mid`; returns `p`'s ftime or -1
  - `pc_tree_rotate(t, xi)` — NNI rotation: given node x at ftime `xi`, transforms `((x,y)u,w)v` into `((w,y)u,x)v`; requires x's parent u and grandparent v to exist and v to be binary; calls `pc_tree_sync` on success; returns 0 on success, -1 if not possible
  - `pc_tree_reroot(t, nid, dist)` — place new root on the branch to node `nid` at distance `dist` from that node (pass `dist<0` for branch midpoint); then calls `pc_tree_sync`
  - `pc_tree_clone(t)` — deep copy of tree including nodes, names, and `m`; `node->q` (`pc_scfg_data_t`) is NOT cloned; caller must `pc_tree_destroy` the result
  - `pc_tree_format(t, &s, &max)` — format tree to Newick; reusable-buffer API (pass `NULL`/`0` first call; `s` and `max` updated in place); returns string length; caller frees `*s`
  - `pc_tree_destroy(t)` — free all nodes and the `pc_tree_t` itself; does NOT free `node->q`; call `pc_scfg_free` first if needed

- **`io.c`** — file I/O (uses `kseq.h` via `KSEQ_INIT(gzFile, gzread)`):
  - `pc_tree_read(fn)` — read a gzip-compressed NHX file; returns `pc_tree_t*` or NULL
  - `pc_msa_read(fn)` — read a gzip-compressed FASTA MSA; validates uniform length; stores plain ASCII in `msa->msa[pos][seq]` layout; returns `pc_msa_t*` or NULL
  - `pc_list_read(o, &n)` — parse comma/space-separated names from string, or from `@file`; returns `char**`

- **`msa.c`** — MSA operations (declared in `phycfg.h`):
  - `pc_msa_infer_rt(msa)` — infer `pc_restype_t` from letter frequencies: ≥50% A/C/G/T → `PC_RT_NT`; ≥80% standard AA letters → `PC_RT_AA`; else `PC_RT_UNKNOWN`
  - `pc_msa_encode(msa, rt)` — set `msa->rt = rt` and `msa->m`; encode ASCII in-place using `kom_nt4_table` (NT) or `kom_aa20_table` (AA); `-`/`.` → `PC_GAP_NT`/`PC_GAP_AA`; if `rt == PC_RT_UNKNOWN`, does nothing
  - `pc_msa_filter(msa, min_cnt)` — in-place column filter (requires prior encode); keeps columns where at least `min_cnt` sequences have a value `< msa->m`; frees dropped rows
  - `pc_msa_select_codon(msa, codon_flag)` — keep only specified codon positions; `codon_flag` bits 0/1/2 select 1st/2nd/3rd positions; operates on encoded CDS (`PC_RT_CODON`) MSA
  - `pc_msa_destroy(msa)` — free all name strings, row arrays, and the struct itself; NULL-safe
  - `pc_msa_t` fields: `len` (alignment length / number of columns), `n_seq`, `rt`, `m` (alphabet size: 4 NT / 20 AA / 256 unknown), `name` (sequence names), `msa` (`uint8_t**`, position-major: `msa[pos][seq]`)
  - Gap/ambiguous constants: `PC_GAP_NT`=5, `PC_GAP_AA`=23 (defined in `phycfg.h`)

- **`model.c`** — substitution model constraints and distance estimation (declared in `phycfg.h`):
  - `pc_model_t` enum: `PC_MD_UNDEF=-1`, `PC_MD_FULL=0` (unconstrained), `PC_MD_REV` (reversible/GTR), `PC_MD_TN93`
  - `pc_model_from_str(model_str)` — parse model name string; recognises `"full"/"FULL"/"."` → `PC_MD_FULL`, `"rev"/"GTR"/"gtr"` → `PC_MD_REV`, `"TN93"/"tn93"` → `PC_MD_TN93`; returns `PC_MD_UNDEF` on unknown input
  - `pc_model_df(md, m)` — degrees of freedom: `PC_MD_FULL` → `m*(m-1)`, `PC_MD_REV` → `(m-1)*(m-2)`, `PC_MD_TN93` → 3; returns -1 for undefined
  - `pc_model_lrt(md0, md1, m, lr)` — P-value under Wilks' theorem: `2*|lr|` as chi-square with `|df(md1)-df(md0)|` d.f.
  - `pc_model_BIC(md0, md1, m, len, lr)` — BIC difference: `2*|lr| - |Δdf|*log(len)` (positive means md0 preferred)
  - `pc_model_matrix(cnt, m, md, tmp)` — apply model constraint: `PC_MD_REV` symmetrises; `PC_MD_TN93` symmetrises then pools transversion off-diagonals by `π_i·π_j·tv` (asserts `m==4`); `PC_MD_FULL` copies unchanged
  - `pc_model_dist_TN93(cnt, &kR, &kY)` — TN93 branch length from 4×4 joint count matrix; writes κ_R and κ_Y; returns expected substitutions per site
  - `pc_model_dist(t, msa, md)` — calls `pc_scfg_alloc` + `pc_scfg_post_cnt` to fill `node->q->jc`, then `pc_model_dist_TN93` per non-root branch; root set to 0.0; asserts `md == PC_MD_TN93`

- **`scfg.c`** — SCFG algorithms and the `scfg` subcommand:
  - **Data layout**: `pc_scfg_data_t` (defined in `phycfg.h`) holds per-node arrays for all positions:
    `p[m×m]` (transition matrix), `jc[m×m]` (posterior joint count), `h[len]` (scaling factor),
    `alpha[len×m]` (α̃), `alpha2[len×m]` (α̃'), `beta[len×m]` (β̃); all in one `calloc` block via flexible array `x[]`
  - `pc_scfg_alloc(t, len)` — allocates `pc_scfg_data_t` at `node->q` for all nodes (skips already-allocated)
  - `pc_scfg_free(t)` — frees `node->q` for all nodes (in `pcpriv.h`)
  - `pc_scfg_init_par(t)` — init `node->q->p`: non-root gets JC matrix from `node->d` (falls back to 1e-3 if `d ≤ 0`); root gets flat `1/m` prior
  - `pc_scfg_inside(t, msa, pos)` — inside pass for column `pos`; writes `q->alpha[pos*m]`, `q->alpha2[pos*m]`, `q->h[pos]`; returns log P(column)
  - `pc_scfg_outside(t, pos)` — outside pass (after inside); writes `q->beta[pos*m]`; β̃(root,a)=q(a)/h_root; β̃(u,b)=(1/h_u)·Σ_a p(b|a)·β̃(par,a)·∏_sib α̃'(sib,a)
  - `pc_scfg_post_cnt(t, msa)` — E-step over all columns: zeros `q->jc`, accumulates normalised posterior branch counts; returns total log likelihood (in `pcpriv.h`)
  - `pc_scfg_em_all(t, msa, ct)` — one EM round: `pc_scfg_post_cnt` then M-step via `pc_model_matrix` + row-normalise into `q->p`
  - `pc_scfg_em1(m, len, ct, xp, yp, up, wp, vp, max_itr, p)` — 1-branch EM for topology `((x,y)u,w)v`; optimises `p[m×m]` for branch `u`; `xp=yp=NULL` for non-NNI case (uses stored `up->q->alpha`); stops at improvement < 1e-6
  - `pc_scfg_nni1(t, msa, ct, max_iter_br)` — NNI with 1-branch EM; runs global inside/outside first; for each eligible internal non-root node evaluates original + 2 NNI rotations; applies best improving move via `pc_tree_rotate`; updates only `up->q->p`; returns log-likelihood improvement
  - `pc_scfg_em5(m, len, ct, xp, yp, up, wp, vp, max_itr, q)` — 5-branch EM for topology `(((x,y)u,w)v,z)p`; `q[5×m×m]` for x/y/u/w/v; each iteration recomputes α̃'(x/y/w), α̃(u), α̃'(u), β̃(u) from stored α̃(x/y/w) and stored β̃(v) without re-running global inside/outside; v's branch is optimised only when `vp->parent != NULL` (v not root)
  - `pc_scfg_nni5(t, msa, ct, max_iter_br)` — NNI with 5-branch EM; runs global inside/outside first; evaluates original + 2 NNI rotations per eligible node; on best move updates x/y/u/w matrices (remapped to post-rotation nodes) and v's matrix (if v has parent) via `pc_tree_rotate`
  - `pc_scfg_model_cmp(t, msa, md0, md1, max_iter_br, diff)` — per-branch log-likelihood ratio `log(P(md0)/P(md1))`; writes to `diff[n_node]` (root entry = 0)

- **`sfunc.c`** — special functions (declared in `pcpriv.h`):
  - `kf_chi2_p(df, x)` — chi-square survival function P(X > x) for X ~ χ²(df); implemented via regularised incomplete gamma: `kf_gammaq(df/2, x/2)`

### Third-party headers

- **`kseq.h`** — macro-based streaming FASTA/Q + kstream parser from [klib](https://github.com/attractivechaos/klib) (MIT). Instantiated in `io.c` via `KSEQ_INIT(gzFile, gzread)`.
- **`khashl.h`** — open-addressing hash table from klib (MIT). Used in `tree.c` for `pc_tree_mark_leaf`.
- **`ketopt.h`** — command-line option parser. Used in `main.c` and `scfg.c`.

## Test Data

`test/CCNE.{mfa,nhx}.gz` — a cyclin E family dataset:
- `.mfa.gz` — multiple sequence alignment in FASTA format (gzip-compressed)
- `.nhx.gz` — phylogenetic tree in Newick/NHX format (gzip-compressed)

Keep these files compressed; phycfg reads gzip input directly via zlib.

## Planned Theory (`tex/scfg.tex`)

The math notes describe the algorithms phycfg is intended to implement, framed as a **Stochastic Context-Free Grammar (SCFG)** on a rooted ordered full binary tree. Key components:

- **Inside algorithm** (`α(v,a)`) — equivalent to Felsenstein's pruning algorithm; computes the probability of observed leaf residues given an ancestral state. Time complexity O(|V|·|Σ|²).
- **Outside algorithm** (`β(u,b)`) — computes the probability of everything outside a subtree given a node state, using the inside values.
- **Branch form** (`η_u(b|a)`) — a factored representation enabling O(|Σ|²) nearest-neighbor interchange (NNI) moves for efficient tree topology search without re-traversing subtrees.
- **Posteriors** (`x_{v,a}`, `x_{u,ab}`) — posterior probability of a residue at a node, and of a substitution `a→b` on a branch, computed from inside/outside values.
- **EM parameter estimation** — updates branch transition probabilities `p_u(b|a)` from sufficient statistics accumulated over alignment columns.
- **Numerical scaling** — scaled versions `α̃` and `β̃` with per-node scale factors `h(v)` to prevent underflow; `h(v)` is chosen so that `Σ_a α̃(v,a) = 1`.
