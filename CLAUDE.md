# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

`main.c` is the entry point and also contains `main_view()`, `main_msaflt()`, and `main_reroot()`. `main_scfg()` lives in `scfg.c`. It dispatches subcommands:
- `view` → `main_view()` at the bottom of `main.c`; accepts `-l STR` (comma/space-separated leaf names or `@file`) to extract and print the minimal induced subtree over those leaves
- `msaflt` → `main_msaflt()`; reads a gzip'd FASTA MSA, infers residue type, encodes, filters columns, and writes decoded FASTA to stdout; accepts `-m INT` (min non-gap/non-ambiguous residues per column, default 1) and `-1`/`-2`/`-3` to select codon positions
- `reroot` → `main_reroot()`; reroots a tree and writes Newick to stdout; by default uses global midpoint rooting; with `-l STR` roots at the midpoint of the branch leading to the LCA of the listed leaves
- `scfg` → `main_scfg()` in `scfg.c`; reads tree and MSA, encodes, matches sequences to leaves; with `-n INT` runs EM then up to INT rounds of NNI topology search, prints final Newick to stdout; `-m INT` sets EM iterations per round (default 100), `-b INT` sets branch EM iterations (default 50)
- `version` → prints `PC_VERSION` from `phycfg.h`

When `kom_verbose >= 3` and the command succeeds, timing/resource info is printed to stderr.

### Library (`libphycfg.a`)

Object files archived: `kommon.o knhx.o tree.o io.o msa.o model.o scfg.o`. Linked with `main.o`.

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
  - `pc_node_t` fields: `n_child`, `ftime` (post-order index), `seq_id` (index into `pc_msa_t::name`, -1 if unmatched), `d` (branch length, -1.0 if absent), `name`, `parent`, `child[]` (flexible array)
  - `pc_tree_t` fields: `n_node`, `m` (alphabet size, set by `pc_tree_match_msa`), `root`, `node` (pointer array in post-order; `node[i]->ftime == i`), `p` (owned `double[n_node×m×m]` transition matrix, NULL until `pc_transmat_init` is called)
  - `pc_tree_parse(str, &end)` — parse Newick/NHX string via `kn_parse`, convert to `pc_tree_t`
  - `pc_tree_expand(root, node)` — post-order traversal; pass `node=NULL` to count, non-NULL to fill
  - `pc_tree_sync(t)` — rebuild `t->node[]` and `ftime` fields from scratch (always frees and reallocates)
  - `pc_tree_mark_leaf(t, n, leaf, mark)` — set `mark[i]=1` for leaves whose names appear in `leaf[0..n-1]`; caller allocates `mark`
  - `pc_tree_reduce(t, mark)` — return a new `pc_tree_t` spanning only leaves with `mark != 0`; suppresses unary internals, accumulates branch lengths; caller must `pc_tree_destroy` result
  - `pc_tree_lca(t, mark)` — return ftime of the LCA of all marked nodes via post-order count propagation; returns -1 if no nodes marked
  - `pc_tree_mid_longest(t, &dist_to_mid)` — find node `p` whose incoming branch contains the diameter midpoint (O(n) post-order); writes distance from `p` to midpoint into `*dist_to_mid`; returns `p`'s ftime or -1
  - `pc_tree_reroot_core(root, p0, dist)` — place new root on the branch to `p0` at distance `dist` from `p0` (pass `dist<0` to use branch midpoint); reverses edges up to old root, suppresses old root if binary; returns new root node; if `p0==root` returns root unchanged; if `dist>p0->d` returns NULL
  - `pc_tree_reroot(t, nid, dist)` — wrapper: calls `pc_tree_reroot_core` then `pc_tree_sync`
  - `pc_tree_clone(t)` — deep copy of tree including nodes, names, `m`, and `p` (if non-NULL); caller must `pc_tree_destroy` the result
  - `pc_tree_format(t, &s, &max)` — format tree to Newick; reusable-buffer API (pass `NULL`/`0` first call; `s` and `max` updated in place); returns string length; caller frees `*s`
  - `pc_tree_destroy(t)` — free all nodes, `t->p`, and the `pc_tree_t` itself

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
  - `pc_model_t` enum: `PC_MD_ERR`, `PC_MD_NULL` (unconstrained), `PC_MD_REV` (reversible/GTR), `PC_MD_TN93`
  - `pc_model_from_str(model_str)` — parse model name string to `pc_model_t`; recognises `"null"/"NULL"`, `"rev"/"GTR"/"gtr"`, `"TN93"/"tn93"`; returns `PC_MD_ERR` on unknown input
  - `pc_model_matrix(cnt, m, md, tmp)` — apply model constraint to raw posterior count matrix `cnt` into `tmp`; `PC_MD_REV` symmetrises; `PC_MD_TN93` symmetrises then pools transversion rates by `π_i·π_j·tv` (nucleotide only, asserts `m==4`); `PC_MD_NULL` copies unchanged
  - `pc_model_dist_TN93(cnt, &kR, &kY)` — compute TN93 branch length from a 4×4 joint count matrix; writes transition rate ratios κ_R and κ_Y; returns branch length scaled by the expected number of substitutions per site
  - `pc_model_dist(t, msa, md)` — compute branch lengths for all branches: calls `pc_transmat_init`, allocates `pc_scfg_buf_t`, runs `pc_scfg_post_cnt` to get posterior joint counts, then calls `pc_model_dist_TN93` per non-root branch and stores result in `node->d`; root branch set to 0.0; currently asserts `md == PC_MD_TN93`

- **`scfg.c`** — SCFG algorithms and the `scfg` subcommand (functions declared in `phycfg.h` where externally visible):
  - `pc_scfg_buf_t` fields: `h` (per-node scaling factor), `*alpha` (α̃), `*alpha2` (α̃'), `*beta` (β̃) — all pointers into a single flat allocation; NOT a flexible-array struct
  - `pc_nni_t` fields: `rotation` (0/1/2), `u` (node ftime), `loglk`, `p[]` (flexible m×m transition matrix) — result of `pc_scfg_em_branch`; caller must free
  - `pc_scfg_buf_new(n_node, m)` — allocates one contiguous block for an array of `n_node` `pc_scfg_buf_t` headers followed by `3*n_node*m` doubles; sets each node's three pointers into the data region
  - `pc_transmat_init(t)` — allocates `t->p` if NULL, then fills `t->p + k*m*m` (m×m, row-major, `p[k*m*m + a*m+b] = P(b|a)`) for each node k; non-root nodes get JC model from `t->node[k]->d` (clamped to ≥1e-3); root node gets flat `1/m` (encodes q(a))
  - `pc_scfg_inside(t, msa, pos, sd)` — inside (Felsenstein) pass for column `pos`; requires binary tree (`assert n_child∈{0,2}`); requires `seq_id≥0` on all leaves; returns `Σ_v log h(v) + log(Σ_a α̃(root,a)·q(a))` = log P(column)
  - `pc_scfg_outside(t, sd)` — outside pass (must follow inside); initializes β̃(root,a) = q(a)/h_root; propagates β̃ downward using α̃' from inside; returns void
  - `pc_scfg_eta(t, sd, eta)` — compute η̃[n_node*m*m] from inside/outside values in `sd`; for non-root u: `η̃(u,b|a) = β̃(par,a)·α̃(u,b)·∏ₖ α̃'(sibₖ,a)`; not defined at root (loop stops before root)
  - `pc_scfg_eta_nni(t, sd, eta)` — compute η̃ for all three NNI rotations; `eta` has shape `(n_node, 3, m, m)`; only written for eligible nodes (internal, non-root); rotation 0 = original `((x,y)u,w)v`, 1 = `((w,y)u,x)v`, 2 = `((x,w)u,y)v`
  - `pc_scfg_em_branch(t, ct, len, eta, u, rotation, max_itr)` — run EM on the m×m transition matrix of branch `u` (seeded from `t->p + u*m*m`) under a given NNI rotation using precomputed `eta[len][n_node*3*m*m]`; returns allocated `pc_nni_t*` (or NULL for leaves/root); caller must free
  - `pc_scfg_post_cnt(t, msa, sd, cnt)` — E-step over all MSA columns: zeros `cnt[n_node*m*m]`, runs inside/outside per column, accumulates normalized posterior branch counts; returns total log likelihood
  - `pc_scfg_em(t, msa, ct, sd)` — one EM round: calls `pc_scfg_post_cnt` then M-step via `pc_model_matrix` + row-normalise `t->p` in-place; returns total log likelihood
  - `pc_scfg_nni(t, msa, ct, max_iter_br)` — one NNI round: computes η̃ for all rotations, runs `pc_scfg_em_branch` for all nodes, applies the best improving move (if any) and rearranges `t->p` to match the new post-order; returns improvement (0 = no move)

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
