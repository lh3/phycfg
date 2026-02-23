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

`main.c` is the entry point and also contains `main_view()` and `main_msaflt()`. `main_scfg()` lives in `scfg.c`. It dispatches subcommands:
- `view` → `main_view()` at the bottom of `main.c`; accepts `-l STR` (comma/space-separated leaf names or `@file`) to extract and print the minimal induced subtree over those leaves
- `msaflt` → `main_msaflt()`; reads a gzip'd FASTA MSA, infers residue type, encodes, filters columns, and writes decoded FASTA to stdout; accepts `-m INT` (min non-gap/non-ambiguous residues per column, default 1) and `-c` (treat as CDS, filter whole codons)
- `scfg` → `main_scfg()` in `scfg.c`; reads tree and MSA, encodes, matches sequences to leaves, runs EM to estimate branch transition matrices, prints per-iteration log likelihood
- `version` → prints `PC_VERSION` from `phycfg.h`

When `kom_verbose >= 3` and the command succeeds, timing/resource info is printed to stderr.

### Library (`libphycfg.a`)

Object files archived: `kommon.o knhx.o tree.o io.o msa.o scfg.o`. Linked with `main.o`.

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
  - `pc_node_t` fields: `n_child`, `ftime` (post-order index), `aux` (caller tag), `seq_id` (index into `pc_msa_t::name`, -1 if unmatched), `d` (branch length, -1.0 if absent), `name`, `parent`, `child[]` (flexible array)
  - `pc_tree_t` fields: `n_node`, `root`, `node` (pointer array in post-order; `node[i]->ftime == i`)
  - `pc_tree_parse(str, &end)` — parse Newick/NHX string via `kn_parse`, convert to `pc_tree_t`
  - `pc_tree_expand(root, node)` — post-order traversal; pass `node=NULL` to count, non-NULL to fill
  - `pc_tree_sync(t)` — rebuild `t->node[]` and `ftime` fields from scratch (always frees and reallocates)
  - `pc_tree_mark_leaf(t, n, leaf)` — set `aux=1` on leaves whose names appear in `leaf[0..n-1]`
  - `pc_tree_reduce(t)` — return a new `pc_tree_t` spanning only leaves with `aux != 0`; suppresses unary internals, accumulates branch lengths; caller must `pc_tree_destroy` result
  - `pc_tree_format(t, &s, &max)` — format tree to Newick; reusable-buffer API (pass `NULL`/`0` first call; `s` and `max` updated in place); returns string length; caller frees `*s`
  - `pc_tree_destroy(t)` — free all nodes and the `pc_tree_t` itself

- **`io.c`** — file I/O (uses `kseq.h` via `KSEQ_INIT(gzFile, gzread)`):
  - `pc_tree_read(fn)` — read a gzip-compressed NHX file; returns `pc_tree_t*` or NULL
  - `pc_msa_read(fn)` — read a gzip-compressed FASTA MSA; validates uniform length; stores plain ASCII in `msa->msa[pos][seq]` layout; returns `pc_msa_t*` or NULL
  - `pc_list_read(o, &n)` — parse comma/space-separated names from string, or from `@file`; returns `char**`

- **`msa.c`** — MSA operations (declared in `phycfg.h`):
  - `pc_msa_infer_rt(msa)` — infer `pc_restype_t` from letter frequencies: ≥50% A/C/G/T → `PC_RT_NT`; ≥80% standard AA letters → `PC_RT_AA`; else `PC_RT_UNKNOWN`
  - `pc_msa_encode(msa, rt)` — set `msa->rt = rt` and `msa->m`; encode ASCII in-place using `kom_nt4_table` (NT) or `kom_aa20_table` (AA); `-`/`.` → `PC_GAP_NT`/`PC_GAP_AA`; if `rt == PC_RT_UNKNOWN`, does nothing
  - `pc_msa_filter(msa, min_cnt, is_cds)` — in-place column filter (requires prior encode); keeps columns where at least `min_cnt` sequences have a value `< msa->m`; with `is_cds`, processes and keeps/discards positions as triplets; frees dropped rows
  - `pc_msa_destroy(msa)` — free all name strings, row arrays, and the struct itself; NULL-safe
  - `pc_msa_t` fields: `len` (alignment length / number of columns), `n_seq`, `rt`, `m` (alphabet size: 4 NT / 20 AA / 256 unknown), `name` (sequence names), `msa` (`uint8_t**`, position-major: `msa[pos][seq]`)
  - Gap/ambiguous constants: `PC_GAP_NT`=5, `PC_GAP_AA`=23 (defined in `phycfg.h`)

- **`scfg.c`** — SCFG algorithms and the `scfg` subcommand (functions declared in `phycfg.h` where externally visible):
  - `pc_scfg_t` fields: `h` (per-node scaling factor), `*alpha` (α̃), `*alpha2` (α̃'), `*beta` (β̃) — all pointers into a single flat allocation; NOT a flexible-array struct
  - `pc_scfg_new(n_node, m)` — allocates one contiguous block for an array of `n_node` `pc_scfg_t` headers followed by `3*n_node*m` doubles; sets each node's three pointers into the data region
  - `pc_mat2d_new(n_row, n_col)` — allocates a 2D `double**` with row pointers and data in one flat block; used for transition matrices
  - `pc_transmat_init(p, m, t)` — fills `p[k]` (m×m, row-major, `p[k][a*m+b] = P(b|a)`) for each node k; non-root nodes get JC model from `t->node[k]->d` (clamped to ≥1e-3); root node gets flat `1/m` for all entries (encodes q(a))
  - `pc_scfg_inside(t, p, msa, pos, sd)` — inside (Felsenstein) pass for column `pos`; requires binary tree (`assert n_child∈{0,2}`); requires `seq_id≥0` on all leaves; returns `Σ_v log h(v) + log(Σ_a α̃(root,a)·q(a))` = log P(column)
  - `pc_scfg_outside(t, p, m, sd)` — outside pass (must follow inside); initializes β̃(root,a) = p[root][a]/h_root = q(a)/h_root; propagates β̃ downward using α̃' from inside; returns void
  - `pc_scfg_em_basic(t, p, msa, sd)` — one EM round over the full MSA: E-step accumulates per-branch normalized sufficient statistics (including root/prior row); M-step renormalizes each transition row in-place; returns total log likelihood `Σ_pos log P(pos)`

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
