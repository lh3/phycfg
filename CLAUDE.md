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
./phycfg view test/TF101005.nhx.gz   # example usage
```

## Architecture

**phycfg** is a C99 command-line tool for phylogenetic tree analysis (version 0.0, early development).

### Command dispatch

`main.c` is the entry point and also contains `main_view()`. It dispatches subcommands:
- `view` → `main_view()` at the bottom of `main.c`; accepts `-l STR` (comma/space-separated leaf names or `@file`) to extract and print the minimal induced subtree over those leaves
- `version` → prints `PC_VERSION` from `phycfg.h`

When `kom_verbose >= 3` and the command succeeds, timing/resource info is printed to stderr.

### Library (`libphycfg.a`)

Object files archived: `kommon.o knhx.o tree.o io.o msa.o`. Linked with `main.o`.

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
  - `pc_node_t` fields: `n_child`, `ftime` (post-order index), `aux` (caller tag), `d` (branch length, -1.0 if absent), `name`, `ptr`, `parent`, `child[]` (flexible array)
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
  - `pc_msa_encode(msa, rt)` — set `msa->rt = rt` then encode ASCII in-place using `kom_nt4_table` (NT) or `kom_aa20_table` (AA); if `rt == PC_RT_UNKNOWN`, does nothing
  - `pc_msa_t` fields: `n_pos`, `n_seq`, `rt`, `name` (sequence names), `msa` (`uint8_t**`, position-major layout)

### Third-party headers

- **`kseq.h`** — macro-based streaming FASTA/Q + kstream parser from [klib](https://github.com/attractivechaos/klib) (MIT). Instantiated in `io.c` via `KSEQ_INIT(gzFile, gzread)`.
- **`khashl.h`** — open-addressing hash table from klib (MIT). Used in `tree.c` for `pc_tree_mark_leaf`.
- **`ketopt.h`** — command-line option parser. Used in `main.c`.

## Test Data

`test/{CCNE,TF101005}.{mfa,nhx}.gz` — a transcription factor family dataset:
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
