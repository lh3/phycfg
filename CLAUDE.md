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

`main.c` is the entry point. It dispatches subcommands to dedicated `main_<command>()` functions:
- `view` → `main_view()` in `view.c`; accepts `-l STR` (comma/space-separated leaf names or `@file`) to extract and print the minimal induced subtree over those leaves
- `version` → prints `PC_VERSION` from `phycfg.h`

When `kom_verbose >= 3` and the command succeeds, timing/resource info is printed to stderr.

### Library (`libphycfg.a`)

Two object files are archived into a static library, then linked with `view.o` and `main.o`:

- **`kommon.c`/`kommon.h`** — general-purpose utilities:
  - Memory macros: `kom_malloc`, `kom_calloc`, `kom_realloc`, `kom_grow` (dynamic array growth with 1.5× expansion)
  - String: `kstring_t` (length-tracked), `kom_strdup`, `kom_sprintf_lite` (limited format: `%d`, `%ld`, `%u`, `%s`, `%c`)
  - Numerics: `kom_parse_num` (parses K/M/G suffixes), `kom_splitmix64` (fast RNG), `kom_u64todbl`
  - Timing/RSS: `kom_realtime()`, `kom_cputime()`, `kom_peakrss()`, `kom_percent_cpu()` (cross-platform: POSIX + Win32)
  - DNA: `kom_nt4_table`, `kom_comp_table`, `kom_revcomp()`
  - Globals: `kom_verbose` (default 3), `kom_dbg`, `kom_dbg_flag`

- **`knhx.c`/`knhx.h`** — Newick/NHX phylogenetic tree parser:
  - `kn_parse(nhx, &n, &error)` — parses a Newick/NHX string into a flat array of `knhx1_t` nodes; error bits defined as `KNERR_*` constants
  - `kn_destroy(n, a)` — frees a node array (name strings, child index arrays, the array itself)
  - `kn_format(node, root, &s)` — formats a node array back to Newick string into a `kstring_t`
  - `kn_extract_marked(a0, n0, &n_out)` — extracts the minimal induced subtree spanning all leaves with `aux != 0`; unary internals are suppressed and their branch lengths accumulated into the surviving child; returns a freshly allocated array (caller must `kn_destroy` it) or NULL if no leaf is marked
  - `knhx1_t` fields: `parent`, `n` (child count), `child[]`, `name`, `d` (branch length, -1.0 if absent), `aux` (caller-defined tag; used by `kn_extract_marked` to identify marked leaves)
  - **Gotcha**: `kn_parse` does not initialise `child` for leaf nodes (only internal nodes get a `calloc`'d child array). Always guard with `if (a[i].n > 0)` before freeing or dereferencing `child`; `kn_destroy` handles this correctly.
  - Can be built standalone with `-DKNHX_MAIN` for testing

### Third-party headers

- **`kseq.h`** — macro-based streaming FASTA/Q parser from [klib](https://github.com/attractivechaos/klib) (MIT). Instantiated via `KSEQ_INIT(type_t, __read)` macro. Included in `view.c` but not yet used.

## Test Data

`test/TF101005.{mfa,nhx}.gz` — a transcription factor family dataset:
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
