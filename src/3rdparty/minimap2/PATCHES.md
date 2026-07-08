# Vendored minimap2 — local modifications

This is a fork of [minimap2](https://github.com/lh3/minimap2) **2.26 (r1175)**,
vendored so panmap can call its aligner as a library. `MM_VERSION` is set to
`2.26-panmap` so the provenance shows up in `minimap2 --version` and BAM `@PG`
lines.

Local changes relative to upstream `v2.26`, in full:

## 1. Expose internal functions for library use

panmap's aligner wrapper (`src/mm_align.c`) drives minimap2's mapping stages
directly, so these were changed from `static` to external linkage:

- `map.c` — `collect_seed_hits`, `collect_seed_hits_heap`, `chain_post`,
  `align_regs`.

## 2. Bug fix: guard NULL `r->p` in `mm_align1`

- `align.c` — a gap-fill segment can produce `n_cigar == 0` while not
  z-dropped, leaving `r->p` unallocated; the non-z-dropped branch then
  dereferenced it. Guarded the dereference. See the comment at the call site.
