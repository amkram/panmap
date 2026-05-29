# Multi-PanMAT screening for panmap (`--meta`)

**Problem.** Obelisk (and viroid) references can't be a single PanMAN — their clades share no
alignable coordinate frame (a global MSA of all 5,139 obelisks is 12,366 cols / 92% gaps). So you
need one PanMAN *per clade* (175 here). But screening a sample against 175 separate panmap indexes
re-seeds the reads 175× — read-seeding is ~95% of `--meta` runtime — so cost scales ~175×.

**Solution (implemented here, no changes to the scorer / mgsr.cpp).** The MGSR lite index is a flat
node list + global `seedHashes` + per-node `seedDeltaIndices`; `MgsrLiteTree::initialize` rebuilds
the tree purely from `parentIndex`, needing only that node[0] is a root. So multiple clades can live
in ONE index as a forest under a synthetic super-root. Two small standalone tools do this:

- **`merge_indexes OUT.idx IN1.idx ...`** (`src/merge_indexes.cpp`) — merges N single-tree MGSR lite
  `.idx` files into one combined `LiteIndex`: prepends a zero-seed `__superroot__` at node 0, appends
  each clade's nodes (parentIndex offset by the clade's node base; clade root → super-root),
  concatenates `seedHashes`, and shifts each clade's `seedDeltaIndices` by its seed base. Node ids are
  namespaced `<clade>::<id>` (internal ids like `node_1` collide across clades; `initialize` exit(1)s
  on duplicates). `seedChange*`/`nodeChangeOffsets` left empty (lite uses the seedInfos path).
- **`combine_panmans OUT.panman IN1.panman ...`** (`src/combine_panmans.cpp`) — merges N single-tree
  `.panman` into one multi-PanMAT `TreeGroup` (disconnected, 0 complex mutations) with the same id
  namespacing; needed only to satisfy panmap's positional `<panman>` arg (the meta path uses `-i`).

**Usage**
```bash
# one-time: build a single-tree lite index per clade, then merge + combine
for p in panmans/clade_*.panman; do c=$(basename $p .panman); \
  panmap --meta --index-mgsr idx/$c.idx --stop index $p; done
merge_indexes all.idx idx/clade_*.idx
combine_panmans all.panman panmans/clade_*.panman
# screen any sample against ALL clades in ONE read-seeding pass:
panmap --meta -i all.idx all.panman reads.fq -o out      # -> out.mgsr.abundance.out (clade::node)
```

**Verified (macOS arm64, 1 thread, real SRR11790905 = 1.62M reads):**
- Equivalence: a clade's abundances in the merged run are byte-identical to its standalone run;
  unrelated clades get ~0 (no cross-talk); no duplicate-id crash.
- Scale: ONE merged pass over **174 clades = 20.0 s**, vs a single-clade run 19.0 s → 174 separate
  runs ≈ **55 min**. **~165× speedup**; read-seeding (~14 s) is paid once. Merged index = 16.9 MB,
  9,275 nodes, 963k seeds; merge step itself is ~0.1 s.

**`panmap --meta -i merged.idx` uses ALL trees** (the merged index is a forest under the
super-root), verified by all 174 clades receiving placements on a real run — not just the first.

**Why not `--index-mgsr` directly on a combined multi-tree PanMAN?** `--index-mgsr` indexes a
single tree. Previously it silently used `trees[0]` on a multi-tree PanMAN (a first-clade-only,
silently-wrong index); it now **errors and directs to `merge_indexes`**. Building one index across
trees by re-indexing a combined PanMAN is unreliable anyway — panman multi-tree *re-processing* is
incomplete (e.g. `panmanUtils -f` on a multi-tree file does not emit all tips), and re-indexing a
tree reloaded from a combined PanMAN can crash. Indexing each clade's own single-tree PanMAN
(robust) and merging the indexes avoids this entirely.

**Caveats.** Covers the buildable clades indexed; small/singleton clades are better served by a
flat all-genomes detection index (the panmap MGSR indexer rejects some tiny/degenerate panmans:
`firstScalarCoord != 0`, same class as `clade_045`'s `getCoordFromScalar`). Standard (non-`--meta`)
placement remains single-tree.
