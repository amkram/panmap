# Troubleshooting

## Installation and build

**Build fails looking for `zlib.h`, Cap'n Proto, or htslib headers.**
Build inside the conda environment, which supplies the include and library paths
the build needs:

```bash
conda env create -f environment.yml && conda activate panmap
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
```

Building without conda? Install the system packages listed under
[Installation → Without conda](installation.md#without-conda-system-packages).

**`conda install` can't find the package.** Include both channels, in order:
`conda install -c conda-forge -c bioconda panmap`.

## Running panmap

**Results look stale after changing seed parameters.** The index is cached next
to the PanMAN and reused on later runs. After changing `-k/-s/-l` or other
index-affecting options, force a rebuild:

```bash
panmap ref.panman reads.fq --reindex -o sample
```

**Placement lands on an unexpected node.** By default, any run past the
placement stage restricts placement to leaf nodes (`--force-leaf` is auto-enabled
unless `--stop place`). To place on internal nodes, stop at placement:

```bash
panmap ref.panman reads.fq --stop place -o sample
```

If placement is ambiguous, `--refine` re-ranks the top candidates by alignment
score (slower but more accurate).

**The VCF is empty.** For a sample that matches its closest reference exactly,
zero variants is correct. Otherwise, check that reads actually reach the genotype
stage (`--stop genotype` or the default `consensus`) and that coverage is
sufficient.

**Only some outputs were written.** panmap stops at the stage given by `--stop`
(`index → place → align → genotype → consensus`). To get a BAM and VCF, run at
least `--stop genotype`; for a consensus, use the default.

**Out of memory or very slow on a large PanMAN.** Placement scales with the
number of nodes. Give panmap more threads (`-t`), and note that the metagenomic
EM (`--meta`) is heavier than single-sample placement. See
[System requirements](installation.md) for the ~8 GB baseline.

## Getting more detail

- `panmap --help-all` lists every option (see the [CLI Reference](cli-reference.md)).
- `-v/--verbose` prints per-stage detail; `-q/--quiet` limits output to errors.
- Reference outputs for the bundled examples are in
  [`examples/expected/`](https://github.com/amkram/panmap/tree/main/examples/expected).

Report other issues at
[github.com/amkram/panmap/issues](https://github.com/amkram/panmap/issues).
