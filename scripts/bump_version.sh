#!/usr/bin/env bash
# Bump the panmap version. Runtime version comes from CMakeLists.txt via
# src/version.h.in, so main.cpp needs no edits. The bioconda recipe lives in
# bioconda-recipes (not this repo) and is bumped there once a release is tagged.
#
# Usage: scripts/bump_version.sh X.Y.Z

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 X.Y.Z" >&2
    exit 1
fi

new_version="$1"

if ! [[ "$new_version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Version must be X.Y.Z (got: $new_version)" >&2
    exit 1
fi

repo_root="$(cd "$(dirname "$0")/.." && pwd)"

# sed -i differs between GNU and BSD/macOS; pass a backup suffix that works on both, then delete it.
sed -i.bak -E "s/^(project\(panmap VERSION )[0-9]+\.[0-9]+\.[0-9]+(\))/\1${new_version}\2/" \
    "${repo_root}/CMakeLists.txt"

# Illustrative BioContainer tags in the docs (keep the version part current; the
# --N build suffix is bioconda's and is left as-is).
sed -i.bak -E "s#(quay\.io/biocontainers/panmap:)[0-9]+\.[0-9]+\.[0-9]+(--)#\1${new_version}\2#g" \
    "${repo_root}/README.md" \
    "${repo_root}/docs/installation.md"

rm -f "${repo_root}/CMakeLists.txt.bak" "${repo_root}/README.md.bak" "${repo_root}/docs/installation.md.bak"

echo "Bumped to ${new_version}."
echo "After tagging v${new_version}, update the bioconda recipe's version + sha256 in bioconda-recipes/recipes/panmap/meta.yaml."
