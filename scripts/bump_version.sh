#!/usr/bin/env bash
# Bump the panmap version in all files that track it.
# The runtime version is generated from CMakeLists.txt via src/version.h.in,
# so main.cpp never needs editing. Recipe files still embed the version
# because bioconda fetches a tarball whose URL contains it.
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

# sed -i differs between GNU and BSD/macOS; pass a backup suffix (portable on both) then remove it.
sed -i.bak -E "s/^(project\(panmap VERSION )[0-9]+\.[0-9]+\.[0-9]+(\))/\1${new_version}\2/" \
    "${repo_root}/CMakeLists.txt" \
    "${repo_root}/recipe/CMakeLists.txt"

sed -i.bak -E "s/^(\{% set version = \")[0-9]+\.[0-9]+\.[0-9]+(\" %\})/\1${new_version}\2/" \
    "${repo_root}/recipe/meta.yaml"

sed -i.bak -E "s/^(  version: \")[0-9]+\.[0-9]+\.[0-9]+(\")/\1${new_version}\2/" \
    "${repo_root}/recipe/recipe.yaml"

rm -f "${repo_root}/CMakeLists.txt.bak" "${repo_root}/recipe/CMakeLists.txt.bak" \
      "${repo_root}/recipe/meta.yaml.bak" "${repo_root}/recipe/recipe.yaml.bak"

echo "Bumped to ${new_version}."
echo "Remember to update the tarball sha256 in recipe/meta.yaml and recipe/recipe.yaml."
