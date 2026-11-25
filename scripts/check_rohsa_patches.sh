#!/usr/bin/env bash
# Check that rohsa/upstream/ROHSA is the upstream base plus the patch series in
# docs/upstream-prs/rohsa/patches. Exit 0: identical. 1: they differ or a patch
# does not apply. 3: the base commit is not in this clone (needs full history).
#
# The base is the first vendored copy (ab4f299), identical to upstream ROHSA
# e36bc49. The compiled binary src/ROHSA that commit tracked is a build product
# here and is left out. To change the vendored Fortran, add a patch to the end
# of the series and rerun this script.
set -euo pipefail

BASE=ab4f299b7e470f6f556f07f5c9e9ade4c3769052
TREE=rohsa/upstream/ROHSA
PATCHES=docs/upstream-prs/rohsa/patches

root=$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)
if ! git -C "$root" cat-file -e "$BASE^{commit}" 2>/dev/null; then
  echo "base commit $BASE not found; fetch the full history" >&2
  exit 3
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/rohsa-patches.XXXXXX")
trap 'rm -rf "$work"' EXIT
mkdir "$work/base" "$work/vendored"

git -C "$root" archive "$BASE:$TREE" | tar -x -C "$work/base"
rm -f "$work/base/src/ROHSA"
for p in "$root/$PATCHES"/*.patch; do
  if ! (cd "$work/base" && git apply --whitespace=nowarn "$p"); then
    echo "patch does not apply: ${p#"$root"/}" >&2
    exit 1
  fi
done

# Tracked files only, so a local build (src/ROHSA, *.o, *.mod) is not compared.
git -C "$root" ls-files -z -- "$TREE" | tar -C "$root" --null -T - -cf - |
  tar -x -C "$work/vendored" --strip-components=3

if ! diff -ruN "$work/base" "$work/vendored"; then
  echo "base plus patches differs from $TREE (diff above: - is base plus patches, + is the vendored tree)" >&2
  exit 1
fi
echo "ok: base plus $(ls "$root/$PATCHES"/*.patch | wc -l) patches equals $TREE"
