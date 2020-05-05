#!/usr/bin/env bash
set -eu -o pipefail
if [[ ${1##*.} = gz ]]; then gzip -cd "$1" else cat "$1" fi | \
paste - - - - | \
sort -k1,1 | \
tr "\t" "\n" | \
if [[ ${1##*.} = gz ]]; then gzip else cat fi \
>"$1.$$.tmp"
mv -v "$1.$$.tmp" "$1"