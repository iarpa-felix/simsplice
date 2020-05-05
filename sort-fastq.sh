#!/usr/bin/env bash
set -eu -o pipefail
if [[ ${1##*.} = gz ]]; then gzip -cd "$1" else cat "$1" fi | \
paste - - - - | \
sort -k1,1 | \
tr "\t" "\n" | \
if [[ ${1##*.} = gz ]]; then gzip else cat fi \
>"${2-$1}.$$.tmp"
mv -v "${2-$1}.$$"{.tmp,}
