#!/usr/bin/env bash
set -eu -o pipefail
if [[ ${1##*.} = gz ]]; then gzip -cd "$1" else cat "$1" fi | \
paste - - - - | \
sort -k1,1 | \
tr "\t" "\n" | \
perl -MDigest::SHA=sha256_hex -ple 's/^([\@\+])(.+)$/$1.sha256_hex($2)/e' | \
if [[ ${1##*.} = gz ]]; then gzip else cat fi \
>"${2-$1}.$$.tmp"
mv -v "${2-$1}.$$"{.tmp,}
