#!/usr/bin/env bash
set -eu -o pipefail
if [[ ${1##*.} = gz ]]; then gzip -cd "$1"; else cat "$1"; fi | \
  perl -MDigest::SHA=sha256_hex -ple 'BEGIN{$s=shift} s/^(\@)(\S+).*$/$1.sha256_hex("$2.$s")/e if ($.-1) % 4 == 0; s/^\+.*$/+/ if ($.-1) % 4 == 2' "${3-}" | \
paste - - - - | \
sort -k1,1 | \
tr "\t" "\n" | \
if [[ ${1##*.} = gz ]]; then pigz||gzip; else cat; fi >"${2-$1}.$$.tmp"
mv -v "${2-$1}"{.$$.tmp,}
