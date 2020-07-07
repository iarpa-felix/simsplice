#!/usr/bin/env bash
set -eu -o pipefail
in1=$1
in2=$2
out1=${3-$1}
out2=${4-$2}
paste <(if [[ ${in1##*.} = gz ]]; then gzip -cd "$in1"; else cat "$in1"; fi |paste - - - -) \
      <(if [[ ${in2##*.} = gz ]]; then gzip -cd "$in2"; else cat "$in2"; fi |paste - - - -) | \
shuf | \
perl -F\\t -lane 'print for @F[0..3]; print STDERR for @F[4..7]' \
  > >(if [[ ${out1##*.} = gz ]]; then pigz||gzip; else cat; fi >"$out1.$$.tmp") \
  2> >(if [[ ${out2##*.} = gz ]]; then pigz||gzip; else cat; fi >"$out2.$$.tmp")
mv -v "$out1"{.$$.tmp,}
mv -v "$out2"{.$$.tmp,}
