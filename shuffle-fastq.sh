#!/usr/bin/env bash
set -eu -o pipefail
in1=$1
out1=${2-$1}
in2=${3-}
out2=${4-${3-}}
paste <(if [[ ${in1##*.} = gz ]]; then gzip -cd "$in1"; else cat "$in1"; fi |paste - - - -) \
      <(if [[ -n ${in2-} ]] then if [[ ${in2##*.} = gz ]]; then gzip -cd "$in2"; else cat "$in2"; fi |paste - - - -; fi) | \
shuf | \
perl -F\\t -lane 'print for @F[0..3]; print STDERR for @F[4..7] if @F >= 8' \
  > >(if [[ ${out1##*.} = gz ]]; then pigz||gzip; else cat; fi >"$out1.$$.tmp") \
  2> >(if [[ -n ${out2-} ]]; then if [[ ${out2##*.} = gz ]]; then pigz||gzip; else cat; fi >"$out2.$$.tmp"; else cat >&2; fi)
mv -v "$out1"{.$$.tmp,}
if [[ -n ${out2-} ]]; then mv -v "$out2"{.$$.tmp,}; fi
