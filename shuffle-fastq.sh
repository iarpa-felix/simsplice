#!/usr/bin/env bash
set -eu -o pipefail
in1=${1-/dev/stdin}
in2=${2-}
paste <(if [[ ${in1##*.} = gz ]]; then gzip -cd "$in1"; else cat "$in1"; fi |paste - - - -) \
      <(if [[ -n ${in2-} ]] then if [[ ${in2##*.} = gz ]]; then gzip -cd "$in2"; else cat "$in2"; fi |paste - - - -; fi) | \
shuf | \
perl -F\\t -lane 'print for @F[0..3]; print STDERR for @F[4..7] if @F >= 8'
