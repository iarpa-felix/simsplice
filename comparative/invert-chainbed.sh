#!/usr/bin/env bash
exec perl -F\\t -lane 'print join "\t",@F[4..7,0..3,8,10,9,11]' "$@"