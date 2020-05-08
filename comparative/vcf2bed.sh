#!/usr/bin/env bash
perl -F\\t -lane '
next if /^#/; 
die if $F[4]=~/,/; 
if ($F[0] ne $c) {
  $q=0; 
  $r=0; 
  $c=$F[0];
} 
$v++;
print join "\t",$F[0],($F[1]-1),($F[1]-1)+length($F[3]),"vcf_$v",0,"+";
print STDERR join "\t",$F[0],$r+(($F[1]-1)-$q),$r+(($F[1]-1)-$q)+length($F[4]),"vcf_$v",0,"+";
$r+=(($F[1]-1)-$q)+length($F[4]); 
$q=($F[1]-1)+length($F[3]); 
' "$@"
