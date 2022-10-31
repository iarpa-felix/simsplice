#!/usr/bin/env bash
perl -F\\t -lane '
next if /^#/; 
die if $F[4]=~/,/; 
if ($F[0] ne $c) {
  $q=0; 
  $r=0; 
  $c=$F[0];
} 
print join "\t",$F[0],$r+(($F[1]-1)-$q),$r+(($F[1]-1)-$q)+length($F[4]),"+",$F[0],$F[1]-1,($F[1]-1)+length($F[3]),"+",1,0,0,length($F[3]); 
print STDERR join "\t",$F[0],$F[1]-1,($F[1]-1)+length($F[3]),"+",$F[0],$r+(($F[1]-1)-$q),$r+(($F[1]-1)-$q)+length($F[4]),"+",1,0,0,length($F[4]); 
$r+=(($F[1]-1)-$q)+length($F[4]); 
$q=($F[1]-1)+length($F[3]); 
' "$@"
