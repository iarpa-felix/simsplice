#!/usr/bin/env bash
perl -F\\t -lane '
if (/^##contig=<([^>]+)>$/) {
  $l={map {split /=/,$_,2} split /,/,$1}; 
  $c{$l->{ID}}=$l->{length}
} 
next if /^#/; 
die if $F[4]=~/,/; 
if ($F[0] ne $c) {
  print join "\t",$c,$r,$r+($c{$c}-$q),"+",$c,$q,$c{$c},"+",1,0,0,$c{$c}-$q if defined $c; 
  print STDERR join "\t",$c,$q,$c{$c},"+",$c,$r,$r+($c{$c}-$q),"+",1,0,0,$c{$c}-$q if defined $c; 
  $q=0; 
  $r=0; 
  $c=$F[0];
} 
print join "\t",$F[0],$r,$r+(($F[1]-1)-$q),"+",$F[0],$q,($F[1]-1),"+",1,0,0,(($F[1]-1)-$q); 
print STDERR join "\t",$F[0],$q,($F[1]-1),"+",$F[0],$r,$r+(($F[1]-1)-$q),"+",1,0,0,(($F[1]-1)-$q); 
$r+=(($F[1]-1)-$q)+length($F[4]); 
$q=($F[1]-1)+length($F[3]); 
END{
  print join "\t",$c,$r,$r+($c{$c}-$q),"+",$c,$q,$c{$c},"+",1,0,0,$c{$c}-$q;
  print STDERR join "\t",$c,$q,$c{$c},"+",$c,$r,$r+($c{$c}-$q),"+",1,0,0,$c{$c}-$q;
}
' "$@"
