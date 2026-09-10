#!/usr/bin/perl
#
# Usage:  ./swap_chain_names_AB_and_order.pl  PDBfile
#

$|=1;

$molecule=$ARGV[0];

$chain1='A';
$chain2='B';

$i1=0;
$i2=0;
open(IN, "$molecule") || die "Can't open PDB file: $molecule ";
while (<IN>) {
  chop($line=$_);
  if(/^ATOM  / || /^HETATM/) {
    $ch=substr($line,21,1);
    $l=length($line);
    $left=substr($line,0,21);
    $right=substr($line,22,$l-22);
    if($ch eq $chain1) { $i2=$i2+1; $DATA2[$i2]="$left$chain2$right"; }
    if($ch eq $chain2) { $i1=$i1+1; $DATA1[$i1]="$left$chain1$right"; }
#    print "$ch\n";
  } 
  elsif(!/^END/ && !/^TER/) {
    print "$line\n";
  }
}
close(IN);

foreach $n (1 .. $i1) {
  print "$DATA1[$n]\n";
}
print "TER\n";
foreach $n (1 .. $i2) {
  print "$DATA2[$n]\n";
}
print "END\n";

exit;

