#!/bin/env perl
use strict;
my %g; # gene_id => \%ref_gene_ids (or gene_names)
my @prep; # array of [line, original_id]
while (<>) {
 my @t=split(/\t/);
 unless (@t>8) { print $_; next }
 my ($gid)=($t[8]=~m/gene_id "(MSTRG\.\d+)"/);
 my ($rn)=($t[8]=~m/ref_gene_id "([^"]+)/);#or gene_name
 unless ($gid && $rn) { print $_; next }
 push(@prep, [$_, $gid]);
 my $h=$g{$gid};
 if ($h) { $h->{$rn}=1 }
 else { $g{$gid}= { $rn=>1 } }
}
foreach my $d (@prep) {
 my ($line, $gid)=@$d;
 my @rids=keys(%{$g{$gid}});
 my $nid=$gid.'|'.join('|',@rids); #new gene_id
 $line=~s/gene_id "MSTRG\.\d+"/gene_id "$nid"/;
 print $line;
}