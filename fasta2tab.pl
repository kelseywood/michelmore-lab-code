#!/usr/bin/perl
# fasta2tab.pl by kelsey
use strict; use warnings;
use FAlite;

die "usage: fasta2tab.pl <FASTA file>" unless @ARGV == 1;

open(my $genome, "<$ARGV[0]") or die "error reading $ARGV[0] for reading";
my $fasta = new FAlite($genome);
while (my $entry = $fasta->nextEntry) {
    my $def = $entry->def;            # get definition line from object
	my $seq = $entry->seq;            # get sequence from object
	print "$def\t$seq\n";
}