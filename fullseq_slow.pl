#!/usr/bin/perl
# script.pl by kelsey
use strict; use warnings;
use FAlite;
use List::MoreUtils qw(uniq);

#gets full sequence from a multiple sequence alignment of an HMMsearch
#if msa is in stockholm please convert to fasta first

die "usage: fullseq.pl <multiple seq alignment> <FASTA file with full seqs>" unless @ARGV == 2;
my $msa = $ARGV[0];
my $genome = $ARGV[1];

#get names from MSA
my @names;

open(my $in, "$msa") or die "error reading $msa for reading";

while (<$in>){
    chomp;
    if ($_ =~ /\# STOCKHOLM/){die "MSA must be in FASTA format"};
    if (m/>/){                      #if the line has a > in it then its a description
        my $des = $_;
        my ($name) = $des =~ m/(\w+)/;  #extract the gene name from description
        push (@names, $name);
    }
}

close $in;

@names = uniq(@names);

open($in, "$genome") or die "error reading $genome for reading";
my $fasta = new FAlite($in);

# typical reading loop
while (my $entry = $fasta->nextEntry) {
    my $def = $entry->def;            # get definition line from object
	my $seq = $entry->seq;            # get sequence from object
	
    foreach my $name (@names){
        if ($def =~ /$name/){
            print "$entry"
        }
    }
}