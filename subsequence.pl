#!/usr/bin/perl
# script.pl by kelsey
use strict; use warnings;
use FAlite;

#prints substrings of fasta sequences

die "usage: fastasub.pl <FASTA file> <start> [end (optional)]" unless @ARGV >= 2;
my $fasta = $ARGV[0];
my $start = $ARGV[1];
my $end = 0;
if (defined $ARGV[2]){
    $end = $ARGV[2];
}

open(my $in, "$fasta") or die "error reading $fasta for reading";

while (<$in>){
    chomp;
    if (m/>/){                      # if the line has a > in it then its a description
        my $des = $_;
        print $des,"\n";
    }
    else{
        my $seq = $_;
        if ($end == 0){
            $end = length($seq);
        }
        print (substr($seq,($start-1),$end),"\n");
    }
}

close $in;