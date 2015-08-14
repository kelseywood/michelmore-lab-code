#!/usr/bin/perl
# randoseq.pl by kelsey
use strict; use warnings;

#gets random sequences from a FASTA file
#input = fasta and number of sequences you want to pull out

die "usage: randoseq.pl <FASTA file> <number of sequences>" unless @ARGV == 2;
my $fasta = $ARGV[0];
my $num = $ARGV[1];

#count number of sequences in FASTA

my $total = 0;
my $count = 0;

open(IN, "$fasta") or die "error reading $fasta for reading";

while (<IN>){
    chomp;
    if (m/>/){
        $total++;
    }
}

close IN;

open(IN, "$fasta") or die "error reading $fasta for reading";

#generate array of random numbers
my @rando = {};
for (my $i = 0; $i < $num; $i++){
    push @rando, int(rand($total));
}

while (<IN>){
    chomp;
    if (m/>/){
        $count++;
        my $line =$_;
        if ($count ~~ @rando){
            print "$line\n";
            my $nextline = <IN>; #get second row of sequence
            print "$nextline";
        }
    }
}

close IN;