#!/usr/bin/perl

use strict; use warnings;

open(FASTQ, $ARGV[0]);

my $fasta = $ARGV[0];

$fasta =~ s{\.[^.]+$}{};

open(FASTA, ">$fasta.FASTA");

my $stat = "";

while(my $header = <FASTQ>) {

	my $sequence = <FASTQ>;
	my $header2 = <FASTQ>;
	my $qual = <FASTQ>;

	$header =~ s/@/>/;

	print FASTA $header;
	print FASTA $sequence;

} #end while






