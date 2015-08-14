#!/usr/bin/perl
use strict;
use warnings;

my $sam = $ARGV[0];

if ($sam eq "") {

	print "Error, not sam input was entered. Please enter a SAM file for analysis\n";
	die "Missing arguments\n";

} #end if

open(FILTER, ">$sam.filter.sam");
open(FASTQ, ">$sam.unMapped.fasta");

open(SAM, $sam);

while (my $samLine = <SAM>) {

	chomp($samLine);

	my @alignment = split("\t", $samLine);

	if (defined($alignment[9]) ) {

		if ($alignment[2] eq "\*") {

			print FASTQ ">", $alignment[0],"/1\n";
			print FASTQ $alignment[9],"\n";

			print FILTER $samLine, "\n";

		} #end if non mapped

	} #end if defined

} #end while

exit;

