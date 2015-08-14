#!/usr/bin/perl
use strict;
use warnings;

if ($ARGV[0] eq "") {

	print "Error, not sam input was entered. Please enter a file with the SAM files to analyzed\n";
	die "Missing arguments\n";

} #end if

open(FILES, $ARGV[0]);

while (my $sam = <FILES>) {

	chomp($sam);

	print "Analyzing ",$sam,"\n";

	open(FILTER, ">$sam.filter.sam");
	open(FASTQ1, ">$sam.1.fastq");
	open(FASTQ2, ">$sam.2.fastq");

	my @forwardQuals = ("73","83","97","113","117");

	my @reverseQuals = ("133","145","153","163","177");

	open(SAM, $sam);

	while (my $samLine = <SAM>) {

		chomp($samLine);

		my @alignment = split("\t", $samLine);

		foreach my $qual1 (@forwardQuals) {

			if ($qual1 == $alignment[1]) {

				print FASTQ1 "@", $alignment[0],"/1\n";
				print FASTQ1 $alignment[9],"\n";
				print FASTQ1 "+\n";
				print FASTQ1 $alignment[10],"\n";

				print FILTER $samLine, "\n";

			} #end if

		} #end foreach

		foreach my $qual2 (@reverseQuals) {

			if ($qual2 == $alignment[1]) {

				print FASTQ2 "@", $alignment[0],"/2\n";
				print FASTQ2 $alignment[9],"\n";
				print FASTQ2 "+\n";
				print FASTQ2 $alignment[10],"\n";

				print FILTER $samLine, "\n";

			} #end if

		} #end foreach

	} #end while

	print "Compressing files for ",$sam,"\n";

	system("gzip $sam.1.fastq");
	system("gzip $sam.2.fastq");
	system("gzip $sam.filter.sam");

} #end while files

exit;

