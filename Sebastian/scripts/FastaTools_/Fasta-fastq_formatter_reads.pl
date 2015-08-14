#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[0])) { die "Please provide a fasta file (one-liner sequence) to convert\n"}

#open DP File
open(INPUT, $ARGV[0]) || die "Cannot open file \"$ARGV[0]\""; # Can't open fasta file to convert

my $file = $ARGV[0];

$file =~ s{\.[^.]+$}{};

open(RESULTS, ">$file.fake.fastq") || die "Cannot open file results file";

while (my $line1 = <INPUT>) {

	if(substr($line1,0,1) eq ">") {

		chomp ($line1);

		my $headerline = substr($line1,1,100);

		my $line2 = <INPUT>;

		chomp ($line2);

		my $length = length($line2);

		my @quality = ();

	   	for(my $i = 0; $i < $length; ++$i) {
			push (@quality, "I");
		} #end for

		if ($length >= 25) {

			print RESULTS "@", $headerline , "\n";
		        print RESULTS $line2, "\n";
		       	print RESULTS "+\n";
			print RESULTS @quality, "\n";

		} #end if

	} #end if line check

} #end while

close(INPUT);
close(RESULTS);

exit;
