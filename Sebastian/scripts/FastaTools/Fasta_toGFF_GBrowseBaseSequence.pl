#!/usr/bin/perl

use strict;
use warnings;

if(!defined($ARGV[0])) {die "Missing argument, please provided a fasta file"}

open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0].gff3");

my $sequenceName;
my $sequenceSize = 0;

while(my $line = <FASTAINPUT>) {

	chomp($line);
	$line =~ s/\r//;

	if( $line =~ /^>/) {

		if (!defined($sequenceName)) {
			my @fastaHeader = split(" ", $line);
		        $sequenceName = $fastaHeader[0];
			$sequenceName =~ s/>//;
			print OUTPUT $sequenceName, "\t";
		} else {

			print OUTPUT ".\tdna\t1\t", $sequenceSize, "\t.\t.\t.\tID=",$sequenceName,";Name=",$sequenceName,";\n";
			my @fastaHeader = split(" ", $line);
		        $sequenceName = $fastaHeader[0];
			$sequenceName =~ s/>//;
			print OUTPUT $sequenceName, "\t";
			$sequenceSize = 0;

		}

	} else {

		$sequenceSize += length($line);

	} #end else

} #end while

print OUTPUT ".\tdna\t1\t", $sequenceSize, "\t.\t.\t.\tID=",$sequenceName,";Name=",$sequenceName,";\n";

close(FASTAINPUT);
close(OUTPUT);

exit;
