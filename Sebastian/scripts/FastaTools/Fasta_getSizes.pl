#!/usr/bin/perl

use strict;
use warnings;

open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0]_sequenceSizes.txt");

my $sequenceName;
my $sequenceSize = 0;

while(my $line = <FASTAINPUT>) {

	chomp($line);
	$line =~ s/\r//;

	if( $line =~ /^>/) {

		if (!defined($sequenceName)) {
			my @header = split(" ", $line);
		        $sequenceName = $header[0];
			$sequenceName =~ s/>//;
			print OUTPUT $sequenceName, "\t";
		} else {

			print OUTPUT $sequenceSize, "\n";
			my @header = split(" ", $line);
		        $sequenceName = $header[0];
			$sequenceName =~ s/>//;
			print OUTPUT $sequenceName, "\t";
			$sequenceSize = 0;

		}

	} else {

		$sequenceSize += length($line);

	} #end else

} #end while

print OUTPUT $sequenceSize;

close(FASTAINPUT);
close(OUTPUT);

exit;

