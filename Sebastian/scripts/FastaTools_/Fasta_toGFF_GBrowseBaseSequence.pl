#!/usr/bin/perl

use strict;
use warnings;

open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0].gff");

my $sequenceName;
my $sequenceSize = 0;

while(my $line = <FASTAINPUT>) {

	chomp($line);
	$line =~ s/\r//;

	if( $line =~ /^>/) {

		if (!defined($sequenceName)) {
		        $sequenceName = substr($line,1,200);
			print OUTPUT $sequenceName, "\t";
		} else {

			print OUTPUT ".\tdna\t1\t", $sequenceSize, "\t.\t.\t.\tID=",$sequenceName,";Name=",$sequenceName,";\n";
			$sequenceName = substr($line,1,200);
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
