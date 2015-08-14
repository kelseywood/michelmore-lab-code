#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[0])) {
	print "missing arguments, provide list of ranges\n";
	die;
} #end if no arguments


open(RANGES, $ARGV[0]);

open(POSITIONS, ">$ARGV[0].positions.txt");

while(my $line = <RANGES>) {

	chomp($line);

        my @rangeInfo = split("\t", $line);

        for(my $i = $rangeInfo[1]; $i <= $rangeInfo[2]; ++$i) {

        	print POSITIONS $rangeInfo[0]. "\t", $i, "\n";

	} # end printing positions

} #end reading ranges

close(RANGES);

exit;