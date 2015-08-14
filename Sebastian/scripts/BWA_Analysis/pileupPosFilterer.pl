#!/usr/bin/perl

open(INPUT, $ARGV[0]);

open(RESULTS, ">$ARGV[0].depth$ARGV[1].filtered.txt");

while ($line = <INPUT>) {

	chop($line);

        @pos = split("\t", $line);

        if ($pos[4] > 1 && $pos[7] > $ARGV[1]) {

        	print RESULTS $line, "\n";

        } #end if

} #end while

close(INPUT);

close(RESULTS);
