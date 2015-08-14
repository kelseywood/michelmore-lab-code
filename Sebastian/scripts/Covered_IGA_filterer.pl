#!/usr/bin/perl

$min_depth = $ARGV[1];

open(INPUT, $ARGV[0]);

open(RESULTS, ">$ARGV[0]_filtered_$ARGV[1]");

while ($line = <INPUT>) {

	chop ($line);

        @data = split("\t", $line);

        if ($data[3] >= $min_depth && $data[9] >= $min_depth) {

        	print RESULTS $line;
		print RESULTS "\n";

        } #end if

} #end while

close(INPUT);
close(RESULTS);

exit;