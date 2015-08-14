#!/usr/bin/perl

$min_depth = $ARGV[1];

open(INPUT, $ARGV[0]);

open(RESULTS, ">$ARGV[0]_filtered_$ARGV[1]");

while ($line = <INPUT>) {

	chop ($line);

        @data = split("\t", $line);

        $covered = 0;

        if($data[11] >= $min_depth) {++$covered;}
        if($data[17] >= $min_depth) {++$covered;}
        if(($data[4]-$data[11]-$data[17]) >= $min_depth) {++$covered;}

        if ($covered == 3) {

        	print RESULTS $line;
		print RESULTS "\n";

        } #end if

} #end while

close(INPUT);
close(RESULTS);

exit;
