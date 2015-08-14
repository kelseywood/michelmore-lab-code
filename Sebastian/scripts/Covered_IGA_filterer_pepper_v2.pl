#!/usr/bin/perl

$min_depth = $ARGV[1];

open(INPUT, $ARGV[0]);

open(RESULTS1, ">Sequence_covered_CM334-EJ_$min_depth.txt");
open(RESULTS2, ">Sequence_covered_CM334-Maor_$min_depth.txt");
open(RESULTS3, ">Sequence_covered_EJ-Maor_$min_depth.txt");

while ($line = <INPUT>) {

	chop ($line);

        @data = split("\t", $line);

        $CM334 = $data[10];
        $Maor = $data[16];
        $EJ = $data[4]-$data[10]-$data[16];


        if($CM334 >= $min_depth && $EJ >= $min_depth) {
        	print RESULTS1 $line;
		print RESULTS1 "\n";
        } #end if

        if($CM334 >= $min_depth && $Maor >= $min_depth) {
        	print RESULTS2 $line;
		print RESULTS2 "\n";
        } #end if

        if($EJ >= $min_depth && $Maor >= $min_depth) {
        	print RESULTS3 $line;
		print RESULTS3 "\n";
        } #end if

} #end while

close(INPUT);
close(RESULTS);

exit;