#!/usr/bin/perl

open(INPUT, $ARGV[0]);

while ($line = <INPUT>) {

	@data = split("\t", $line);

	$code = shift(@data);

        open (OUTPUT, ">$code.txt");

        foreach $SPP(@data) {

		@SPPData = split("-", $SPP);

        	print OUTPUT $SPPData[0], "\t";

                @SPPRange = split("_", $SPPData[1]);

		print OUTPUT $SPPRange[0], "\t", $SPPRange[1], "\n";

        } #end foreach

        close(OUTPUT);

} #end while

close (INPUT);

exit;

