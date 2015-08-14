#!/usr/bin/perl

open(INPUT, $ARGV[0]);

$headerLine = <INPUT>;

chop($headerLine);

my @headers = split("\t", $headerLine);

for($i=0; $i < 24; ++$i) {
     	open($headers[$i], ">$headers[$i].txt");
} #end for

while ($line = <INPUT>) {

	chop($line);

	@SPPS = split("\t", $line);

        for($i=0; $i < 24; ++$i) {

	        $SPP = @SPPS[$i];
        	@SPP_2 = split("-", $SPP);

	        @SPP_range = split("_", @SPP_2[1]);

        	print {$headers[$i]} ($SPP_2[0], "\t", $SPP_range[0], "\t", $SPP_range[1], "\n");

      	} #end for

} #end while

close(INPUT);

for($i=0; $i < 24; ++$i) {
     	close($header[$i]);
} #end for

exit;