#!/usr/bin/perl

open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0]_shortName.txt");

while($line = <FASTAINPUT>) {

	chomp($line);

        @header = split(" ", $line);

        $sequence = <FASTAINPUT>;

        chomp($sequence);

        print OUTPUT $header[0], "\n";
        print OUTPUT $sequence, "\n";

} #end while

close(FASTAINPUT);
close(OUTPUT);

exit;
