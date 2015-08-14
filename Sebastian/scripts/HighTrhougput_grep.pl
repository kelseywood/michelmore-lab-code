#!/usr/bin/perl
use strict; use warnings;

open(INPUT1,$ARGV[0]) || die "Cannot open list of string to search"; # Table where to select to
open(RESULTS,">$ARGV[2]")|| die "Cannot open the Results file"; # output file

my $forExtract = ();

while(my $selected = <INPUT1>) {

	chomp($selected);

	system("grep \"$selected\" $ARGV[1] > temp.txt");

	open(TEMP, "temp.txt");

	while(my $hit = <TEMP>) {

		print RESULTS $hit;

	} #end while

} #end while

close(INPUT1);
close(RESULTS);

exit;
