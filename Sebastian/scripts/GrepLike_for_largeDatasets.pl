#!/usr/bin/perl

open(INPUT1,$ARGV[0]) || die "Cannot open file \"$ARGV[0], table to select\""; # Table where to select to
open(INPUT2,$ARGV[1]) || die "Cannot open file \"$ARGV[1], list to extrax\""; # List of name to extract
open(RESULTS,">$ARGV[2]")|| die "Cannot open the Results file"; # output file

my %forExtract = ();

while($selected = <INPUT2>) {

	chomp($selected);

        @toSelect = split("\t", $selected);

	$forExtract {$toSelect[0]} = 1;

} #end while

close(INPUT2);

while ($line = <INPUT1>) {

	chomp($line);

	@data = split("\t", $line);

	if(defined($forExtract{$data[0]})) {
		print RESULTS $line, "\n";
	} #end if

} #end while

close(INPUT1);
close(RESULTS);

exit;

#end
