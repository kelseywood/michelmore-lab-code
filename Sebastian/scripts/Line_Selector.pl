#!/usr/bin/perl
use strict;
use warnings;

if( !defined($ARGV[2])) { die "not enough arguments\n" }


open(INPUT1,$ARGV[0]) || die "Cannot open file \"$ARGV[0]\""; # Table where to select to
open(INPUT2,$ARGV[1]) || die "Cannot open file \"$ARGV[1]\""; # Line indeces to extract
open(RESULTS,">$ARGV[2]")|| die "Cannot open the Results file"; # output file

my %forExtract = ();

while(my $selected = <INPUT2>) {

	chomp($selected);

        my @toSelect = split("\t", $selected);

	$forExtract {$toSelect[0]} = 1;

} #end while

close(INPUT2);

my $lineNum = 1;

while (my $line = <INPUT1>) {

	if(defined($forExtract{$lineNum})) {
		print RESULTS $line;
	} #end if

	++$lineNum;
	
} #end while

close(INPUT1);
close(RESULTS);

exit;

#end
