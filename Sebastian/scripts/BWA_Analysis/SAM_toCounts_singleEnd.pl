#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[0])) {

	print "Error, not sam input was entered. Please enter a SAM file for analysis\n";
	die "Missing arguments\n";

} #end if

my $sam = $ARGV[0];

open(SAM, $sam) or die "Can't open $sam\n";
open(COUNTS, ">$sam.counts.single.txt");

my %references;
my $mapReads;
my $totalReads;

while (my $samLine = <SAM>) {

	if($samLine =~ /^@/) {next}

	chomp($samLine);

	my @alignment = split("\t", $samLine);

	if ($alignment[2] ne "\*" && defined($alignment[9])) {

		++$references{$alignment[2]};
		++$mapReads

	} #end of map/unmap

	++$totalReads;

} #end while
close(SAM);

print "Number of map reads:	", $mapReads, "\n";
print "Number of unmap reads:	", $totalReads-$mapReads, "\n";

foreach my $reference (sort {$a cmp $b} keys %references) {

	print COUNTS $reference , "\t", $references{$reference}, "\n";

} #end counting references

exit;













