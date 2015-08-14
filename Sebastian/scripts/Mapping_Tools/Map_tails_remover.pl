#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[1])) {
	print "missing arguments\n";
	print "map file and tail threshold required\n";
	
	die;
} #end if no arguments


my $map = $ARGV[0];
my $threshold = $ARGV[1];
open(MAP, $map) or die "Can't open file $map\n";

my @mapOrder;

while (my $positions = <MAP>) {
	
	chomp($positions);

	push(@mapOrder, $positions);

} #end reading map

close(MAP);

# Check the beggining
while(@mapOrder) {

	my $first = shift (@mapOrder);

	my @firstMarker = split("\t", $first);
	my @secondMarker = split("\t", $mapOrder[0]);

	if( ($secondMarker[2] - $firstMarker[2]) < $threshold) {
		unshift(@mapOrder, $first);

		last;
	} #end if

} #end while

# Check the end
while(@mapOrder) {

	my $last = pop (@mapOrder);

	my $nMarker = scalar(@mapOrder);

	my @lastMarker = split("\t", $last);
	my @almostLastMarker = split("\t", $mapOrder[$nMarker-1]);

	if( ($lastMarker[2] - $almostLastMarker[2]) < $threshold) {
		push(@mapOrder, $last);

		last;
	} #end if

} #end while

# Print de-tailed map and recalculate distances (if first marker still 0, there's no difference)

open(DETAILED, ">$map.withoutTails.map");

my $first = shift (@mapOrder);
my @firstMarker = split("\t", $first);

print DETAILED $firstMarker[0], "\t", $firstMarker[1], "\t", "0", "\t", $firstMarker[3], "\n";

foreach my $marker (@mapOrder) {

	my @markerData = split("\t", $marker);

	print DETAILED $markerData[0], "\t", $markerData[1], "\t", $markerData[2]-$firstMarker[2], "\t", $markerData[3], "\n";

} #end printing map

close(DETAILED);

exit;

