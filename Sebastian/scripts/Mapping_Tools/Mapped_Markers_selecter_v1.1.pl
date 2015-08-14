#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[2])) {
	print "missing arguments\n";
	print "Loc file, Map file and output file required\n";
	
	die;
} #end if no arguments

## Diferences from version 0
# Add argument for output file name
# Add lots of comments
# security checkpoints

my %markers;

my $loc = $ARGV[0];

my $map = $ARGV[1];

my $output = $ARGV[2];


open(LOC, $loc) or die "Can't open file $loc\n";
open(MAP, $map) or die "Can't open file $map\n";
open(ORDERED, ">$output") or die "Can't open file $output\n";

# Reading in loci information

my $header = "true";
while(my $line = <LOC>) {

	if($header eq "true") {

		print ORDERED $line;

		if($line =~ /^;/) {$header = "false"}

	} else {

		my @data = split("\t", $line);

		$markers{$data[0]} = $line;

	} #end if header

} #end reading genotypes

close(LOC);

# Reading map information and printing new loc file

while (my $positions = <MAP>) {

	my @info = split("\t", $positions);

	if(defined($markers{$info[1]})) {

		print ORDERED $markers{$info[1]};

	} else {

		print "Genotype for marker $info[1] not found\n";

	} #end genotype defined

} #end reading map

close(MAP);
close(ORDERED);

exit;

