#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[1])) {
	print "missing arguments\n";
	print "Loc file and Map file required\n";
	
	die;
} #end if no arguments

my %markers;

my $loc = $ARGV[0];

my $map = $ARGV[1];


open(LOC, $loc) or die "Can't open file $loc\n";
open(MAP, $map) or die "Can't open file $map\n";
open(ORDERED, ">$map.loc");

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

