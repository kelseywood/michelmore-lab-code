#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[2])) {
	print "missing arguments\n";
	print "Loc file, LG file and prefix required\n";
	
	die;
} #end if no arguments

my %markers;

open(LOC, $ARGV[0]);

my $header = "true";
my @individualNames;

while(my $line = <LOC>) {

	if($header eq "true") {

		if($line =~ /^;/) {
			$header = "false";
                        @individualNames = split("\t", $line);
                        shift(@individualNames);
		} #end if

	} else {

		my @data = split("\t", $line);

		$markers{$data[0]} = $line;

	} #end if header

} #end reading genotypes

close(LOC);

open(LGS, $ARGV[1]);

my %LGs;

while (my $positions = <LGS>) {

	my @info = split("\t", $positions);

	push(@{$LGs{$info[0]}}, $info[1]);

} #end reading map

close(LGS);

foreach my $LG (keys %LGs) {

	open(RESULTS, ">$ARGV[2].LG$LG.loc");

	print RESULTS ";";

	foreach my $individual( @individualNames ) {

		print RESULTS "\t", $individual;

	} #end foreach

	foreach my $marker (@{$LGs{$LG}}) {

		if(defined($markers{$marker})) {

			print RESULTS $markers{$marker};

		} #end if haplotype defined

	} #end foreach marker

	close(RESULTS);

} #en foreach LG

exit;

