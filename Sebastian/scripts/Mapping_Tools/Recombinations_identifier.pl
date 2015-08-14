#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[1])) {
	print "Requires a loci file and a map file\n";
	print "missing arguments\n";
	die;
} #end if no arguments


open(LOC, $ARGV[0]);

my $header = "true";
my %markers;
my @individualNames;

while(my $line = <LOC>) {

	if($header eq "true") {

 		if($line =~ /^;/) {

			chomp($line);

        		$header = "false";

			$header =~ s/\t/\ti/g;

                        @individualNames = split("\t", $line);
                        shift(@individualNames);
                } #end if it's the individual names

	} else {

        	chomp($line);

		my @data = split("\t", $line);

                my $markername = shift(@data);

		$markers{$markername} = \@data;

	} #end if header

} #end reading genotypes

close(LOC);

my %map;

open(MAP, $ARGV[1]);

while (my $pos = <MAP>) {

	chomp($pos);

	my @info = split("\t", $pos);

	if(!defined($info[1])){next}

	push (@{$map{$info[0]}}, $info[1]);

} #end loading the map

my %crossOvers;

foreach my $lg ( keys %map ) {

	my @positions = @{$map{$lg}};

	for (my $j = 0; $j < (scalar(@individualNames)); ++$j) {

	        my @transitions;

	        if($markers{$positions[0]}[$j] eq $markers{$positions[1]}[$j] || $markers{$positions[0]}[$j] eq "-" || $markers{$positions[1]}[$j] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	        if($markers{$positions[1]}[$j] eq $markers{$positions[2]}[$j] || $markers{$positions[1]}[$j] eq "-" || $markers{$positions[2]}[$j] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	        if($markers{$positions[2]}[$j] eq $markers{$positions[3]}[$j] || $markers{$positions[2]}[$j] eq "-" || $markers{$positions[3]}[$j] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	        if($markers{$positions[3]}[$j] eq $markers{$positions[4]}[$j] || $markers{$positions[3]}[$j] eq "-" || $markers{$positions[4]}[$j] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	        if($markers{$positions[4]}[$j] eq $markers{$positions[5]}[$j] || $markers{$positions[4]}[$j] eq "-" || $markers{$positions[5]}[$j] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }

		for (my $i = 1; $i < (scalar(@positions)-2); ++$i) {

			if($transitions[0] == 0 && $transitions[1] == 1 && $transitions[2] == 0 && $transitions[3] == 0 && $transitions[4] == 0 ) {

				++$crossOvers{$lg}{$individualNames[$j]};

			} #end if crossover found

		} #end for indviduals

        } #end for haplotype

} #end foreach individuals

open(RESULTS, ">$ARGV[1].crossovers.txt");

print RESULTS ";";
foreach my $ind (@individualNames) { print RESULTS "\t", $ind}
print RESULTS "\n";

foreach my $lgs (sort { $a <=> $b} keys %crossOvers) {

	print RESULTS $lgs;

	foreach my $individual (@individualNames) {

		if (!defined($crossOvers{$lgs}{$individual}) ) { print RESULTS "\t", 0 }

		else { print RESULTS "\t", $crossOvers{$lgs}{$individual} }

	} #end foreach individual

	print RESULTS "\n";

} #end foreach

exit;















