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

	my @transitions;

	for (my $i = 1; $i < (scalar(@positions)-2); ++$i) {

		for (my $j = 0; $j < (scalar(@individualNames)); ++$j) {

#			print $positions[$i], "\t", $j, "\n";

			if( ($markers{$positions[$i]}[$j] ne $markers{$positions[$i+1]}[$j]) && ($markers{$positions[$i]}[$j] ne "-") && ($markers{$positions[$i+1]}[$j] ne "-") ) {

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
















