#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[0])) {
	print "missing arguments\n";
	die;
} #end if no arguments


open(LOC, $ARGV[0]);

my $header = "true";
my @markers;
my %individuals;
my @individualNames;

while(my $line = <LOC>) {

	if($header eq "true") {

 		if($line =~ /^;/) {
        		$header = "false";

                        @individualNames = split("\t", $line);
                        shift(@individualNames);
                } #end if it's the individual names

	} else {

        	chomp($line);

                $line =~ s/(H|h)/-/g;

		my @data = split("\t", $line);

                my $markername = shift(@data);

                push(@markers, $markername);

                my $i = 0;

                foreach my $datapoint (@data) {

                	push(@{$individuals{$individualNames[$i]}}, $datapoint);
                        ++$i;

                } #end loading hash

	} #end if header

} #end reading genotypes

close(LOC);

my %doubleRecombinantsInd;
my %doubleRecombinantsMarkers;

foreach my $individual( keys %individuals ) {

	my @haplotype = @{$individuals{$individual}};

	my @transitions;

	if($haplotype[0] eq $haplotype[1] || $haplotype[0] eq "-" || $haplotype[1] eq "-" || $haplotype[0] eq "H" || $haplotype[1] eq "H" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[1] eq $haplotype[2] || $haplotype[1] eq "-" || $haplotype[2] eq "-" || $haplotype[1] eq "H" || $haplotype[2] eq "H" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[2] eq $haplotype[3] || $haplotype[2] eq "-" || $haplotype[3] eq "-" || $haplotype[2] eq "H" || $haplotype[3] eq "H" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[3] eq $haplotype[4] || $haplotype[3] eq "-" || $haplotype[4] eq "-" || $haplotype[3] eq "H" || $haplotype[4] eq "H" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[4] eq $haplotype[5] || $haplotype[4] eq "-" || $haplotype[5] eq "-" || $haplotype[4] eq "H" || $haplotype[5] eq "H" ) { push(@transitions, 0) } else { push(@transitions, 1) }

	for (my $i = 1; $i < (scalar(@haplotype)-6); ++$i) {

        	shift(@transitions);

                if($haplotype[$i + 4] eq $haplotype[$i + 5] || $haplotype[$i + 4] eq "-" || $haplotype[$i + 5] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }

#                print @transitions, "\n";

                if( (eval join("+", @transitions)) > 1 ) {

			if ( $transitions[0] == 0 && $transitions[1] == 1 && $transitions[2] == 1 && $transitions[3] == 0 && $transitions[4] == 0 ) {

                           ++$doubleRecombinantsInd{$individual};
                             ++$doubleRecombinantsMarkers{$markers[$i + 2]};
                             $individuals{$individual}[$i + 2] = "-";
                             $transitions[1] = 0;
                             $transitions[2] = 0;

                        } elsif ( $transitions[0] == 0 && $transitions[1] == 0 && $transitions[2] == 1 && $transitions[3] == 1 && $transitions[4] == 0) {

                             ++$doubleRecombinantsInd{$individual};
                             ++$doubleRecombinantsMarkers{$markers[$i + 3]};

                             $individuals{$individual}[$i + 3] = "-";
                             $transitions[2] = 0;
                             $transitions[3] = 0;

                        } elsif ( $transitions[0] == 0 && $transitions[1] == 1 && $transitions[2] == 1 && $transitions[3] == 1 && $transitions[4] == 0)  {

                             ++$doubleRecombinantsInd{$individual};
                             ++$doubleRecombinantsMarkers{$markers[$i + 2]};
                             ++$doubleRecombinantsMarkers{$markers[$i + 3]};

                             $individuals{$individual}[$i + 2] = "-";
                             $individuals{$individual}[$i + 3] = "-";
                             $transitions[1] = 0;
                             $transitions[2] = 0;
                             $transitions[3] = 0;

                        } elsif ( $transitions[0] == 0 && $transitions[1] == 1 && $transitions[2] == 1 && $transitions[3] == 1 && $transitions[4] == 1)  {

                             ++$doubleRecombinantsInd{$individual};
                             ++$doubleRecombinantsMarkers{$markers[$i + 2]};

                             $individuals{$individual}[$i + 2] = "-";
                             $transitions[1] = 0;

                        } #end checking conditions

                } #there is a double recombination

        } #end for haplotype

} #end foreach individuals

open(INDIVIDUALCOUNTS, ">$ARGV[0].countsPerIndividual.loc");

foreach my $individual( keys %doubleRecombinantsInd ) {

	print INDIVIDUALCOUNTS $individual, "\t", $doubleRecombinantsInd{$individual}, "\n";

} #end foreach

open(MARKERCOUNTS, ">$ARGV[0].countsPerMarker.loc");

foreach my $marker( keys %doubleRecombinantsMarkers ) {

	print MARKERCOUNTS $marker, "\t", $doubleRecombinantsMarkers{$marker}, "\n";

} #end foreach


exit;