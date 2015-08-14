#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[1])) {
	die "missing arguments\n";
} #end if no arguments


open(LOC, $ARGV[0]) or die "Can't open $ARGV[0]";

my $header = "true";
my @markers;
my %individuals;
my @individualNames;

my %indInfo;
my %markerInfo;

while(my $line = <LOC>) {

	if($header eq "true") {

 		if($line =~ /^;/) {
        		$header = "false";

                        $line =~ s/\;//;

                        @individualNames = split("\t", $line);
                } #end if it's the individual names

	} else {

        	chomp($line);

		my @data = split("\t", $line);

                my $markername = shift(@data);

                push(@markers, $markername);

                my $i = 0;

                foreach my $datapoint (@data) {

                	push(@{$individuals{$individualNames[$i]}}, $datapoint);
                        ++$i;

			if($datapoint eq "-") {++$indInfo{$individualNames[$i]}{'missing'};++$markerInfo{$markername}{'missing'}}


                } #end loading hash

	} #end if header

} #end reading genotypes

close(LOC);

foreach my $individual( keys %individuals ) {

	my @haplotype = @{$individuals{$individual}};

	my @transitions;

	if($haplotype[0] eq $haplotype[1] || $haplotype[0] eq "-" || $haplotype[1] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[1] eq $haplotype[2] || $haplotype[1] eq "-" || $haplotype[2] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[2] eq $haplotype[3] || $haplotype[2] eq "-" || $haplotype[3] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[3] eq $haplotype[4] || $haplotype[3] eq "-" || $haplotype[4] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }
	if($haplotype[4] eq $haplotype[5] || $haplotype[4] eq "-" || $haplotype[5] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }

	for (my $i = 1; $i < (scalar(@haplotype)-6); ++$i) {

        	shift(@transitions);

                if($haplotype[$i + 4] eq $haplotype[$i + 5] || $haplotype[$i + 4] eq "-" || $haplotype[$i + 5] eq "-" ) { push(@transitions, 0) } else { push(@transitions, 1) }

#                print @transitions, "\n";

                if( (eval join("+", @transitions)) > 1 ) {

			if ( $transitions[0] == 0 && $transitions[1] == 1 && $transitions[2] == 1 && $transitions[3] == 0 && $transitions[4] == 0 ) {

                             ++$indInfo{$individual}{'DR'};
                             ++$markerInfo{$markers[$i + 2]}{'DR'};
                             $individuals{$individual}[$i + 2] = "-";
                             $transitions[1] = 0;
                             $transitions[2] = 0;

                        } elsif ( $transitions[0] == 0 && $transitions[1] == 0 && $transitions[2] == 1 && $transitions[3] == 1 && $transitions[4] == 0) {

                             ++$indInfo{$individual}{'DR'};
                             ++$markerInfo{$markers[$i + 3]}{'DR'};

                             $individuals{$individual}[$i + 3] = "-";
                             $transitions[2] = 0;
                             $transitions[3] = 0;

                        } elsif ( $transitions[0] == 0 && $transitions[1] == 1 && $transitions[2] == 1 && $transitions[3] == 1 && $transitions[4] == 0)  {

                             ++$indInfo{$individual}{'DR'};
                             ++$markerInfo{$markers[$i + 2]}{'DR'};
                             ++$markerInfo{$markers[$i + 3]}{'DR'};

                             $individuals{$individual}[$i + 2] = "-";
                             $individuals{$individual}[$i + 3] = "-";
                             $transitions[1] = 0;
                             $transitions[2] = 0;
                             $transitions[3] = 0;

                        } elsif ( $transitions[0] == 0 && $transitions[1] == 1 && $transitions[2] == 1 && $transitions[3] == 1 && $transitions[4] == 1)  {

                             ++$indInfo{$individual}{'DR'};
                             ++$markerInfo{$markers[$i + 2]}{'DR'};

                             $individuals{$individual}[$i + 2] = "-";
                             $transitions[1] = 0;

                        } #end checking conditions

                } #there is a double recombination

        } #end for haplotype

} #end foreach individuals

my $sumDRInd = 0;
my $sumMissingInd = 0;

my $maxDRInd = 0;
my $maxMissingInd = 0;

open(INDIVIDUALCOUNTS, ">$ARGV[0].countsPerIndividual.loc");

foreach my $individual( keys %indInfo ) {

	if(!defined($indInfo{$individual}{'missing'})) {$indInfo{$individual}{'missing'} = 0}

	if(!defined($indInfo{$individual}{'DR'})) {$indInfo{$individual}{'DR'} = 0}

	print INDIVIDUALCOUNTS $individual, "\t", $indInfo{$individual}{'missing'}, "\t", $indInfo{$individual}{'DR'}, "\n";

	$sumDRInd += $indInfo{$individual}{'DR'};
	$sumMissingInd += $indInfo{$individual}{'missing'};

	if($indInfo{$individual}{'DR'} > $maxDRInd) { $maxDRInd = $indInfo{$individual}{'DR'}}
	if($indInfo{$individual}{'missing'} > $maxMissingInd) { $maxMissingInd = $indInfo{$individual}{'missing'}}

} #end foreach

my $averageDRInd = $sumDRInd/scalar(keys %indInfo);
my $averageMissingInd = $sumMissingInd/scalar(keys %indInfo);

my $sumDRMarker = 0;
my $sumMisingMarker = 0;

my $maxDRMarker = 0;
my $maxMissingMarker = 0;

open(MARKERCOUNTS, ">$ARGV[0].countsPerMarker.loc");

foreach my $marker( keys %markerInfo ) {

	if(!defined($markerInfo{$marker}{'missing'})) {$markerInfo{$marker}{'missing'} = 0}

	if(!defined($markerInfo{$marker}{'DR'})) {$markerInfo{$marker}{'DR'} = 0}

	print MARKERCOUNTS $marker, "\t", $markerInfo{$marker}{'missing'}, "\t", $markerInfo{$marker}{'DR'}, "\n";

	$sumDRMarker += $markerInfo{$marker}{'DR'};
	$sumMisingMarker += $markerInfo{$marker}{'missing'};

	if($markerInfo{$marker}{'DR'} > $maxDRMarker) { $maxDRMarker = $markerInfo{$marker}{'DR'}}
	if($markerInfo{$marker}{'missing'} > $maxMissingMarker) { $maxMissingMarker = $markerInfo{$marker}{'missing'}}

} #end foreach

my $averageDRMarker = $sumDRMarker/scalar(keys %markerInfo);
my $averageMissingMarker = $sumMisingMarker/scalar(keys %markerInfo);


open(MAP, $ARGV[1]) or die "Can't open $ARGV[1]";

my %positions;

my $prevPos = 0;
my $maxGap = 0;

my $numMarkers;

while(my $mapLine = <MAP>) {

	chomp($mapLine);

	my @mapInfo = split("\t", $mapLine);

	++$positions{$mapInfo[2]};

	if( ($mapInfo[2]-$prevPos) > $maxGap) {$maxGap = $mapInfo[2]-$prevPos}

	++$numMarkers;

	$prevPos = $mapInfo[2];

} #end while

my @markersPerBin =  sort {$positions{$b} <=> $positions{$a}} keys %positions;

my @geneticPositions =  sort {$b <=> $a} keys %positions;

my $numGeneticPositions = scalar(@geneticPositions);

my $size = shift(@geneticPositions);

my $maxMarkersPerBin = $positions{shift(@markersPerBin)};

open(STATS, ">$ARGV[1].stats.txt");

print STATS "Total number of markers = $numMarkers\n";
print STATS "Total size in cM = $size\n";
print STATS "Total number of genetic locations = $numGeneticPositions\n";
print STATS "Max bin size = $maxMarkersPerBin\n";
print STATS "Max gap size = $maxGap\n\n";
print STATS "Average number of double recombinants per individual = ", sprintf("%.3f", $averageDRInd), " \n";
print STATS "Average number of double recombinants per marker = ", sprintf("%.3f", $averageDRMarker), " \n";
print STATS "Maximun number of double recombinants per individual = $maxDRInd \n";
print STATS "Maximun number of double recombinants per marker = $maxDRMarker \n\n";
print STATS "Average number of missings per individual = ", sprintf("%.3f", $averageMissingInd), " \n";
print STATS "Average number of missings per marker = ", sprintf("%.3f", $averageMissingMarker), " \n";
print STATS "Maximun number of missings per individual = $maxMissingInd \n";
print STATS "Maximun number of missings per marker = $maxMissingMarker \n";

close(STATS);

exit;