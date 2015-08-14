#!/usr/bin/perl
#######################################################
#                                                     #
#                                       #
#                                                     #
#              Sebastian Reyes-Chin-Wo                #
#               sreyesch@ucdavis.edu                  #
#                                                     #
#######################################################
use strict;
use warnings;

if( !(defined($ARGV[4])) || (defined($ARGV[5]))) {

	print "Missing Arguments or extra arguments.\n";
	print "Please input the info2 file (only one file), a spacing distance to merge regions and the source for the gff file.\n";
	die "Missing Arguments";

} #end if

my $spacing = $ARGV[1];

my $source = $ARGV[2];

my $minAlignment = $ARGV[3];

my $minIdentity = $ARGV[4];

open(INFO2, $ARGV[0]);
open(MATCHEDREGIONS, ">$ARGV[0].minAlignment$minAlignment.minIdentity$minIdentity.matched_regions.gff");

my %matches;

while(my $line = <INFO2>) {

	chomp($line);

	my @data = split("\t", $line);

        if ($data[2] ne "no_hits_found") {

	        my @lengths = split("/", $data[13]);

	        my $covQuery = $data[6]/$lengths[0];

	        my $covSubject = $data[6]/$lengths[1];

	        my @sizes = split("/", $data[13]);

	        my $start;

	        my $end;

	        if ( ($data[12] - $data[11]) < 0 ) {

	                $start = $data[12];
	                $end = $data[11];

	        } else {

	                $start = $data[11];
	                $end = $data[12];

	        } #end else

	        if ( ($minAlignment <= $data[6]) && ($minIdentity <= $data[4]) ) {

	                push( @{$matches{$data[1]}}, [$start, $end, $data[0]]);

	        } #end if for good match

           } #end if it is a hit

} #end while

foreach my $query (sort keys %matches) {

	my $startRegion = 0;

	my $endRegion = 0;

	my $hitCount = 0;

	my $gapCount = 0;

	my @hits = ();

	my @AoA = @{$matches{$query}};

	my @sorted = sort {$a->[0] <=> $b->[0]} @AoA;

	for (@sorted) {

		my @match = @$_;

		if ( ($startRegion == 0) || ($endRegion == 0) ) {

			$startRegion = $match[0];

			$endRegion = $match[1];

			++$hitCount;

			push(@hits, $match[2]);

		} elsif ( $endRegion > $match[0] ) {

			#some kind of overlap

			if ( ($startRegion <= $match[0]) || ($endRegion < $match[1]) ) {

				$endRegion = $match[1];

				++$hitCount;

				push(@hits, $match[2]);

			} elsif ( ($startRegion <= $match[0]) || ($endRegion >= $match[1]) ) {

				++$hitCount;

				push(@hits, $match[2]);

			}

		} elsif ( $endRegion > ($match[0] - $spacing) ) {

			#gap found

			++$gapCount;

			++$hitCount;

			push(@hits, $match[2]);

			$endRegion = $match[1];

		} else {

			if( ($endRegion-$startRegion) != 0 ) {

				print MATCHEDREGIONS $query, "\t", $source, "\t", "BlastMatches", "\t", $startRegion, "\t", $endRegion, "\t", ".", "\t", ".", "\t", ".", "\t";

				print MATCHEDREGIONS "ID=", $query,"-",$startRegion, "_", $endRegion,";";
				print MATCHEDREGIONS "Name=", $query,"-",$startRegion, "_", $endRegion,";";
				print MATCHEDREGIONS "NumHits=", $hitCount, ";", "Gaps=", $gapCount,";", "Hits=";

				print MATCHEDREGIONS $hits[0];

				shift(@hits);

				foreach my $hit (@hits) {

					print MATCHEDREGIONS  "/", $hit;

				} #end foreach

				print MATCHEDREGIONS ";\n";

			} #end if

			$startRegion = $match[0];

			$endRegion = $match[1];

			$hitCount = 1;

			$gapCount = 0;

			@hits = ();

			push(@hits, $match[2]);

		} #end else

	} #end for

	if( ($endRegion-$startRegion) != 0 ) {

		print MATCHEDREGIONS $query, "\t", $source, "\t", "BlastMatches", "\t", $startRegion, "\t", $endRegion, "\t", ".", "\t", ".", "\t", ".", "\t";

		print MATCHEDREGIONS "ID=", $query,"-",$startRegion, "_", $endRegion,";";
		print MATCHEDREGIONS "Name=", $query,"-",$startRegion, "_", $endRegion,";";
		print MATCHEDREGIONS "NumHits=", $hitCount, ";", "Gaps=", $gapCount,";", "Hits=";

		print MATCHEDREGIONS $hits[0];

		shift(@hits);

		foreach my $hit (@hits) {

			print MATCHEDREGIONS "/", $hit;

		} #end foreach

		print MATCHEDREGIONS ";\n";

	} #end if

}#end foreach

exit;
