#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[0])) {

	print "Error, not sam input was entered. Please enter a SAM file for analysis\n";
	die "Missing arguments\n";

} #end if

my $sam = $ARGV[0];

open(COUNTSPAIRED, ">$sam.counts.paired.txt");
open(COUNTSSINGLE, ">$sam.counts.nonPaired.txt");
open(EQUALLY, ">$sam.multipleEquallyScoredMappings.txt");

print COUNTSPAIRED "#Reference\tRawPairedReadCounts\n";
print COUNTSSINGLE "#Reference\tRawSingleReadCounts\n";
print EQUALLY "#ReadName\tReference1\tReference2\n";

my %reads;
my $unmapReads;

open(SAM, $sam);
while (my $samLine = <SAM>) {

	if($samLine =~ /^@/) {next}

	chomp($samLine);

	my @alignment = split("\t", $samLine);

	if ($alignment[2] ne "\*" && defined($alignment[9])) {

		++$reads{$alignment[0]}{'reference'}{$alignment[2]};
		$reads{$alignment[0]}{'score'}{$alignment[2]} += $alignment[4];

	} else {

		++$unmapReads;

	} #end of map/unmap

} #end while
close(SAM);

print "Number of map pairs of reads:	", scalar(keys %reads), "\n";
print "Number of unmap pairs of reads:	", $unmapReads/2, "\n";


my %references;

my %singleEndMappedReferences;

foreach my $read (keys %reads) {

	if ( (scalar (keys %{$reads{$read}{'reference'}})) == 1) {

		my @ref = keys %{$reads{$read}{'reference'}};

		if( $reads{$read}{'reference'}{$ref[0]} > 1) { ++$references{$ref[0]} }
		else { ++$singleEndMappedReferences{$ref[0]} } #split single mappings

	} elsif ( (scalar (keys %{$reads{$read}{'reference'}})) == 2) {

		my @ref = keys %{$reads{$read}{'reference'}};

		if ( $reads{$read}{'reference'}{$ref[0]} > $reads{$read}{'reference'}{$ref[1]} ) {

			++$references{$ref[0]};

		} elsif ( $reads{$read}{'reference'}{$ref[0]} < $reads{$read}{'reference'}{$ref[1]} ) {

			++$references{$ref[1]};

		} elsif ( $reads{$read}{'reference'}{$ref[0]} == $reads{$read}{'reference'}{$ref[1]} ) {


			if ( $reads{$read}{'score'}{$ref[0]} > $reads{$read}{'score'}{$ref[1]} ) {

				if( $reads{$read}{'reference'}{$ref[0]} > 1) { ++$references{$ref[0]} }
				else { ++$singleEndMappedReferences{$ref[0]} } #split single mappings

			} elsif ( $reads{$read}{'score'}{$ref[0]} < $reads{$read}{'score'}{$ref[1]} ) {

				if( $reads{$read}{'reference'}{$ref[1]} > 1) { ++$references{$ref[1]} }
				else { ++$singleEndMappedReferences{$ref[1]} } #split single mappings

			} elsif ( $reads{$read}{'score'}{$ref[0]} == $reads{$read}{'score'}{$ref[1]} ) {

				print EQUALLY $read, "\t", $ref[0], "\t", $ref[1], "\n";

				if( $reads{$read}{'reference'}{$ref[0]} > 1) {
					$references{$ref[0]} += 0.5;
					$references{$ref[1]} += 0.5;
				} else {
					$singleEndMappedReferences{$ref[0]} += 0.5;
					$singleEndMappedReferences{$ref[1]} += 0.5;
				} #split single mappings

			} #end elsif


		} #end elsif

	} #end multiple references

} #end counting references

foreach my $reference (sort {$a cmp $b} keys %references) {

	print COUNTSPAIRED $reference , "\t", $references{$reference}, "\n";

} #end counting references

foreach my $reference (sort {$a cmp $b} keys %singleEndMappedReferences) {

	print COUNTSSINGLE $reference , "\t", $singleEndMappedReferences{$reference}, "\n";

} #end counting references

exit;













