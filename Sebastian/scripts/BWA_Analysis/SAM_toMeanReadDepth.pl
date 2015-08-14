#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[0])) {

	print "Error, not sam input was entered. Please enter a SAM file for analysis\n";
	die "Missing arguments\n";

} #end if

my $sam = $ARGV[0];
my $phred = $ARGV[1];
my $sequenceSizes = $ARGV[2];
my $prefix = $ARGV[3];

open(COUNTS, ">$prefix.meanDepth.txt");

print COUNTS "#Reference\tRawReadCounts\tTotalReadLength\tMeanDepth\n";

my %reads;
my $unmapReads = 0;

open(SAM, $sam);
while (my $samLine = <SAM>) {

	if($samLine =~ /^@/) {next}

	chomp($samLine);

	my @alignment = split("\t", $samLine);

	if ($alignment[2] ne "\*" && defined($alignment[9]) && $alignment[4] >= $phred) {

		if(defined($reads{$alignment[0]}{'flag'}{$alignment[2]})) {
			if($reads{$alignment[0]}{'flag'}{$alignment[2]} == $alignment[1]) {next}
		} #end if defined

		++$reads{$alignment[0]}{'reference'}{$alignment[2]};
		$reads{$alignment[0]}{'score'}{$alignment[2]} += $alignment[4];
		$reads{$alignment[0]}{'length'}{$alignment[2]} += length($alignment[9]);
		$reads{$alignment[0]}{'flag'}{$alignment[2]} = $alignment[1];

	} else {

		++$unmapReads;

	} #end of map/unmap

} #end while
close(SAM);

#print "Number of map pairs of reads:	", scalar(keys %reads), "\n";
#print "Number of unmap pairs of reads:	", $unmapReads/2, "\n";

print STDERR "Done loading mapping information, will do coverage calculation\n\n";

my %references;

foreach my $read (keys %reads) {

	my @ref = keys %{$reads{$read}{'reference'}};

	if ( (scalar (keys %{$reads{$read}{'reference'}})) == 1) {

		$references{$ref[0]}{'length'} +=  $reads{$read}{'length'}{$ref[0]};
		++$references{$ref[0]}{'count'};

	} elsif ( (scalar (keys %{$reads{$read}{'reference'}})) == 2) {

		if ( $reads{$read}{'reference'}{$ref[0]} > $reads{$read}{'reference'}{$ref[1]} ) {

			$references{$ref[0]}{'length'} +=  $reads{$read}{'length'}{$ref[0]};
			++$references{$ref[0]}{'count'};

		} elsif ( $reads{$read}{'reference'}{$ref[0]} < $reads{$read}{'reference'}{$ref[1]} ) {

			$references{$ref[1]}{'length'} +=  $reads{$read}{'length'}{$ref[1]};
			++$references{$ref[1]}{'count'};

		} elsif ( $reads{$read}{'reference'}{$ref[0]} == $reads{$read}{'reference'}{$ref[1]} ) {

			if ( $reads{$read}{'score'}{$ref[0]} > $reads{$read}{'score'}{$ref[1]} ) {

				$references{$ref[0]}{'length'} +=  $reads{$read}{'length'}{$ref[0]};
				++$references{$ref[0]}{'count'};

			} elsif ( $reads{$read}{'score'}{$ref[0]} < $reads{$read}{'score'}{$ref[1]} ) {

				$references{$ref[1]}{'length'} +=  $reads{$read}{'length'}{$ref[1]};
				++$references{$ref[1]}{'count'};

			} elsif ( $reads{$read}{'score'}{$ref[0]} == $reads{$read}{'score'}{$ref[1]} ) {

					$references{$ref[0]}{'length'} +=  $reads{$read}{'length'}{$ref[0]}/2;
					$references{$ref[0]}{'count'} += 0.5;
					$references{$ref[1]}{'length'} +=  $reads{$read}{'length'}{$ref[1]}/2;
					$references{$ref[1]}{'count'} += 0.5;

			} #end elsif


		} #end elsif

	} #end multiple references

} #end counting references

print STDERR "Done calculating coverages, will load size information\n\n";

my %sizes;

while( my $sizeLine = <SIZES>) {

	chomp($sizeLine);

	


} #end while leading sizes

print STDERR "Done loading size information, will print data\n\n";

foreach my $reference (sort {$a cmp $b} keys %references) {

	print COUNTS $reference , "\t", $references{$reference}{'count'}, "\t", $references{$reference}{'length'}, "\n";

} #end counting references

#foreach my $reference (sort {$a cmp $b} keys %singleEndMappedReferences) {
#
#	print COUNTSSINGLE $reference , "\t", $singleEndMappedReferences{$reference}, "\n";
#
#} #end counting references

exit;













