#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[0])) {

	print "Error, not gsnap input was entered. Please enter a gsnap output file for analysis\n";
	die "Missing arguments\n";

} #end if

open(GSNAP, $ARGV[0]);
open(SAM, ">$ARGV[0].sam");

my $firstRead = "YES";

while(my $firstLine = <GSNAP>) {

	chomp($firstLine);

	my @firstArray = split("\t", $firstLine);

	my $readName = $firstArray[3];

	my $seq = $firstArray[0];

	my $qual = $firstArray[2];

	my $alignmentLine = <GSNAP>;

	my $flag;

	if($firstRead eq "YES") { $flag = '99'; $firstRead = "NO"}
	else { $flag = '147'; $firstRead = "YES"}

	chomp($alignmentLine);

	my @secondArray = split ("\t", $alignmentLine);

	my @alignCoordinates = split("\\.\\.", $secondArray[1]);

	my $alignSeq = substr($seq, $alignCoordinates[0], ($alignCoordinates[1]-$alignCoordinates[0]+1) );

	my $alignQual = substr($qual, $alignCoordinates[0], ($alignCoordinates[1]-$alignCoordinates[0]+1) );

	my @location = split(":", $secondArray[2]);

	my $reference = $location[0];

	$reference =~ s/\+|\-//;

	my @coordinates = split("\\.\\.", $location[1]);

	my $end = $coordinates[0]-$coordinates[1];

	my @alignInfo = split(",", $secondArray[4]);

	my @mapq = split(":", $alignInfo[2]);

	print SAM join("\t", ($readName, $flag, $reference, "0", $mapq[1], $coordinates[0],"=", $coordinates[0], $end, $alignSeq, $alignQual)), "\n";

	while(my $alignmentPieceLine = <GSNAP>) {

		if($alignmentPieceLine =~ m/^$/) {last}

		chomp($alignmentPieceLine);


	my @secondExonArray = split ("\t", $alignmentPieceLine);

	my @exonAlignCoordinates = split("\\.\\.", $secondExonArray[1]);

	my $exonAlignSeq = substr($seq, $exonAlignCoordinates[0], ($exonAlignCoordinates[1]-$exonAlignCoordinates[0]+1) );

	my $exonAlignQual = substr($qual, $exonAlignCoordinates[0], ($exonAlignCoordinates[1]-$exonAlignCoordinates[0]+1) );

	my @exonLocation = split(":", $secondExonArray[2]);

	my @exonCoordinates = split("\\.\\.", $exonLocation[1]);

	my $exonEnd = $exonCoordinates[0]-$exonCoordinates[1];

		print SAM join("\t", ($readName, $flag, $reference, "0", $mapq[1], $exonCoordinates[0],"=", $exonCoordinates[0], $exonEnd, $exonAlignSeq, $exonAlignQual)), "\n";

	} #end alignment

} #end while file

close(GSNAP);
close(SAM);

exit;
